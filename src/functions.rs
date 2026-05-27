use cryiorust::{cbf::Cbf, frame::{Array, Frame}, poni::{DetectorConfig, Poni}};
use integrustio::integrator::{Cake, Integrator};
use core::f64;
use std::{cmp::Ordering, fs::File, io::{self, Write}, path::{Path, PathBuf}, sync::Arc, vec};
use glob::glob;
pub struct ImagePoni{
    pub poni:Poni,
    pub cbf:Cbf,
}
impl ImagePoni {
    pub fn build(ponifile:&Path, cbffile: &Path, dc: Option<Arc<DetectorConfig>>) -> ImagePoni{
        let poni = Poni::open(ponifile, dc).unwrap() ;
        let cbf = Cbf::open(cbffile).unwrap();
        //let dirname = cbffile.parent().unwrap();
        ImagePoni { poni, cbf }
    }
    pub fn integrate(self, tthmin:f64, tthmax:f64, tthbins:usize, chimin:f64, chimax: f64, 
                chibins: usize, pfactor: f64)->Cake{
        let mut i = Integrator::new();
        i.set_poni(self.poni);
        i.set_radial_bins(tthbins);
        i.set_azimuthal_bins(chibins);
        i.set_polarization(pfactor);
        i.set_azimuthal_range(Some((chimin,chimax)));
        i.set_radial_range(Some((tthmin,tthmax)));
        i.init(self.cbf.array());
        let (cake, _) = i.integrate_cake(self.cbf.array()).unwrap();
        cake
    }

    pub fn get_cake(self,tthmin:f64, tthmax:f64, tthbins:usize, chimin:f64, chimax: f64, 
                chibins: usize, pfactor: f64, cakedir:&String)->Cake{
        let mut fname = self.cbf.name().to_string();
        fname.push_str(".edf");
        let cake = self.integrate(tthmin, tthmax, tthbins, chimin, chimax, chibins, pfactor);
        
        if cakedir != ""{
            
            let cakefile = format!("{}/{}",cakedir,fname);
            cake.store(cakefile, None).unwrap();
        }
        cake
    }
}


pub(crate) struct MultiFile{
    ilist :Vec<ImagePoni>,
    tthmin: f64,
    tthmax: f64,
    tthbins: usize,
    chimin: f64,
    chimax: f64,
    chibins: usize,
    pfactor:f64
}

impl MultiFile{

    pub fn build(cbfdir:&String, ponidir: &String, tthmin:f64, tthmax:f64, tthbins:usize, chimin: f64, chimax: f64, 
                chibins:usize, pfactor:f64) -> MultiFile {
                    let mut dc = None;
                    let pattern = format!("{cbfdir}/*.cbf");
                    let cbffiles = glob(&pattern).unwrap();
                    let mut cbffile: Option<Arc<PathBuf>>;
                    
                    let mut ilist:Vec<ImagePoni> = Vec::new();
                    for fresult in cbffiles{
                        cbffile = Some(Arc::new(fresult.unwrap()));
                        let cbfclone = cbffile.clone();
                        
                        let basename = cbfclone.unwrap().file_name().unwrap().to_str().unwrap().replace(".cbf", "");
                        let ponifiles = glob(&format!("{ponidir}/*.poni")).unwrap();
                        for presult in ponifiles{
                            let ponifile = presult.unwrap();
                            let pbasefile = ponifile.file_name().unwrap().to_str().unwrap().replace(".poni","");
                            
                            if basename == pbasefile{
                                
                                println!("{:?}, {:?}", cbffile.clone().unwrap(),ponifile);
                                let ip = match dc{
                                    None => {let ip = ImagePoni::build(&ponifile,&cbffile.clone().unwrap(), None);
                                    dc = ip.poni.detector_config.clone();ip},
                                    Some(ref dc) => ImagePoni::build(&ponifile,&cbffile.clone().unwrap(), Some(dc.clone())),
                                };
                                ilist.push(ip);
                                break;
                            }
                        }
                    }
                    if ilist.len() < 2{
                        let nitems = ilist.len();
                        panic!("build function requires at least 2 pairs of ponis and cbfs, found {nitems}. cbf and poni files must match in the base name")
                    }
                    //let mut it = cbffiles.into_iter();
                    
                    MultiFile { ilist,  tthmin, tthmax, tthbins, chimin, chimax, chibins, pfactor }
                }
    pub fn integrate_all(self, cakedir: &String)->Vec<Cake>{
        let mut cakes :Vec<Cake> = Vec::new();
        let mut count = 0;
        println!("integrating images");
        for ip in self.ilist{
            print!("{count}, ");
            io::stdout().flush().unwrap();
            cakes.push(ip.get_cake(self.tthmin, self.tthmax, self.tthbins, self.chimin, self.chimax, self.chibins, self.pfactor, cakedir));
            count += 1;
        };
        cakes
    }

    pub fn average_cakes(self, medianfilter:f64, cakedir: &String, avdir: &String){
        let cakes = self.integrate_all(cakedir);
        println!("\naveraging cakes");
        //let mut medianvec: Vec<f64> = Vec::new();
        let c0 = &cakes[0];
        let rpos = &c0.radial_positions;
        let azpos = &c0.azimuthal_positions;
        let chisize = c0.cake.dim1();
        let radsize = c0.cake.dim2();
        
        let mut avvec : Vec<f64> = vec![0.;c0.cake.len()];
        let mut vec1d : Vec<f64> = vec![0.; radsize];
        let mut div1d : Vec<f64> = vec![0.; radsize];
        let mut sigma: Vec<f64> = vec![0.;radsize];
        let datalen = c0.cake.data().len();
        println!("dim1: {chisize}");
        println!("dim2: {radsize}");
        println!("array size {datalen}");
        for i in 0..datalen{
            let index1d = i%radsize;
            let mut atemp: Vec<f64> = Vec::new();
            for c in &cakes{
                let item = c.cake.data()[i];
                if item > 0.{
                    atemp.push(item);
                }
            }
            let mut div: f64 = 0.;
            let mut intensity :f64 = 0.;
            let med = getmedian(&atemp);
            for val in atemp{
                if !(val > med*medianfilter) & !(val < med/medianfilter){
                    intensity += val;
                    div += 1.;
                }
            }
            if div > 0.1 {
                intensity = intensity/div;
            }
            else {
                intensity = 0.;
            };
            avvec[i] = intensity;
            
            if intensity > 0.{
                vec1d[index1d] += intensity;
                div1d[index1d] += 1.;
            }
        };
        println!("averaging 1d");
        for j in 0..vec1d.len(){
            let intj = vec1d[j];
            vec1d[j] = intj/div1d[j];
            sigma[j] = intj.powf(0.5)/div1d[j]; //approximation of error, square root intensity divide by number of values
        }
        let a:Array = Array::with_data(chisize,radsize, avvec);
        println!("{vec1d:?}");
        

        
        let mut newcake:Cake = Default::default(); 
        newcake.cake = a; 
        newcake.radial_positions = rpos.clone();
        newcake.azimuthal_positions = azpos.clone();
        newcake.radial.intensity = vec1d;   
        newcake.radial.sigma = sigma;
        newcake.radial.positions = rpos.clone();
        let fnameav = format!("{avdir}/avcake.edf");
        newcake.store(fnameav, None).unwrap();
    }
}

fn save1d(fname:String, vec1d: Vec<f64>){
    let mut outstring = String::new();
    for item in vec1d{
        outstring = outstring + &String::from(format!("{item}\n"));
    }
    let mut file = File::create(fname).unwrap();
    file.write(outstring.as_bytes()).unwrap();    
}

fn getmedian(medvec:&Vec<f64>)->f64{
    //let svec = sortvec(medvec);
    let mut svec = medvec.clone();
    svec.sort_by(cmpf64);
    let vlen = svec.len();
    let pos : usize = vlen/2;
    if vlen == 0 {
        return 0.;
    }
    if vlen%2 == 0 {
        return (svec[pos-1] + svec[pos])/2.;
    }
    return svec[pos]
}

fn cmpf64(a:&f64,b:&f64)->Ordering{
    if a > b    {
        return Ordering::Greater
    }
    else if  a < b{
        return  Ordering::Less;
    }
    return Ordering::Equal;
    }
    


#[cfg(test)]
mod tests{
    use super::*;
    fn floatcompare(f1:f64, f2:f64)->bool{
        let b = (f1 > f2 - 0.01) & (f1 < f2+0.01);
        b
    }
    #[test]
    fn medtest(){
        let v = vec![1.3, 5.7, 2.2, 9.8, 7.3];
        let v2 = vec![53.4, 6.9, 8.2, 1.1, 17.3, 24.1];
        let m1 = getmedian(&v);
        let m2 = getmedian(&v2);
        let b1 = floatcompare(m1, 5.7);
        let b2 = floatcompare(m2,12.75);
        println!("{m1}");
        println!("{m2}");
        assert_eq!(b1, true);
        assert_eq!(b2, true);
    }
}