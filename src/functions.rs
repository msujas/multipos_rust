use cryiorust::{cbf::Cbf, edf::{self, Edf}, frame::{ Array, Frame, HeaderEntry::Float}, poni::{DetectorConfig, Poni}};
use fluosubtraction_rust::functions::fluosub_curvefit;
use integrustio::integrator::{Cake, Integrator};
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator,  ParallelIterator};
use core::f64;
use std::{cmp::Ordering, fs::{File, create_dir}, io::{self, Write}, path::{Path, PathBuf}, sync::Arc, vec};
use glob::{ glob};

pub struct ImagePoni{
    pub poni:Poni,
    pub cbf:Cbf,
}
impl ImagePoni {
    pub fn build(ponifile:&Path, cbffile: &Path, dc: Option<Arc<DetectorConfig>>, mask: Option<&Array>) -> ImagePoni{
        let poni = Poni::open(ponifile, dc).unwrap() ;
        let mut cbf = Cbf::open(cbffile).unwrap();

        let flux = match cbf.header().get("# Flux "){
            Some(Float(f64)) => f64.clone(),
            _ => panic!("couldn\'t find flux for {cbffile:?}"),
        };
        if let Some(mask) = mask{
            for (i , m) in cbf.array_mut().data_mut().iter_mut().zip(mask.data().iter()){
                *i = *i/ flux;
                if *m > 0. {
                    *i = -1.;
                }
            }
        }
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
        i.set_solid_angle(true, true);
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


pub struct MultiFile{
    ilist :Vec<ImagePoni>,
    tthmin: f64,
    tthmax: f64,
    tthbins: usize,
    chimin: f64,
    chimax: f64,
    chibins: usize,
    pfactor:f64,

}

impl MultiFile{

    pub fn build(cbfdir:&String, ponidir: &String, tthmin:f64, tthmax:f64, tthbins:usize, chimin: f64, chimax: f64, 
                chibins:usize, pfactor:f64, maskfile: Option<&Path>, maskdir: Option<String>) -> MultiFile {
                    let mut dc = None;
                    let pattern = format!("{cbfdir}/*.cbf");
                    let cbffiles = glob(&pattern).unwrap();
                    let mut cbffile: Option<Arc<PathBuf>>;
                    let binding: edf::Edf;
                    let mask: Option<&Array> = match maskfile {
                        None => None,
                        Some(f) => {binding = edf::Edf::open(f).unwrap();
                                          Some(binding.array())},
                    };
                    let mut usedmask : Option<&Array> = mask.clone();
                    let mut ilist:Vec<ImagePoni> = Vec::new();
                    let mut mbinding:Edf;
                    for fresult in cbffiles{
                        cbffile = Some(Arc::new(fresult.unwrap()));
                        let cbfclone = cbffile.clone();
                        
                        let basename = cbfclone.unwrap().file_name().unwrap().to_str().unwrap().replace(".cbf", "");
                        let ponifiles = glob(&format!("{ponidir}/*.poni")).unwrap();
                        if let Some(md) = &maskdir{
                            let maskfiles = glob(&format!("{md}/*.edf")).unwrap();
                            for mresult in maskfiles{
                                let m = mresult.unwrap();
                                let mut breakloop:bool = false;
                                let mbase = m.file_name().unwrap().to_str().unwrap().replace(".edf", "");
                                
                                usedmask = match basename.find(&mbase){
                                    None => mask,
                                    Some(_s) => {mbinding = Edf::open(m).unwrap();
                                    breakloop=true;
                                    Some(mbinding.array())},
                                };
                                if breakloop{
                                    println!("{basename}.cbf using mask {mbase}.edf");
                                    break;
                                }
                                
                                
                            }
                        }
                        for presult in ponifiles{
                            let ponifile = presult.unwrap();
                            let pbasefile = ponifile.file_name().unwrap().to_str().unwrap().replace(".poni","");
                            
                            if basename == pbasefile{
                                
                                println!("{:?}, {:?}", cbffile.clone().unwrap(),ponifile);
                                let ip = match dc{
                                    None => {let ip = ImagePoni::build(&ponifile,&cbffile.clone().unwrap(), None, usedmask);
                                    dc = ip.poni.detector_config.clone();ip},
                                    Some(ref dc) => ImagePoni::build(&ponifile,&cbffile.clone().unwrap(),
                                     Some(dc.clone()), usedmask),
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
        //let mut cakes :Vec<Cake> = Vec::new(); //vec![Default::default(); self.ilist.len()];
        
        println!("integrating images");

        let cakes: Vec<Cake> = self.ilist.into_par_iter()
        .enumerate()
        .map(|(i,ip)|{
            print!("{i}, ");
            io::stdout().flush().unwrap();
            ip.get_cake(self.tthmin, self.tthmax, 
                    self.tthbins, self.chimin, self.chimax, self.chibins, self.pfactor, cakedir)
        }).collect();
     
        cakes
    }

    pub fn average_cakes(self, medianfilter:f64, cakedir: &String, avdir: &String, cakemaskfile: Option<String>)->Cake{
        let cakes = self.integrate_all(cakedir);
        println!("\naveraging cakes");
        let c0 = &cakes[0];
        let rpos = &c0.radial_positions;
        let azpos = &c0.azimuthal_positions;
        let chisize = c0.cake.dim1();
        let radsize = c0.cake.dim2();
        let _ = create_dir(&avdir);

        let mut avvec : Vec<f64> = vec![0.;c0.cake.len()];
        let mut vec1d : Vec<f64> = vec![0.; radsize];
        let mut div1d : Vec<f64> = vec![0.; radsize];
        let mut sigma: Vec<f64> = vec![0.;radsize];
        let datalen = c0.cake.data().len();
        println!("dim1: {chisize}");
        println!("dim2: {radsize}");
        println!("array size {datalen}");

        let tmp :Edf;
        let cakemask = match cakemaskfile{
            None => None,
            Some(s)  if Edf::open(s.clone()).unwrap().array().data().len() == datalen =>{ 
                tmp = Edf::open(s).unwrap();
                Some(tmp.array())}
            _ => {println!("mismatch in cake mask and data length. Ignoring mask");
                None}  
        };

        for i in 0..datalen{
            let index1d = i%radsize;
            let mut atemp: Vec<f64> = Vec::new();
            for c in &cakes{
                let item = c.cake.data()[i];
                if item > 0.{
                    if let Some(cakemask) = cakemask {
                        if cakemask.data()[i] > 0.{
                        continue;}
                    }
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
        //println!("{vec1d:?}");
        
        let fname1d  = format!("{avdir}/av1d_2.xye");
        save1d(fname1d, rpos, &vec1d, Some(&sigma));
        let mut newcake:Cake = Default::default(); 
        newcake.cake = a; 
        newcake.radial_positions = rpos.clone();
        newcake.azimuthal_positions = azpos.clone();
        newcake.radial.intensity = vec1d;   
        newcake.radial.sigma = sigma;
        newcake.radial.positions = rpos.clone();
        let fnameav = format!("{avdir}/avcake.edf");
        newcake.store(fnameav, None).unwrap();
        
        let av1d_alt = cakeav(&cakes, cakemask, medianfilter);
        let fname1d_alt = format!("{avdir}/av1d.xy");
        save1d(fname1d_alt, rpos, &av1d_alt, None);

        newcake
    }

    pub fn integrate_fluosub(self, medianfilter:f64, cakedir: &String, avdir: &String, cakemaskfile: Option<String>, 
                        fluo_k0:f64, tthindex:usize)->Cake{
        let pfactor = self.pfactor;
        let cake = self.average_cakes(medianfilter, cakedir, avdir, cakemaskfile);
        let newcake = fluosub_curvefit(fluo_k0, cake, pfactor, tthindex);
        let fsdir = format!("{avdir}fluoSub");
        let _ = create_dir(&fsdir);
        let fname = format!("{}/avcake.edf", &fsdir);
        let fname1d = format!("{}/avcake.xye", &fsdir);
        newcake.store(fname, None).unwrap();
        let av1d = &newcake.radial.intensity;
        let tth = &newcake.radial.positions.to_vec();
        let sigma = &newcake.radial.sigma;
        save1d(fname1d, tth, av1d, Some(sigma));
        newcake
    }
}



fn cakeav(cakelist: &Vec<Cake>, cakemask: Option<&Array>, medianfilter:f64)-> Vec<f64>{

    let c0 = &cakelist[0];
    let chisize = c0.cake.dim1();
    let tthsize= c0.cake.dim2();
    let mut av1d : Vec<f64> = vec![0.; tthsize];
    let mut divvec: Vec<f64> = vec![0.;tthsize];
    let mut index: usize;

    for c in cakelist{
        for i in 0..tthsize{
            let mut vtemp : Vec<f64> = Vec::new();
            for j in 0..chisize{
                index = i + j*tthsize;
                let value = c.cake.data()[index];
                if c.cake.data()[index] > 0.{
                    if let Some(cakemask)=cakemask{
                        if cakemask.data()[index] > 0.1{
                            continue;
                        }
                    }
                    vtemp.push(value);
                }
            }
            let med = getmedian(&vtemp);
            let mut intensity = 0.;
            let mut div = 0.;
            for item in vtemp{
                if (item < med*medianfilter) & (item > med/medianfilter){
                    intensity += item;
                    div += 1.;
                }
            }
            if div > 0. {
                av1d[i] += intensity/div;
                divvec[i] += 1.;
            }    
        }
    }
    for (x, d) in av1d.iter_mut().zip(divvec.iter_mut()){
        if *d > 0. {
        *x = *x/ *d;
        }
    }
av1d
}

fn save1d(fname:String, tthrange: &Vec<f64>, vec1d: &Vec<f64>, sigma : Option<&Vec<f64>>){
    let mut outstring = String::new();
    //for (x,y ) in  tthrange.iter().zip(vec1d.iter()){
    let mut x:f64;
    let mut y:f64;
    let mut e:f64;
    let dosig:bool = match sigma  {
        None => false,
        Some(_s) => true
    };
    for i in 0..tthrange.len(){
        x = tthrange[i];
        y=vec1d[i];
        outstring = outstring + &String::from(format!("{x:.6} {y:.6}"));
        if dosig{
            e = sigma.unwrap()[i];
            outstring = outstring + &String::from(format!(" {e:.6}"));
            }
        outstring = outstring + &String::from("\n");
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