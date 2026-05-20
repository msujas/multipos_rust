use cryiorust::{cbf::Cbf, frame::{Array, Frame}, poni::Poni};
use integrustio::integrator::{Cake, Diffractogram, Integrable, IntegrationType, Integrator, Pattern, PatternType};
use core::f64;
use std::{cmp::Ordering, ops::Mul, path::Path, vec};
use glob::glob;
pub struct ImagePoni{
    pub poni:Poni,
    pub cbf:Cbf,
}
impl ImagePoni {
    pub fn build(ponipath:&Path, cbfpath: &Path) -> ImagePoni{
        let poni = Poni::read_file(ponipath).unwrap() ;
        let cbf = Cbf::read_file(cbfpath).unwrap();
        ImagePoni { poni, cbf }
    }
    pub fn integrate(self, tthmin:f64, tthmax:f64, tthbins:usize, chimin:f64, chimax: f64, 
                chibins: usize, pfactor: f64)->Diffractogram{
        let ranges = [tthmin,tthmax];
        let chirange = [chimin,chimax];
        let data = Integrable{
            array: self.cbf.array(),
            radial_range: &ranges,
            azimuthal_range: &chirange,
            integration_type: IntegrationType::Cake
        };
        let mut i = Integrator::new();
        i.set_poni(self.poni);
        i.set_radial_bins(tthbins);
        i.set_azimuthal_bins(chibins);
        i.set_polarization(pfactor);
        i.init(self.cbf.array());
        let d = i.integrate(&data).unwrap();
        d
    }

    pub fn get_cake(self,tthmin:f64, tthmax:f64, tthbins:usize, chimin:f64, chimax: f64, 
                chibins: usize, pfactor: f64)->Cake{
        let d = self.integrate(tthmin, tthmax, tthbins, chimin, chimax, chibins, pfactor);
        let p: integrustio::integrator::PatternType = d.data;
        let cake = match p{
            PatternType::Cake(c) => c,
            _=>panic!("some issue with getting cake image"),
        };
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
    pfactor:f64
}

impl MultiFile{

    pub fn build(cbfdir:&String, ponidir: &String, tthmin:f64, tthmax:f64, tthbins:usize, chimin: f64, chimax: f64, 
                chibins:usize, pfactor:f64) -> MultiFile {
                    let pattern = format!("{cbfdir}/*.cbf");
                    let cbffiles = glob(&pattern).unwrap();
                    let mut ilist:Vec<ImagePoni> = Vec::new();
                    for fresult in cbffiles{
                        let cbffile = fresult.unwrap();
                        let basename = cbffile.file_name().unwrap().to_str().unwrap().replace(".cbf", "");
                        let ponifiles = glob(&format!("{ponidir}/*.poni")).unwrap();
                        for presult in ponifiles{
                            let ponifile = presult.unwrap();
                            let pbasefile = ponifile.file_name().unwrap().to_str().unwrap().replace(".poni","");
                            
                            if basename == pbasefile{
                                println!("{cbffile:?}, {ponifile:?}");
                                ilist.push(ImagePoni::build(&ponifile, &cbffile));
                                break;
                            }
                        }
                    }
                    if ilist.len() < 2{
                        let nitems = ilist.len();
                        panic!("build function requires at least 2 pairs of ponis and cbfs, found {nitems}. cbf and poni files must match in the base name")
                    }
                    MultiFile { ilist, tthmin, tthmax, tthbins, chimin, chimax, chibins, pfactor }
                }
    pub fn integrate_all(self)->Vec<Cake>{
        let mut cakes :Vec<Cake> = Vec::new();
        for ip in self.ilist{
            cakes.push(ip.get_cake(self.tthmin, self.tthmax, self.tthbins, self.chimin, self.chimax, self.chibins, self.pfactor));
        };
        cakes
    }

    pub fn average_cakes(self, medianfilter:f64){
        let cakes = self.integrate_all();
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
        for i in 0..c0.cake.data().len(){
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
            intensity = intensity/div;
            avvec[i] = intensity;
            
            if intensity > 0.{
                vec1d[index1d] += intensity;
                div1d[index1d] += 1.;
            }
        };

        for j in 0..vec1d.len(){
            let intj = vec1d[j];
            vec1d[j] = intj/div1d[j];
            sigma[j] = intj.powf(0.5)/div1d[j]; //approximation of error, square root intensity divide by number of values
        }
        let a:Array = Array::with_data(chisize,radsize, avvec);
        /*
        let newp:Pattern = Pattern{
            positions:rpos.clone(),
            intensity: vec1d,
            sigma: sigma,
            new_line: '\n'
        };
        
        
        let newcake:Cake = Cake { radial_positions:rpos.clone(), azimuthal_positions: azpos.clone(), cake: a, radial: newp };
         */
    }
}


fn getmedian(medvec:&Vec<f64>)->f64{
    //let svec = sortvec(medvec);
    let mut svec = medvec.clone();
    svec.sort_by(cmpf64);
    let vlen = svec.len();
    let pos : usize = vlen/2;
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