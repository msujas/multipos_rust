use cryiorust::{cbf::Cbf, frame::{Frame}, poni::Poni};
use integrustio::integrator::{Cake, Diffractogram, Integrable, IntegrationType, Integrator, PatternType};
use core::f64;
use std::{cmp::Ordering, path::Path, vec};

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
    pub fn integrate_all(self)->Vec<Cake>{
        let mut cakes :Vec<Cake> = Vec::new();
        for ip in self.ilist{
            cakes.push(ip.get_cake(self.tthmin, self.tthmax, self.tthbins, self.chimin, self.chimax, self.chibins, self.pfactor));
        };
        cakes
    }

    pub fn average_cakes(self, medianfilter:i32){
        let cakes = self.integrate_all();
        let mut medianvec: Vec<f64> = Vec::new();
        for cake in cakes{
            let d = cake.cake.data();
        }
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
    
fn floatcompare(f1:f64, f2:f64)->bool{
    let b = (f1 > f2 - 0.01) & (f1 < f2+0.01);
    b
}

#[cfg(test)]
mod tests{
    use super::*;

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