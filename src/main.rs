use clap::Parser;
use cryiorust::cbf;
use integrustio::integrator::PatternType;
use multipos_rust::{ImagePoni, MultiFile};
use std::path::Path;
use crate::params::Params;
use glob::glob;
use std::ffi::OsStr;

mod params;

fn main() {
    
    let ap = Params::parse();
    let tthmin = ap.tthmin;
    let tthmax = ap.tthmax;
    let tthbins = ap.tthbins;
    let chimin = ap.chimin;
    let chimax = ap.chimax;
    let chibins = ap.chibins;
    let pfactor = ap.pfactor;
    let cbfdir = &ap.cbfdir;
    let ponidir = &ap.ponidir;
     


    let cbfpath = Path::new("D:/beamlineData/April2026/multipositions3/Si/Si_dty260.00_dtz117.50_001_0001p.cbf");
    let ponipath = Path::new("D:/beamlineData/April2026/multipositions3/Si/poni/Si_dty260.00_dtz117.50_001_0001p.poni");
    let ip = ImagePoni::build(ponipath, cbfpath);
    let d = ip.integrate(0.75, 58.,5000, 2.,358., 356,0.85);
    let p: integrustio::integrator::PatternType = d.data;
    let cbfbase = cbfpath.file_name().unwrap().to_str().unwrap().replace(".cbf", "");
    let ponibase = ponipath.file_name().unwrap().to_str().unwrap().replace(".poni", "");
    println!("{cbfbase}");
    println!("{ponibase}");
    let b = cbfbase == ponibase ;
    println!("{b}");

    let mf = MultiFile::build(cbfdir, ponidir, tthmin, tthmax, tthbins, chimin, chimax, chibins, pfactor);
    p.store("D:/beamlineData/April2026/multipositions3/cake.edf", None, None).unwrap();
    let cake = match p{
        PatternType::Cake(c) => c,
        _=>panic!("some issue with getting cake image"),
    };
    let d1 = cake.cake.dim1();
    let d2 = cake.cake.dim2();
    print!("d1: {d1}, d2 {d2}");
}
