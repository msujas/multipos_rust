use std::process::exit;
use std::{path::Path, time::Instant};
use clap::Parser;
use multiposrust::MultiFile;
use multiposrust::imagereader::ImageFormat;

use crate::params::ParamsFF;

mod params;


fn main(){
    let ap = ParamsFF::parse();
    let cbfdir = ap.cbfdir;
    let tthmin = ap.tthmin;
    let tthmax = ap.tthmax;
    let tthbins = ap.tthbins;
    let chimin = ap.chimin;
    let chimax = ap.chimax;
    let pfactor = ap.pfactor;
    let ponidir = ap.ponidir;
    let maskfileo = ap.maskfile;
    let maskdir = ap.maskdir;
    let ffmin = ap.ffmin;
    let ffmax = ap.ffmax;
    let ponipattern = &ap.ponipattern;
    let ymotor = &ap.ymotor;
    let zmotor = &ap.zmotor;
    let saveponis = ap.saveponis;
    let unit = ap.unit;
    let fluxentry = ap.fluxentry;
    let fileextension = ap.fileextension;

    let imageformat = match ImageFormat::fromextension(fileextension){
        Ok(i) => i,
        Err(_e) => {eprintln!("Exiting"); exit(2)},
    };

    let chibins = 50; //not used for calculating the flat field, but needed to make MultiFile object
    let tmp: String;
    let maskfile = match maskfileo{
        None => None,
        Some(f) => {tmp=f;
                    Some(Path::new(&tmp))}
    };
    let now = Instant::now();
    let mf = match MultiFile::buildinterpolate(&cbfdir, &ponidir, imageformat, tthmin, tthmax, tthbins, chimin, chimax, 
        chibins, pfactor,maskfile, maskdir, ponipattern, ymotor, zmotor, saveponis, Some(&unit), None, fluxentry){
            Ok(m) => m,
            Err(_e) => {eprintln!("error. Exiting"); exit(1)}
         };
    mf.calculateflatfield(&cbfdir,ffmin,ffmax);
    let elapsed = now.elapsed();
    println!("program took {} s", elapsed.as_secs());
}