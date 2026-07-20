use std::path::Path;
use std::process::exit;
use std::{time::Instant};
use clap::Parser;
use multiposrust::MultiFile;
use multiposrust::params::Params;


fn main() {
    let now = Instant::now();
    
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
    let savecakes = ap.savecakes;
    let subdir = ap.cakesubdir;
    let maskfilestr = ap.maskfile;
    let cakemaskfile = ap.cakemaskfile;
    let maskdir = ap.maskdir;
    let ymotor = ap.ymotor;
    let zmotor = ap.zmotor;
    let ponipattern = ap.ponipattern;
    let fluosub = ap.fluosub;
    let fluok0 = ap.fluok0;
    let saveponis = ap.saveponis;
    let unit = ap.unit;

    let tmp:String;
    let maskfile: Option<&Path> = match maskfilestr {
        None => None,
        Some(s) => {tmp = s;
                        Some(Path::new(&tmp))},
    };

    
    let cakedir = match savecakes {
        true => Some(&format!("{cbfdir}/{subdir}")),
        false => None,        
    };
    
    let mf = match MultiFile::buildinterpolate(cbfdir, ponidir, tthmin, tthmax, tthbins, chimin, chimax, chibins, 
        pfactor, maskfile,maskdir, &ponipattern, &ymotor, &zmotor, saveponis, Some(&unit)){
            Ok(m) => m,
            Err(_e) => {eprintln!("exiting"); exit(1)}
        };
    let e1 = now.elapsed();
    println!("loading cbfs and ponis took {} s", e1.as_secs());
    let avdir = format!("{cbfdir}/{subdir}");
    if fluosub{
        mf.integrate_fluosub(4., cakedir, &avdir, cakemaskfile, fluok0, tthbins*96/100);
    }
    else {
        mf.average_cakes(4., cakedir, &avdir, cakemaskfile);
    };
    let elapsed =  now.elapsed();
    println!("");
    println!("program took {} s", elapsed.as_secs());
}
