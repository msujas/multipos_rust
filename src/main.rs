use std::path::Path;
use std::{fs::create_dir, time::Instant};
use clap::Parser;
use crate::functions::MultiFile;
use crate::params::Params;

mod functions;
mod params;

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
    let tmp:String;
    let maskfile = match maskfilestr {
        None => None,
        Some(s) => {tmp = s;
                        Some(Path::new(&tmp))},
    };

    let mut cakedir: String = String::from("");

    if savecakes {
        cakedir.push_str(&format!("{cbfdir}/{subdir}"));
        let _ = create_dir(&cakedir);
    }

    let mf = MultiFile::build(cbfdir, ponidir, tthmin, tthmax, tthbins, chimin, chimax, chibins, pfactor, maskfile);
    let e1 = now.elapsed();
    println!("loading cbfs and ponis took {} s", e1.as_secs());
    let avdir = format!("{cbfdir}/{subdir}");
    mf.average_cakes(4., &cakedir, &avdir);
    let elapsed =  now.elapsed();
    println!("");
    println!("program took {} s", elapsed.as_secs());
}
