use std::{path::Path, time::Instant};

use crate::params::ParamsFF;
use clap::Parser;
use multipos_rust::functions::MultiFile;

mod functions;
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

    let chibins = 50;
    let tmp: String;
    let maskfile = match maskfileo{
        None => None,
        Some(f) => {tmp=f;
                    Some(Path::new(&tmp))}
    };
    let now = Instant::now();
    let mf = MultiFile::build(&cbfdir, &ponidir, tthmin, tthmax, tthbins, chimin, chimax, chibins, pfactor, maskfile, maskdir);
    mf.calculateflatfield();
    let elapsed = now.elapsed();
    println!("program took {}", elapsed.as_secs());
}