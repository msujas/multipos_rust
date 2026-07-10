use std::{path::Path, time::Instant};
use clap::Parser;
use multiposrust::MultiFile;

use multiposrust::params::ParamsFF;

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

    let chibins = 50; //not used for calculating the flat field, but needed to make MultiFile object
    let tmp: String;
    let maskfile = match maskfileo{
        None => None,
        Some(f) => {tmp=f;
                    Some(Path::new(&tmp))}
    };
    let now = Instant::now();
    //let mf = MultiFile::build(&cbfdir, &ponidir, tthmin, tthmax, tthbins, chimin, chimax, chibins, pfactor, maskfile, maskdir);
    let mf = MultiFile::buildinterpolate(&cbfdir, &ponidir, tthmin, tthmax, tthbins, chimin, chimax, chibins, pfactor,
         maskfile, maskdir, ponipattern, ymotor, zmotor, saveponis);
    mf.calculateflatfield(&cbfdir,ffmin,ffmax);
    let elapsed = now.elapsed();
    println!("program took {} s", elapsed.as_secs());
}