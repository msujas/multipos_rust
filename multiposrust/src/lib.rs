use chrono::{DateTime, Datelike, Local};
use cryiorust::{cbf::Cbf, edf::{self, Edf}, frame::{ Array, Frame, Header, HeaderEntry::{self, Float}}, 
poni::{DetectorConfig, Poni}};
use fluosubtraction_rust::functions::fluosub_curvefit;
use integrustio::{ geometry::{IntoGeometry, Units::{self}}, integrator::{Cake, Integrator, KEY_BUBBLE_MADE}};
use rayon::iter::{IndexedParallelIterator, IntoParallelIterator,  ParallelIterator};
use core::f64;
use std::{borrow::Cow,  f64::consts::PI, fs::{File, create_dir}, io::{self, BufWriter,  Write}, path::{Path, PathBuf}, 
sync::Arc, vec};
use glob::{ glob};
use functions::{save1d, cakeav, getmedian, closestindexordered,getyz};

use crate::{functions::yzcompare, poniinterpolator::{Interpolators, PoniList}};

#[derive(Debug)]
pub struct BuildError;

pub mod params;
mod functions;
pub mod poniinterpolator;
pub struct ImagePoni{
    pub poni:Poni,
    pub cbf:Cbf,
}

fn buildip(poni:Poni, cbffile: &Path,mask: Option<&Array>) -> Result<ImagePoni, BuildError>{
        let mut cbf = Cbf::open(cbffile).expect(&format!("couldn't open file: {cbffile:?}"));
        let flux = match cbf.header().get("# Flux "){
            Some(Float(f64)) => f64.clone(),
            _ => {eprintln!("couldn\'t find flux for {cbffile:?}");
                return Err(BuildError)},
        };
        if let Some(mask) = mask{
            for (i , m) in cbf.array_mut().data_mut().iter_mut().zip(mask.data().iter()){
                *i = *i/ flux;
                if *m > 0. {
                    *i = -1.;
                }
            }
        }
        Ok(ImagePoni { poni, cbf })
    }

impl ImagePoni {
    pub fn build(ponifile:&Path, cbffile: &Path, dc: Option<Arc<DetectorConfig>>, mask: Option<&Array>) -> Result<ImagePoni, BuildError>{
        let poni = Poni::open(ponifile, dc).unwrap() ;
        buildip(poni, cbffile, mask)
    }

    pub fn buildfromponi(poni:Poni, cbffile: &Path,mask: Option<&Array>) -> Result<ImagePoni, BuildError>{
        buildip(poni, cbffile, mask)
    }

    pub fn integrate(self, tthmin:f64, tthmax:f64, tthbins:usize, chimin:f64, chimax: f64, 
                chibins: usize, pfactor: f64, units: &str)->Cake{
        let units = strtounits(Some(units));
        let mut i = Integrator::new();
        i.set_poni(self.poni);
        i.set_radial_bins(tthbins);
        i.set_units(units);
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
                chibins: usize, pfactor: f64, cakedir:Option<&String>, units:&str )->Cake{
        let mut fname = self.cbf.name().to_string();
        fname.push_str(".edf");
        let cake = self.integrate(tthmin, tthmax, tthbins, chimin, chimax, chibins, pfactor, units);
        
        if let Some(cakedir) = cakedir{
            let cakefile = format!("{}/{}",cakedir,fname);
            println!("saving individual cake to {}",&cakefile);
            cake.store(&cakefile, None).unwrap();
        }
        cake
    }
}

fn strtounits(unitstr:Option<&str>)->Units{
    let units = match unitstr{
        None => Units::TwoTheta,
        Some(s) => {match s{
            "TwoTheta" => Units::TwoTheta,
            "2theta" => Units::TwoTheta,
            "QA" => Units::QA,
            "Qnm" => Units::Qnm,
            _ => {println!("couldn't interpret unit string, defaulting to 2theta");Units::TwoTheta},
        }}
    };
    units
}

fn optiontostr(unito:Option<&str>)-> String{
    let units = match unito {
        None => "TwoTheta",
        Some(s) => {match s{
            "twotheta" | "2theta" | "TwoTheta" | "2Theta" => "TwoTheta",
            "QA"|"Qa"|"qa" => "QA",
            "Qnm"|"qnm" => "Qnm",
            _ => {println!("couldn't interpret unit string, defaulting to TwoTheta");"TwoTheta"}

        }}
    };
    String::from(units)
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
    units:String,

}



impl MultiFile{

    pub fn build(cbfdir:&String, ponidir: &String, tthmin:f64, tthmax:f64, tthbins:usize, chimin: f64, chimax: f64, 
                chibins:usize, pfactor:f64, maskfile: Option<&Path>, maskdir: Option<String>, unit:Option<&str>, ymotor:&String,
            zmotor:&String) -> Result<MultiFile,BuildError> {
                    //let units = strtoUnits(units);
                    let units = optiontostr(unit);
                    let mut dc: Option<Arc<DetectorConfig>> = None;
                    let pattern = format!("{cbfdir}/*.cbf");
                    let cbffiles = glob(&pattern).unwrap();
                    let mut cbffile: Arc<PathBuf>;
                    let binding: edf::Edf;
                    let mask: Option<&Array> = match maskfile {
                        None => None,
                        Some(f) => {binding = edf::Edf::open(f)
                                          .expect(format!("couldn't open or find mask file {f:?}").as_str());
                                          Some(binding.array())},
                    };
                    let mut usedmask : Option<&Array> = mask.clone();
                    let mut ilist:Vec<ImagePoni> = Vec::new();
                    let mut mbinding:Edf;
                    let mut getdetconf=true;
                    for fresult in cbffiles{
                        cbffile = Arc::new(fresult.unwrap());
                        let cbfclone = cbffile.clone();
                        
                        let basename = cbfclone.file_name().unwrap().to_str().unwrap().replace(".cbf", "");
                        let ponifiles = glob(&format!("{ponidir}/*.poni")).unwrap();
                        if let Some(md) = &maskdir{
                            let maskfiles = glob(&format!("{md}/*.edf")).unwrap();
                            for mresult in maskfiles{
                                let m = mresult.unwrap();
                                let mut breakloop:bool = false;
                                let mbase = m.file_name().unwrap().to_str().unwrap().replace(".edf", "");
                                usedmask = match yzcompare(&cbffile.clone(), &m, ymotor, zmotor){
                                    false => mask,
                                    true => {mbinding = Edf::open(m).unwrap();
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
                            if yzcompare(&cbffile.clone(), &ponifile, ymotor, zmotor){
                                
                                println!("{:?}, {:?}", cbffile.clone(),ponifile);

                                if getdetconf{
                                    let ptemp = Poni::open(&ponifile, None)
                                    .expect(&format!("couldn't read poni file {:?}",&ponifile));
                                    dc = ptemp.detector_config;
                                    getdetconf = false;
                                }
                                let ip = match ImagePoni::build(&ponifile,&cbffile.clone(), dc.clone(),
                                 usedmask){
                                    Ok(ip) => ip,
                                    Err(e) => return Err(e),
                                };
                                ilist.push(ip);
                                break;
                            }
                        }
                    }
                    if ilist.len() < 2{
                        let nitems = ilist.len();
                        eprintln!("build function requires at least 2 pairs of ponis and cbfs, found {nitems}. cbf and poni files must match in the base name");
                        return Err(BuildError)
                    };
                    //let mut it = cbffiles.into_iter();
                    
                    Ok(MultiFile { ilist,  tthmin, tthmax, tthbins, chimin, chimax, chibins, pfactor, units })
                }

    pub fn buildinterpolate(cbfdir:&String, ponidir: &String,tthmin:f64,tthmax:f64,tthbins:usize,chimin:f64, chimax:f64, 
        chibins:usize,pfactor:f64, maskfile: Option<&Path>, maskdir: Option<String>, ponipattern:&String, ymotor:&String, 
        zmotor:&String,saveponis:bool, unit:Option<&str>)-> Result<MultiFile, BuildError>{
        let units = optiontostr(unit);
        let plist = PoniList::build(ponidir, ponipattern, ymotor, zmotor);
        let p0 = plist.ponilist[0].poni.clone();

        let t = plist.gettriangulations();
        let interp = Interpolators::build(&t);

        let cbfpattern = format!("{cbfdir}/*.cbf");
        let mut ilist: Vec<ImagePoni> = Vec::new();
        let binding: Edf;
        let mask: Option<&Array> = match maskfile {
            None => None,
            Some(f) => {binding = edf::Edf::open(f).expect(&format!("couldn't find mask file {f:?}"));
                                Some(binding.array())},
        };
        let cbffiles = glob(&cbfpattern).unwrap();
        let mut usedmask = mask.clone();
        let mut etmp: Edf;
        let outponidir = format!("{}/savedponis",cbfdir);
        if !Path::new(&outponidir).exists() & saveponis{
            std::fs::create_dir(&outponidir).unwrap();
        }

        for fresult in cbffiles{
            let f: PathBuf = fresult.unwrap();
            let fstring = String::from(f.to_str().unwrap());
            let (yo,zo) = getyz(&f, ymotor, zmotor);
            let y = match yo{
                None =>  {eprintln!("couldn't find y value for file {}", &fstring);return Err(BuildError)},
                Some(val) => val,
            };
            let z = match zo{
                None => {eprintln!("couldn't find z value for file {}", &fstring);return Err(BuildError)},
                Some(val) => val,
            };
            
            let poni = match PoniList::interpolatexy( y, z, p0.clone(), &interp){
                Err(_e) => return Err(BuildError),
                Ok(p) => p,
            };
           
            if let Some(ref md) = maskdir{
                let mfiles = glob(&format!("{md}/*.edf")).unwrap();
                usedmask = mask.clone();
                for mresult in mfiles{
                    let m = mresult.unwrap();
                    if yzcompare(&f, &m, ymotor, zmotor){
                        println!("using mask {:?} for {:?}",&m, &f);
                        etmp = Edf::open(&m).expect(&format!("couldn't open mask file {:?}",&m));
                        usedmask = Some(etmp.array());
                        break;
                    }
                }
            }
            if saveponis{
                let cbfbase = Path::new(&fstring).file_name()
                .unwrap();
                let outponistr = String::from(cbfbase.to_str().unwrap()).replace(".cbf", ".poni");
                let mut file = File::create(format!("{outponidir}/{}", outponistr))
                .expect(&format!("couldn't create file: {}/{}",&outponidir,&outponistr));
                file.write_all( poni.to_string().as_bytes()).unwrap();
            }
            let ip = match ImagePoni::buildfromponi(poni, &f, usedmask){
                Ok(i) => i,
                Err(e) => return Err(e),
            };
            ilist.push(ip);
        }
        if ilist.len() < 2{
            eprintln!("coudn't find any files in {cbfdir}");
            return Err(BuildError)
        }
        Ok(MultiFile{ilist, tthmin,tthmax, tthbins, chimin,chimax,chibins,pfactor, units})
    }
    pub fn integrate_all(self, cakedir: Option<&String>)->Vec<Cake>{
        //let mut cakes :Vec<Cake> = Vec::new(); //vec![Default::default(); self.ilist.len()];
        
        println!("integrating images");

        let cakes: Vec<Cake> = self.ilist.into_par_iter()
        .enumerate()
        .map(|(i,ip)|{
            print!("{i}, ");
            io::stdout().flush().unwrap();
            ip.get_cake(self.tthmin, self.tthmax, self.tthbins, self.chimin, self.chimax, 
                self.chibins, self.pfactor, cakedir, &self.units)
        }).collect();
     
        cakes
    }

    pub fn average_cakes(self, medianfilter:f64, cakedir: Option<&String>, avdir: &String, cakemaskfile: Option<String>)->Cake{
        let cakes = self.integrate_all(cakedir);
        println!("\naveraging cakes");
        let c0 = &cakes[0];
        let rpos = &c0.radial_positions;
        let azpos = &c0.azimuthal_positions;
        let chisize = c0.cake.dim1();
        let radsize = c0.cake.dim2();
        let _ = create_dir(&avdir);
        if let Some(cd) = cakedir{
            let _ = create_dir(&cd);
        }


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
        println!("saving cake to {}",&fnameav);
        newcake.store(&fnameav, None).unwrap();
        
        let av1d_alt = cakeav(&cakes, cakemask, medianfilter);
        let fname1d_alt = format!("{avdir}/av1d.xy");
        save1d(fname1d_alt, rpos, &av1d_alt, None);

        newcake
    }

    pub fn integrate_fluosub(self, medianfilter:f64, cakedir: Option<&String>, avdir: &String, cakemaskfile: Option<String>, 
                        fluo_k0:f64, tthindex:usize)->Cake{
        let pfactor = self.pfactor;
        if let Some(cd) = cakedir{
            let _ = create_dir(&cd);
        };
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

    pub fn calculateflatfield(&self, outdir:&String, ffmin: f64, ffmax: f64){
        let tthmin = self.tthmin;
        let tthmax = self.tthmax;
        let tthbins = self.tthbins;
        let chimin = self.chimin;
        let chimax = self.chimax;
        let tthspacing = (tthmax - tthmin)/(tthbins as f64 - 1.);
        let mut tthrange : Vec<f64> = Vec::new();
        for i in 0..tthbins{
            tthrange.push(tthmin + tthspacing * i as f64);
        };
        let mut tthvalues: Vec<f64> = vec![0.; tthbins];
        let mut tthdiv: Vec<f64> = vec![0.; tthbins];
        let mut tthindexvec : Vec<Vec<i32>> = Vec::new();
        let deg = 180./PI;
        println!("calculating 1d pattern");
        let cbf0 = &self.ilist[0].cbf;
        let dim1 = cbf0.dim1();
        let dim2 = cbf0.dim2();
        let scale = 1e7;
        println!("dim1: {dim1}, dim2: {dim2}");
        for (n,ip) in self.ilist.iter().enumerate(){
            tthindexvec.push(Vec::new());
            let geo = ip.poni.geometry(&Units::TwoTheta, self.pfactor, 0., 0.);
            
            let data = ip.cbf.array().data();    
            
            print!("{n}, ");
            io::stdout().flush().unwrap();
            for (i,pix) in data.iter().enumerate(){
                let y = i / dim2;
                let x =  i % dim2;
                //let (tth, chi) = geo.compute_tth_chi(y as f64, x as f64);
                let pd = geo.compute_pixel(y, x);
                let pol = pd.polar;
                let sa = pd.sa;
                let tth = pd.tth;
                let chi = pd.chi;
                let tthdeg = tth*deg;
                let chideg = chi*deg;
                if (tthdeg < tthmin - tthspacing/2.) | (tthdeg > tthmax + tthspacing/2.) | (chideg < chimin) | (chideg > chimax) | (*pix <= 0.){
                    tthindexvec[n].push(-1);
                    continue;
                }
                let tthindex = closestindexordered(&tthrange, tthdeg);
                
                tthvalues[tthindex] += *pix * scale/(pol * sa);
                tthdiv[tthindex] += 1.;
                tthindexvec[n].push(tthindex as i32);
            }
        }
        let mut tthav: Vec<f64> = Vec::new();
        for (val, div) in tthvalues.iter().zip(tthdiv.iter()){
            tthav.push(val/div);
        }
        let out1d = format!("{outdir}/ff1d.xy");
        save1d(out1d, &tthrange, &tthav, None);
        println!("\npattern calculated");

        let cbfsize = self.ilist[0].cbf.array().data().len();
        let dim1 = self.ilist[0].cbf.dim1();
        let dim2 = self.ilist[0].cbf.dim2();

        let mut gainsum : Vec<f64> = vec![0.;cbfsize];
        let mut gaindiv : Vec<f64> = vec![0.;cbfsize];

        //println!("{tthindexvec:?}");
        for (ip, tthiv) in self.ilist.iter().zip(tthindexvec.iter()){
            let geo = ip.poni.geometry(&Units::TwoTheta, self.pfactor, 0., 0.);
            let data = ip.cbf.array().data();
            for (i, (pix, tthi)) in data.iter().zip(tthiv).enumerate(){
                let y = i / dim2;
                let x = i % dim2;
                //let (tth, chi) = geo.compute_tth_chi(y as f64, x as f64);
                let pd = geo.compute_pixel(y, x);
                let pol = pd.polar;
                let sa = pd.sa;
                if *tthi < 0 {
                    continue;
                }

                let gain = (pix * scale/(pol * sa)) / tthav[*tthi as usize];
                gainsum[i] += gain;
                gaindiv[i] += 1.;
            }
        }

        let mut flatfield : Vec<f64> = Vec::new();
        for (g,d) in gainsum.iter().zip(gaindiv.iter()){
            if *d <= 0. {
                flatfield.push(-1.);
                continue;
            }
            let value = g/d;
            
            if (value < ffmin) | (value > ffmax){
                flatfield.push(-1.);
                continue;
            }
            
            flatfield.push(value);
             
        }
        let mut header = Header::new();
        header.insert(Cow::Borrowed(KEY_BUBBLE_MADE), HeaderEntry::Number(1)); //tells bubble not to integrate
        let flatfieldarray = Array::with_data(dim1, dim2, flatfield);
        let dt : DateTime<Local> = Local::now();
        let dtstring = format!("{:04}{:02}{:02}",dt.year(), dt.month(), dt.day());
        let fname = format!("{}/{}_p{:.2}_flatfield.edf", outdir, dtstring, self.pfactor);
        let mut writer = BufWriter::new(File::create(&fname).unwrap());
        Edf::save_array(&flatfieldarray, &mut header, &mut writer, edf::DataType::F32).unwrap();
        println!("flat-field image saved to {}",&fname);
    
    }
}



#[cfg(test)]
mod tests{
    use super::*;

    #[test]
    fn mfbuildtest(){
        let cbfdir = String::from("D:\\beamlineData\\April2026\\multipos3_tilt\\C60");
        let ponidir = String::from("D:\\beamlineData\\April2026\\multipos3_tilt\\C60\\ponis_rp");
        let tthmin = 0.75;
        let tthmax = 68.;
        let tthbins:usize = 5000;
        let chimin = 220.;
        let chimax = 320.;
        let chibins:usize = 101;
        let pfactor = 0.85;
        let maskfile = None;
        let maskdir = None;
        let unit = Some("TwoTheta");
        let ymotor = String::from("dty");
        let zmotor = String::from("dtz");
        let _mf = MultiFile::build(&cbfdir, &ponidir, tthmin, tthmax, tthbins, chimin, chimax, chibins, 
            pfactor, maskfile, maskdir, unit, &ymotor, &zmotor).unwrap();

    }

}