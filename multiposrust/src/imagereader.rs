use std::path::Path;
use cryiorust::{cbf::Cbf, edf::Edf, eiger::Eiger, frame::{Array, Frame, HeaderEntry::Float}};


#[derive(Debug)]
enum ImageFormat {
    Cbf,
    Edf,
    Eiger,    
}


impl ImageFormat{
    fn fromfilename(filename:&Path)->Result<ImageFormat,ImageReadError >{
        let ext = filename.extension().unwrap().to_str().unwrap();
        match ext{
            "cbf" => return Ok(ImageFormat::Cbf),
            "edf" => return Ok(ImageFormat::Edf),
            "hdf5" => return Ok(ImageFormat::Eiger),
            _ => return Err(ImageReadError)
        }
    }
}


#[derive(Debug)]
pub struct ImageReadError;


pub struct ImageFlux{
    pub namestem: String,
    pub array: Array,
    pub flux: f64,
}

impl ImageFlux{
    pub fn readimage(imagefile: &Path)-> Result<ImageFlux, ImageReadError>{
        let extension = ImageFormat::fromfilename(imagefile)?;
        let c: Cbf;
        let e: Edf;
        let ei: Eiger;
        let dim1: usize;
        let dim2: usize;
        let fluxheaderstring = "# Flux ";
        let name = String::from(imagefile.file_stem().unwrap().to_str().unwrap());
        let fluxo: Option<&cryiorust::frame::HeaderEntry>;
        let im = match extension{
            ImageFormat::Cbf => {c = match Cbf::open(imagefile){Err(_e) => {eprintln!("couldn't open image: {imagefile:?} as cbf"); 
                                                                                            return Err(ImageReadError)},
                                                                        Ok(cb) => cb,};

                                dim1 = c.dim1();
                                dim2 = c.dim2();
                                fluxo = c.header().get(fluxheaderstring);
                                c.array().data().clone()},
            ImageFormat::Edf => {e = match Edf::open(imagefile) { Err(_e) => {eprintln!("couldn't open image {imagefile:?} as edf");
                                                                                                return Err(ImageReadError)},
                                                                    Ok(ed)=> ed};
                                dim1 = e.dim1();
                                dim2 = e.dim2();
                                fluxo = e.header().get(fluxheaderstring);
                                e.array().data().clone()},
            ImageFormat::Eiger => {ei = match Eiger::open(imagefile){Err(_e) =>{ eprintln!("couldn't open {imagefile:?} as hdf5");
                                                                                                return Err(ImageReadError)},
                                                                            Ok(ei) => ei,};
                                    dim1 = ei.dim1();
                                    dim2 = ei.dim2();
                                    fluxo = ei.header().get(fluxheaderstring);
                                    ei.array().data().clone()},
        };
        let flux = match fluxo {
            Some(Float(f)) => f.clone(),
            _ => return Err(ImageReadError),
            
        };

        let array = Array::with_data(dim1, dim2, im);
        Ok(ImageFlux{namestem: name, array, flux})
    }
}