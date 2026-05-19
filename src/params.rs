use std::path::Path;
use clap::{Parser};

/// program for processing X-ray scattering from multiple detector positions
#[derive(Parser,Debug)]
#[command(version, about, long_about = None)]
pub(crate) struct Params{
    /// minimum 2theta
    #[arg(short, long)]
    pub tthmin : f64,
    /// maximum 2theta
    #[arg(short='x', long)]
    pub tthmax : f64,
    /// number of 2theta bins
    #[arg(short='b', long, default_value_t = 5000)]
    pub tthbins : usize,
    /// minimum chi
    #[arg(short='c', long, default_value_t = 2.)]
    pub chimin : f64,
    /// maximum chi
    #[arg(short='a', long, default_value_t= 358.)]
    pub chimax : f64,
    /// number of chi bins
    #[arg(short='i', long, default_value_t= 356)]
    pub chibins : usize,
    /// polarization factor
    #[arg(short, long, default_value_t = 0.85)]
    pub pfactor: f64,
    /// cbf directory
    #[arg(short='d', long, default_value  = ".")]
    pub cbfdir: String,
    /// poni directory
    #[arg(short='o', long)]
    pub ponidir: String,  

}