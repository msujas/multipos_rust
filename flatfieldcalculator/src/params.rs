use clap::{Parser};

/// program for calculating flat field given cbfs from multiple detector positions
#[derive(Parser,Debug)]
#[command(version, about, long_about = None)]
pub struct ParamsFF{
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
    #[arg(long, default_value_t= 358.)]
    pub chimax : f64,
    /// polarization factor
    #[arg(short, long, default_value_t = 0.85)]
    pub pfactor: f64,
    /// cbf directory
    #[arg(short='d', long, default_value  = ".")]
    pub cbfdir: String,
    /// poni directory - only need extreme detector positions and the rest will be interpolated
    #[arg(long)]
    pub ponidir: String,
    /// mask file path (optional)
    #[arg(short, long, default_value=None)]
    pub maskfile: Option<String>,
    /// directory for individual masks (optional, matched with y and z positions)
    #[arg(long, default_value=None)]
    pub maskdir: Option<String>,
    /// minimum allowed flat field value
    #[arg(long, default_value_t=0.7)]
    pub ffmin: f64,
    /// maximum allowed flat field value
    #[arg(long, default_value_t=1.5)]
    pub ffmax: f64,
    /// string pattern used to find poni files in directory (must include asterix)
    #[arg(long, default_value="*.poni")]
    pub ponipattern: String,
    /// ymotor name used to find detector y position in file name (format ..._<ymotor>yyy.yy_<zmotor>zzz.zz_...)
    #[arg(long, default_value="dty")]
    pub ymotor: String,
    /// z motor name used to find detector z position in file name (format ..._<ymotor>yyy.yy_<zmotor>zzz.zz_...)
    #[arg(long, default_value="dtz")]
    pub zmotor: String,
    /// save individual ponis
    #[arg(long, default_value_t=false)]
    pub saveponis: bool,
    /// integration units, default 2theta. Options TwoTheta/2Theta/2theta/twotheta, QA/qa, Qnm/qnm. 
    /// Will default to TwoTheta if invalid
    #[arg(short, long, default_value="TwoTheta")]
    pub unit: String,
}