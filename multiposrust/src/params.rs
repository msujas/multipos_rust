use clap::{Parser};

/// program for processing X-ray scattering from multiple detector positions
#[derive(Parser,Debug)]
#[command(version, about, long_about = None)]
pub struct Params{
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
    #[arg(short='i', long, default_value_t= 357)]
    pub chibins : usize,
    /// polarization factor
    #[arg(short, long, default_value_t = 0.85)]
    pub pfactor: f64,
    /// cbf directory
    #[arg(short='d', long, default_value  = ".")]
    pub cbfdir: String,
    /// poni directory - only need most extreme positions and the rest will be interpolated
    #[arg(short='o', long)]
    pub ponidir: String,
    /// save individual cakes or not
    #[arg(short, long)]
    pub savecakes: bool,
    /// subdirectory to store cake file
    #[arg(short='u', long, default_value="cakes")]
    pub cakesubdir: String,
    /// mask file path (optional)
    #[arg(short, long, default_value=None)]
    pub maskfile: Option<String>,
    /// cake mask path (optional)
    #[arg(short='k', long, default_value=None)]
    pub cakemaskfile: Option<String>,
    /// directory for individual masks (matched with y and z positions)
    #[arg(long, default_value=None)]
    pub maskdir: Option<String>,
    /// string pattern used to find poni files in directory (must include asterix)
    #[arg(long, default_value="*.poni")]
    pub ponipattern: String,
    /// ymotor name used to find detector y position in file name (format ..._<ymotor>yyy.yy_<zmotor>zzz.zz_...)
    #[arg(long, default_value="dty")]
    pub ymotor: String,
    /// z motor name used to find detector z position in file name (format ..._<ymotor>yyy.yy_<zmotor>zzz.zz_...)
    #[arg(long, default_value="dtz")]
    pub zmotor: String,
    /// do fluo subtraction or not
    #[arg(short, long, default_value=None)]
    pub fluosub: bool,
    /// fluok starting value
    #[arg(long, default_value_t=1.)]
    pub fluok0: f64,
    /// save individual ponis
    #[arg(short, long, default_value=None)]
    pub saveponis: bool,
    /// integration unit. Options TwoTheta/2Theta/2theta/twotheta, QA/qa, Qnm/qnm. 
    /// Will default to TwoTheta if invalid
    #[arg(short, long, default_value="TwoTheta")]
    pub unit: String,

}

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
    #[arg(short='a', long, default_value_t= 358.)]
    pub chimax : f64,
    /// polarization factor
    #[arg(short, long, default_value_t = 0.85)]
    pub pfactor: f64,
    /// cbf directory
    #[arg(short='d', long, default_value  = ".")]
    pub cbfdir: String,
    /// poni directory - only need extreme detector positions and the rest will be interpolated
    #[arg(short='o', long)]
    pub ponidir: String,
    /// mask file path (optional)
    #[arg(short, long, default_value=None)]
    pub maskfile: Option<String>,
    /// directory for individual masks (matched with y and z positions)
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
    #[arg(short, long, default_value=None)]
    pub saveponis: bool,
    /// integration units, default 2theta. Options TwoTheta/2Theta/2theta/twotheta, QA/qa, Qnm/qnm. 
    /// Will default to TwoTheta if invalid
    #[arg(short, long, default_value="TwoTheta")]
    pub unit: String,
}