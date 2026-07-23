use clap::{ Parser};

#[derive(Parser,Debug)]
#[command(version,  about = concat!("program for processing X-ray scattering from multiple detector positions.\n",
            "Images and poni files must contain detector positions in name separated by _.\n",
            "e.g. <name>_<ymotor><ypos>_<zmotor><zpos>.cbf or <name>_<ymotor><ypos>_<zmotor><zpos>_<image number>.cbf"), 
long_about=None)]
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
    #[arg(long, default_value_t= 358.)]
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
    /// poni directory - only need corner (convex hull) positions and the rest will be interpolated
    #[arg(long)]
    pub ponidir: String,
    /// save individual cakes or not
    #[arg(long)]
    pub savecakes: bool,
    /// subdirectory to store cake file
    #[arg(long, default_value="cakes")]
    pub cakesubdir: String,
    /// mask file path (optional)
    #[arg(short, long, default_value=None)]
    pub maskfile: Option<String>,
    /// cake mask path (optional)
    #[arg(short='k', long, default_value=None)]
    pub cakemaskfile: Option<String>,
    /// directory for individual masks (optional, matched with y and z positions)
    #[arg(long, default_value=None)]
    pub maskdir: Option<String>,
    /// string pattern used to find poni files in directory (must include asterix)
    #[arg(long, default_value="*.poni")]
    pub ponipattern: String,
    #[arg(long, default_value="dty", help = concat!("ymotor name used to find detector y position in file name\n",
    "(format ..._<ymotor>yyy.yy_<zmotor>zzz.zz_...)"))]
    pub ymotor: String,
    #[arg(long, default_value="dtz", help = concat!("z motor name used to find detector z position in file name\n",
    "(format ..._<ymotor>yyy.yy_<zmotor>zzz.zz_...)"))]
    pub zmotor: String,
    /// do fluo subtraction or not
    #[arg(short, long, default_value_t=false)]
    pub fluosub: bool,
    /// fluok starting value
    #[arg(long, default_value_t=1.)]
    pub fluok0: f64,
    /// save individual ponis
    #[arg(long, default_value_t=false)]
    pub saveponis: bool,
    #[arg(short, long, default_value="TwoTheta", 
    help = "integration unit. Options TwoTheta/2Theta/2theta/twotheta, QA/qa, Qnm/qnm.\nWill default to TwoTheta if invalid")]
    pub unit: String,

}