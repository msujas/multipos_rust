use integrustio::integrator::PatternType;
use multipos_rust::{ImagePoni};
use std::path::Path;

fn main() {
    let cbfpath = Path::new("D:/beamlineData/April2026/multipositions2/Si_dtx0/Si_dty261.34_dtz121.50_001_0001p.cbf");
    let ponipath = Path::new("D:/beamlineData/April2026/multipositions2/Si_dty261.34_dtz121.50_MD.poni");
    let ip = ImagePoni::build(ponipath, cbfpath);
    let d = ip.integrate(0.75, 58.,5000, 2.,358., 356,0.85);
    let p: integrustio::integrator::PatternType = d.data;

    p.store("D:/beamlineData/April2026/multipositions2/Si_dtx0/cake.edf", None, None).unwrap();
    let cake = match p{
        PatternType::Cake(c) => c,
        _=>panic!("some issue with getting cake image"),
    };
    
}
