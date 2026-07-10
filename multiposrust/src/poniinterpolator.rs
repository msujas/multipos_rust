use std::{ sync::Arc};

use cryiorust::poni::{DetectorConfig, Poni};
use glob::glob;
use spade::{DelaunayTriangulation, HasPosition,  Point2, Triangulation};

pub fn getyz(fname:&String, ymotor:&String, zmotor:&String)->(Option<f64>,Option<f64>){

        let fsplit = fname.split("_");
        let mut yo: Option<f64> = None;
        let mut zo: Option<f64> = None;
        for item in fsplit{
            if item.contains(ymotor){
                let y = item.replace(ymotor, "").parse::<f64>().unwrap();
                yo = Some(y);
            }
            if item.contains(zmotor){
                let z = item.replace(zmotor, "").parse::<f64>().unwrap();
                zo = Some(z);
            }
        }
        (yo,zo)   
}

pub struct PointHeight{
    pub position: Point2<f64>,
    pub height: f64,
}


impl HasPosition for PointHeight {
    type Scalar = f64;

    fn position(&self) -> Point2<f64> {
        self.position
    }
}
#[derive(Clone)]
pub struct PoniYZ{
    pub y: f64,
    pub z: f64,
    pub poni: Poni,
}

impl PoniYZ{
    fn fromfile(ponifile:&String, y:f64, z:f64, dc:Option<Arc<DetectorConfig>>)-> PoniYZ{
        let poni = Poni::open(ponifile, dc).unwrap();
        PoniYZ { y, z, poni }
    }

    fn fromparams(y:f64,z:f64, poni1:f64,poni2:f64, distance:f64,rot1:f64,rot2:f64, rot3:f64, wavelength:f64, dc:Option<Arc<DetectorConfig>>, version:f32, pixel1: f64, pixel2:f64 )-> PoniYZ{
        let mut poni = Poni::new();
        poni.poni1 = poni1;
        poni.poni2 = poni2;
        poni.distance = distance;
        poni.rot1 = rot1;
        poni.rot2 = rot2;
        poni.rot3 = rot3;
        poni.wavelength = wavelength;
        poni.detector_config = dc;
        poni.version = version;
        poni.pixel1 = pixel1;
        poni.pixel2 = pixel2;

        PoniYZ { y, z, poni }
    }

    fn new()-> PoniYZ{
        
        PoniYZ { y:0.,z: 0.,poni: Poni::new() }
    }
}

#[derive(Clone)]
pub struct PoniList{
    pub ponilist: Vec<PoniYZ>
}

impl PoniList{
    pub fn build(ponidir:&String, ponipattern: &String, ymotor: &String, zmotor: &String)-> PoniList{
        let mut plist: Vec<PoniYZ> = Vec::new();
        let mut dc: Option<Arc<DetectorConfig>> = None;
        let searchpattern = format!("{ponidir}/{ponipattern}");
        let ponifiles = glob(&searchpattern).unwrap();
        for fresult in ponifiles{
            let fname = fresult.as_ref().unwrap();
            println!("{fname:?}");
            let f = String::from(fresult.as_ref().unwrap().file_name().unwrap().to_str().unwrap());
            let (yo,zo) = getyz(&f, ymotor, zmotor);
            let mut pyz = PoniYZ::new();
            match yo {
                None => {panic!("couldn't find y value for {f}");}
                Some(y) => {pyz.y = y;}
            }
            match zo {
                None => {panic!("couldn't find z value for {f}");}
                Some(z) => {pyz.z = z;}
            }
            let poni = match dc{
                None => {let p = Poni::open(fname, None).unwrap(); dc = p.detector_config.clone();
                p},
                Some(ref dc) => Poni::open(fname, Some(dc.clone())).unwrap()
            };
            pyz.poni = poni;
            plist.push(pyz);
        }
        PoniList { ponilist:plist}
    }

    pub fn getinterpolators(self)->(DelaunayTriangulation::<PointHeight>, DelaunayTriangulation::<PointHeight>,DelaunayTriangulation::<PointHeight>,
    DelaunayTriangulation::<PointHeight>,DelaunayTriangulation::<PointHeight>,DelaunayTriangulation::<PointHeight>,){
        let mut tponi1 =  DelaunayTriangulation::<PointHeight>::new();
        let mut tponi2 = DelaunayTriangulation::<PointHeight>::new();
        let mut trot1 = DelaunayTriangulation::<PointHeight>::new();
        let mut trot2 = DelaunayTriangulation::<PointHeight>::new();
        let mut trot3 = DelaunayTriangulation::<PointHeight>::new();
        let mut tdist = DelaunayTriangulation::<PointHeight>::new();
        for pyz in self.ponilist{
            tponi1.insert(PointHeight{position:Point2::new( pyz.y, pyz.z), height: pyz.poni.poni1}).unwrap();
            tponi2.insert(PointHeight{position:Point2::new( pyz.y, pyz.z), height: pyz.poni.poni2}).unwrap();
            trot1.insert(PointHeight{position:Point2::new( pyz.y, pyz.z), height: pyz.poni.rot1}).unwrap();
            trot2.insert(PointHeight{position:Point2::new( pyz.y, pyz.z), height: pyz.poni.rot2}).unwrap();
            tdist.insert(PointHeight{position:Point2::new( pyz.y, pyz.z), height: pyz.poni.distance}).unwrap();
            trot3.insert(PointHeight{position:Point2::new( pyz.y, pyz.z), height: pyz.poni.rot3}).unwrap();
        }
        /*
        let nnponi1: NaturalNeighbor<'static, DelaunayTriangulation<PointHeight>> = tponi1.natural_neighbor();
        let nnponi2: NaturalNeighbor<'static, DelaunayTriangulation<PointHeight>> = tponi2.natural_neighbor();
        let nnrot1: NaturalNeighbor<'static, DelaunayTriangulation<PointHeight>> = trot1.natural_neighbor();
        let nnrot2: NaturalNeighbor<'static, DelaunayTriangulation<PointHeight>> = trot2.natural_neighbor();
        let nnrot3: NaturalNeighbor<'static, DelaunayTriangulation<PointHeight>> = trot3.natural_neighbor();
        let nndist: NaturalNeighbor<'static, DelaunayTriangulation<PointHeight>> = tdist.natural_neighbor();
        */
        (tponi1, tponi2,tdist,trot1,trot2,trot3)
        
    }
    fn interpolateponi(self,y:f64,z:f64)-> Poni{
        let p0 = self.ponilist[0].poni.clone();
        let dc = p0.detector_config.clone();
        let version = p0.version;
        let pixel1 = p0.pixel1;
        let pixel2 = p0.pixel2;
        let wavelength = self.ponilist[0].poni.wavelength;
        let (tponi1, tponi2,tdist,trot1,trot2, trot3) = self.getinterpolators();
        let nnponi1  = tponi1.natural_neighbor();
        let nnponi2 = tponi2.natural_neighbor();
        let nnrot1  = trot1.natural_neighbor();
        let nnrot2  = trot2.natural_neighbor();
        let nnrot3  = trot3.natural_neighbor();
        let nndist  = tdist.natural_neighbor();
        let poni1 = nnponi1.interpolate(|v| v.data().height, Point2 ::new(y,z)).unwrap();
        let poni2 = nnponi2.interpolate(|v| v.data().height, Point2 ::new(y,z)).unwrap();
        let dist = nndist.interpolate(|v| v.data().height, Point2 ::new(y,z)).unwrap();
        let rot1 = nnrot1.interpolate(|v| v.data().height, Point2 ::new(y,z)).unwrap();
        let rot2 = nnrot2.interpolate(|v| v.data().height, Point2 ::new(y,z)).unwrap();
        let rot3 = nnrot3.interpolate(|v| v.data().height, Point2 ::new(y,z)).unwrap();
        let mut poni = Poni::new();
        poni.poni1=poni1;
        poni.poni2 = poni2;
        poni.distance = dist;
        poni.rot1 = rot1;
        poni.rot2 = rot2;
        poni.rot3 = rot3;
        poni.detector_config = dc;
        poni.wavelength = wavelength;
        poni.version = version;
        poni.pixel1 = pixel1;
        poni.pixel2 = pixel2;
        poni
    }
}


#[cfg(test)]
mod test {
    use glob::glob;

use crate::poniinterpolator::getyz;

    #[test]
    fn getyztest(){
        let fname= String::from("emptyCap_dty138.79_dtz108.00_001_0001p.poni");
        let ymotor = String::from("dty");
        let zmotor = String::from("dtz");
        let (yo,zo) = getyz(&fname, &ymotor, &zmotor);
        let y = yo.unwrap();
        let z = zo.unwrap();
        assert_eq!(y, 138.79);
        assert_eq!(z, 108.00);
    }

}