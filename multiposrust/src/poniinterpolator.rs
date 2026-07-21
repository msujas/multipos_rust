use std::{ sync::Arc};
use crate::functions::getyz;
use cryiorust::poni::{DetectorConfig, Poni};
use glob::glob;
use spade::{DelaunayTriangulation, HasPosition, NaturalNeighbor, Point2, Triangulation};


#[derive(Debug)]
pub struct InterpolationError;

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
#[derive(Clone, Debug)]
pub struct PoniYZ{
    pub y: f64,
    pub z: f64,
    pub poni: Poni,
}

impl PoniYZ{
    pub fn fromfile(ponifile:&String, y:f64, z:f64, dc:Option<Arc<DetectorConfig>>)-> PoniYZ{
        let poni = Poni::open(ponifile, dc).unwrap();
        PoniYZ { y, z, poni }
    }

    pub fn fromparams(y:f64,z:f64, poni1:f64,poni2:f64, distance:f64,rot1:f64,rot2:f64, rot3:f64, wavelength:f64, dc:Option<Arc<DetectorConfig>>, version:f32, pixel1: f64, pixel2:f64 )-> PoniYZ{
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

    pub fn new()-> PoniYZ{
        
        PoniYZ { y:0.,z: 0.,poni: Poni::new() }
    }
}

#[derive(Clone, Debug)]
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
            let (yo,zo) = getyz(fname, ymotor, zmotor);
            let mut pyz = PoniYZ::new();
            match yo {
                None => {panic!("couldn't find y value for {:?}", fname);}
                Some(y) => {pyz.y = y;}
            }
            match zo {
                None => {panic!("couldn't find z value for {:?}",fname);}
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
        if plist.len() == 0{
            panic!("didn't find any poni files at {ponidir} with pattern {ponipattern}");
        }
        PoniList { ponilist:plist}
    }

    pub fn gettriangulations(self)->Triangluators{
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
        let nnponi1: NaturalNeighbor<'_, DelaunayTriangulation<PointHeight>> = tponi1.natural_neighbor();
        let nnponi2: NaturalNeighbor<'_, DelaunayTriangulation<PointHeight>> = tponi2.natural_neighbor();
        let nnrot1: NaturalNeighbor<'_, DelaunayTriangulation<PointHeight>> = trot1.natural_neighbor();
        let nnrot2: NaturalNeighbor<'_, DelaunayTriangulation<PointHeight>> = trot2.natural_neighbor();
        let nnrot3: NaturalNeighbor<'_, DelaunayTriangulation<PointHeight>> = trot3.natural_neighbor();
        let nndist: NaturalNeighbor<'_, DelaunayTriangulation<PointHeight>> = tdist.natural_neighbor();
        */
        //[nnponi1,nnponi2, nndist, nnrot1,nnrot2, nnrot3]
        //(tponi1, tponi2,tdist,trot1,trot2,trot3)
        Triangluators { tponi1, tponi2, tdist, trot1, trot2, trot3 }
        
    }

    
    pub fn interpolatexy( y:f64, z:f64, poni0:Poni, interp: &Interpolators)->Result<Poni, InterpolationError>{
        let interpolationerrormessage = "couldn't interpolate poni value, likely out of range. Maybe ponis for some corner positions are missing";
        let poni1 = match interp.nnponi1.interpolate(|v| v.data().height, Point2::new(y,z)){
            Some(p) => p,
            None => {eprintln!("\ny: {y}, z: {z},\n{interpolationerrormessage}");
            return Err(InterpolationError)}
        };
        let poni2 = interp.nnponi2.interpolate(|v| v.data().height, Point2::new(y,z))
        .expect(interpolationerrormessage);
        let dist = interp.nndist.interpolate(|v| v.data().height, Point2::new(y,z))
        .expect(interpolationerrormessage);
        let rot1 = interp.nnrot1.interpolate(|v| v.data().height, Point2::new(y,z))
        .expect(interpolationerrormessage);
        let rot2 = interp.nnrot2.interpolate(|v| v.data().height, Point2::new(y,z))
        .expect(interpolationerrormessage);
        let rot3 = interp.nnrot3.interpolate(|v| v.data().height, Point2::new(y,z))
        .expect(interpolationerrormessage);

        let pversion = poni0.version;
        let pix1 = poni0.pixel1;
        let pix2 = poni0.pixel2;
        let wavelength = poni0.wavelength;
        let dc = poni0.detector_config;
        let mut poni = Poni::new();
        poni.poni1 = poni1;
        poni.poni2 = poni2;
        poni.distance = dist;
        poni.rot1 = rot1;
        poni.rot2= rot2;
        poni.rot3 = rot3;
        poni.detector_config = dc.clone();
        poni.version = pversion;
        poni.wavelength = wavelength;
        poni.pixel1 = pix1;
        poni.pixel2 = pix2;

        Ok(poni)

     }

    pub fn interpolateponi(self,y:f64,z:f64)-> Poni{
        let p0 = self.ponilist[0].poni.clone();
        let dc = p0.detector_config.clone();
        let version = p0.version;
        let pixel1 = p0.pixel1;
        let pixel2 = p0.pixel2;
        let wavelength = self.ponilist[0].poni.wavelength;
        let t = self.gettriangulations();
        let interp = Interpolators::build(&t);
        let poni1 = interp.nnponi1.interpolate(|v| v.data().height, Point2 ::new(y,z)).unwrap();
        let poni2 = interp.nnponi2.interpolate(|v| v.data().height, Point2 ::new(y,z)).unwrap();
        let dist = interp.nndist.interpolate(|v| v.data().height, Point2 ::new(y,z)).unwrap();
        let rot1 = interp.nnrot1.interpolate(|v| v.data().height, Point2 ::new(y,z)).unwrap();
        let rot2 = interp.nnrot2.interpolate(|v| v.data().height, Point2 ::new(y,z)).unwrap();
        let rot3 = interp.nnrot3.interpolate(|v| v.data().height, Point2 ::new(y,z)).unwrap();
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

pub struct Triangluators{
    tponi1: DelaunayTriangulation::<PointHeight>,
    tponi2: DelaunayTriangulation::<PointHeight>,
    tdist: DelaunayTriangulation::<PointHeight>,
    trot1: DelaunayTriangulation::<PointHeight>,
    trot2: DelaunayTriangulation::<PointHeight>,
    trot3: DelaunayTriangulation::<PointHeight>,
}


impl Triangluators{
    pub fn build(ponilist: PoniList)-> Triangluators{
        let mut tponi1 =  DelaunayTriangulation::<PointHeight>::new();
        let mut tponi2 = DelaunayTriangulation::<PointHeight>::new();
        let mut trot1 = DelaunayTriangulation::<PointHeight>::new();
        let mut trot2 = DelaunayTriangulation::<PointHeight>::new();
        let mut trot3 = DelaunayTriangulation::<PointHeight>::new();
        let mut tdist = DelaunayTriangulation::<PointHeight>::new();
        for pyz in ponilist.ponilist{
            tponi1.insert(PointHeight{position:Point2::new( pyz.y, pyz.z), height: pyz.poni.poni1}).unwrap();
            tponi2.insert(PointHeight{position:Point2::new( pyz.y, pyz.z), height: pyz.poni.poni2}).unwrap();
            trot1.insert(PointHeight{position:Point2::new( pyz.y, pyz.z), height: pyz.poni.rot1}).unwrap();
            trot2.insert(PointHeight{position:Point2::new( pyz.y, pyz.z), height: pyz.poni.rot2}).unwrap();
            tdist.insert(PointHeight{position:Point2::new( pyz.y, pyz.z), height: pyz.poni.distance}).unwrap();
            trot3.insert(PointHeight{position:Point2::new( pyz.y, pyz.z), height: pyz.poni.rot3}).unwrap();
        };
        Triangluators { tponi1, tponi2, tdist, trot1, trot2, trot3 }
    }
}

pub struct Interpolators<'a>{
    nnponi1: NaturalNeighbor<'a, DelaunayTriangulation<PointHeight>> ,
    nnponi2: NaturalNeighbor<'a, DelaunayTriangulation<PointHeight>> ,
    nndist: NaturalNeighbor<'a, DelaunayTriangulation<PointHeight>>,
    nnrot1: NaturalNeighbor<'a, DelaunayTriangulation<PointHeight>> ,
    nnrot2: NaturalNeighbor<'a, DelaunayTriangulation<PointHeight>> ,
    nnrot3: NaturalNeighbor<'a, DelaunayTriangulation<PointHeight>> ,
}

impl Interpolators<'_>{
    pub fn build<'a>(t: &'a Triangluators)->Interpolators<'a>{
        let nnponi1 = t.tponi1.natural_neighbor();
        let nnponi2 = t.tponi2.natural_neighbor();
        let nndist = t.tdist.natural_neighbor();
        let nnrot1 = t.trot1.natural_neighbor();
        let nnrot2 = t.trot2.natural_neighbor();
        let nnrot3 = t.trot3.natural_neighbor();
        Interpolators { nnponi1, nnponi2, nndist, nnrot1, nnrot2, nnrot3 }
    }
}

