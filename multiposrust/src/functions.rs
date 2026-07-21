use std::{cmp::Ordering, fs::File, io::Write, path::Path};
use cryiorust::{frame::Array};
use integrustio::integrator::Cake;

pub fn getyz(fname:&Path, ymotor:&String, zmotor:&String)->(Option<f64>,Option<f64>){
        let f = fname.file_stem().unwrap().to_str().unwrap();
        let fsplit = f.split("_");
        let mut yo: Option<f64> = None;
        let mut zo: Option<f64> = None;
        for item in fsplit{
            if item.contains(ymotor){
                let ystring = item.replace(ymotor, "");
                let y = ystring.parse::<f64>()
                .expect(&format!("couldn't convert y-string: {ystring} to float"));
                yo = Some(y);
            }
            if item.contains(zmotor){
                let zstring = item.replace(zmotor, "");
                let z = zstring.parse::<f64>()
                .expect(&format!("couldn't convert z-string: {zstring} to float"));
                zo = Some(z);
            }
        }
        (yo,zo)   
}

pub fn yzcompare(file1:&Path, file2: &Path, ymotor:&String, zmotor:&String)->bool{
    let tolerance = 0.00001; //using a tolerance for float comparison
    let (y1,z1) = getyz(file1, ymotor, zmotor);
    let (y2,z2) = getyz(file2, ymotor, zmotor);
    let y1 = match y1{
        None => return false,
        Some(y) => y,
    };
    let y2 = match y2{
        None => return false,
        Some(y) => y,
    };
    let z1 = match z1{
        None => return false,
        Some(z) => z,
    };
    let z2 = match z2{
        None => return false,
        Some(z) => z,
    };
    ((y1 - y2).abs() < tolerance) & ((z1-z2).abs() < tolerance)
}

pub fn cakeav(cakelist: &Vec<Cake>, cakemask: Option<&Array>, medianfilter:f64)-> Vec<f64>{

    let c0 = &cakelist[0];
    let chisize = c0.cake.dim1();
    let tthsize= c0.cake.dim2();
    let mut av1d : Vec<f64> = vec![0.; tthsize];
    let mut divvec: Vec<f64> = vec![0.;tthsize];
    let mut index: usize;

    for c in cakelist{
        for i in 0..tthsize{
            let mut vtemp : Vec<f64> = Vec::new();
            for j in 0..chisize{
                index = i + j*tthsize;
                let value = c.cake.data()[index];
                if c.cake.data()[index] > 0.{
                    if let Some(cakemask)=cakemask{
                        if cakemask.data()[index] > 0.1{
                            continue;
                        }
                    }
                    vtemp.push(value);
                }
            }
            let med = getmedian(&vtemp);
            let mut intensity = 0.;
            let mut div = 0.;
            for item in vtemp{
                if (item < med*medianfilter) & (item > med/medianfilter){
                    intensity += item;
                    div += 1.;
                }
            }
            if div > 0. {
                av1d[i] += intensity/div;
                divvec[i] += 1.;
            }    
        }
    }
    for (x, d) in av1d.iter_mut().zip(divvec.iter_mut()){
        if *d > 0. {
        *x = *x/ *d;
        }
    }
av1d
}

pub fn save1d(fname:String, tthrange: &Vec<f64>, vec1d: &Vec<f64>, sigma : Option<&Vec<f64>>){
    let mut outstring = String::new();
    //for (x,y ) in  tthrange.iter().zip(vec1d.iter()){
    let mut x:f64;
    let mut y:f64;
    let mut e:f64;
    let dosig:bool = match sigma  {
        None => false,
        Some(_s) => true
    };
    for i in 0..tthrange.len(){
        x = tthrange[i];
        y=vec1d[i];
        outstring = outstring + &String::from(format!("{x:.6} {y:.6}"));
        if dosig{
            e = sigma.unwrap()[i];
            outstring = outstring + &String::from(format!(" {e:.6}"));
            }
        outstring = outstring + &String::from("\n");
        }
    println!("saving 1d pattern to {}", &fname);
    let mut file = File::create(&fname).expect(&format!("error creating file {:?}",&fname));
    file.write(outstring.as_bytes()).unwrap();    
}

pub fn getmedian(medvec:&Vec<f64>)->f64{
    //let svec = sortvec(medvec);
    let mut svec = medvec.clone();
    svec.sort_by(cmpf64);
    let vlen = svec.len();
    let pos : usize = vlen/2;
    if vlen == 0 {
        return 0.;
    }
    if vlen%2 == 0 {
        return (svec[pos-1] + svec[pos])/2.;
    }
    return svec[pos]
}

pub fn cmpf64(a:&f64,b:&f64)->Ordering{
    if a > b    {
        return Ordering::Greater
    }
    else if  a < b{
        return  Ordering::Less;
    }
    return Ordering::Equal;
    }
    
/*
pub fn closestindex(avec: &Vec<f64>, value: f64)->usize{
    let mut minindex: usize = 0;
    let mut mindiffsq: f64 = (value - avec[0]).powi(2);
    for (i,v) in avec.iter().enumerate(){
        let diffsq = (value-v).powi(2);
        if diffsq < mindiffsq{
            mindiffsq = diffsq;
            minindex = i;
        }
    }
    minindex
}
 */

/// a function for finding the index of the closest value in an ordered vector. 
/// The vector must be ordered already or it will give the wrong value
pub fn closestindexordered(orderedvec: &Vec<f64>, value: f64)->usize{
    let size = orderedvec.len();
    let mut i0:usize = 0;
    let mut i = size/2;
    let mut iend = size-1;

    let mut midvalue = orderedvec[i];
    loop {
        if value < midvalue{
            iend = i;
            i = (i + i0)/2;
            
        }
        else {
            i0 = i;
            i = (iend + i )/2;
        }
        midvalue = orderedvec[i];
        if (iend - i0) <= 1{
            let de = (value - orderedvec[iend]).powi(2);
            let d0 = (value - orderedvec[i0]).powi(2);
            if d0 < de{
                i = i0;
            }
            else {
                i = iend;
            }
            break;
        }
    };
    i
}




#[cfg(test)]
mod tests{
    use super::*;
    fn floatcompare(f1:f64, f2:f64)->bool{
        let b = (f1 > f2 - 0.01) & (f1 < f2+0.01);
        b
    }
    #[test]
    fn medtest(){
        let v = vec![1.3, 5.7, 2.2, 9.8, 7.3];
        let v2 = vec![53.4, 6.9, 8.2, 1.1, 17.3, 24.1];
        let m1 = getmedian(&v);
        let m2 = getmedian(&v2);
        let b1 = floatcompare(m1, 5.7);
        let b2 = floatcompare(m2,12.75);
        println!("{m1}");
        println!("{m2}");
        assert_eq!(b1, true);
        assert_eq!(b2, true);
    }

    /*
    #[test]
    fn indextest(){
        let i = closestindex(&vec![1., 3., 7., 5.2], 5.3);
        assert_eq!(i,3)
    }
     */

    #[test]
    fn indextest2(){
        let mut avec:Vec<f64> = Vec::new();
        let space = 0.5;
        let start = 0.5;
        for n in 0..1000{
            avec.push(start + space*(n as f64));
        };
        let i = closestindexordered(&avec, 426.);
        println!("{}",avec[851]);
        println!("{}",avec[i]);
        println!("{i}");
        assert_eq!(i, 851);
        let i2 = closestindexordered(&avec, -5.);
        assert_eq!(i2,0);

        let i3 = closestindexordered(&avec, 52.3);
        println!("{i3}, {}", avec[i3]);
        assert_eq!(i3,104);
        assert!((avec[i3] - 52.3).powi(2) < 0.25);

    }  
        #[test]
        fn getyztest(){
        let fname= String::from("emptyCap_dty138.79_dtz108.00_001_0001p.poni");
        let fpath = Path::new(&fname);
        let fpath2 = Path::new("D:\\data\\July2026\\dty138.79_dtz108.00.poni");
        let ymotor = String::from("dty");
        let zmotor = String::from("dtz");
        let (yo,zo) = getyz(fpath, &ymotor, &zmotor);
        let y = yo.unwrap();
        let z = zo.unwrap();
        assert_eq!(y, 138.79);
        assert_eq!(z, 108.00);
        let (yo,zo) = getyz(fpath2, &ymotor, &zmotor);
        let y = yo.unwrap();
        let z = zo.unwrap();
        assert_eq!(y, 138.79);
        assert_eq!(z, 108.00);
    }

    #[test]
    fn yzcomparetest(){
        let f1 = Path::new("emptyCap_dty138.79_dtz108.00_001_0001p.cbf");
        let f2 = Path::new("dty138.79_dtz108.00_001_0001p.edf");
        let ymotor = String::from("dty");
        let zmotor = String::from("dtz");
        assert!(yzcompare(f1, f2, &ymotor, &zmotor));
    }

}