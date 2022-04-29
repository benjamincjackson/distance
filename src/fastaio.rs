use std::io;
use std::fs::File;
use std::io::BufReader;
use std::io::BufRead;
use std::io::{Error, ErrorKind};

#[path = "encoding.rs"]
mod encoding;
use encoding::*;

pub struct FastaRecord {
    pub id: String,
    pub description: String,
    pub seq: String,
}

pub fn read_fasta(filename: &String) -> io::Result<()> {
    let f = File::open(filename)?;
    let reader = BufReader::new(f);

    let mut first = true;
    let mut id = String::with_capacity(80);
    let mut description = String::with_capacity(80);
    let mut seq = String::new();
    let mut l = String::new();

    // let mut counter = 0;

    for line in reader.lines() {
        l = line.unwrap();
        if first {
            if l.starts_with(">") {
                description = l.strip_prefix(">").unwrap().to_string();
                id = description.split_whitespace().collect::<Vec<&str>>()[0].to_string();
                first = false
            } else {
                return Err(Error::new(ErrorKind::Other, "badly formatted fasta file"));
            }
            continue
        }
        if l.starts_with(">") {
            let FR: FastaRecord = FastaRecord { id: id.clone(), description: description.clone(), seq: seq.clone() };
            // counter+=1;
            // println!("{}", counter);
            println!(">{}", FR.id);
            // println!("{}", FR.description);
            println!("{}", FR.seq);
            description = l.strip_prefix(">").unwrap().to_string();
            id = description.split_whitespace().collect::<Vec<&str>>()[0].to_string();
            seq = "".to_string()
        } else {
            seq.push_str(&l);
        }
    }

    let FR: FastaRecord = FastaRecord { id: id.clone(), description: description.clone(), seq: seq.clone() };
    // println!("{}", counter);
    println!(">{}", FR.id);
    // println!("{}", FR.description);
    println!("{}", FR.seq);

    Ok(())
}

pub fn align_dims(filename: &String) -> io::Result<(usize, usize)> {
    let f = File::open(filename)?;
    let reader = BufReader::new(f);

    let mut n = 0;
    let mut w = 0;
    let mut l = String::new();

    for line in reader.lines() {
        l = line.unwrap();
        
        if l.chars().nth(0).unwrap() == '>' {
			n += 1;
		}

        if n == 1 && l.chars().nth(0).unwrap() != '>' {
            w += l.chars().count();
        }
    }

    Ok((w, n))
}

pub fn populate_array(filename: &String) -> io::Result<(Vec<Vec<u8>>, Vec<String>)> {
    
    let a = encoding_array();

    let (w, n) = align_dims(filename).unwrap();

    let mut byte_array: Vec<Vec<u8>> = vec![vec![0; w]; n];
    let mut fasta_ids: Vec<String> = vec!["".to_string(); n];

    let mut n:usize = 0;
    let mut counter:usize = 0;
    let mut l = String::new();
    let mut id = String::with_capacity(80);
    let mut description = String::with_capacity(80);
    let mut first = true;

    let f = File::open(filename)?;
    let reader = BufReader::new(f);

    for line in reader.lines() {
        l = line.unwrap();

        if first {
            description = l.strip_prefix(">").unwrap().to_string();
            id = description.split_whitespace().collect::<Vec<&str>>()[0].to_string();
            fasta_ids[n as usize] = id;
            first = false;
            continue
        } else if l.chars().nth(0).unwrap() == '>' {
            n += 1;
            counter = 0;
            description = l.strip_prefix(">").unwrap().to_string();
            id = description.split_whitespace().collect::<Vec<&str>>()[0].to_string();
            fasta_ids[n] = id;
            continue
		}
        for nuc in l.bytes() {
            byte_array[n][counter] = a[nuc as usize];
            counter += 1;
        }
    }

   Ok((byte_array, fasta_ids))
}

pub fn fasta_consensus(filename: &String) -> io::Result<String> {

    let mut consensus = String::new();

    let a = encoding_array();

    let f = File::open(filename)?;
    let reader = BufReader::new(f);

    let mut firstline: bool = true;
    let mut firstrec: bool = true;
    let mut l = String::new();
    let mut counter = 0;
    let mut w = 0;
  
    let mut a: [u8; 256] = [5; 256];
    a[136] = 1; // => 'A'
    a[72] = 2; // => 'G'
    a[40] = 3; // => 'C'
    a[24] = 4; // => 'T'

    for line in reader.lines() {
        l = line.unwrap();

        if firstline {
            if l.starts_with(">") {
                firstline = false
            } else {
                return Err(Error::new(ErrorKind::Other, "badly formatted fasta file"));
            }
            continue
        } else if l.chars().nth(0).unwrap() == '>' {
            if firstrec {
                w = counter;
                firstrec = false;
                let mut counts: Vec<Vec<u8>> = vec![vec![0;5]; w];
            }
            counter = 0;
            continue
		}
        for nuc in l.bytes() {
            counter += 1;
        }
    }

    Ok(consensus)
}