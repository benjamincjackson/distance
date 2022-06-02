use std::fs::File;
use std::io;
use std::io::BufRead;
use std::io::BufReader;
use std::io::{Error, ErrorKind};

#[path = "encoding.rs"]
mod encoding;
use encoding::*;

pub struct FastaRecord {
    pub id: String,
    pub description: String,
    pub seq: String,
}

#[derive(Clone)]
pub struct EncodedFastaRecord {
    pub id: String,
    pub description: String,
    pub seq: Vec<u8>,
    pub count_A: usize,
    pub count_T: usize,
    pub count_G: usize,
    pub count_C: usize,
}

impl EncodedFastaRecord {
    fn new() -> EncodedFastaRecord {
        EncodedFastaRecord {
            id: String::new(),
            description: String::new(),
            seq: vec![0; 0],
            count_A: 0,
            count_T: 0,
            count_C: 0,
            count_G: 0,
        }
    }
    fn newknownwidth(w: usize) -> EncodedFastaRecord {
        EncodedFastaRecord {
            id: String::new(),
            description: String::new(),
            seq: vec![0; w],
            count_A: 0,
            count_T: 0,
            count_C: 0,
            count_G: 0,
        }
    }
    fn count_bases(&mut self) {
        for i in 0..self.seq.len() {
            match self.seq[i] {
                136 => self.count_A += 1,
                24 => self.count_T += 1,
                72 => self.count_G += 1,
                40 => self.count_C += 1,
                _ => continue,
            }
        }
    }
}

pub fn align_dims(filename: &String) -> io::Result<(usize, usize)> {
    let f = File::open(filename)?;
    let reader = BufReader::new(f);

    let mut n = 0;
    let mut w = 0;
    let mut l = String::new();

    for line in reader.lines() {
        l = line.unwrap();

        if l.starts_with('>') {
            n += 1;
        }

        if n == 1 && !l.starts_with('>') {
            w += l.chars().count();
        }
    }

    Ok((w, n))
}

pub fn align_width(filename: &String) -> io::Result<usize> {
    let f = File::open(filename)?;
    let reader = BufReader::new(f);

    let mut w = 0;
    let mut n = 0;
    let mut l = String::new();

    for line in reader.lines() {
        l = line.unwrap();

        if l.starts_with('>') {
            n += 1;
        }

        if n == 2 {
            break;
        }

        if n == 1 && !l.starts_with('>') {
            w += l.chars().count();
        }
    }

    Ok(w)
}

pub fn populate_array(filename: &String) -> io::Result<(Vec<Vec<u8>>, Vec<String>)> {
    let a = encoding_array();

    let (w, n) = align_dims(filename).unwrap();

    let mut byte_array: Vec<Vec<u8>> = vec![vec![0; w]; n];
    let mut fasta_ids: Vec<String> = vec!["".to_string(); n];

    let mut ncounter = 0;
    let mut counter: usize = 0;
    let mut l = String::new();
    let mut id = String::with_capacity(80);
    let mut description = String::with_capacity(80);
    let mut first = true;

    let f = File::open(filename)?;
    let reader = BufReader::new(f);

    for line in reader.lines() {
        l = line.unwrap();

        if first {
            description = l.strip_prefix('>').unwrap().to_string();
            id = description.split_whitespace().collect::<Vec<&str>>()[0].to_string();
            fasta_ids[ncounter as usize] = id;
            first = false;
            continue;
        } else if l.starts_with('>') {
            ncounter += 1;
            counter = 0;
            description = l.strip_prefix('>').unwrap().to_string();
            id = description.split_whitespace().collect::<Vec<&str>>()[0].to_string();
            fasta_ids[ncounter] = id;
            continue;
        }
        for nuc in l.bytes() {
            byte_array[ncounter][counter] = a[nuc as usize];
            counter += 1;
        }
    }

    Ok((byte_array, fasta_ids))
}

pub fn populate_struct_array(filename: &String) -> io::Result<Vec<EncodedFastaRecord>> {
    let v = encoding_array();

    let (w, n) = align_dims(filename).unwrap();

    let mut structs: Vec<EncodedFastaRecord> = vec![EncodedFastaRecord::newknownwidth(w); n];

    let mut ncounter = 0;
    let mut nuccounter: usize = 0;
    let mut l = String::new();
    let mut first = true;

    let f = File::open(filename)?;
    let reader = BufReader::new(f);

    for line in reader.lines() {
        l = line.unwrap();

        if first {
            structs[ncounter].description = l.strip_prefix('>').unwrap().to_string();
            structs[ncounter].id = structs[ncounter]
                .description
                .split_whitespace()
                .collect::<Vec<&str>>()[0]
                .to_string();
            first = false;
            continue;
        } else if l.starts_with('>') {
            ncounter += 1;
            nuccounter = 0;
            structs[ncounter].description = l.strip_prefix('>').unwrap().to_string();
            structs[ncounter].id = structs[ncounter]
                .description
                .split_whitespace()
                .collect::<Vec<&str>>()[0]
                .to_string();
            continue;
        }
        for nuc in l.bytes() {
            structs[ncounter].seq[nuccounter] = v[nuc as usize];
            nuccounter += 1;
        }
    }

    for i in 0..structs.len() {
        structs[i].count_bases();
    }

    Ok(structs)
}

pub fn fasta_consensus(filename: &String) -> io::Result<String> {
    let w = align_width(filename)?;
    let mut counts: Vec<Vec<usize>> = vec![vec![0; 17]; w];

    let mut lookup: [usize; 256] = [17; 256];
    lookup['A' as usize] = 0;
    lookup['a' as usize] = 0;
    lookup['G' as usize] = 1;
    lookup['g' as usize] = 1;
    lookup['C' as usize] = 2;
    lookup['c' as usize] = 2;
    lookup['T' as usize] = 3;
    lookup['t' as usize] = 3;
    lookup['R' as usize] = 4;
    lookup['r' as usize] = 4;
    lookup['M' as usize] = 5;
    lookup['m' as usize] = 5;
    lookup['W' as usize] = 6;
    lookup['w' as usize] = 6;
    lookup['S' as usize] = 7;
    lookup['s' as usize] = 7;
    lookup['K' as usize] = 8;
    lookup['k' as usize] = 8;
    lookup['Y' as usize] = 9;
    lookup['y' as usize] = 9;
    lookup['V' as usize] = 10;
    lookup['v' as usize] = 10;
    lookup['H' as usize] = 11;
    lookup['h' as usize] = 11;
    lookup['D' as usize] = 12;
    lookup['d' as usize] = 12;
    lookup['B' as usize] = 13;
    lookup['b' as usize] = 13;
    lookup['N' as usize] = 14;
    lookup['n' as usize] = 14;
    lookup['-' as usize] = 15;
    lookup['?' as usize] = 16;

    let f = File::open(filename)?;
    let reader = BufReader::new(f);

    let mut firstline = true;
    let mut l = String::new();
    let mut nuccounter = 0;

    for line in reader.lines() {
        l = line.unwrap();

        if firstline {
            if l.starts_with('>') {
                firstline = false
            } else {
                return Err(Error::new(ErrorKind::Other, "badly formatted fasta file"));
            }
            continue;
        } else if l.starts_with('>') {
            if nuccounter != w {
                return Err(Error::new(
                    ErrorKind::Other,
                    "different length sequences, is this an alignment?",
                ));
            }
            nuccounter = 0;
            continue;
        }
        for nuc in l.bytes() {
            counts[nuccounter][lookup[nuc as usize]] += 1;
            nuccounter += 1;
            if nuccounter > w {
                return Err(Error::new(
                    ErrorKind::Other,
                    "different length sequences, is this an alignment?",
                ));
            }
        }
    }

    let back_translate: [char; 17] = [
        'A', 'G', 'C', 'T', 'R', 'M', 'W', 'S', 'K', 'Y', 'V', 'H', 'D', 'B', 'N', '-', '?',
    ];

    let mut consensus = String::new();

    let mut maxidx: usize = 0;
    let mut maxval: usize = 0;

    for array in counts {
        maxidx = 0;
        maxval = 0;
        for (i, val) in array.iter().enumerate() {
            if val > &maxval {
                maxval = *val;
                maxidx = i;
            }
        }
        consensus.push(back_translate[maxidx])
    }

    Ok(consensus)
}
