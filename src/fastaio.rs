use std::fs::File;
use std::io;
use std::io::BufRead;
use std::io::BufReader;
use std::io::{Error, ErrorKind};

#[path = "encoding.rs"]
mod encoding;
use encoding::*;

// One (unencoded) fasta record. Currently unused
pub struct FastaRecord {
    pub id: String,
    pub description: String,
    pub seq: String,
}
impl FastaRecord {
    // Encode a plain fasta record
    pub fn encode(&self) -> EncodedFastaRecord {
        let mut EFR = EncodedFastaRecord::newknownwidth(self.seq.len());
        EFR.id = self.id.clone();
        EFR.description = self.description.clone();

        let a = encoding_array();
        for i in 0..self.seq.len() {
            EFR.seq[i] = a[self.seq.as_bytes()[i] as usize];
        }

        EFR
    }
}

// One encoded fasta record.
#[derive(Clone)]
pub struct EncodedFastaRecord {
    pub id: String,
    pub description: String,
    pub seq: Vec<u8>,
    pub count_A: usize, // the base contents are needed to calculate tn93 distance
    pub count_T: usize,
    pub count_G: usize,
    pub count_C: usize,
    pub differences: Vec<usize>, // differences from the consensus - makes for fast snp-distances for low-diversity datasets
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
            differences: vec![0; 0],
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
            differences: vec![0; 0],
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
    pub fn get_differences(&mut self, other: &EncodedFastaRecord) {
        self.differences.clear();
        for i in 0..self.seq.len() {
            if (self.seq[i] & other.seq[i]) < 16 {
                self.differences.push(i);
            }
        }
    }
}

// How wide is the alignment
pub fn align_width(filename: &str) -> io::Result<usize> {
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

// Read and encode input files to a vector (of vectors of structs)
// First dimension is file, second dimension is sequence record.
pub fn populate_struct_array(files: &Vec<&str>, measure: &str) -> io::Result<Vec<Vec<EncodedFastaRecord>>> {

    let mut widths: Vec<usize> = vec![0; files.len()];
    let mut structs_vec: Vec<Vec<EncodedFastaRecord>> = vec![];
    let ea = encoding_array();

    for (i, file) in files.into_iter().enumerate()  {
        let w = align_width(file)
                    .unwrap();
        widths.push(w);
        if i == 1 {
            if widths[0] != widths[1] {
                panic!("Files have different widths");
            }
        }

        let mut efr = EncodedFastaRecord::newknownwidth(w);
        let mut structs: Vec<EncodedFastaRecord> = vec![];
    
        let mut nuccounter: usize = 0;
        let mut l = String::new();
        let mut first = true;
    
        let f = File::open(file)?;
        let reader = BufReader::new(f);
    
        for line in reader.lines() {
            l = line.unwrap();
    
            if first {
                efr.description = l.strip_prefix('>').unwrap().to_string();
                efr.id = efr
                    .description
                    .split_whitespace()
                    .collect::<Vec<&str>>()[0]
                    .to_string();
                first = false;
                continue;
            } else if l.starts_with('>') {
                nuccounter = 0;
                
                structs.push(efr.clone());
                efr.description = l.strip_prefix('>').unwrap().to_string();
                efr.id = efr
                    .description
                    .split_whitespace()
                    .collect::<Vec<&str>>()[0]
                    .to_string();
                
                continue;
            }
            for nuc in l.bytes() {
                efr.seq[nuccounter] = ea[nuc as usize];
                nuccounter += 1;
            }
        }
    
        structs.push(efr.clone());

        structs_vec.push(structs);
    }

    // Need to do some extra work depending on which distance measure is used.
    match measure {
        // For the fast snp-distance, need to calculate the consensus then get the differences from 
        // it for each record (in each file)
        "n" => {
            let consensus = consensus(&structs_vec);
            for i in 0..structs_vec.len() {
                for j in 0..structs_vec[i].len() {
                    structs_vec[i][j].get_differences(&consensus);
                }
            }
        }
        // For Tamura and Nei (1993), need to calculate the base content of each record.
        "tn93" => {
            for i in 0..structs_vec.len() {
                for j in 0..structs_vec[i].len() {
                    structs_vec[i][j].count_bases();
                }
            }
        }
        _ => (),
    }

    Ok(structs_vec)
}

// Calculate the (ATGC) consensus sequence from all the input data held in memory
pub fn consensus(efras: &Vec<Vec<EncodedFastaRecord>>) -> EncodedFastaRecord {

    // Alignment width
    let w = efras[0][0].seq.len();
    // Counts of ATGC for each alignment column
    let mut counts: Vec<Vec<usize>> = vec![vec![0; 4]; w];

    // for looking up the encoded bases. Non-ATGC characters map to A, which is fine
    // as we don't really care what the consensus sequence actually is, just that it
    // enables us to reduce the computational burden
    let mut lookup: [usize; 256] = [0; 256];
    lookup[136] = 0;
    lookup[72] = 1;
    lookup[40] = 2;
    lookup[24] = 3;

    // For every record, at every alignment column, count which base occurs
    for efra in efras {
        for record in efra {
            for i in 0..record.seq.len() {
                let nuc = record.seq[i];
                counts[i][lookup[nuc as usize]] += 1;
            }
        }
    }

    // For back-translating the counts to (encoded) nucleotides
    let back_translate: [u8; 4] = [
        136, 72, 40, 24,
    ];

    let mut consensus: Vec<u8> = vec![0; w];

    let mut maxidx: usize = 0;
    let mut maxval: usize = 0;

    for (i, array) in counts.into_iter().enumerate() {
        maxidx = 0;
        maxval = 0;
        for (i, val) in array.iter().enumerate() {
            if val > &maxval {
                maxval = *val;
                maxidx = i;
            }
        }
        consensus[i] = back_translate[maxidx]
    }

    let mut EFR = EncodedFastaRecord::newknownwidth(w);
    EFR.seq = consensus;
    
    EFR
}
