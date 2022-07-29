use std::fs::File;
use std::io;
use crossbeam_channel::{Sender};
use bio::io::fasta;
use bio::io::fasta::{Record};

#[path = "encoding.rs"]
mod encoding;
use encoding::*;

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
    pub idx: usize,
}

impl EncodedFastaRecord {
    pub fn new() -> EncodedFastaRecord {
        EncodedFastaRecord {
            id: String::new(),
            description: String::new(),
            seq: vec![0; 0],
            count_A: 0,
            count_T: 0,
            count_C: 0,
            count_G: 0,
            differences: vec![0; 0],
            idx: 0,
        }
    }
    pub fn new_known_width(w: usize) -> EncodedFastaRecord {
        EncodedFastaRecord {
            id: String::new(),
            description: String::new(),
            seq: vec![0; w],
            count_A: 0,
            count_T: 0,
            count_C: 0,
            count_G: 0,
            differences: vec![0; 0],
            idx: 0,
        }
    }
    pub fn count_bases(&mut self) {
        self.count_A = 0;
        self.count_T = 0;
        self.count_G = 0;
        self.count_C = 0;
        let mut counting = [0; 256];
        for i in 0..self.seq.len() {
            counting[self.seq[i] as usize] += 1
        }
        self.count_A = counting[136];
        self.count_T = counting[24];
        self.count_G = counting[72];
        self.count_C = counting[40];
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

#[derive(Clone)]
pub struct Records {
    pub records: Vec<EncodedFastaRecord>,
    pub idx: usize,
}

pub fn encode(record: &Record) -> Result<EncodedFastaRecord, String> {
    
    let ea  = encoding_array();
    let mut efr = EncodedFastaRecord::new_known_width(record.seq().len());

    efr.id = record.id().to_string();
    match record.desc() {
        Some(desc) => efr.description = desc.to_string(),
        None => (),
    }
    
    for (i, nuc) in record.seq().iter().enumerate() {
        if ea[*nuc as usize] == 0 {
            let mut message = "invalid nucleotide character in record: ".to_string();
            message.push(*nuc as char);
            return Err(message)
        }
        efr.seq[i] = ea[*nuc as usize]
    }

    Ok(efr)
}

pub fn encode_count_bases(record: &Record) -> Result<EncodedFastaRecord, String> {
    
    let ea  = encoding_array();
    let mut efr = EncodedFastaRecord::new_known_width(record.seq().len());

    let mut counting = [0; 256];

    efr.id = record.id().to_string();
    match record.desc() {
        Some(desc) => efr.description = desc.to_string(),
        None => (),
    }
    
    for (i, nuc) in record.seq().iter().enumerate() {
        if ea[*nuc as usize] == 0 {
            let mut message = "invalid nucleotide character in record: ".to_string();
            message.push(*nuc as char);
            return Err(message)
        }
        efr.seq[i] = ea[*nuc as usize];
        counting[*nuc as usize] += 1;
    }

    efr.count_A = counting['A' as usize];
    efr.count_T = counting['T' as usize];
    efr.count_G = counting['G' as usize];
    efr.count_C = counting['C' as usize];

    Ok(efr)
}

pub fn encode_get_differences(record: &Record, other: &EncodedFastaRecord) -> Result<EncodedFastaRecord, String> {
    
    let ea  = encoding_array();
    let mut efr = EncodedFastaRecord::new_known_width(record.seq().len());

    efr.id = record.id().to_string();
    match record.desc() {
        Some(desc) => efr.description = desc.to_string(),
        None => (),
    }
    
    for (i, nuc) in record.seq().iter().enumerate() {
        if ea[*nuc as usize] == 0 {
            let mut message = "invalid nucleotide character in record: ".to_string();
            message.push(*nuc as char);
            return Err(message)
        }
        efr.seq[i] = ea[*nuc as usize];
        if (efr.seq[i] & other.seq[i]) < 16 {
            efr.differences.push(i)
        }
    }

    Ok(efr)
}

// Load the records in a list of fasta files into a vector of vector of records in memory.
pub fn load_fastas(files: &Vec<String>, measure: &str) -> io::Result<Vec<Vec<EncodedFastaRecord>>> {
    
    let mut widths: Vec<usize> = vec![0; files.len()];
    let mut structs_vec: Vec<Vec<EncodedFastaRecord>> = vec![];

    for (i, file) in files.into_iter().enumerate() {

        let mut structs: Vec<EncodedFastaRecord> = vec![];

        let f = File::open(file).unwrap();
        let reader = fasta::Reader::new(f);

        let mut first = true;

        for r in reader.records() {
            
            let record = r.expect("Error during fasta record parsing");
            let efr = encode(&record).unwrap();

            if first {
                widths.push(efr.seq.len());
                if i == 1 {
                    if widths[0] != widths[1] {
                        panic!("Files have different widths");
                    }
                }
                first = false
            }

            structs.push(efr.clone());
        }
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

// Stream the records in a fasta file by passing them down a channel. Avoids loading the whole file into memory.
pub fn stream_fasta(file: &str, loaded: &Vec<Vec<EncodedFastaRecord>>, measure: &str, batchsize: usize, channel: Sender<Records>) {

    let mut c = EncodedFastaRecord::new();
    if measure == "n" {
       c = consensus(&loaded);
    }

    let w = loaded[0][0].seq.len();

    let f = File::open(file).unwrap();
    let reader = fasta::Reader::new(f);

    let mut idx_counter = 0;
    let mut batch_counter = 0;
    let mut record_vec: Vec<EncodedFastaRecord> = vec![];

    for r in reader.records() {
        
        let record = r.expect("Error during fasta record parsing");
        if record.seq().len() != w {
            panic!("streamed alignment is not the same width as loaded alignment")
        }

        let mut efr = EncodedFastaRecord::new_known_width(w);

        match measure {
            "n" => {
                efr = encode_get_differences(&record, &c).unwrap()
            }
            "tn93" => {
                efr = encode_count_bases(&record).unwrap()
            }
            _ => efr = encode(&record).unwrap()
        }

        record_vec.push(efr);
        batch_counter += 1;

        if batch_counter == batchsize {
            channel
                .send(Records{
                        records: record_vec.clone(), 
                        idx: idx_counter,
                    })
                .unwrap();
            
            idx_counter += 1;
            batch_counter = 0;
            record_vec.clear();
        }
    }

    // send the last batch
    if record_vec.len() > 0 {
        channel
            .send(Records{
                records: record_vec.clone(), 
                idx: idx_counter,
            })
        .unwrap()
    }

    drop(channel)
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

    let mut efr = EncodedFastaRecord::new_known_width(w);
    efr.seq = consensus;
    
    efr
}
