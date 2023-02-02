use std::io;
use std::io::{Error, ErrorKind};
use crossbeam_channel::{Sender};
use bio::io::fasta;
use bio::io::fasta::{Record};

#[path = "encoding.rs"]
mod encoding;
use encoding::*;

// One encoded fasta record.
#[derive(Clone, Debug,  PartialEq)]
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
            if (self.seq[i] < 240) && (self.seq[i] != other.seq[i]) { // any difference apart from N/-/? is relevant here, not just certain nucleotide differences, because of the triangularity of query vs consensus, target vs consensus, query vs target.
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

fn encode_get_differences(record: &Record, other: &EncodedFastaRecord) -> Result<EncodedFastaRecord, String> {
    
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
        if (efr.seq[i] < 240) && (efr.seq[i] != other.seq[i]) { // any difference apart from N/-/? is relevant here, not just certain nucleotide differences, because of the triangularity of query vs consensus, target vs consensus, query vs target.
            efr.differences.push(i)
        }
    }

    Ok(efr)
}

// Load the records in a list of fasta files into a vector of records in memory.
pub fn load_fasta<T: io::Read>(input: T) -> io::Result<Vec<EncodedFastaRecord>> {
    
    let mut width: usize = 0;
    let mut records: Vec<EncodedFastaRecord> = vec![];

    let reader = fasta::Reader::new(input);

    let mut first = true;

    for r in reader.records() {
        
        let record = r.expect("Error during fasta record parsing");
        let efr = encode(&record).unwrap();

        if first {
            width = efr.seq.len();
            first = false;
        } else if efr.seq.len() != width {
            return Err(Error::new(ErrorKind::Other, "Different length sequences in alignment"))
        }

        records.push(efr);
    }

    Ok(records)
}

// Stream the records in a fasta file by passing them down a channel. Avoids loading the whole file into memory.
pub fn stream_fasta<T: io::Read>(stream: T, loaded: &Vec<Vec<EncodedFastaRecord>>, measure: &str, consen: Option<EncodedFastaRecord>, batchsize: usize, channel: Sender<Records>) {

    let w = loaded[0][0].seq.len();

    let reader = fasta::Reader::new(stream);

    let mut idx_counter = 0;
    let mut batch_counter = 0;
    let mut record_vec: Vec<EncodedFastaRecord> = vec![];

    let mut consensus = EncodedFastaRecord::new_known_width(w);
    match measure {
        "n" => {
            consensus = consen.unwrap()
        }
        _ => ()
    }

    for r in reader.records() {
        
        let record = r.expect("Error during fasta record parsing");
        if record.seq().len() != w {
            panic!("streamed alignment is not the same width as loaded alignment")
        }

        let mut efr = EncodedFastaRecord::new_known_width(w);

        match measure {
            "tn93" => {
                efr = encode_count_bases(&record).unwrap()
            }
            "n" => {
                efr = encode_get_differences(&record, &consensus).unwrap()
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
        136, 72, 40, 24, // [A, G, C, T]
    ];

    let mut consensus: Vec<u8> = vec![0; w];

    for (i, array) in counts.into_iter().enumerate() {
        let mut maxidx = 0;
        let mut maxval = 0;
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

#[cfg(test)]
mod tests {
    use super::*;
    use bio::io::fasta::{Reader};
    use crossbeam_channel::bounded;

    const FASTA: &[u8] = b">target
ATGATGATGATGCCC
";

    const OTHER: &[u8] = b">target
ATTATTATGATGCCC
";

    fn read(fasta: &[u8]) -> Record {
        let reader = Reader::new(fasta);
        let record = reader.records().next().unwrap().unwrap();
        record
    }

    #[test]
    fn test_count_bases() {
        let mut record = encode(&read(FASTA)).unwrap();
        record.count_bases();

        assert_eq!(record.count_A, 4);
        assert_eq!(record.count_T, 4);
        assert_eq!(record.count_C, 3);
        assert_eq!(record.count_G, 4);
    }

    #[test]
    fn test_get_differences() {
        let mut record = encode(&read(FASTA)).unwrap();
        let other = encode(&read(OTHER)).unwrap();

        record.get_differences(&other);

        assert_eq!(record.differences, vec![2,5]);
    }

    #[test]
    fn test_encode() {
        let temp = read(FASTA);
        let record = encode(&temp).unwrap();

        let desired_result: Vec<u8> = vec![136, 24, 72, 136, 24, 72, 136, 24, 72, 136, 24, 72, 40, 40, 40];

        assert_eq!(record.seq, desired_result);
    }

    #[test]
    fn test_encode_count_bases() {
        let temp = read(FASTA);
        let record = encode_count_bases(&temp).unwrap();

        assert_eq!(record.count_A, 4);
        assert_eq!(record.count_T, 4);
        assert_eq!(record.count_C, 3);
        assert_eq!(record.count_G, 4);
    }

    #[test]
    fn test_encode_get_differences() {
        let temp1 = read(OTHER);
        let other = encode(&temp1).unwrap();

        let temp2 = read(FASTA);
        let record = encode_get_differences(&temp2, &other).unwrap();            
        
        assert_eq!(record.differences, vec![2,5]);
    }

    #[test]
    fn test_load_alignment() {
        let record_vec = load_fasta(FASTA).unwrap();
        
        assert_eq!(record_vec.len(), 1);
        assert_eq!(record_vec[0].seq, vec![136, 24, 72, 136, 24, 72, 136, 24, 72, 136, 24, 72, 40, 40, 40]);
    }

    #[test]
    fn test_consensus() {
        let temp1 = read(OTHER);
        let other = encode(&temp1).unwrap();
        let temp2 = read(FASTA);
        let record = encode(&temp2).unwrap();  

        let mut loaded = vec![vec![record.clone(), other.clone()]];
        let mut c = consensus(&loaded);
        
        assert_eq!(c.seq, vec![136, 24, 72, 136, 24, 72, 136, 24, 72, 136, 24, 72, 40, 40, 40]);

        loaded = vec![vec![record.clone(), record.clone()]];
        c = consensus(&loaded);
        
        assert_eq!(c.seq, vec![136, 24, 72, 136, 24, 72, 136, 24, 72, 136, 24, 72, 40, 40, 40]);

        loaded = vec![vec![other.clone(), other.clone()]];
        c = consensus(&loaded);
        
        assert_eq!(c.seq, vec![136, 24, 24, 136, 24, 24, 136, 24, 72, 136, 24, 72, 40, 40, 40]);
    }

    #[test]
    fn test_stream_alignment() {
        let temp1 = read(OTHER);
        let other = encode(&temp1).unwrap();
        let temp2 = read(FASTA);
        let record = encode(&temp2).unwrap();  

        let mut loaded = vec![vec![record, other]];
        let c = consensus(&loaded);

        let (sx, rx) = bounded(1);

        stream_fasta(FASTA, &loaded, "raw", None, 1, sx.clone());
        assert_eq!(rx.recv().unwrap().records[0], loaded[0][0]);
        assert!(rx.is_empty());

        loaded[0][0].count_bases();
        stream_fasta(FASTA, &loaded, "tn93", None, 1, sx.clone());
        assert_eq!(rx.recv().unwrap().records[0], loaded[0][0]);
        assert!(rx.is_empty());

        loaded[0][1].get_differences(&c);
        stream_fasta(OTHER, &loaded, "n", Some(c), 1, sx.clone());
        assert_eq!(rx.recv().unwrap().records[0], loaded[0][1]);
        assert!(rx.is_empty());        
    }

}