use bio::io::fasta;
use bio::io::fasta::Record;
use crossbeam_channel::Sender;
use std::io;

use crate::{DistanceError, Result};

#[path = "encoding.rs"]
mod encoding;
use encoding::*;

// One encoded fasta record.
#[derive(Clone, Debug, PartialEq)]
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
            seq: Vec::new(),
            count_A: 0,
            count_T: 0,
            count_C: 0,
            count_G: 0,
            differences: Vec::new(),
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
            differences: Vec::new(),
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
            if (self.seq[i] < 240) && (self.seq[i] != other.seq[i]) {
                // any difference apart from N/-/? is relevant here, not just certain nucleotide differences, because of the triangularity of query vs consensus, target vs consensus, query vs target.
                self.differences.push(i);
            }
        }
    }
}
impl Default for EncodedFastaRecord {
    fn default() -> Self {
            Self::new()
        }
    }

#[derive(Clone)]
pub struct Records {
    pub records: Vec<EncodedFastaRecord>,
    pub idx: usize,
}

fn err_message_invalid_nuc(r: String, c: char) -> String {
    format!("Invalid nucleotide character in record '{}': '{}'", r, c)
}

fn err_message_different_length_seqs(w1: usize, w2: usize) -> String {
    format!("Different length sequences in alignment(s): {} vs {}", w1, w2)
}

fn err_message_empty_fasta_file() -> String {
    "Empty FASTA file".to_string()
}

pub fn encode(record: &Record) -> Result<EncodedFastaRecord> {
    let encoding_array = encoding_array();
    let mut encoded_record = EncodedFastaRecord::new_known_width(record.seq().len());

    encoded_record.id = record.id().to_string();
    if let Some(desc) = record.desc() {
        encoded_record.description = desc.to_string()
    }

    for (i, nuc) in record.seq().iter().enumerate() {
        if encoding_array[*nuc as usize] == 0 {
            return Err(DistanceError::Message(err_message_invalid_nuc(record.id().to_string(), *nuc as char)));
        }
        encoded_record.seq[i] = encoding_array[*nuc as usize]
    }

    Ok(encoded_record)
}

pub fn encode_count_bases(record: &Record) -> Result<EncodedFastaRecord> {
    let encoding_array = encoding_array();
    let mut encoded_record = EncodedFastaRecord::new_known_width(record.seq().len());

    let mut counting = [0; 256];

    encoded_record.id = record.id().to_string();
    if let Some(desc) = record.desc() {
        encoded_record.description = desc.to_string()
    }

    for (i, nuc) in record.seq().iter().enumerate() {
        if encoding_array[*nuc as usize] == 0 {
            return Err(DistanceError::Message(err_message_invalid_nuc(record.id().to_string(), *nuc as char)));
        }
        encoded_record.seq[i] = encoding_array[*nuc as usize];
        counting[*nuc as usize] += 1;
    }

    encoded_record.count_A = counting['A' as usize];
    encoded_record.count_T = counting['T' as usize];
    encoded_record.count_G = counting['G' as usize];
    encoded_record.count_C = counting['C' as usize];

    Ok(encoded_record)
}

fn encode_get_differences(
    record: &Record,
    other: &EncodedFastaRecord,
) -> Result<EncodedFastaRecord> {
    let encoding_array = encoding_array();
    let mut encoded_record = EncodedFastaRecord::new_known_width(record.seq().len());

    encoded_record.id = record.id().to_string();
    if let Some(desc) = record.desc() {
        encoded_record.description = desc.to_string()
    }

    for (i, nuc) in record.seq().iter().enumerate() {
        if encoding_array[*nuc as usize] == 0 {
            return Err(DistanceError::Message(err_message_invalid_nuc(record.id().to_string(), *nuc as char)));
        }
        encoded_record.seq[i] = encoding_array[*nuc as usize];
        if (encoded_record.seq[i] < 240) && (encoded_record.seq[i] != other.seq[i]) {
            // any difference apart from N/-/? is relevant here, not just certain nucleotide differences, because of the triangularity of query vs consensus, target vs consensus, query vs target.
            encoded_record.differences.push(i)
        }
    }

    Ok(encoded_record)
}

// Load the records in a list of fasta files into a vector of records in memory.
fn load_fasta<T: io::Read>(input: T) -> Result<Vec<EncodedFastaRecord>> {
    let mut width: usize = 0;
    let mut records: Vec<EncodedFastaRecord> = vec![];
    let mut first = true;

    let reader = fasta::Reader::new(input);

    for r in reader.records() {
        let record = r?;
        let encoded_record = encode(&record)?;

        if first {
            width = encoded_record.seq.len();
            first = false;
        } else if encoded_record.seq.len() != width {
            return Err(DistanceError::Message(err_message_different_length_seqs(encoded_record.seq.len(), width)));
        }

        records.push(encoded_record);
    }

    if records.is_empty() {
        return Err(DistanceError::Message(err_message_empty_fasta_file()));
    }

    Ok(records)
}

pub fn load_fastas<T: io::Read>(inputs: Vec<T>) -> Result<Vec<Vec<EncodedFastaRecord>>> {
    let mut loaded = vec![];
    for (counter, file) in inputs.into_iter().enumerate() {
        loaded.push(load_fasta(file)?);
        if counter == 1 && loaded[0][0].seq.len() != loaded[1][0].seq.len() {
            return Err(DistanceError::Message(err_message_different_length_seqs(loaded[0][0].seq.len(), loaded[1][0].seq.len())));
        }
    }

    Ok(loaded)
}

// Stream the records in a fasta file by passing them down a channel. Avoids loading the whole file into memory.
pub fn stream_fasta<T: io::Read>(
    stream: T,
    loaded: &[Vec<EncodedFastaRecord>],
    measure: &str,
    consen: Option<EncodedFastaRecord>,
    batchsize: usize,
    channel: Sender<Records>,
) -> Result<()> {
    let width: usize = loaded[0][0].seq.len();

    let reader = fasta::Reader::new(stream);

    let mut idx_counter: usize = 0;
    let mut batch_counter: usize = 0;
    let mut record_counter: usize = 0;
    let mut record_vec: Vec<EncodedFastaRecord> = vec![];

    let mut consensus = EncodedFastaRecord::new_known_width(width);
    if measure == "n" {
        match consen {
            None => return Err(DistanceError::Message(
                "Expected a consensus sequence to be generated when the distance measure is n".to_string(),
            )),
            Some(encoded_record) => consensus = encoded_record,
        }
    }

    for r in reader.records() {
        let record = r?;
        record_counter += 1;

        if record.seq().len() != width {
            return Err(DistanceError::Message(err_message_different_length_seqs(record.seq().len(), width)));
        }

        let encoded_record: EncodedFastaRecord = match measure {
            "tn93" => encode_count_bases(&record)?,
            "n" => encode_get_differences(&record, &consensus)?,
            _ => encode(&record)?,
        };

        record_vec.push(encoded_record);
        batch_counter += 1;

        if batch_counter == batchsize {
            channel.send(Records {
                records: record_vec.clone(),
                idx: idx_counter,
            })?;

            idx_counter += 1;
            batch_counter = 0;
            record_vec.clear();
        }
    }

    // send the last batch
    if !record_vec.is_empty() {
        channel.send(Records {
            records: record_vec.clone(),
            idx: idx_counter,
        })?;
    }

    drop(channel);

    if record_counter == 0 {
        return Err(DistanceError::Message(err_message_empty_fasta_file()))
    }

    Ok(())
}

// Calculate the (ATGC) consensus sequence from all the input data held in memory
pub fn consensus(encoded_records_vecs: &Vec<Vec<EncodedFastaRecord>>) -> EncodedFastaRecord {
    // Alignment width
    let width: usize = encoded_records_vecs[0][0].seq.len();
    // Counts of ATGC for each alignment column
    let mut counts: Vec<Vec<usize>> = vec![vec![0; 4]; width];

    // for looking up the encoded bases. Non-ATGC characters map to A, which is fine
    // as we don't really care what the consensus sequence actually is, just that it
    // enables us to reduce the computational burden
    let mut lookup: [usize; 256] = [0; 256];
    lookup[136] = 0;
    lookup[72] = 1;
    lookup[40] = 2;
    lookup[24] = 3;

    // For every record, at every alignment column, count which base occurs
    for encoded_records in encoded_records_vecs {
        for encoded_record in encoded_records {
            for (i, nuc) in encoded_record.seq.iter().enumerate() {
                counts[i][lookup[*nuc as usize]] += 1;
            }
        }
    }

    // For back-translating the counts to (encoded) nucleotides
    let back_translate: [u8; 4] = [
        136, 72, 40, 24, // [A, G, C, T]
    ];

    let mut consensus: Vec<u8> = vec![0; width];

    for (i, array) in counts.into_iter().enumerate() {
        let mut maxidx: usize = 0;
        let mut maxval: usize = 0;
        for (i, val) in array.iter().enumerate() {
            if val > &maxval {
                maxval = *val;
                maxidx = i;
            }
        }
        consensus[i] = back_translate[maxidx]
    }

    let mut encoded_record = EncodedFastaRecord::new_known_width(width);
    encoded_record.seq = consensus;

    encoded_record
}

#[cfg(test)]
mod tests {
    use super::*;
    use bio::io::fasta::Reader;
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

        assert_eq!(record.differences, vec![2, 5]);
    }

    #[test]
    fn test_encode() {
        let temp = read(FASTA);
        let record = encode(&temp).unwrap();

        let desired_result: Vec<u8> = vec![
            136, 24, 72, 136, 24, 72, 136, 24, 72, 136, 24, 72, 40, 40, 40,
        ];

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

        assert_eq!(record.differences, vec![2, 5]);
    }

    #[test]
    fn test_load_alignment() {
        let record_vec = load_fasta(FASTA).unwrap();

        assert_eq!(record_vec.len(), 1);
        assert_eq!(
            record_vec[0].seq,
            vec![136, 24, 72, 136, 24, 72, 136, 24, 72, 136, 24, 72, 40, 40, 40]
        );
    }

    #[test]
    fn test_consensus() {
        let temp1 = read(OTHER);
        let other = encode(&temp1).unwrap();
        let temp2 = read(FASTA);
        let record = encode(&temp2).unwrap();

        let mut loaded = vec![vec![record.clone(), other.clone()]];
        let mut c = consensus(&loaded);

        assert_eq!(
            c.seq,
            vec![136, 24, 72, 136, 24, 72, 136, 24, 72, 136, 24, 72, 40, 40, 40]
        );

        loaded = vec![vec![record.clone(), record.clone()]];
        c = consensus(&loaded);

        assert_eq!(
            c.seq,
            vec![136, 24, 72, 136, 24, 72, 136, 24, 72, 136, 24, 72, 40, 40, 40]
        );

        loaded = vec![vec![other.clone(), other.clone()]];
        c = consensus(&loaded);

        assert_eq!(
            c.seq,
            vec![136, 24, 24, 136, 24, 24, 136, 24, 72, 136, 24, 72, 40, 40, 40]
        );
    }

    #[test]
    fn test_stream_alignment() -> Result<()> {
        let temp1 = read(OTHER);
        let other = encode(&temp1).unwrap();
        let temp2 = read(FASTA);
        let record = encode(&temp2).unwrap();

        let mut loaded = vec![vec![record, other]];
        let c = consensus(&loaded);

        let (sx, rx) = bounded(1);

        let _ = stream_fasta(FASTA, &loaded, "raw", None, 1, sx.clone())?;

        assert_eq!(rx.recv().unwrap().records[0], loaded[0][0]);
        assert!(rx.is_empty());

        loaded[0][0].count_bases();
        let _ = stream_fasta(FASTA, &loaded, "tn93", None, 1, sx.clone())?;
        assert_eq!(rx.recv().unwrap().records[0], loaded[0][0]);
        assert!(rx.is_empty());

        loaded[0][1].get_differences(&c);
        let _ = stream_fasta(OTHER, &loaded, "n", Some(c), 1, sx.clone())?;
        assert_eq!(rx.recv().unwrap().records[0], loaded[0][1]);
        assert!(rx.is_empty());

        Ok(())
    }
}
