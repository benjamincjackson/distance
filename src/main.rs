use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::Write;
use std::io::BufWriter;
use crossbeam_channel::{bounded, Sender, Receiver};
use std::thread;
use clap::Parser;

mod fastaio;
use crate::fastaio::*;

mod distance;
use crate::distance::tn93;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
	// input alignment
	#[clap(short, long)]
	input: String,
	#[clap(short, long)]
    output: String
}

struct Pair<'a> {
    seq1: &'a EncodedFastaRecord,
    seq2: &'a EncodedFastaRecord,
    idx: usize
}

struct Result {
    id1: String,
    id2: String,
    dist: f64,
    idx: usize
}

fn main() -> io::Result<()> {
	let args = Args::parse();

    let efra = Box::new(populate_struct_array(&args.input)).unwrap(); 

    let (pair_sender, pair_receiver) = bounded(50);
    let (distance_sender, distance_receiver) = bounded(50);
    
    let future_1 = generate_pairs(efra, pair_sender);
    let future_2 = gather_write(&args.output, distance_receiver);

    Ok(())
}

async fn generate_pairs<'a>(sequences: Vec<EncodedFastaRecord>, sender: Sender<Pair<'a>>) {
    let mut counter: usize = 0;
    for i in 0..sequences.len()-1 {
        for j in i+1..sequences.len() {
            sender.send(Pair{seq1: &sequences[i], seq2: &sequences[j], idx: counter.clone()}).unwrap();
            counter += 1;
        }
    }
}

async fn gather_write(filename: &str, rx: Receiver<Result>) -> io::Result<()> {

    let f = File::create(filename)?;
    let mut buf = BufWriter::new(f);
    writeln!(buf, "sequence1\tsequence2\tdistance")?;

    let mut m: HashMap<usize,Result> = HashMap::new();

    let mut counter: usize = 0;

    for r in rx {

        m.insert(r.idx, r);

        if m.contains_key(&counter) {
            let r = m.remove(&counter).unwrap();
            writeln!(buf, "{}\t{}\t{}", &r.id1, &r.id2, r.dist)?;
            counter += 1;
        }

    }

    while m.len() > 0 {
        let r = m.remove(&counter).unwrap();
        writeln!(buf, "{}\t{}\t{}", &r.id1, &r.id2, r.dist)?;
        counter += 1;
    }

    Ok(())
}
