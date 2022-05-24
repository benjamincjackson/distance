use clap::Parser;
use crossbeam_channel::{bounded, Sender, Receiver};
use crossbeam_utils::sync::WaitGroup;
use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::{Write, BufWriter};
use std::thread;

mod fastaio;
use crate::fastaio::*;

mod distance;
use crate::distance::*;

#[derive(Parser, Debug, Clone)]
#[clap(author, version, about, long_about = None)]
struct Args {
    #[clap(short, long, default_value = "raw", possible_values = ["raw", "n", "jc69", "tn93"], help = "which distance measure to use")]
    measure: String,
    #[clap(short, long, help = "how many threads to spin up for pairwise comparisons")]
    threads: usize,
    // input alignment
	#[clap(short, long, help = "input alignment in fasta format")]
	input: String,
	#[clap(short, long, help = "output file in tab-separated-value format")]
    output: String,
}

#[derive(Clone)]
struct Pair {
    seq1: EncodedFastaRecord,
    seq2: EncodedFastaRecord,
    idx: usize
}
impl Pair {
    fn distance(&self, measure: &String) -> FloatInt {
        match measure.as_str() {
            "n" => snp(&self.seq1, &self.seq2),
            "raw" => raw(&self.seq1, &self.seq2),
            "jc69" => jc69(&self.seq1, &self.seq2),
            "tn93" => tn93(&self.seq1, &self.seq2),
            // TO DO - learn Rust's error handling properly and apply it here
            _ => panic!("unknown distance measure")
        }
    }
}

#[derive(Clone)]
struct Distance {
    id1: String,
    id2: String,
    dist: FloatInt,
    idx: usize
}


fn main() -> io::Result<()> {
	let args = Args::parse();

    let efra = populate_struct_array(&args.input)?;

    let (pair_sender, pair_receiver) = bounded(50);
    let (distance_sender, distance_receiver) = bounded(50);

    let wg_dist = WaitGroup::new();

    thread::spawn( {
        move || {
            generate_pairs(efra, pair_sender);
        }
    });

    let write = thread::spawn( {
        move || {
            gather_write(&args.output, distance_receiver);
        }
    });

    let mut workers = Vec::new();
    for i in 0..args.threads {
        workers.push((pair_receiver.clone(), distance_sender.clone()))
    }

    for worker in workers {
        let wg_dist = wg_dist.clone();
        let measure = args.measure.clone();
        thread::spawn( move || {
            for message in worker.0.iter() {
                let d = message.distance(&measure);
                worker.1.send(Distance{id1: message.seq1.id.clone(), id2: message.seq2.id.clone(), dist: d, idx: message.idx});
            }
            drop(wg_dist);
       });
    }

    wg_dist.wait();
    drop(distance_sender);

    write.join().unwrap();

    Ok(())
}

fn generate_pairs(sequences: Vec<EncodedFastaRecord>, sender: Sender<Pair>) {
    let mut counter: usize = 0;
    for i in 0..sequences.len()-1 {
        for j in i+1..sequences.len() {
            sender.send(Pair{seq1: sequences[i].clone(), seq2: sequences[j].clone(), idx: counter.clone()}).unwrap();
            counter += 1;
        }
    }
    drop(sender);
}

fn gather_write(filename: &str, rx: Receiver<Distance>) -> io::Result<()> {

    let f = File::create(filename)?;
    let mut buf = BufWriter::new(f);
    writeln!(buf, "sequence1\tsequence2\tdistance")?;

    let mut m: HashMap<usize,Distance> = HashMap::new();

    let mut counter: usize = 0;

    for r in rx.iter() {

        m.insert(r.idx, r);

        if m.contains_key(&counter) {
            let r = m.remove(&counter).unwrap();
            match r.dist {
                FloatInt::Int(d) => writeln!(buf, "{}\t{}\t{}", &r.id1, &r.id2, d)?,
                FloatInt::Float(d) => writeln!(buf, "{}\t{}\t{}", &r.id1, &r.id2, d)?,
            }
            counter += 1;
        }

    }

    while m.len() > 0 {
        let r = m.remove(&counter).unwrap();
        match r.dist {
            FloatInt::Int(d) => writeln!(buf, "{}\t{}\t{}", &r.id1, &r.id2, d)?,
            FloatInt::Float(d) => writeln!(buf, "{}\t{}\t{}", &r.id1, &r.id2, d)?,
        }
        counter += 1;
    }

    Ok(())
}
