use async_std;
use clap::Parser;
use crossbeam_channel::{bounded, Sender, Receiver};
use crossbeam_utils::sync::WaitGroup;
use futures::executor::block_on;
use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::Write;
use std::io::BufWriter;
use std::thread;

mod fastaio;
use crate::fastaio::*;

mod distance;
use crate::distance::tn93;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    #[clap(short, long)]
    threads: usize,
    // input alignment
	#[clap(short, long)]
	input: String,
	#[clap(short, long)]
    output: String,
}

#[derive(Clone)]
struct Pair {
    seq1: EncodedFastaRecord,
    seq2: EncodedFastaRecord,
    idx: usize
}

#[derive(Clone)]
struct Result {
    id1: String,
    id2: String,
    dist: f64,
    idx: usize
}

fn main() -> io::Result<()> {
	let args = Args::parse();

    // let efra = Box::new(populate_struct_array(&args.input)).unwrap(); 
    let efra = populate_struct_array(&args.input)?;

    let (pair_sender, pair_receiver) = bounded(50);
    let (distance_sender, distance_receiver) = bounded(50);

    let wg_write = WaitGroup::new();
    let wg_dist = WaitGroup::new();

    thread::spawn( {
        move || {
            generate_pairs(efra, pair_sender);
        }
    });

    thread::spawn( {
        let wg_write = wg_write.clone();
        move || {
            gather_write(&args.output, distance_receiver);
            drop(wg_write);
        }
    });

    let mut workers = Vec::new();
    for i in 0..args.threads {
        workers.push((pair_receiver.clone(), distance_sender.clone()))
    }

    for worker in workers {
        let wg_dist = wg_dist.clone();
        thread::spawn( move || {
            for message in worker.0.iter() {
                let d = tn93(&message.seq1, &message.seq2);
                worker.1.send(Result{id1: message.seq1.id.clone(), id2: message.seq2.id.clone(), dist: d, idx: message.idx});
            }
            drop(wg_dist);
       });
    }

    wg_dist.wait();
    drop(distance_sender);

    wg_write.wait();

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

fn gather_write(filename: &str, rx: Receiver<Result>) -> io::Result<()> {

    let f = File::create(filename)?;
    let mut buf = BufWriter::new(f);
    writeln!(buf, "sequence1\tsequence2\tdistance")?;

    let mut m: HashMap<usize,Result> = HashMap::new();

    let mut counter: usize = 0;

    for r in rx.iter() {

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
