use clap::{Command, Arg};
use crossbeam_channel::{bounded, Receiver, Sender};
use crossbeam_utils::sync::WaitGroup;
use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::{BufWriter, Write};
use std::thread;

mod fastaio;
use crate::fastaio::*;

mod distance;
use crate::distance::*;

#[derive(Clone)]
struct Pair {
    seq1: EncodedFastaRecord,
    seq2: EncodedFastaRecord,
    idx: usize,
}
impl Pair {
    fn distance(&self, measure: &str) -> FloatInt {
        match measure {
            "raw" => raw(&self.seq1, &self.seq2),
            "n" => snp(&self.seq1, &self.seq2),
            "jc69" => jc69(&self.seq1, &self.seq2),
            "k80" => k80(&self.seq1, &self.seq2),
            "tn93" => tn93(&self.seq1, &self.seq2),
            // should never get this far because the options are defined in the clap cli:
            _ => panic!("unknown distance measure"),
        }
    }
}

#[derive(Clone)]
struct Distance {
    id1: String,
    id2: String,
    dist: FloatInt,
    idx: usize,
}

fn main() -> io::Result<()> {

    let m = Command::new("distance")
        .version("0.1.0")
        .arg(Arg::new("threads")
            .short('t')
            .long("threads")
            .takes_value(true)
            .help("how many threads to spin up for pairwise comparisons"))
        .arg(Arg::new("input")
            .short('i')
            .long("input")
            .help("input alignment(s) in fasta format")
            .takes_value(true)
            .multiple_values(true)
            .max_values(2)
            .required(true))
        .arg(Arg::new("measure")
            .short('m')
            .long("measure")
            .takes_value(true)
            .default_value("raw")
            .possible_values(["n", "raw", "jc69", "k80", "tn93"])
            .help("which distance measure to use"))
        .arg(Arg::new("output")
            .short('o')
            .long("output")
            .takes_value(true)
            .help("output file in tab-separated-value format")
            .required(true))
        .get_matches();

    let inputs: Vec<&str> = m.values_of("input").unwrap().collect();

    let (pair_sender, pair_receiver) = bounded(50);
    let (distance_sender, distance_receiver) = bounded(50);

    let wg_dist = WaitGroup::new();

    let measure = m
        .value_of("measure")
        .unwrap()
        .to_owned();

    if inputs.len() == 1 {
        let mut efra = populate_struct_array(inputs[0])?;
        thread::spawn({
            move || {
                generate_pairs_square(efra, pair_sender);
            }
        });
    } else {
        let mut efra1 = populate_struct_array(inputs[0])?;
        let mut efra2 = populate_struct_array(inputs[1])?;
        thread::spawn({
            move || {
                generate_pairs_rect(efra1, efra2, pair_sender);
            }
        });
    }

    let output = m
        .value_of("output")
        .unwrap()
        .to_owned();

    let write = thread::spawn({
        move || {
            gather_write(&output, distance_receiver);
        }
    });

    let threads = m
        .value_of("threads")
        .unwrap()
        .parse::<usize>()
        .unwrap();

    let mut workers = Vec::new();
    for _i in 0..threads {
        workers.push((pair_receiver.clone(), distance_sender.clone()))
    }    

    for worker in workers {
        let wg_dist = wg_dist.clone();
        let measure = measure.clone();
        thread::spawn(move || {
            for message in worker.0.iter() {
                let d = message.distance(&measure);
                worker.1.send(Distance {
                    id1: message.seq1.id.clone(),
                    id2: message.seq2.id.clone(),
                    dist: d,
                    idx: message.idx,
                });
            }
            drop(wg_dist);
        });
    }

    wg_dist.wait();
    drop(distance_sender);

    write.join().unwrap();

    Ok(())
}

fn generate_pairs_square(sequences: Vec<EncodedFastaRecord>, sender: Sender<Pair>) {
    let mut counter: usize = 0;
    for i in 0..sequences.len() - 1 {
        for j in i + 1..sequences.len() {
            sender
                .send(Pair {
                    seq1: sequences[i].clone(),
                    seq2: sequences[j].clone(),
                    idx: counter,
                })
                .unwrap();
            counter += 1;
        }
    }
    drop(sender);
}

fn generate_pairs_rect(sequences1: Vec<EncodedFastaRecord>, sequences2: Vec<EncodedFastaRecord>, sender: Sender<Pair>) {
    let mut counter: usize = 0;
    for i in 0..sequences1.len() {
        for j in 0..sequences2.len() {
            sender
                .send(Pair {
                    seq1: sequences1[i].clone(),
                    seq2: sequences2[j].clone(),
                    idx: counter,
                })
                .unwrap();
            counter += 1;
        }
    }
    drop(sender);
}

fn gather_write(filename: &str, rx: Receiver<Distance>) -> io::Result<()> {
    let f = File::create(filename)?;
    let mut buf = BufWriter::new(f);
    writeln!(buf, "sequence1\tsequence2\tdistance")?;

    let mut m: HashMap<usize, Distance> = HashMap::new();

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

    while !m.is_empty() {
        let r = m.remove(&counter).unwrap();
        match r.dist {
            FloatInt::Int(d) => writeln!(buf, "{}\t{}\t{}", &r.id1, &r.id2, d)?,
            FloatInt::Float(d) => writeln!(buf, "{}\t{}\t{}", &r.id1, &r.id2, d)?,
        }
        counter += 1;
    }

    Ok(())
}
