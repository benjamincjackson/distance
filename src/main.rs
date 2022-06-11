use clap::{Command, Arg};
use crossbeam_channel::{bounded, Receiver, Sender};
use crossbeam_utils::sync::WaitGroup;
use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::{BufWriter, Write};
use std::sync::Arc;
use std::thread;

mod fastaio;
use crate::fastaio::*;

mod distance;
use crate::distance::*;

// A struct for passing the location of one pairwise comparison down a channel (between threads)
struct Pair {
    seq1_idx: usize,
    seq2_idx: usize,
    idx_counter: usize,
}

// A struct for passing one pair's distance down a channel (between threads)
struct Distance {
    id1: String,
    id2: String,
    dist: FloatInt,
    idx: usize,
}

fn main() -> io::Result<()> {

    // Define the command-line interface
    let m = Command::new("distance")
        .version("0.1.0")
        .arg(Arg::new("threads")
            .short('t')
            .long("threads")
            .takes_value(true)
            .default_value("1")
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
            .possible_values(["n", "n_high", "raw", "jc69", "k80", "tn93"])
            .help("which distance measure to use"))
        .arg(Arg::new("output")
            .short('o')
            .long("output")
            .takes_value(true)
            .help("output file in tab-separated-value format")
            .required(true))
        .get_matches();

    // One or two input fasta file names
    let inputs: Vec<&str> = m.values_of("input").unwrap().collect();

    // mpmc channels for passing pairs of sequences and their distances between threads
    let (pair_sender, pair_receiver) = bounded(50);
    let (distance_sender, distance_receiver) = bounded(50);

    // We will use this waitgroup to make sure all possible pair-distances are calculated before 
    // dropping the sending end of the distance channel
    let wg_dist = WaitGroup::new();

    // Which distance measure to use
    let measure = m
        .value_of("measure")
        .unwrap()
        .to_owned();

    // The aligned sequence data
    let efras = populate_struct_array(&inputs, &measure)
                .unwrap();

    // The sample sizes (for generating the pairs)
    let mut ns = vec![];
    for efra in &efras {
        ns.push(efra.len());
    }

    // If there is one input file, generate all pairwise comparisons within the alignment,
    // else there are two input files, so generate all pairwise comparisons between the alignments.
    // We spin up a thread do to this and move on to the next part of the program
    if inputs.len() == 1 {
        thread::spawn({
            move || {
                generate_pairs_square(ns[0], pair_sender);
            }
        });
    } else {
        thread::spawn({
            move || {
                generate_pairs_rect(ns[0], ns[1], pair_sender);
            }
        });
    }

    // The output file name
    let output = m
        .value_of("output")
        .unwrap()
        .to_owned();

    // We spin up a thread to write the output as it arrives down the distance channel.
    let write = thread::spawn({
        move || {
            gather_write(&output, distance_receiver);
        }
    });

    // How many additional threads to use for calculating distances - need at least 1.
    let mut threads = m
        .value_of("threads")
        .unwrap()
        .parse::<usize>()
        .unwrap();

    if threads < 1 {
        let threads = 1;
    }

    // A vector of receiver/sender tuples, cloned to share between each thread in threads.
    let mut workers = Vec::new();
    for _i in 0..threads {
        workers.push((pair_receiver.clone(), distance_sender.clone()))
    }

    // Wrap the data in an Arc so that we can share it between threads.
    let arc = Arc::new(efras);

    // Spin up the threads to do the actual work
    for worker in workers {
        let wg_dist = wg_dist.clone();
        let measure = measure.clone();
        let arc = arc.clone();
        thread::spawn(move || {
            for message in worker.0.iter() {
                // calculate the distance
                let d = match measure.as_str() {
                    "raw" => raw(&arc[0][message.seq1_idx], &arc[arc.len()-1][message.seq2_idx]),
                    "n" => snp2(&arc[0][message.seq1_idx], &arc[arc.len()-1][message.seq2_idx]),
                    "n_high" => snp(&arc[0][message.seq1_idx], &arc[arc.len()-1][message.seq2_idx]),
                    "jc69" => jc69(&arc[0][message.seq1_idx], &arc[arc.len()-1][message.seq2_idx]),
                    "k80" => k80(&arc[0][message.seq1_idx], &arc[arc.len()-1][message.seq2_idx]),
                    "tn93" => tn93(&arc[0][message.seq1_idx], &arc[arc.len()-1][message.seq2_idx]),
                    // should never get this far because the options are defined in the cli:
                    _ => panic!("unknown distance measure"),
                };
                // send the distance the the writer
                worker.1.send(Distance {
                    id1: arc[0][message.seq1_idx].id.clone(),
                    id2: arc[arc.len()-1][message.seq2_idx].id.clone(),
                    dist: d,
                    idx: message.idx_counter,
                })
                .unwrap();
            }
            // when the pair channel is empty (and all distances are calculated) we can drop the cloned distance waitgroup (for this thread)
            drop(wg_dist);
        });
    }

    // When all the distances have been calculated, we can drop the sending end of the distance channel
    wg_dist.wait();
    drop(distance_sender);

    // Joins when all the pairwise comparisons have been written, and then we're done.
    write.join().unwrap();

    Ok(())
}

// Given the sample size of a single alignment, generate all possible pairwise comparisons
// within it, and pass them down a channel.
fn generate_pairs_square(n: usize, sender: Sender<Pair>) {
    // we send counter down the channel in order to later retain input order in the output
    let mut counter: usize = 0;
    for i in 0..n - 1 {
        for j in i + 1..n {
            sender
                .send(Pair {
                    seq1_idx: i,
                    seq2_idx: j,
                    idx_counter: counter,
                })
                .unwrap();
            counter += 1;
        }
    }
    drop(sender);
}

// Given the sample size of two alignments, generate all possible pairwise comparisons
// between them, and pass them down a channel.
fn generate_pairs_rect(n1: usize, n2: usize, sender: Sender<Pair>) {
    let mut counter: usize = 0;
    for i in 0..n1 {
        for j in 0..n2 {
            sender
                .send(Pair {
                    seq1_idx: i,
                    seq2_idx: j,
                    idx_counter: counter,
                })
                .unwrap();
            counter += 1;
        }
    }
    drop(sender);
}

// Write the distances as they arrive. Uses a hashmap whos keys are indices to write the results in the
// order they are produced by generate_pairs_*()
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
