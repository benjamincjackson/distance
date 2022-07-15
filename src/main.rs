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
#[derive(Clone)]
struct Pair {
    seq1_idx: usize,
    seq2_idx: usize,
}

#[derive(Clone)]
struct Pairs {
    pairs: Vec<Pair>,
    idx: usize,
}

// A struct for passing one pair's distance down a channel (between threads)
#[derive(Clone)]
struct Distance {
    id1: String,
    id2: String,
    dist: FloatInt,
}

#[derive(Clone)]
struct Distances {
    distances: Vec<Distance>,
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
            .help("How many threads to spin up for pairwise comparisons"))
        .arg(Arg::new("input")
            .short('i')
            .long("input")
            .help("Input alignment file(s) in fasta format. Loaded into memory")
            .takes_value(true)
            .multiple_values(true)
            .max_values(2)
            .required_unless_present("licenses"))
        .arg(Arg::new("stream")
            .short('s')
            .long("stream")
            .takes_value(true)
            .multiple_values(false)
            .help("Input alignment file in fasta format. Streamed from disk. Requires exactly one file also be specifed to -i"))
        .arg(Arg::new("measure")
            .short('m')
            .long("measure")
            .takes_value(true)
            .default_value("raw")
            .possible_values(["n", "n_high", "raw", "jc69", "k80", "tn93"])
            .help("Which distance measure to use"))
        .arg(Arg::new("output")
            .short('o')
            .long("output")
            .takes_value(true)
            .help("Output file in tab-separated-value format")
            .required_unless_present("licenses"))
        .arg(Arg::new("batchsize")
            .long("batchsize")
            .short('b')
            .default_value("1")
            .help("Try setting this >(>) 1 if you are struggling to get a speedup when adding threads"))
        .arg(Arg::new("licenses")
            .long("licenses")
            .short('l')
            .takes_value(false)
            .help("Print licence information and exit"))
        .get_matches();

    if m.is_present("licenses") {
        println!("{}", licences());
        std::process::exit(0);
    }

    // One or two input fasta file names
    let inputs: Vec<&str> = m.values_of("input").unwrap().collect();

    // mpmc channel for passing pairs of distances between threads
    let (distances_sender, distances_receiver) = bounded(50);

    // We will use this waitgroup to make sure all possible pair-distances are calculated before 
    // dropping the sending end of the distance channel
    let wg_dist = WaitGroup::new();

    // Which distance measure to use
    let measure = m
        .value_of("measure")
        .unwrap()
        .to_owned();

    // batch size - to tune the workload per message so that threads aren't fighting over pairs_receiver's lock as often
    let batchsize = m
        .value_of("batchsize")
        .unwrap()
        .parse::<usize>()
        .unwrap(); 

    // The aligned sequence data
    let efras = load_fastas(&inputs, &measure)
        .unwrap();

    // The sample sizes (for generating the pairs)
    let mut ns = vec![];
    for efra in &efras {
        ns.push(efra.len());
    }

    // The output file name
    let output = m
        .value_of("output")
        .unwrap()
        .to_owned();

    // We spin up a thread to write the output as it arrives down the distance channel.
    let write = thread::spawn({
        move || {
            gather_write(&output, distances_receiver).unwrap();
        }
    });

    // How many additional threads to use for calculating distances - need at least 1.
    let mut threads = m
        .value_of("threads")
        .unwrap()
        .parse::<usize>()
        .unwrap();

    if threads < 1 {
        threads = 1;
    }

    // Wrap the loaded data in an Arc so that we can share it between threads.
    let arc = Arc::new(efras);

    // if a stream is specified, 
    if m.values_of("stream").is_some() {

        // then the user must specify exactly one input file
        if inputs.len() != 1 {
            eprintln!("If you stream one file, you must also provide exactly one other file to -i/--input");
            std::process::exit(1);
        }
        
        // and then we do streamy stuff...
                
        let stream = m
            .value_of("stream")
            .unwrap()
            .to_owned();
        
        // an mpmc channel for passing the streamed fasta records between threads
        let (record_sender, record_receiver) = bounded(50);

        // we spin up a thread to read the streamed input
        thread::spawn({
            let measure = measure.clone();
            let arc = arc.clone();
            move || {
                stream_fasta(&stream, &arc, &measure, batchsize, record_sender);
            }
        });

        // A vector of receiver/sender tuples, cloned to share between each thread in threads.
        let mut workers = Vec::new();
        for _i in 0..threads {
            workers.push((record_receiver.clone(), distances_sender.clone()))
        }

        // Spin up the threads that do the distance-calculating
        for worker in workers {
            let wg_dist = wg_dist.clone();
            let measure = measure.clone();
            let arc = arc.clone();       
            thread::spawn(move || {
                let mut distances: Vec<Distance> = vec![];
                for message in worker.0.iter() {
                    for target in message.records {
                        for query in &arc[0] {
                            let d = match measure.as_str() {
                                "raw" => raw(query, &target),
                                "n" => snp2(query, &target),
                                "n_high" => snp(query, &target),
                                "jc69" => jc69(query, &target),
                                "k80" => k80(query, &target),
                                "tn93" => tn93(query, &target),
                                // should never get this far because the options are defined in the cli:
                                _ => panic!("unknown distance measure"),
                            };
                            // add the ids and the distance to this batch's temporary vector
                            distances.push(Distance{
                                id1: target.id.clone(),
                                id2: query.id.clone(),
                                dist: d,
                            });                        
                        }
                    }
                    // send this batch of distances to the writer
                    worker.1
                        .send(Distances{
                            distances: distances.clone(),
                            idx: message.idx,
                            })
                        .unwrap();
                    // clear the vector ready for the next batch
                    distances.clear();
                }

                // when the pair channel is empty (and all distances are calculated) we can drop the cloned distance waitgroup (for this thread)
                drop(wg_dist);
            });
        }
    } else {
        // otherwise all the data is in memory, so we can do things a little differently

        // an mpmc channel for passing pairs of sequences between threads
        let (pairs_sender, pairs_receiver) = bounded(50);

        // If there is one input file, generate all pairwise comparisons within the alignment,
        // else there are two input files, so generate all pairwise comparisons between the alignments.
        // We spin up a thread do to this and move on to the next part of the program
        if inputs.len() == 1 {
            thread::spawn({
                move || {
                    generate_pairs_square(ns[0], batchsize, pairs_sender);
                }
            });
        } else {
            thread::spawn({
                move || {
                    generate_pairs_rect(ns[0], ns[1], batchsize, pairs_sender);
                }
            });
        }

        // A vector of receiver/sender tuples, cloned to share between each thread in threads.
        let mut workers = Vec::new();
        for _i in 0..threads {
            workers.push((pairs_receiver.clone(), distances_sender.clone()))
        }

        // Spin up the threads that do the distance-calculating
        for worker in workers {
            let wg_dist = wg_dist.clone();
            let measure = measure.clone();
            let arc = arc.clone();
            thread::spawn(move || {
                let mut distances: Vec<Distance> = vec![];
                // for each batch
                for message in worker.0.iter() {
                    // for each pair in this batch
                    for pair in message.pairs {
                        // calculate the distance
                        let d = match measure.as_str() {
                            "raw" => raw(&arc[0][pair.seq1_idx], &arc[arc.len()-1][pair.seq2_idx]),
                            "n" => snp2(&arc[0][pair.seq1_idx], &arc[arc.len()-1][pair.seq2_idx]),
                            "n_high" => snp(&arc[0][pair.seq1_idx], &arc[arc.len()-1][pair.seq2_idx]),
                            "jc69" => jc69(&arc[0][pair.seq1_idx], &arc[arc.len()-1][pair.seq2_idx]),
                            "k80" => k80(&arc[0][pair.seq1_idx], &arc[arc.len()-1][pair.seq2_idx]),
                            "tn93" => tn93(&arc[0][pair.seq1_idx], &arc[arc.len()-1][pair.seq2_idx]),
                            // should never get this far because the options are defined in the cli:
                            _ => panic!("unknown distance measure"),
                        };
                        // add the ids and the distance to this batch's temporary vector
                        distances.push(Distance{
                            id1: arc[0][pair.seq1_idx].id.clone(),
                            id2: arc[arc.len()-1][pair.seq2_idx].id.clone(),
                            dist: d,
                        });
                    }

                    // send this batch of distances to the writer
                    worker.1
                        .send(Distances{
                            distances: distances.clone(),
                            idx: message.idx,
                        })
                    .unwrap();

                    // clear the vector ready for the next batch
                    distances.clear();
                }

                // when the pair channel is empty (and all distances are calculated) we can drop the cloned distance waitgroup (for this thread)
                drop(wg_dist);
            });
        }
    }

    // When all the distances have been calculated, we can drop the sending end of the distance channel
    wg_dist.wait();
    drop(distances_sender);

    // Joins when all the pairwise comparisons have been written, and then we're done.
    write.join().unwrap();

    Ok(())
}

// Given the sample size of a single alignment, generate all possible pairwise comparisons
// within it, and pass them down a channel.
fn generate_pairs_square(n: usize, size: usize, sender: Sender<Pairs>) {
    
    // this counter is used to send the correct-sized batch
    let mut size_counter: usize = 0;

    // this counter is sent down the channel in order to later retain input order in the output
    let mut idx_counter: usize = 0;

    let mut pair_vec: Vec<Pair> = vec![];

    for i in 0..n - 1 {
        for j in i + 1..n {

            pair_vec
                .push(Pair {
                    seq1_idx: i,
                    seq2_idx: j,
                });

            size_counter += 1;

            // when we reach the batch size, we send this batch of pairs
            if size_counter == size {

                sender
                    .send( Pairs {
                        pairs: pair_vec.clone(),
                        idx: idx_counter,
                    })
                    .unwrap();

                size_counter = 0;
                idx_counter += 1;
                pair_vec.clear();
            }
        }
    }

    // send the last batch
    if pair_vec.len() > 0 {
        sender
            .send( Pairs {
                pairs: pair_vec.clone(),
                idx: idx_counter,
            })
        .unwrap();
    }
    
    drop(sender);
}

// Given the sample size of two alignments, generate all possible pairwise comparisons
// between them, and pass them down a channel.
fn generate_pairs_rect(n1: usize, n2: usize, size: usize, sender: Sender<Pairs>) {

    // this counter is used to send the correct-sized batch
    let mut size_counter: usize = 0;

    // this counter is sent down the channel in order to later retain input order in the output
    let mut idx_counter: usize = 0;

    let mut pair_vec: Vec<Pair> = vec![];
    
    for i in 0..n1 {
        for j in 0..n2 {
            pair_vec
                .push(Pair {
                    seq1_idx: i,
                    seq2_idx: j,
                });

            size_counter += 1;

            // when we reach the batch size, we send this batch of pairs
            if size_counter == size {

                sender
                    .send( Pairs {
                        pairs: pair_vec.clone(),
                        idx: idx_counter,
                    })
                    .unwrap();

                size_counter = 0;
                idx_counter += 1;
                pair_vec.clear();
            }
        }
    }

    // send the last batch
    if pair_vec.len() > 0 {
        sender
            .send( Pairs {
                pairs: pair_vec.clone(),
                idx: idx_counter,
            })
        .unwrap();
    }

    drop(sender);
}

// Write the distances as they arrive. Uses a hashmap whos keys are indices to write the results in the
// order they are produced by generate_pairs_*()
fn gather_write(filename: &str, rx: Receiver<Distances>) -> io::Result<()> {
    
    let f = File::create(filename)?;
    let mut buf = BufWriter::new(f);
    writeln!(buf, "sequence1\tsequence2\tdistance")?;

    let mut m: HashMap<usize, Distances> = HashMap::new();

    let mut counter: usize = 0;

    for r in rx.iter() {
        m.insert(r.idx, r);
        while m.contains_key(&counter) {
            let rv = m.remove(&counter).unwrap();
            for result in rv.distances {
                match result.dist {
                    FloatInt::Int(d) => writeln!(buf, "{}\t{}\t{}", &result.id1, &result.id2, d)?,
                    FloatInt::Float(d) => writeln!(buf, "{}\t{}\t{}", &result.id1, &result.id2, d)?,
                }
            }
            counter += 1;            
        }
    }

    Ok(())
}

fn licences() -> String {
    let licences = "
Copyright 2022, Ben Jackson

distance is licensed under the GNU LIBRARY GENERAL PUBLIC LICENSE, Version 2

distance incorporates some of Rust-Bio, which is licensed under the MIT licence:

The MIT License (MIT)

Copyright (c) 2016 Johannes Köster, the Rust-Bio team, Google Inc.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and 
associated documentation files (the \"Software\"), to deal in the Software without restriction, 
including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, 
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial 
portions of the Software.

THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT 
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES 
OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
".to_string();

    licences
}