use clap::{Command, Arg, ArgMatches};
use crossbeam_channel::{bounded, Receiver, Sender};
use crossbeam_utils::sync::WaitGroup;
use std::collections::HashMap;
use std::io;
use std::io::{BufWriter, Write};
use std::fs::File;
use std::sync::Arc;
use std::thread;


mod measures;
use crate::measures::*;

mod fastaio;
use crate::fastaio::*;

// A struct for passing the location of one pairwise comparison down a channel (between threads)
#[derive(Clone, Debug, PartialEq)]
struct Pair {
    seq1_idx: usize,
    seq2_idx: usize,
}

#[derive(Clone, Debug, PartialEq)]
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

pub fn get_cli_arguments() -> ArgMatches {
    
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

        m
}

pub struct Setup {
    inputs: Vec<String>,
    stream: String,
    measure: String,
    threads: usize,
    batchsize: usize,
    output: String,
    distances_channel: (Sender<Distances>, Receiver<Distances>),
    wg_dist: WaitGroup,
    efras: Vec<Vec<EncodedFastaRecord>>,
    ns: Vec<usize>,
}
impl Setup {
    fn new() -> Setup {
        Setup{
            inputs: Vec::new(),
            stream: String::new(),
            measure: String::new(),
            threads: 1,
            batchsize: 1,
            output: String::new(),
            distances_channel: bounded(50),
            wg_dist: WaitGroup::new(),
            efras: vec![vec![EncodedFastaRecord::new();0];0],
            ns: Vec::new(),
        }
    }
}

pub fn set_up(m: &ArgMatches) -> Setup {

    let mut setup = Setup::new();

    // One or two input fasta file names
    setup.inputs = m
        .values_of("input")
        .unwrap()
        .map(|s| String::from(s))
        .collect::<Vec<String>>();

    if m.value_of("stream").is_some() {
        setup.stream = m
        .value_of("stream")
        .unwrap()
        .to_owned();
    }

    if setup.stream.len() != 0 {
        if setup.inputs.len() != 1 {
            eprintln!("If you stream one file, you must also provide exactly one other file to -i/--input");
            std::process::exit(1);
        }
    }

    // Which distance measure to use
    setup.measure = m
        .value_of("measure")
        .unwrap()
        .to_string();

    // batch size - to tune the workload per message so that threads aren't fighting over pairs_receiver's lock as often
    setup.batchsize = m
        .value_of("batchsize")
        .unwrap()
        .parse::<usize>()
        .unwrap(); 

    // The aligned sequence data
    setup.efras = load_fastas(&setup.inputs, &setup.measure)
        .unwrap();

    // The sample sizes (for generating the pairs)
    for efra in &setup.efras {
        setup.ns.push(efra.len());
    }

    // The output file name
    setup.output = m
        .value_of("output")
        .unwrap()
        .to_owned();

    // How many additional threads to use for calculating distances - need at least 1.
    setup.threads = m
        .value_of("threads")
        .unwrap()
        .parse::<usize>()
        .unwrap();

    if setup.threads < 1 {
        setup.threads = 1;
    }

    setup

}

pub fn stream(setup: Setup) {

    let arc = Arc::new(setup.efras);
    
    let (records_sender, records_receiver) = bounded(50);
    let (distances_sender, distances_receiver) = setup.distances_channel;

    let batchsize = setup.batchsize.to_owned();
    let output = setup.output.to_owned();
    let stream = setup.stream.to_owned();
    let measure = setup.measure.to_owned();

    // We spin up a thread to write the output as it arrives down the distance channel.
    let write = thread::spawn({
        move || {
            gather_write(&output, distances_receiver).unwrap();
        }
    });
    
    thread::spawn({
        let measure = measure.clone();
        let arc = arc.clone();
        move || {
            stream_fasta(&stream, &arc, &measure, batchsize, records_sender);
        }
    });

    // A vector of receiver/sender tuples, cloned to share between each thread in threads.
    let mut workers = Vec::new();
    for _i in 0..setup.threads {
        workers.push((records_receiver.clone(), distances_sender.clone()))
    }

    // Spin up the threads that do the distance-calculating
    for worker in workers {
        let wg_dist = setup.wg_dist.clone();
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

    // When all the distances have been calculated, we can drop the sending end of the distance channel
    setup.wg_dist.wait();
    drop(distances_sender);

    // Joins when all the pairwise comparisons have been written, and then we're done.
    write.join().unwrap();

}

pub fn load(setup: Setup) {

    let arc = Arc::new(setup.efras);

    let (pairs_sender, pairs_receiver) = bounded(50);
    let (distances_sender, distances_receiver) = setup.distances_channel;

    let inputs = setup.inputs.to_owned();
    let ns = setup.ns.to_owned();
    let batchsize = setup.batchsize.to_owned();
    let output = setup.output.to_owned();
    let measure = setup.measure.to_owned();

    // We spin up a thread to write the output as it arrives down the distance channel.
    let write = thread::spawn({
        move || {
            gather_write(&output, distances_receiver).unwrap();
        }
    });

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
    for _i in 0..setup.threads {
        workers.push((pairs_receiver.clone(), distances_sender.clone()))
    }

    // Spin up the threads that do the distance-calculating
    for worker in workers {
        let wg_dist = setup.wg_dist.clone();
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

    // When all the distances have been calculated, we can drop the sending end of the distance channel
    setup.wg_dist.wait();
    drop(distances_sender);

    // Joins when all the pairwise comparisons have been written, and then we're done.
    write.join().unwrap();

}

pub fn run(setup: Setup) {
    if setup.stream.len() > 0 {
        stream(setup)
    } else {
        load(setup)
    }
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


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_generate_pairs_square() {
        let n: usize = 4;
        let batch_size = 1;
        let (sx, rx) = bounded(50);

        thread::spawn({
            let sx = sx.clone();
            move || {
                generate_pairs_square(n, batch_size, sx);
            }
        });

        assert_eq!(rx.recv().unwrap(), Pairs{pairs: vec![Pair{seq1_idx: 0, seq2_idx: 1}], idx: 0 });
        assert_eq!(rx.recv().unwrap(), Pairs{pairs: vec![Pair{seq1_idx: 0, seq2_idx: 2}], idx: 1 });
        assert_eq!(rx.recv().unwrap(), Pairs{pairs: vec![Pair{seq1_idx: 0, seq2_idx: 3}], idx: 2 });
        assert_eq!(rx.recv().unwrap(), Pairs{pairs: vec![Pair{seq1_idx: 1, seq2_idx: 2}], idx: 3 });
        assert_eq!(rx.recv().unwrap(), Pairs{pairs: vec![Pair{seq1_idx: 1, seq2_idx: 3}], idx: 4 });
        assert_eq!(rx.recv().unwrap(), Pairs{pairs: vec![Pair{seq1_idx: 2, seq2_idx: 3}], idx: 5 });
        assert!(rx.is_empty());

        let batch_size_2 = 4;

        thread::spawn({
            let sx = sx.clone();
            move || {
                generate_pairs_square(n, batch_size_2, sx);
            }
        });

        assert_eq!(rx.recv().unwrap(), Pairs{pairs: vec![Pair{seq1_idx: 0, seq2_idx: 1}, Pair{seq1_idx: 0, seq2_idx: 2}, Pair{seq1_idx: 0, seq2_idx: 3}, Pair{seq1_idx: 1, seq2_idx: 2}], idx: 0 });
        assert_eq!(rx.recv().unwrap(), Pairs{pairs: vec![Pair{seq1_idx: 1, seq2_idx: 3}, Pair{seq1_idx: 2, seq2_idx: 3}], idx: 1 });
        assert!(rx.is_empty());
    }

    #[test]
    fn test_generate_pairs_rect() {
        let n1: usize = 2;
        let n2: usize = 2;
        let batch_size = 1;
        let (sx, rx) = bounded(50);

        thread::spawn({
            let sx = sx.clone();
            move || {
                generate_pairs_rect(n1, n2, batch_size, sx);
            }
        });

        assert_eq!(rx.recv().unwrap(), Pairs{pairs: vec![Pair{seq1_idx: 0, seq2_idx: 0}], idx: 0 });
        assert_eq!(rx.recv().unwrap(), Pairs{pairs: vec![Pair{seq1_idx: 0, seq2_idx: 1}], idx: 1 });
        assert_eq!(rx.recv().unwrap(), Pairs{pairs: vec![Pair{seq1_idx: 1, seq2_idx: 0}], idx: 2 });
        assert_eq!(rx.recv().unwrap(), Pairs{pairs: vec![Pair{seq1_idx: 1, seq2_idx: 1}], idx: 3 });
        assert!(rx.is_empty());

        let batch_size_2 = 4;

        thread::spawn({
            let sx = sx.clone();
            move || {
                generate_pairs_rect(n1, n2, batch_size_2, sx);
            }
        });

        assert_eq!(rx.recv().unwrap(), Pairs{pairs: vec![Pair{seq1_idx: 0, seq2_idx: 0}, Pair{seq1_idx: 0, seq2_idx: 1}, Pair{seq1_idx: 1, seq2_idx: 0}, Pair{seq1_idx: 1, seq2_idx: 1}], idx: 0 });
        assert!(rx.is_empty());
    }
}