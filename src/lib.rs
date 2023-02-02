use clap::{Command, Arg, ArgMatches, crate_version};
use crossbeam_channel::{bounded, Receiver, Sender};
use crossbeam_utils::sync::WaitGroup;
use std::collections::HashMap;
use std::io;
use std::io::{BufWriter};
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
        .version(crate_version!())
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
            .default_value("stdin"))
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
            .help("Output file in tab-separated-value format. Omit this option to print to stdout"))
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
    stream: Option<Box<dyn io::Read + Send>>,
    measure: String,
    threads: usize,
    batchsize: usize,
    writer: BufWriter<Box<dyn io::Write + Send>>,
    distances_channel: (Sender<Distances>, Receiver<Distances>),
    wg_dist: WaitGroup,
    efras: Vec<Vec<EncodedFastaRecord>>,
    ns: Vec<usize>,
    consensus: Option<EncodedFastaRecord>,
}
impl Setup {
    fn new() -> Setup {
        Setup{
            stream: None,
            measure: String::new(),
            threads: 1,
            batchsize: 1,
            writer: BufWriter::new(Box::new(io::stdout())),
            distances_channel: bounded(100),
            wg_dist: WaitGroup::new(),
            efras: vec![vec![EncodedFastaRecord::new();0];0],
            ns: Vec::new(),
            consensus: None,
        }
    }
}

pub fn set_up(m: &ArgMatches) -> Setup {

    let mut setup = Setup::new();

    // One or two input fasta file names (or stdin)
    let mut inputs: Vec<Box<dyn std::io::Read>> = Vec::new();
    let raw_inputs: Vec<&str> = m.values_of("input").unwrap().into_iter().collect();
    for s in &raw_inputs {
        if s == &"stdin" {
            inputs.push(Box::new(io::stdin()))
        } else {
            inputs.push(Box::new(File::open(s).unwrap()))
        }
    }

    if m.value_of("stream").is_some() {
        if inputs.len() != 1 || (inputs.len() == 1 && raw_inputs[0] == "stdin") {
            eprintln!("If you stream one file from disk, you must also provide exactly one other file to -i/--input");
            std::process::exit(1);
        }
        setup.stream = Some(
            Box::new(
                File::open(
                    m
                        .value_of("stream")
                        .unwrap()
                )
                .unwrap()
            )
        );
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
    for file in inputs {
        setup.efras.push(load_fasta(file).unwrap())
    }

    // Need to do some extra work depending on which distance measure is used.
    match setup.measure.as_str() {
        // For the fast snp-distance, need to calculate the consensus then get the differences from 
        // it for each record (in each file)
        "n" => {
            let consensus = consensus(&setup.efras);
            for i in 0..setup.efras.len() {
                for j in 0..setup.efras[i].len() {
                    setup.efras[i][j].get_differences(&consensus);
                }
            }
            setup.consensus = Some(consensus);
        }
        // For Tamura and Nei (1993), need to calculate the base content of each record.
        "tn93" => {
            for i in 0..setup.efras.len() {
                for j in 0..setup.efras[i].len() {
                    setup.efras[i][j].count_bases();
                }
            }
        }
        _ => (),
    }

    // The sample sizes (for generating the pairs)
    for efra in &setup.efras {
        setup.ns.push(efra.len());
    }

    if m.value_of("output").is_some() {
        setup.writer = BufWriter::new(Box::new(File::create(m.value_of("output").unwrap()).unwrap()))
    }

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
    
    let (records_sender, records_receiver) = bounded(100);
    let (distances_sender, distances_receiver) = setup.distances_channel;

    let batchsize = setup.batchsize.to_owned();
    let measure = setup.measure.to_owned();

    // We spin up a thread to write the output as it arrives down the distance channel.
    let write = thread::spawn({
        move || {
            gather_write(setup.writer, distances_receiver).unwrap();
        }
    });
    
    thread::spawn({
        let measure = measure.clone();
        let arc = arc.clone();
        let s = setup.stream.expect("Stream not specified");
        move || {
            stream_fasta(s, &arc, &measure, setup.consensus, batchsize, records_sender);
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

            // when the target channel has been dropped (by stream_fasta() after it reaches the end of the file) we can drop the cloned distance waitgroup (for this thread)
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

    let (pairs_sender, pairs_receiver) = bounded(100);
    let (distances_sender, distances_receiver) = setup.distances_channel;

    let ns = setup.ns.to_owned();
    let batchsize = setup.batchsize.to_owned();
    let measure = setup.measure.to_owned();

    // We spin up a thread to write the output as it arrives down the distance channel.
    let write = thread::spawn({
        move || {
            gather_write(setup.writer, distances_receiver).unwrap();
        }
    });

    // If there is one input file, generate all pairwise comparisons within the alignment,
    // else there are two input files, so generate all pairwise comparisons between the alignments.
    // We spin up a thread do to this and move on to the next part of the program
    if ns.len() == 1 {
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
    if setup.stream.is_some() {
        stream(setup);
    } else {
        load(setup);
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

// Write the distances as they arrive. Uses a hashmap whose keys are indices to write the results in the
// order they are produced by generate_pairs_*()
fn gather_write<T: io::Write>(mut writer: T, rx: Receiver<Distances>) -> io::Result<()> {
    
    writeln!(writer, "sequence1\tsequence2\tdistance")?;

    let mut m: HashMap<usize, Distances> = HashMap::new();

    let mut counter: usize = 0;

    for r in rx.iter() {
        m.insert(r.idx, r);
        while m.contains_key(&counter) {
            let rv = m.remove(&counter).unwrap();
            for result in rv.distances {
                match result.dist {
                    FloatInt::Int(d) => writeln!(writer, "{}\t{}\t{}", &result.id1, &result.id2, d)?,
                    FloatInt::Float(d) => writeln!(writer, "{}\t{}\t{}", &result.id1, &result.id2, d)?,
                }
            }
            counter += 1;            
        }
    }

    writer.flush()?;

    Ok(())
}


#[cfg(test)]
mod tests {
    use std::io::{BufReader, BufRead, Read};

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

// const LOAD: &[u8] = b">seq1
// ATGATG
// >seq2
// ATGATT
// ";

    // #[test]
    // fn test_load() {

    //     let mut buffer: Vec<u8> = Vec::new();

    //     thread::scope(|s| {
    //         let j = s.spawn(|| {
    //             let reader = BufReader::new(io::stdin());
    //             for b in reader.bytes() {
    //                 let byte = b.unwrap();
    //                 buffer.push(byte)
    //             }
    //         });
    //     });

    //     let setup = Setup{
    //         stream: None,
    //         measure: "n_high".to_string(),
    //         threads: 1,
    //         batchsize: 1,
    //         writer: BufWriter::new(Box::new(io::stdout())),
    //         distances_channel: bounded(100),
    //         wg_dist: WaitGroup::new(),
    //         efras: vec![load_fasta(LOAD).unwrap()],
    //         ns: vec![2],
    //         consensus: None,
    //     };

    //     load(setup);

    //     // print!("{}", String::from_utf8_lossy(&b));
    // }
}