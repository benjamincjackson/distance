use clap::{value_parser, crate_version, Arg, ArgMatches, Command};
use crossbeam_channel::{bounded, Receiver, Sender};
use crossbeam_utils::sync::WaitGroup;
use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::BufWriter;
use std::num::ParseIntError;
use std::sync::Arc;
use std::thread;
use std::thread::JoinHandle;
use thiserror::Error;

mod measures;
use crate::measures::*;
mod fastaio;
use crate::fastaio::*;

type Result<T> = std::result::Result<T, DistanceError>;

#[derive(Debug, Error)]
pub enum DistanceError {
    #[error(transparent)]
    IOError(#[from] io::Error),
    #[error(transparent)]
    ParseIntError(#[from] ParseIntError),
    #[error(transparent)]
    MatchesError(#[from] clap::parser::MatchesError),
    #[error(transparent)]
    ChanSendFastaErr(#[from] crossbeam_channel::SendError<fastaio::Records>),
    #[error(transparent)]
    ChanSendDistErr(#[from] crossbeam_channel::SendError<Distances>),
    #[error(transparent)]
    ChanSendPairErr(#[from] crossbeam_channel::SendError<Pairs>),
    #[error(transparent)]
    ChanRecvErr(#[from] crossbeam_channel::RecvError),
    #[error("")]
    Message(String),
}

// A struct for passing the location of one pairwise comparison down a channel (between threads)
#[derive(Clone, Debug, PartialEq)]
struct Pair {
    seq1_idx: usize,
    seq2_idx: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Pairs {
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
pub struct Distances {
    distances: Vec<Distance>,
    idx: usize,
}

pub fn get_cli_arguments() -> ArgMatches {
    // Define the command-line interface
    Command::new("distance")
        .version(crate_version!())
        .about("Calculate genetic distances within/between fasta-format alignments of DNA sequences")
        .override_usage(
            "All sequences across all input files must be the same length.\n       \
             \n       \
             distance alignment.fasta\n       \
             cat alignment.fasta | distance\n       \
             distance alignment.fasta -o distances.tsv\n       \
             distance -t 8 -m jc69 alignment.fasta -o jc69.tsv\n       \
             distance alignment1.fasta alignment2.fasta > distances2.tsv\n       \
             distance -i smallAlignment.fasta -s bigAlignment.fasta -o distances3.tsv\n       \
             cat bigAlignment.fasta | distance smallAlignment.fasta -s - > distances3.tsv\n       \
             "
       )
        .arg(Arg::new("input")
            .short('i')
            .long("input")
            .help("One or two input alignment files in fasta format. Loaded into memory. This flag can be omitted and the files passed as positional arguments")
            .num_args(0..=2)
            .required(false))
        .arg(Arg::new("input_pos_1")
            .index(1)
            .required(false)
            .hide(true))
        .arg(Arg::new("input_pos_2")
            .index(2)
            .required(false)
            .hide(true))
        .arg(Arg::new("stream")
            .short('s')
            .long("stream")
            .required(false)
            .help("One input alignment file in fasta format. Streamed from disk (or stdin using \"-s -\"). Requires exactly one file also be loaded"))
        .arg(Arg::new("measure")
            .short('m')
            .long("measure")
            .default_value("raw")
            .value_parser(["n", "n_high", "raw", "jc69", "k80", "tn93"])
            .help("Which distance measure to use"))
        .arg(Arg::new("output")
            .short('o')
            .long("output")
            .help("Output file in tab-separated-value format. Omit this option to print to stdout"))
        .arg(Arg::new("threads")
            .short('t')
            .long("threads")
            .value_parser(value_parser!(usize))
            .help("How many threads to spin up for pairwise comparisons. Omitting this option spins up the number of available CPUs"))
        .arg(Arg::new("batchsize")
            .long("batchsize")
            .short('b')
            .default_value("1")
            .value_parser(value_parser!(usize))
            .help("Try setting this >(>) 1 to tune the workload per thread"))
        .arg(Arg::new("licenses")
            .long("licenses")
            .short('l')
            .num_args(0)
            .help("Print licence information and exit"))
        .get_matches()
}

pub struct Setup {
    loaded_fastas: Vec<Vec<EncodedFastaRecord>>,
    streamed_fasta: Option<Box<dyn io::Read + Send>>,
    writer: BufWriter<Box<dyn io::Write + Send>>,
    measure: String,
    n_threads: usize,
    batchsize: usize,
    distances_channel: (Sender<Distances>, Receiver<Distances>),
    distances_waitgroup: WaitGroup,
    sample_sizes: Vec<usize>,
    consensus: Option<EncodedFastaRecord>,
}
impl Setup {
    fn new() -> Setup {
        Setup {
            loaded_fastas: vec![vec![]],
            streamed_fasta: None,
            writer: BufWriter::new(Box::new(io::stdout())),
            measure: String::new(),
            n_threads: 1,
            batchsize: 1,
            distances_channel: bounded(100),
            distances_waitgroup: WaitGroup::new(),
            sample_sizes: Vec::new(),
            consensus: None,
        }
    }
}

pub fn set_up(m: &ArgMatches) -> Result<Setup> {
    let mut setup = Setup::new();

    // One or two input fasta file names (or stdin)

    // Inputs from positional arguments
    let mut pos_inputs: Vec<String> = Vec::new();
    if let Some(ip1) = m.get_one::<String>("input_pos_1") {
        pos_inputs.push(ip1.to_string());
    }
    if let Some(ip2) = m.get_one::<String>("input_pos_2") {
        pos_inputs.push(ip2.to_string());
    }

    // Inputs from -i/--input flag
    let mut flag_inputs: Vec<String> = Vec::new();
    if let Some(fi) = m.get_many::<String>("input") {
        flag_inputs = fi.map(|s| s.into()).collect();
    }

    if !pos_inputs.is_empty() && !flag_inputs.is_empty() {
        return Err(DistanceError::Message("For loading input files, don't use both positional arguments and the -i/--input flag".to_string()));
    }

    let consolidated_inputs: Vec<String> = [flag_inputs, pos_inputs].concat();
    let mut inputs: Vec<Box<dyn std::io::Read>> = Vec::new();

    if consolidated_inputs.is_empty() {
        inputs.push(Box::new(io::stdin()))
    }
    for path in &consolidated_inputs {
        inputs.push(Box::new(File::open(path)?))
    }

    if let Some(stream) = m.get_one::<String>("stream") {
        if consolidated_inputs.len() != 1 {
            return Err(DistanceError::Message("If you stream one file, you must also provide exactly one other file to be loaded".to_string()));
        }
        match stream.as_str() {
            "-" => {
                setup.streamed_fasta = Some(Box::new(io::stdin()))
            },
            _ => {
                setup.streamed_fasta = Some(Box::new(File::open(stream)?))
            }
        }
    }

    // Which distance measure to use
    setup.measure = m.get_one::<String>("measure").unwrap().into();

    // batch size - to tune the workload per message so that threads aren't fighting over pairs_receiver's lock as often
    setup.batchsize = *m.get_one::<usize>("batchsize").unwrap();

    // The aligned sequence data
    setup.loaded_fastas = load_fastas(inputs)?;

    // Need to do some extra work depending on which distance measure is used.
    match setup.measure.as_str() {
        // For the fast snp-distance, need to calculate the consensus then get the differences from
        // it for each record (in each file)
        "n" => {
            let consensus = consensus(&setup.loaded_fastas);
            for i in 0..setup.loaded_fastas.len() {
                for j in 0..setup.loaded_fastas[i].len() {
                    setup.loaded_fastas[i][j].get_differences(&consensus);
                }
            }
            setup.consensus = Some(consensus);
        }
        // For Tamura and Nei (1993), need to calculate the base content of each record.
        "tn93" => {
            for i in 0..setup.loaded_fastas.len() {
                for j in 0..setup.loaded_fastas[i].len() {
                    setup.loaded_fastas[i][j].count_bases();
                }
            }
        }
        _ => (),
    }

    // The sample sizes (for generating the pairs)
    for records_vec in &setup.loaded_fastas {
        setup.sample_sizes.push(records_vec.len());
    }

    if let Some(output) = m.get_one::<String>("output") {
        setup.writer = BufWriter::new(Box::new(File::create(output)?))
    }

    let threads = m.get_one::<usize>("threads");
    match threads {
        Some(n) => {
            if *n < 1 {
                setup.n_threads = 1;
            } else {
                setup.n_threads = *n;
            }
        }
        None => {
            setup.n_threads = num_cpus::get();
        }
    }

    Ok(setup)
}

pub fn stream(setup: Setup) -> Result<()> {
    let loaded_fastas_arc = Arc::new(setup.loaded_fastas);

    let (records_sender, records_receiver) = bounded(100);
    let (distances_sender, distances_receiver) = setup.distances_channel;

    let measure = setup.measure.clone();

    // We spin up a thread to write the output as it arrives down the distance channel.
    let write_joinhandle = thread::spawn({
        move || {
            let result = gather_write(setup.writer, distances_receiver);
            match result {
                Err(e) => Err(e),
                Ok(_) => Ok(()),
            }
        }
    });

    let stream_joinhandle = thread::spawn({
        let measure = measure.clone();
        let loaded_fastas_arc = loaded_fastas_arc.clone();
        let s = setup.streamed_fasta.expect("Stream not specified");
        move || {
            let r = stream_fasta(
                s,
                &loaded_fastas_arc,
                &measure,
                setup.consensus,
                setup.batchsize,
                records_sender,
            );
            match r {
                Err(e) => Err(e),
                Ok(_) => Ok(()),
            }
        }
    });

    // A vector of receiver/sender tuples, cloned to share between each thread in threads.
    let mut workers = Vec::new();
    for _ in 0..setup.n_threads {
        workers.push((records_receiver.clone(), distances_sender.clone()))
    }

    let f = get_distance_function(&measure);

    // Spin up the threads that do the distance-calculating
    for worker in workers {
        let wg_dist = setup.distances_waitgroup.clone();
        let loaded_fastas_arc = loaded_fastas_arc.clone();
        thread::spawn(move || {
            let mut distances: Vec<Distance> = vec![];
            for message in worker.0.iter() {
                for record_2 in message.records.iter() {
                    for record_1 in &loaded_fastas_arc[0] {
                        let d = f(record_1, record_2);
                        // add the ids and the distance to this batch's temporary vector
                        distances.push(Distance {
                            id1: record_1.id.clone(),
                            id2: record_2.id.clone(),
                            dist: d,
                        });
                    }
                }
                // send this batch of distances to the writer
                worker
                    .1
                    .send(Distances {
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

    stream_joinhandle
        .join()
        .unwrap()?;

    // When all the distances have been calculated, we can drop the sending end of the distance channel
    setup.distances_waitgroup.wait();
    drop(distances_sender);

    // Joins when all the pairwise comparisons have been written, and then we're done.
    write_joinhandle
        .join()
        .unwrap()?;

    Ok(())
}

pub fn load(setup: Setup) -> Result<()> {
    let loaded_fastas_arc = Arc::new(setup.loaded_fastas);

    let (pairs_sender, pairs_receiver) = bounded(100);
    let (distances_sender, distances_receiver) = setup.distances_channel;

    let sample_sizes = setup.sample_sizes.clone();
    let measure = setup.measure.clone();

    // We spin up a thread to write the output as it arrives down the distance channel.
    let write_joinhandle = thread::spawn({
        move || {
            let r = gather_write(setup.writer, distances_receiver);
            match r {
                Err(e) => Err(e),
                Ok(_) => Ok(()),
            }
        }
    });

    // If there is one input file, generate all pairwise comparisons within the alignment,
    // else there are two input files, so generate all pairwise comparisons between the alignments.
    // We spin up a thread do to this and move on to the next part of the program
    let make_pairs_joinhandle: JoinHandle<std::prelude::v1::Result<(), DistanceError>> = if sample_sizes.len() == 1 {
    thread::spawn({
        move || {
            let r = generate_pairs_square(sample_sizes[0], setup.batchsize, pairs_sender);
            match r {
                Err(e) => Err(e),
                Ok(_) => Ok(()),
            }
        }
    })
    } else {
        thread::spawn({
            move || {
                let r = generate_pairs_rectangle(sample_sizes[0], sample_sizes[1], setup.batchsize, pairs_sender);
                match r {
                    Err(e) => Err(e),
                    Ok(_) => Ok(()),
                }
            }
        })
    };

    // A vector of receiver/sender tuples, cloned to share between each thread in threads.
    let mut workers = Vec::new();
    for _ in 0..setup.n_threads {
        workers.push((pairs_receiver.clone(), distances_sender.clone()))
    }

    // Which distance function to use
    let f = get_distance_function(&measure);

    // Spin up the threads that do the distance-calculating
    for worker in workers {
        let wg_dist = setup.distances_waitgroup.clone();
        let loaded_fastas_arc = loaded_fastas_arc.clone();
        thread::spawn(move || {
            let mut distances: Vec<Distance> = vec![];
            // for each batch
            for message in worker.0.iter() {
                // for each pair in this batch
                for pair in message.pairs {
                    // calculate the distance
                    let record_1 = &loaded_fastas_arc[0][pair.seq1_idx];
                    let record_2 = &loaded_fastas_arc[loaded_fastas_arc.len() - 1][pair.seq2_idx];
                    let d = f(record_1, record_2);
                    // add the ids and the distance to this batch's temporary vector
                    distances.push(Distance {
                        id1: record_1.id.clone(),
                        id2: record_2.id.clone(),
                        dist: d,
                    });
                }

                // send this batch of distances to the writer
                worker
                    .1
                    .send(Distances {
                        distances: distances.clone(),
                        idx: message.idx,
                    }).unwrap();

                // clear the vector ready for the next batch
                distances.clear();
            }

            // when the pair channel is empty (and all distances are calculated) we can drop the cloned distance waitgroup (for this thread)
            drop(wg_dist);
        });
    }

    make_pairs_joinhandle
        .join()
        .unwrap()?;

    // When all the distances have been calculated, we can drop the sending end of the distance channel
    setup.distances_waitgroup.wait();
    drop(distances_sender);

    // Joins when all the pairwise comparisons have been written, and then we're done.
    write_joinhandle
        .join()
        .unwrap()?;

    Ok(())
}

// Return the correct distance function given the CLI input
fn get_distance_function(s: &str) -> fn(&EncodedFastaRecord, &EncodedFastaRecord) -> FloatInt {
    match s {
        "raw" => raw,
        "n" => snp_consensus,
        "n_high" => snp,
        "jc69" => jc69,
        "k80" => k80,
        "tn93" => tn93,
        // should never get this far because the options are defined in the cli:
        _ => panic!("Unknown distance measure"),
    }
}

pub fn run(setup: Setup) -> Result<()> {
    if setup.streamed_fasta.is_some() {
        stream(setup)?;
    } else {
        load(setup)?;
    }

    Ok(())
}

// Given the sample size of a single alignment, generate all possible pairwise comparisons
// within it, and pass them down a channel.
fn generate_pairs_square(n: usize, size: usize, sender: Sender<Pairs>) -> Result<()> {
    // this counter is used to send the correct-sized batch
    let mut size_counter: usize = 0;

    // this counter is sent down the channel in order to later retain input order in the output
    let mut idx_counter: usize = 0;

    let mut pair_vec: Vec<Pair> = Vec::new();

    for i in 0..n - 1 {
        for j in i + 1..n {
            pair_vec.push(Pair {
                seq1_idx: i,
                seq2_idx: j,
            });

            size_counter += 1;

            // when we reach the batch size, we send this batch of pairs
            if size_counter == size {
                sender
                    .send(Pairs {
                        pairs: pair_vec.clone(),
                        idx: idx_counter,
                    })?;

                size_counter = 0;
                idx_counter += 1;
                pair_vec.clear();
            }
        }
    }

    // send the last batch
    if !pair_vec.is_empty() {
        sender
            .send(Pairs {
                pairs: pair_vec.clone(),
                idx: idx_counter,
            })?;
    }

    drop(sender);

    Ok(())
}

// Given the sample size of two alignments, generate all possible pairwise comparisons
// between them, and pass them down a channel.
fn generate_pairs_rectangle(n1: usize, n2: usize, size: usize, sender: Sender<Pairs>) -> Result<()> {
    // this counter is used to send the correct-sized batch
    let mut size_counter: usize = 0;

    // this counter is sent down the channel in order to later retain input order in the output
    let mut idx_counter: usize = 0;

    let mut pair_vec: Vec<Pair> = Vec::new();

    for i in 0..n1 {
        for j in 0..n2 {
            pair_vec.push(Pair {
                seq1_idx: i,
                seq2_idx: j,
            });

            size_counter += 1;

            // when we reach the batch size, we send this batch of pairs
            if size_counter == size {
                sender
                    .send(Pairs {
                        pairs: pair_vec.clone(),
                        idx: idx_counter,
                    })?;

                size_counter = 0;
                idx_counter += 1;
                pair_vec.clear();
            }
        }
    }

    // send the last batch
    if !pair_vec.is_empty() {
        sender
            .send(Pairs {
                pairs: pair_vec.clone(),
                idx: idx_counter,
            })?;
    }

    drop(sender);

    Ok(())
}

fn handle_broken_pipe(r: std::result::Result<(), io::Error>) -> Result<()> {
    if let Err(e) = r {
        if e.kind() == io::ErrorKind::BrokenPipe {
            std::process::exit(0);
        } else {
            return Err(DistanceError::IOError(e))
        }
    }

    Ok(())
}

// Write the distances as they arrive. Uses a hashmap whose keys are indices to write the results in the
// order they are produced by generate_pairs_*()
fn gather_write<T: io::Write>(mut writer: T, rx: Receiver<Distances>) -> Result<()> {
    let r = writeln!(writer, "sequence1\tsequence2\tdistance");
    handle_broken_pipe(r)?;

    let mut m: HashMap<usize, Distances> = HashMap::new();

    let mut counter: usize = 0;

    for r in rx.iter() {
        m.insert(r.idx, r);
        while m.contains_key(&counter) {
            let rv = m.remove(&counter).unwrap();
            for result in rv.distances {
                match result.dist {
                    FloatInt::Int(d) => {
                        let r = writeln!(writer, "{}\t{}\t{}", &result.id1, &result.id2, d);
                        handle_broken_pipe(r)?;
                    }
                    FloatInt::Float(d) => {
                        let r = writeln!(writer, "{}\t{}\t{:.12}", &result.id1, &result.id2, d);
                        handle_broken_pipe(r)?;
                    }
                }
            }
            counter += 1;
        }
    }

    let r = writer.flush();
    handle_broken_pipe(r)?;
    
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_generate_pairs_square() -> Result<()> {
        let n: usize = 4;
        let batch_size = 1;
        let (sx, rx) = bounded(50);

        let jh = thread::spawn({
            let sx = sx.clone();
            move || {
                let r = generate_pairs_square(n, batch_size, sx);
                match r {
                    Err(e) => Err(e),
                    Ok(_) => Ok(()),
                }
            }
        });

        assert_eq!(
            rx.recv().unwrap(),
            Pairs {
                pairs: vec![Pair {
                    seq1_idx: 0,
                    seq2_idx: 1
                }],
                idx: 0
            }
        );
        assert_eq!(
            rx.recv().unwrap(),
            Pairs {
                pairs: vec![Pair {
                    seq1_idx: 0,
                    seq2_idx: 2
                }],
                idx: 1
            }
        );
        assert_eq!(
            rx.recv().unwrap(),
            Pairs {
                pairs: vec![Pair {
                    seq1_idx: 0,
                    seq2_idx: 3
                }],
                idx: 2
            }
        );
        assert_eq!(
            rx.recv().unwrap(),
            Pairs {
                pairs: vec![Pair {
                    seq1_idx: 1,
                    seq2_idx: 2
                }],
                idx: 3
            }
        );
        assert_eq!(
            rx.recv().unwrap(),
            Pairs {
                pairs: vec![Pair {
                    seq1_idx: 1,
                    seq2_idx: 3
                }],
                idx: 4
            }
        );
        assert_eq!(
            rx.recv().unwrap(),
            Pairs {
                pairs: vec![Pair {
                    seq1_idx: 2,
                    seq2_idx: 3
                }],
                idx: 5
            }
        );
        assert!(rx.is_empty());

        jh
            .join()
            .unwrap()?;

        let batch_size_2 = 4;

        let jh = thread::spawn({
            let sx = sx.clone();
            move || {
                let r = generate_pairs_square(n, batch_size_2, sx);
                match r {
                    Err(e) => Err(e),
                    Ok(_) => Ok(()),
                }
            }
        });

        assert_eq!(
            rx.recv().unwrap(),
            Pairs {
                pairs: vec![
                    Pair {
                        seq1_idx: 0,
                        seq2_idx: 1
                    },
                    Pair {
                        seq1_idx: 0,
                        seq2_idx: 2
                    },
                    Pair {
                        seq1_idx: 0,
                        seq2_idx: 3
                    },
                    Pair {
                        seq1_idx: 1,
                        seq2_idx: 2
                    }
                ],
                idx: 0
            }
        );
        assert_eq!(
            rx.recv().unwrap(),
            Pairs {
                pairs: vec![
                    Pair {
                        seq1_idx: 1,
                        seq2_idx: 3
                    },
                    Pair {
                        seq1_idx: 2,
                        seq2_idx: 3
                    }
                ],
                idx: 1
            }
        );
        assert!(rx.is_empty());

    jh
        .join()
        .unwrap()?;

    Ok(())
    }

    #[test]
    fn test_generate_pairs_rectangle() -> Result<()> {
        let n1: usize = 2;
        let n2: usize = 2;
        let batch_size = 1;
        let (sx, rx) = bounded(50);

        let jh = thread::spawn({
            let sx = sx.clone();
            move || {
                let r = generate_pairs_rectangle(n1, n2, batch_size, sx);
                match r {
                    Err(e) => Err(e),
                    Ok(_) => Ok(()),
                }
            }
        });

        jh
            .join()
            .unwrap()?;

        assert_eq!(
            rx.recv().unwrap(),
            Pairs {
                pairs: vec![Pair {
                    seq1_idx: 0,
                    seq2_idx: 0
                }],
                idx: 0
            }
        );
        assert_eq!(
            rx.recv().unwrap(),
            Pairs {
                pairs: vec![Pair {
                    seq1_idx: 0,
                    seq2_idx: 1
                }],
                idx: 1
            }
        );
        assert_eq!(
            rx.recv().unwrap(),
            Pairs {
                pairs: vec![Pair {
                    seq1_idx: 1,
                    seq2_idx: 0
                }],
                idx: 2
            }
        );
        assert_eq!(
            rx.recv().unwrap(),
            Pairs {
                pairs: vec![Pair {
                    seq1_idx: 1,
                    seq2_idx: 1
                }],
                idx: 3
            }
        );
        assert!(rx.is_empty());

        let batch_size_2 = 4;

        let jh = thread::spawn({
            let sx = sx.clone();
            move || {
                let r = generate_pairs_rectangle(n1, n2, batch_size_2, sx);
                match r {
                    Err(e) => Err(e),
                    Ok(_) => Ok(()),
                }
            }
        });

        assert_eq!(
            rx.recv().unwrap(),
            Pairs {
                pairs: vec![
                    Pair {
                        seq1_idx: 0,
                        seq2_idx: 0
                    },
                    Pair {
                        seq1_idx: 0,
                        seq2_idx: 1
                    },
                    Pair {
                        seq1_idx: 1,
                        seq2_idx: 0
                    },
                    Pair {
                        seq1_idx: 1,
                        seq2_idx: 1
                    }
                ],
                idx: 0
            }
        );
        assert!(rx.is_empty());

    
    jh
        .join()
        .unwrap()?;

    Ok(())
    }

    const FASTA_1: &[u8] = b">seq1
ATGATG
>seq2
ATGATC
";

    const FASTA_2: &[u8] = b">seqA
ATGATG
";

    use std::io::{pipe, Read};

    #[test]
    fn test_integration_1() {
        // Test with m = "n", batch size 1
        let (mut reader, writer) = pipe().unwrap();
        let mut setup = Setup::new();
        setup.streamed_fasta = None;
        setup.measure = "n".to_string();
        setup.n_threads = 1;
        setup.batchsize = 1;
        setup.writer = BufWriter::new(Box::new(writer));
        setup.loaded_fastas = load_fastas(vec![FASTA_1]).unwrap();
        setup.sample_sizes = setup.loaded_fastas.iter().map(|r| r.len()).collect();
        let c = consensus(&setup.loaded_fastas);
        for i in 0..setup.loaded_fastas.len() {
            for j in 0..setup.loaded_fastas[i].len() {
                setup.loaded_fastas[i][j].get_differences(&c);
            }
        }
        setup.consensus = Some(c);

        let expected = "sequence1\tsequence2\tdistance
seq1\tseq2\t1
".to_string();
        let result = run(setup);
        assert!(result.is_ok());
        let mut output = String::new();
        let _ = reader.read_to_string(&mut output).unwrap();
        assert_eq!(expected, output);

        // Test with batch size 2
        let (mut reader, writer) = pipe().unwrap();
        let mut setup = Setup::new();
        setup.streamed_fasta = None;
        setup.measure = "n".to_string();
        setup.n_threads = 1;
        setup.batchsize = 2;
        setup.writer = BufWriter::new(Box::new(writer));
        setup.loaded_fastas = load_fastas(vec![FASTA_1]).unwrap();
        setup.sample_sizes = setup.loaded_fastas.iter().map(|r| r.len()).collect();
        let c = consensus(&setup.loaded_fastas);
        for i in 0..setup.loaded_fastas.len() {
            for j in 0..setup.loaded_fastas[i].len() {
                setup.loaded_fastas[i][j].get_differences(&c);
            }
        }
        setup.consensus = Some(c);

        let expected = "sequence1\tsequence2\tdistance
seq1\tseq2\t1
".to_string();
        let result = run(setup);
        assert!(result.is_ok());
        let mut output = String::new();
        let _ = reader.read_to_string(&mut output).unwrap();
        assert_eq!(expected, output);

        // Test with m = "n", batch size 2, threads = 2
        let (mut reader, writer) = pipe().unwrap();
        let mut setup = Setup::new();
        setup.streamed_fasta = None;
        setup.measure = "n".to_string();
        setup.n_threads = 2;
        setup.batchsize = 2;
        setup.writer = BufWriter::new(Box::new(writer));
        setup.loaded_fastas = load_fastas(vec![FASTA_1]).unwrap();
        setup.sample_sizes = setup.loaded_fastas.iter().map(|r| r.len()).collect();
        let c = consensus(&setup.loaded_fastas);
        for i in 0..setup.loaded_fastas.len() {
            for j in 0..setup.loaded_fastas[i].len() {
                setup.loaded_fastas[i][j].get_differences(&c);
            }
        }
        setup.consensus = Some(c);

        let expected = "sequence1\tsequence2\tdistance
seq1\tseq2\t1
".to_string();
        let result = run(setup);
        assert!(result.is_ok());
        let mut output = String::new();
        let _ = reader.read_to_string(&mut output).unwrap();
        assert_eq!(expected, output);

    }

    #[test]
    fn test_integration_2() {
        // Test with m = "n_high", stream
        let (mut reader, writer) = pipe().unwrap();
        let mut setup = Setup::new();
        setup.streamed_fasta = Some(Box::new(FASTA_2));
        setup.measure = "n_high".to_string();
        setup.n_threads = 1;
        setup.batchsize = 1;
        setup.writer = BufWriter::new(Box::new(writer));
        setup.loaded_fastas = load_fastas(vec![FASTA_1]).unwrap();
        setup.sample_sizes = setup.loaded_fastas.iter().map(|r| r.len()).collect();

        let expected = "sequence1\tsequence2\tdistance
seq1\tseqA\t0
seq2\tseqA\t1
".to_string();
        let result = run(setup);
        assert!(result.is_ok());
        let mut output = String::new();
        let _ = reader.read_to_string(&mut output).unwrap();
        assert_eq!(expected, output);

        // Test with m = "n_high", stream, batch size 2
        let (mut reader, writer) = pipe().unwrap();
        let mut setup = Setup::new();
        setup.streamed_fasta = Some(Box::new(FASTA_2));
        setup.measure = "n_high".to_string();
        setup.n_threads = 1;
        setup.batchsize = 2;
        setup.writer = BufWriter::new(Box::new(writer));
        setup.loaded_fastas = load_fastas(vec![FASTA_1]).unwrap();
        setup.sample_sizes = setup.loaded_fastas.iter().map(|r| r.len()).collect();

        let expected = "sequence1\tsequence2\tdistance
seq1\tseqA\t0
seq2\tseqA\t1
".to_string();
        let result = run(setup);
        assert!(result.is_ok());
        let mut output = String::new();
        let _ = reader.read_to_string(&mut output).unwrap();
        assert_eq!(expected, output);

        // Test with m = "n_high", stream, batch size 2, threads = 2
        let (mut reader, writer) = pipe().unwrap();
        let mut setup = Setup::new();
        setup.streamed_fasta = Some(Box::new(FASTA_2));
        setup.measure = "n_high".to_string();
        setup.n_threads = 2;
        setup.batchsize = 2;
        setup.writer = BufWriter::new(Box::new(writer));
        setup.loaded_fastas = load_fastas(vec![FASTA_1]).unwrap();
        setup.sample_sizes = setup.loaded_fastas.iter().map(|r| r.len()).collect();

        let expected = "sequence1\tsequence2\tdistance
seq1\tseqA\t0
seq2\tseqA\t1
".to_string();
        let result = run(setup);
        assert!(result.is_ok());
        let mut output = String::new();
        let _ = reader.read_to_string(&mut output).unwrap();
        assert_eq!(expected, output);
    }

    #[test]
    fn test_integration_3() {
        // Test with m = "n_high", two loaded inputs
        let (mut reader, writer) = pipe().unwrap();
        let mut setup = Setup::new();
        setup.streamed_fasta = None;
        setup.measure = "n_high".to_string();
        setup.n_threads = 1;
        setup.batchsize = 1;
        setup.writer = BufWriter::new(Box::new(writer));
        setup.loaded_fastas = load_fastas(vec![FASTA_1, FASTA_2]).unwrap();
        setup.sample_sizes = setup.loaded_fastas.iter().map(|r| r.len()).collect();

        let expected = "sequence1\tsequence2\tdistance
seq1\tseqA\t0
seq2\tseqA\t1
".to_string();
        let result = run(setup);
        assert!(result.is_ok());
        let mut output = String::new();
        let _ = reader.read_to_string(&mut output).unwrap();
        assert_eq!(expected, output);

        // Test with m = "n_high", two loaded inputs, batchsize 2
        let (mut reader, writer) = pipe().unwrap();
        let mut setup = Setup::new();
        setup.streamed_fasta = None;
        setup.measure = "n_high".to_string();
        setup.n_threads = 1;
        setup.batchsize = 2;
        setup.writer = BufWriter::new(Box::new(writer));
        setup.loaded_fastas = load_fastas(vec![FASTA_1, FASTA_2]).unwrap();
        setup.sample_sizes = setup.loaded_fastas.iter().map(|r| r.len()).collect();

        let expected = "sequence1\tsequence2\tdistance
seq1\tseqA\t0
seq2\tseqA\t1
".to_string();
        let result = run(setup);
        assert!(result.is_ok());
        let mut output = String::new();
        let _ = reader.read_to_string(&mut output).unwrap();
        assert_eq!(expected, output);

        // Test with m = "n_high", two loaded inputs, batchsize 2, threads 2
        let (mut reader, writer) = pipe().unwrap();
        let mut setup = Setup::new();
        setup.streamed_fasta = None;
        setup.measure = "n_high".to_string();
        setup.n_threads = 2;
        setup.batchsize = 2;
        setup.writer = BufWriter::new(Box::new(writer));
        setup.loaded_fastas = load_fastas(vec![FASTA_1, FASTA_2]).unwrap();
        setup.sample_sizes = setup.loaded_fastas.iter().map(|r| r.len()).collect();

        let expected = "sequence1\tsequence2\tdistance
seq1\tseqA\t0
seq2\tseqA\t1
".to_string();
        let result = run(setup);
        assert!(result.is_ok());
        let mut output = String::new();
        let _ = reader.read_to_string(&mut output).unwrap();
        assert_eq!(expected, output);

        // Test with m = "n_high", two loaded inputs, reverse order
        let (mut reader, writer) = pipe().unwrap();
        let mut setup = Setup::new();
        setup.streamed_fasta = None;
        setup.measure = "n_high".to_string();
        setup.n_threads = 1;
        setup.batchsize = 1;
        setup.writer = BufWriter::new(Box::new(writer));
        setup.loaded_fastas = load_fastas(vec![FASTA_2, FASTA_1]).unwrap();
        setup.sample_sizes = setup.loaded_fastas.iter().map(|r| r.len()).collect();

        let expected = "sequence1\tsequence2\tdistance
seqA\tseq1\t0
seqA\tseq2\t1
".to_string();
        let result = run(setup);
        assert!(result.is_ok());
        let mut output = String::new();
        let _ = reader.read_to_string(&mut output).unwrap();
        assert_eq!(expected, output);
    }
}