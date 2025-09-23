use clap::value_parser;
use clap::{crate_version, Arg, ArgMatches, Command};
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

mod measures;
use crate::measures::*;
mod fastaio;
use crate::fastaio::*;

type Result<T> = std::result::Result<T, DistanceError>;

#[derive(Debug)]
pub enum DistanceError {
    IOError(io::Error),
    ParseIntError(ParseIntError),
    ChanSendFastaErr(crossbeam_channel::SendError<fastaio::Records>),
    ChanSendDistErr(crossbeam_channel::SendError<Distances>),
    ChanSendPairErr(crossbeam_channel::SendError<Pairs>),
    ChanRecvErr(crossbeam_channel::RecvError),
    Message(String),
}
impl From<io::Error> for DistanceError {
    fn from(err: io::Error) -> Self {
        DistanceError::IOError(err)
    }
}
impl From<crossbeam_channel::SendError<fastaio::Records>> for DistanceError {
    fn from(err: crossbeam_channel::SendError<fastaio::Records>) -> Self {
        DistanceError::ChanSendFastaErr(err)
    }
}
impl From<crossbeam_channel::SendError<Distances>> for DistanceError {
    fn from(err: crossbeam_channel::SendError<Distances>) -> Self {
        DistanceError::ChanSendDistErr(err)
    }
}
impl From<crossbeam_channel::SendError<Pairs>> for DistanceError {
    fn from(err: crossbeam_channel::SendError<Pairs>) -> Self {
        DistanceError::ChanSendPairErr(err)
    }
}
impl From<ParseIntError> for DistanceError {
    fn from(err: ParseIntError) -> Self {
        DistanceError::ParseIntError(err)
    }
}

fn err_message_stream_input_count() -> DistanceError {
    DistanceError::Message("If you stream one file, you must also provide exactly one other file to be loaded".to_string())
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
        .override_usage(
            "\n       \
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
            .default_value("1")
            .value_parser(value_parser!(usize))
            .help("How many threads to spin up for pairwise comparisons"))
        .arg(Arg::new("batchsize")
            .long("batchsize")
            .short('b')
            .default_value("1")
            .value_parser(value_parser!(usize))
            .help("Try setting this >(>) 1 if you are struggling to get a speedup when adding threads"))
        .arg(Arg::new("licenses")
            .long("licenses")
            .short('l')
            .num_args(0)
            .help("Print licence information and exit"))
        .get_matches()
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
        Setup {
            stream: None,
            measure: String::new(),
            threads: 1,
            batchsize: 1,
            writer: BufWriter::new(Box::new(io::stdout())),
            distances_channel: bounded(100),
            wg_dist: WaitGroup::new(),
            efras: vec![vec![]],
            ns: Vec::new(),
            consensus: None,
        }
    }
}

pub fn set_up(m: &ArgMatches) -> Result<Setup> {
    let mut setup = Setup::new();

    // One or two input fasta file names (or stdin)

    // Inputs from positional arguments
    let mut pos_inputs: Vec<String> = vec![];
    if let Some(ip1) = m.get_one::<String>("input_pos_1") {
        pos_inputs.push(ip1.to_string());
    }
    if let Some(ip2) = m.get_one::<String>("input_pos_2") {
        pos_inputs.push(ip2.to_string());
    }

    // Inputs from -i/--input flag
    let mut flag_inputs: Vec<String> = vec![];
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
    for s in &consolidated_inputs {
        inputs.push(Box::new(File::open(s)?))
    }

    if let Some(s) = m.get_one::<String>("stream") {
        if consolidated_inputs.len() != 1 {
            return Err(err_message_stream_input_count());
        }
        match s.as_str() {
            "-" => {
                setup.stream = Some(Box::new(io::stdin()))
            },
            _ => {
                setup.stream = Some(Box::new(File::open(s)?))
            }
        }
    }

    // Which distance measure to use
    setup.measure = m.get_one::<String>("measure").unwrap().into();

    // batch size - to tune the workload per message so that threads aren't fighting over pairs_receiver's lock as often
    setup.batchsize = *m.get_one::<usize>("batchsize").unwrap();

    // The aligned sequence data
    setup.efras = load_fastas(inputs)?;

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

    if let Some(output) = m.get_one::<String>("output") {
        setup.writer = BufWriter::new(Box::new(File::create(output)?))
    }

    // How many additional threads to use for calculating distances - need at least 1.
    setup.threads = *m.get_one::<usize>("threads").unwrap();

    if setup.threads < 1 {
        setup.threads = 1;
    }

    Ok(setup)
}

pub fn stream(setup: Setup) -> Result<()> {
    let arc = Arc::new(setup.efras);

    let (records_sender, records_receiver) = bounded(100);
    let (distances_sender, distances_receiver) = setup.distances_channel;

    let batchsize = setup.batchsize.to_owned();
    let measure = setup.measure.to_owned();

    // We spin up a thread to write the output as it arrives down the distance channel.
    let write = thread::spawn({
        move || {
            let result = gather_write(setup.writer, distances_receiver);
            match result {
                Err(e) => Err(e),
                Ok(_) => Ok(()),
            }
        }
    });

    let stream = thread::spawn({
        let measure = measure.clone();
        let arc = arc.clone();
        let s = setup.stream.expect("Stream not specified");
        move || {
            let r = stream_fasta(
                s,
                &arc,
                &measure,
                setup.consensus,
                batchsize,
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
    for _i in 0..setup.threads {
        workers.push((records_receiver.clone(), distances_sender.clone()))
    }

    let f = get_distance_function(measure.as_str());

    // Spin up the threads that do the distance-calculating
    for worker in workers {
        let wg_dist = setup.wg_dist.clone();
        let arc = arc.clone();
        thread::spawn(move || {
            let mut distances: Vec<Distance> = vec![];
            for message in worker.0.iter() {
                for target in message.records {
                    for query in &arc[0] {
                        let d = f(query, &target);
                        // add the ids and the distance to this batch's temporary vector
                        distances.push(Distance {
                            id1: target.id.clone(),
                            id2: query.id.clone(),
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

    stream
        .join()
        .unwrap()?;

    // When all the distances have been calculated, we can drop the sending end of the distance channel
    setup.wg_dist.wait();
    drop(distances_sender);

    // Joins when all the pairwise comparisons have been written, and then we're done.
    write
        .join()
        .unwrap()?;

    Ok(())
}

pub fn load(setup: Setup) -> Result<()> {
    let arc = Arc::new(setup.efras);

    let (pairs_sender, pairs_receiver) = bounded(100);
    let (distances_sender, distances_receiver) = setup.distances_channel;

    let ns = setup.ns.to_owned();
    let batchsize = setup.batchsize.to_owned();
    let measure = setup.measure.to_owned();

    // We spin up a thread to write the output as it arrives down the distance channel.
    let write = thread::spawn({
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

    let make_pairs: JoinHandle<std::prelude::v1::Result<(), DistanceError>> = if ns.len() == 1 {
    thread::spawn({
        move || {
            let r = generate_pairs_square(ns[0], batchsize, pairs_sender);
            match r {
                Err(e) => Err(e),
                Ok(_) => Ok(()),
            }
        }
    })
    } else {
        thread::spawn({
            move || {
                let r = generate_pairs_rect(ns[0], ns[1], batchsize, pairs_sender);
                match r {
                    Err(e) => Err(e),
                    Ok(_) => Ok(()),
                }
            }
        })
    };

    // A vector of receiver/sender tuples, cloned to share between each thread in threads.
    let mut workers = Vec::new();
    for _i in 0..setup.threads {
        workers.push((pairs_receiver.clone(), distances_sender.clone()))
    }

    // Which distance function to use
    let f = get_distance_function(measure.as_str());

    // Spin up the threads that do the distance-calculating
    for worker in workers {
        let wg_dist = setup.wg_dist.clone();
        let arc = arc.clone();
        thread::spawn(move || {
            let mut distances: Vec<Distance> = vec![];
            // for each batch
            for message in worker.0.iter() {
                // for each pair in this batch
                for pair in message.pairs {
                    // calculate the distance
                    let d = f(&arc[0][pair.seq1_idx], &arc[arc.len() - 1][pair.seq2_idx]);
                    // add the ids and the distance to this batch's temporary vector
                    distances.push(Distance {
                        id1: arc[0][pair.seq1_idx].id.clone(),
                        id2: arc[arc.len() - 1][pair.seq2_idx].id.clone(),
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

    make_pairs
        .join()
        .unwrap()?;

    // When all the distances have been calculated, we can drop the sending end of the distance channel
    setup.wg_dist.wait();
    drop(distances_sender);

    // Joins when all the pairwise comparisons have been written, and then we're done.
    write
        .join()
        .unwrap()?;

    Ok(())
}

// Return the correct distance function given the CLI input
fn get_distance_function(s: &str) -> fn(&EncodedFastaRecord, &EncodedFastaRecord) -> FloatInt {
    match s {
        "raw" => raw,
        "n" => snp2,
        "n_high" => snp,
        "jc69" => jc69,
        "k80" => k80,
        "tn93" => tn93,
        // should never get this far because the options are defined in the cli:
        _ => panic!("Unknown distance measure"),
    }
}

pub fn run(setup: Setup) -> Result<()> {
    if setup.stream.is_some() {
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

    let mut pair_vec: Vec<Pair> = vec![];

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
fn generate_pairs_rect(n1: usize, n2: usize, size: usize, sender: Sender<Pairs>) -> Result<()> {
    // this counter is used to send the correct-sized batch
    let mut size_counter: usize = 0;

    // this counter is sent down the channel in order to later retain input order in the output
    let mut idx_counter: usize = 0;

    let mut pair_vec: Vec<Pair> = vec![];

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
    fn test_generate_pairs_rect() -> Result<()> {
        let n1: usize = 2;
        let n2: usize = 2;
        let batch_size = 1;
        let (sx, rx) = bounded(50);

        let jh = thread::spawn({
            let sx = sx.clone();
            move || {
                let r = generate_pairs_rect(n1, n2, batch_size, sx);
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
                let r = generate_pairs_rect(n1, n2, batch_size_2, sx);
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
