use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::Write;
use std::io::BufWriter;
use crossbeam_channel::bounded;
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

struct Result {
    id1: String,
    id2: String,
    dist: f64,
    idx: usize
}

fn main() -> io::Result<()> {
	let args = Args::parse();

//	 let (w, n) = align_dims(&args.input).unwrap();
	
//	 println!("{} {}", w, n);

//	let ba = populate_array(&args.input);

    let efra = Box::new(populate_struct_array(&args.input)).unwrap(); 
    // let efra = populate_struct_array(&args.input)?; 

    let (pair_sender, pair_receiver) = bounded(50);
    // let (distance_sender, distance_receiver) = bounded(50);

    thread::spawn( move || {
    
        let mut d: f64 = 0.0;
        let mut idx: usize = 0;

        for i in 0..efra.len()-1 {
            for j in i+1..efra.len() {
                d = tn93(&efra[i], &efra[j]);
                let r = Result{id1: efra[i].id.clone(), id2: efra[j].id.clone(), dist: d.clone(), idx: idx.clone()};
                pair_sender.send(r).unwrap();
                idx += 1;
            }
        }
   });


   let f = File::create(&args.output)?;
   let mut buf = BufWriter::new(f);
   writeln!(buf, "sequence1\tsequence2\tdistance")?;

   for r in pair_receiver {
        writeln!(buf, "{}\t{}\t{}", &r.id1, &r.id2, r.dist)?;
   }

   buf.flush()?;

   Ok(())
}

// fn gather_write(filename: &str, rx: Receiver<Result>) -> io::Result<()> {

//     let f = File::create(filename)?;
//     let mut buf = BufWriter::new(f);
//     writeln!(buf, "sequence1\tsequence2\tdistance")?;

//     let mut m: HashMap<usize,Result> = HashMap::new();

//     let mut counter: usize = 0;

//     for r in rx {

//         m.insert(r.idx, r);

//         if m.contains_key(&counter) {
//             let r = m.remove(&counter).unwrap();
//             writeln!(buf, "{}\t{}\t{}", &r.id1, &r.id2, r.dist)?;
//             counter += 1;
//         }

//     }

//     while m.len() > 0 {
//         let r = m.remove(&counter).unwrap();
//         writeln!(buf, "{}\t{}\t{}", &r.id1, &r.id2, r.dist)?;
//         counter += 1;
//     }

//     Ok(())
// }
