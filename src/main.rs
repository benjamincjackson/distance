use std::io;
use std::io::Write;
use std::io::BufWriter;
use std::fs::File;
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
    output: String,
}

fn main() -> io::Result<()> {
	let args = Args::parse();

//	 let (w, n) = align_dims(&args.input).unwrap();
	
//	 println!("{} {}", w, n);

//	let ba = populate_array(&args.input);

   let efra = populate_struct_array(&args.input)?;

   let mut d: f64 = 0.0;

   let f = File::create(&args.output)?;
   let mut buf = BufWriter::new(f);
   writeln!(buf, "sequence1,sequence2,distance")?;

   for i in 0..efra.len()-1 {
	   for j in i+1..efra.len() {
            d = tn93(&efra[i], &efra[j]);
            writeln!(buf, "{}\t{}\t{}", &efra[i].id, &efra[j].id, d)?;
	   }
   }

   buf.flush()?;

   Ok(())
}
