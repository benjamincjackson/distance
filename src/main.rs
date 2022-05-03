mod fastaio;
use crate::fastaio::*;
use clap::Parser;

mod distance;
use crate::distance::tn93;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
	// input alignment
	#[clap(short, long)]
	input: String,
}

fn main() {
	let args = Args::parse();

//	 let (w, n) = align_dims(&args.input).unwrap();
	
//	 println!("{} {}", w, n);

//	let ba = populate_array(&args.input);

   let efra = populate_struct_array(&args.input).unwrap();

   let mut d: f64 = 0.0;

   println!("sequence1\tsequence2\tdistance");

   for i in 0..efra.len()-1 {
	   for j in i+1..efra.len() {
			d = tn93(&efra[i], &efra[j]);
			println!("{}\t{}\t{}", &efra[i].id, &efra[j].id, d)
	   }
   }

}
