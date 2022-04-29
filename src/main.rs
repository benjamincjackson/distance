mod fastaio;
use crate::fastaio::*;
use clap::Parser;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    // input alignment
    #[clap(short, long)]
    input: String,
}

fn main() {
    let args = Args::parse();

    let (w, n) = align_dims(&args.input).unwrap();
    
    println!("{} {}", w, n);

   let ba = populate_array(&args.input);
}
