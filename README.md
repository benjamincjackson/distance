# distance

A command-line program to calculate pairwise evolutionary distances within or between alignments of DNA sequences in fasta format.

This program is a work in progress and may contain bugs. It was written primarily as a vehicle for learning Rust.

## Installation

You will need to [install Rust](https://www.rust-lang.org/tools/install) first.

Then you should be able to install distance by running:

```
cargo install --git https://github.com/benjamincjackson/distance
```

and check it worked:
```
distance -h # the binary has been installed somewhere in your $PATH
```

or clone the repository and build it:

```
git clone https://github.com/benjamincjackson/distance.git
cd distance
cargo build --release
```
and check it worked:
```
./target/release/distance --version # the binary gets built in the repo's directory
```

## Usage

You can calculate all pairwise distances within a single alignmet by providing a single input file, like:

`distance -i alignment.fasta -o distances.tsv`

or all pairwise comparisons between two files like:

`distance -i alignment1.fasta alignment2.fasta -o distances2.tsv`

All alignments are read into memory.

The output is a tab-separated-value file with three columns:

```
> head distances.tsv
sequence1	sequence2	distance
seq1	seq2	0.0007199917715226112
seq1	seq3	0.00188795825895922
seq1	seq4	0.0019496632399858206
seq1	seq5	0.0028537046587588104
seq1	seq6	0.0028618715950624093
seq1	seq7	0.0017487313125771498
seq1	seq8	0.002784940691077875
seq1	seq9	0.0007885624164295265
seq1	seq10	0.001381740301910256
```

Different distance measures are available. These are `n` - the total number of nucleotide differences, `raw` - the number of nucleotide differences _per site_, `jc69` - Jukes and Cantor's ([1969](https://books.google.co.uk/books?id=FDHLBAAAQBAJ&lpg=PA21&ots=bmgnXDW6mB&dq=jukes%20cantor%201969&lr&pg=PA34#v=onepage&q=jukes%20cantor%201969&f=false)) evolutionary distance, `k80` - Kimura's ([1980](https://doi.org/10.1007/bf01731581)) evolutionary distance, and `tn93` - Tamura and Nei's ([1993](https://doi.org/10.1093/oxfordjournals.molbev.a040023)) evolutionary distance.

Use the `-t` option to use spin up multiple threads for pairwise comparisons (in addition to a thread for i/o), e.g.:

```
distance -t 8 -m jc69 -i aligned.fasta -o jc69.tsv
```

`❯ distance -h`:

```
distance 0.1.0

USAGE:
    distance [OPTIONS] --input <input>... --output <output>

OPTIONS:
    -h, --help                 Print help information
    -i, --input <input>...     input alignment(s) in fasta format
    -m, --measure <measure>    which distance measure to use [default: raw] [possible values: n,
                               raw, jc69, k80, tn93]
    -o, --output <output>      output file in tab-separated-value format
    -t, --threads <threads>    how many threads to spin up for pairwise comparisons
    -V, --version              Print version information
```

