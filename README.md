# distance

A command-line program to calculate pairwise genetic distances within or between alignments of DNA sequences in fasta format.

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

You can calculate all pairwise distances within a single alignment by providing a single input file, like:

`distance -i alignment.fasta -o distances.tsv`

or reading from `stdin` and writing to `stdout`:

`cat alignment.fasta | distance`

or all pairwise comparisons between two alignments like:

`distance -i alignment1.fasta alignment2.fasta -o distances2.tsv`

All alignments provided to `-i` (or when you read from `stdin`) are read into memory.

You can also calculate all pairwise distances between one alignment in memory and one alignment streamed from disk, using the `-s / --stream` flag, like:

`distance -i aSmallAlignment.fasta -s aBigAlignment.fasta -o distances3.tsv`

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

Different distance measures are available. These are `n` - the total number of nucleotide differences, `n_high` - should give exactly the same answer as `n` but might be faster for very high diversity datasets (?), or if comparing a small number of sequences in `-i` with a large number of sequences in `-s`, `raw` - the number of nucleotide differences _per site_, `jc69` - Jukes and Cantor's ([1969](https://books.google.co.uk/books?id=FDHLBAAAQBAJ&lpg=PA21&ots=bmgnXDW6mB&dq=jukes%20cantor%201969&lr&pg=PA34#v=onepage&q=jukes%20cantor%201969&f=false)) evolutionary distance, `k80` - Kimura's ([1980](https://doi.org/10.1007/bf01731581)) evolutionary distance, and `tn93` - Tamura and Nei's ([1993](https://doi.org/10.1093/oxfordjournals.molbev.a040023)) evolutionary distance.

Use the `-t` option to use spin up multiple threads for pairwise comparisons (in addition to a thread for i/o), e.g.:

```
distance -t 8 -m jc69 -i aligned.fasta -o jc69.tsv
```

and the `-b` option to tune the workload per thread, which may result in extra efficiency:

```
distance -t 8 -b 1000 -m jc69 -i aligned.fasta -o jc69.tsv
```

## Help

```
> distance -h
distance 0.1.2

USAGE:
    distance [OPTIONS]

OPTIONS:
    -b, --batchsize <batchsize>    Try setting this >(>) 1 if you are struggling to get a speedup
                                   when adding threads [default: 1]
    -h, --help                     Print help information
    -i, --input <input>...         Input alignment file(s) in fasta format. Loaded into memory
                                   [default: stdin]
    -l, --licenses                 Print licence information and exit
    -m, --measure <measure>        Which distance measure to use [default: raw] [possible values: n,
                                   n_high, raw, jc69, k80, tn93]
    -o, --output <output>          Output file in tab-separated-value format. Omit this option to
                                   print to stdout
    -s, --stream <stream>          Input alignment file in fasta format. Streamed from disk.
                                   Requires exactly one file also be specifed to -i
    -t, --threads <threads>        How many threads to spin up for pairwise comparisons [default: 1]
    -V, --version                  Print version information
```

## Acknowledgements

This program incorporates [rust-bio](https://rust-bio.github.io/). This program makes use of the [bitwise coding scheme for nucleotides](http://ape-package.ird.fr/misc/BitLevelCodingScheme.html) by Emmanuel Paradis, as used in ape ([Paradis, 2004](https://doi.org/10.1093/bioinformatics/btg412)). Equation (7) in Tamura and Nei ([1993](https://doi.org/10.1093/oxfordjournals.molbev.a040023)) is also rearranged according to ape's source code.
