# distance

A command-line program to calculate pairwise genetic distances within or between alignments of DNA sequences in fasta format.

## Installation

You will need to [install Rust](https://www.rust-lang.org/tools/install) first.

Then you should be able to install distance by running:

```
cargo install --git https://github.com/benjamincjackson/distance
```

and the binary will be installed somewhere in your `$PATH`:

```
distance -h
```

or clone the repository and build it:

```
git clone https://github.com/benjamincjackson/distance.git
cd distance
cargo build --release
```

and the binary gets built in the repo's directory:

```
./target/release/distance --version
```

## Usage

You can calculate all pairwise distances within a single alignment by providing a single input file, like:

```
distance alignment.fasta > distances.tsv
```

or, equivalently:

```
distance -i alignment.fasta -o distances.tsv
```

or, equivalently, but reading from `stdin` and printing to `stdout`:

```
cat alignment.fasta | distance
```

Or you can calculate all pairwise comparisons between two alignments like:

```
distance alignment1.fasta alignment2.fasta > distances2.tsv
```

or, equivalently:

```
distance -i alignment1.fasta alignment2.fasta -o distances2.tsv
```

Alignments read from `stdin` or provided as positional files (or provided to `-i`, which is equivalent) are read into memory by default (but see below).

You can also calculate all pairwise distances between one alignment in memory and one alignment streamed from disk, using the `-s / --stream` flag, like:

```
distance -i smallAlignment.fasta -s bigAlignment.fasta -o distances3.tsv
```

You can also force `-s` to accept input from `stdin` if you pass it `-`, like:

```
cat bigAlignment.fasta | distance smallAlignment.fasta -s - > distances3.tsv
```

The output is a tab-separated-value file with three columns:

```
> head distances.tsv
sequence1	sequence2	distance
seq1	seq2	0.000719991771
seq1	seq3	0.001887958258
seq1	seq4	0.001949663239
seq1	seq5	0.002853704658
seq1	seq6	0.002861871595
seq1	seq7	0.001748731312
seq1	seq8	0.002784940691
seq1	seq9	0.000788562416
seq1	seq10	0.001381740301
```

Different distance measures are available. These are `n` - the total number of nucleotide differences, `n_high` - should give exactly the same answer as `n` but might be faster for very high diversity datasets (?), or if comparing a small number of sequences in `-i` with a large number of sequences in `-s`, `raw` - the number of nucleotide differences _per site_, `jc69` - Jukes and Cantor's ([1969](https://books.google.co.uk/books?id=FDHLBAAAQBAJ&lpg=PA21&ots=bmgnXDW6mB&dq=jukes%20cantor%201969&lr&pg=PA34#v=onepage&q=jukes%20cantor%201969&f=false)) evolutionary distance, `k80` - Kimura's ([1980](https://doi.org/10.1007/bf01731581)) evolutionary distance, and `tn93` - Tamura and Nei's ([1993](https://doi.org/10.1093/oxfordjournals.molbev.a040023)) evolutionary distance.

Use the `-t` option to use spin up multiple threads for pairwise comparisons (in addition to a thread for i/o), e.g.:

```
distance -t 8 -m jc69 -i alignment.fasta -o jc69.tsv
```

and the `-b` option to tune the workload per thread, which may result in extra efficiency:

```
distance -t 8 -b 1000 -m jc69 -i alignment.fasta -o jc69.tsv
```

## Help

```
> distance -h
Usage:
       distance alignment.fasta
       cat alignment.fasta | distance
       distance alignment.fasta -o distances.tsv
       distance -t 8 -m jc69 alignment.fasta -o jc69.tsv
       distance alignment1.fasta alignment2.fasta > distances2.tsv
       distance -i smallAlignment.fasta -s bigAlignment.fasta -o distances3.tsv
       cat bigAlignment.fasta | distance smallAlignment.fasta -s - > distances3.tsv

Options:
  -i, --input [<input>...]     One or two input alignment files in fasta format. Loaded into memory. This flag can be omitted and the files passed as positional arguments
  -s, --stream <stream>        One input alignment file in fasta format. Streamed from disk (or stdin using "-s -"). Requires exactly one file also be loaded
  -m, --measure <measure>      Which distance measure to use [default: raw] [possible values: n, n_high, raw, jc69, k80, tn93]
  -o, --output <output>        Output file in tab-separated-value format. Omit this option to print to stdout
  -t, --threads <threads>      How many threads to spin up for pairwise comparisons [default: 1]
  -b, --batchsize <batchsize>  Try setting this >(>) 1 if you are struggling to get a speedup when adding threads [default: 1]
  -l, --licenses               Print licence information and exit
  -h, --help                   Print help
  -V, --version                Print version
```

## Acknowledgements

This program incorporates [rust-bio](https://rust-bio.github.io/). This program makes use of the [bitwise coding scheme for nucleotides](http://ape-package.ird.fr/misc/BitLevelCodingScheme.html) by Emmanuel Paradis, as used in ape ([Paradis, 2004](https://doi.org/10.1093/bioinformatics/btg412)). Equation (7) in Tamura and Nei ([1993](https://doi.org/10.1093/oxfordjournals.molbev.a040023)) is also rearranged according to ape's source code.
