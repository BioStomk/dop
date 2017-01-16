# dop
A program for drawing dot plot by finding k-mer matches of two DNA sequences

## Installation

```sh
$ git clone https://github.com/stomk/dop.git
$ cd dop
$ cd src && make
```

## Getting started
Move into the example directory

```sh
$ cd example
```

and type

```
$ ../dop.py example.fasta
```

This will give example.fasta.20.png (20 indicates default k-mer size (20 bp)).

<img width=600 src="https://github.com/stomk/dop/blob/images/example/example.fasta.20.png", alt="default">

```
$ ../dop.py example.fasta -c
```

Will give a output file of the same name as above, with two-color dots (red dots represent forward matches and green dots reverse matches).

<img width=600 src="https://github.com/stomk/dop/blob/images/example/example.fasta.20c.png", alt="-c">

You can change the k-mer size with `-k` option.

```
$ ../dop.py example.fasta -k 40
```

<img width=600 src="https://github.com/stomk/dop/blob/images/example/example.fasta.40.png", alt="-k 40">

For more usage, type `../dop.py -h`.
