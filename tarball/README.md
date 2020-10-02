# Distribution Sensitive Bucketing (DSB)
Author(s): Chengze Shen, Mihir Mongia, Arash Gholami Davoodi, Guillaume Marcais, Hosein Mohimani

## What It Does
DSB is a C++ based program to find overlaps among sequences and alignments of queries in a given genome. There are two input files, reference and query in _fasta_ format. Reference could be a set of reads or a reference genome. For each query sequence, the program outputs the reference sequences that overlap with it. In cases where query and reference are the same files, the program discards self-overlaps.

The goal is to find as many as true overlapping sequences and alignments while minimizing false positives.

## Requirements
### Linux
We tested our program on **Ubuntu 18.04** with **g++ 7.5.0** and above and ISO standard **-std=c++11**.
### macOS
We also tested our program on **macOS 10.14.6** with **g++ 4.2.1** and above, ISO standard **-std=c++11** and **Apple LLVM 10.0.1**.

## Installation
1. Make sure the correct version of C++ compiler has been installed.
2. _**By default**_, _**use**_ `make` _**to generate everything (DSBMain, DataGeneration).**_
3. There are two binary executables you could generate.
    * To generate **DSBMain** (for DSB), please `make main`.
    * To generate **DataGeneration** (for generating simulation data), please `make gen`.
4. To remove installation, please `make clean`. This will remove everything generated after your initial download, except data files you generated using **DataGeneration**.

## How To Run
### **DSBMain**
```bash
./DSBMain -q [query file] -r [reference file] -i [insertion rate] -d [deletion rate] -m [mutation rate] -a [add threshold] -k [kill threshold] -o [name] {-vh}

-h              Print this block of information.
-v              Verbose mode.
-q [path/to/q]  Path to the query file.
-r [path/to/r]  Path to the reference file.
-i [0<=i<1]     Insertion rate.
-d [0<=d<1]     Deletion rate.
-m [0<=e<0.5]   Mutation rate when neither insertion/deletion happens.
-a [a>0]        Threshold for a node to be considered a bucket.
-k [k>0]        Threshold for a node to be pruned.
-o [name]       Specify output file name.
```
To view the full helper message with command line, please use `./DSBMain -h`.

##### _Example 1_
Starting with query **data/pacbio_reads_5000.fasta** and reference **data/ecoli_genome_full.fasta**, we use **12\% insertion rate, 2\% deletion rate, and 1\% mutation rate** for the PacBio sequencing data:
```bash
./DSBMain -q data/pacbio_reads_5000.fasta -r data/ecoli_genome_full.fasta -i 0.12 -d 0.02 -m 0.01 -a 25000 -k 250000000
```

To search reads against reads, we use the following command:
```bash
./DSBMain -q data/pacbio_reads_5000.fasta -r data/pacbio_reads_5000.fasta -i 0.12 -d 0.02 -m 0.01 -a 25000 -k 250000000
```

Without specifying the output name, the program will print the results to a default file named **output.txt**.
### **DataGeneration**
```bash
./DataGeneration -i [insertion rate] -d [deletion rate] -m [mutation rate] -n [number of sequences] -s [initial length of a sequence] -p [path] -vh

-h              Print this block of information.
-v              Verbose mode.
-i [0<=i<1]     Insertion rate.
-d [0<=d<1]     Deletion rate.
-m [0<=e<0.5]   Mutation rate when neither insertion/deletion happens.
-n [n>0]        Number of sequence pairs to generate.
-s [s>0]        Initial length of a sequence pair.
-p [path]       Where generated data will be written to. Default: data/
```
To view the full helper message with command line, please use `./DataGeneration -h`.

##### _Example 2_
We generate **1000 sequence pairs** (each sequence pair contains a reference r and a query q) with **insertion, deletion, and mutation rate of 5%** and length of approximately **100 bp** using the following command:
```bash
./DataGeneration -i 0.05 -d 0.05 -m 0.05 -n 1000 -s 100 -p data
```
The command will generate 2 files under **data/** directory, namely:
* q_1000_100_0.05_0.05.fasta
* r_1000_100_0.05_0.05.fasta

where **q_1000_100_0.05_0.05.fasta** corresponds to the queries in the sequence pairs and **r_1000_100_0.05_0.05.fasta** corresponds to the references.
