# Distribution Sensitive Bucketing (DSB)
Author(s): Chengze Shen, Mihir Mongia, Arash Gholami Davoodi, Guillaume Marcais, Hosein Mohimani

## What It Does
DSB is a C++ based program to find overlaps among sequences and alignments of queries in a given genome. There are two input files: _X_ and _Y_ data in _.fasta_ format. The program will output, for every sequence in _Y_, the sequences in _X_ that overlaps with it. In case of a single file (_X_=_Y_), the program will avoid outputting self overlaps.

The goal is to find as many as true overlapping sequences and alignments while minimizing false positives.

## Requirements
The program is compiled with **g++ 4.2.1** and above, ISO standard **-std=c++11**. It was mainly tested on a MacBook Pro with **Apple LLVM 10.0.1**.

## Installation
1. Make sure the correct version of C++ compiler has been installed.
2. _**By default**_, _**use**_ `make` _**to generate everything (DSBMain, DataGeneration).**_
3. There are two binary executables you could generate.
    * To generate **DSBMain** (for DSB), please `make main`.
    * To generate **DataGeneration** (for generating simulation data), please `make gen`.
4. To remove installation, please `make clean`. This will remove everything generated after your initial download, except data files you generated using **DataGeneration**.

## How To Run
#### **DSBMain**
```bash
./DSBMain -x [X file] -y [Y file] -i [insertion rate] -d [deletion rate] -m [mutation rate] -a [threshold 1] -k [threshold 2] -o [name] {-vh}

-h              Print this block of information.
-v              Verbose mode.
-x [path/to/x]  Path to X data file.
-y [path/to/y]  Path to Y data file.
-i [0<=i<1]     Insertion rate.
-d [0<=d<1]     Deletion rate.
-m [0<=e<0.5]   Mutation rate when neither insertion/deletion happens.
-a [a>0]        Threshold for a node to be considered a bucket.
-k [k>0]        Threshold for a node to be pruned.
-o [name]       Specify output file name.
```
To view the full helper message with command line, please use `./DSBMain -h`.

##### _Example 1_
We have data files **data/sample_X.txt** and **data/sample_Y.txt**. We assume that these files represent DNA sequencing data with **5\% insertion rate, 5\% deletion rate, and 5\% mutation rate**. We want our data structure to have **threshold 1=1000** and **threshold 2=1000000**. Then, we need to run the program with the following command:
```bash
./DSBMain -x data/sample_X.txt -y data/sample_Y.txt -i 0.05 -d 0.05 -m 0.05 -a 1000 -k 1000000
```

Without specifying the output name, the program will print the results to a default file named **output.txt**.
#### **DataGeneration**
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
We want to generate data files under the **data/** directory. We want the **insertion rate, deletion rate, and mutation rate all to be 5\%**. We want to generate **1000 sequence pairs** (each sequence pair contains two sequences: one for X and one for Y. They are generated with the intended insertion/deletion/mutation rates), each to be approximately **100 bp in length**. Then, we need to run the program with the following command:
```bash
./DataGeneration -i 0.05 -d 0.05 -m 0.05 -n 1000 -s 100 -p data
```
The command will generate 4 files under **data/** directory, namely:
* X1000_100_0.05_0.05.txt
* Y1000_100_0.05_0.05.txt
* s_1000_100_0.05_0.05.fasta
* q_1000_100_0.05_0.05.fasta

The first two files are simply sequences extracted from _.fasta_ files without header lines that can be used for other purposes.
