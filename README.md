# TF-IDF method for LGT detection.
This is the source code of TF-IDF program. This program is specific for LGT detection. 

NOTE: This version of the TF-IDF program has been modified to fix a memory leak that was occuring during the final stage when the results were being written to file. The way in which the program is called on the command line has also changed and the log messages has been made more informative. The results produced by this version should be identical as the core algorithm is unchanged.

## Input Files
The input should be seqeuence part of a FASTA file. e.g. there are 5 sequences for LGT detection, the input file should be (minus the '#' comments which are added for this example):

ACCCGGGGTTTTCAAA # seq1
ACCCTTGGGGCCCAAT # seq2
ACCCGGGTTCCAAAAA # seq3
ACGGTTGGGGCCCAAT # seq4
ACCCTTGGGGCCCCCT # seq5

The header of each sequence is not required.

The group information should be provided in a separate file. e.g. there are 5 sequences in a dataset, and the first 2 sequences are in one group, the rest 3 sequences are in another group. Then the group information should be:

2     # amount of groups
1 2   # sequence IDs in group 1
3 4 5 # sequence IDs in group 2


## To Compile
The program is written in C++. It can be compiled by GCC 4.4.7 (has also been tested with 9.4.0). If you have any question, please contact me by email y.cong@uq.edu.au.

`g++ -O3 -o tf-idf tf-idf.cpp`


## To Run
An example of how to run this program from command line:

`tf-idf seqfile.txt 40 0.05 < speciesInfo.txt`

Where:
- `seqfile.txt` is the file with your sequences (in the format shown above)
- `40` is the *k*-mer size to use (40 is the default for the original version of program; we suggest you keep this set at 40)
- `0.05` is the significant level used for the significant test
- `speciesInfo.txt` is the group information file (in the format shown above)


