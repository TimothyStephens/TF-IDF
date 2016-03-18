# TF-IDF method for LGT detection.
This is the source code of TF-IDF program. This is program is specific for LGT detection. The input should be seqeuence part of a FASTA file. e.g. there are 5 sequences for LGT detection, the input file should be:

ACCCGGGGTTTTCAAA # seq1

ACCCTTGGGGCCCAAT # seq2

ACCCGGGTTCCAAAAA # seq3

ACGGTTGGGGCCCAAT # seq4

ACCCTTGGGGCCCCCT # seq5

The header of each sequence is not required.

The group information should be provided in a separate file. e.g. there are 5 sequences in a dataset, and the first 2 sequences are in one group, the rest 8 sequences are in another group. Then the group information should be:

2     # amount of groups

1 2   # sequence IDs in group 1

3 4 5 # sequence IDs in group 2



The program is written by C++. It can be compiled by GCC 4.4.7.
