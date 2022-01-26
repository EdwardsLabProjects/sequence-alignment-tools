## Introduction

The Primer Match suite of tools is designed to find and count exact and near exact matches of short oligonucleotide sequences in large genomic databases. All matches to the oligos will be output, the tools guarantee a complete enumeration of all matches consistent with the search options. Substitutions, insertions and deletions can be prohibited in the start, end, 5' or 3' bases of oligos. Many options for constraining acceptable alignments and input/output formats are provided. The tools automatically optimize the sequence search strategies to match the search parameters.

## Installation

See the binary [releases](https://github.com/EdwardsLabProjects/sequence-alignment-tools) in this repository. 
    
## Programs

### primer_match
[primer_match](primer_match.html) finds and counts exact and near exact instances of short DNA sequences, usually primers, in a (much) larger DNA sequence database such as the human genome. 

### pcr_match
[pcr_match](pcr_match.html) finds pairs of short DNA sequences, usually primers, in a (much) larger DNA sequence database such as the human genome. 

### compress_seq
[compress_seq](compress_seq.html) reformats multi-FASTA sequence databases for efficient searching by primer_match and pcr_match. 

## Examples

* Determine the set of PCR primers that occur exactly once in the human genome

  1. Normalize the genome sequence database using compress_seq

     ```
     % compress_seq -i genome.fasta -n true
     ```

  2. Count primer occurrences in the genome that match with at most one string edit, but in which the 5' most 7 bases must match exactly

     ```
     % primer_match -i genome.fasta -P primers.txt -r -k 1 -5 7 -c -a
     ```

* Find all occurrences of PCR primer pairs matching exactly, oriented towards each other, with maximum amplicon length 1000, and extract the amplicons in fasta format.

     ```
     % compress_seq -i genome.fasta -n true
     % pcr_match -i genome.fasta -P primers.txt -M 1000 -A ">%i /len=%l /Ns=%N /edits=(%>e,%<e)\n%@\n"
     ```

