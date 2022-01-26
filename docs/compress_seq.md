## Name

compress_seq - Normalize and compress a multi-FASTA sequence database

## Synopsis

```
compress_seq -i fasta_sequence_database [ options ]
```

## Description

compress_seq takes a multi-FASTA sequence database and splits it into separate sequence, header and index files. It can output the sequence data in a variety of forms, depending of the command line options given, including a DNA optimized normalized form, and a bit compressed form. It can also add an arbitrary end of sequence character to each sequence entry, and force each sequence character to uppercase.

If compress_seq is used to pre-process a sequence database for primer_match, the best performance will result from the use of the -n true option, to index and normalize the sequence database.

By default, compress_seq splits the multi-FASTA sequence database file db into three files: db.seq, db.hdr and db.idb. In order to make compress_seq more suitable for use in computational analysis pipelines, it only re-constructs the sequence database component files if the file system timestamps indicate that it is necessary. Note that the -F option can be used to force each component to be re-made.

## Command Line Options

-i fasta_sequence_database

Name of the multi-Fasta sequence database file to process. Required.

-I ( true | false )

Write fasta index file in binary format. Default: true.

-n ( true | false )

Create a normalized version of the sequence data. This creates additional files with suffixes .sqn and .tbl. Default: false.

-D ( true | false )

Optimize the normalized sequence data for DNA sequence. Default: true.

-z ( true | false )

Create a bit compressed normalized version of the sequence data. This creates additional files with suffixes .sqz and .tbz. Default: false.

-u ( true | false )

Force the sequence data to uppercase characters. Default: false.

-e ( true | false )

Add an end of sequence character after each entry from the multi-FASTA sequence database file. This can help ensure that text matching algorithms cannot find a match that straddles two FASTA entries. Default: true.

-E eos

Use eos as the end of sequence character. eos is the ascii code for the desired character, it may be specified as a decimal, octal or hexadecimal number. Default: 12 (newline).

-S ( true | false )

Insert end of sequence character before initial sequence entry. Default: false.

-F ( true | false )

Force each component of the compressed sequence database to be regenerated, even if the file timestamps indicate that this isn't necessary. Default: false.

-C ( true | false )

Cleanup unnecessary temporary files. Default: true.

-B

Use buffered standard I/O rather than mmap to stream through the sequence database. On some platforms, where the use of mmap is somewhat unpredictable, this option may make it possible to run compress_seq reliably. 

-v 

Output the release tag of the binary.

-h 

Command-line help.

## See Also

[primer_match](primer_match.html), [pcr_match](pcr_match.html)

## Author

[Nathan Edwards](http://edwardslab.bcmb.georgetown.edu/)
