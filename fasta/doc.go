/*
Package fasta provides routines for reading and writing FASTA files. Routines
are also provided to read and write aligned fasta files.

The format used is the one described by NCBI:
http://blast.ncbi.nlm.nih.gov/blastcgihelp.shtmlj

By default, sequences are checked to make sure they contain only valid
characters: a-z, A-Z, * and -. All lowercases letters are translated to their
upper case equivalent.
*/
package fasta
