/*
Package pdb provides some minimal support for extracting information from PDB
files. Currently, this information includes the primary amino acid sequence for
each chain in the PDB file, as well as the start and stop residue indices of
each chain's ATOM records. Information not currently included: the 3
dimensional coordinates found in the ATOM records.

Anything in a PDB file that isn't protein related is ignored.

The 'Entry' type defined here is used in several places throughout bcbgo.
*/
package pdb
