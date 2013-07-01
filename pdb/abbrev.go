package pdb

import (
	"fmt"

	"github.com/TuftsBCB/seq"
)

var aminoMap = map[string]seq.Residue{
	"UNK": 'X',
	"ALA": 'A', "ARG": 'R', "ASN": 'N', "ASP": 'D', "CYS": 'C',
	"GLU": 'E', "GLN": 'Q', "GLY": 'G', "HIS": 'H', "ILE": 'I',
	"LEU": 'L', "LYS": 'K', "MET": 'M', "PHE": 'F', "PRO": 'P',
	"SER": 'S', "THR": 'T', "TRP": 'W', "TYR": 'Y', "VAL": 'V',

	// Residues that are not classified as "missing" and are not
	// listed in MODRES.
	"ASX": 'X', "GLX": 'X', "DLE": 'X',

	// these are nucleotides. grrr.
	"DOP": 'X',
	"8OG": 'X',

	// misc
	"NH2": 'X',
}

var deoxyMap = map[string]seq.Residue{
	"DA": 'A', "DC": 'C', "DG": 'G', "DT": 'T', "DI": 'I', "DU": 'U',
}

var riboMap = map[string]seq.Residue{
	"A": 'A', "C": 'C', "G": 'G', "U": 'U', "I": 'I', "T": 'T',

	// wtf?
	"N": 'X',
}

type SequenceType int

const (
	SeqProtein = iota
	SeqDeoxy
	SeqRibo
)

func (typ SequenceType) String() string {
	switch typ {
	case SeqProtein:
		return "Protein"
	case SeqDeoxy:
		return "Deoxy"
	case SeqRibo:
		return "Ribo"
	}
	panic(fmt.Sprintf("Unknown sequence type: %d", typ))
}

func getAbbrev(abbrev string) (seq.Residue, error) {
	typ := getAbbrevType(abbrev)
	switch typ {
	case SeqProtein:
		return getAmino(abbrev), nil
	case SeqDeoxy:
		return getDeoxy(abbrev), nil
	case SeqRibo:
		return getRibo(abbrev), nil
	}
	panic(fmt.Sprintf("Unknown sequence type: %d", typ))
}

func getAbbrevType(abbrev string) SequenceType {
	switch len(abbrev) {
	case 3:
		return SeqProtein
	case 2:
		return SeqDeoxy
	case 1:
		return SeqRibo
	}
	panic(fmt.Sprintf("Unknown abbreviation type: %s (length: %d)",
		abbrev, len(abbrev)))
}

func getAmino(threeAbbrev string) seq.Residue {
	if v, ok := aminoMap[threeAbbrev]; ok {
		return v
	}
	return 'X'
}

func getDeoxy(twoAbbrev string) seq.Residue {
	if v, ok := deoxyMap[twoAbbrev]; ok {
		return v
	}
	return 'X'
}

func getRibo(oneAbbrev string) seq.Residue {
	if v, ok := riboMap[oneAbbrev]; ok {
		return v
	}
	return 'X'
}
