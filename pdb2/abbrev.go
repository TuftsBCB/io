package pdb2

import (
	"fmt"

	"github.com/BurntSushi/bcbgo/seq"
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
		return getAmino(abbrev)
	case SeqDeoxy:
		return getDeoxy(abbrev)
	case SeqRibo:
		return getRibo(abbrev)
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

func getAmino(threeAbbrev string) (seq.Residue, error) {
	if v, ok := aminoMap[threeAbbrev]; ok {
		return v, nil
	}
	return 'X', fmt.Errorf("Unknown amino residue: %s", threeAbbrev)
}

func getDeoxy(twoAbbrev string) (seq.Residue, error) {
	if v, ok := deoxyMap[twoAbbrev]; ok {
		return v, nil
	}
	return 'X', fmt.Errorf("Unknown deoxy residue: %s", twoAbbrev)
}

func getRibo(oneAbbrev string) (seq.Residue, error) {
	if v, ok := riboMap[oneAbbrev]; ok {
		return v, nil
	}
	return 'X', fmt.Errorf("Unknown ribo residue: %s", oneAbbrev)
}
