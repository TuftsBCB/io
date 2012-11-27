package pdb2

import (
	"fmt"

	"github.com/BurntSushi/bcbgo/seq"
)

var aminoMap = map[string]seq.Residue{
	"ALA": 'A', "ARG": 'R', "ASN": 'N', "ASP": 'D', "CYS": 'C',
	"GLU": 'E', "GLN": 'Q', "GLY": 'G', "HIS": 'H', "ILE": 'I',
	"LEU": 'L', "LYS": 'K', "MET": 'M', "PHE": 'F', "PRO": 'P',
	"SER": 'S', "THR": 'T', "TRP": 'W', "TYR": 'Y', "VAL": 'V',
	"SEC": 'U', "PYL": 'O',
	"UNK": 'X', "ACE": 'X', "NH2": 'X',
	"ASX": 'X', "GLX": 'X',
	"MSE": 'M', "CSA": 'C', "LLP": 'K', "CSW": 'C', "STE": 'X',
}

var deoxyMap = map[string]seq.Residue{
	"DA": 'A', "DC": 'C', "DG": 'T', "DI": 'I',
}

var riboMap = map[string]seq.Residue{
	"A": 'A', "C": 'C', "G": 'G', "U": 'U', "I": 'I',
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

func getAbbrev(abbrev string) seq.Residue {
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
