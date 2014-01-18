package pdbx

import "github.com/TuftsBCB/seq"

func aminoCodonToLetter(codon string) seq.Residue {
	switch codon {
	case "ALA":
		return 'A'
	case "ARG":
		return 'R'
	case "ASN":
		return 'N'
	case "ASP":
		return 'D'
	case "CYS":
		return 'C'
	case "GLU":
		return 'E'
	case "GLN":
		return 'Q'
	case "GLY":
		return 'G'
	case "HIS":
		return 'H'
	case "ILE":
		return 'I'
	case "LEU":
		return 'L'
	case "LYS":
		return 'K'
	case "MET":
		return 'M'
	case "PHE":
		return 'F'
	case "PRO":
		return 'P'
	case "SER":
		return 'S'
	case "THR":
		return 'T'
	case "TRP":
		return 'W'
	case "TYR":
		return 'Y'
	case "VAL":
		return 'V'
	case "UNK":
		return 'X'
	default:
		panic(sf("Unknown amino acid codon: %s", codon))
	}
}
