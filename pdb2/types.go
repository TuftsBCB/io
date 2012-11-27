package pdb2

import (
	"github.com/BurntSushi/bcbgo/seq"
)

type Entry struct {
	Path   string
	IdCode string
	Chains []Chain
}

type Chain struct {
	Entry    *Entry
	Ident    byte
	SeqType  SequenceType
	Sequence []seq.Residue
	Models   []Model
}

type Model struct {
	Entry    *Entry
	Chain    *Chain
	Num      int
	Residues []Residue
}

type Residue struct {
	Name        seq.Residue
	SequenceNum int

	// This is only populated when parsing XML files.
	// It's an index into Chain.Sequence (the SEQRES sequence).
	// This should be -1 when it's unavailable. In which case, we have
	// to guess the correspondence between SEQRES and ATOM with SequenceNum.
	ResidueIndex int
	Atoms        []Atom
}

type Atom struct {
	Name    string
	X, Y, Z float64
}

// SequenceCaAtoms returns a slice of all Ca atoms for the chain in
// correspondence with the sequence in SEQRES. (Only the first model is used.)
// Note that a slice of pointers is returned, since not all residues
// necessarily correspond to a alpha-carbon ATOM.
//
// If the ResidueIndex is not present in 'Residue' (i.e., read from a regular
// PDB file), then this correspondence is subject to corruption. SequenceCaAtoms
// does its best not to produce any false positives, but it's still possible.
//
// When ResidueIndex is present (in the PDB XML files), then each ATOM record
// precisely corresponds to a residue in SEQRES, and therefore the
// correspondence will be correct. (But could still contain gaps, for residues
// that did not crystallize.)
//
// In sum, a list of atom pointers is returned with length equal to the number
// of residues in the SEQRES record for this chain. Some pointers may be nil.
func (c Chain) SequenceCaAtoms() []*[3]float64 {
	cas := make([]*[3]float64, len(c.Sequence))
	for _, r := range c.Models[0].Residues {
		ca := r.Ca()
		if ca == nil {
			continue
		}

		// If we have a residue index, then we know precisely which residue
		// we want.
		if r.ResidueIndex != -1 {
			cas[r.ResidueIndex-1] = ca
			continue
		}

		// On the other hand, the SequenceNum (in ATOM records) is crap. It's
		// basically up to the author of the particular PDB entry. i.e., it
		// doesn't even have to be a number! (Then why the fuck is it called
		// a residue *sequence* number?)
		sn := r.SequenceNum - 1
		if sn >= 0 && sn < len(c.Sequence) && c.Sequence[sn] == r.Name {
			cas[sn] = ca
			continue
		}
	}
	return cas
}

// CaAtoms returns all alpha-carbon atoms in the chain. If there is more than
// one model, only the first model is used.
func (c Chain) CaAtoms() [][3]float64 {
	cas := make([][3]float64, 0, len(c.Models[0].Residues))
	for _, r := range c.Models[0].Residues {
		if ca := r.Ca(); ca != nil {
			cas = append(cas, *ca)
		}
	}
	return cas
}

// Ca returns the alpha-carbon atom in this residue.
// If one does not exist, nil is returned.
func (r Residue) Ca() *[3]float64 {
	for _, atom := range r.Atoms {
		if atom.Name == "CA" {
			coords := atom.Coords()
			return &coords
		}
	}
	return nil
}

func (a Atom) Coords() [3]float64 {
	return [3]float64{a.X, a.Y, a.Z}
}
