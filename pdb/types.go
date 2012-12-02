package pdb

import (
	"fmt"

	"github.com/BurntSushi/bcbgo/seq"
)

type Entry struct {
	Path   string
	IdCode string
	Chains []*Chain
}

type Chain struct {
	Entry    *Entry
	Ident    byte
	SeqType  SequenceType
	Sequence []seq.Residue
	Models   []*Model
	Missing  []*Residue
}

type Model struct {
	Entry    *Entry
	Chain    *Chain
	Num      int
	Residues []*Residue
}

type Residue struct {
	Name          seq.Residue
	SequenceNum   int
	InsertionCode byte
	Atoms         []Atom
}

type Atom struct {
	Name string
	Het  bool
	Coords
}

type Coords struct {
	X, Y, Z float64
}

// Chain returns a chain with the given identifier.
// If such a chain does not exist, nil is returned.
func (entry *Entry) Chain(ident byte) *Chain {
	for i := range entry.Chains {
		if entry.Chains[i].Ident == ident {
			return entry.Chains[i]
		}
	}
	return nil
}

// OneChain returns a single chain in the PDB file. If there is more than one
// chain, OneChain will panic. This is convenient when you expect a PDB file to
// have only a single chain, but don't know the name.
func (entry *Entry) OneChain() *Chain {
	if len(entry.Chains) != 1 {
		panic(fmt.Sprintf("OneChain can only be called on PDB entries with "+
			"ONE chain. But the '%s' PDB entry has %d chains.",
			entry.Path, len(entry.Chains)))
	}
	return entry.Chains[0]
}

// IsProtein returns true if the chain consists of amino acids.
func (c Chain) IsProtein() bool {
	return c.SeqType == SeqProtein && len(c.Models) > 0
}

// SequenceCaAtomSlice attempts to extract a contiguous slice of alpha-carbon
// ATOM records based on *residue* index. Namely, if a contiguous slice cannot
// be found, nil is returned. If there is more than one model, the first model
// is used.
func (c Chain) SequenceCaAtomSlice(start, end int) []Coords {
	return c.Models[0].SequenceCaAtomSlice(start, end)
}

// SequenceCaAtoms returns a slice of all Ca atoms for the chain in
// correspondence with the sequence in SEQRES (automatically using the first
// model).
//
// See Model.SequenceCaAtoms for the deets.
func (c Chain) SequenceCaAtoms() []*Coords {
	return c.Models[0].SequenceCaAtoms()
}

// SequenceAtoms returns a slice of all residues for the chain in
// correspondence with the sequence in SEQRES (automatically using the first
// model). Namely, the mapping is sparse, since not all SEQRES residues have
// an ATOM record.
//
// See Model.SequenceCaAtoms for the deets.
func (c Chain) SequenceAtoms() []*Residue {
	return c.Models[0].SequenceAtoms()
}

// CaAtoms returns all alpha-carbon atoms in the chain. If there is more than
// one model, only the first model is used.
func (c Chain) CaAtoms() []Coords {
	return c.Models[0].CaAtoms()
}

// SequenceCaAtomSlice attempts to extract a contiguous slice of alpha-carbon
// ATOM records based on *residue* index. Namely, if a contiguous slice cannot
// be found, nil is returned.
func (m Model) SequenceCaAtomSlice(start, end int) []Coords {
	residues := m.SequenceCaAtoms()
	atoms := make([]Coords, end-start)
	for i, cai := 0, start; cai < end; i, cai = i+1, cai+1 {
		if residues[cai] == nil {
			return nil
		}
		atoms[i] = *residues[cai]
	}
	return atoms
}

// SequenceCaAtoms returns a slice of all Ca atoms for the model in
// correspondence with the sequence in SEQRES.
// Note that a slice of pointers is returned, since not all residues
// necessarily correspond to a alpha-carbon ATOM.
//
// This method can proceed in one of two ways. First, if "REMARK 465" is
// present in the PDB file, it will be used to determine the positions of the
// holes in the sequence (i.e., residues in SEQRES without an ATOM record).
// This method is generally reliable, since REMARK 465 lists all residues
// in SEQRES that don't have an ATOM record. This will fail if there are any
// unreported missing residues.
//
// If "REMARK 465" is absent, then we have to rely on the order of ATOM records
// to correspond to a residue index in the SEQRES sequence. This will fail
// with an error if there are any unreported missing residues.
//
// Generally, false positives are limited by returning errors if corruption
// is detected. However, false positives can be returned in pathological cases
// (like long strings of low complexity regions or UNKNOWN amino acids), but
// they are rare. Probably on the order of a handful in the entire PDB.
//
// In sum, a list of atom pointers is returned with length equal to the number
// of residues in the SEQRES record for this model. Some pointers may be nil.
func (m Model) SequenceCaAtoms() []*Coords {
	mapping := m.SequenceAtoms()
	cas := make([]*Coords, len(mapping))
	for i := range mapping {
		if mapping[i] != nil {
			cas[i] = mapping[i].Ca()
		}
	}
	return cas
}

// SequenceAtoms is just like SequenceCaAtoms, except it returns the residues
// instead of the alpha-carbon coordinates directly. The advantage here is to
// get a mapping that isn't limited by the presence of alpha-carbon atoms.
//
// See SequenceCaAtoms for the deets.
func (m Model) SequenceAtoms() []*Residue {
	if len(m.Chain.Missing) > 0 {
		return m.seqAtomsChunksMerge()
	}
	return m.seqAtomsGuess()
}

// CaAtoms returns all alpha-carbon atoms in the model.
// This includes multiple alpha-carbon atoms belonging to the same residue.
// It does not include HETATMs.
func (m Model) CaAtoms() []Coords {
	cas := make([]Coords, 0, len(m.Residues))
	for _, r := range m.Residues {
		for _, atom := range r.Atoms {
			if atom.Name == "CA" && !atom.Het {
				cas = append(cas, atom.Coords)
			}
		}
	}
	return cas
}

// Ca returns the alpha-carbon atom in this residue.
// If one does not exist, nil is returned.
func (r Residue) Ca() *Coords {
	for _, atom := range r.Atoms {
		if atom.Name == "CA" {
			return &atom.Coords
		}
	}
	return nil
}

func (coords Coords) String() string {
	return fmt.Sprintf("%0.3f %0.3f %0.3f", coords.X, coords.Y, coords.Z)
}
