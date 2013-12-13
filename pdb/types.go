package pdb

import (
	"fmt"

	"github.com/TuftsBCB/seq"
	"github.com/TuftsBCB/structure"
)

type Entry struct {
	Path   string
	IdCode string
	Chains []*Chain

	// SCOP is set whenever we see an identifier that looks like a
	// SCOP id. We use this to determine how to satisfy the Bower interface,
	// so that each entry has a unique ID.
	// Similarly for CATH.
	Scop string
	Cath string
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
	structure.Coords
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
//
// IsProtein also returns true if there are no SEQRES records.
func (c Chain) IsProtein() bool {
	return c.SeqType == -1 || (c.SeqType == SeqProtein && len(c.Models) > 0)
}

// SequenceCaAtomSlice attempts to extract a contiguous slice of alpha-carbon
// ATOM records based on *residue* index. Namely, if a contiguous slice cannot
// be found, nil is returned. If there is more than one model, the first model
// is used.
func (c Chain) SequenceCaAtomSlice(start, end int) []structure.Coords {
	return c.Models[0].SequenceCaAtomSlice(start, end)
}

// SequenceCaAtoms returns a slice of all Ca atoms for the chain in
// correspondence with the sequence in SEQRES (automatically using the first
// model).
//
// See Model.SequenceCaAtoms for the deets.
func (c Chain) SequenceCaAtoms() []*structure.Coords {
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
func (c Chain) CaAtoms() []structure.Coords {
	return c.Models[0].CaAtoms()
}

// AsSequence returns the chain as a sequence with an appropriate name.
// (e.g., 1tcfA)
func (c *Chain) AsSequence() seq.Sequence {
	name := fmt.Sprintf("%s%c", c.Entry.IdCode, c.Ident)
	return seq.Sequence{name, c.Sequence}
}

// SequenceCaAtomSlice attempts to extract a contiguous slice of alpha-carbon
// ATOM records based on *residue* index. Namely, if a contiguous slice cannot
// be found, nil is returned.
func (m Model) SequenceCaAtomSlice(start, end int) []structure.Coords {
	residues := m.SequenceCaAtoms()
	if start < 0 || start >= end || end > len(residues) {
		panic(fmt.Sprintf(
			"Invalid range [%d, %d). Must be in [%d, %d).",
			start, end, 0, len(residues)))
	}
	atoms := make([]structure.Coords, end-start)
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
func (m Model) SequenceCaAtoms() []*structure.Coords {
	mapping := m.SequenceAtoms()
	cas := make([]*structure.Coords, len(mapping))
	for i := range mapping {
		if mapping[i] != nil {
			if coords, ok := mapping[i].Ca(); ok {
				cas[i] = &coords
			}
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
func (m Model) CaAtoms() []structure.Coords {
	cas := make([]structure.Coords, 0, len(m.Residues))
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
func (r Residue) Ca() (structure.Coords, bool) {
	for _, atom := range r.Atoms {
		if atom.Name == "CA" {
			return atom.Coords, true
		}
	}
	return structure.Coords{}, false
}
