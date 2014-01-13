package pdbx

import (
	"github.com/TuftsBCB/seq"
	"github.com/TuftsBCB/structure"

	"github.com/BurntSushi/cif"
)

// Entry corresponds to a single data block inside a PDBx/mmCIF formatted file.
// Usually, a file contains only one data block.
//
// Note that an Entry does *not* encapsulate all information in a PDBx/mmCIF
// file. In fact, it's only a small subset of the information mostly related
// to the structures and sequences in the file. For this reason, access to the
// raw CIF data block from the file is exposed in this type.
type Entry struct {
	// The underlying CIF file. This provides raw access to attributes of a
	// PDBx file that are not captured by the types in this package.
	CIF *cif.DataBlock

	// The four letter PDB identifier corresponding to this entry.
	// It is stored in lowercase.
	IdCode string

	// A list of all structures (chains) in this entry.
	// Note that each structure may contain multiple models.
	Chains []*Chain
}

// Chain corresponds to a single unique unit of structure within this PDB file.
type Chain struct {
	// A pointer back to its corresponding PDB entry.
	Entry *Entry

	// The entity that this structure corresponds to.
	EntityId int

	// A single character that uniquely identifies this structure in an entry.
	Ident byte

	// The type of this structure, corresponding precisely to the value
	// in the _entity_poly.type field.
	Type string

	// The sequence as specified in the _entity_poly_seq category (NOT the
	// _entity_poly.pdbx_seq_one_letter_code data item).
	Seq []seq.Residue

	// All models corresponding to this structure. There may be zero models,
	// although in the vast majority of cases, there is precisely one (but
	// there may be more, particularly for NMR structures).
	Models []*Model
}

// Model corresponds to a single structural model in a PDB file. There may be
// more than one model for each chain, although in practice, there is usually
// only one.
type Model struct {
	// A pointer back to its corresponding PDB entry.
	Entry *Entry

	// A pointer back to its corresponding chain.
	Chain *Chain

	// The model number (usually starting at 1).
	Num int

	// The ATOM records grouped by sequence residue.
	// Note that this slice has exactly the same size as Chain.Seq.
	// Since not every residue has a corresponding ATOM record, some
	// values in this slice may be nil.
	Sites []*Site

	// For convenience, the coordinates of each alpha-carbon atom found in
	// the ATOM records for this structure are stored here. Like Sites, this
	// is also guaranteed to have the same length as Chain.Seq, and since
	// not every reside has an ATOM record, some values may be nil.
	//
	// Note that if more than one alpha-carbon ATOM record exists for a
	// single residue, the first set of coordinates is used. (However, it is
	// not omitted from the set of atoms in each value in Sites.)
	AlphaCarbons []*structure.Coords
}

// Site corresponds to the set of ATOM records for a particular residue in
// a structure's sequence.
type Site struct {
	// A one-letter code corresponding to the residue for this ATOM record.
	Residue seq.Residue

	// The set of atoms (in no particular order) for this site.
	Atoms []Atom
}

// Atom represents a single ATOM record.
type Atom struct {
	Name string
	Het  bool
	structure.Coords
}
