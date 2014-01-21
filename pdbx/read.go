package pdbx

import (
	"fmt"
	"io"
	"os"

	"github.com/BurntSushi/cif"
	"github.com/TuftsBCB/seq"
	"github.com/TuftsBCB/structure"
)

var (
	ef = fmt.Errorf
	lf = func(f string, v ...interface{}) { fmt.Fprintf(os.Stderr, f, v...) }
	sf = fmt.Sprintf
)

// Read reads exactly one PDB entry from the reader given. If there are 0
// entries or more than 1 entry, an error is returned.
//
// An error is also returned if the reader could not be interpreted as a valid
// PDBx/mmCIF file (which must be a valid CIF file).
//
// This function is useful for reading standard PDBx/mmCIF files obtained from
// the PDB. (But in general, a PDBx/mmCIF file may contain more than one entry.)
func Read(r io.Reader) (*Entry, error) {
	entries, err := ReadAll(r)
	if err != nil {
		return nil, err
	} else if len(entries) != 1 {
		return nil, ef("Expected one PDB entry but got %d.", len(entries))
	}
	return entries[0], nil
}

// ReadAll reads all PDB entries from the reader provided. If you're reading
// PDBx/mmCIF files from the PDB, then use the Read function which guarantees
// the reader given only has a single entry.
//
// An error is returned if the reader could not be interpreted as a valid
// PDBx/mmCIF file (which must be a valid CIF file).
func ReadAll(r io.Reader) ([]*Entry, error) {
	cf, err := cif.Read(r)
	if err != nil {
		return nil, err
	}
	var entries []*Entry
	for _, block := range cf.Blocks {
		e, err := ReadCIFDataBlock(block)
		if err != nil {
			return nil, err
		}
		entries = append(entries, e)
	}
	return entries, nil
}

// ReadCIFDataBlock converts a PDBx/mmCIF data block to a PDB entry.
// It is exposed in the public interface so that clients can freely mix
// Entry objects and their corresponding underlying data blocks.
//
// An error is returned if the data block given does not correspond to a valid
// PDBx/mmCIF entry.
func ReadCIFDataBlock(b *cif.DataBlock) (*Entry, error) {
	e := &Entry{CIF: b}
	if err := e.read(b); err != nil {
		return nil, err
	}
	return e, nil
}

func (e *Entry) read(b *cif.DataBlock) error {
	if id, ok := b.Items["entry.id"]; !ok {
		return ef("Could not find ID code (in 'entry.id').")
	} else {
		e.Id = id.String()
	}
	e.Title = value(b, "struct.title").String()
	e.Descriptor = value(b, "struct.pdbx_descriptor").String()

	if err := e.readEntities(b); err != nil {
		return err
	}
	if err := e.readSequences(b); err != nil {
		return err
	}
	if err := e.readChains(b); err != nil {
		return err
	}
	if err := e.readAtomSites(b); err != nil {
		return err
	}
	return nil
}

func (e *Entry) readEntities(b *cif.DataBlock) error {
	loop := asLoop(b, "entity.id", "entity.type", "entity.formula_weight",
		"entity.pdbx_number_of_molecules")
	ids := loop[0].Strings()
	if ids == nil {
		return ef("There are no entities in this PDBx/mmCIF file.")
	}

	e.Entities = make(map[byte]*Entity, 2)
	types := loop[1].Strings()
	weights := loop[2].Floats()
	mols := loop[3].Floats()
	for i := range ids {
		ent := &Entity{
			Entry:  e,
			Id:     ids[i][0],
			Seq:    make([]seq.Residue, 0, 50),
			Chains: make(map[byte]*Chain, 2),
		}
		if types != nil {
			ent.Type = types[i]
		}
		if weights != nil {
			ent.FormulaWeight = weights[i]
		}
		if mols != nil {
			ent.NumOfMolecules = mols[i]
		}
		e.Entities[ent.Id] = ent
	}
	return nil
}

func (e *Entry) readSequences(b *cif.DataBlock) error {
	loop := asLoop(b, "entity_poly_seq.entity_id", "entity_poly_seq.mon_id")
	eids := loop[0].Strings()
	residues := loop[1].Strings()
	for i := range eids {
		ent := e.Entities[eids[i][0]]
		ent.Seq = append(ent.Seq, aminoCodonToLetter(residues[i]))
	}
	return nil
}

func (e *Entry) readChains(b *cif.DataBlock) error {
	loop := asLoop(b, "struct_asym.id", "struct_asym.entity_id")
	chainids := loop[0].Strings()
	eids := loop[1].Strings()
	if chainids == nil || eids == nil {
		return ef("Could not find chain attributes 'struct_asym.id' and " +
			"'struct_asym.entity_id'.")
	}
	numModelsTag := "pdbx_nmr_ensemble.conformers_submitted_total_number"
	numModels := value(b, numModelsTag).Int()
	if numModels == 0 {
		numModels = 1
	}
	for i := range chainids {
		ent := e.Entities[eids[i][0]]
		chain := &Chain{Entity: ent, Id: chainids[i][0]}
		ent.Chains[chain.Id] = chain
	}
	return nil
}

func (e *Entry) readAtomSites(b *cif.DataBlock) error {
	loop := asLoop(b, "atom_site.group_pdb", "atom_site.label_atom_id",
		"atom_site.label_asym_id", "atom_site.label_entity_id",
		"atom_site.label_seq_id", "atom_site.cartn_x", "atom_site.cartn_y",
		"atom_site.cartn_z", "atom_site.pdbx_pdb_model_num",
		"atom_site.label_comp_id")
	groups, atoms := loop[0].Strings(), loop[1].Strings()
	chainids, eids := loop[2].Strings(), loop[3].Strings()
	seqids, modelids := loop[4].Ints(), loop[8].Ints()
	xs, ys, zs := loop[5].Floats(), loop[6].Floats(), loop[7].Floats()
	comps := loop[9].Strings()
	if groups == nil || atoms == nil || chainids == nil || eids == nil ||
		seqids == nil || modelids == nil || comps == nil ||
		xs == nil || ys == nil || zs == nil {
		return ef("The given PDBx/mmCIF data has no ATOM/HETATM records.")
	}

	var curEntity *Entity
	var curChain *Chain
	var curModel *Model
	var sitei int
	for i := range groups {
		eid, cid := eids[i][0], chainids[i][0]
		mid, sid := modelids[i]-1, seqids[i]-1

		curEntity = e.Entities[eid]
		curChain = curEntity.Chains[cid]
		if mid == len(curChain.Models) {
			curChain.Models = append(curChain.Models, &Model{
				Chain:           curChain,
				Num:             mid + 1,
				Sites:           make([]Site, 0, 100),
				AlphaCarbonsSeq: make([]*structure.Coords, len(curEntity.Seq)),
				AlphaCarbons:    make([]structure.Coords, 0, 100),
			})
		}
		curModel = curChain.Models[mid]

		sitei = len(curModel.Sites) - 1
		if len(curModel.Sites) == 0 || curModel.Sites[sitei].SeqIndex != sid {
			curModel.Sites = append(curModel.Sites, Site{
				SeqIndex: sid,
				Comp:     comps[i],
				Atoms:    make([]Atom, 0, 5),
			})
			sitei++
		}
		atom := Atom{
			Name:   atoms[i],
			Het:    groups[i] == "HETATM",
			Coords: structure.Coords{X: xs[i], Y: ys[i], Z: zs[i]},
		}
		curModel.Sites[sitei].Atoms = append(curModel.Sites[sitei].Atoms, atom)
		if atoms[i] == "CA" {
			if sid >= 0 {
				curModel.AlphaCarbonsSeq[sid] = &atom.Coords
			}
			curModel.AlphaCarbons = append(curModel.AlphaCarbons, atom.Coords)
		}
	}
	return nil
}

// value returns the data value tagged by "key". If it does not exist, then
// an empty string is returned (wrapped in a cif.Value).
func value(b *cif.DataBlock, key string) cif.Value {
	if v, ok := b.Items[key]; ok {
		return v
	}
	return cif.AsValue("")
}

// asLoop retrieves the Loop containing the data tag "key". If a loop does
// not exist, then one is created with a single row with columns corresponding
// to "key" and each of the tags in "others". If the tag in "key" or any
// tag in "others" does not exist, an empty string is used for its value.
//
// The purpose of this function is to abstract over whether some data set in
// a PDBx/CIF file is represented as a loop or not. For example, if a PDBx file
// has only one entity, then the "entity.*" tags are not in a loop. But if there
// is more than one entity, they are declared as a loop.
func asLoop(b *cif.DataBlock, key string, others ...string) []cif.ValueLoop {
	tags := append([]string{key}, others...)
	asColumns := func(loop *cif.Loop) []cif.ValueLoop {
		vloop := make([]cif.ValueLoop, len(tags))
		for i, tag := range tags {
			vloop[i] = loop.Get(tag)
		}
		return vloop
	}

	if loop, ok := b.Loops[key]; ok {
		return asColumns(loop)
	}
	loop := &cif.Loop{
		Columns: make(map[string]int, len(tags)),
		Values:  make([]cif.ValueLoop, len(tags)),
	}
	for i, tag := range tags {
		loop.Columns[tag] = i
		switch v := value(b, tag).Raw().(type) {
		case string:
			loop.Values[i] = cif.AsValues([]string{v})
		case int:
			loop.Values[i] = cif.AsValues([]int{v})
		case float64:
			loop.Values[i] = cif.AsValues([]float64{v})
		default:
			panic(sf("Unknown value type %T for %s.", v, tag))
		}
	}
	return asColumns(loop)
}
