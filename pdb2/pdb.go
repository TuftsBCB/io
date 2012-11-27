package pdb2

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"path"
	"strconv"

	"github.com/BurntSushi/bcbgo/seq"
)

type pdbParser struct {
	entry    *Entry
	curModel int
	line     []byte
	modified map[modification]string
	seqres   map[byte][]string
}

type modification struct {
	chainIdent byte
	from       string
}

func ReadPDB(fp string) (*Entry, error) {
	var reader io.Reader
	var err error

	reader, err = os.Open(fp)
	if err != nil {
		return nil, err
	}

	// If the file is gzipped, use the gzip decompressor.
	if path.Ext(fp) == ".gz" {
		reader, err = gzip.NewReader(reader)
		if err != nil {
			return nil, err
		}
	}

	entry := &Entry{
		Path:   fp,
		Chains: make([]Chain, 0),
	}

	// Now traverse each line, and process it according to the record name.
	// Note that it is imperative that we preserve the order of ATOM records
	// as we read them. We are currently trying to replicate Fragbag, and this
	// is what Fragbag does. (A more stable approach would require more
	// information from the PDB file; like differentiating models, since
	// sorting on ATOM serial number isn't good enough.)
	breader := bufio.NewReaderSize(reader, 1000)
	parser := pdbParser{
		entry:    entry,
		curModel: 1,
		line:     nil,
		modified: make(map[modification]string, 5),
		seqres:   make(map[byte][]string, 2),
	}
	for {
		// We ignore 'isPrefix' here, since we never care about lines longer
		// than 1000 characters, which is the size of our buffer.
		line, _, err := breader.ReadLine()
		if err == io.EOF && len(line) == 0 {
			break
		} else if err != io.EOF && err != nil {
			return nil, err
		}
		parser.line = line
		if err := parser.parseLine(); err != nil {
			return nil, err
		}
	}

	// We've got to back and translate SEQRES residues to their single-letter
	// abbreviations.
	// We didn't do this before, because the MODRES records typically come
	// afterwards. Ug.
	for i := range parser.entry.Chains {
		chain := &parser.entry.Chains[i]
		for _, r := range parser.seqres[chain.Ident] {
			newr := parser.residueAbbrev(chain.Ident, r)
			chain.Sequence = append(chain.Sequence, newr)
		}
	}

	// If we didn't pick up any chains, this probably isn't a valid PDB file.
	if len(entry.Chains) == 0 {
		return nil, fmt.Errorf("The file '%s' does not appear to be a valid "+
			"PDB file.", entry.Path)
	}

	// If we couldn't find an Id code, inspect the base name of the file path.
	if len(entry.IdCode) == 0 {
		name := path.Base(entry.Path)
		switch {
		case len(name) >= 7 && name[0:3] == "pdb":
			entry.IdCode = name[3:7]
		case len(name) == 7: // cath
			entry.IdCode = name[0:4]
		}
	}

	return entry, nil
}

func (p *pdbParser) parseLine() error {
	var err error

	switch p.cols(1, 6) {
	case "HEADER":
		p.entry.IdCode = p.cols(63, 66)
	case "MODEL":
		p.curModel, err = p.atoi(11, 14)
		if err != nil {
			return err
		}
	case "SEQRES":
		p.parseSeqres()
	case "MODRES":
		mod := modification{p.at(17), p.cols(13, 15)}
		p.modified[mod] = p.cols(25, 27)
	case "ATOM":
		if err := p.parseAtom(); err != nil {
			return err
		}
	}
	return nil
}

func (p pdbParser) parseSeqres() {
	chain := p.getChain(p.at(12))
	if _, ok := p.seqres[chain.Ident]; !ok {
		p.seqres[chain.Ident] = make([]string, 0, 25)
	}
	for c := 20; c <= 68; c += 4 {
		res := p.cols(c, c+2)
		if len(res) == 0 {
			break
		}
		if chain.SeqType == -1 {
			// Label the sequence type if we haven't yet.
			// This basically just looks at the length of the unabbreviated
			// residue. len of 3 = protein, len of 2 = deoxy, len of 1 = ribo.
			chain.SeqType = getAbbrevType(res)
		}
		p.seqres[chain.Ident] = append(p.seqres[chain.Ident], res)
	}
}

func (p pdbParser) parseAtom() error {
	ident := p.at(22)
	res := p.cols(18, 20)
	seqNum, err := p.atoi(23, 26)
	if err != nil {
		return err
	}

	residue := p.getResidue(ident, res, seqNum)
	atom := Atom{
		Name: p.cols(13, 16),
	}

	atom.X, err = p.atof(31, 38)
	if err != nil {
		return err
	}
	atom.Y, err = p.atof(39, 46)
	if err != nil {
		return err
	}
	atom.Z, err = p.atof(47, 54)
	if err != nil {
		return err
	}

	residue.Atoms = append(residue.Atoms, atom)
	return nil
}

func (p pdbParser) getChain(ident byte) *Chain {
	for i, chain := range p.entry.Chains {
		if chain.Ident == ident {
			return &p.entry.Chains[i]
		}
	}
	chain := Chain{
		Entry:    p.entry,
		Ident:    ident,
		SeqType:  -1,
		Sequence: make([]seq.Residue, 0, 25),
		Models:   make([]Model, 0, 1),
	}
	p.entry.Chains = append(p.entry.Chains, chain)
	return &p.entry.Chains[len(p.entry.Chains)-1]
}

func (p pdbParser) getModel(ident byte) *Model {
	chain := p.getChain(ident)
	for i, model := range chain.Models {
		if model.Num == p.curModel {
			return &chain.Models[i]
		}
	}
	model := Model{
		Entry:    p.entry,
		Chain:    chain,
		Num:      p.curModel,
		Residues: make([]Residue, 0, 25),
	}
	chain.Models = append(chain.Models, model)
	return &chain.Models[len(chain.Models)-1]
}

func (p pdbParser) getResidue(ident byte, res string, seqNum int) *Residue {
	model := p.getModel(ident)
	for i, residue := range model.Residues {
		if residue.SequenceNum == seqNum {
			return &model.Residues[i]
		}
	}
	residue := Residue{
		Name:         p.residueAbbrev(model.Chain.Ident, res),
		SequenceNum:  seqNum,
		ResidueIndex: -1,
		Atoms:        make([]Atom, 0, 4),
	}
	model.Residues = append(model.Residues, residue)
	return &model.Residues[len(model.Residues)-1]
}

func (p pdbParser) atoi(start, end int) (int, error) {
	return strconv.Atoi(p.cols(start, end))
}

func (p pdbParser) atof(start, end int) (float64, error) {
	return strconv.ParseFloat(p.cols(start, end), 64)
}

func (p pdbParser) cols(start, end int) string {
	rs, re := start-1, end
	if rs >= len(p.line) || rs < 0 {
		return ""
	}
	if re > len(p.line) || re < 0 || re < rs {
		return ""
	}
	return string(bytes.TrimSpace(p.line[rs:re]))
}

func (p pdbParser) at(column int) byte {
	i := column - 1
	if i < 0 || i >= len(p.line) {
		return 0
	}
	return p.line[i]
}

func (p pdbParser) residueAbbrev(chainIdent byte, long string) seq.Residue {
	if res, ok := p.modified[modification{chainIdent, long}]; ok {
		return getAbbrev(res)
	}
	return getAbbrev(long)
}
