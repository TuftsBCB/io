package pdb

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"path"
	"strconv"
	"strings"

	"github.com/TuftsBCB/seq"
)

type pdbParser struct {
	entry    *Entry
	curModel int
	line     []byte
	modified map[string]string
	seqres   map[byte][]string

	missing map[byte][]missingResidue
	r465    bool

	processed map[seen]bool
	lastSeen  seen
}

type missingResidue struct {
	residue string
	seqNum  int
	insCode byte
}

type seen struct {
	chain byte
	model int
}

func ReadPDB(fp string) (*Entry, error) {
	var reader io.Reader
	var err error

	f, err := os.Open(fp)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	reader = f

	// If the file is gzipped, use the gzip decompressor.
	if path.Ext(fp) == ".gz" {
		reader, err = gzip.NewReader(reader)
		if err != nil {
			return nil, err
		}
	}

	return Read(reader, fp)
}

func Read(reader io.Reader, fpath string) (*Entry, error) {
	entry := &Entry{
		Path:   fpath,
		Chains: make([]*Chain, 0),
	}

	// Now traverse each line, and process it according to the record name.
	// Note that it is imperative that we preserve the order of ATOM records
	// as we read them. We are currently trying to replicate Fragbag, and this
	// is what Fragbag does. (A more stable approach would require more
	// information from the PDB file; like differentiating models, since
	// sorting on ATOM serial number isn't good enough.)
	breader := bufio.NewReaderSize(reader, 1000)
	parser := pdbParser{
		entry:     entry,
		curModel:  1,
		line:      nil,
		modified:  make(map[string]string),
		seqres:    make(map[byte][]string),
		missing:   make(map[byte][]missingResidue),
		r465:      false,
		processed: make(map[seen]bool),
		lastSeen:  seen{0, 0},
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
	for _, chain := range parser.entry.Chains {
		for _, r := range parser.seqres[chain.Ident] {
			newr, err := parser.residueAbbrev(r)
			if err != nil {
				return nil, fmt.Errorf("PDB (SEQRES) '%s': %s",
					parser.entry.Path, err)
			}
			chain.Sequence = append(chain.Sequence, newr)
		}
	}

	// We also need to translate the missing residues.
	for _, chain := range parser.entry.Chains {
		for _, r := range parser.missing[chain.Ident] {
			abbrev, err := parser.residueAbbrev(r.residue)
			if err != nil {
				return nil, fmt.Errorf("PDB (REMARK 465) '%s': %s",
					parser.entry.Path, err)
			}
			newr := &Residue{
				Name:          abbrev,
				SequenceNum:   r.seqNum,
				InsertionCode: r.insCode,
				Atoms:         nil,
			}
			chain.Missing = append(chain.Missing, newr)
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
		case len(name) >= 6: // SCOP
			entry.scop = stripExt(name)
			entry.IdCode = name[1:5]
		case len(name) == 7: // cath
			entry.cath = stripExt(name)
			entry.IdCode = name[0:4]
		}
	}

	return entry, nil
}

func stripExt(s string) string {
	i := strings.Index(s, ".")
	if i == -1 {
		return s
	}
	return s[0:i]
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
		residue := p.cols(25, 27)
		if len(residue) == 0 {
			residue = "UNK"
		}
		p.modified[p.cols(13, 15)] = residue
	case "HET":
		residue := p.cols(8, 10)

		// Only add this if there isn't a MODRES for this residue.
		if _, ok := p.modified[residue]; !ok {
			p.modified[residue] = "UNK"
		}
	case "ATOM":
		fallthrough
	case "HETATM":
		if err := p.parseAtom(); err != nil {
			return err
		}
	case "TER":
		// Whatever ATOM record we saw last will have its {chain, model}
		// added to the 'processed' map. This effectively disallows any further
		// additions to this particular {chain, model}.
		// This cuts off HETATM's after "TER", but I'm fine with that for now.
		// fmt.Printf("%c %d\n", p.lastSeen.chain, p.lastSeen.model)
		p.processed[p.lastSeen] = true
	case "REMARK":
		num, err := p.atoi(8, 10)
		if err != nil {
			// Some people like to put in their own remarks and not conform
			// to the SPEC. Why............
			return nil
		}
		switch num {
		case 465:
			if err := p.parseRemark465(); err != nil {
				return err
			}
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

func (p *pdbParser) parseAtom() error {
	ident := p.at(22)

	// If we've already processed this {chain, model}, then don't look
	// for more ATOMs.
	if p.processed[seen{ident, p.curModel}] {
		return nil
	}

	res := p.cols(18, 20)
	if res == "HOH" { // ignore water HETATMs
		return nil
	}

	insCode := p.at(27)
	if insCode == ' ' {
		insCode = 0
	}

	seqNum, err := p.atoi(23, 26)
	if err != nil {
		return err
	}

	// If this is a HETATM, we often see weird residues, so we'll suppress an
	// error if we don't recognize the residue. The residue will be set to 'X'.
	flexible := p.cols(1, 6) == "HETATM"
	residue, err := p.getResidue(ident, res, seqNum, insCode, flexible)
	if err != nil {
		return err
	}
	atom := Atom{
		Name:   p.cols(13, 16),
		Het:    p.cols(1, 6) == "HETATM",
		Coords: Coords{},
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

	p.lastSeen = seen{ident, p.curModel}
	residue.Atoms = append(residue.Atoms, atom)
	return nil
}

func (p *pdbParser) parseRemark465() error {
	residue := p.cols(16, 18)

	if residue == "RES" && p.at(20) == 'C' && p.cols(22, 27) == "SSSEQI" {
		p.r465 = true
		return nil
	}
	if !p.r465 {
		return nil
	}
	if p.r465 && (len(residue) == 0 || len(residue) > 3) {
		p.r465 = false
		return nil
	}

	chainIdent := p.at(20)

	seqNum, err := p.atoi(22, 26)
	if err != nil {
		return err
	}

	insCode := p.at(27)
	if insCode == ' ' {
		insCode = 0
	}

	if _, ok := p.missing[chainIdent]; !ok {
		p.missing[chainIdent] = make([]missingResidue, 0)
	}
	missing := missingResidue{residue, seqNum, insCode}
	p.missing[chainIdent] = append(p.missing[chainIdent], missing)
	return nil
}

func (p *pdbParser) getChain(ident byte) *Chain {
	for i, chain := range p.entry.Chains {
		if chain.Ident == ident {
			return p.entry.Chains[i]
		}
	}
	chain := &Chain{
		Entry:    p.entry,
		Ident:    ident,
		SeqType:  -1,
		Sequence: make([]seq.Residue, 0, 25),
		Models:   make([]*Model, 0),
		Missing:  make([]*Residue, 0),
	}
	p.entry.Chains = append(p.entry.Chains, chain)
	return chain
}

func (p pdbParser) getModel(ident byte) *Model {
	chain := p.getChain(ident)
	for i, model := range chain.Models {
		if model.Num == p.curModel {
			return chain.Models[i]
		}
	}

	model := &Model{
		Entry:    p.entry,
		Chain:    chain,
		Num:      p.curModel,
		Residues: make([]*Residue, 0, 25),
	}
	chain.Models = append(chain.Models, model)
	return model
}

func (p pdbParser) getResidue(ident byte,
	res string, seqNum int, insCode byte, flexible bool) (*Residue, error) {

	model := p.getModel(ident)
	for i := range model.Residues {
		if model.Residues[i].SequenceNum == seqNum &&
			model.Residues[i].InsertionCode == insCode {

			return model.Residues[i], nil
		}
	}
	abbrev, err := p.residueAbbrev(res)
	if !flexible && err != nil {
		return nil, fmt.Errorf("PDB (ATOM) '%s': %s", p.entry.Path, err)
	}

	residue := &Residue{
		Name:          abbrev,
		SequenceNum:   seqNum,
		InsertionCode: insCode,
		Atoms:         make([]Atom, 0, 2),
	}
	model.Residues = append(model.Residues, residue)
	return residue, nil
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

func (p pdbParser) residueAbbrev(long string) (seq.Residue, error) {
	if res, ok := p.modified[long]; ok {
		return getAbbrev(res)
	}
	abbrev, err := getAbbrev(long)

	// If we couldn't find a residue, check the missing residues. If this
	// residue shows up there, let's just ignore it and return 'X'.
	// If it isn't in the missing residues, then quit with an error.
	if err != nil {
		for _, missing := range p.missing {
			for _, residue := range missing {
				if residue.residue == long {
					return 'X', nil
				}
			}
		}
		return 'X', err
	}
	return abbrev, nil
}
