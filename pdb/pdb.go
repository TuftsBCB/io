package pdb

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"path"
	"strconv"
	"strings"

	"github.com/BurntSushi/bcbgo/apps/matt"
)

// AminoThreeToOne is a map from three letter amino acids to their
// corresponding single letter representation.
var AminoThreeToOne = map[string]byte{
	"ALA": 'A', "ARG": 'R', "ASN": 'N', "ASP": 'D', "CYS": 'C',
	"GLU": 'E', "GLN": 'Q', "GLY": 'G', "HIS": 'H', "ILE": 'I',
	"LEU": 'L', "LYS": 'K', "MET": 'M', "PHE": 'F', "PRO": 'P',
	"SER": 'S', "THR": 'T', "TRP": 'W', "TYR": 'Y', "VAL": 'V',
	"SEC": 'U', "PYL": 'O',
	"UNK": 'X', "ACE": 'X', "NH2": 'X',
	"ASX": 'X', "GLX": 'X',
}

// AminoOneToThree is the reverse of AminoThreeToOne. It is created in
// this packages 'init' function.
var AminoOneToThree = map[byte]string{}

func init() {
	// Create a reverse map of AminoThreeToOne.
	for k, v := range AminoThreeToOne {
		AminoOneToThree[v] = k
	}
}

// Entry represents all information known about a particular PDB file (that
// has been implemented in this package).
//
// Currently, a PDB entry is simply a file path and a map of protein chains.
type Entry struct {
	Path           string
	IdCode         string
	Classification string
	Chains         []*Chain
}

// New creates a new PDB Entry from a file. If the file cannot be read, or there
// is an error parsing the PDB file, an error is returned.
//
// If the file name ends with ".gz", gzip decompression will be used.
func New(fileName string) (*Entry, error) {
	var reader io.Reader
	var err error

	reader, err = os.Open(fileName)
	if err != nil {
		return nil, err
	}

	// If the file is gzipped, use the gzip decompressor.
	if path.Ext(fileName) == ".gz" {
		reader, err = gzip.NewReader(reader)
		if err != nil {
			return nil, err
		}
	}

	entry := &Entry{
		Path:   fileName,
		Chains: make([]*Chain, 0),
	}

	// Now traverse each line, and process it according to the record name.
	// Note that it is imperative that we preserve the order of ATOM records
	// as we read them. We are currently trying to replicate Fragbag, and this
	// is what Fragbag does. (A more stable approach would require more
	// information from the PDB file; like differentiating models, since
	// sorting on ATOM serial number isn't good enough.)
	breader := bufio.NewReaderSize(reader, 1000)
	for {
		// We ignore 'isPrefix' here, since we never care about lines longer
		// than 1000 characters, which is the size of our buffer.
		line, _, err := breader.ReadLine()
		if err == io.EOF {
			break
		} else if err != nil {
			return nil, err
		}

		// The record name is always in the fix six columns.
		switch strings.TrimSpace(string(line[0:6])) {
		case "HEADER":
			if err := entry.parseHeader(line); err != nil {
				return nil, err
			}
		case "SEQRES":
			entry.parseSeqres(line)
		case "ATOM":
			entry.parseAtom(line)
		}
	}

	// If we didn't pick up any chains, this probably isn't a valid PDB file.
	if len(entry.Chains) == 0 {
		return nil, fmt.Errorf("The file '%s' does not appear to be a valid "+
			"PDB file.", fileName)
	}

	// If we couldn't find an Id code, inspect the base name of the file path.
	if len(entry.IdCode) == 0 {
		name := path.Base(fileName)
		switch {
		case len(name) >= 7 && name[0:3] == "pdb":
			entry.IdCode = name[3:7]
		case len(name) == 7: // cath
			entry.IdCode = name[0:4]
		}
	}

	return entry, nil
}

// Chain looks for the chain with identifier ident and returns it. 'nil' is
// returned if the chain could not be found.
func (e *Entry) Chain(ident byte) *Chain {
	for _, chain := range e.Chains {
		if chain.Ident == ident {
			return chain
		}
	}
	return nil
}

// OneChain returns a single chain in the PDB file. If there is more than one
// chain, OneChain will panic. This is convenient when you expect a PDB file to
// have only a single chain, but don't know the name.
func (e *Entry) OneChain() *Chain {
	if len(e.Chains) != 1 {
		panic(fmt.Sprintf("OneChain can only be called on PDB entries with "+
			"ONE chain. But the '%s' PDB entry has %d chains.",
			e.Path, len(e.Chains)))
	}
	return e.Chains[0]
}

// PDBArg is a convenience method for creating a PDBArg that can be used in
// the 'matt' package. It sets 'Location' in PDBArg to 'Path' from 'Entry'.
func (e *Entry) PDBArg() matt.PDBArg {
	return matt.PDBArg{Location: e.Path}
}

// Name returns the base name of the path of this PDB entry.
func (e *Entry) Name() string {
	return path.Base(e.Path)
}

// String returns a list of all chains, their residue start/stop indices,
// and the amino acid sequence.
func (e *Entry) String() string {
	lines := make([]string, 0)
	for _, chain := range e.Chains {
		lines = append(lines, chain.String())
	}
	return strings.Join(lines, "\n")
}

// getOrMakeChain looks for a chain in the 'Chains' slice corresponding to the
// chain indentifier. If one exists, it is returned. If one doesn't exist,
// it is created, memory is allocated and it is returned.
func (e *Entry) getOrMakeChain(ident byte) *Chain {
	if ident == ' ' {
		ident = '_'
	}

	chain := e.Chain(ident)
	if chain != nil {
		return chain
	}
	newChain := &Chain{
		Entry:            e,
		Ident:            ident,
		Sequence:         make([]byte, 0, 10),
		AtomResidueStart: 0,
		AtomResidueEnd:   0,
		CaAtoms:          make(Atoms, 0, 30),
	}
	e.Chains = append(e.Chains, newChain)
	return newChain
}

// parseHeader loads the "idCode" and "classification" fields from the
// header record.
//
// If the fields are already filled, then we've seen a second header record
// and therefore report an error.
func (e *Entry) parseHeader(line []byte) error {
	if len(e.Classification) > 0 || len(e.IdCode) > 0 {
		return fmt.Errorf("More than one HEADER record was found.")
	}
	e.Classification = strings.TrimSpace(string(line[10:50]))
	e.IdCode = strings.TrimSpace(string(line[62:66]))
	return nil
}

// parseSeqres loads all pertinent information from SEQRES records in a PDB
// file. In particular, amino acid resides are read and added to the chain's
// "Sequence" field. If a residue isn't a valid amino acid, it is simply
// ignored.
//
// N.B. This assumes that the SEQRES records are in order in the PDB file.
func (e *Entry) parseSeqres(line []byte) {
	chain := e.getOrMakeChain(line[11])

	// Residues are in columns 19-21, 23-25, 27-29, ..., 67-69
	for i := 19; i <= 67; i += 4 {
		end := i + 3

		// If we're passed the end of this line, quit.
		if end >= len(line) {
			break
		}

		// Get the residue. If it's not in our sequence map, skip it.
		residue := strings.TrimSpace(string(line[i:end]))
		if single, ok := AminoThreeToOne[residue]; ok {
			chain.Sequence = append(chain.Sequence, single)
		}
	}
}

// parseAtom loads all pertinent information from ATOM records in a PDB file.
// Currently, this only includes deducing the amino acid residue start and
// stop indices. (Note that the length of the range is not necessarily
// equivalent to the length of the amino acid sequence found in the SEQRES
// records.)
//
// ATOM records without a valid amino acid residue in columns 18-20 are ignored.
func (e *Entry) parseAtom(line []byte) {
	chain := e.getOrMakeChain(line[21])

	// An ATOM record is only processed if it corresponds to an amino acid
	// residue. (Which is in columns 17-19.)
	residue := strings.TrimSpace(string(line[17:20]))
	if _, ok := AminoThreeToOne[residue]; !ok {
		// Sanity check. I'm pretty sure that only amino acids have three
		// letter abbreviations.
		if len(residue) == 3 {
			panic(fmt.Sprintf("The residue '%s' found in PDB file '%s' has "+
				"length 3, but is not in my amino acid map.",
				residue, e.Path))
		}
		return
	}

	// The residue sequence number is in columns 22-25. Grab it, trim it,
	// and look for an integer.
	snum := strings.TrimSpace(string(line[22:26]))
	inum := int(0)
	if num, err := strconv.ParseInt(snum, 10, 32); err == nil {
		inum = int(num)
		switch {
		case chain.AtomResidueStart == 0 || inum < chain.AtomResidueStart:
			chain.AtomResidueStart = inum
		case chain.AtomResidueEnd == 0 || inum > chain.AtomResidueEnd:
			chain.AtomResidueEnd = inum
		}
	}

	// Build an Atom value. We need the serial number from columns 6-10,
	// the atom name from columns 12-15, the amino acid residue from
	// columns 17-19 (we already have that: 'residue'), the residue sequence
	// number from columns 22-25 (already have that too: 'inum'), and the
	// three dimension coordinates in columns 30-37 (x), 38-45 (y), and
	// 46-53 (z).
	atom := Atom{
		Name:       strings.TrimSpace(string(line[12:16])),
		Residue:    residue,
		ResidueInd: inum,
		Coords:     [3]float64{},
	}

	serialStr := strings.TrimSpace(string(line[6:11]))
	if serial64, err := strconv.ParseInt(serialStr, 10, 32); err == nil {
		atom.Serial = int(serial64)
	}

	xstr := strings.TrimSpace(string(line[30:38]))
	ystr := strings.TrimSpace(string(line[38:46]))
	zstr := strings.TrimSpace(string(line[46:54]))
	if x64, err := strconv.ParseFloat(xstr, 64); err == nil {
		atom.Coords[0] = x64
	}
	if y64, err := strconv.ParseFloat(ystr, 64); err == nil {
		atom.Coords[1] = y64
	}
	if z64, err := strconv.ParseFloat(zstr, 64); err == nil {
		atom.Coords[2] = z64
	}

	// Now add our atom to the chain.
	chain.Atoms = append(chain.Atoms, atom)
	if atom.Name == "CA" {
		chain.CaAtoms = append(chain.CaAtoms, atom)
	}
}

// Chain represents a protein chain or subunit in a PDB file. Each chain has
// its own identifier, amino acid sequence (if its a protein sequence), and
// the start and stop residue indices of the ATOM coordinates.
//
// It also contains a slice of all carbon-alpha ATOM records corresponding
// to an amino acid.
type Chain struct {
	Entry                            *Entry
	Ident                            byte
	Sequence                         []byte
	AtomResidueStart, AtomResidueEnd int
	Atoms                            Atoms
	CaAtoms                          Atoms
}

// PDBArg is a convenience method for creating a PDBArg for a specific
// chain that can be used in the 'matt' package.
func (c *Chain) PDBArg() matt.PDBArg {
	return matt.PDBArg{
		Location: c.Entry.Path,
		Chain:    c.Ident,
	}
}

// ValidProtein returns true when there are ATOM records corresponding to
// a protein backbone.
func (c *Chain) ValidProtein() bool {
	return c.AtomResidueStart != c.AtomResidueEnd
}

// String returns a FASTA-like formatted string of this chain and all of its
// related information.
func (c *Chain) String() string {
	return strings.TrimSpace(
		fmt.Sprintf("> Chain %c (%d, %d) :: length %d\n%s",
			c.Ident, c.AtomResidueStart, c.AtomResidueEnd,
			len(c.Sequence), string(c.Sequence)))
}

// Atom contains information about an ATOM record, including the serial
// number, the residue (and residue sequence number), the atom name and the
// three dimensional coordinates.
type Atom struct {
	Serial     int
	Name       string
	ResidueInd int
	Residue    string

	// Coords is a triple where the first element is X, the second is Y and
	// the third is Z.
	Coords [3]float64
}

func (a Atom) String() string {
	return fmt.Sprintf("(%d, %s, %d, %s, [%0.4f %0.4f %0.4f])",
		a.Serial, a.Name, a.ResidueInd, a.Residue,
		a.Coords[0], a.Coords[1], a.Coords[2])
}

// Atoms names a slice of Atom for sorting.
type Atoms []Atom

func (as Atoms) String() string {
	lines := make([]string, len(as))
	for i, atom := range as {
		lines[i] = atom.String()
	}
	return strings.Join(lines, "\n")
}
