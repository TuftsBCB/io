package pdbx

import (
	"bytes"
	"compress/gzip"
	"flag"
	"io"
	"io/ioutil"
	"log"
	"os"
	"strings"
	"testing"

	"github.com/TuftsBCB/io/pdb"
)

var (
	flagCifFile = ""
	flagPdbFile = ""
)

func init() {
	flag.StringVar(&flagCifFile, "cif-file", flagCifFile,
		"When set, tests using a PDBx file will be run with the file given.")
	flag.StringVar(&flagPdbFile, "pdb-file", flagPdbFile,
		"When set, tests using a PDB file will be run with the file given.")
	flag.Parse()
}

func TestRead(t *testing.T) {
	if len(flagCifFile) == 0 {
		return
	}
	e := openCif(flagCifFile)
	lf("Entry: %s\n", e.Id)
	lf("Title: %s\n", e.Title)
	for _, ent := range e.Entities {
		lf("  Entity: %c, %s, %f, %f\n", ent.Id, ent.Type, ent.FormulaWeight,
			ent.NumOfMolecules)
		lf("    %s\n", ent.Seq)
		for _, chain := range ent.Chains {
			lf("    Chain: %c (models: %d)\n", chain.Id, len(chain.Models))
			lf("      Atom:    %#v\n", chain.Models[0].Sites[0].Atoms[0])
			lf("      SeqAtom: %#v\n", chain.Models[0].AlphaCarbons[13])
		}
		lf("\n")
	}
}

func BenchmarkReadCIF(b *testing.B) {
	if len(flagCifFile) == 0 {
		return
	}

	r := gzReader(flagCifFile)
	bs, err := ioutil.ReadAll(r)
	if err != nil {
		log.Fatal(err)
	}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, err := Read(bytes.NewReader(bs)); err != nil {
			log.Fatal(err)
		}
	}
}

func BenchmarkReadPDB(b *testing.B) {
	if len(flagPdbFile) == 0 {
		return
	}

	r := gzReader(flagPdbFile)
	bs, err := ioutil.ReadAll(r)
	if err != nil {
		log.Fatal(err)
	}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if e, err := pdb.Read(bytes.NewReader(bs), flagPdbFile); err != nil {
			log.Fatal(err)
		} else {
			for _, chain := range e.Chains {
				chain.SequenceCaAtoms()
			}
		}
	}
}

// openCif opens a PDBx/mmCIF file, accounting for gzip compression.
// If there is an error, the current test fails.
func openCif(fp string) *Entry {
	e, err := Read(gzReader(fp))
	if err != nil {
		log.Fatal(err)
	}
	return e
}

// gzReader returns a reader for the given file.
// If there is a ".gz" extension, the reader is wrapped in a gzip reader.
// If there was a problem opening the file, the current test/benchmark fails.
func gzReader(fp string) io.Reader {
	var f io.Reader
	f, err := os.Open(fp)
	if err != nil {
		log.Fatal(err)
	}

	if strings.HasSuffix(fp, ".gz") {
		f, err = gzip.NewReader(f)
		if err != nil {
			log.Fatal(err)
		}
	}
	return f
}
