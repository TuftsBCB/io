package pdb

import (
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/BurntSushi/bcbgo/io/fasta"
	"github.com/BurntSushi/bcbgo/seq"
)

// look at:
// 1aco
// 3cfs

var flagPdbDir = ""

func init() {
	log.SetFlags(0)

	flag.StringVar(&flagPdbDir, "pdb-dir", flagPdbDir,
		"Directory to use to find PDB files (recursively).")
	flag.Parse()
}

func TestReadPDB(t *testing.T) {
	entry := readPDB()
	if len(entry.Chains) != 2 {
		t.Fatalf("Expected 2 chains but got %d.", len(entry.Chains))
	}
}

// This tests the correspondence between residues in SEQRES records with
// residues in ATOM records. In particular, the correspondence should match
// residues classified as disordered by the PDB in the "pdb_disorder.fasta"
// file, which contains all of the "disordered" sequences from:
// http://www.pdb.org/pdb/files/ss_dis.txt
//
// A disordered sequence only contains two residues: '-' for a residue with
// an ATOM record and 'X' for a residue without an ATOM record.
//
// N.B. The "pdb_disorder.fasta" file isn't included in the repository because
// it's >50MB. If you need it, contact me or generate it from ss_dis.txt.
func TestSequenceCa(t *testing.T) {
	if len(flagPdbDir) == 0 {
		return
	}

	disordered := readDisordered()
	if disordered == nil {
		// Don't fail the test if we can't get the file.
		return
	}

	// So I think some of the disordered classifications made by the PDB
	// in pdb_disorder.fasta are wrong. I suspect they generated their
	// correspondence by actual sequence alignment, rather than looking at
	// the REMARK 465 residues.
	//
	// The following lists the chains I believe my method gets right that
	// the PDB gets wrong. They were manually added by burntsushi.
	//
	// Thus, we don't test these.
	pdbWrong := []string{
		// I think they are just off by one here.
		"2CFAB", "3ANZA", "4APCA", "3ARCC", "3ARCD", "3ARCF", "3ARCJ", "2B66E",
		"1B8MB", "2BFEA", "2BFLB",

		// The PDB says these are wholly disordered, by the actual file gives
		// ATOM records for many of the residues. So...
		"2A1DA", "2A1DB", "2A1DD", "2A1DE", "2A1DF", "2A1DH", "2AI8B", "3AI8B",
		"3AI8A",
	}
	pdbWrongMap := make(map[string]bool)
	for _, idcode := range pdbWrong {
		pdbWrongMap[idcode] = true
	}

	// These are corrupt PDB entries. Meaning that they don't fulfill promises
	// that they state they will fill.
	// i.e., Not including a missing residue in REMARK 465.
	corrupt := []string{
		"4A18A", "4A18C", "4A18L", "4A19A", "4A19C", "4A19Q", "4A1BA", "4A1BC",
		"4A1DA", "4A1DC", "2AAZA", "2AAZB", "2AAZC", "2AAZD", "2AAZE", "2AAZF",
		"2AAZG", "2AAZH", "2AAZI", "2AAZJ", "2AAZK", "2AAZL", "2AAZM", "2AAZN",
		"2AAZO", "2AAZP", "2BEQA", "2BEQB", "2BEQC", "2BEQD", "2BEQE", "2BEQF",
		"2BETA", "2BETB", "2BETC", "2BETD", "2BETE",

		// order of residues in REMARK 465 does not match SEQRES.
		"2C24B",

		// These are interesting. Tons of residues are marked
		// UNKNOWN, and it's therefore difficult to get a good alignment.
		"4A1SA", "4A1SB", "2AXTX", "3BZ1Y", "3BZ2Y",
	}
	corruptMap := make(map[string]bool)
	for _, idcode := range corrupt {
		corruptMap[idcode] = true
	}

	walkDir(flagPdbDir, func(path string) error {
		entry, err := ReadPDB(path)
		if err != nil {
			t.Fatalf("%s\n", err)
		}

		for _, chain := range entry.Chains {
			key := chainKey(chain)
			if pdbWrongMap[key] || corruptMap[key] {
				continue
			}
			testDisordered(t, chain, disordered)
		}
		return nil
	})
}

// Use 'go test --pdb-dir blah'.
func TestManyReadPDB(t *testing.T) {
	if len(flagPdbDir) == 0 {
		return
	}

	walkDir(flagPdbDir, func(path string) error {
		_, err := ReadPDB(path)
		if err != nil {
			t.Fatalf("%s\n", err)
		}
		// for _, chain := range entry.Chains {
		// chain.SequenceCaAtoms()
		// }
		return nil
	})
}

func ExampleSeqresCas() {
	entry := readPDB()
	cas, err := entry.Chains[0].SequenceCaAtoms()
	assert(err)
	if ca := cas[12]; ca != nil {
		fmt.Println(ca)
	}

	// Output:
	// -13.956 15.047 -5.126
}

func ExampleMissingResidues() {
	entry := readPDB()
	last := entry.Chain('A').Missing[12]
	fmt.Printf("%c %d\n", last.Name, last.SequenceNum)

	// Output:
	// S 211
}

func ExampleReadPDB() {
	entry := readPDB()
	fmt.Printf("%s\n", entry.Chains[1].Sequence)

	res := entry.Chains[0].Models[1].Residues[0]
	atom := res.Atoms[1]
	fmt.Printf("%s %c %0.3f %0.3f %0.3f\n",
		atom.Name, res.Name, atom.X, atom.Y, atom.Z)

	// Output:
	// AYIGPYL
	// CA S -18.866 9.770 -5.303
}

func BenchmarkReadPDB(b *testing.B) {
	for i := 0; i < b.N; i++ {
		readPDB()
	}
}

// chain is a chain from a PDB file, and disordered is a map of all chains
// from the PDB to sequences that specify disordered residues.
// (See TestSequenceCa.)
func testDisordered(t *testing.T,
	chain *Chain, disordered map[string]seq.Sequence) {

	key := chainKey(chain)
	answerSeq, ok := disordered[key]
	if !ok {
		log.Printf("Could not find answer for '%s'.", key)
		return
	}
	answer := answerSeq.Residues

	mapping, err := chain.SequenceAtoms()
	if err != nil {
		log.Printf("%s\n", err)
		return
	}

	guess := make([]seq.Residue, len(mapping))
	for i, residue := range mapping {
		if residue == nil {
			guess[i] = 'X'
		} else {
			guess[i] = '-'
		}
	}

	wrongo := func() {
		log.Println(key)
		log.Println("Guess:")
		log.Printf("%s\n\n", guess)
		log.Println("Answer:")
		log.Printf("%s\n\n", answer)
	}

	if len(guess) != len(answer) {
		wrongo()
		t.Fatalf("Lengths do not match. %d != %d", len(guess), len(answer))
	}
	for i := 0; i < len(guess); i++ {
		if guess[i] != answer[i] {
			wrongo()
			t.Fatalf("Residues at position %d are different. '%c' != '%c'.",
				i, guess[i], answer[i])
		}
	}
}

func readPDB() *Entry {
	entry, err := ReadPDB("pdb2nmb.ent.gz")
	assert(err)
	return entry
}

// This file isn't in the repo, so if it isn't there, we just return nil.
// The caller must not fail with nil returned.
func readDisordered() map[string]seq.Sequence {
	f, err := os.Open("pdb_disordered.fasta")
	if err != nil {
		return nil
	}

	r := fasta.NewReader(f)
	r.TrustSequences = true
	m := make(map[string]seq.Sequence, 10000)
	for {
		s, err := r.Read()
		if err == io.EOF {
			break
		}
		assert(err)

		m[strings.ToUpper(s.Name)] = s
	}
	return m
}

func chainKey(chain *Chain) string {
	key := fmt.Sprintf("%s%c", chain.Entry.IdCode, chain.Ident)
	return strings.ToUpper(key)
}

func walkDir(dir string, f func(path string) error) {
	walk := func(path string, info os.FileInfo, err error) error {
		if err != nil || info.IsDir() {
			return nil
		}

		suffix := func(s string) bool { return strings.HasSuffix(path, s) }
		if !suffix(".pdb") && !suffix(".ent.gz") {
			return nil
		}

		return f(path)
	}
	filepath.Walk(dir, walk)
}

func assert(err error) {
	if err != nil {
		log.Fatalln(err)
	}
}
