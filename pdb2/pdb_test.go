package pdb2

import (
	"fmt"
	"log"
	"testing"
)

func TestReadPDB(t *testing.T) {
	entry := readPDB()
	if len(entry.Chains) != 2 {
		t.Fatalf("Expected 2 chains but got %d.", len(entry.Chains))
	}
}

func ExampleSeqresCas() {
	entry := readPDB()
	for _, ca := range entry.Chains[0].SequenceCaAtoms() {
		fmt.Printf("%v\n", ca)
	}

	// Output:
	// Bah
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

func readPDB() *Entry {
	entry, err := ReadPDB("pdb2nmb.ent.gz")
	assert(err)
	return entry
}

func assert(err error) {
	if err != nil {
		log.Fatalln(err)
	}
}
