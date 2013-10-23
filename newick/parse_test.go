package newick

import (
	"testing"
)

func TestParser(t *testing.T) {
	v := sample("(A,B,(X,Y)C)ROOT;(A,B,C)ROOT;")
	// v, _ := os.Open("/media/Hummingbird/data/bio/mattbench/astral_true.tree")

	r := NewReader(v)
	trees, err := r.ReadAll()
	if err != nil {
		t.Fatal(err)
	}

	for i := range trees {
		ef("Tree %d:\n", i)
		ef("%s", trees[i])
		ef("-------------------------\n")
	}
}
