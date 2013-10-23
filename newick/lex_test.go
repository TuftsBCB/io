package newick

import (
	"bytes"
	"fmt"
	"io"
	"os"
	"testing"
)

func ef(format string, v ...interface{}) {
	fmt.Fprintf(os.Stderr, format, v...)
}

func sample(s string) io.Reader {
	return bytes.NewReader([]byte(s))
}

func TestLexer(t *testing.T) {
	// v := sample("(A,B,(C,D)E)F;")
	// v := sample("(,,(,));")
	// v := sample("(:0.1,:0.2,(:0.3,:0.4):0.5);")
	// v := sample("(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);")
	// v := sample("((d1qbea_:0.597492,d1dwna_:0.632208):0.162939," +
	// "(d1gav0_:0.526213,(d1unaa_:0.457107,d2iznb1:0.523093):0.043387);")
	// v := sample("(B)A;")
	v := sample("((X,Y)C)ROOT;")

	// v, _ := os.Open("/media/Hummingbird/data/bio/mattbench/astral_true.tree")
	lx := lex(v)
	for {
		item := lx.nextItem()
		if item.typ == itemEOF {
			break
		} else if item.typ == itemError {
			t.Fatal(item.val)
		}
		// ef("%s\n", item)
	}
}
