package msa

import (
	"bytes"
	"fmt"
	"testing"

	"github.com/TuftsBCB/seq"
)

func TestFasta(t *testing.T) {
	tests := []*bytes.Buffer{
		makeBuffer(inpAlignFasta),
	}
	answers := [][]string{
		alignFasta,
	}
	for i := 0; i < len(tests); i++ {
		computed, err := Read(tests[i])
		if err != nil {
			t.Fatalf("%s", err)
		}
		answer := makeMSA(makeSeqs(answers[i]))
		testEqualAlign(t, computed, answer)
	}
}

func TestA2M(t *testing.T) {
	tests := []*bytes.Buffer{
		makeBuffer(inpAlignA2M),
	}
	answers := [][]string{
		alignA2M,
	}
	for i := 0; i < len(tests); i++ {
		computed, err := Read(tests[i])
		if err != nil {
			t.Fatalf("%s", err)
		}
		answer := makeMSA(makeSeqs(answers[i]))
		testEqualAlign(t, computed, answer)
	}
}

func TestA3M(t *testing.T) {
	tests := []*bytes.Buffer{
		makeBuffer(inpAlignA3M),
	}
	answers := [][]string{
		alignA3M,
	}
	for i := 0; i < len(tests); i++ {
		computed, err := Read(tests[i])
		if err != nil {
			t.Fatalf("%s", err)
		}
		answer := makeMSA(makeSeqs(answers[i]))
		testEqualAlign(t, computed, answer)
	}
}

func TestReaderError(t *testing.T) {
	_, err := Read(bytes.NewBuffer(testBadAlignedInput))
	if err == nil {
		t.Fatalf("Expected an error for sequences of unequal length.")
	}
}

func BenchmarkReader(b *testing.B) {
	for i := 0; i < b.N; i++ {
		Read(bytes.NewBuffer(testAlignedInput))
	}
}

func BenchmarkReaderTrusted(b *testing.B) {
	for i := 0; i < b.N; i++ {
		ReadTrusted(bytes.NewBuffer(testAlignedInput))
	}
}

func testEqualAlign(t *testing.T, computed, answer seq.MSA) {
	if computed.Len() != answer.Len() {
		t.Fatalf("Lengths of MSAs differ: %d != %d",
			computed.Len(), answer.Len())
	}

	scomputed := makeStrings(computed.Entries)
	sanswer := makeStrings(answer.Entries)
	if len(scomputed) != len(sanswer) {
		t.Fatalf("\nLengths of entries in MSAs differ: %d != %d",
			len(scomputed), len(sanswer))
	}
	for i := 0; i < len(scomputed); i++ {
		c, a := scomputed[i], sanswer[i]
		if c != a {
			t.Fatalf("\nComputed sequence in MSA is\n\n%s\n\n"+
				"but answer is\n\n%s", c, a)
		}
	}
}

func makeBuffer(strs []string) *bytes.Buffer {
	buf := new(bytes.Buffer)
	for _, str := range strs {
		buf.WriteString(str + "\n")
	}
	return buf
}

func makeMSA(seqs []seq.Sequence) seq.MSA {
	msa := seq.NewMSA()
	msa.AddSlice(seqs)
	return msa
}

func makeSeqs(strs []string) []seq.Sequence {
	seqs := make([]seq.Sequence, len(strs))
	for i, str := range strs {
		seqs[i] = seq.Sequence{
			Name:     fmt.Sprintf("%d", i),
			Residues: []seq.Residue(str),
		}
	}
	return seqs
}

func makeStrings(seqs []seq.Sequence) []string {
	strs := make([]string, len(seqs))
	for i, s := range seqs {
		strs[i] = fmt.Sprintf("%s", s.Residues)
	}
	return strs
}
