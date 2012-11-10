package fasta

import (
	"bytes"
	"flag"
	"fmt"
	"io"
	"io/ioutil"
	"log"
	"os"
	"testing"
)

var (
	flagFastaFile = ""
	flagFastaOut  = ""
)

func init() {
	flag.StringVar(&flagFastaFile, "fasta", flagFastaFile,
		"The fasta file to use for benchmarks.")
	flag.StringVar(&flagFastaOut, "fasta-out", flagFastaOut,
		"The fasta file write to (for debugging).")
	flag.Parse()

	log.SetFlags(0)
}

func TestReadAll(t *testing.T) {
	r := NewReader(bytes.NewBuffer(testFastaInput))
	all, err := r.ReadAll()
	if err != nil {
		t.Fatalf("%s", err)
	}
	testLastEntry(t, all[len(all)-1])
}

func TestRead(t *testing.T) {
	var last, entry Entry
	var err error

	r := NewReader(bytes.NewBuffer(testFastaInput))
	for {
		entry, err = r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("%s", err)
		}
		last = entry
	}

	testLastEntry(t, last)
}

func TestReadWrite(t *testing.T) {
	if len(flagFastaFile) == 0 {
		t.Fatalf("Please set the '--fasta path/to/file.fasta' flag.")
	}

	f, err := os.Open(flagFastaFile)
	if err != nil {
		t.Fatalf("%s", err)
	}

	contents, err := ioutil.ReadAll(f)
	if err != nil {
		t.Fatalf("%s", err)
	}

	entries, err := NewReader(bytes.NewBuffer(contents)).ReadAll()
	if err != nil {
		t.Fatalf("%s", err)
	}

	buf := new(bytes.Buffer)
	if err := NewWriter(buf).WriteAll(entries); err != nil {
		t.Fatalf("%s", err)
	}

	if len(flagFastaOut) != 0 {
		bs := buf.Bytes()
		if err := ioutil.WriteFile(flagFastaOut, bs, 0666); err != nil {
			t.Fatalf("%s", err)
		}
	}

	testBytesEqual(t, contents, buf.Bytes())
}

func testBytesEqual(t *testing.T, bs1, bs2 []byte) {
	if len(bs1) != len(bs2) {
		t.Fatalf("Lengths not equal: %d != %d", len(bs1), len(bs2))
	}
	for i := 0; i < len(bs1); i++ {
		if bs1[i] != bs2[i] {
			t.Fatalf("Byte %d not equal: %c != %c", i, bs1[i], bs2[i])
		}
	}
}

func testLastEntry(t *testing.T, last Entry) {
	// Check the value of the last sequence.
	answer := "MQFSTVASIAAIAAVASAASNITTATVTEESTTLVTITSCEDHVCSETVSPALVSTATVTVN" +
		"DVIT*YTTWCPLPTTEAPKNTTSPAPTEKPTEKPTEKPTQQGSSTQTVTSYTGAAVKALPAAGALLAG" +
		"AAALLL*"
	ours := string(last.Sequence)
	if answer != ours {
		t.Fatalf("The last sequence should be\n%s\nbut we got\n%s",
			answer, ours)
	}

	answer = "YDR134C YDR134C SGDID:S000002541, Chr IV from 721481-721071, " +
		"Genome Release 64-1-1, reverse complement, " +
		"pseudogene, \"Hypothetical protein\""
	ours = last.Header
	if answer != ours {
		t.Fatalf("The last header should be\n%s\nbut we got\n%s",
			answer, ours)
	}
}

func ExampleRead() {
	var last Entry

	if len(flagFastaFile) == 0 {
		log.Fatalf("Please set the '--fasta path/to/file.fasta' flag.")
	}

	f, err := os.Open(flagFastaFile)
	if err != nil {
		log.Fatalf("%s", err)
	}

	r := NewReader(f)
	for {
		entry, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatalf("%s", err)
		}
		last = entry
	}
	fmt.Printf("%d\n", 0*len(last.Header))
	// Output: 0
}
