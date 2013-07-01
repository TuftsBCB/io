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

	"github.com/TuftsBCB/seq"
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
	var last, entry seq.Sequence
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
	entries, err := NewReader(bytes.NewBuffer(testFastaInput)).ReadAll()
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

	testBytesEqual(t, testFastaInput, buf.Bytes())
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

func testLastEntry(t *testing.T, last seq.Sequence) {
	// Check the value of the last sequence.
	answer := "MQFSTVASIAAIAAVASAASNITTATVTEESTTLVTITSCEDHVCSETVSPALVSTATVTVN" +
		"DVITYTTWCPLPTTEAPKNTTSPAPTEKPTEKPTEKPTQQGSSTQTVTSYTGAAVKALPAAGALLAG" +
		"AAALLL"
	ours := fmt.Sprintf("%s", last.Residues)
	if answer != ours {
		t.Fatalf("The last sequence should be\n%s\nbut we got\n%s",
			answer, ours)
	}

	answer = "YDR134C YDR134C SGDID:S000002541"
	ours = last.Name
	if answer != ours {
		t.Fatalf("The last header should be\n%s\nbut we got\n%s",
			answer, ours)
	}
}

func BenchmarkRead(b *testing.B) {
	if len(flagFastaFile) == 0 {
		log.Fatalf("Please set the '--fasta path/to/file.fasta' flag.")
	}

	f, err := os.Open(flagFastaFile)
	if err != nil {
		log.Fatalf("%s", err)
	}

	for i := 0; i < b.N; i++ {
		_, err := f.Seek(0, os.SEEK_SET)
		if err != nil {
			log.Fatalf("%s", err)
		}

		r := NewReader(f)
		for {
			_, err := r.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				log.Fatalf("%s", err)
			}
		}
	}
}

func BenchmarkReadTrusted(b *testing.B) {
	if len(flagFastaFile) == 0 {
		log.Fatalf("Please set the '--fasta path/to/file.fasta' flag.")
	}

	f, err := os.Open(flagFastaFile)
	if err != nil {
		log.Fatalf("%s", err)
	}

	for i := 0; i < b.N; i++ {
		_, err := f.Seek(0, os.SEEK_SET)
		if err != nil {
			log.Fatalf("%s", err)
		}

		r := NewReader(f)
		r.TrustSequences = true
		for {
			_, err := r.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				log.Fatalf("%s", err)
			}
		}
	}
}
