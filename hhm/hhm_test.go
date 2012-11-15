package hhm

import (
	"flag"
	"log"
	"os"
	"testing"
)

var (
	flagReadFile  = ""
	flagWriteFile = ""
)

func init() {
	flag.StringVar(&flagReadFile, "hhm", flagReadFile,
		"The HHM file to use for tests.")
	flag.StringVar(&flagWriteFile, "out", flagWriteFile,
		"The HHM file to write to.")
	flag.Parse()

	log.SetFlags(0)
}

func TestReadWrite(t *testing.T) {
	r, w := getFiles()

	hhm, err := Read(r)
	if err != nil {
		t.Fatalf("%s", err)
	}

	if err := Write(w, hhm); err != nil {
		t.Fatalf("%s", err)
	}
}

func BenchmarkReadWrite(b *testing.B) {
	for i := 0; i < b.N; i++ {
		r, w := getFiles()

		hhm, err := Read(r)
		if err != nil {
			log.Fatalf("%s", err)
		}

		if err := Write(w, hhm); err != nil {
			log.Fatalf("%s", err)
		}
	}
}

func getFiles() (*os.File, *os.File) {
	if len(flagReadFile) == 0 {
		log.Fatalf("Please set the '--hhm path/to/file.hhm' flag.")
	}
	if len(flagWriteFile) == 0 {
		log.Fatalf("Please set the '--out path/to/file.hhm' flag.")
	}

	r, err := os.Open(flagReadFile)
	if err != nil {
		log.Fatalf("%s", err)
	}

	w, err := os.Create(flagWriteFile)
	if err != nil {
		log.Fatalf("%s", err)
	}
	return r, w
}
