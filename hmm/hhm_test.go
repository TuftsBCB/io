package hhm

import (
	"flag"
	"fmt"
	"log"
	"os"
	"strings"
	"testing"
)

var (
	flagReadFile   = ""
	flagWriteFile  = ""
	flagSliceStart = -1
	flagSliceEnd   = -1
)

func init() {
	flag.StringVar(&flagReadFile, "hhm", flagReadFile,
		"The HHM file to use for tests.")
	flag.StringVar(&flagWriteFile, "out", flagWriteFile,
		"The HHM file to write to.")
	flag.IntVar(&flagSliceStart, "s", flagSliceStart,
		"The start of the HHM slice.")
	flag.IntVar(&flagSliceEnd, "e", flagSliceEnd,
		"The start of the HHM slice.")
	flag.Parse()

	log.SetFlags(0)
}

func TestSlice(t *testing.T) {
	r, _ := getFiles()
	s, e := getSlice()

	wSuffix := fmt.Sprintf("_%d-%d.test.hhm", s, e)
	w, err := os.Create(strings.Replace(flagReadFile, ".hhm", wSuffix, -1))
	if err != nil {
		log.Fatalf("%s", err)
	}

	hhm, err := Read(r)
	if err != nil {
		t.Fatalf("%s", err)
	}
	if err := Write(w, hhm.Slice(s, e)); err != nil {
		t.Fatalf("%s", err)
	}
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

func getSlice() (int, int) {
	if flagSliceStart == -1 {
		log.Fatalf("Please set the '-s INDEX' flag.")
	}
	if flagSliceEnd == -1 {
		log.Fatalf("Please set the '-e INDEX' flag.")
	}
	return flagSliceStart, flagSliceEnd
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
