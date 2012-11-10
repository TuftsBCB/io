// Example sequential shows how to read and output information about PDB files
// in sequence.
package main

import (
	"flag"
	"fmt"
	"os"
	"path"

	"github.com/BurntSushi/bcbgo/pdb"
)

func main() {
	if flag.NArg() < 1 {
		usage()
	}
	for _, pdbfile := range flag.Args() {
		// pdb.New will return an error if 'pdbfile' could not be read.
		// It will automatically decompress gzipped files that end with
		// a '.gz' extension.
		entry, err := pdb.New(pdbfile)
		if err != nil {
			fmt.Println(pdbfile)
			fmt.Println(err)
			return
		}
		fmt.Printf("%s\n", pdbfile)
		fmt.Printf("%s\n--------------------------\n", entry)
	}
	fmt.Println("")
}

func init() {
	flag.Usage = usage
	flag.Parse()
}

func usage() {
	fmt.Fprintf(os.Stderr, "Usage: %s pdb-file [ pdb-file ... ]\n",
		path.Base(os.Args[0]))
	flag.PrintDefaults()
	fmt.Fprintf(os.Stderr, "\nex. './%s ../../../data/samples/*.pdb'\n",
		path.Base(os.Args[0]))
	os.Exit(1)
}
