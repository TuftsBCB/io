// Example parallel shows how to read and output information about PDB files
// in parallel. It works by spawning several workers (by default, the number
// of workers equals the number of CPUs available) that each read and parse
// PDB files.
package main

import (
	"flag"
	"fmt"
	"os"
	"path"
	"runtime"

	"github.com/BurntSushi/bcbgo/pdb"
)

var (
	// flagWorkers controls how many worker goroutines are spawned to process
	// PDB files. By default, it is set to the number of CPUs.
	flagWorkers int
)

func init() {
	// For this example's purpose, force parallelism.
	runtime.GOMAXPROCS(runtime.NumCPU())

	flag.IntVar(&flagWorkers, "workers", runtime.NumCPU(),
		"The number of workers to use to process PDB files. This is "+
			"limited by the maximum allowable open file descriptors.")
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

// entryOrError values are sent on channels. They correspond to the results
// returned by 'pdb.New'.
type entryOrErr struct {
	file  string
	entry *pdb.Entry
	err   error
}

// pdbWorker is meant to be executed as a goroutine. It continually picks off
// PDB files from the pdbfiles channel, reads them and sends the results
// (including an error if one occurred) via the results channel.
func pdbWorker(pdbfiles chan string, results chan entryOrErr) {
	for pdbfile := range pdbfiles {
		// pdb.New will return an error if 'pdbfile' could not be read.
		// It will automatically decompress gzipped files that end with
		// a '.gz' extension.
		entry, err := pdb.New(pdbfile)
		results <- entryOrErr{
			file:  pdbfile,
			entry: entry,
			err:   err,
		}
	}
}

// collect is meant to be run as a single goroutine that reads values sent on
// the 'results' channel. Errors are printed to stderr and valid entry
// information is echoed to stdout.
//
// collect will stop on its own after it has processed N results, where N is
// the number of arguments passed to this command.
func collect(results chan entryOrErr, done chan struct{}) {
	i := 1
	total := flag.NArg()
	for pdb := range results {
		if pdb.err != nil {
			fmt.Fprintf(os.Stderr, "%s error: %s\n", pdb.file, pdb.err)
		} else {
			fmt.Printf("%s\n", pdb.file)
			fmt.Printf("%s\n--------------------------\n", pdb.entry)
		}

		i += 1
		if i > total {
			break
		}
	}
	fmt.Println("")
	done <- struct{}{}
}

func main() {
	if flag.NArg() < 1 {
		usage()
	}

	// Setup buffered channels through which to pass PDB file locations and
	// PDB entry (or error) values.
	pdbfiles := make(chan string, 100)
	results := make(chan entryOrErr, 100)

	// done is a ping channel that allows 'main' to wait until 'collect' has
	// finished executing.
	done := make(chan struct{}, 0)

	// collect picks off all the results and prints them.
	go collect(results, done)

	// Start the workers.
	for i := 0; i < flagWorkers; i++ {
		go pdbWorker(pdbfiles, results)
	}

	// Finally read through each argument passed to the command, and send the
	// file location over the 'pdbfiles' channel. One of the workers that we
	// just started will receive it, process it and send it to the 'results'
	// channel. (Which is processed by the 'collect' goroutine.)
	for _, pdbfile := range flag.Args() {
		pdbfiles <- pdbfile
	}

	// Wait for 'collect' to finish.
	<-done
}
