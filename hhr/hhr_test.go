package hhr

import (
	"flag"
	"fmt"
	"log"
	"os"
)

var (
	flagReadFile = ""
)

func init() {
	flag.StringVar(&flagReadFile, "hhr", flagReadFile,
		"The HHR file to use for tests.")
	flag.Parse()

	log.SetFlags(0)
}

func ExampleRead() {
	r := getFile()

	hhr, err := Read(r)
	if err != nil {
		log.Fatalf("%s", err)
	}

	hit := hhr.Hits[9]
	fmt.Println(hit.Num)
	fmt.Println(hit.Name)
	fmt.Printf("%0.3f\n", hit.Prob)
	fmt.Println(hit.EValue)
	fmt.Println(hit.PValue)
	fmt.Println(hit.ViterbiScore)
	fmt.Println(hit.SSScore)
	fmt.Println(hit.NumAlignedCols)
	fmt.Println(hit.QueryStart)
	fmt.Println(hit.QueryEnd)
	fmt.Println(hit.TemplateStart)
	fmt.Println(hit.TemplateEnd)
	fmt.Println(hit.NumTemplateCols)
	fmt.Printf("%s\n", hit.Aligned.QSeq)
	fmt.Printf("%s\n", hit.Aligned.TSeq)
	// Output:
	// 10
	// 2fxaD
	// 0.571
	// 0.24
	// 0.00011
	// 34.3
	// 0
	// 47
	// 106
	// 154
	// 46
	// 94
	// 207
	// IGNSAFELLLEVAKSGEKGINTMDLAQVTGQDPRSVTGRIKKINH--LLTS 
	// LNINEHHILWIAY--QLNGASISEIAKFGVMHVSTAFNFSKKLEERGYLRF 
}

func getFile() *os.File {
	if len(flagReadFile) == 0 {
		log.Fatalf("Please set the '--hhr path/to/file.hhr' flag.")
	}

	r, err := os.Open(flagReadFile)
	if err != nil {
		log.Fatalf("%s", err)
	}
	return r
}
