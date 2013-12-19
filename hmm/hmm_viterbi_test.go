package hmm

import (
	"fmt"
	"log"
	"os"

	"github.com/TuftsBCB/seq"
)

func ExampleViterbiScore() {
	query := "IVEGQDAEVGLSPWQVMLFRKSPQELLCGASLISDRWVLTAAHCLLYPPWDKNFTVDDLLVR" +
		"IGKHSRTRYERKVEKISMLDKIYIHPRYNWKENLDRDIALLKLKRPIELSDYIHPVCLPDKQTAAKL" +
		"LHAGFKGRVTGWGNRRETWTTSVAEVQPSVLQVVNLPLVERPVCKASTRIRITDNMFCAGYKPGEGK" +
		"RGDACEGDSGGPFVMKSPYNNRWYQMGIVSWGEGCDRDGKYGFYTHVFRLKKWIQKVIDRLGS"
	squery := seq.NewSequenceString("query", query)

	hmmf, err := os.Open("sermam.hmm")
	if err != nil {
		log.Fatal(err)
	}
	defer hmmf.Close()

	profile, err := ReadHMM(hmmf)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Println(profile.HMM.ViterbiScore(squery))
	// Output:
	// NotWorking
}
