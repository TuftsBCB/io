package hhm

import (
	"github.com/BurntSushi/bcbgo/seq"
)

var hmmScale = 1000.0

type HHM struct {
	Meta      HHMMeta
	Secondary HHMSecondary
	MSA       seq.MSA
	HMM       *seq.HMM
}

type HHMMeta struct {
	// Corresponds to the HHsearch version for which this format was
	// first introduced.
	FormatVersion string

	// The name of the HMM and an optional description.
	Name string

	// Family if the sequence is from SCOP or PFAM.
	Fam string

	// Command that was used to generate the file.
	Com string

	// Diversity of the alignment, calculated as the exp of the negative
	// entropy averaged over all columns of the alignment.
	Neff seq.Prob

	// ???
	File string

	// Number of match states and number of columns. Doesn't appear to be
	// in machine readable format.
	Leng string

	// Pre-filter details? Not machine readable.
	Filt string

	// EVD parameters. (Not used in HHsuite 2+, I think.)
	EvdLambda, EvdMu float64

	// Whether the HMM has pseudo count correction.
	Pct bool

	// Date file was generated.
	// e.g., 'Sat Nov 10 21:31:12 2012'
	Date string

	// ???
	Desc string
}

// An HHMSecondary represents secondary structure information that *may* be
// in an HHM file. Any or all members of an HHMSecondary may be nil.
type HHMSecondary struct {
	SSdssp    *seq.Sequence
	SAdssp    *seq.Sequence
	SSpred    *seq.Sequence
	SSconf    *seq.Sequence
	Consensus *seq.Sequence
}
