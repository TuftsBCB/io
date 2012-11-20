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

// Slice will create a new collection of secondary structure sequences where
// each sequence has been sliced with the given start/end parameters.
func (ss HHMSecondary) Slice(start, end int) HHMSecondary {
	var ssdssp, sadssp, sspred, ssconf, ssconsensus *seq.Sequence

	if ss.SSdssp != nil {
		t := ss.SSdssp.Slice(start, end)
		ssdssp = &t
	}
	if ss.SAdssp != nil {
		t := ss.SAdssp.Slice(start, end)
		sadssp = &t
	}
	if ss.SSpred != nil {
		t := ss.SSpred.Slice(start, end)
		sspred = &t
	}
	if ss.SSconf != nil {
		t := ss.SSconf.Slice(start, end)
		ssconf = &t
	}
	if ss.Consensus != nil {
		t := ss.Consensus.Slice(start, end)
		ssconsensus = &t
	}
	return HHMSecondary{
		SSdssp:    ssdssp,
		SAdssp:    sadssp,
		SSpred:    sspred,
		SSconf:    ssconf,
		Consensus: ssconsensus,
	}
}

// Slice dices up an entire HHM file. The slice indices should be in terms of
// the number of match/delete states in the underlying HMM.
// All secondary structure annotations are also sliced.
// The multiple sequence alignment is also sliced.
// The NEFF for the HHM is also re-computed as the average of all NeffM scores
// in each HMM column.
func (hhm *HHM) Slice(start, end int) *HHM {
	hmm := hhm.HMM.Slice(start, end)

	meta := hhm.Meta
	meta.Neff = 0
	for _, node := range hmm.Nodes {
		meta.Neff += node.NeffM
	}
	meta.Neff /= seq.Prob(len(hmm.Nodes))

	return &HHM{
		Meta:      meta,
		Secondary: hhm.Secondary.Slice(start, end),
		MSA:       hhm.MSA.Slice(start, end),
		HMM:       hmm,
	}
}
