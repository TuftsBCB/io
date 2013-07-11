package pdb

import (
	"fmt"

	"github.com/TuftsBCB/structure"
)

// RMSD is a convenience function for computing the RMSD between two sets of
// residues, where each set is take from a chain of a PDB entry. Note that RMSD
// is only computed using carbon-alpha atoms.
//
// Each set of atoms to be used is specified by a four-tuple: a PDB entry file,
// a chain identifier, and the start and end residue numbers to use as a range.
// (Where the range is inclusive.)
//
// An error will be returned if: chainId{1,2} does not correspond to a chain
// in entry{1,2}. The ranges specified by start{1,2}-end{1,2} are not valid.
// The ranges specified by start{1,2}-end{1,2} do not correspond to precisely
// the same number of carbon-alpha atoms.
func RMSD(entry1 *Entry, chainId1 byte, start1, end1 int,
	entry2 *Entry, chainId2 byte, start2, end2 int) (float64, error) {

	chain1 := entry1.Chain(chainId1)
	if chain1 == nil {
		return 0.0, fmt.Errorf("The chain '%c' could not be found in '%s'.",
			chainId1, entry1.Path)
	}
	chain2 := entry2.Chain(chainId2)
	if chain2 == nil {
		return 0.0, fmt.Errorf("The chain '%c' could not be found in '%s'.",
			chainId2, entry2.Path)
	}
	return RMSDChains(chain1, start1, end1, chain2, start2, end2)
}

// RMSDChains is the same as RMSD, except it uses *Chain values directly.
func RMSDChains(chain1 *Chain, start1, end1 int,
	chain2 *Chain, start2, end2 int) (float64, error) {

	// In order to fetch the appropriate carbon-alpha atoms, we need to
	// traverse each chain's carbon-alpha atom slice and pick only the carbon
	// alpha atoms with residue indices in the range specified.
	struct1 := chain1.SequenceCaAtomSlice(start1-1, end1)
	struct2 := chain2.SequenceCaAtomSlice(start2-1, end2)

	// Verify that neither of the atom sets is 0.
	if struct1 == nil || len(struct1) == 0 {
		return 0.0, fmt.Errorf("The range '%d-%d' (for chain %c in %s) does "+
			"not correspond to any carbon-alpha ATOM records.",
			start1, end1, chain1.Ident, chain1.Entry.Path)
	}
	if struct2 == nil || len(struct2) == 0 {
		return 0.0, fmt.Errorf("The range '%d-%d' (for chain %c in %s) does "+
			"not correspond to any carbon-alpha ATOM records.",
			start2, end2, chain2.Ident, chain2.Entry.Path)
	}

	// If we don't have the same number of atoms from each chain, we can't
	// compute RMSD.
	if len(struct1) != len(struct2) {
		return 0.0, fmt.Errorf("The range '%d-%d' (%d ATOM records for chain "+
			"%c in %s) does not correspond to the same number of carbon-alpha "+
			"atoms as the range '%d-%d' (%d ATOM records for chain %c in %s). "+
			"It is possible that the PDB file does not contain a carbon-alpha "+
			"atom for every residue index in the ranges.",
			start1, end1, len(struct1), chain1.Ident, chain1.Entry.Path,
			start2, end2, len(struct2), chain2.Ident, chain2.Entry.Path)
	}

	// We're good to go...
	return structure.RMSD(struct1, struct2), nil
}
