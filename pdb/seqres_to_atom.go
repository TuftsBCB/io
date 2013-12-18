package pdb

import (
	"unicode"

	"github.com/TuftsBCB/seq"
)

func (m Model) seqAtomsGuess() []*Residue {
	seqres := m.Chain.Sequence
	mapping := make([]*Residue, len(seqres))
	if len(seqres) != len(m.Residues) {
		// This is a last ditch effort. Use the ATOM sequence number as an
		// index into the SEQRES residues.
		for _, r := range m.Residues {
			si := r.SequenceNum - 1
			if si >= 0 && si < len(seqres) && seqres[si] == r.Name {
				mapping[si] = r
			}
		}
		return mapping
	}

	for i, r := range m.Residues {
		if i < len(seqres) && seqres[i] == r.Name {
			mapping[i] = r
		}
	}
	return mapping
}

func (m Model) seqAtomsAlign() []*Residue {
	seqres := m.Chain.Sequence
	atomResidues := make([]seq.Residue, len(m.Residues))
	for i, r := range m.Residues {
		atomResidues[i] = r.Name
	}
	aligned := seq.NeedlemanWunsch(seqres, atomResidues, seq.SubstBlosum62)

	mapped := make([]*Residue, len(seqres))
	atomi := 0
	for i, r := range aligned.B {
		if r == '-' {
			mapped[i] = nil
		} else {
			mapped[i] = m.Residues[atomi]
			atomi++
		}
	}
	return mapped
}

// Attempts to accomplish the same thing as seqAtomsWithMissing, but instead
// of mapping one residue at a time, we map *chunks* at a time. That is,
// a chunk is any group of contiguous residues.
func (m Model) seqAtomsChunksMerge() []*Residue {
	seqres := m.Chain.Sequence

	// Check to make sure that the total number of missing residues, plus the
	// total number of ATOM record residues equal the total number of residues
	// in the SEQRES records. Otherwise, the merge will fail.
	if len(seqres) != len(m.Chain.Missing)+len(m.Residues) {
		return m.seqAtomsAlign()
	}

	result := make([]*Residue, len(seqres))
	mchunks := chunk(m.Chain.Missing, m.Residues)
	fchunks := chunk(m.Residues, m.Chain.Missing)

	// If the PDB file is corrupted, a merge will fail.
	// So we fall back to alignment.
	if ok := merge(result, 0, seqres, nil, mchunks, fchunks); !ok {
		return m.seqAtomsAlign()
	}

	// X out any residues that don't have ATOM records.
	for i := range result {
		if result[i].Atoms == nil {
			result[i] = nil
		}
	}
	return result
}

func merge(result []*Residue, end int, seqres []seq.Residue,
	chunkToMerge []*Residue, mchunks, fchunks [][]*Residue) bool {

	if chunkToMerge != nil {
		i := 0
		for ; i < len(chunkToMerge); i++ {
			if chunkToMerge[i].Name != seqres[end+i] {
				return false
			}
			result[end+i] = chunkToMerge[i]
		}
		end += i
	}
	switch {
	case len(mchunks) == 0 && len(fchunks) == 0:
		return true
	case len(mchunks) == 0:
		return merge(result, end, seqres, fchunks[0], nil, fchunks[1:])
	case len(fchunks) == 0:
		return merge(result, end, seqres, mchunks[0], mchunks[1:], nil)
	}

	// This is a little weird. We really want to alternate chunks, since
	// a chunk is split presumably where there are holes.
	// But, we don't force it.
	// So we express a preference for it by trying the alternate first, and
	// fallback to similar chunk second.
	switch {
	case len(chunkToMerge) > 0 && chunkToMerge[0].Atoms == nil:
		// The current chunk is from missing, so prefer filled.
		return merge(result, end, seqres, fchunks[0], mchunks, fchunks[1:]) ||
			merge(result, end, seqres, mchunks[0], mchunks[1:], fchunks)
	case len(chunkToMerge) > 0 && chunkToMerge[0].Atoms != nil:
		// The current chunk is from filled, so prefer missing.
		return merge(result, end, seqres, mchunks[0], mchunks[1:], fchunks) ||
			merge(result, end, seqres, fchunks[0], mchunks, fchunks[1:])
	}

	// doesn't matter what we do here.
	// (this is the beginning)
	return merge(result, end, seqres, mchunks[0], mchunks[1:], fchunks) ||
		merge(result, end, seqres, fchunks[0], mchunks, fchunks[1:])
}

// chunk splits up residues in a list of contiguous segments.
func chunk(residues []*Residue, other []*Residue) [][]*Residue {
	if len(residues) == 0 {
		return nil
	}

	chunks := make([][]*Residue, 0)

	cur := make([]*Residue, 1)
	last := residues[0]
	cur[0] = last
	for i := 1; i < len(residues); i++ {
		r := residues[i]

		if last.isContiguous(r, other) {
			cur = append(cur, r)
			last = r
			continue
		}

		chunks = append(chunks, cur)
		cur = make([]*Residue, 1)
		cur[0] = r
		last = r
	}
	chunks = append(chunks, cur)

	return chunks
}

func (a Residue) less(b Residue) bool {
	return a.SequenceNum < b.SequenceNum ||
		(a.SequenceNum == b.SequenceNum && a.InsertionCode < b.InsertionCode)
}

func (a Residue) equals(b Residue) bool {
	return a.SequenceNum == b.SequenceNum && a.InsertionCode == b.InsertionCode
}

func (a *Residue) isContiguous(b *Residue, other []*Residue) bool {
	// If any residue in "other" is next to a or b, then we cannot claim
	// contiguity.
	for i := range other {
		if a.SequenceNum == other[i].SequenceNum {
			return false
		}
		if b.SequenceNum == other[i].SequenceNum {
			return false
		}
	}
	return a.isNext(b)
}

func (a *Residue) isNext(b *Residue) bool {
	asn, bsn := a.SequenceNum, b.SequenceNum

	// We make the assumption that any two residues with the same sequence
	// number are always contiguous.
	// Note that this has problems. Consider the case with two ATOM records
	// with sequence numbers 36 and 37, but a missing residue with sequence
	// number 36 and insertion code A. So it should go 36 -> 36A -> 37.
	// To remedy this, the caller must make sure that all missing (or filled)
	// residues don't have the same sequence number.
	if asn == bsn {
		return true
	}
	if asn+1 == bsn || asn-1 == bsn {
		return true
	}
	return false
}

type residues []*Residue

func (rs residues) String() string {
	bs := make([]byte, len(rs))
	for i, r := range rs {
		switch {
		case r == nil:
			bs[i] = '-'
		case r.Atoms == nil:
			bs[i] = byte(unicode.ToLower(rune(r.Name)))
		default:
			bs[i] = byte(r.Name)
		}
	}
	return string(bs)
}

type gappedResidues []*Residue

func (rs gappedResidues) String() string {
	bs := make([]byte, len(rs))
	for i, r := range rs {
		if r == nil || r.Atoms == nil {
			bs[i] = '-'
		} else {
			bs[i] = byte(r.Name)
		}
	}
	return string(bs)
}
