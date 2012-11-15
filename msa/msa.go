package msa

import (
	"fmt"
	"io"

	"github.com/BurntSushi/bcbgo/io/fasta"
	"github.com/BurntSushi/bcbgo/seq"
)

func translateA2M(b byte) (seq.Residue, bool) {
	switch {
	case b >= 'a' && b <= 'z':
		return seq.Residue(b), true
	case b >= 'A' && b <= 'Z':
		return seq.Residue(b), true
	case b == '*':
		return 0, true
	case b == '-':
		return '-', true
	case b == '.':
		return '.', true
	}
	return 0, false
}

// Read will read a single MSA from the input, where the input can be formatted
// in FASTA, A2M or A3M formats. Sequences are read until io.EOF.
func Read(reader io.Reader) (seq.MSA, error) {
	r := fasta.NewReader(reader)
	r.TrustSequences = false
	return read(r)
}

// ReadTrusted will read a single MSA from trusted input, where the input can
// be formatted/ in FASTA, A2M or A3M formats. Sequences are read until io.EOF.
//
// "Trust" in this context means that the input doesn't contain any illegal
// characters in the sequence. Trusting the input should be faster.
func ReadTrusted(reader io.Reader) (seq.MSA, error) {
	r := fasta.NewReader(reader)
	r.TrustSequences = true
	return read(r)
}

func read(r *fasta.Reader) (seq.MSA, error) {
	msa := seq.NewMSA()
	for {
		s, err := readSequence(r)
		if err == io.EOF {
			break
		}
		if err != nil {
			return seq.MSA{}, err
		}

		msa.Add(s)
		if len(msa.Entries) > 1 {
			// We can't use 's' directly, because a sequence added to an MSA
			// may be modified if it isn't already in A2M format.
			lastEntry := msa.Entries[len(msa.Entries)-1]
			if lastEntry.Len() != msa.Entries[0].Len() {
				return seq.MSA{},
					fmt.Errorf("Sequence '%s' has length %d, but other "+
						"sequences have length %d.",
						s.Name, lastEntry, msa.Entries[0].Len())
			}
		}
	}
	return msa, nil
}

func readSequence(r *fasta.Reader) (s seq.Sequence, err error) {
	s, err = r.ReadSequence(translateA2M) // A2M encompasses FASTA/A3M
	if !s.IsNull() {
		return s, nil
	}
	if err == io.EOF {
		return seq.Sequence{}, err
	}
	if err != nil {
		return seq.Sequence{}, err
	}
	panic("unreachable")
}

type formatSeq func(row int) seq.Sequence

// WriteFasta writes a multiple sequence alignment to the output in aligned
// FASTA format. Aligned FASTA format uses upper case characters to indicate
// matches, lower case characters to indicate insertions, and '-' characters to
// indicate deletions/insertions.
func WriteFasta(w io.Writer, msa seq.MSA) error {
	formatter := func(row int) seq.Sequence {
		return msa.GetFasta(row)
	}
	return write(w, msa, formatter)
}

// WriteA2M writes a multiple sequence alignment to the output in
// A2M format. A2M format uses upper case characters to indicate
// matches, lower case and '.' characters to indicate insertions, and '-'
// characters to indicate deletions.
func WriteA2M(w io.Writer, msa seq.MSA) error {
	formatter := func(row int) seq.Sequence {
		return msa.GetA2M(row)
	}
	return write(w, msa, formatter)
}

// WriteA3M writes a multiple sequence alignment to the output in
// A3M format. A3M format uses upper case characters to indicate
// matches, lower case characters to indicate insertions, and '-'
// characters to indicate deletions.
//
// A3M format is a more compact way to write an MSA than FASTA or A2M.
func WriteA3M(w io.Writer, msa seq.MSA) error {
	formatter := func(row int) seq.Sequence {
		return msa.GetA3M(row)
	}
	return write(w, msa, formatter)
}

func write(writer io.Writer, msa seq.MSA, formatter formatSeq) error {
	w := fasta.NewWriter(writer)
	w.Asterisk = false
	for row := range msa.Entries {
		if err := w.Write(formatter(row)); err != nil {
			return err
		}
	}
	return w.Flush()
}
