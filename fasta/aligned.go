package fasta

import (
	"fmt"
	"io"

	"github.com/BurntSushi/bcbgo/seq"
)

// An AlignedReader reads a MSA from aligned FASTA encoded input.
//
// The TrustSequence field embedded in the Reader may be used.
type AlignedReader struct {
	// See the exported fields of Reader for options.
	*Reader
	seqLen int // set after the first read
}

// NewAlignedReader creates a new AlignedReader that is ready to read an
// MSA from some io.Reader.
func NewAlignedReader(r io.Reader) *AlignedReader {
	return &AlignedReader{
		Reader: NewReader(r),
		seqLen: -1,
	}
}

// Read will read all entries in the aligned FASTA input and return them as
// an MSA (multiple sequence alignment).
// If an error is encountered, processing is stopped, and the error is
// returned.
// All Entry's have the same sequence length, otherwise an error occurs.
func (r *AlignedReader) Read() (seq.MSA, error) {
	seqs := make([]seq.Sequence, 0, 100)
	for {
		s, err := r.read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return seq.MSA{}, err
		}
		seqs = append(seqs, s)
	}
	msa := seq.NewMSA()
	msa.AddSlice(seqs)
	return msa, nil
}

// read will read the next entry in the aligned FASTA input.
// The format roughly corresponds to that described by NCBI:
// http://blast.ncbi.nlm.nih.gov/blastcgihelp.shtml
//
// The aligned format follows immediately from the format of a regular FASTA
// file: all sequences must be the same length, '-' indicate gaps, and the
// n'th letter of any sequence is the n'th column in the alignment.
//
// See (*Reader).Read for more details.
func (r *AlignedReader) read() (s seq.Sequence, err error) {
	s, err = r.ReadSequence(TranslateAligned)
	if !seqIsNull(s) {
		if r.seqLen == -1 {
			r.seqLen = len(s.Residues)
			return s, nil
		} else {
			if r.seqLen != len(s.Residues) {
				return seq.Sequence{},
					fmt.Errorf("Sequence '%s' has length %d, but other "+
						"sequences have length %d.",
						s.Name, len(s.Residues), r.seqLen)
			}
			return s, nil
		}
	}
	if err == io.EOF {
		return seq.Sequence{}, err
	}
	if err != nil {
		return seq.Sequence{}, fmt.Errorf("Error on line %d: %s", r.line, err)
	}
	panic("unreachable")
}

// TranslateAligned tests each residue in a sequence in an alignment.
// The only characters allowed are: a-z, A-Z, ., - and *.
// If a '*' character is found, it is deleted. No other translation is
// performed.
func TranslateAligned(b byte) (seq.Residue, bool) {
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

// An AlignedWriter writes entries to an aligned FASTA encoded file.
//
// See the exported fields of Writer for options that can be set.
type AlignedWriter struct {
	*Writer
	seqLen int
}

// NewAlignedWriter createa a new aligned FASTA writer that can write an MSA
// to an io.Writer.
func NewAlignedWriter(w io.Writer) *AlignedWriter {
	return &AlignedWriter{
		Writer: NewWriter(w),
		seqLen: -1,
	}
}

// Flush writes any buffered data to the underlying io.Writer.
func (w *AlignedWriter) Flush() error {
	return w.Writer.Flush()
}

// write writes a single alignment entry to the underlying io.Writer.
//
// An error is returned if the length of the sequence is not the same length
// as other sequences that have already been written.
//
// You may need to call Flush in order for the changes to be written.
//
// XXX: Currently, the sequence is not checked. Should it be?
func (w *AlignedWriter) write(msa seq.MSA, row int) error {
	s := msa.GetFasta(row)
	if w.seqLen == -1 {
		w.seqLen = len(s.Residues)
	} else if w.seqLen != len(s.Residues) {
		return fmt.Errorf("Sequence '%s' has length %d, but other sequences "+
			"have length %d.", s.Name, len(s.Residues), w.seqLen)
	}
	out := fmt.Sprintf("%s\n", SequenceFasta(s, w.Columns))
	_, err := w.buf.WriteString(out)
	return err
}

// Write writes an MSA of aligned FASTA entries to the underyling
// io.Writer, and calls Flush.
func (w *AlignedWriter) Write(msa seq.MSA) error {
	for row := range msa.Entries {
		if err := w.write(msa, row); err != nil {
			return err
		}
	}
	return w.Flush()
}
