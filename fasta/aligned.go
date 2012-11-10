package fasta

import (
	"fmt"
	"io"
)

type AlignedReader struct {
	// See the exported fields of Reader for options.
	*Reader
	seqLen int // set after the first read
}

func NewAlignedReader(r io.Reader) *AlignedReader {
	return &AlignedReader{
		Reader: NewReader(r),
		seqLen: -1,
	}
}

// ReadAll will read all entries in the aligned FASTA input and return them as
// a slice.
// If an error is encountered, processing is stopped, and the error is
// returned.
// All Entry's have the same sequence length, otherwise an error occurs.
func (r *AlignedReader) ReadAll() ([]Entry, error) {
	entries := make([]Entry, 0, 100)
	for {
		entry, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}
		entries = append(entries, entry)
	}
	return entries, nil
}

// Read will read the next entry in the aligned FASTA input.
// The format roughly corresponds to that described by NCBI:
// http://blast.ncbi.nlm.nih.gov/blastcgihelp.shtml
//
// The aligned format follows immediate from the format of a regular FASTA
// file: all sequences must be the same length, '-' indicate gaps, and the
// n'th letter of any sequence is the n'th column in the alignment.
//
// See (*Reader).Read for more details.
func (r *AlignedReader) Read() (entry Entry, err error) {
	entry, err = r.ReadEntry(TranslateNormal)
	if !entry.isNull() {
		if r.seqLen == -1 {
			r.seqLen = len(entry.Sequence)
			return entry, nil
		} else {
			if r.seqLen != len(entry.Sequence) {
				return Entry{},
					fmt.Errorf("Sequence '%s' has length %d, but other "+
						"sequences have length %d.",
						entry.Header, len(entry.Sequence), r.seqLen)
			}
			return entry, nil
		}
	}
	if err == io.EOF {
		return Entry{}, err
	}
	if err != nil {
		return Entry{}, fmt.Errorf("Error on line %d: %s", r.line, err)
	}
	panic("unreachable")
}

// An AlignedWriter writes entries to an aligned FASTA encoded file.
//
// See the exported fields of Writer for options that can be set.
type AlignedWriter struct {
	*Writer
	seqLen int
}

// NewAlignedWriter createa a new aligned FASTA writer that can write FASTA
// entries to an io.Writer.
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

// Write writes a single alignment entry to the underlying io.Writer.
//
// An error is returned if the length of the sequence is not the same length
// as other sequences that have already been written.
//
// You may need to call Flush in order for the changes to be written.
//
// XXX: Currently, the sequence is not checked. Should it be?
func (w *AlignedWriter) Write(entry Entry) error {
	if w.seqLen == -1 {
		w.seqLen = len(entry.Sequence)
	} else if w.seqLen != len(entry.Sequence) {
		return fmt.Errorf("Sequence '%s' has length %d, but other sequences "+
			"have length %d.", entry.Header, len(entry.Sequence), w.seqLen)
	}
	s := fmt.Sprintf("%s\n", entry.StringCols(w.Columns))
	_, err := w.buf.WriteString(s)
	return err
}

// WriteAll writes a slice of aligned FASTA entries to the underyling
// io.Writer, and calls Flush.
func (w *AlignedWriter) WriteAll(entries []Entry) error {
	for _, entry := range entries {
		if err := w.Write(entry); err != nil {
			return err
		}
	}
	return w.Flush()
}
