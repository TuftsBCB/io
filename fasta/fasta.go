package fasta

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"strings"
	"unicode"

	"github.com/BurntSushi/bcbgo/seq"
)

// SequenceFasta returns the FASTA string corresponding to a sequence with the
// sequence wrapped at the number of columns given.
//
// If cols is <= 0, then no wrapping is done.
func SequenceFasta(s seq.Sequence, cols int) string {
	if cols <= 0 {
		return fmt.Sprintf(">%s\n%s", s.Name, s.Residues)
	}

	wrapped := SequenceString(s, cols)
	return fmt.Sprintf(">%s\n%s", s.Name, strings.Join(wrapped, "\n"))
}

// SequenceStrings chops up one long sequence into multiple strings
// based on the number of columns provided.
//
// cols must be greater than 0.
func SequenceString(s seq.Sequence, cols int) []string {
	if cols <= 0 {
		panic("cols must be greater than 0")
	}

	wrapped := make([]string, 1+((len(s.Residues)-1)/cols))
	for i := range wrapped {
		start := cols * i
		end := start + cols
		if end > len(s.Residues) {
			end = len(s.Residues)
		}
		wrapped[i] = fmt.Sprintf("%s", s.Residues[start:end])
	}
	return wrapped
}

// A Reader reads entries from FASTA encoded input.
//
// If TrustSequences is true, then sequence data will not be checked to make
// sure that it conforms to the NCBI spec. (See the Read method for details.)
// By default, TrustSequences is false.
type Reader struct {
	// When set to true, the sequences will not be checked for errors.
	// If you trust the data, this may improve performance.
	// This may be set at any time.
	TrustSequences bool
	buf            *bufio.Reader
	line           int
	nextHeader     []byte
}

// NewReader creates a new Reader that is ready to read sequences from some
// io.Reader.
func NewReader(r io.Reader) *Reader {
	return &Reader{
		TrustSequences: false,
		buf:            bufio.NewReader(r),
		line:           1,
		nextHeader:     nil,
	}
}

// ReadAll will read all entries in the FASTA input and return them as a slice.
// If an error is encountered, processing is stopped, and the error is
// returned.
func (r *Reader) ReadAll() ([]seq.Sequence, error) {
	seqs := make([]seq.Sequence, 0, 100)
	for {
		s, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}
		seqs = append(seqs, s)
	}
	return seqs, nil
}

// Read will read the next entry in the FASTA input.
// The format roughly corresponds to that described by NCBI:
// http://blast.ncbi.nlm.nih.gov/blastcgihelp.shtml
//
// In particular, the only characters allowed in the sequence section
// are a-z, A-Z, * and -. Any other character will result in an error.
//
// All lower case letters in the sequence section are translated to upper case.
//
// Blank lines, leading and trailing whitespace are always ignored (regardless
// of where they are).
//
// No distinction is made between DNA/RNA or amino acid sequences. (Currently.)
//
// It is NOT safe to call this function from multiple goroutines.
//
// If the underlying reader is seekable, it is OK to use its seek operation
// provided that you call (*Reader).SeekerReset before the next time Read is
// called. If you don't, the behavior is undefined. Moreover, seeking will
// result in erroneous line numbers in error messages. Finally, you MUST seek
// to a location that corresponds precisely to an entry boundary. i.e., the
// file pointer should be at a '>' character.
func (r *Reader) Read() (s seq.Sequence, err error) {
	s, err = r.ReadSequence(TranslateNormal)
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

// SeekerReset will reset the internal state of Reader to allow Read to be
// called at arbitrary entry boundaries in the input.
//
// See the comments for Read for more details.
func (r *Reader) SeekerReset() {
	r.nextHeader = nil
}

// ReadSequence is exported for use in other packages that read FASTA-like
// files.
//
// The 'translate' function is used when sequences are checked for valid
// characters.
//
// If you're just reading FASTA files, this method SHOULD NOT be used.
func (r *Reader) ReadSequence(translate Translator) (seq.Sequence, error) {
	s := seq.Sequence{}
	seenHeader := false

	// Before entering the main loop, we have to check to see if we've
	// already read this entry's header.
	if r.nextHeader != nil {
		s.Name = trimHeader(r.nextHeader)
		r.nextHeader = nil
		seenHeader = true
	}
	for {
		line, err := r.buf.ReadBytes('\n')
		if err == io.EOF {
			if len(line) == 0 {
				return s, io.EOF
			}
		} else if err != nil {
			return seq.Sequence{}, fmt.Errorf("Error on line %d: %s",
				r.line, err)
		}
		line = bytes.TrimSpace(line)

		// If it's empty, increment the counter and skip ahead.
		if len(line) == 0 {
			r.line++
			continue
		}

		// If we haven't seen the header yet, this better be it.
		if !seenHeader {
			if line[0] != '>' {
				return seq.Sequence{},
					fmt.Errorf("Expected '>', got '%c' on line %d.",
						line[0], r.line)
			}

			// Trim the '>' and load this line into the header.
			s.Name = trimHeader(line)
			seenHeader = true

			r.line++
			continue
		} else if line[0] == '>' {
			// This means we've begun reading the next entry.
			// So slap this line into 'nextHeader' and return the current entry.
			r.nextHeader = line

			r.line++
			return s, nil
		}

		// Finally, time to start reading the sequence.
		// If we trust the sequences, then we can just append this line
		// willy nilly. Otherwise we've got to check each character.
		if s.Residues == nil {
			s.Residues = make([]seq.Residue, 0, 50)
		}
		if r.TrustSequences {
			for _, b := range line {
				s.Residues = append(s.Residues, seq.Residue(b))
			}
		} else {
			for _, b := range line {
				bNew, ok := translate(b)
				if !ok {
					return seq.Sequence{},
						fmt.Errorf("Invalid character '%c' on line %d.",
							b, r.line)
				}

				// If the zero byte is returned from translate, then we
				// don't keep this residue around.
				if bNew > 0 {
					s.Residues = append(s.Residues, bNew)
				}
			}
		}
		r.line++
	}
	panic("unreachable")
}

// A Translator is a function that accepts a single character, checks whether
// it's valid, and optionally maps it to a new character. Additionally, if the
// zero byte is returned, then the character should not be included in the
// final sequence.
//
// Translators are ONLY applicable to developers writing their own parsers for
// FASTA-like files. They should not be used to read regular FASTA files.
type Translator func(b byte) (seq.Residue, bool)

// TranslateNormal is the default translator for regular (NOT aligned) FASTA
// files.
func TranslateNormal(b byte) (seq.Residue, bool) {
	switch {
	case b >= 'a' && b <= 'z':
		return seq.Residue(unicode.ToTitle(rune(b))), true
	case b >= 'A' && b <= 'Z':
		return seq.Residue(b), true
	case b == '*':
		return 0, true
	case b == '-':
		return '-', true
	}
	return 0, false
}

func trimHeader(line []byte) string {
	return string(bytes.TrimSpace(bytes.TrimLeft(line, ">")))
}

// A Writer writes entries to a FASTA encoded file.
//
// The 'Columns' corresponds to the number of columns at which a sequence is
// wrapped. If it's <= 0, then no wrapping will be used.
//
// The header text is never wrapped.
type Writer struct {
	// The number of columns to wrap a sequence at. By default, this
	// is set to 60. A value <= 0 will result in no wrapping.
	Columns int

	// Whether to a '*' at the end of each sequence.
	// By default, this is false.
	Asterisk bool

	buf *bufio.Writer
}

// NewWriter createa a new FASTA writer that can write FASTA entries to
// an io.Writer.
func NewWriter(w io.Writer) *Writer {
	return &Writer{
		Columns:  60,
		Asterisk: false,
		buf:      bufio.NewWriter(w),
	}
}

// Flush writes any buffered data to the underlying io.Writer.
func (w *Writer) Flush() error {
	return w.buf.Flush()
}

// Write writes a single FASTA entry to the underlying io.Writer.
//
// You may need to call Flush in order for the changes to be written.
//
// XXX: Currently, the sequence is not checked. Should it be?
func (w *Writer) Write(s seq.Sequence) error {
	var out string
	if w.Asterisk {
		out = fmt.Sprintf("%s*\n", SequenceFasta(s, w.Columns))
	} else {
		out = fmt.Sprintf("%s\n", SequenceFasta(s, w.Columns))
	}
	_, err := w.buf.WriteString(out)
	return err
}

// WriteAll writes a slice of FASTA entries to the underyling io.Writer, and
// calls Flush.
func (w *Writer) WriteAll(seqs []seq.Sequence) error {
	for _, s := range seqs {
		if err := w.Write(s); err != nil {
			return err
		}
	}
	return w.Flush()
}
