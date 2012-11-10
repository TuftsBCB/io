package fasta

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"strings"
	"unicode"
)

// An Entry corresponds to an entry in a FASTA file. That is, it is a tuple
// of a single line header and a sequence over multiple lines concatenated
// into a single slice of bytes.
//
// N.B. The header part of a FASTA entry can also be in a standard format.
// As of right now, there are no facilities in this package to parse those
// formats, but they are forthcoming.
type Entry struct {
	Header   string
	Sequence []byte
}

// Output is a string in FASTA format, with the sequence wrapped at 60
// columns.
func (e Entry) String() string {
	return e.StringCols(60)
}

// StringCols returns the FASTA string corresponding to this entry with the
// sequence wrapped at the number of columns given.
//
// If cols is <= 0, then no wrapping is done.
func (e Entry) StringCols(cols int) string {
	if cols <= 0 {
		return fmt.Sprintf(">%s\n%s", e.Header, string(e.Sequence))
	}

	wrapped := make([]string, 1+((len(e.Sequence)-1)/cols))
	for i := range wrapped {
		start := cols * i
		end := start + cols
		if end > len(e.Sequence) {
			end = len(e.Sequence)
		}
		wrapped[i] = string(e.Sequence[start:end])
	}
	return fmt.Sprintf(">%s\n%s", e.Header, strings.Join(wrapped, "\n"))
}

func (e Entry) isNull() bool {
	return len(e.Header) == 0 && e.Sequence == nil
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
func (r *Reader) ReadAll() ([]Entry, error) {
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
func (r *Reader) Read() (entry Entry, err error) {
	entry, err = r.ReadEntry(TranslateNormal)
	if !entry.isNull() {
		return entry, nil
	}
	if err == io.EOF {
		return Entry{}, err
	}
	if err != nil {
		return Entry{}, fmt.Errorf("Error on line %d: %s", r.line, err)
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

// ReadEntry is exported for use in other packages that read FASTA-like files.
//
// The 'translate' function is used when sequences are checked for valid
// characters.
//
// If you're just reading FASTA files, this method SHOULD NOT be used.
func (r *Reader) ReadEntry(translate Translator) (Entry, error) {
	entry := Entry{}
	seenHeader := false

	// Before entering the main loop, we have to check to see if we've
	// already read this entry's header.
	if r.nextHeader != nil {
		entry.Header = trimHeader(r.nextHeader)
		r.nextHeader = nil
		seenHeader = true
	}
	for {
		line, err := r.buf.ReadBytes('\n')
		if err == io.EOF {
			if len(line) == 0 {
				return entry, io.EOF
			}
		} else if err != nil {
			return Entry{}, err
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
				return Entry{}, fmt.Errorf("Expected '>', got '%c'.", line[0])
			}

			// Trim the '>' and load this line into the header.
			entry.Header = trimHeader(line)
			seenHeader = true

			r.line++
			continue
		} else if line[0] == '>' {
			// This means we've begun reading the next entry.
			// So slap this line into 'nextHeader' and return the current entry.
			r.nextHeader = line

			r.line++
			return entry, nil
		}

		// Finally, time to start reading the sequence.
		// If we trust the sequences, then we can just append this line
		// willy nilly. Otherwise we've got to check each character.
		if entry.Sequence == nil {
			entry.Sequence = make([]byte, 0, 50)
		}
		if !r.TrustSequences {
			for i, b := range line {
				bNew, ok := translate(b)
				if !ok {
					return Entry{},
						fmt.Errorf("Invalid character '%c' on line %d.",
							b, r.line)
				}
				line[i] = bNew
			}
		}

		// Careful. Since we're using 'line...' notation, the bytes here
		// are *copied*. Which is what we want, lest we pin our sequences
		// so that the garbage collector can't free them.
		entry.Sequence = append(entry.Sequence, line...)

		r.line++
	}
	panic("unreachable")
}

// A Translator is a function that accepts a single character, checks whether
// it's valid, and optionally maps it to a new character.
//
// Translators are ONLY applicable to developers writing their own parsers for
// FASTA-like files. They should not be used to read regular FASTA files.
type Translator func(b byte) (byte, bool)

// TranslateNormal is the default translator for regular (and aligned) FASTA
// files.
func TranslateNormal(b byte) (byte, bool) {
	switch {
	case b >= 'a' && b <= 'z':
		b = byte(unicode.ToTitle(rune(b)))
	case b >= 'A' && b <= 'Z':
	case b == '*':
	case b == '-':
	default:
		return 0, false
	}
	return b, true
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
	buf     *bufio.Writer
}

// NewWriter createa a new FASTA writer that can write FASTA entries to
// an io.Writer.
func NewWriter(w io.Writer) *Writer {
	return &Writer{
		Columns: 60,
		buf:     bufio.NewWriter(w),
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
func (w *Writer) Write(entry Entry) error {
	s := fmt.Sprintf("%s\n", entry.StringCols(w.Columns))
	_, err := w.buf.WriteString(s)
	return err
}

// WriteAll writes a slice of FASTA entries to the underyling io.Writer, and
// calls Flush.
func (w *Writer) WriteAll(entries []Entry) error {
	for _, entry := range entries {
		if err := w.Write(entry); err != nil {
			return err
		}
	}
	return w.Flush()
}
