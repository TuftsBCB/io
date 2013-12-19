package msa

import (
	"bufio"
	"bytes"
	"fmt"
	"io"

	"github.com/TuftsBCB/seq"
)

// ReadStockholm reads an MSA from a Stockholm formatted file. Note that
// features are completely ignored. This reader only checks for the Stockholm
// header (and version), and then slurps up the sequence data into an MSA.
func ReadStockholm(r io.Reader) (seq.MSA, error) {
	return readStockholm(r, false)
}

// ReadStockholmTrusted is the same as ReadStockholm, except it does not check
// if each residue is valid. This may be faster.
func ReadStockholmTrusted(r io.Reader) (seq.MSA, error) {
	return readStockholm(r, true)
}

// WriteStockholm writes the given MSA to the writer in the Stockholm format.
// This does not write any features. It only creates a minimal valid Stockholm
// file with the header (and version) along with the sequences (names and
// residues).
func WriteStockholm(w io.Writer, msa seq.MSA) error {
	var err error
	pf := func(format string, v ...interface{}) {
		if err != nil {
			return
		}
		_, err = fmt.Fprintf(w, format, v...)
	}
	pf("# STOCKHOLM 1.0\n")
	for row := 0; row < len(msa.Entries) && err == nil; row++ {
		s := msa.GetA2M(row)
		pf("%s %s\n", s.Name, s.Residues)
	}
	pf("//\n")
	return err
}

func readStockholm(r io.Reader, trusted bool) (seq.MSA, error) {
	msa := seq.NewMSA()
	ef := fmt.Errorf

	scanner := bufio.NewScanner(r)
	if scanner.Scan() {
		first := bytes.ToLower(bytes.Trim(scanner.Bytes(), " #"))
		if !bytes.Equal([]byte("stockholm 1.0"), first) {
			fmt.Printf("%s\n", first)
			return seq.MSA{}, ef("First line does not contain 'STOCKHOLM 1.0'.")
		}
	}
	for scanner.Scan() {
		line := bytes.TrimSpace(scanner.Bytes())
		if line[0] == '#' {
			continue
		}
		if line[0] == '/' && line[1] == '/' { // alignment done, says the spec
			break
		}

		pieces := bytes.Fields(line)
		residues, err := asResidues(pieces[len(pieces)-1], trusted)
		if err != nil {
			return seq.MSA{}, err
		}

		s := seq.Sequence{
			Name:     string(concat(pieces[0 : len(pieces)-1])),
			Residues: residues,
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
	if err := scanner.Err(); err != nil {
		return seq.MSA{}, err
	}
	return msa, nil
}

func asResidues(brs []byte, trusted bool) ([]seq.Residue, error) {
	rs := make([]seq.Residue, 0, len(brs))
	for _, b := range brs {
		if trusted {
			rs = append(rs, seq.Residue(b))
		} else {
			bNew, ok := translateStockholm(b)
			if !ok {
				return nil, fmt.Errorf("Invalid Stockholm residue '%c'.", b)
			} else if bNew > 0 {
				rs = append(rs, bNew)
			}
		}
	}
	return rs, nil
}

func concat(bs [][]byte) []byte {
	var ret []byte
	for _, b := range bs {
		ret = append(ret, b...)
	}
	return ret
}

func translateStockholm(b byte) (seq.Residue, bool) {
	switch {
	case b >= 'a' && b <= 'z':
		return seq.Residue(b), true
	case b >= 'A' && b <= 'Z':
		return seq.Residue(b), true
	case b == '-':
		return '-', true
	case b == '.':
		return '.', true
	}
	return 0, false
}
