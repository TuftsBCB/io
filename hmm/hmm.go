package hmm

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
)

// WriteHMM writes an hmm file that can be read by HMMER.
//
// Currently, this function is broken. Do not use.
func WriteHMM(w io.Writer, hmm *HMM) error {
	buf := bufio.NewWriter(w)
	if err := writeMeta(buf, hmm.Meta, hmm.HMM); err != nil {
		return err
	}
	if err := writeHMM(buf, hmm.HMM); err != nil {
		return err
	}
	if _, err := buf.WriteString("//\n"); err != nil {
		return err
	}
	return buf.Flush()
}

// ReadHMM reads an hmm file produced by HMMER.
//
// Currently, this function is broken. Do not use.
func ReadHMM(r io.Reader) (*HMM, error) {
	// An hmm file as two logical sections: 1) Meta data and 2) the HMM.
	bmeta, bhmm := new(bytes.Buffer), new(bytes.Buffer)

	scanner := bufio.NewScanner(r)
	parsingMeta := true
MAIN:
	for scanner.Scan() {
		line := trim(scanner.Bytes())

		// First check if we should do a mode change.
		// i.e., a line starting with 'HMM' means switching from parsing meta
		// data to parsing the actual HMM.
		// Finally, a line starting with '\\' means STOP.
		switch {
		case hasPrefix(line, "HMM"):
			parsingMeta = false
		case hasPrefix(line, "//"):
			break MAIN
		}

		// Now add the bytes to the appropriate buffer based on the mode.
		var err error
		if parsingMeta {
			_, err = bmeta.Write(append(line, '\n'))
		} else {
			_, err = bhmm.Write(append(line, '\n'))
		}
		if err != nil {
			return nil, fmt.Errorf("Error updating buffer for hmm: %s", err)
		}
	}
	if err := scanner.Err(); err != nil {
		return nil, fmt.Errorf("Error reading hmm: %s", err)
	}

	meta, err := readMeta(bmeta)
	if err != nil {
		return nil, fmt.Errorf("Error reading meta data from hmm: %s", err)
	}

	hmm, err := readHMM(bhmm)
	if err != nil {
		return nil, fmt.Errorf("Error reading HMM data from hmm: %s", err)
	}

	return &HMM{
		Meta: meta,
		HMM:  hmm,
	}, nil
}
