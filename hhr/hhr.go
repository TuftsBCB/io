package hhr

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"strconv"
	"unicode"

	"github.com/BurntSushi/bcbgo/seq"
)

type HHR struct {
	Query        string
	MatchColumns int
	NumSeqs      string
	Neff         seq.Prob
	SearchedHMMs int
	Date         string
	Command      string
	Hits         []Hit
}

type Hit struct {
	Num             int
	Name            string
	Prob            float64
	EValue          float64
	PValue          float64
	ViterbiScore    float64
	SSScore         float64
	NumAlignedCols  int
	QueryStart      int
	QueryEnd        int
	TemplateStart   int
	TemplateEnd     int
	NumTemplateCols int
}

func Read(r io.Reader) (*HHR, error) {
	bmeta, bhits := new(bytes.Buffer), new(bytes.Buffer)
	mode := 1 // 1 for meta, 2 for hits

	buf := bufio.NewReader(r)
MAIN:
	for {
		line, err := buf.ReadBytes('\n')
		if err == io.EOF && len(line) == 0 {
			break
		}
		if err != nil && err != io.EOF {
			return nil, fmt.Errorf("Error reading hhm: %s", err)
		}
		line = trim(line)

		// Unfortunately, HHR files really aren't meant to be parsed. At all.
		// They are totally human readable, so I expect this parser is pretty
		// brittle.
		// There are two portions that we care about in this file: the meta
		// information and the hit list. (The rest we don't parse, which are
		// alignments.) The meta information is always first. We don't
		// transition into the second part (the hit list) until we see
		// the header line starting with " No Hit". After that, a hit is on
		// each subsequent line. Once an empty line is reached, the hit list
		// is complete and we quit.
		switch {
		case mode == 1 && len(line) == 0: // skip empty lines in meta data
			continue MAIN
		case mode == 1 && hasPrefix(line, "No Hit"): // hit list begins
			mode = 2
			continue MAIN
		case mode == 2 && len(line) == 0: // done with hit list
			break MAIN
		}
		switch mode {
		case 1: // meta data
			_, err = bmeta.Write(append(line, '\n'))
		case 2: // hit list
			_, err = bhits.Write(append(line, '\n'))
		default:
			panic(fmt.Sprintf("BUG: Unknown mode: %d", mode))
		}
		if err != nil {
			return nil, fmt.Errorf("Error updating buffer for hhr: %s", err)
		}
	}

	hhr, err := readMeta(bmeta)
	if err != nil {
		return nil, fmt.Errorf("Error reading meta data from hhr: %s", err)
	}

	hhr.Hits, err = readHits(bhits)
	if err != nil {
		return nil, fmt.Errorf("Error reading hit list from hhr: %s", err)
	}

	return hhr, nil
}

func readMeta(buf *bytes.Buffer) (*HHR, error) {
	hhr := &HHR{}
	for {
		line, err := buf.ReadBytes('\n')
		if err == io.EOF && len(line) == 0 {
			break
		}
		if err != nil && err != io.EOF {
			return nil, err
		}
		line = trim(line)

		switch {
		case hasPrefix(line, "Query"):
			hhr.Query = str(line[5:])
		case hasPrefix(line, "Match_columns"):
			hhr.MatchColumns, err = strconv.Atoi(str(line[13:]))
			if err != nil {
				return nil, err
			}
		case hasPrefix(line, "No_of_seqs"):
			hhr.NumSeqs = str(line[10:])
		case hasPrefix(line, "Neff"):
			f, err := strconv.ParseFloat(str(line[4:]), 64)
			if err != nil {
				return nil, err
			}
			hhr.Neff = seq.Prob(f)
		case hasPrefix(line, "Searched_HMMs"):
			hhr.SearchedHMMs, err = strconv.Atoi(str(line[13:]))
			if err != nil {
				return nil, err
			}
		case hasPrefix(line, "Date"):
			hhr.Date = str(line[4:])
		case hasPrefix(line, "Command"):
			hhr.Command = str(line[7:])
		}
	}
	return hhr, nil
}

func readHits(buf *bytes.Buffer) (hits []Hit, err error) {
	hits = make([]Hit, 0, 10)
	for {
		line, err := buf.ReadBytes('\n')
		if err == io.EOF && len(line) == 0 {
			break
		}
		if err != nil && err != io.EOF {
			return nil, err
		}
		line = trim(line)

		if len(line) == 0 {
			panic(fmt.Sprintf("BUG: Empty line in hit list."))
		}
		hit := Hit{}

		// Grrrr. Variable length fields with fixed width fields without
		// a delimiter.
		numRest := bytes.SplitN(line, []byte{' '}, 2)

		hit.Num, err = strconv.Atoi(str(numRest[0]))
		if err != nil {
			return nil, err
		}

		// Thankfully, the hit name is the only variable length field.
		hit.Name = str(numRest[1][0:32])

		// So now we can split the rest by whitespace.
		// Except for when there isn't any whitespace to delimit columns!
		// Oh my. *facepalm*
		delim := func(r rune) bool {
			return unicode.IsSpace(r) || r == '('
		}
		rest := bytes.FieldsFunc(numRest[1][32:], delim)

		hit.Prob, err = readFloat(rest[0])
		if err != nil {
			return nil, err
		}
		hit.Prob /= 100.0

		hit.EValue, err = readFloat(rest[1])
		if err != nil {
			return nil, err
		}

		hit.PValue, err = readFloat(rest[2])
		if err != nil {
			return nil, err
		}

		hit.ViterbiScore, err = readFloat(rest[3])
		if err != nil {
			return nil, err
		}

		hit.SSScore, err = readFloat(rest[4])
		if err != nil {
			return nil, err
		}

		hit.NumAlignedCols, err = strconv.Atoi(str(rest[5]))
		if err != nil {
			return nil, err
		}

		// query/template range look like '{start}-{end}' where '{...}' is
		// an integer.
		queryRange := bytes.Split(rest[6], []byte{'-'})
		templateRange := bytes.Split(rest[7], []byte{'-'})

		hit.QueryStart, err = strconv.Atoi(str(queryRange[0]))
		if err != nil {
			return nil, err
		}

		hit.QueryEnd, err = strconv.Atoi(str(queryRange[1]))
		if err != nil {
			return nil, err
		}

		hit.TemplateStart, err = strconv.Atoi(str(templateRange[0]))
		if err != nil {
			return nil, err
		}

		hit.TemplateEnd, err = strconv.Atoi(str(templateRange[1]))
		if err != nil {
			return nil, err
		}

		numPart := rest[8][1 : len(rest[8])-1] // i.e., remove parens in '(52)'.
		hit.NumTemplateCols, err = strconv.Atoi(str(numPart))
		if err != nil {
			return nil, err
		}

		hits = append(hits, hit)
	}
	return hits, nil
}

func hasPrefix(bs []byte, prefix string) bool {
	return bytes.HasPrefix(bs, []byte(prefix))
}

func trim(bs []byte) []byte {
	return bytes.TrimSpace(bs)
}

func readFloat(bs []byte) (float64, error) {
	f, err := strconv.ParseFloat(str(bs), 64)
	if err != nil {
		return 0, err
	}
	return f, nil
}

func str(bs []byte) string {
	return string(bytes.TrimSpace(bs))
}
