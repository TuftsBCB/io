package hhr

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"log"
	"strconv"
	"strings"
	"unicode"

	"github.com/BurntSushi/bcbgo/seq"
)

var _ = log.Println

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
	Aligned         Alignment
}

type Alignment struct {
	QSeq, QConsensus, QDssp, QPred, QConf []seq.Residue
	TSeq, TConsensus, TDssp, TPred, TConf []seq.Residue
}

func Read(r io.Reader) (*HHR, error) {
	mode := 1 // 1 for meta, 2 for hits, 3 for alignments
	bmeta := new(bytes.Buffer)
	bhits := new(bytes.Buffer)
	balign := new(bytes.Buffer)

	buf := bufio.NewReader(r)
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
			continue
		case mode == 1 && hasPrefix(line, "No Hit"): // hit list begins
			mode = 2
			continue
		case mode == 2 && len(line) == 0: // done with hit list, alignment time
			mode = 3
			continue
		}
		switch mode {
		case 1: // meta data
			_, err = bmeta.Write(append(line, '\n'))
		case 2: // hit list
			_, err = bhits.Write(append(line, '\n'))
		case 3: // alignments
			_, err = balign.Write(append(line, '\n'))
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

	if err = readAlignments(balign, hhr.Query, hhr.Hits); err != nil {
		return nil, fmt.Errorf("Error reading alignments from hhr: %s", err)
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

		numPart := rest[8][0 : len(rest[8])-1] // i.e., remove parens in '(52)'.
		hit.NumTemplateCols, err = strconv.Atoi(str(numPart))
		if err != nil {
			return nil, err
		}

		hits = append(hits, hit)
	}
	return hits, nil
}

func readAlignments(buf *bytes.Buffer, queryName string, hits []Hit) error {
	hi := -1 // the current hit index
	for {
		line, err := buf.ReadBytes('\n')
		if err == io.EOF && len(line) == 0 {
			break
		}
		if err != nil && err != io.EOF {
			return err
		}
		line = trim(line)
		if len(line) == 0 {
			continue
		}

		switch {
		case hasPrefix(line, "No "):
			hitno, err := strconv.Atoi(str(line[3:]))
			if err != nil {
				return err
			}
			hi = -1
			for i := range hits {
				if hits[i].Num == hitno {
					hi = i
					break
				}
			}
			if hi == -1 {
				return fmt.Errorf("Invalid hit number: %d\n", hi)
			}
			hits[hi].Aligned = Alignment{
				QSeq:  make([]seq.Residue, 0),
				QDssp: make([]seq.Residue, 0),
				QPred: make([]seq.Residue, 0),
				QConf: make([]seq.Residue, 0),
				TSeq:  make([]seq.Residue, 0),
				TDssp: make([]seq.Residue, 0),
				TPred: make([]seq.Residue, 0),
				TConf: make([]seq.Residue, 0),
			}
		case hasPrefix(line, ">"):
			if !strings.HasPrefix(str(line[1:]), hits[hi].Name) {
				return fmt.Errorf("Hit names don't match: '%s' != '%s'.",
					str(line[1:]), hits[hi].Name)
			}
		case hasPrefix(line, "Q"): // query part of alignment
			rest := line[2:]
			aligned := &hits[hi].Aligned
			switch {
			case strings.HasPrefix(queryName, str(rest[0:14])):
				aligned.QSeq = append(aligned.QSeq, getSeq(rest)...)
			case hasPrefix(rest, "Consensus"):
				aligned.QConsensus = append(aligned.QConsensus, getSeq(rest)...)
			case hasPrefix(rest, "ss_dssp"):
				aligned.QDssp = append(aligned.QDssp, getSeq(rest)...)
			case hasPrefix(rest, "ss_pred"):
				aligned.QPred = append(aligned.QPred, getSeq(rest)...)
			case hasPrefix(rest, "ss_conf"):
				aligned.QConf = append(aligned.QConf, getSeq(rest)...)
			}
		case hasPrefix(line, "T"): // template part of alignment
			rest := line[2:]
			aligned := &hits[hi].Aligned
			switch {
			case strings.HasPrefix(hits[hi].Name, str(rest[0:14])):
				aligned.TSeq = append(aligned.TSeq, getSeq(rest)...)
			case hasPrefix(rest, "Consensus"):
				aligned.TConsensus = append(aligned.TConsensus, getSeq(rest)...)
			case hasPrefix(rest, "ss_dssp"):
				aligned.TDssp = append(aligned.TDssp, getSeq(rest)...)
			case hasPrefix(rest, "ss_pred"):
				aligned.TPred = append(aligned.TPred, getSeq(rest)...)
			case hasPrefix(rest, "ss_conf"):
				aligned.TConf = append(aligned.TConf, getSeq(rest)...)
			}
		}
	}
	return nil
}

func hasPrefix(bs []byte, prefix string) bool {
	return bytes.HasPrefix(bs, []byte(prefix))
}

func trim(bs []byte) []byte {
	return bytes.TrimSpace(bs)
}

func getSeq(line []byte) []seq.Residue {
	fs := bytes.Fields(line[17:])
	var fseq []byte
	if len(fs) == 1 {
		fseq = fs[0]
	} else {
		fseq = fs[1]
	}
	rs := make([]seq.Residue, len(fseq))
	for i, r := range fseq {
		rs[i] = seq.Residue(r)
	}
	return rs
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
