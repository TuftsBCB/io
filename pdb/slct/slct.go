package slct

import (
	"bufio"
	"fmt"
	"io"
	"strconv"
	"strings"
)

type Reader struct {
	*bufio.Reader
}

type Entry struct {
	Threshold int
	ChainID   string
	NAA       int
	Res       float64
	Rfac      float64
	Method    string
	NSid      int
	NBck      int
	NNAA      int
	NHlx      int
	NBta      int
	Compound  string
}

func NewReader(r io.Reader) *Reader {
	return &Reader{bufio.NewReader(r)}
}

func (r *Reader) Read() (*Entry, error) {
	for {
		line, err := r.ReadString('\n')
		if err != nil && err != io.EOF {
			return nil, err
		}
		line = strings.TrimSpace(line)
		if len(line) == 0 && err == io.EOF {
			return nil, io.EOF
		}

		if len(line) == 0 || line[0] == '#' {
			continue
		}
		return toEntry(strings.Fields(line))
	}
}

func (r *Reader) ReadAll() ([]*Entry, error) {
	var entries []*Entry
	for {
		entry, err := r.Read()
		if err == io.EOF {
			return entries, nil
		}
		if err != nil {
			return nil, err
		}
		entries = append(entries, entry)
	}
	return entries, nil
}

func toEntry(record []string) (*Entry, error) {
	if len(record) < 12 {
		return nil, fmt.Errorf(
			"Record has fewer than 12 fields. Not a valid PDB select 25 file.")
	}

	var err error
	atoi := strconv.Atoi
	atof := func(s string) (float64, error) { return strconv.ParseFloat(s, 64) }
	e := new(Entry)

	e.ChainID = record[1]
	e.Method = record[5]
	e.Compound = strings.Join(record[11:], " ")
	if e.Compound[0] == '\'' && e.Compound[len(e.Compound)-1] == '\'' {
		e.Compound = strings.TrimSpace(e.Compound[1 : len(e.Compound)-1])
	}

	if e.Threshold, err = atoi(record[0]); err != nil {
		return nil, err
	}
	if e.NAA, err = atoi(record[2]); err != nil {
		return nil, err
	}
	if e.Res, err = atof(record[3]); err != nil {
		return nil, err
	}
	if e.Rfac, err = atof(record[4]); err != nil {
		return nil, err
	}
	if e.NSid, err = atoi(record[6]); err != nil {
		return nil, err
	}
	if e.NBck, err = atoi(record[7]); err != nil {
		return nil, err
	}
	if e.NNAA, err = atoi(record[8]); err != nil {
		return nil, err
	}
	if e.NHlx, err = atoi(record[9]); err != nil {
		return nil, err
	}
	if e.NBta, err = atoi(record[10]); err != nil {
		return nil, err
	}
	return e, nil
}
