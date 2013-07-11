package pdb

import (
	"fmt"
	"strings"

	"github.com/TuftsBCB/structure"
)

func (e *Entry) Id() string {
	return strings.ToLower(e.IdCode)
}

func (e *Entry) Data() string {
	return ""
}

func (e *Entry) Atoms() [][]structure.Coords {
	chunks := make([][]structure.Coords, len(e.Chains))
	for i, chain := range e.Chains {
		chunks[i] = chain.CaAtoms()
	}
	return chunks
}

func (chain *Chain) Id() string {
	switch {
	case len(chain.Entry.scop) > 0:
		return chain.Entry.scop
	case len(chain.Entry.cath) > 0:
		return chain.Entry.cath
	}
	return fmt.Sprintf("%s%c", strings.ToLower(chain.Entry.IdCode), chain.Ident)
}

func (chain *Chain) Data() string {
	return ""
}

func (chain *Chain) Atoms() [][]structure.Coords {
	return [][]structure.Coords{chain.CaAtoms()}
}
