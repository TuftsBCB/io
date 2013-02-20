package pdb

import (
	"fmt"
	"strings"
)

func (e *Entry) IdString() string {
	return strings.ToLower(e.IdCode)
}

func (e *Entry) AtomChunks() [][]Coords {
	chunks := make([][]Coords, len(e.Chains))
	for i, chain := range e.Chains {
		chunks[i] = chain.CaAtoms()
	}
	return chunks
}

func (chain *Chain) IdString() string {
	switch {
	case len(chain.Entry.scop) > 0:
		return chain.Entry.scop
	case len(chain.Entry.cath) > 0:
		return chain.Entry.cath
	}
	return fmt.Sprintf("%s%c", strings.ToLower(chain.Entry.IdCode), chain.Ident)
}

func (chain *Chain) AtomChunks() [][]Coords {
	return [][]Coords{chain.CaAtoms()}
}
