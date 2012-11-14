package hhm

import (
	"io"

	"github.com/BurntSushi/bcbgo/seq"
)

type HHM struct {
	HMM *seq.HMM
}

func Read(r io.Reader) HHM {
	return HHM{}
}
