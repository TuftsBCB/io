package hhm

import (
	"bufio"
	"fmt"
	"io"
	"strings"

	"github.com/BurntSushi/bcbgo/io/fasta"
	"github.com/BurntSushi/bcbgo/io/msa"
	"github.com/BurntSushi/bcbgo/seq"
)

func Write(w io.Writer, hhm *HHM) error {
	buf := bufio.NewWriter(w)
	if err := writeMeta(buf, hhm); err != nil {
		return err
	}
	if _, err := buf.WriteString("SEQ\n"); err != nil {
		return err
	}
	if err := writeSecondary(buf, hhm); err != nil {
		return err
	}
	if err := writeMSA(buf, hhm); err != nil {
		return err
	}
	if _, err := buf.WriteString("#\n"); err != nil {
		return err
	}
	if err := writeHMM(buf, hhm); err != nil {
		return err
	}
	if _, err := buf.WriteString("//\n"); err != nil {
		return err
	}
	return buf.Flush()
}

func writeMeta(buf *bufio.Writer, hhm *HHM) (err error) {
	defer func() {
		if r := recover(); r != nil {
			var ok bool
			if err, ok = r.(error); ok {
				return
			}
			panic(r)
		}
	}()
	w := func(format string, v ...interface{}) {
		s := fmt.Sprintf(format+"\n", v...)
		if _, err = buf.WriteString(s); err != nil {
			panic(err)
		}
	}

	meta := hhm.Meta
	if len(meta.FormatVersion) > 0 {
		w("%s", meta.FormatVersion)
	}
	if len(meta.Name) > 0 {
		w("NAME  %s", meta.Name)
	}
	if len(meta.Fam) > 0 {
		w("FAM   %s", meta.Fam)
	}
	if len(meta.File) > 0 {
		w("FILE  %s", meta.File)
	}
	if len(meta.Com) > 0 {
		w("COM   %s", meta.Com)
	}
	if len(meta.Date) > 0 {
		w("DATE  %s", meta.Date)
	}
	if len(meta.Leng) > 0 {
		// Not sure if this is right. ???
		m := len(hhm.HMM.Nodes)
		w("LENG  %d match states, %d columns in multiple alignment", m, m)
	}
	if len(meta.Filt) > 0 {
		w("FILT  %s", meta.Filt)
	}
	w("NEFF  %f", meta.Neff)
	if meta.EvdLambda != 0 || meta.EvdMu != 0 {
		w("EVD   %0.4f  %0.4f", meta.EvdLambda, meta.EvdMu)
	}
	if meta.Pct {
		w("PCT   True")
	} else {
		w("PCT   False")
	}
	if len(meta.Desc) > 0 {
		w("DESC  %s", meta.Desc)
	}
	return nil
}

func writeSecondary(buf *bufio.Writer, hhm *HHM) error {
	ss := hhm.Secondary
	towrite := make([]seq.Sequence, 0, 5)
	if ss.SSdssp != nil {
		towrite = append(towrite, *ss.SSdssp)
	}
	if ss.SAdssp != nil {
		towrite = append(towrite, *ss.SAdssp)
	}
	if ss.SSpred != nil {
		towrite = append(towrite, *ss.SSpred)
	}
	if ss.SSconf != nil {
		towrite = append(towrite, *ss.SSconf)
	}
	if ss.Consensus != nil {
		towrite = append(towrite, *ss.Consensus)
	}
	w := fasta.NewWriter(buf)
	w.Asterisk = false
	w.Columns = 0
	return w.WriteAll(towrite)
}

func writeMSA(buf *bufio.Writer, hhm *HHM) error {
	return msa.WriteA3M(buf, hhm.MSA)
}

func writeHMM(buf *bufio.Writer, hhm *HHM) (err error) {
	defer func() {
		if r := recover(); r != nil {
			var ok bool
			if err, ok = r.(error); ok {
				return
			}
			panic(r)
		}
	}()
	w := func(format string, v ...interface{}) {
		s := fmt.Sprintf(format, v...)
		if _, err = buf.WriteString(s); err != nil {
			panic(err)
		}
	}

	hmm := hhm.HMM
	w("NULL   ")
	must(writeEmissions(buf, hmm.Alphabet, hmm.Null))
	w("\n")

	w("HMM    ")
	must(writeAlphabet(buf, hmm.Alphabet))
	w("\n")
	w("       M->M\tM->I\tM->D\tI->M\tI->I\tD->M\tD->D\tNeff\tNeff_I\tNeff_D\n")
	w("       0\t*\t0\t*\t*\t*\t*\t*\t*\t*\n")

	for _, node := range hmm.Nodes {
		w("%c %d  ", node.Residue, node.NodeNum)
		must(writeEmissions(buf, hmm.Alphabet, node.MatEmit))
		w("  %d\n", node.NodeNum)
		w("   ")
		must(writeTransitions(buf, node.Transitions))
		w("\t")
		must(writeDiversity(buf, node.NeffM, node.NeffI, node.NeffD))
		w("\n\n")
	}

	return nil
}

func writeAlphabet(buf *bufio.Writer, alphabet []seq.Residue) error {
	alpha := make([]string, len(alphabet))
	for i, residue := range alphabet {
		alpha[i] = fmt.Sprintf("%c", residue)
	}
	_, err := buf.WriteString(strings.Join(alpha, "\t"))
	return err
}

func writeEmissions(
	buf *bufio.Writer, alphabet []seq.Residue, ep seq.EProbs) error {

	probs := make([]string, len(alphabet))
	for i, residue := range alphabet {
		probs[i] = probStr(ep[residue])
	}
	_, err := buf.WriteString(strings.Join(probs, "\t"))
	return err
}

func writeTransitions(buf *bufio.Writer, tp seq.TProbs) error {
	towrite := []string{
		probStr(tp.MM),
		probStr(tp.MI),
		probStr(tp.MD),
		probStr(tp.IM),
		probStr(tp.II),
		probStr(tp.DM),
		probStr(tp.DD),
	}
	_, err := buf.WriteString(strings.Join(towrite, "\t"))
	return err
}

func writeDiversity(buf *bufio.Writer, m, i, d seq.Prob) error {
	towrite := []string{
		neffStr(m),
		neffStr(i),
		neffStr(d),
	}
	_, err := buf.WriteString(strings.Join(towrite, "\t"))
	return err
}

func probStr(p seq.Prob) string {
	if p.IsMin() {
		return "*"
	}
	scaled := int(-seq.Prob(hmmScale) * p)
	return fmt.Sprintf("%d", scaled)
}

func neffStr(p seq.Prob) string {
	scaled := int(seq.Prob(hmmScale) * p)
	return fmt.Sprintf("%d", scaled)
}

func must(err error) {
	if err != nil {
		panic(err)
	}
}
