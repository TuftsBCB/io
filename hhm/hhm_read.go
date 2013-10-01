package hhm

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"strconv"
	"strings"

	"github.com/TuftsBCB/io/fasta"
	"github.com/TuftsBCB/seq"
)

func Read(r io.Reader) (*HHM, error) {
	// An hhm file as four logical sections: 1) Meta data, 2) secondary
	// structure info (optional), 3) A2M formatted MSA and 4) the HMM.
	// We group 2+3 together, and store each of the three portions in their
	// own buffer. We then parse each separately.
	bmeta, bseq, bhmm := new(bytes.Buffer), new(bytes.Buffer), new(bytes.Buffer)
	mode := 1 // 1 for meta, 2 for sequence and 3 for hmm

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

		// First check if we should do a mode change.
		// i.e., a line starting with 'SEQ' means mode 1 -> 2, and
		// a line starting with '#' means mode 2 -> 3.
		// Finally, a line starting with '\\' means STOP.
		switch {
		case hasPrefix(line, "SEQ"):
			mode = 2
			continue MAIN
		case hasPrefix(line, "#"):
			mode = 3
			continue MAIN
		case hasPrefix(line, "//"):
			break MAIN
		}

		// Now add the bytes to the appropriate buffer based on the mode.
		switch mode {
		case 1: // meta data
			_, err = bmeta.Write(append(line, '\n'))
		case 2: // sequences
			_, err = bseq.Write(append(line, '\n'))
		case 3: // hmm
			_, err = bhmm.Write(append(line, '\n'))
		default:
			panic(fmt.Sprintf("BUG: Unknown mode: %d", mode))
		}
		if err != nil {
			return nil, fmt.Errorf("Error updating buffer for hhm: %s", err)
		}
	}

	meta, err := readMeta(bmeta)
	if err != nil {
		return nil, fmt.Errorf("Error reading meta data from hhm: %s", err)
	}

	ss, msa, err := readSeqs(bseq)
	if err != nil {
		return nil, fmt.Errorf("Error reading sequence data from hhm: %s", err)
	}

	hmm, err := readHMM(bhmm)
	if err != nil {
		return nil, fmt.Errorf("Error reading HMM data from hhm: %s", err)
	}

	return &HHM{
		Meta:      meta,
		Secondary: ss,
		MSA:       msa,
		HMM:       hmm,
	}, nil
}

func readMeta(buf *bytes.Buffer) (HHMMeta, error) {
	meta := HHMMeta{}
	for {
		line, err := buf.ReadBytes('\n')
		if err == io.EOF && len(line) == 0 {
			break
		}
		if err != nil && err != io.EOF {
			return HHMMeta{}, err
		}

		line = trim(line)
		switch {
		case hasPrefix(line, "HH"):
			meta.FormatVersion = str(line)
		case hasPrefix(line, "NAME"):
			meta.Name = str(line[4:])
		case hasPrefix(line, "FAM"):
			meta.Fam = str(line[3:])
		case hasPrefix(line, "FILE"):
			meta.File = str(line[4:])
		case hasPrefix(line, "LENG"):
			meta.Leng = str(line[4:])
		case hasPrefix(line, "FILT"):
			meta.Filt = str(line[4:])
		case hasPrefix(line, "NEFF"):
			// You'd think we could use readNeff here, but does the HHM
			// format store all Neff values equally? NOOOOOOOOOOOOOOOOOOOO.
			f, err := strconv.ParseFloat(str(line[4:]), 64)
			if err != nil {
				return HHMMeta{}, err
			}
			meta.Neff = seq.Prob(f)
		case hasPrefix(line, "EVD"):
			fields := bytes.Fields(bytes.TrimSpace(line[3:]))
			if len(fields) != 2 {
				return HHMMeta{}, fmt.Errorf("Invalid EVD format: '%s'", line)
			}

			lambda, err := strconv.ParseFloat(string(fields[0]), 64)
			if err != nil {
				return HHMMeta{}, fmt.Errorf("Error EVD lambda '%s': %s",
					string(fields[0]), err)
			}
			meta.EvdLambda = lambda

			mu, err := strconv.ParseFloat(string(fields[1]), 64)
			if err != nil {
				return HHMMeta{}, fmt.Errorf("Error EVD mu '%s': %s",
					string(fields[1]), err)
			}
			meta.EvdMu = mu
		case hasPrefix(line, "PCT"):
			meta.Pct = true
		case hasPrefix(line, "DESC"):
			meta.Desc = str(line[4:])
		case hasPrefix(line, "COM"):
			meta.Com = str(line[3:])
		case hasPrefix(line, "DATE"):
			meta.Date = str(line[4:])
		}
	}
	return meta, nil
}

func readSeqs(buf *bytes.Buffer) (HHMSecondary, seq.MSA, error) {
	// Remember, the sequence portion of an HHM file actually has two parts.
	// The first part is optional and contains secondary structure information.
	// These SS sequences can be identified by special sequence headers:
	// "ss_dssp", "sa_dssp", "ss_pred", "ss_conf", and "Consensus".
	// If a sequence doesn't contain a special header, then that signifies that
	// we should start reading the MSA, which comes after the SS information.
	doneSS := false

	ss := HHMSecondary{}
	msa := seq.NewMSA()

	reader := fasta.NewReader(buf)
	reader.TrustSequences = true
	seqs, err := reader.ReadAll()
	if err != nil {
		return HHMSecondary{}, seq.MSA{}, err
	}

	for _, s := range seqs {
		s := s
		if !doneSS {
			switch {
			case strings.HasPrefix(s.Name, "ss_dssp"):
				ss.SSdssp = &s
			case strings.HasPrefix(s.Name, "sa_dssp"):
				ss.SAdssp = &s
			case strings.HasPrefix(s.Name, "ss_pred"):
				ss.SSpred = &s
			case strings.HasPrefix(s.Name, "ss_conf"):
				ss.SSconf = &s
			case strings.HasPrefix(s.Name, "Consensus"):
				ss.Consensus = &s
			default:
				doneSS = true
			}
		}
		if doneSS {
			msa.Add(s)
		}
	}
	return ss, msa, nil
}

func readHMM(buf *bytes.Buffer) (hmm *seq.HMM, err error) {
	var nullFields []string
	hmm = new(seq.HMM)
	for {
		line, err := buf.ReadBytes('\n')
		if err == io.EOF && len(line) == 0 {
			break
		}
		if err != nil && err != io.EOF {
			panic(err)
		}
		line = trim(line)

		if len(line) == 0 {
			continue
		}
		switch {
		case hasPrefix(line, "NULL"):
			// We can't read the NULL emissions yet, because we don't have
			// an alphabet. (Which we'll get on the next line.)
			// We'll slurp this into a seq.EProbs value in a little bit, after
			// we get an alphabet.
			nullFields = strings.Fields(str(line[4:]))
		case hasPrefix(line, "HMM"):
			// We slurp up three lines here. The first is the alphabet
			// (the current line). The second is the ordering of transition
			// probabilities. And the third are transition probabilities for
			// the begin state. We ignore the second two.
			if _, err := demandLine(buf); err != nil {
				return nil, fmt.Errorf("%s (expected transition ordering)", err)
			}
			if _, err := demandLine(buf); err != nil {
				return nil, fmt.Errorf("%s (expected start transitions)", err)
			}

			// Get the ordering of the alphabet.
			hmm.Alphabet = make([]seq.Residue, 0, 20)
			residues := bytes.Split(trim(line[3:]), []byte{'\t'})
			for _, residue := range residues {
				hmm.Alphabet = append(hmm.Alphabet, seq.Residue(residue[0]))
			}

			// Remember those null probabilities? Well, we have an alphabet now.
			ep, err := readEmissions(hmm.Alphabet, nullFields)
			if err != nil {
				return nil, fmt.Errorf("Could not read NULL emissions '%s': %s",
					strings.Join(nullFields, " "), err)
			}
			hmm.Null = ep
		default: // finally, reading a node in the HMM
			// Each node in the HMM is made up of two lines.
			// The first line starts with the amino acid in the reference
			// sequence, followed by the node number, followed by N match
			// emission probabilities (where N = len(alphabet)), and finally
			// followed by the node number again. (What the *fuck* is up with
			// that? Seriously.)
			//
			// The second line is made up of 7 transition probabilities,
			// followed by 3 diversity (the 'neff' stuff) scores.
			//
			// Also, each field is separated by spaces OR tabs. Lovely, eh?
			line2, err := demandLine(buf)
			if err != nil {
				return nil, fmt.Errorf("%s (expected transition probs)", err)
			}
			fields1 := strings.Fields(string(line))
			fields2 := strings.Fields(string(line2))
			node := seq.HMMNode{
				Residue: seq.Residue(fields1[0][0]),
			}

			node.NodeNum, err = strconv.Atoi(fields1[1])
			if err != nil {
				return nil, fmt.Errorf("Could not parse node number '%s': %s",
					fields1[1], err)
			}

			node.MatEmit, err = readEmissions(hmm.Alphabet, fields1[2:])
			if err != nil {
				return nil, fmt.Errorf("Could not read emissions '%s': %s",
					strings.Join(fields1[2:], " "), err)
			}

			node.InsEmit = seq.NewEProbs(hmm.Alphabet)
			for residue := range node.InsEmit {
				node.InsEmit[residue] = hmm.Null[residue]
			}

			node.Transitions, err = readTransitions(fields2)
			if err != nil {
				return nil, fmt.Errorf("Could not read transitions '%s': %s",
					strings.Join(fields2, " "), err)
			}

			node.NeffM, node.NeffI, node.NeffD, err = readDiversity(fields2[7:])
			if err != nil {
				return nil, fmt.Errorf("Could not read diversity '%s': %s",
					strings.Join(fields2[7:], " "), err)
			}

			hmm.Nodes = append(hmm.Nodes, node)
		}
	}
	return hmm, nil
}

func readEmissions(alphabet []seq.Residue, flds []string) (seq.EProbs, error) {
	var err error

	ep := seq.NewEProbs(alphabet)
	for i := 0; i < len(alphabet); i++ {
		if ep[alphabet[i]], err = readProb(flds[i]); err != nil {
			return nil, err
		}
	}
	return ep, nil
}

func readTransitions(fields []string) (tp seq.TProbs, err error) {
	if tp.MM, err = readProb(fields[0]); err != nil {
		return
	}
	if tp.MI, err = readProb(fields[1]); err != nil {
		return
	}
	if tp.MD, err = readProb(fields[2]); err != nil {
		return
	}
	if tp.IM, err = readProb(fields[3]); err != nil {
		return
	}
	if tp.II, err = readProb(fields[4]); err != nil {
		return
	}
	if tp.DM, err = readProb(fields[5]); err != nil {
		return
	}
	if tp.DD, err = readProb(fields[6]); err != nil {
		return
	}
	return
}

func readDiversity(fields []string) (m, i, d seq.Prob, err error) {
	if m, err = readNeff(fields[0]); err != nil {
		return
	}
	if i, err = readNeff(fields[1]); err != nil {
		return
	}
	if d, err = readNeff(fields[2]); err != nil {
		return
	}
	return
}

func demandLine(buf *bytes.Buffer) ([]byte, error) {
	line, err := buf.ReadBytes('\n')
	if err == io.EOF {
		return nil, fmt.Errorf("Unexpected EOF when reading HMM in hhm.")
	}
	if err != nil {
		return nil, err
	}
	return trim(line), nil
}

// readProb reads a probability (transition or emissions) from an hhm file and
// returns a Prob value in log_2 form.
func readProb(fstr string) (seq.Prob, error) {
	f, err := seq.NewProb(fstr)
	if err != nil {
		return f, fmt.Errorf("Error reading probability '%s': %s", fstr, err)
	}
	if f.IsMin() {
		return f, nil
	}
	return -f / seq.Prob(hmmScale), nil
}

// readNeff reads a diversity value.
func readNeff(fstr string) (seq.Prob, error) {
	f, err := seq.NewProb(fstr)
	if err != nil {
		return f, fmt.Errorf("Error reading neff '%s': %s", fstr, err)
	}
	if f.IsMin() {
		return f, nil
	}
	return f / seq.Prob(hmmScale), nil
}

func hasPrefix(bs []byte, prefix string) bool {
	return bytes.HasPrefix(bs, []byte(prefix))
}

func trim(bs []byte) []byte {
	return bytes.TrimSpace(bs)
}

func str(bs []byte) string {
	return string(bytes.TrimSpace(bs))
}
