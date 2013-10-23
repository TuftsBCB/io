package newick

import (
	"fmt"
	"io"
	"strconv"
	"strings"
)

// Reader corresponds to the state necessary to read trees from Newick
// formatted input.
type Reader struct {
	*lexer
}

// NewReader returns a reader ready for reading trees from `r`.
func NewReader(r io.Reader) *Reader {
	return &Reader{lex(r)}
}

// ReadAll returns all of the Newick trees in the source input. The first
// error that occurs is returned with no trees. The error is never `io.EOF`.
func (lx *Reader) ReadAll() ([]*Tree, error) {
	trees := make([]*Tree, 0)
	for {
		tree, err := lx.ReadTree()
		if err == io.EOF {
			break
		} else if err != nil {
			return nil, err
		}
		trees = append(trees, tree)
	}
	return trees, nil
}

// ReadTree reads a single tree from the source input. If the end of the
// input is reached, then a nil `Tree` is returned with `io.EOF` as the error.
func (lx *Reader) ReadTree() (*Tree, error) {
	item := lx.nextItem()
	parent := &Tree{}
	if item.typ == itemTerminal {
		return parent, nil
	} else if item.typ == itemEOF {
		return nil, io.EOF
	}

	if err := lx.parse(parent, item); err != nil {
		return nil, err
	}

	item = lx.nextItem()
	if item.typ != itemTerminal {
		return nil, expectErr(item, fmt.Sprintf("a terminal '%c'", terminal))
	}
	return parent, nil
}

func (lx *Reader) parse(parent *Tree, next item) error {
	switch next.typ {
	case itemSubtree:
		if err := setLabelLength(parent, next.val); err != nil {
			return errf(next.line, err.Error())
		}
		return nil
	case itemDescendentsStart:
		// good to go!
	default:
		return expectErr(next, "a descendent list or a subtree")
	}

	// If we're here, then we're starting a descendent list.
	// Now we should expected one or more subtrees or descendent lists.
TOKENS:
	for {
		item := lx.nextItem()
		switch item.typ {
		case itemSubtree:
			child := &Tree{}
			if err := setLabelLength(child, item.val); err != nil {
				return errf(item.line, err.Error())
			}
			parent.Children = append(parent.Children, *child)
		case itemDescendentsStart:
			child := &Tree{}
			if err := lx.parse(child, item); err != nil {
				return err
			}
			parent.Children = append(parent.Children, *child)
		case itemDescendentsEnd:
			break TOKENS
		default:
			return expectErr(item, "a descendent list or a subtree")
		}
	}

	// After a descendent list is done, we should always expect a subtree.
	item := lx.nextItem()
	if item.typ != itemSubtree {
		return expectErr(item, "a subtree")
	}
	if err := setLabelLength(parent, item.val); err != nil {
		return errf(item.line, err.Error())
	}
	return nil
}

func setLabelLength(t *Tree, label string) error {
	label = strings.TrimSpace(label)
	if len(label) == 0 {
		return nil
	}
	pieces := strings.SplitN(label, ":", 2)
	t.Label = pieces[0]

	if len(pieces) == 2 {
		length, err := strconv.ParseFloat(pieces[1], 64)
		if err != nil {
			return fmt.Errorf("Invalid branch length: %s", err)
		}
		t.Length = &length
	}
	return nil
}

func expectErr(item item, expected string) error {
	if item.typ == itemError {
		return errf(item.line, item.val)
	}
	return errf(item.line, "Unexpected %s, expected %s.", item.typ, expected)
}

func errf(line int, format string, v ...interface{}) error {
	return fmt.Errorf("Error on line %d: %s", line, fmt.Sprintf(format, v...))
}
