package newick

import (
	"bytes"
	"fmt"
	"strings"
)

// Tree corresponds to any value representable in a Newick format. Each
// tree value corresponds to a single node.
type Tree struct {
	// All children of this node, which may be empty.
	Children []Tree

	// The label of this node. If it's empty, then this node does
	// not have a name.
	Label string

	// The branch length of this node corresponding to the distance between
	// it and its parent node. If it's `nil`, then no distance exists.
	Length *float64
}

// String recursively converts a tree to a string, with whitespace indenting
// to indicate depth.
func (tree *Tree) String() string {
	buf := new(bytes.Buffer)
	pf := func(format string, v ...interface{}) {
		fmt.Fprintf(buf, format, v...)
	}

	var out func(t *Tree, depth int)
	out = func(t *Tree, depth int) {
		name, length := t.Label, ""
		if len(name) == 0 {
			name = "N/A"
		}
		if t.Length != nil {
			length = fmt.Sprintf(" (%f)", *t.Length)
		}
		pf("%s%s%s\n", strings.Repeat("  ", depth), name, length)
		for _, child := range t.Children {
			out(&child, depth+1)
		}
	}
	out(tree, 0)
	return buf.String()
}
