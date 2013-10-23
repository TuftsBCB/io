package newick

import (
	"bufio"
	"fmt"
	"io"
	"strings"
	"unicode/utf8"
)

type itemType int

const (
	itemError itemType = iota
	itemEOF
	itemTerminal
	itemDescendentsStart
	itemDescendentsEnd
	itemSubtree
)

const (
	eof           = 0
	terminal      = ';'
	descDelimiter = ','
	descStart     = '('
	descEnd       = ')'
	quoteStart    = '\''
	quoteEnd      = '\''
	lengthStart   = ':'
)

const unquoteBanned = " ()[]':;,"

type stateFn func(lx *lexer) stateFn

type lexer struct {
	input io.Reader
	buf   string
	start int
	pos   int
	width int
	line  int
	state stateFn
	items chan item
}

type item struct {
	typ  itemType
	val  string
	line int
}

func (lx *lexer) nextItem() item {
	for {
		select {
		case item := <-lx.items:
			return item
		default:
			lx.state = lx.state(lx)
		}
	}
	panic("not reached")
}

func lex(input io.Reader) *lexer {
	lx := &lexer{
		input: bufio.NewReader(input),
		buf:   "",
		state: lexDescendents,
		line:  1,
		items: make(chan item, 10),
	}
	return lx
}

func (lx *lexer) current() string {
	return lx.buf[lx.start:lx.pos]
}

func (lx *lexer) emit(typ itemType) {
	lx.items <- item{typ, lx.current(), lx.line}
	lx.buf = lx.buf[lx.pos:]
	lx.start, lx.pos = 0, 0
}

func (lx *lexer) next() (r rune) {
	if lx.pos >= len(lx.buf) {
		buf := make([]byte, 4096)
		n, err := lx.input.Read(buf)
		if err != nil || lx.pos >= len(lx.buf)+n {
			lx.width = 0
			return eof
		}
		lx.buf += string(buf[0:n])
	}

	if lx.buf[lx.pos] == '\n' {
		lx.line++
	}
	r, lx.width = utf8.DecodeRuneInString(lx.buf[lx.pos:])
	lx.pos += lx.width
	return r
}

// ignore skips over the pending input before this point.
func (lx *lexer) ignore() {
	lx.start = lx.pos
}

// backup steps back one rune. Can be called only once per call of next.
func (lx *lexer) backup() {
	lx.pos -= lx.width
	if lx.pos < len(lx.buf) && lx.buf[lx.pos] == '\n' {
		lx.line--
	}
}

// accept consumes the next rune if it's equal to `valid`.
func (lx *lexer) accept(valid rune) bool {
	if lx.next() == valid {
		return true
	}
	lx.backup()
	return false
}

// peek returns but does not consume the next rune in the input.
func (lx *lexer) peek() rune {
	r := lx.next()
	lx.backup()
	return r
}

// errorf stops all lexing by emitting an error and returning `nil`.
// Note that any value that is a character is escaped if it's a special
// character (new lines, tabs, etc.).
func (lx *lexer) errorf(format string, values ...interface{}) stateFn {
	for i, value := range values {
		if v, ok := value.(rune); ok {
			values[i] = escapeSpecial(v)
		}
	}
	lx.items <- item{
		itemError,
		fmt.Sprintf(format, values...),
		lx.line,
	}
	return nil
}

func lexDescendents(lx *lexer) stateFn {
	r := lx.next()
	if isBlank(r) || isNL(r) {
		return lexSkip(lx, lexDescendents)
	}

	switch r {
	case descStart:
		lx.ignore()
		lx.emit(itemDescendentsStart)
		return lexSubtreeStart
	case eof:
		if lx.pos > lx.start {
			return lx.errorf("Unexpected EOF.")
		}
		lx.emit(itemEOF)
		return nil
	}
	lx.backup()
	lx.ignore()
	return lexLabel
}

func lexSubtreeStart(lx *lexer) stateFn {
	r := lx.next()
	if isBlank(r) || isNL(r) {
		return lexSkip(lx, lexSubtreeStart)
	} else if isSubtreeEnd(r) {
		lx.backup()
		lx.ignore()
		return lexSubtreeEnd
	} else if r == descStart {
		lx.backup()
		return lexDescendents
	}
	lx.backup()
	lx.ignore()
	return lexLabelStart
}

func lexSubtreeEnd(lx *lexer) stateFn {
	lx.emit(itemSubtree)
	r := lx.next()
	switch r {
	case descDelimiter:
		lx.ignore()
		return lexDescendents
	case descEnd:
		lx.ignore()
		lx.emit(itemDescendentsEnd)
		return lexLabelStart
	case terminal:
		lx.ignore()
		lx.emit(itemTerminal)
		return lexDescendents
	case eof:
		lx.ignore()
		lx.emit(itemTerminal)
		lx.emit(itemEOF)
		return nil
	}
	return lx.errorf("Expected end of subtree ('%s', '%s' or '%s') but got "+
		"'%s' instead.", descDelimiter, descEnd, terminal, r)
}

func lexLabelStart(lx *lexer) stateFn {
	r := lx.next()
	if isBlank(r) || isNL(r) {
		return lexSkip(lx, lexSubtreeStart)
	}
	lx.backup()
	return lexLabel
}

func lexLabel(lx *lexer) stateFn {
	r := lx.next()
	if r == lengthStart {
		return lexLengthStart
	} else if isSubtreeEnd(r) {
		lx.backup()
		return lexSubtreeEnd
	} else if strings.ContainsRune(unquoteBanned, r) || isNL(r) {
		return lx.errorf("Found '%s' in an unquoted label, which may not "+
			"contain the following characters: '%s'.", r, unquoteBanned)
	}
	return lexLabel
}

func lexLengthStart(lx *lexer) stateFn {
	r := lx.next()
	if r == '-' {
		return lexLengthNum
	}
	lx.backup()
	return lexLengthNum
}

func lexLengthNum(lx *lexer) stateFn {
	r := lx.next()
	if isSubtreeEnd(r) {
		lx.backup()
		return lexSubtreeEnd
	} else if r == '.' {
		return lexLengthDecimal
	} else if isDigit(r) {
		return lexLengthNum
	}
	return lx.errorf("Expected a '.', a digit or the end of a subtree, but "+
		"got '%s' instead.", r)
}

func lexLengthDecimal(lx *lexer) stateFn {
	r := lx.next()
	if isSubtreeEnd(r) {
		lx.backup()
		return lexSubtreeEnd
	} else if isDigit(r) {
		return lexLengthDecimal
	}
	return lx.errorf("Expected a digit or the end of a subtree, but "+
		"got '%s' instead.", r)
}

// lexSkip ignores all slurped input and moves on to the next state.
func lexSkip(lx *lexer, nextState stateFn) stateFn {
	return func(lx *lexer) stateFn {
		lx.ignore()
		return nextState
	}
}

func isSubtreeEnd(r rune) bool {
	return r == descDelimiter || r == descEnd || r == terminal || r == eof
}

func isBlank(r rune) bool {
	return r == '\t' || r == ' '
}

func isNL(r rune) bool {
	return r == '\n' || r == '\r'
}

func isDigit(r rune) bool {
	return r >= '0' && r <= '9'
}

func (itype itemType) String() string {
	switch itype {
	case itemError:
		return "Error"
	case itemEOF:
		return "EOF"
	case itemTerminal:
		return "Terminal"
	case itemDescendentsStart:
		return "Descendents (start)"
	case itemDescendentsEnd:
		return "Descendents (end)"
	case itemSubtree:
		return "Subtree"
	}
	panic(fmt.Sprintf("BUG: Unknown type '%s'.", itype))
}

func (item item) String() string {
	return fmt.Sprintf("(%s, %s)", item.typ.String(), item.val)
}

func escapeSpecial(c rune) string {
	switch c {
	case '\n':
		return "\\n"
	}
	return string(c)
}
