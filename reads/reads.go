package reads

import (
	"strings"
	"fmt"
	"genomics/utils"
	"bufio"
)

type Read struct {
	Name	string
	Guff	string
	Nts     []byte
	Quality []byte
}

type ReadMsg struct {
	Read
	End	bool
}

func (r *Read) Output(w *bufio.Writer) {
	fmt.Fprintf(w, "@%s\n", r.Guff)
	fmt.Fprintf(w, "%s\n", string(r.Nts))
	fmt.Fprintf(w, "+%s\n", r.Guff)
	fmt.Fprintf(w, "%s\n", string(r.Quality))
}

func ParseFastq(fname string, output chan ReadMsg) error {
	var msg ReadMsg
	var destination *[]byte
	var post bool
	var ret error

	utils.Lines(fname, func(line string, err error) bool {
		if err != nil {
			ret = err
			return false
		}
		if line[0] == '@' {
			fields := strings.Split(line, " ")
			msg.Read.Name = fields[0][1:]
			msg.Read.Guff = line[1:]
			destination = &msg.Read.Nts
		} else if line[0] == '+' {
			destination = &msg.Read.Quality
			post = true
		} else {
			*destination = []byte(line)
			if post {
				output <- msg
				post = false
			}
		}
		return true
	})
	if ret != nil {
		return ret
	}
	output <- ReadMsg{Read{"", "", nil, nil}, true}
	return nil
}
