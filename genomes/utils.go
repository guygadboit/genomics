package genomes

func ReverseComplement(nts []byte) []byte {
	ret := make([]byte, len(nts))

	for i := 0; i < len(nts); i++ {

		var nt byte
		switch nts[i] {
		case 'G':
			nt = 'C'
		case 'C':
			nt = 'G'
		case 'A':
			nt = 'T'
		case 'T':
			nt = 'A'
		}

		j := len(nts) - i - 1
		ret[j] = nt
	}
	return ret
}
