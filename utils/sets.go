package utils

func ToSet[S comparable](s []S) map[S]bool {
	ret := make(map[S]bool)
	for _, item := range s {
		ret[item] = true
	}
	return ret
}

func FromSet[T comparable](s map[T]bool) []T {
	ret := make([]T, 0, len(s))
	for k, _ := range s {
		ret = append(ret, k)
	}
	return ret
}

// a gets b unioned in with it
func Union[T comparable](a map[T]bool, b map[T]bool) {
	for k, _ := range b {
		a[k] = true
	}
}

// Intersection of two sets.
func Intersection[T comparable](a map[T]bool, b map[T]bool) map[T]bool {
	ret := make(map[T]bool)
	for k, _ := range a {
		if b[k] {
			ret[k] = true
		}
	}
	return ret
}

// Return the set a-b
func Difference[T comparable](a map[T]bool, b map[T]bool) map[T]bool {
	ret := make(map[T]bool)
	for k, _ := range a {
		if !b[k] {
			ret[k] = true
		}
	}
	return ret
}
