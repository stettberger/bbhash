// Package bbhash implements the BBHash algorithm for minimal perfect hash functions.
package bbhash

import (
	"fmt"
)

const (
	// Heuristic: 32 levels should be enough for even very large key sets
	initialLevels = 32

	// Maximum number of attempts (level) at making a perfect hash function.
	// Per the paper, each successive level exponentially reduces the
	// probability of collision.
	maxLevel uint = 200
)

// BBHash represents a minimal perfect hash for a set of keys.
type BBHash struct {
	bits     []*bitVector
	ranks    []uint64
	saltHash uint64 // precomputed hash of the salt
	gamma    float64
	revIndex []uint64 // reverse index: only used for reverse mapping
}

// NewSequential creates a new BBHash for the given keys. The keys must be unique.
// This creates the BBHash in a single goroutine.
// The gamma parameter is the expansion factor for the bit vector; the paper recommends
// a value of 2.0. The larger the value the more memory will be consumed by the BBHash.
// The salt parameter is used to salt the hash function. Depending on your use case,
// you may use a cryptographic- or a pseudo-random number for the salt.
func NewSequential(gamma float64, salt uint64, keys []uint64) (*BBHash, error) {
	if gamma <= 1.0 {
		gamma = 2.0
	}
	bb := &BBHash{
		bits:     make([]*bitVector, 0, initialLevels),
		saltHash: saltHash(salt),
		gamma:    gamma,
	}
	if err := bb.compute(keys); err != nil {
		return nil, err
	}
	return bb, nil
}

// Find returns a unique index representing the key in the minimal hash set.
//
// The return value is meaningful ONLY for keys in the original key set
// (provided at the time of construction of the minimal hash set).
//
// If the key is in the original key set, the return value is guaranteed to be
// in the range [1, len(keys)].
//
// If the key is not in the original key set, two things can happen:
// 1. The return value is 0, representing that the key was not in the original key set.
// 2. The return value is in the expected range [1, len(keys)], but is a false positive.
func (bb *BBHash) Find(key uint64) uint64 {
	for lvl, bv := range bb.bits {
		i := hash(bb.saltHash, uint64(lvl), key) % bv.Size()
		if bv.IsSet(i) {
			return bb.ranks[lvl] + bv.Rank(i)
		}
	}
	return 0
}

// compute computes the minimal perfect hash for the given keys.
func (bb *BBHash) compute(keys []uint64) error {
	sz := uint64(len(keys))
	ld := newLevelData(sz, bb.gamma)

	// loop exits when keys == nil, i.e., when there are no more keys to re-hash
	for lvl := uint(0); keys != nil; lvl++ {
		sz = ld.current.Size()
		// precompute the level hash to speed up the key hashing
		lvlHash := levelHash(bb.saltHash, uint64(lvl))

		// find colliding keys and possible bit vector positions for non-colliding keys
		for _, k := range keys {
			i := keyHash(lvlHash, k) % sz
			if ld.current.IsSet(i) {
				// found one or more collisions at index i
				ld.collisions.Set(i)
				continue
			}
			ld.current.Set(i)
		}

		// remove bit vector position assignments for colliding keys and add them to the redo set
		for _, k := range keys {
			i := keyHash(lvlHash, k) % sz
			if ld.collisions.IsSet(i) {
				ld.redo = append(ld.redo, k)
				// unset the bit since there was a collision
				ld.current.Unset(i)
			}
		}

		// save the current bit vector for the current level
		bb.bits = append(bb.bits, ld.current)
		// move to next level and compute the set of keys to re-hash (that had collisions)
		keys = ld.nextLevel()

		if lvl > maxLevel {
			return fmt.Errorf("can't find minimal perfect hash after %d tries", lvl)
		}
	}
	bb.computeLevelRanks()
	return nil
}

// lvlData holds intermediate results for the current level
type lvlData struct {
	current    *bitVector // bit vector for current level  : A in the paper
	collisions *bitVector // collisions at current level   : C in the paper
	redo       []uint64   // keys to re-hash at next level : F in the paper
	gamma      float64
}

func newLevelData(sz uint64, gamma float64) *lvlData {
	words := words(sz, gamma)
	return &lvlData{
		gamma:      gamma,
		current:    newBitVector(words),
		collisions: newBitVector(words),
		redo:       make([]uint64, 0, sz/2), // heuristic: only 1/2 of the keys will collide
	}
}

// nextLevel moves to the next level and returns the keys that need to be re-hashed.
func (ld *lvlData) nextLevel() []uint64 {
	remainingKeys := ld.redo
	numKeys := uint64(len(remainingKeys))
	if numKeys == 0 {
		return nil
	}
	// Reset redo set for the next level, reusing the slice for the next level's keys.
	// It is possible to reuse this slice since we are adding new keys to the redo slice
	// after the corresponding keys have been processed from the current level's "keys" slice.
	// That is, the redo slice is always a subset of the keys slice, and the redo keys' indices
	// are trailing the indices of the keys slice.
	ld.redo = ld.redo[:0:numKeys]
	// number of words for the next level's bit vector
	words := words(numKeys, ld.gamma)
	ld.current = newBitVector(words)
	ld.collisions.Reset(words)
	return remainingKeys
}

// computeLevelRanks computes the total rank of each level.
// The total rank is the rank for all levels up to and including the current level.
func (bb *BBHash) computeLevelRanks() {
	// Initializing the rank to 1, since the 0 index is reserved for not-found.
	var rank uint64 = 1
	bb.ranks = make([]uint64, len(bb.bits))
	for l, bv := range bb.bits {
		bb.ranks[l] = rank
		rank += bv.OnesCount()
	}
}
