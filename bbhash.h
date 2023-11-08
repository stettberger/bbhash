#ifndef __BBHASH_H
#define __BBHASH_H

#include <stdint.h>

#define CACHELINE_SZ 64
#define BCVEC_BITS_PER_ENTRY (CACHELINE_SZ / 2 * 8)

#define ARRAY_SIZE(a) (sizeof((a))/sizeof(*(a)))

struct bbhash_colvec {
  uint64_t v[CACHELINE_SZ / 2 / sizeof(uint64_t)];
  uint64_t c[CACHELINE_SZ / 2 / sizeof(uint64_t)];
};


struct bbhash_level {
    union {
        // The collision-vector representation is only used during
        // construction.
        struct bbhash_colvec *col_vec;
        uint64_t *v;
    };
    uint64_t mmap_size;

    uint64_t  elems;
    uint64_t  bits;
    uint64_t  rank;

};

struct bbhash {
    struct bbhash_level **levels;
    int level_count;
};

struct bbhash BBHashCompute(double gamma, uintptr_t _keys, unsigned long key_count);

uint64_t * BBHashGetLevel(struct bbhash, int lvl, uint64_t *bits, uint64_t *rank);

void BBHashFree(struct bbhash);

#endif
