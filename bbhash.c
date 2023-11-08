#include "bbhash.h"
#include <sys/mman.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <sys/param.h>
#include <linux/perf_event.h>

#include <string.h>
#include <sys/syscall.h>
#include <sys/ioctl.h>

#define die(msg) do { perror(msg); exit(EXIT_FAILURE); } while(0)
#define likely(x)    __builtin_expect (!!(x), 1)
#define unlikely(x)  __builtin_expect (!!(x), 0)


const uint64_t m = 0x880355f21e6d1965;

// mix is a compression function for fast hashing.
static inline uint64_t mix(uint64_t h) {
  h ^= h >> 23;
  h *= 0x2127599bf4325c37;
  h ^= h >> 47;
  return h;
}

// levelHash returns the hash of the given level.
static inline uint64_t levelHash(uint64_t level) { return mix(level) * m; }

// keyHash returns the hash of a key given a level hash.
static inline uint64_t keyHash(uint64_t levelHash, uint64_t key) {
  uint64_t h = levelHash;
  h ^= mix(key);
  h *= m;
  h = mix(h);
  return h;
}

// hash returns the hash of the current level and key.
static inline uint64_t hash(uint64_t level, uint64_t key) {
  return keyHash(levelHash(level), key);
}

// return word %  p;
static inline uint64_t fastrange64(uint64_t word, uint64_t p) {
    // This "range"  yields a different result than %
    // return (uint64_t)(((__uint128_t)word * (__uint128_t)p) >> 64);
    return word % p;
}

static struct bbhash_level *bbhash_level_alloc(uint64_t bits) {
    struct bbhash_level *lvl = malloc(sizeof(struct bbhash_level));
    lvl->bits      = (bits + 63) / 64 * 64;
    lvl->elems     = (bits + BCVEC_BITS_PER_ENTRY - 1) / BCVEC_BITS_PER_ENTRY;
    // This also zeros the bits.
    if (bits < ((1UL<<20) * 8)) {
    alloc_small:
        size_t alloc_size = lvl->elems * sizeof(struct bbhash_colvec);
        int rc = posix_memalign((void**)&lvl->col_vec, 
                                sizeof(struct bbhash_colvec),
                                alloc_size);
        if (rc != 0) die("posix_memalign");
        memset(lvl->col_vec, 0, alloc_size);
        lvl->mmap_size = 0;
    } else {
        size_t size = lvl->elems*sizeof(struct bbhash_colvec);
        // Round up to next page size
        size += (1UL<<21) - 1;
        size &= ~((1UL<<21) - 1);
        lvl->mmap_size = size;
        lvl->col_vec = mmap(0, size, PROT_READ|PROT_WRITE,
                         MAP_ANON|MAP_PRIVATE|MAP_POPULATE|MAP_HUGETLB,
                         -1, 0);
        if (lvl->col_vec == MAP_FAILED) {
            perror("mmap(MAP_HUGETLB) failed, not enough huge pages?");
            goto alloc_small;
        }
    }
     
    return lvl;
}

static void bbhash_level_free(struct bbhash_level *lvl) {
    if (lvl->mmap_size > 0 ) {
        int rc = munmap(lvl->col_vec, lvl->mmap_size);
        if (rc != 0) die("munmap");
    } else {
        free(lvl->col_vec);
    }
    free(lvl);
}

#define destruct_bit_position(bit) \
    uint64_t bvec, entry, offset, mask;          \
    bvec   = bit / BCVEC_BITS_PER_ENTRY;          \
    offset = bit % BCVEC_BITS_PER_ENTRY;           \
    entry  = offset / 64;                           \
    offset = offset % 64;                             \
    mask = (1UL << offset)

static inline _Bool bbhash_colvec_get(struct bbhash_level * lvl, uint64_t bit) {
    destruct_bit_position(bit);
    
    if (lvl->col_vec[bvec].v[entry] & mask)
        return true;
    return false;
}

static inline void bbhash_colvec_prefetch(struct bbhash_level * lvl, uint64_t bit, bool rw) {
    destruct_bit_position(bit);
    
    __builtin_prefetch(&(lvl->col_vec[bvec].v[entry]), rw, 1);
}


// Returns true on collisison
static inline _Bool bbhash_colvec_inc(struct bbhash_level *lvl, uint64_t bit) {
    destruct_bit_position(bit);

    uint64_t col = lvl->col_vec[bvec].v[entry] & mask;

    if (0) {
        if (!col) {
            lvl->col_vec[bvec].v[entry] |= mask;
            return false;
        } else {
            lvl->col_vec[bvec].c[entry] |= mask;
            return true;
        }
    } else {
        // This variant is actually slower, although being branchless
        lvl->col_vec[bvec].v[entry] |= mask;
        lvl->col_vec[bvec].c[entry] |= col;
        return col != 0;
    }

    
}

static inline _Bool bbhash_get(struct bbhash_level *lvl, uint64_t bit) {
    if (lvl->v[bit / 64] & (1UL << (bit % 64)))
        return true;
    return false;
}

static inline void bbhash_prefetch(struct bbhash_level *lvl,
                                          uint64_t bit, bool rw) {
    __builtin_prefetch(&(lvl->v[bit / 64]), rw, 1);
}

static struct bbhash _BBHashCompute(double gamma, uintptr_t _keys,
                            unsigned long key_count) {
    struct bbhash ret = {0};
    
    uint64_t *keys = (uint64_t *)_keys;
    uint64_t keys_remaining = key_count;

    uint64_t *scratch = malloc(key_count     * sizeof(uint64_t));
    uint64_t  redo_idx, recheck_idx;

    #define STRIDE 128
    uint64_t *hashes = malloc(sizeof(uint64_t)*STRIDE);

    uint64_t depth_sum = 0;

    uint64_t rank = 1;

    for (int level = 0; ; level ++) {
        // Create a new vector to our result bbhash
        struct bbhash_level *lvl = bbhash_level_alloc(keys_remaining * gamma);
        ret.levels = realloc(ret.levels, sizeof(struct bbhash_level *) * ++ret.level_count);
        ret.levels[level] = lvl;

        // Precalculate the level hash to speed up the hash calculation
        uint64_t lvl_hash = levelHash(level);

        // Reset scratch pointers. The redoc vector grows from the
        // beginning of scratch, while the recheck vector grows from
        // the end.
        redo_idx = 0; // redo index grows up
        recheck_idx = key_count - 1; // recheck index grows down

        for (unsigned long i = 0; i < keys_remaining;) {
            int count = MIN(STRIDE, keys_remaining - i);
            for (int ii = 0; ii < count; ii++) {
                uint64_t h = keyHash(lvl_hash, keys[i+ii]);
                h = fastrange64(h, lvl->bits);
                hashes[ii] = h;
                bbhash_colvec_prefetch(lvl, h, true);
            }
            for (int ii = 0; ii < count; ii++) {
                // Set the bits in the collision vector
                _Bool collision = bbhash_colvec_inc(lvl, hashes[ii]);

                // Partition keys array into keys that we surely have to
                // redo and those that might be ok.
                if (collision) {
                    scratch[redo_idx++] = keys[i+ii];
                } else {
                    scratch[recheck_idx--] = keys[i+ii];
                }
            }
            // Move loop forward
            i += count;
        }

        // This loops combines both collision bit vectors into a
        // single bit-vector. Furthermore, we compact the previously
        // interleaved bit-vector into a single dense bit-vector.
        uint64_t bits = 0, collisions = 0;
        for (unsigned i = 0; i < lvl->elems; i++) {
            struct bbhash_colvec *e = &lvl->col_vec[i];
            for (unsigned ii = 0; ii < ARRAY_SIZE(e->v); ii++){
                uint64_t bloom = e->v[ii] & ~(e->c[ii]);
                lvl->v[i * ARRAY_SIZE(e->v) + ii] = bloom;
                bits += __builtin_popcountll(bloom);
                //collisions += __builtin_popcountll(lvl->col_vec[i].c[ii]);
            }
        }
        
        lvl->rank = rank;
        rank += bits;

        // This loop iterates over all keys that are scheduled for
        // rechecking as they did not provoke a collision during their
        // bit-set operation. For collisiding bits, this loop will
        // move all keys to the redo vector that set the col_vec->v
        // bit.

        // IMPORTANT: In this loop the bit-vector is already
        // compressed, hence we have to use bbhash_{prefetch, get}.
        for (unsigned long i = recheck_idx+1; i < key_count;) {
            int count = MIN(STRIDE, key_count - i);
            for (int ii = 0; ii < count; ii++) {
                uint64_t h = keyHash(lvl_hash, scratch[i+ii]);
                h = fastrange64(h, lvl->bits);
                hashes[ii] = h;
                bbhash_prefetch(lvl, h, false);
            }
            for (int ii = 0; ii < count; ii++) {
                if (!bbhash_get(lvl, hashes[ii])) {
                    scratch[redo_idx++] = scratch[i+ii];
                }
            }
            i += count;
        }

        printf("lvl %d: %lu / %lu (%d redo)\n", level, bits, lvl->bits, redo_idx);

        depth_sum += bits * (level + 1);

        if (redo_idx == 0) {
            break;
        }

        keys_remaining = redo_idx;
        keys = scratch;
    }

    printf("average depth: %.2f\n", depth_sum/(double)key_count);

    free(scratch);
    free(hashes);

    return ret;
}

void BBHashFree(struct bbhash hh) {
    printf("Free a bbhash with %d levels\n", hh.level_count);
    for (unsigned i = 0; i < hh.level_count; i++) {
        bbhash_level_free(hh.levels[i]);
    }
}

uint64_t * BBHashGetLevel(struct bbhash hh, int i, uint64_t *bits, uint64_t *rank) {
    *bits = hh.levels[i]->bits;
    *rank = hh.levels[i]->rank;
    return hh.levels[i]->v;
}

////////////////////////////////////////////////////////////////
// Measurement code
////////////////////////////////////////////////////////////////
#include <time.h>

inline double time_delta(struct timespec *const ts) {
    struct timespec N;
    if (clock_gettime(CLOCK_REALTIME, &N)) {
        exit(1);
    }
    double ret = ((double)N.tv_sec - ts->tv_sec) +
                 (double)(N.tv_nsec - ts->tv_nsec) / 1e9;
    *ts = N;
    return ret;
}

// When reading from the perf descriptor, the kernel returns an
// event record in the following format (if PERF_FORMAT_GROUP |
// PERF_FORMAT_ID are enabled).
// Example (with id0=100, id1=200): {.nr = 2, .values = {{41433, 200}, {42342314, 100}}}
typedef uint64_t perf_event_id; // For readability only
struct read_format {
    uint64_t nr;
    struct {
        uint64_t value;
        perf_event_id id; // PERF_FORMAT_ID
    } values[];
};

// Structure to hold a perf group
struct perf_handle {
    int group_fd;   // First perf_event fd that we create
    int nevents;    // Number of registered events
    size_t rf_size; // How large is the read_format buffer (derived from nevents)
    struct read_format *rf; // heap-allocated buffer for the read event
};

// Syscall wrapper for perf_event_open(2), as glibc does not have one
int sys_perf_event_open(struct perf_event_attr *attr,
                    pid_t pid, int cpu, int group_fd,
                    unsigned long flags) {
    return syscall(__NR_perf_event_open, attr, pid, cpu, group_fd, flags);
}

// Add a perf event of given (type, config) with default config. If p
// is not yet initialized (p->group_fd <=0), the perf_event becomes
// the group leader. The function returns an id that can be used in
// combination with perf_event_get.
perf_event_id perf_event_add(struct perf_handle *p, int type, int config) {
    // SOL_IF
    struct perf_event_attr attr;

    memset(&attr, 0, sizeof(struct perf_event_attr));
    attr.type = type;
    attr.size = sizeof(struct perf_event_attr);
    attr.config = config;
    attr.disabled = 1;
    attr.exclude_kernel = 1;
    attr.exclude_hv = 1;
    attr.read_format = PERF_FORMAT_GROUP | PERF_FORMAT_ID;
    int fd = sys_perf_event_open(&attr, 0, -1,
                             p->group_fd > 0 ? p->group_fd : -1,
                             0);
    if (fd < 0) die("perf_event_open");
    if (p->group_fd <= 0)
        p->group_fd = fd;

    p->nevents ++;

    perf_event_id id;
    if (ioctl(fd, PERF_EVENT_IOC_ID, &id) < 0)
        die("perf/IOC_ID");
    return id;
    // SOL_ELSE
    // FIXME: Create event with perf_event_open
    // FIXME: Get perf_event_id with PERF_EVENT_IOC_ID
    return -1;
    // SOL_END
}

// Resets and starts the perf measurement
void perf_event_start(struct perf_handle *p) {
    // SOL_IF
    // Reset and enable the event group
    ioctl(p->group_fd, PERF_EVENT_IOC_RESET,  PERF_IOC_FLAG_GROUP);
    ioctl(p->group_fd, PERF_EVENT_IOC_ENABLE, PERF_IOC_FLAG_GROUP);
    // SOL_ELSE
    // FIXME: PERF_EVENT_IOC_{RESET, ENABLE}
    // SOL_END
}

// Stops the perf measurement and reads out the event
void perf_event_stop(struct perf_handle *p) {
    // SOL_IF
    // Stop the tracing for the whole event group
    ioctl(p->group_fd, PERF_EVENT_IOC_DISABLE, PERF_IOC_FLAG_GROUP);

    // Allocate a read_format buffer if not done yet.
    if (p->rf == NULL) {
        p->rf_size = sizeof(uint64_t) + 2 * p->nevents * sizeof(uint64_t);
        p->rf = malloc(p->rf_size);
    }

    // get the event from the kernel. Our buffer should be sized exactly righ
    if (read(p->group_fd, p->rf, p->rf_size) < 0)
        die("read");
    // SOL_ELSE
    // FIXME: PERF_EVENT_IOC_DISABLE
    // FIXME: Read event from the group_fd into an allocated buffer
    // SOL_END
}


// After the measurement, this helper extracts the event counter for
// the given perf_event_id (which was returned by perf_event_add)
uint64_t perf_event_get(struct perf_handle *p, perf_event_id id) {
    // SOL_IF
    for (unsigned i = 0; i < p->rf->nr; i++) {
        if (p->rf->values[i].id == id) {
            return p->rf->values[i].value;
        }
    }
    // SOL_END
    return -1;
}



struct bbhash BBHashCompute(double gamma, uintptr_t keys,
                            unsigned long key_count) {

    struct timespec start;
    // Create and initialize a new perf handle
    struct perf_handle p;
    memset(&p, 0, sizeof(p));

    // Create three new perf events that we want to monitor for our
    // matrix multiplication algorithms
    perf_event_id id_instrs =
        perf_event_add(&p, PERF_TYPE_HARDWARE, PERF_COUNT_HW_INSTRUCTIONS);
    perf_event_id id_cycles =
        perf_event_add(&p, PERF_TYPE_HARDWARE, PERF_COUNT_HW_CPU_CYCLES);
    perf_event_id id_cache_miss = 
        perf_event_add(&p, PERF_TYPE_HARDWARE, PERF_COUNT_HW_CACHE_MISSES);
    perf_event_id id_cache_refs =
        perf_event_add(&p, PERF_TYPE_HARDWARE, PERF_COUNT_HW_CACHE_REFERENCES);
    // perf_event_id id_branch_miss =
    //     perf_event_add(&p, PERF_TYPE_HARDWARE, PERF_COUNT_HW_BRANCH_MISSES);
    // perf_event_id id_branchs =
    //     perf_event_add(&p, PERF_TYPE_HARDWARE, PERF_COUNT_HW_BRANCH_INSTRUCTIONS);

    perf_event_start(&p);
    time_delta(&start);

    struct bbhash ret = _BBHashCompute(gamma, keys, key_count);

    double delta = time_delta(&start);
    perf_event_stop(&p);

    double instrs = perf_event_get(&p, id_instrs);
    double cycles = perf_event_get(&p, id_cycles);
    double misses = perf_event_get(&p, id_cache_miss);
    double refs = perf_event_get(&p, id_cache_refs);

    //double branch_miss = perf_event_get(&p, id_branch_miss)
    //    / (double) perf_event_get(&p, id_branchs);



    printf("C version: %.2f s, cache-miss ratio: %f%%, ins/cycle: %.2f, mem-BW: %.2f MiB/s\n", delta,
           misses / refs *100,
           instrs/cycles,
           misses * 64 / (1 << 20)
        );
    printf("Per Key: %.2f cyc., %.3f D$-miss, %.3f D$-acc\n",
           instrs / key_count,
           misses / key_count,
           refs   / key_count);
    
    return ret;
}
