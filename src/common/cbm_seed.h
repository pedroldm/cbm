/*
 * Deterministic seed derivation shared by ENS (C), ILS and CBMLKH (C++).
 *
 * Every random decision of a run is a pure function of the run's seed and a
 * logical index (iteration number, trajectory index, sub-problem hash), never
 * of the clock, the PID or the order in which threads finish. The Python
 * runner mirrors cbm_derive_seed() so the derived seeds can be recomputed.
 */
#ifndef CBM_SEED_H
#define CBM_SEED_H

#include <stdint.h>

static inline uint64_t cbm_splitmix64(uint64_t x)
{
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

/* Seed of sub-stream `index` of `seed`, in [1, 2^31 - 1]: valid for LKH's SEED
 * (unsigned) and Linkern's -s (int, where 0 would mean "seed from the clock"). */
static inline uint32_t cbm_derive_seed(uint64_t seed, uint64_t index)
{
    uint64_t h = cbm_splitmix64(cbm_splitmix64(seed) ^ index);
    return (uint32_t)(h % 2147483647ULL) + 1u;
}

/* Small portable PRNG (splitmix64 stream): unlike rand(), its sequence does not
 * depend on the C library, so a seed replays identically on another machine. */
typedef struct { uint64_t state; } cbm_rng;

static inline uint64_t cbm_rng_next(cbm_rng *rng)
{
    rng->state += 0x9E3779B97F4A7C15ULL;
    uint64_t z = rng->state;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

#endif /* CBM_SEED_H */
