/*
 * Copyright 2026 ValeEng
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#ifndef SPODY_RANDOM_H
#define SPODY_RANDOM_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================
 * Reproducible random numbers for Monte Carlo sampling
 *
 * Generator: Philox4x64-10, a counter-based generator (Salmon,
 * Moraes, Dror, Shaw, "Parallel random numbers: as easy as 1, 2, 3",
 * SC'11, 2011). A block of four 64-bit words is a pure function of a
 * 256-bit counter and a 128-bit key, so any word of any stream is
 * computed directly, with no state carried from the words before it
 * and no overlap between streams by construction. The same algorithm
 * is numpy.random.Philox, std::philox4x64 (C++26) and Random123's
 * philox4x64_10; this implementation is written from the paper and
 * checked against the Random123 known-answer vectors.
 *
 * Stream convention (one stream per Monte Carlo case and per
 * dispersed quantity):
 *
 *     word k of stream (seed, id, sub, dom) = Philox(ctr = {k / 4, id, dom, 0},
 *                                                    key = {seed, sub})[k % 4]
 *
 * so case `id` draws the same numbers whatever the number of cases,
 * the thread count or the order the cases run in, and each dispersed
 * quantity reads its own substream `sub`: adding, removing or
 * reordering one quantity never changes the draws of the others
 * (common random numbers between two configurations). Substream 0 is
 * the initial state; a parameter takes spody_random_substream_id of
 * its name. Two quantities of one run must NEVER share a substream
 * (they would draw identical numbers): callers check it, and any new
 * kind of dispersed quantity needs an id no other can take.
 *
 * The third counter word is the stream DOMAIN: 0 for the draws made
 * once per case (initial state, parameters), 1 for the process noise
 * drawn along the trajectory. Two streams of different domains never
 * share a counter, whatever their substreams, so the separation holds
 * by construction (no hash can collide across it), and adding process
 * noise leaves every domain-0 draw unchanged.
 *
 * Normal deviates: the inverse of the standard normal CDF applied to
 * one uniform (AS241, see spody_normal_quantile). One uniform gives
 * one normal, monotonically, which is what stratified designs (Latin
 * hypercube, scrambled Sobol) need and Box-Muller does not give.
 *
 * Thread safety: no global state. A SpodyRandomStream belongs to its
 * caller (one per case); spody_philox4x64 and spody_normal_quantile
 * are pure.
 * ============================================================ */

/* One Philox4x64-10 block: out = Philox(ctr, key). */
void spody_philox4x64(const uint64_t ctr[4], const uint64_t key[2],
                      uint64_t out[4]);

/* Sequential reader over one stream: refills a block every 4 words. */
typedef struct {
    uint64_t key[2];
    uint64_t ctr[4];     /* ctr[0] = index of the NEXT block to compute */
    uint64_t block[4];
    int      next;       /* next unused word of block; 4 = refill first */
} SpodyRandomStream;

/* Stream domains: the third counter word (see the convention above). */
enum {
    SPODY_RANDOM_DOMAIN_DRAW          = 0,
    SPODY_RANDOM_DOMAIN_PROCESS_NOISE = 1
};

/* Position `s` at word 0 of stream (seed, id, sub) of domain
 * SPODY_RANDOM_DOMAIN_DRAW. */
void spody_random_stream_init(SpodyRandomStream *s, uint64_t seed,
                              uint64_t id, uint64_t sub);

/* Position `s` at word 0 of stream (seed, id, sub) of domain `dom`. */
void spody_random_stream_init_domain(SpodyRandomStream *s, uint64_t seed,
                                     uint64_t id, uint64_t sub, uint64_t dom);

/* Substream id of a named quantity: the 64-bit FNV-1a hash
 * (Fowler-Noll-Vo) of the NUL-terminated `name`, e.g. the batch
 * target path "spacecraft.drag.Cd". Stable across versions as long as
 * the name is. Not collision-free by construction: the caller rejects
 * a run where two quantities (or a quantity and the reserved 0) get
 * the same id. */
uint64_t spody_random_substream_id(const char *name);

/* Next 64-bit word of the stream. */
uint64_t spody_random_next_u64(SpodyRandomStream *s);

/* Map a 64-bit word to the open interval (0, 1): the top 52 bits k
 * give u = (k + 1/2) * 2^-52, the odd multiples of 2^-53. Every value
 * is exact in a double, 0 and 1 never occur, and u and 1 - u are both
 * on the grid, so the extreme normal deviate is -/+8.2095 sigma
 * exactly symmetric. (53 bits with the half step would put the top
 * value at 1 - 2^-54, which rounds to 1.0.) */
double spody_random_u64_to_open01(uint64_t x);

/* Standard normal quantile z = Phi^-1(p) for p in (0, 1): Wichura,
 * "Algorithm AS 241: The percentage points of the normal
 * distribution", Applied Statistics 37(3), 477-484, 1988 (PPND16,
 * about 16 significant digits). Returns NaN for p <= 0, p >= 1 or a
 * NaN argument: the stream never produces those, so a NaN downstream
 * is a caller bug, not a tail value. */
double spody_normal_quantile(double p);

/* Next standard normal deviate of the stream: one word, one deviate. */
double spody_random_next_normal(SpodyRandomStream *s);

/* First-order Gauss-Markov (Ornstein-Uhlenbeck) process
 *
 *     dx/dt = -x / tau + w(t),   stationary variance sigma^2,
 *     E[x(t) x(t + d)] = sigma^2 exp(-|d| / tau)
 *
 * (Uhlenbeck & Ornstein, Phys. Rev. 36, 823, 1930) sampled EXACTLY at
 * the n strictly increasing times t[0..n-1] (Gillespie, "Exact
 * numerical simulation of the Ornstein-Uhlenbeck process and its
 * integral", Phys. Rev. E 54(2), 2084, 1996):
 *
 *     x[0]   = sigma z[0]                      (stationary start)
 *     x[j+1] = phi x[j] + sigma sqrt(1 - phi^2) z[j+1],
 *     phi    = exp(-(t[j+1] - t[j]) / tau)
 *
 * with z[j] the next n normal deviates of `s`, in order. The node
 * values have exactly the statistics of the continuous process for
 * any spacing.
 *
 * `scale` (NULL = constant sigma) gives a time-varying sigma_j =
 * sigma scale[j] at node j, in both terms above: the variance then
 * relaxes towards sigma_j^2 with time constant tau / 2, the
 * activity-scaled stochastic error of Wright (AGI, "Real-time
 * estimation of local atmospheric density"). With scale NULL the
 * result is bit-identical to the constant-sigma recursion.
 *
 * Returns 0, or -1 (x untouched) for n = 0, sigma < 0, a scale < 0,
 * tau <= 0, non-finite inputs or non-increasing times. */
int spody_gauss_markov_nodes(SpodyRandomStream *s, double sigma, const double *scale,
                             double tau, const double *t, size_t n, double *x);

#ifdef __cplusplus
}
#endif

#endif /* SPODY_RANDOM_H */
