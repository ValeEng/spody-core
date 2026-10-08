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
#include <math.h>
#include "spody_random.h"

/* Philox4x64-10 (Salmon, Moraes, Dror, Shaw, SC'11, 2011): the two
 * round multipliers and the two Weyl increments added to the key
 * between rounds (golden ratio and sqrt(3) - 1 as 64-bit fractions).
 * Private to this generator, so they live here, not in spody_const.h. */
static const uint64_t philox_m0 = 0xD2E7470EE14C6C93ULL;
static const uint64_t philox_m1 = 0xCA5A826395121157ULL;
static const uint64_t philox_w0 = 0x9E3779B97F4A7C15ULL;
static const uint64_t philox_w1 = 0xBB67AE8584CAA73BULL;
enum { PHILOX_ROUNDS = 10 };

/* 64-bit FNV-1a (Fowler, Noll, Vo): offset basis and prime, for the
 * substream id of a named quantity. */
static const uint64_t fnv1a64_offset = 0xCBF29CE484222325ULL;
static const uint64_t fnv1a64_prime  = 0x00000100000001B3ULL;

/* Step of the open uniform grid (k + 1/2) * 2^-52, see
 * spody_random_u64_to_open01. */
static const double two_pow_m52 = 2.220446049250313080847263336181640625e-16;

#if defined(_MSC_VER) && !defined(SPODY_RANDOM_PORTABLE_MUL)
#include <intrin.h>
#endif
#if !defined(SPODY_RANDOM_PORTABLE_MUL) && defined(__SIZEOF_INT128__)
__extension__ typedef unsigned __int128 uint128;   /* GCC/Clang extension */
#endif

/* ---- 64 x 64 -> 128-bit product ---------------------------------------
 *
 * Each Philox round multiplies two 64-bit words and keeps BOTH halves
 * of the 128-bit product: the high half depends on every bit of the
 * inputs and does the mixing. Standard C has no 128-bit integer, so the
 * product comes from the compiler where it has one (GCC/Clang
 * __int128, MSVC _umul128 on x64 and __umulh on ARM64), otherwise from
 * four 32 x 32 -> 64 partial products recombined by hand. Defining
 * SPODY_RANDOM_PORTABLE_MUL forces the by-hand path, so a test can
 * check that all paths give the same bits. */
static void mulhilo64(uint64_t a, uint64_t b, uint64_t *hi, uint64_t *lo)
{
#if !defined(SPODY_RANDOM_PORTABLE_MUL) && defined(__SIZEOF_INT128__)
    const uint128 p = (uint128)a * b;
    *hi = (uint64_t)(p >> 64);
    *lo = (uint64_t)p;
#elif !defined(SPODY_RANDOM_PORTABLE_MUL) && defined(_MSC_VER) && defined(_M_X64)
    *lo = _umul128(a, b, hi);
#elif !defined(SPODY_RANDOM_PORTABLE_MUL) && defined(_MSC_VER) && defined(_M_ARM64)
    *lo = a * b;
    *hi = __umulh(a, b);
#else
    const uint64_t a_lo = (uint32_t)a, a_hi = a >> 32;
    const uint64_t b_lo = (uint32_t)b, b_hi = b >> 32;
    const uint64_t p0 = a_lo * b_lo, p1 = a_lo * b_hi;
    const uint64_t p2 = a_hi * b_lo, p3 = a_hi * b_hi;
    /* middle column: carry out of the low word plus the two cross
     * terms' low halves; at most 3 * (2^32 - 1), no overflow */
    const uint64_t mid = (p0 >> 32) + (uint32_t)p1 + (uint32_t)p2;
    *lo = (mid << 32) | (uint32_t)p0;
    *hi = p3 + (p1 >> 32) + (p2 >> 32) + (mid >> 32);
#endif
}

void spody_philox4x64(const uint64_t ctr[4], const uint64_t key[2],
                      uint64_t out[4])
{
    uint64_t c0 = ctr[0], c1 = ctr[1], c2 = ctr[2], c3 = ctr[3];
    uint64_t k0 = key[0], k1 = key[1];
    for (int r = 0; r < PHILOX_ROUNDS; ++r) {
        if (r > 0) {             /* key schedule: Weyl increments */
            k0 += philox_w0;
            k1 += philox_w1;
        }
        uint64_t hi0, lo0, hi1, lo1;
        mulhilo64(philox_m0, c0, &hi0, &lo0);
        mulhilo64(philox_m1, c2, &hi1, &lo1);
        c0 = hi1 ^ c1 ^ k0;
        c1 = lo1;
        c2 = hi0 ^ c3 ^ k1;
        c3 = lo0;
    }
    out[0] = c0; out[1] = c1; out[2] = c2; out[3] = c3;
}

void spody_random_stream_init_domain(SpodyRandomStream *s, uint64_t seed,
                                     uint64_t id, uint64_t sub, uint64_t dom)
{
    s->key[0] = seed; s->key[1] = sub;
    s->ctr[0] = 0;    s->ctr[1] = id; s->ctr[2] = dom; s->ctr[3] = 0;
    s->block[0] = s->block[1] = s->block[2] = s->block[3] = 0;
    s->next = 4;
}

void spody_random_stream_init(SpodyRandomStream *s, uint64_t seed,
                              uint64_t id, uint64_t sub)
{
    spody_random_stream_init_domain(s, seed, id, sub, SPODY_RANDOM_DOMAIN_DRAW);
}

int spody_gauss_markov_nodes(SpodyRandomStream *s, double sigma, const double *scale,
                             double tau, const double *t, size_t n, double *x)
{
    if (n == 0 || !(sigma >= 0.0) || !(tau > 0.0) || !isfinite(sigma)
        || !isfinite(tau) || !isfinite(t[0]))
        return -1;
    for (size_t j = 0; j < n; ++j) {
        if (j > 0 && (!(t[j] > t[j - 1]) || !isfinite(t[j]))) return -1;
        if (scale && (!(scale[j] >= 0.0) || !isfinite(scale[j]))) return -1;
    }
    /* Gillespie (1996); 1 - phi^2 = -expm1(-2 dt / tau) keeps
     * full precision when dt << tau. */
    x[0] = (scale ? sigma * scale[0] : sigma) * spody_random_next_normal(s);
    for (size_t j = 1; j < n; ++j) {
        const double dt  = t[j] - t[j - 1];
        const double phi = exp(-dt / tau);
        const double sj  = scale ? sigma * scale[j] : sigma;
        x[j] = phi * x[j - 1]
             + sj * sqrt(-expm1(-2.0 * dt / tau)) * spody_random_next_normal(s);
    }
    return 0;
}

uint64_t spody_random_next_u64(SpodyRandomStream *s)
{
    if (s->next == 4) {
        spody_philox4x64(s->ctr, s->key, s->block);
        ++s->ctr[0];
        s->next = 0;
    }
    return s->block[s->next++];
}

uint64_t spody_random_substream_id(const char *name)
{
    uint64_t h = fnv1a64_offset;
    for (const unsigned char *p = (const unsigned char *)name; *p; ++p) {
        h ^= *p;
        h *= fnv1a64_prime;
    }
    return h;
}

double spody_random_u64_to_open01(uint64_t x)
{
    return ((double)(x >> 12) + 0.5) * two_pow_m52;
}

/* ---- AS241 (PPND16) ---------------------------------------------------
 *
 * Three rational approximations of degree 7/7: the central one in
 * r = 0.180625 - q^2 for |p - 1/2| <= 0.425, two tail ones in
 * r = sqrt(-ln(min(p, 1-p))) split at r = 5. Coefficients and splits
 * as published by Wichura (1988); they are one table, kept here like
 * the RK45 tableau in spody_integrators.c. Checked against
 * scipy.special.ndtri to 1.0e-15 relative over p in [1e-300, 1-1e-6]
 * plus the two grid extremes. */
static const double as241_a[8] = {
    3.3871328727963666080,     1.3314166789178437745e+2,
    1.9715909503065514427e+3,  1.3731693765509461125e+4,
    4.5921953931549871457e+4,  6.7265770927008700853e+4,
    3.3430575583588128105e+4,  2.5090809287301226727e+3 };
static const double as241_b[8] = {
    1.0,                       4.2313330701600911252e+1,
    6.8718700749205790830e+2,  5.3941960214247511077e+3,
    2.1213794301586595867e+4,  3.9307895800092710610e+4,
    2.8729085735721942674e+4,  5.2264952788528545610e+3 };
static const double as241_c[8] = {
    1.42343711074968357734,    4.63033784615654529590,
    5.76949722146069140550,    3.64784832476320460504,
    1.27045825245236838258,    2.41780725177450611770e-1,
    2.27238449892691845833e-2, 7.74545014278341407640e-4 };
static const double as241_d[8] = {
    1.0,                       2.05319162663775882187,
    1.67638483018380384940,    6.89767334985100004550e-1,
    1.48103976427480074590e-1, 1.51986665636164571966e-2,
    5.47593808499534494600e-4, 1.05075007164441684324e-9 };
static const double as241_e[8] = {
    6.65790464350110377720,    5.46378491116411436990,
    1.78482653991729133580,    2.96560571828504891230e-1,
    2.65321895265761230930e-2, 1.24266094738807843860e-3,
    2.71155556874348757815e-5, 2.01033439929228813265e-7 };
static const double as241_f[8] = {
    1.0,                       5.99832206555887937690e-1,
    1.36929880922735805310e-1, 1.48753612908506148525e-2,
    7.86869131145613259100e-4, 1.84631831751005468180e-5,
    1.42151175831644588870e-7, 2.04426310338993978564e-15 };
static const double as241_split_central = 0.425;
static const double as241_const_central = 0.180625;   /* 0.425^2 */
static const double as241_split_tail    = 5.0;
static const double as241_const_tail    = 1.6;

/* Horner evaluation of the degree-7 polynomial c[0] + c[1] x + ... */
static double poly7(const double c[8], double x)
{
    double s = c[7];
    for (int i = 6; i >= 0; --i) s = s * x + c[i];
    return s;
}

double spody_normal_quantile(double p)
{
    if (!(p > 0.0 && p < 1.0)) return NAN;      /* also catches NaN */
    const double q = p - 0.5;
    if (fabs(q) <= as241_split_central) {
        const double r = as241_const_central - q * q;
        return q * poly7(as241_a, r) / poly7(as241_b, r);
    }
    double r = sqrt(-log(q < 0.0 ? p : 1.0 - p));
    double z;
    if (r <= as241_split_tail) {
        r -= as241_const_tail;
        z = poly7(as241_c, r) / poly7(as241_d, r);
    } else {
        r -= as241_split_tail;
        z = poly7(as241_e, r) / poly7(as241_f, r);
    }
    return q < 0.0 ? -z : z;
}

double spody_random_next_normal(SpodyRandomStream *s)
{
    return spody_normal_quantile(
        spody_random_u64_to_open01(spody_random_next_u64(s)));
}
