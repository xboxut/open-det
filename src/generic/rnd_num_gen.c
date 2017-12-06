
#include <stdio.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>


#include "rnd_num_gen.h"
/***********************************************************************************
*DONNEES GENERALES
***********************************************************************************/

inline OD_rng *init_rng(const OD_rng * rng)
{
	OD_rng *ret=malloc(sizeof(OD_rng));

	ret->name=rng->name;
	ret->state_sze=rng->state_sze;


	ret->seeds=malloc(ret->state_sze);
	ret->states=malloc(ret->state_sze);
	ret->set=rng->set;
	ret->get32=rng->get32;
	ret->get64=rng->get64;
	ret->get_unif_double=rng->get_unif_double;
	ret->get_fast_unif_double=rng->get_fast_unif_double;
	ret->get_unif_float=rng->get_unif_float;
	ret->jump=rng->jump;
	return ret;
}

inline void seed_rng(OD_rng *rng, uint32_t seed)
{

    rng->set(rng->seeds,rng->states,seed);

}


inline uint64_t rng_uniform64(OD_rng *rng)
{
return rng->get64(rng->states);
}

inline uint32_t rng_uniform(OD_rng * rng)
{
return rng->get32(rng->states);
}
inline double drng_uniformf(OD_rng *rng)
{
	return rng->get_fast_unif_double(rng->states);
}
inline double drng_uniform(OD_rng *rng)
{
	return rng->get_unif_double(rng->states);
}

inline void rng_jump(OD_rng *rng)
{
   rng->jump(rng->states);
}

/***********************************************************************************
*DONNEES RELATIVES A XOSHIFT128
***********************************************************************************/




/****************************************************************
uint32_t xorshift128_(void *state)
-----------------------------------------------------------------
utilité:
Routine permettant de génerer des séquences de nombres pseudo-
aléatoires. Ce générateur est relativement rapide *4 par rapport
à Mersenne twister GSL.

-----------------------------------------------------------------
Provenance:

Xorshift RNGs". Marsaglia. Journal of Statistical Software 8 (14)
------------------------------------------------------------------La routine:
* La routine à une période de 2^128-1
* Elle réussi les tests die hard,
* Elle loupe MaxtrixRank et Linearcomp de BIGcrush

*****************************************************************/
static inline uint32_t xorshift128_(void *state)
{
	/* Algorithm "xor128" from p. 5 of Marsaglia, "Xorshift RNGs" */
	uint32_t *st=(uint32_t *)state;
	uint32_t s, t = st[3];
	t ^= t << 11;
	t ^= t >> 8;
	st[3] = st[2]; st[2] = st[1]; st[1] = s = st[0];
	t ^= s;
	t ^= s >> 19;
	st[0] = t;
	return t;
}


float xorshift128_uniform_(void *state)
{
	return (float)(0.5+ (signed)xorshift128_(state) *0.2328306e-9);
}

const OD_rng *xorshift125=NULL;

/***********************************************************************************
*DONNEES RELATIVES A XORWOW
***********************************************************************************/


/* The state array must be initialized to not be all zero in the first four words */
static uint32_t xorwow_(void *st)
{
	uint32_t *state=(uint32_t *)st;
	/* Algorithm "xorwow" from p. 5 of Marsaglia, "Xorshift RNGs" */
	uint32_t s, t = state[3];
	t ^= t >> 2;
	t ^= t << 1;
	state[3] = state[2]; state[2] = state[1]; state[1] = s = state[0];
	t ^= s;
	t ^= s << 4;
	state[0] = t;
	return t + (state[4] += 362437);
}

const OD_rng *xorwow=NULL;
/***********************************************************************************
*DONNEES RELATIVES A XORSHIFT128+
***********************************************************************************/

static uint64_t xorshift128plusv1(void *st) {

	uint64_t *s=(uint64_t *)st;

	uint64_t x = s[0];
	uint64_t const y = s[1];
	s[0] = y;
	x ^= x << 23; // a
	s[1] = x ^ y ^ (x >> 17) ^ (y >> 26); // b, c
	return s[1] + y;
}

/*  Written in 2014-2016 by Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>. */

/* This generator has been replaced by xoroshiro128plus, which is
   significantly faster and has better statistical properties.

   It might be nonetheless useful for languages in which low-level rotate
   instructions are not available. Due to the relatively short period it
   is acceptable only for applications with a mild amount of parallelism;
   otherwise, use a xorshift1024* generator.

   Note that the lowest bit of this generator is an LFSR of degree 128;
   thus, it will fail linearity tests. The next bit can be described by an
   LFSR of degree 8256, but in the long run it will fail linearity tests,
   too. The other bits needs a much higher degree to be represented as
   LFSRs.

   We suggest to use a sign test to extract a random Boolean value, and
   right shifts to extract subsets of bits.

   The state must be seeded so that it is not everywhere zero. If you have
   a 64-bit seed, we suggest to seed a splitmix64 generator and use its
   output to fill s.

   A previous version of this generator was adding the two halves of the
   newly computed state. This version adds the two halves of the *current*
   state (as xoroshiro128plus does), which improves speed due to better
   internal parallelization from the CPU. The resulting streams are off by
   one step. */


static uint64_t xorshift128plusv2(void *st)
{
	uint64_t *s=(uint64_t *)st;
	uint64_t s1 = s[0];
	const uint64_t s0 = s[1];
	const uint64_t result = s0 + s1;
	s[0] = s0;
	s1 ^= s1 << 23; // a
	s[1] = s1 ^ s0 ^ (s1 >> 18) ^ (s0 >> 5); // b, c
	return result;
}

/* This is the jump function for the generator. It is equivalent
   to 2^64 calls to next(); it can be used to generate 2^64
   non-overlapping subsequences for parallel computations. */

static void jumpxorshift128plusv2_(void *st) {
	static const uint64_t JUMP[] = { 0x8a5cd789635d2dff, 0x121fd2155c472f96 };
	uint64_t *s=(uint64_t *)st;
	uint64_t s0 = 0;
	uint64_t s1 = 0;
	for(int i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
		for(int b = 0; b < 64; b++) {
			if (JUMP[i] & UINT64_C(1) << b) {
				s0 ^= s[0];
				s1 ^= s[1];
			}
			xorshift128plusv2(st) ;
		}

	s[0] = s0;
	s[1] = s1;
}

const OD_rng *xorshift128pv1=NULL;
const OD_rng *xorshift128pv2=NULL;

/***********************************************************************************
*DONNEES RELATIVES A XOROSHIRO128
***********************************************************************************/



/****************************************************************
uint64_t xoroshiro128(void* st)
-----------------------------------------------------------------
utilité:
Routine permettant de génerer des séquences de nombres pseudo-
aléatoires.

La fonction retourne un uint64_t. Pour tirer un nombre aléatoire
double depuis ce générateur, se référer aux routines
xoroshio128d précise mais contenant une multiplication ou
 à xoroshiro128fastd, qui permet de convertir un uint64_t via un
 shift de bit vers la mantisse.

-----------------------------------------------------------------
Provenance:

Sebastiano Vigna et David Blackman;
http://xoroshiro.di.unimi.it/

Written in 2016 by David Blackman and Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>.
 This is the successor to xorshift128+. It is the fastest full-period
   generator passing BigCrush without systematic failures, but due to the
   relatively short period it is acceptable only for applications with a
   mild amount of parallelism; otherwise, use a xorshift1024* generator.

   Beside passing BigCrush, this generator passes the PractRand test suite
   up to (and included) 16TB, with the exception of binary rank tests, as
   the lowest bit of this generator is an LFSR of degree 128. The next bit
   can be described by an LFSR of degree 8256, but in the long run it will
   fail linearity tests, too. The other bits needs a much higher degree to
   be represented as LFSRs.

   We suggest to use a sign test to extract a random Boolean value, and
   right shifts to extract subsets of bits.

   Note that the generator uses a simulated rotate operation, which most C
   compilers will turn into a single instruction. In Java, you can use
   Long.rotateLeft(). In languages that do not make low-level rotation
   instructions accessible xorshift128+ could be faster.

   The state must be seeded so that it is not everywhere zero. If you have
   a 64-bit seed, we suggest to seed a splitmix64 generator and use its
   output to fill s.
------------------------------------------------------------------

La routine prend en paramètre un pointeur ver l'état du générateur de nombre
aléatoire. Pour éviter tout problème, au démarrage, l'état du générateur
ne doit pas etre nul.
*****************************************************************/


static inline uint64_t rotl(const uint64_t x, int k) {
	return (x << k) | (x >> (64 - k));
}

static inline uint64_t xoroshiro128_(void* st) {

	uint64_t *s=(uint64_t *)st;
	const uint64_t s0 = s[0];
	uint64_t s1 = s[1];
	const uint64_t result = s0 + s1;

	s1 ^= s0;
	s[0] = rotl(s0, 55) ^ s1 ^ (s1 << 14); // a, b
	s[1] = rotl(s1, 36); // c

	return result;
}

static  double xoroshiro128d(void* st)
{
	return (xoroshiro128_(st) >> 11) * (1. / (UINT64_C(1) << 53));
}

static  double xoroshiro128fastd(void *st)
{
		union { uint64_t i; double d; } u = { .i = UINT64_C(0x3FF) << 52 | xoroshiro128_(st) >> 12 };
       return u.d - 1.0;
}

static  float xoroshiro128f(void *st)
{
		union { uint64_t i; double d; } u = { .i = UINT32_C(0x3F) << 23 | xoroshiro128_(st) >> (64-23) };
       return u.d - 1.0;
}

static  uint32_t xoroshiro32(void *st) {
    union { uint64_t value; struct { uint32_t high, low; }; } converter;
    converter.value=xoroshiro128_(st);
    return converter.high;
}
/* This is the jump function for the generator. It is equivalent
   to 2^64 calls to next(); it can be used to generate 2^64
   non-overlapping subsequences for parallel computations. */

void jumpxoroshiro128(void* st) {
	static const uint64_t JUMP[] = { 0xbeac0467eba5facb, 0xd86b048b86aa9922 };
    uint64_t *s=(uint64_t *)st;
	uint64_t s0 = 0;
	uint64_t s1 = 0;
	for(int i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
		for(int b = 0; b < 64; b++) {
			if (JUMP[i] & UINT64_C(1) << b) {
				s0 ^= s[0];
				s1 ^= s[1];
			}
			xoroshiro128_(st);
		}

	s[0] = s0;
	s[1] = s1;
}

static void xoroshiro128_set(void * st,void *se, uint32_t seed)
{

// Le set de xoroshiro utilise le générateur de nombre aléatoire mixsplit64 de Vigna.

    uint64_t *s=(uint64_t *)st;
    uint64_t *see=(uint64_t *)se;
    uint64_t x=(uint64_t)seed;
    uint64_t z = (x += 0x9e3779b97f4a7c15);
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
	z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
	s[0]=z ^ (z >> 31);


    x += 0x9e3779b97f4a7c15;
        x += 0x9e3779b97f4a7c15;
            x += 0x9e3779b97f4a7c15;

    z = (x += 0x9e3779b97f4a7c15);
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
	z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
	s[1]=z ^ (z >> 31);

    see[0]=s[0];
    see[1]=s[1];

}

/*
typedef struct{

const char *name;

uint8_t warm_up;

void *seeds;


uint64_t limit;

void *states;


 void (*set)(uint32_t seed);
 uint32_t (*get)();
 uint64_t (*get)();
 double (*get_unif_double)();
 double (*get_quick_unif_double)();
 float (*get_unif_float)();
 void(*jump)();

}OD_rng;
*/
static const OD_rng ran_xoshi128 = {
  "xoroshiro128",                   /* name */
  0,                           /* WARM_UP */
  2*sizeof(uint64_t),							/* 2 etat de 64 bit*/
  NULL,                            /* Sauvegarde de graines */
  NULL,								/* etat courant*/
  &xoroshiro128_set,						/*set*/
  &xoroshiro32,						/* nombre aleat uint32*/
  &xoroshiro128_,
  &xoroshiro128d,						/* nombre aleat double*/
  &xoroshiro128fastd,						/* nombre aleat double rapide*/
	&xoroshiro128f,						/* nombre aleat flottant->> a tester*/
  &jumpxoroshiro128						/* fonction jump*/
};
const OD_rng *xoroshiro128 = &ran_xoshi128;


/*
 static inline double to_double(uint64_t x) {
       const union { uint64_t i; double d; } u = { .i = UINT64_C(0x3FF) << 52 | x >> 12 };
       return u.d - 1.0;
    }
*/
/*
uint32_t lower_32_bits(uint64_t value) {
    union { uint64_t value; struct { uint32_t high, low; }; } converter;
    converter.value = value;
    return converter.low;
}

*/
