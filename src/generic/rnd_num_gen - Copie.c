

#include <stdint.h>





/* routine permettant de generer un nombre aleatoire selon la methode xorshift 128 de
Marsaglia.

Voir: "Xorshift RNGs". Marsaglia. Journal of Statistical Software 8 (14).

*/
//seed
uint32_t xorshft_x=20, xorshft_y=546216547, xorshft_z=687954679, xorshft_w=1478529;

inline uint32_t xorshift128(void) {
    uint32_t t = xorshft_x ^ (xorshft_x << 11);
    xorshft_x = xorshft_y; xorshft_y = xorshft_z; xorshft_z = xorshft_w;
    return xorshft_w = xorshft_w ^ (xorshft_w >> 19) ^ t ^ (t >> 8);
}

float xorshift128_uniform(void) {
    uint32_t t = xorshft_x ^ (xorshft_x << 11);
    xorshft_x = xorshft_y; xorshft_y = xorshft_z; xorshft_z = xorshft_w;
    xorshft_w = xorshft_w ^ (xorshft_w >> 19) ^ t ^ (t >> 8);
    return (float)(0.5+ (signed)xorshft_w *0.2328306e-9);
}



/*

// The state array must be initialized to not be all zero
uint32_t xorshift128(uint32_t state[static 4])
{
	// Algorithm "xor128" from p. 5 of Marsaglia, "Xorshift RNGs"
	uint32_t s, t = state[3];
	t ^= t << 11;
	t ^= t >> 8;
	state[3] = state[2]; state[2] = state[1]; state[1] = s = state[0];
	t ^= s;
	t ^= s >> 19;
	state[0] = t;
	return t;
}

*/
