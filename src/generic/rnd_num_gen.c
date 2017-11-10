

#include <stdint.h>
/* routine permettant de generer un nombre aleatoire selon la methode xorshift 128 de
Marsaglia.

Voir: "Xorshift RNGs". Marsaglia. Journal of Statistical Software 8 (14).

*/
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
