




/****************************************************************
//DONNEE
OD_rng_num
-----------------------------------------------------------------
utilit�:

Structure permettant de g�rer de multiple g�n�rateurs de nombre
al�atoire. Cela permet de s'affranchir de la GSL et de proposer
des g�n�rateur de nombre al�atoires rapides.
-----------------------------------------------------------------
Provenance:
CEA Cadarache, DEN/DER/SPEx/LPE
G. de Izarra
------------------------------------------------------------------



*/

typedef struct{

const char *name;// nom du g�n�rateur

uint8_t warm_up;// s il y a eu un warm up ou non.
uint8_t state_sze; // taille des �tats du g�n�rateur
void *seeds;// les graines utilis�es, en cas de besoin.




void *states;// etat courant du g�n�rateur

// pointeur sur les fonctions du g�n�rateur
 void (*set)(void * st,void *se, uint32_t seed);
 uint32_t (*get32)(void *st);
 uint64_t (*get64)(void *st);
 double (*get_unif_double)(void *st);
 double (*get_fast_unif_double)(void *st);
 float (*get_unif_float)(void *st);
 void(*jump)(void *st);

}OD_rng;


/****************************************************************
//DONNEE
Constantes de type OD_rng
-----------------------------------------------------------------
utilit�:

Ces structures contiennent tout le n�cessaire pour d�finir des
g�n�rateur de nombre al�atoire � l'aide de la routine init_rng.
-----------------------------------------------------------------
Provenance:
CEA Cadarache, DEN/DER/SPEx/LPE
G. de Izarra
------------------------------------------------------------------



*/

extern const OD_rng * xorwow; // p�riode de 2^192-2^32
//extern const OD_rng * xorshift64s;
extern const OD_rng * xorshift128;
extern const OD_rng * xorshift128pv1;
extern const OD_rng * xorshift128pv2;
extern const OD_rng * xoroshiro128;
//extern const OD_rng * xorshift1024p;



 OD_rng *init_rng(const OD_rng * rng);
 void seed_rng(OD_rng *rng, uint32_t seed);
 uint64_t rng_uniform64(OD_rng *rng);
 uint32_t rng_uniform(OD_rng * rng);
 double drng_uniformf(OD_rng *rng);
 double drng_uniform(OD_rng *rng);
 void rng_jump(OD_rng *rng);

