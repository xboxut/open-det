



#define OD_SURF_NONE		0
#define OD_SURF_XPLANE		1
#define OD_SURF_YPLANE		2
#define OD_SURF_ZPLANE		3	
#define OD_SURF_PLANE		4

#define OD_SURF_SPHERE		5

#define OD_SURF_ZCYLIND		6
#define OD_SURF_XCYLIND		7
#define OD_SURF_YCYLIND		8
#define OD_SURF_CYLIND		9

/****************************************************************
//DONNEE
OD_rng_num
-----------------------------------------------------------------
utilité:

Structure permettant la gestion de la géométrie.
Pour l'instant, on s'intéresse aux probleme 0D ou 1D, mais par 
la suite, on aimerait étendre les possibilités du code
en simulant des géométries 3 dimensions. 
La description des géométries se fait en assemblant des surfaces 
ensembles. Par soucis de simplicité, l'ensemble des surface générées
par le code sera stockée dans la liste chainée base geom.
-----------------------------------------------------------------
Provenance:
CEA Cadarache, DEN/DER/SPEx/LPE
G. de Izarra 
------------------------------------------------------------------



*/



typedef struct{

int id;//identifier utile notamment pour le binding du code avec une interface lua ou autre suivi couant a travers la suface etc...

char type; //type de surface définie: plan sphere, cylindre ou autre
double coordinates[9];//coordonnées  de la surface 

int mat_id; // identifier d'un materiaux pour des intéractions surfaciques. Utile par exemple pour un modèle simple d'optique avec absorption reflexion.
int voltage;// Utile pour calculer le champ électrique 

OD_geom_surf *next;



}OD_geom_surf;




extern OD_geom_surf * base_geom_surf;



/* routine pour vider la liste chainée des surfaces*/
int free_base_geom_surf();

