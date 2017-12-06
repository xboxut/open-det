
/****************************************************************
Définition des particules pour une simulation monte carlo des traces
d'ionisation.




*/
#define PARTICLE_ATOM                   0
#define PARTICLE_MOLECULE               1
#define PARTICLE_ION                    2
#define PARTICLE_ION_MOLECULE           3
#define PARTICLE_ELECTRON               4


typedef struct{

// données spatiales
OD_vect3d position;//definition de la position de la particule.
OD_vect3d direction;// vecteur  directeur de la particule.

// donnée internes
double kin_energy; // énergie cinétique en eV
double M; // Masse en uma

double creation_time;// temps de création -> utile par exemple pour suivre la durée de vie des états excités

char type;// type de particule considéré

// Définition du type d'ion ou d'atome
double Z;
double *Zmol; // definition des molecules -> a faire
double charge_state;

int level_ref[4];

}pic_particle;
