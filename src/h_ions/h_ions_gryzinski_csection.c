
#include <math.h>


/****************************************************************
Gryz_Eexch_ie_csection()
-----------------------------------------------------------------
utilité:
Routine pour calculer la section efficace d'echange d énergie
entre un ion et un électron du cortège électronique d'une cible
de Gryzinski.

Les hypothèses effectuées sont :
- des collisions binaires
- une fonction de distribution des vitesses des électrons estimée
 via une "analyse asymptotique".
- un traitement classique des collisions.
-----------------------------------------------------------------
Provenance:

classical theory of atomic collision 1: theory of inelastic
collision.
de Michal Gryzinski
Phys. Rev. 138, 336-358
------------------------------------------------------------------
La fonction prend en paramètre:
Vp      la vitesse du projectile (en m/s)
q       la charge électrique du projectile (en charge elementaire)
dE      l'énergie échangée (en ev)
Ubind  l'énergie de liaison des électrons de la couche considérée
		(en ev)
*****************************************************************/

double Gryz_Eexch_ie_csection(double Vp,double q, double dE, double Ubind)
{


//const double v0=2.179e6;// en m/s (V0=e²/(hbar*4*pi*epsilon 0)
const double sigma0=q*q*6.56e-14; // en cm² eV² (sigma0=pi e^4Zq^2)

double ve=sqrt(Ubind*2.0/9.10938356e-31); // Definition de la vitesse électronique en m/s
double epsilon=0.5*9.10938356e-31*(ve*ve)/1.602e-19; //energie cinetique de l'électron en eV.



double fv=(ve/Vp)*(ve/Vp)*pow(Vp*Vp/(Vp*Vp+ve*ve),3.0/2.0);
//printf("fv %e\n",fv);
double DEmax=epsilon*4.0*pow(Vp/ve,2.0)*(1.0+ve/Vp);
//printf("DEmax %e\n",DEmax);
double G=fv*(Vp*Vp/(Vp*Vp+ve*ve)*dE/epsilon+4.0/3.0*log(2.7+Vp/ve))*(1.0-pow(dE/DEmax,1+Vp*Vp/(ve*ve)));
//printf("G %e %e\n",G,(1.0-pow(dE/DEmax,1+Vp*Vp/(ve*ve))));

return sigma0/pow(dE,3.0)*G;

}



