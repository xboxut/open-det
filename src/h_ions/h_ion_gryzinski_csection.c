
#include <math.h>

#include "./../generic/constant.h"

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

La routine retourne la section efficace en cm²/ev

Note: les sections efficaces de Grynzinski ne prennent en compte qu'
un seul électron de la couche électronique des atomes cibles.
Il est nécessaire de sommer la contribution de l'ensemble des
electron de l'atome pour avoir un résultat correct.
Généralement, seul les électrons de la couche externe comptent
car leur énergie de liaison est un ou deux ordre de grandeur
plus faible que celle des couches pus profondes.
Or, la section efficace varie comme 1/dE^3 avec dE l'énergie
transférée à l'électron-> pour certain e-, du fait de l'e de liaison
importante, la section efficace est négligeable.

L'amplitude des sections efficaces de gryzinski peut être loin
de la valeur exacte. Selon la méthode de calcul du ralentissement
des ions énergétiques, il peut être judicieux de normaliser
l'intégrale de la section efficace de transfert d'énergie
par le pouvoir d'arret électronique afin d'être sur de conserver
l'énergie.
*****************************************************************/

double Gryz_Eexch_ie_csection(double Vp,double Mp,double q, double dE, double Ubind)
{

//const double v0=2.179e6;// en m/s (V0=e²/(hbar*4*pi*epsilon 0)
const double sigma0=q*q*6.56e-14; // en cm² eV² (sigma0=pi e^4Zq^2)

double ve=sqrt(Ubind*OD_e*2.0/OD_me); // Definition de la vitesse électronique en m/s

double epsilon=Ubind; //energie cinetique de l'électron en eV.
//double epsilon=0.5*OD_uma*Mp*Vp*Vp/OD_e;


double fv=(ve/Vp)*(ve/Vp)*pow(Vp*Vp/(Vp*Vp+ve*ve),3.0/2.0);
//printf("fv %e\n",fv);
double DEmax=epsilon*4.0*pow(Vp/ve,2.0)*(1.0+ve/Vp);
//printf("Demax  %e vp/ve %e\n",DEmax,Vp/ve);
//system("pause");

//printf("DEmax %e\n",DEmax);
double G=fv*(Vp*Vp/(Vp*Vp+ve*ve)*(dE+Ubind)/epsilon+4.0/3.0*log(2.7+Vp/ve))*(1.0-pow((dE+Ubind)/DEmax,1+Vp*Vp/(ve*ve)));
//printf("G %e %e\n",G,(1.0-pow(dE/DEmax,1+Vp*Vp/(ve*ve))));
if(	(dE+Ubind)/DEmax>1)
{
	 	G = 0.0;
}

return sigma0/pow(dE+Ubind,3.0)*G;

}




/****************************************************************
Gryz_Stop_pow_i_csection()
-----------------------------------------------------------------
utilité:
Routine pour calculer la section efficace de pouvoir d'arret
entre un ion et un atome de Gryzinski. Seul les processus d'excitation
et d'ionisation (collision inelastiques) sont couverts par ce modèle.
D'autre processus comme le ralentissement coulombien ne sont pas
pris en compte.

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
double ne,double Ep,double Eexc,double Vp,double Ubind)
La fonction prend en paramètre:
ne      Le nombre d'electron susceptible susceptible de participer
        à la collision
Eexv    L'énergie d'excitation du niveau ou l'énergie d'ionisation.
Vp      la vitesse du projectile (en m/s)
Ubind  l'énergie de liaison des électrons de la couche considérée
		(en ev)

La routine retourne la section efficace en cm².ev

Note: les sections efficaces de Grynzinski ne prennent en compte qu'
un seul électron de la couche électronique des atomes cibles.
Il est nécessaire de prendre en compte la contribution de l'ensemble des
electrons de l'atome pour avoir un résultat correct.

L'amplitude des sections efficaces de gryzinski peut être loin
de la valeur exacte. Selon la méthode de calcul du ralentissement
des ions énergétiques, il peut être judicieux de normaliser
l'intégrale de la section efficace de transfert d'énergie
par le pouvoir d'arret électronique afin d'être sur de conserver
l'énergie.
*****************************************************************/

double Gryz_Stop_pow_i_csection(double ne,double q,double Ep,double Eexc,double Vp,double Ubind)
{

//const double v0=2.179e6;// en m/s (V0=e²/(hbar*4*pi*epsilon 0)
const double sigma0=q*q*6.56e-14; // en cm² eV² (sigma0=pi e^4Zq^2)

double ve=sqrt(Ubind*2.0/9.10938356e-31); // Definition de la vitesse électronique en m/s
double epsilon=0.5*9.10938356e-31*(ve*ve)/1.602e-19; //energie cinetique de l'électron en eV.



double fv=(ve/Vp)*(ve/Vp)*pow(Vp*Vp/(Vp*Vp+ve*ve),3.0/2.0);
//printf("fv %e\n",fv);
double DEmax=epsilon*4.0*pow(Vp/ve,2.0)*(1.0+ve/Vp);

//y=fv*(vq^2/(vq^2+ve^2)*Uexc/Ekn*log(Demax/Uexc)+4/3*(1-Uexc/Demax)*log(2.7+vq/ve))*(1-(Uexc/Demax)^(1+vq^2/ve^2));

double Gs=fv*(Vp*Vp/(Vp*Vp+ve*ve)*Eexc/Ep*log(DEmax/Eexc)+4.0/3.0*(1.0-Eexc/DEmax)*log(2.7+Vp/ve)*(1.0-pow(Eexc/DEmax,1.0+Vp*Vp/(ve*ve))));

return ne*sigma0/Eexc*Gs;
}
