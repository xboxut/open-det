
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "./../../h_ions/h_ion_chargestate.h"
#include "./../../h_ions/h_ion_ziegler_stoppow.h"
#include "./../../h_ions/h_ion_srivastava_stoppow.h"
#include "./../../generic/constant.h"

double puissanceArret(double V, double Zt, double Zp)  //v en m/s
{
	//constantes
    const double V0=2.19e6; // vitesse de bohr en m/s
	const double epsilon0=8.854187e-12;
	double N = 2.50e+025; //densite du milieu en 1/m3
	//double A  = 39.;     //nombre d'atome du gaz

	//la distribution electronique dans les atomes cibles
	double f = pow(Zt,1.0/3.0);

	//etat de charge BETZ
	//double q = 0.5*V*pow(Zp, 0.45)/V0;

	//etat de charge Schiwietz pour les  projectiles ayant un z compris entre 1 et 92 pour des cibles GAZEUSES avec un Z entre 1 et 54.
    double q=q_Schiwietz_Grande(V, Zp, Zt);


	//condition
	double ki = 2.0*V0*q/V;

	//la puissance d'arret
	double s = 1.0/pow(4*M_PI*epsilon0,2.0)*1.0/(V*V0) * (4.0*M_PI*pow(1.6e-19,4.0)*N/(92.0*1.66e-27*9.1e-31))*q*q*f*(3.0*pow(ki,-1.0/3.0)+(1.0/ki));
	//double s =M_PI*2.5027e+025*pow(6.62e-34/(2.*M_PI),2.0)/(95.0*1.66e-27*9.1e-31)*pow(Zp,0.9)*f*(3.0*pow(Zp,-0.15)+pow(Zp,-0.45))*V/V0;
 	return s;
}


int main()
{

    FILE* fich=fopen("comparaison_StoppingPOW.txt","w");

   double dE=100;
   double V=sqrt((dE*1000*1.6e-19*2)/(95*OD_uma));
   double Skush;
   double Ssrvi;
   double E=100;
    for(int i=1;i<1000;i++)
    {
        E+=dE;
        V=sqrt((E*1000*1.6e-19*2)/(92*OD_uma));
        Skush=puissanceArret(V, 18.0, 41.0);
       // Ssrvi=1e3*OD_e*1.6637*Ziegler_elec_stoppow(V,92.0,41.0,18.0)
        Ssrvi=1000*Srivastava_elec_stoppow(V,92.0,41.0,q_Schiwietz_Grande(V, 41.0, 18.0),40.0,18.0);

        printf("%e %e %e\n",i*dE,Skush*(92.0*1.66e-27)/(100*Ssrvi),Ssrvi);///(95.0*1.66e-27)
    }

    fclose(fich);
    return 0;
}
