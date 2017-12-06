

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "./../h_ions/h_ion_chargestate.h"
#include "./../h_ions/h_ion_ziegler_stoppow.h"
#include "./../h_ions/h_ion_srivastava_stoppow.h"
#include "./../generic/constant.h"
#include "hion_spectrum.h"

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
	double s = 1.0/pow(4*M_PI*epsilon0,2.0)*1.0/(V*V0) * (4.0*M_PI*pow(1.6e-19,4.0)*N/(95.0*1.66e-27*9.1e-31))*q*q*f*(3.0*pow(ki,-1.0/3.0)+(1.0/ki));
	//double s =M_PI*2.5027e+025*pow(6.62e-34/(2.*M_PI),2.0)/(95.0*1.66e-27*9.1e-31)*pow(Zp,0.9)*f*(3.0*pow(Zp,-0.15)+pow(Zp,-0.45))*V/V0;
 	return s;
}


// Ziegler_elec_stoppow(double V,double Mp,double Zp,double Zt)
//double E,double L, double E0, double RE0,double n,double rrate

double * simple_transport_hion_gas( double E0, double RE0, double n,double L,double rr, double Mp, double Zp, double Zt, int nE )
{

    double *tab=malloc(sizeof(double)*nE*2);

    double dV=sqrt((E0)*OD_e*2.0/(Mp*OD_uma))/nE;
    double R=0;

    double V1,V2;
    double S1,S2;
    double f1,f2;
// calcul du range
FILE *fich=fopen("comparaison.txt","w");
    for( int i= 1;i<nE;i++)
    {
        V1=dV*(double)i;
        S1=1.e6*OD_e/601.0*100.0*Ziegler_elec_stoppow(V1,Mp,Zp,Zt)*1./(Mp*OD_uma)+OD_e*2.50e4*100.0*Ziegler_nuc_stoppow(V1,Mp,Zp, 40.0,Zt)*1./(Mp*OD_uma);

        V2=dV*(double)(i+1);
        S2=1.e6*OD_e/601.0*100.0*Ziegler_elec_stoppow(V2,Mp,Zp,Zt)*1./(Mp*OD_uma)+OD_e*2.50e4*100.0*Ziegler_nuc_stoppow(V2,Mp,Zp, 40.0,Zt)*1./(Mp*OD_uma);
//printf("V1=%e S1=%e V2=%e S2=%e Par2=%e \n",V1,S1,V2,S2, S2/puissanceArret(V2, Zt,Zp));
        fprintf(fich,"%e %e %e %e %e\n", V2*V2*0.5*(Mp*OD_uma)/(1.6e-19), V2,1.e6*OD_e/601.0*100.0*Ziegler_elec_stoppow(V2,Mp,Zp,Zt)*1./(Mp*OD_uma),OD_e*2.50e4*100.0*Ziegler_nuc_stoppow(V2,Mp,Zp, 40.0,Zt)*1./(Mp*OD_uma),puissanceArret(V2, Zt,Zp) );

        if(isnan(S1) || isnan(S2))
        {
             S1=OD_e*2.50e4*100.0*Ziegler_nuc_stoppow(V1,Mp,Zp, 40.0,Zt)*1./(Mp*OD_uma);
             S2=OD_e*2.50e4*100.0*Ziegler_nuc_stoppow(V2,Mp,Zp, 40.0,Zt)*1./(Mp*OD_uma);
            R+=(V1/S1+V2/S2)*0.5*(V2-V1);
      //  printf("S1=%e\n",S1);
        }
        else
        {
        R+=(V1/S1+V2/S2)*0.5*(V2-V1);
        }

    }
 printf("R=%e",R);
   fclose(fich);
/*
	for(int i = 0; i<nv; i++)
	{
	 	//printf("nv %d\n",i);
		H[i] = 0;

		for(int j = i; j<nv-1; j++)
		{

		 	S1 = puissanceArret( (double)j*dv);

			if(j==0)
			double h1=0;
			else{


			f1 = spectre(i*dv);
					double h1 = f1/S1;


			S2 = puissanceArret((double) (j+1)*dv);
			f2 = spectre((double) (j+1)*dv);
						double h2 =f2/S2;

								H[i]+= (h1 + h2)*0.5*dv;
				}

		}

	}
	*/

	fich=fopen("fonction_distri.txt","w");

    for( int i= 1;i<nE;i++)
    {
            tab[i*2]=dV*(double)(i);
            tab[i*2+1]=0;
            printf("%d\n",i);
        for(int j = i; j<nE-1; j++)
		{

                V1=dV*(double)j;
                f1=hion_analytic_spec_dep_o(0.5*(Mp*OD_uma)*V1*V1/(OD_e),L, E0, RE0,n,1)*Mp*OD_uma*V1/OD_e;

                S1=1.e6*OD_e/601.0*100.0*Ziegler_elec_stoppow(V1,Mp,Zp,Zt)*1./(Mp*OD_uma)+OD_e*2.50e4*100.0*Ziegler_nuc_stoppow(V1,Mp,Zp, 40.0,Zt)*1./(Mp*OD_uma);


                V2=dV*(double)(j+1);
                f2=hion_analytic_spec_dep_o(0.5*(Mp*OD_uma)*V2*V2/(OD_e),L, E0, RE0,n,1)*Mp*OD_uma*V2/OD_e;
                S2=1.e6*OD_e/601.0*100.0*Ziegler_elec_stoppow(V2,Mp,Zp,Zt)*1./(Mp*OD_uma)+OD_e*2.50e4*100.0*Ziegler_nuc_stoppow(V2,Mp,Zp, 40.0,Zt)*1./(Mp*OD_uma);
                if(isnan(S1) || isnan(S2))
                {
                S1=OD_e*2.50e4*100.0*Ziegler_nuc_stoppow(V1,Mp,Zp, 40.0,Zt)*1./(Mp*OD_uma);
                S2=OD_e*2.50e4*100.0*Ziegler_nuc_stoppow(V2,Mp,Zp, 40.0,Zt)*1./(Mp*OD_uma);
                }
                tab[i*2+1]+=(f1/S1 + f2/S2)*0.5*(V2-V1);
		}
        tab[i*2+1]/=R;

        fprintf(fich,"%e %e  %e %e\r\n",tab[i*2],tab[i*2+1],0.5*(Mp*OD_uma)*tab[i*2]*tab[i*2]/(OD_e),hion_analytic_spec_dep_o(0.5*(Mp*OD_uma)*tab[i*2]*tab[i*2]/(OD_e),L, E0, RE0,n,1));

    }
    double sum=0;
    for( int i= 1;i<nE-1;i++)
    {
        V1=dV*(double)i;
        f1=hion_analytic_spec_dep_o(0.5*(Mp*OD_uma)*V1*V1/(OD_e),L, E0, RE0,n,1)*Mp*OD_uma*V1/OD_e;

        V2=dV*(double)(i+1);
        f2=hion_analytic_spec_dep_o(0.5*(Mp*OD_uma)*V2*V2/(OD_e),L, E0, RE0,n,1)*Mp*OD_uma*V2/OD_e;

//        sum+=(f1+f2)*0.5*(0.5*(Mp*OD_uma)*V2*V2/(OD_e)-0.5*(Mp*OD_uma)*V1*V1/(OD_e));
            sum+=(f1+f2)*0.5*(V2-V1);

    }

    printf("integrale %e",sum);
    fclose(fich);



    return tab;
}

