


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "./../h_ions/h_ion_srivastava_stoppow.h"
#include "./../generic/constant.h"




double Srivastava_elec_stoppow(double V,double Mp,double Zp,double eff_charge,double Mt,double Zt)
{
    double S=0;

//    B=v/c -> c v de la lumiere
//    x=2zV0/V z charge ionique
   // double b=V/OD_c;
    double X=2.0*eff_charge*OD_v0/V;
  //  printf("X %e Eff charge %e V0  %e V %e\n",X,eff_charge,OD_v0,V );
    double I=0;// a calculer
    double fz=0;

    // calcul de f(z)
    if(Zt<=45)
    fz=0.28*pow(Zt,2.0/3.0);
    else
    fz=pow(Zt,1.0/3.0);


    I=(Zt-2.0)/Zt*log(13.6*pow((Zt-2.0)/(2.717*fz),2.0))+2.0*log(13.6*Zt*Zt);
    I=exp(I);

    //printf("Valeur de I %e\n",I);
    if(X>1.0)
    {

        if(V>=0.5*Zt*OD_v0*X)
        {

            S=63.65*eff_charge*eff_charge*Zt/(Mt*(V*V*1e-12))*log10(11.39*(V*V*1e-12)/(I*X));
        }
        else if(0.5*Zt*OD_v0*X>V && V>=0.5*Zt*OD_v0*pow(X,1./3.))
        {
    //        printf("ttttttt\n");

            S=13.79*eff_charge*eff_charge/(Mt*(V*V*1e-12))*(3.0*(Zt-2.0)+3.0*(Zt-2.0)*log(2.0*fz*V/((Zt-2.0)*OD_v0))+6.0*log(2.0*V/(Zt*OD_v0))+2.0*fz*V/(OD_v0*X)-Zt*log(X));
        }
        else if(0.5*Zt*OD_v0*pow(X,1./3.)>V)
        {
              printf("ttttttt\n");
            S=12.68*fz*eff_charge*eff_charge/(Mt*V*1e-6)*(3.0*pow(X,-1.0/3.0)+1/X);

        }


    }
    else
    {

        if(V>=0.5*Zt*OD_v0)
        {
            S=63.65*eff_charge*eff_charge*Zt/(Mt*(V*V*1e-12))*log10(11.39*(V*V*1e-12)/I);

        }
        else
        {
            S=60.6*fz*eff_charge*eff_charge/(Mt*V*1e-6);

        }

    }


 return S;
}
