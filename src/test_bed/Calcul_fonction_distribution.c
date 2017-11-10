

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "./../h_ions/h_ion_chargestate.h"
#include "./../h_ions/h_ion_ziegler_stoppow.h"
#include "./../h_ions/h_ion_srivastava_stoppow.h"
#include "./../generic/constant.h"
#include "./../simple_models/hion_spectrum.h"
#include "./../simple_models/hion_gastransport.h"

double Gryz_Eexch_ie_csection(double Vp,double Mp,double q, double dE, double Ubind);
int main()
{

double E0=68.0e6;
double Mp=139;
double Zp=57;
double R=7.74e-6;
//double E0=99.01e6;
//double Mp=95;
//double Zp=42;
//double R=10.02e-6;

double L=3e-6;
double n=1.65;
double Zt=18.0;
double Mt=40.0;


 printf("Chargement des donnees P aret proton\n");

    setup_Ziegler_elec_stoppow();

double * tab=simple_transport_hion_gas( E0, R, n,L,1.0, Mp, Zp, Zt, 2000);


FILE *fich=fopen("spectre_e.txt","wb");

 double dV=sqrt((E0)*OD_e*2.0/(Mp*OD_uma))/2000;
 double V1=0;
double V2=0;
 double q=0;
 double Se1=0;
 double Se2=0;

 double integ=0;
 double nbe=0;
 double Emean=0;
for(int j=0;j<4000;j++)
{
    integ=0;
for(int i=1;i<2000-1;i++)
{
V1=dV*(double)i;
q= q_Schiwietz_Grande(V1, Zp, Zt);
Se1=6.0*Gryz_Eexch_ie_csection(V1,Mp,q,(double)j, 15.7);
Se1+=2.0*Gryz_Eexch_ie_csection(V1,Mp,q,(double)j, 29.3);

V2=dV*(double)(i+1);
q= q_Schiwietz_Grande(V2, Zp, Zt);
Se2=6.0*Gryz_Eexch_ie_csection(V2,Mp,q, (double)j, 15.7);
Se2+=2.0*Gryz_Eexch_ie_csection(V2,Mp,q, (double)j, 29.3);

integ+=(V1*tab[i*2+1]*Se1*1e-4+V2*tab[(i+1)*2+1]*Se2*1e-4)*0.5*(V2-V1);
}
printf("%d %e\n",j,integ);
fprintf(fich,"%e %e\n",(double)j,integ*2.67e25);
nbe+=integ*2.67e25;
Emean+=integ*2.67e25*(double)j;

}

printf("nbe %e Emean %e",nbe,Emean/nbe);
fclose(fich);
  free_Ziegler_elec_stoppow();

return 0;
}
