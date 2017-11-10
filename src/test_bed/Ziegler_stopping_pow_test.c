

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "./../h_ions/h_ion_chargestate.h"
#include "./../h_ions/h_ion_ziegler_stoppow.h"
#include "./../h_ions/h_ion_srivastava_stoppow.h"
#include "./../generic/constant.h"
int main()
{
    double V=0;
	printf("*********************************************************\n");

	printf("*********Essai des puissances d'arret de Ziegler*********\n\n");

    printf("Chargement des donnees P aret proton\n");

    setup_Ziegler_elec_stoppow();
    system("echo \%cd\%");


	printf("Creation de puissance d'arret: niobium dans argon\n");
    FILE *fich=fopen("NB_in_AR_Ziegler.txt","w");
	for(int i=10;i<100000;i++)
    {
        V=sqrt(2.0*(double)i*1e3*OD_e/(92.0*OD_uma));
    fprintf(fich,"%d %e %e\n",i, Ziegler_elec_stoppow(V,92.0,41.0,18.0), Ziegler_nuc_stoppow(V,92,41, 40,18));

    }
    fclose(fich);
	printf("                               niobium dans xenon\n");

	fich=fopen("NB_in_XE_Ziegler.txt","w");
	for(int i=10;i<100000;i++)
    {
        V=sqrt(2.0*(double)i*1e3*OD_e/(92.0*OD_uma));
    fprintf(fich,"%d %e %e\n",i, Ziegler_elec_stoppow(V,92.0,41.0,54.0), Ziegler_nuc_stoppow(V,92,41, 40,54));

    }
    fclose(fich);


    free_Ziegler_elec_stoppow();


    printf("*********************************************************\n");
	printf("**Essai des puissances d'arret de SRIVASTAVA ET MUKHERJI**\n\n");



    printf("Creation de puissance d'arret: niobium dans argon\n");
    fich=fopen("NB_in_AR_Srivastava.txt","w");
	for(int i=10;i<100000;i+=10)
    {
        V=sqrt(2.0*(double)i*1e3*OD_e/(92.0*OD_uma));
    fprintf(fich,"%d %e\n",i, 1000.0*Srivastava_elec_stoppow( V,92.0,41.0,q_Schiwietz_Grande(V, 41.0, 18.0),40.0,18.0));

    }
     fclose(fich);

     fich=fopen("NB_in_XE_Srivastava.txt","w");
	for(int i=10;i<100000;i+=10)
    {
        V=sqrt(2.0*(double)i*1e3*OD_e/(92.0*OD_uma));
    fprintf(fich,"%d %e\n",i, 1000.0*Srivastava_elec_stoppow( V,92.0,41.0,q_Schiwietz_Grande(V, 41.0, 54.0),131.0,54.0));

    }
     fclose(fich);

return 0;
}
