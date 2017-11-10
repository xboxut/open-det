
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "./../generic/constant.h"
#include "./../generic/file_loader.h"
#include "./../generic/interpolation.h"

/****************************************************************
Ziegler_nuc_stoppow()
-----------------------------------------------------------------
utilité:

Routine permettant de calculer la puissance d'arret nucléaire
d'un projectile P sur un milieu d'éléments Zt. La routine retourne
la puissance d'arret en ev/(10^15 atome/cm²)
-----------------------------------------------------------------
Provenance:

(nouvelle version)
SRIM The Stopping and Range of Ions in matter.
J.F. Ziegker M.D. Ziegler J.P. Biersack

(ancienne version)
The electronic and nuclear stopping of energetic ions
J.F. Ziegler
Apllied physics Letters 31, 544(1977)
------------------------------------------------------------------
La fonction prend en paramètre:

V la vitesse du projectile (en m.s)
Mp la masse du projectile (en uma)
Zp la charge du projectile (en charge élémentaire)
Mt la masse de la cible (en uma)
Zt la charge de la cible (en charge élémentaire)

La fonction retourne la puissance d'arret nucléaire en
eV/10^15 atome/cm².

Les resultats fournis sont issus de fit sur des données expérimentales.
Par rapport a la dernière version de Srim, les résultats a faible énergies
sont correct mais ceux a haute énergies sont légèrement sous estimés 20 30 %.
*****************************************************************/
double Ziegler_nuc_stoppow(double V,double Mp,double Zp, double Mt,double Zt)
{
//calcul de l energie du projectile
double Ep=Mp*OD_uma*0.5*V*V/OD_e/1e3; // Energie du projectile en keV

double e=32.53*Mt*Ep/(Zp*Zt*(Mp+Mt)*sqrt(pow(Zp,2.0/3.0)+pow(Zt,2.0/3.0)));
//printf("e;%e\n",e);
double S=0;// puissance d'arret en unite LSS reduite

/*
Routine provenant de :
The electronic and nuclear stopping of energetic ions
------------- IL A ETE CHOISI DE SUPPRIMER CETTE ROUTINE
---- CAR ELLE DONNE LES MEMES Sn QUE CELLE DE SRIM
MAIS AVEC DES ECARTS IMPORTANTS A BASSE ENERGIE!!!
if(e<0.01)
S=1.593*sqrt(e);
else if(e>=0.1 && e<=10.0)
S=1.7*sqrt(e)*log(e+exp(1))/(1+6.8*e+3.4*pow(e,3.0/2.0));
else if(e>10.0)
S=log(0.47*e)/(2.0*e);
*/
if(e<=30)
S=log(1.0+1.1383*e)/(2*(e+0.01321*pow(e,0.21226)+0.19593*pow(e,0.5)));
else
S=log(0.47*e)/(2.0*e);



//conversion en eV/(10^15 a/cm²)
// Correction par rapport au papier de Ziegler
//printf("mult %e \n",(8.462*Zp*Zt*Mp/((Mp+Mt)*sqrt(pow(Zp,2./3.)+pow(Zt,2./3.)))));
S*=(8.462*Zp*Zt*Mp/((Mp+Mt)*sqrt(pow(Zp,2./3.)+pow(Zt,2./3.))));

return S;
}



double **stoppow_Ziegler_proton=NULL; // definition du pointeur qui contiendra les données.
int *stoppow_Ziegler_proton_sze=NULL;// taille de chacun des tableau de donnees
double *stoppow_Ziegler_proton_Z=NULL;
double stoppow_Ziegler_proton_targetnb=0;// nombre de données chargée en mémoire
double setup_Ziegler_elec_stoppow()
{
	//POUR L INSTANT, 16 élements sont integres a la base de donnee.
	stoppow_Ziegler_proton_targetnb=16;
	stoppow_Ziegler_proton_Z=malloc(sizeof(double)*stoppow_Ziegler_proton_targetnb);
	stoppow_Ziegler_proton_sze=malloc(sizeof(int)*stoppow_Ziegler_proton_targetnb);
	stoppow_Ziegler_proton=malloc(sizeof(double*)*stoppow_Ziegler_proton_targetnb);

	//pour l'instant, chargement manuel.
	stoppow_Ziegler_proton_Z[0]=1.0;
	stoppow_Ziegler_proton[0]=univ_loader_xytox(&stoppow_Ziegler_proton_sze[0],"./data/proton_spower/PROT_1.txt");

	stoppow_Ziegler_proton_Z[1]=2.0;
	stoppow_Ziegler_proton[1]=univ_loader_xytox(&stoppow_Ziegler_proton_sze[1],"./data/proton_spower/PROT_2.txt");

	stoppow_Ziegler_proton_Z[2]=6.0;
	stoppow_Ziegler_proton[2]=univ_loader_xytox(&stoppow_Ziegler_proton_sze[2],"./data/proton_spower/PROT_6.txt");

	stoppow_Ziegler_proton_Z[3]=7.0;
	stoppow_Ziegler_proton[3]=univ_loader_xytox(&stoppow_Ziegler_proton_sze[3],"./data/proton_spower/PROT_7.txt");

	stoppow_Ziegler_proton_Z[4]=8.0;
	stoppow_Ziegler_proton[4]=univ_loader_xytox(&stoppow_Ziegler_proton_sze[4],"./data/proton_spower/PROT_8.txt");

	stoppow_Ziegler_proton_Z[5]=10.0;
	stoppow_Ziegler_proton[5]=univ_loader_xytox(&stoppow_Ziegler_proton_sze[5],"./data/proton_spower/PROT_10.txt");

	stoppow_Ziegler_proton_Z[6]=13.0;
	stoppow_Ziegler_proton[6]=univ_loader_xytox(&stoppow_Ziegler_proton_sze[6],"./data/proton_spower/PROT_13.txt");

	stoppow_Ziegler_proton_Z[7]=14.0;
	stoppow_Ziegler_proton[7]=univ_loader_xytox(&stoppow_Ziegler_proton_sze[7],"./data/proton_spower/PROT_14.txt");

	stoppow_Ziegler_proton_Z[8]=18.0;
	stoppow_Ziegler_proton[8]=univ_loader_xytox(&stoppow_Ziegler_proton_sze[8],"./data/proton_spower/PROT_18.txt");

	stoppow_Ziegler_proton_Z[9]=21.0;
	stoppow_Ziegler_proton[9]=univ_loader_xytox(&stoppow_Ziegler_proton_sze[9],"./data/proton_spower/PROT_21.txt");

	stoppow_Ziegler_proton_Z[10]=22.0;
	stoppow_Ziegler_proton[10]=univ_loader_xytox(&stoppow_Ziegler_proton_sze[10],"./data/proton_spower/PROT_22.txt");

	stoppow_Ziegler_proton_Z[11]=26.0;
	stoppow_Ziegler_proton[11]=univ_loader_xytox(&stoppow_Ziegler_proton_sze[11],"./data/proton_spower/PROT_26.txt");

	stoppow_Ziegler_proton_Z[12]=36.0;
	stoppow_Ziegler_proton[12]=univ_loader_xytox(&stoppow_Ziegler_proton_sze[12],"./data/proton_spower/PROT_36.txt");

	stoppow_Ziegler_proton_Z[13]=54.0;
	stoppow_Ziegler_proton[13]=univ_loader_xytox(&stoppow_Ziegler_proton_sze[13],"./data/proton_spower/PROT_54.txt");

	stoppow_Ziegler_proton_Z[14]=79.0;
	stoppow_Ziegler_proton[14]=univ_loader_xytox(&stoppow_Ziegler_proton_sze[14],"./data/proton_spower/PROT_79.txt");

	stoppow_Ziegler_proton_Z[15]=92.0;
	stoppow_Ziegler_proton[15]=univ_loader_xytox(&stoppow_Ziegler_proton_sze[15],"./data/proton_spower/PROT_92.txt");


	return 1;
}

/****************************************************************
Ziegler_elec_stoppow()
-----------------------------------------------------------------
utilité:

Routine permettant de calculer la puissance d'arret électronique
d'un projectile P sur un milieu d'éléments Zt. La routine retourne
la puissance d'arret en kev/(mg/cm²).

Il est nécessaire de charger les données concernant lepouvoir d'arret
des protons de PSTAR pour être capable d'obtenir des résultats
-----------------------------------------------------------------
Provenance:

The electronic and nuclear stopping of energetic ions
J.F. Ziegler
Apllied physics Letters 31, 544(1977)
------------------------------------------------------------------
La fonction prend en paramètre:

V la vitesse du projectile (en m.s)
Mp la masse du projectile (en uma)
Zp la charge du projectile (en charge élémentaire)
Mt la masse de la cible (en uma)
Zt la charge de la cible (en charge élémentaire)

La fonction retourne la puissance d'arret électronique en
keV/mg/cm².

Les resultats fournis sont issus de fit sur des données expérimentales.
En conjonction avec la puissance d'arret nucleaire de Ziegler,
Ziegler_nuc_stoppow, l'écart entre la puissance d'arret
experimentale et le résultat de ces fonctions est de 4.5%.
*****************************************************************/
double Ziegler_elec_stoppow(double V,double Mp,double Zp,double Zt)
{

	double s=0;// p arret reduite
	double sp;//p arret proton
	double E;//energie du projectile, en MeV
	double v1=0;
	double v2=0;
	int p_spow_ind=-1;


	// Première étape, calcul de la puissance d'arret électronique réduite
	v1=0.886*V/OD_v0*pow(Zp,-2.0/3.0);
	v2=v1+0.0378*sin(0.5*M_PI*v1);
	s=1-exp(-v2)*(1.034-0.1777*exp(-0.08114*Zp));
	s*=s;
    //printf("s reduit %e\n",s);
	//Deuxieme etape, recherche de la puissance d'arret electronique des protons correspondant.
	E=0.5*OD_uma*V*V/(OD_e*1.0e6); // calcul de l energie du projectile en MeV;


   /* printf("V1 %e V2: %e, terme1 %e terme 2 %e \n",v1,v2,exp(-v2),(1.034-0.1777*exp(-0.08114*Zp)));
    printf("sred %e\n",s);*/
	// recherche de l'existence de la puissance d'arret demande.
	if(stoppow_Ziegler_proton==NULL) return NAN; // pas de données chargées

	for(int i=0;i<stoppow_Ziegler_proton_targetnb;i++)
    {

		if(stoppow_Ziegler_proton_Z[i]==Zt){p_spow_ind=i;break;}
    }

	if(p_spow_ind==-1) return NAN; // pas de donnee existante

        sp=linear_interp_xy_x(stoppow_Ziegler_proton[p_spow_ind],stoppow_Ziegler_proton_sze[p_spow_ind],E);
     //   printf("E %e sp: %e\n",E,sp);

	return sp*s*Zp*Zp;
}



/****************************************************************
free_Ziegler_elec_stoppow()
-----------------------------------------------------------------
utilité:

Routine permettant de décharger les données concernant le pouvoir
d'arrêt des protons dans divers milieux (PSTAR).

-----------------------------------------------------------------
Provenance:

Les données déchargées proviennent à l'heure actuelle de:

*****************************************************************/
int free_Ziegler_elec_stoppow()
{

	if(stoppow_Ziegler_proton==NULL)return 0;

	for (int i=0;i<stoppow_Ziegler_proton_targetnb;i++)
	free(stoppow_Ziegler_proton[i]);
	free(stoppow_Ziegler_proton);
	stoppow_Ziegler_proton=NULL;

	free(stoppow_Ziegler_proton_Z);
	free(stoppow_Ziegler_proton_sze);
	stoppow_Ziegler_proton_targetnb=0;

	return 1;
}
