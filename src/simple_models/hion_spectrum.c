


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#ifdef USE_GSL
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#endif // USE_GSL

#include "./../generic/OD_vect3d.h"


/****************************************************************
hion_analytic_spec_dep_o()
-----------------------------------------------------------------
utilit�:

Cette routine permet de calculer le spectre d'ions �nerg�tiques
en sortie d'un d�pot de profondeur L.


-----------------------------------------------------------------
Provenance:
CEA Cadarache, DEN/DER/SPEx/LPE
G. de Izarra et A. Amamra
------------------------------------------------------------------
La fonction prend en param�tre:

E   l'�nergie d'inter�t (l'unit� est libre mais doit etre similaire
    a celle de E0)
L   L'�paisseur du d�pot (en cm).
E0  L'�nergie initiale de l'ion.
RE0 Le range d'un ion � E0 dans le mat�riau consid�r� (en cm).
n   une puissance pour la loi de mod�ration empirique E=E0(1-x/R)^n
rrate le taux de r�action � l'origine des ions �nerg�tiques en (cm^3)


*****************************************************************/
double hion_analytic_spec_dep_o(double E,double L, double E0, double RE0,double n,double rrate)
{
    double Lmax=(1.0-pow(E/E0,(1.0/n)))*RE0;
    double S=0;
    if(Lmax>L)
    Lmax=L;

//    y=lmax^2/2*R*E0^(-1/n)*1/n*E^(1/n-1)/(R*(1-(E/E0)^(1/n)))^2;
    S=Lmax*Lmax/2.0*RE0*pow(E0,-1.0/n)*1.0/n*pow(E,1.0/n-1.0)/(pow((RE0*(1.0-pow(E/E0,1.0/n))),2.0))*1.0/L*rrate;

 return S;
}


#ifdef USE_GSL
/****************************************************************
hion_mc0_spec_dep_o()
-----------------------------------------------------------------
utilit�:

La fonction calcule un spectre en �nergie et en angle des ions en
sortie d'un d�p�t semi infini d'�paisseur L.

-----------------------------------------------------------------
Provenance:
CEA Cadarache, DEN/DER/SPEx/LPE
G. de Izarra et A. Amamra
------------------------------------------------------------------
La fonction prend en param�tre:

r       Un g�n�rateur de nombre al�atoire.
L       L'�paisseur du d�pot (en cm).
E0      L'�nergie initiale de l'ion.
RE0     Le range d'un ion � E0 dans le mat�riau consid�r� (en cm).
n       une puissance pour la loi de mod�ration empirique E=E0(1-x/R)^n
rrate   le taux de r�action � l'origine des ions �nerg�tiques en (cm^3)
dE      la taille d'un bin en �nergie.
dang    La taille d'un bin angulaire. Attention, le code fournit le cos
        de l'angle � la normale de la surface de sortie(z). Car on consid�re
        que l'�mission est homog�ne dans le plan xy (milieu semi infini).


*****************************************************************/
double ** hion_mc0_spec_dep_o(gsl_rng * r,int N,double E0, double L, double RE0,double n,double rrate,double dE,double dang,double *eprob )
{

    int nE, nang;
    double **tabular;

    double z0,x,y,z; // z0 position ou l'ion est g�n�r� dans le depot. x,y,z, vec directeur
    double norm,k; // normalisation du vec directeur de l ion + k -> distance
    double E;// var pour stocker l energie d'un ion en sortie de d�pot

    int Eind=0,angind=0;// indcies pour remplir le tableau;
    int cnt=0; // compteur d'ion sortant du d�pot

    // initialisation du tableau de r�sultat.
    nE=(int)ceil(E0/dE);
    nang=(int)ceil(1/dang);

    tabular=malloc(sizeof(double *)*nE);
    for(int i=0;i<nE;i++)
    {
        tabular[i]=malloc(sizeof(double)*nang);
        for(int j=0;j<nang;j++)
        tabular[i][j]=0;

    }



    for(int i=0;i<N;i++)
    {


    z0=L*gsl_rng_uniform(r);

    x=(2.0*gsl_rng_uniform(r)-1); // Emission sur 2pi sr
    y=(2.0*gsl_rng_uniform(r)-1);
    z=gsl_rng_uniform(r);

    norm=sqrt(x*x+y*y+z*z);

    z/=norm;
    x/=norm;
    y/=norm;


            if(z!=0)
            k=(L-z0)/z;
            else
                k=100000000000000.0; // arret du fragment de fission dans le d�pot-> trajectoire // � la sortie du d�pot

            E=E0*pow((1.0-k/RE0),n);
            if(E>0)
            {
                //Calcul du bin d energie et d angle que l on va incr�menter
                Eind=(int)floor(E/(dE+0.000001));
                angind=(int)floor(z/(dang+0.000001));
                tabular[Eind][angind]+=1.0;
                cnt++;
            }

    }

    // Normalisation du spectre
    for(int i=0;i<nE;i++)
    {
        for(int j=0;j<nang;j++)
        {
            tabular[i][j]*=rrate/(double)N;
        }

    }

    // Si l'utilisateur le demande, calcul de la proba totale de sortie de l'ion.
    if(eprob!=NULL)
        *eprob=(double) cnt/(double) N;


    return tabular;
}


/****************************************************************
hion_analytic_spec_dep_o()
-----------------------------------------------------------------
utilit�:

Cette routine tire al�atoirement un ion �nerg�tique g�n�r�
dans un d�pot semi-infini d'�paisseur L et le propage jusqu'� ce
qu'il s'arrete ou bien atteigne la surface du d�pot.

Le mod�le utilise des lois de mod�ration ou:
- la trajectoire de l'ion est consid�r� rectiligne.
-La perte d'�nergie est continue et suit une loi semi empirique.
- le milieu est semi infini d'une �paisseur L inf�rieur � R le range
  d'un ion dans ce milieu.


-----------------------------------------------------------------
Provenance:
CEA Cadarache, DEN/DER/SPEx/LPE
G. de Izarra et A. Amamra
------------------------------------------------------------------
La fonction prend en param�tre:

r   un g�n�rateur de nombre aleatoire de la GSL
E   l'�nergie d'inter�t (l'unit� est libre mais doit etre similaire
    a celle de E0)
L   L'�paisseur du d�pot (en cm).
E0  L'�nergie initiale de l'ion.
RE0 Le range d'un ion � E0 dans le mat�riau consid�r� (en cm).
n   une puissance pour la loi de mod�ration empirique E=E0(1-x/R)^n

La fonction retourne l'�nergie de  l'ion en sortie de d�p�t
ou bien -1.0 si ce dernier est stopp� dans le d�pot.
Pour r�cuperer la direction de propagation de l'ion, on peut
passer en param�tre un pointeur sur un vecteur. Si l'argument est
mis a NULL, rien n'est copi�.
*****************************************************************/
double hion_mc1_spec_dep_o(gsl_rng * r,double E0, double L, double RE0, double n, OD_vect3d *dir)
{

    double z0;//profondeur a laquelle un fragment de fission est cree.
    double x,y,z; //vecteur directeur de ff
    double norm;
    double k=0; // distance parcourue dans le depot
    double E;
    z0=L*gsl_rng_uniform(r);

    x=(2.0*gsl_rng_uniform(r)-1); // Emission sur 2pi sr
    y=(2.0*gsl_rng_uniform(r)-1);
    z=gsl_rng_uniform(r);

    norm=sqrt(x*x+y*y+z*z);

    z/=norm;
    x/=norm;
    y/=norm;


            if(z!=0)
            k=(L-z0)/z;
            else
                k=100000000000000.0; // arret du fragment de fission dans le d�pot-> trajectoire // � la sortie du d�pot

            E=E0*pow((1.0-k/RE0),n);

    if(dir!=NULL)
    {
        dir->x=x;
        dir->y=z;
        dir->z=z;

    }

    if(E>0)
    {
        return E;
    }
    else
    {
     return -1.0;
    }
}
#endif



int test_hion_spectrum()
{

    //Donn�es pour les fragments de fission lourds


    double R=10.4e-6*100.0;
    double E0=68e6;
    double n=1.45;
    double rrate=1.0;
    double L=3.0e-6*100.0;
    double **res;
    OD_vect3d vect_dir;

    double E;
    FILE *fich=NULL;


    printf("**********MODELE DE TRANSPORT D ION EN SORTIE DE DEPOT ACTIF ***\n\n\n");

    printf("Cas des fragments de fission lourds les plus probables:\n");
    printf("E0: %e\nREO: %e\nn: %e\nL= %e\nrrate: %e\n\n",E0,R,n,L,rrate);
    printf("Cas des fragments de fission lourds les plus probables:\n");

    printf("*******Modele analytique*****\n");
    printf("creation du fichier 'TESThion_analytic_spec_dep_o.txt'\n");

    fich=fopen("TESThion_analytic_spec_dep_o.txt","w");
    if(fich==NULL)
    {
        printf("IMPOSSIBLE DE CREER LE FICHIER\nSORTIE DU TEST\n");
        return 0;
    }

    for(int i=1;i<500;i++)
        fprintf(fich,"%e %e\r\n",(double)i*68.e6/500.0,hion_analytic_spec_dep_o((double)i*68.e6/500.0,L, E0,  R,n, rrate));

    fclose(fich);

    #ifdef USE_GSL
    printf("*******Modele monte_carlo one shot*****\n");

    printf("creation du fichier 'TESThion_mc1_spec_dep_o.txt'\n");
    fich=fopen("TESThion_mc1_spec_dep_o.txt","w");
    if(fich==NULL)
    {
        printf("IMPOSSIBLE DE CREER LE FICHIER\nSORTIE DU TEST\n");
        return 0;
    }


    printf("Tirage de 30000 fragments de fission lourds\n");

    // init du gene de nombre aleatoire

        gsl_rng * r = gsl_rng_alloc (gsl_rng_mt19937);
        gsl_rng_set(r,time(NULL));

    for(int i=0;i<30000;i++)
    {


    E=hion_mc1_spec_dep_o(r,E0,L,R, n, &vect_dir);
    if(i%10==0)
    printf("%d E:%e x:%e y:%e z:%e\n",i, E,vect_dir.x,vect_dir.y,vect_dir.z);
    fprintf(fich,"%e %e %e %e\r\n", E,vect_dir.x,vect_dir.y,vect_dir.z);

    }

    fclose(fich);


     printf("*******Modele monte_carlo spectre de sortie total*****\n");

        printf("creation du fichier 'TESThion_mc0_spec_dep_o.txt'\n");
        printf("construction du tableau de l'energie de sortie.\n pas d'energie: %e, pas d'angle %e\n",0.5e6,0.005);
    fich=fopen("TESThion_mc0_spec_dep_o.txt","w");
    if(fich==NULL)
    {
        printf("IMPOSSIBLE DE CREER LE FICHIER\nSORTIE DU TEST\n");
        return 0;
    }

    res=hion_mc0_spec_dep_o(r,1e7,E0, L, R,n,1,0.5e6,0.005,NULL );

    for(int i=0;i<(int)ceil(E0/0.5e6);i++)
    {

        for(int j=0;j<(int)ceil(1/0.005);j++)
        {
            fprintf(fich,"%e %e %e\r\n",(double)i*0.5e6,(double)j*0.005,res[i][j]);
        }
    }

    fclose(fich);
    #endif
    printf("fin du test des modeles concernant le spectre de sortie des fragments de fission.\n");
    return 1;
}

