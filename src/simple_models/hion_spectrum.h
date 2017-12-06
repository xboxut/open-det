




double hion_analytic_spec_dep_o(double E,double L, double E0, double RE0,double n,double rrate);



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
double ** hion_mc0_spec_dep_o(gsl_rng * r,int N,double E0, double L, double RE0,double n,double rrate,double dE,double dang,double *eprob );


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
double hion_mc1_spec_dep_o(gsl_rng * r,double E0, double L, double RE0, double n, OD_vect3d *dir);


#endif // USE_GSL


int test_hion_spectrum();
