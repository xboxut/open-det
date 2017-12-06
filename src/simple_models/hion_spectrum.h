




double hion_analytic_spec_dep_o(double E,double L, double E0, double RE0,double n,double rrate);



#ifdef USE_GSL

/****************************************************************
hion_mc0_spec_dep_o()
-----------------------------------------------------------------
utilité:

La fonction calcule un spectre en énergie et en angle des ions en
sortie d'un dépôt semi infini d'épaisseur L.

-----------------------------------------------------------------
Provenance:
CEA Cadarache, DEN/DER/SPEx/LPE
G. de Izarra et A. Amamra
------------------------------------------------------------------
La fonction prend en paramètre:

r       Un générateur de nombre aléatoire.
L       L'épaisseur du dépot (en cm).
E0      L'énergie initiale de l'ion.
RE0     Le range d'un ion à E0 dans le matériau considéré (en cm).
n       une puissance pour la loi de modération empirique E=E0(1-x/R)^n
rrate   le taux de réaction à l'origine des ions énergétiques en (cm^3)
dE      la taille d'un bin en énergie.
dang    La taille d'un bin angulaire. Attention, le code fournit le cos
        de l'angle à la normale de la surface de sortie(z). Car on considère
        que l'émission est homogène dans le plan xy (milieu semi infini).


*****************************************************************/
double ** hion_mc0_spec_dep_o(gsl_rng * r,int N,double E0, double L, double RE0,double n,double rrate,double dE,double dang,double *eprob );


/****************************************************************
hion_analytic_spec_dep_o()
-----------------------------------------------------------------
utilité:

Cette routine tire aléatoirement un ion énergétique généré
dans un dépot semi-infini d'épaisseur L et le propage jusqu'à ce
qu'il s'arrete ou bien atteigne la surface du dépot.

Le modèle utilise des lois de modération ou:
- la trajectoire de l'ion est considéré rectiligne.
-La perte d'énergie est continue et suit une loi semi empirique.
- le milieu est semi infini d'une épaisseur L inférieur à R le range
  d'un ion dans ce milieu.


-----------------------------------------------------------------
Provenance:
CEA Cadarache, DEN/DER/SPEx/LPE
G. de Izarra et A. Amamra
------------------------------------------------------------------
La fonction prend en paramètre:

r   un générateur de nombre aleatoire de la GSL
E   l'énergie d'interêt (l'unité est libre mais doit etre similaire
    a celle de E0)
L   L'épaisseur du dépot (en cm).
E0  L'énergie initiale de l'ion.
RE0 Le range d'un ion à E0 dans le matériau considéré (en cm).
n   une puissance pour la loi de modération empirique E=E0(1-x/R)^n

La fonction retourne l'énergie de  l'ion en sortie de dépôt
ou bien -1.0 si ce dernier est stoppé dans le dépot.
Pour récuperer la direction de propagation de l'ion, on peut
passer en paramètre un pointeur sur un vecteur. Si l'argument est
mis a NULL, rien n'est copié.
*****************************************************************/
double hion_mc1_spec_dep_o(gsl_rng * r,double E0, double L, double RE0, double n, OD_vect3d *dir);


#endif // USE_GSL


int test_hion_spectrum();
