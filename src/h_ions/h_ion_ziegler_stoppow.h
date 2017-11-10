





/****************************************************************
Ziegler_nuc_stoppow()
-----------------------------------------------------------------
utilité:

Routine permettant de calculer la puissance d'arret nucléaire
d'un projectile P sur un milieu d'éléments Zp. La routine retourne
la puissance d'arret en ev/(10^15 atome/cm²)
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

La fonction retourne la puissance d'arret nucléaire en
eV/10^15 atome/cm².

Les resultats fournis sont issus de fit sur des données expérimentales.
En conjonction avec la puissance d'arret électronique de Ziegler,
Ziegler_elec_stoppow, l'écart relation entre la puissance d'arret
experimentale et le résultat de ces fonctions est de 4.5%.
*****************************************************************/
double Ziegler_nuc_stoppow(double V,double Mp,double Zp, double Mt,double Zt);




double setup_Ziegler_elec_stoppow();

double Ziegler_elec_stoppow(double V,double Mp,double Zp,double Zt);

int free_Ziegler_elec_stoppow();
