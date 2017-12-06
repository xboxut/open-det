


/****************************************************************
q_Schiwietz_Grande()
-----------------------------------------------------------------
utilité:

Routine de calcul d'une charge moyenne provenant du papier:
-----------------------------------------------------------------
Provenance:
Improvised charge-state formulas
G. Schiwietz P.L. Grande
Nuclear intruments and method in physic research B 175 177 (2001)

------------------------------------------------------------------
La fonction prend en paramètre:
vp la vitesse du projectile en cm/s
Zp la charge nucleaire du projectile
Zt la charge nucleaire de la cible

ATTENTION:
Elle est valide pour les  projectiles ayant un z compris entre 1
et 92 pour des cibles GAZEUSES avec un Z entre  1 et 54.
*****************************************************************/
double q_Schiwietz_Grande(double vp, double Zp, double Zt);
