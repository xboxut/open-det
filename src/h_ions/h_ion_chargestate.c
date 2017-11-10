

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "./../generic/constant.h"

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
vp la vitesse du projectile en m/s
Zp la charge nucleaire du projectile
Zt la charge nucleaire de la cible

ATTENTION:
Elle est valide pour les  projectiles ayant un z compris entre 1
et 92 pour des cibles GAZEUSES avec un Z entre  1 et 54.
*****************************************************************/
double q_Schiwietz_Grande(double vp, double Zp, double Zt)
{


    double q;
    double x=vp/OD_v0*pow(Zp,-0.52)*pow(Zt,0.03-0.017*pow(Zp,-0.52)*vp/OD_v0);

    x=pow(x,1.0+0.4/Zp);
    q=Zp*(376.0*x+pow(x,6.0))/(1428.-1206.*pow(x,0.5)+690.0*x+pow(x,6.0));
    return q;
}
