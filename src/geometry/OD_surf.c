

#include <stdio.h>
#include <stdlib.h>

#include "OD_surf.h"



/*...............................................................
//DONNEE
OD_geom_surf * base_geom_surf
-----------------------------------------------------------------
utilité:

Définition d'un point d'ancrage pour la liste chainée
de surface que l'on peut déclarer dans un jeu de données.

*/

OD_geom_surf * base_geom_surf=NULL;






/****************************************************************
int free_base_geom_surf()
-----------------------------------------------------------------
utilité:
Routine permettant de supprimer l'ensemble des surfaces générées
lors d'une simulation. Elle retourne le nombre de surface supprimée
.
Par défaut, la routine part du point d'ancrage "base_geom_surf"
de la liste chainée de surface
-----------------------------------------------------------------
Provenance:
CEA Cadarache, DEN/DER/SPEx/LPE
G. de Izarra 
------------------------------------------------------------------


*****************************************************************/
int free_base_geom_surf()
{
	OD_geom_surf * current=base_geom_surf;
	OD_geom_surf * next=NULL;
	int surf_nb=0;
	
	while(current!=NULL)
	{
		next=current->next;
		free(current);
		current=next;
		surf_nb++
	}
		return surf_nb;
}