

#include <stdio.h>
#include <stdlib.h>

#include "OD_surf.h"



/*...............................................................
//DONNEE
OD_geom_surf * base_geom_surf
-----------------------------------------------------------------
utilit�:

D�finition d'un point d'ancrage pour la liste chain�e
de surface que l'on peut d�clarer dans un jeu de donn�es.

*/

OD_geom_surf * base_geom_surf=NULL;






/****************************************************************
int free_base_geom_surf()
-----------------------------------------------------------------
utilit�:
Routine permettant de supprimer l'ensemble des surfaces g�n�r�es
lors d'une simulation. Elle retourne le nombre de surface supprim�e
.
Par d�faut, la routine part du point d'ancrage "base_geom_surf"
de la liste chain�e de surface
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