
#include <math.h>
#include <stdio.h>

#include "interpolation.h"




double linear_interp_xy_x(double *data,int data_sze,double x)
{
	int i=0;
	for(;i<data_sze;i+=2)
    {
		if(x<data[i])break;
   // printf("I %d   %e \n",i,data[i]);
    }
	// Si l'interpolation n'est pas dans les bornes demandées -> on extrapole pas on sort avec un message d erreur.
	if(i==0 || i==data_sze-2)//
	{
		     return NAN;
	}


	return (data[i+1]-data[i-1])/(data[i]-data[i-2])*(x-data[i-2])+data[i-1];
}
