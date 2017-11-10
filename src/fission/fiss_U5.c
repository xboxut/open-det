
#include <stdlib.h>
#include <math.h>





/****************************************************************
fiss_U5_get_rnd_ff()
-----------------------------------------------------------------
utilité:

Routine pour le calcul d'une masse et d'une énergie d'un fragment
de fission emis par un atome d'U 235.

la masse du fragment est tiré de manière aléatoire, puis l'energie
cinétique est calculée en considérant un nombre de neutrons moyen
émis de x et une énergie embarquée de xx MeV.
-----------------------------------------------------------------
Provenance:
Rendement de la fission de l'U5: JEFF3.1.1
 A. Koning, R. Forrest, M. Kellett, R. Mills, H. Henriksson, Y. Rugama,
 The JEFF-3.1 Nuclear Data Library, JEFF Report 21, OECD/NEA, Paris,
 France, 2006, ISBN 92-64-02314-3.

Nombre moyen de neutron :  IAEA-CRP-STD
S. A. Badikov, C. Zhenpeng, A. D. Carlson, E. V. Gai, G. M. Hale,
F.-J. Hambsh, H. M. Hofmann, T. Kawano, N. M. Larson, V. G. Pronyaev,
 D.L. Smith, Soo-Youl Ho, S. Tagesen, H. K. Vonach, H.L. Nichols,
 IAEA CRP "International Evaluation of Neutron Cross-Section Standards",
 IAEA Scientific and Technical Report STI/PUB/1292, November 2007,
 International Atomic Energy Agency, Vienna, Austria, ISBN 90-0-100807-4.
http://www-nds.iaea.org/standards/

Energie des ff:
Kaye and Laby  table of physical and chemical constant
http://www.kayelaby.npl.co.uk/atomic_and_nuclear_physics/4_7/4_7_1.html
------------------------------------------------------------------
La fonction prend en paramètre:
urnd un nombre aleatoire genere selon une distribution uniforme
mass un pointeur vers un double qui contiendra la masse du fragment
E un pointeur vers un double qui contiendra l'énergie du fragment.


La fonction stocke la masse, l'énergie et le Z des fragments de
fission dans les variables pointées.

Elle retourne 0 en cas de problème et 1 si le calcul est un succès.

*****************************************************************/
int fiss_U5_get_rnd_ff(double urnd,double *mass, double * E, double *Z)
{
	//static const m_emit_neutron=2.43;
	// Energie totale emportée par les ff 169.1 MeV
	static const double cumpdf[119]={8.555000e-006, 1.275500e-005, 6.675500e-005, 9.167550e-004,
	9.300950e-004, 9.308242e-004, 9.314348e-004, 9.366358e-004, 9.367619e-004, 9.375509e-004,
	9.376773e-004, 9.381773e-004, 9.381774e-004, 9.381774e-004, 9.381775e-004, 9.381778e-004, 9.381791e-004, 9.381839e-004,
	9.382004e-004, 9.382550e-004, 9.384251e-004, 9.389273e-004, 9.406260e-004, 9.447743e-004, 9.686837e-004, 1.011117e-003, 1.113359e-003, 1.356929e-003, 1.994788e-003, 2.985319e-003,
	4.627599e-003, 7.406698e-003, 1.239778e-002, 1.903152e-002, 2.884564e-002, 4.151780e-002, 5.900393e-002, 8.283769e-002, 1.119289e-001, 1.412922e-001, 1.712741e-001, 2.029436e-001, 2.354340e-001,
	2.681598e-001, 2.996900e-001, 3.293825e-001, 3.582464e-001, 3.890511e-001, 4.202980e-001, 4.461324e-001, 4.675666e-001, 4.830821e-001, 4.924539e-001, 4.971886e-001, 4.992424e-001, 4.999405e-001,
	5.002269e-001, 5.003712e-001, 5.004983e-001, 5.005967e-001, 5.006559e-001, 5.007357e-001, 5.008002e-001, 5.008569e-001, 5.009372e-001, 5.009987e-001, 5.010667e-001, 5.011419e-001, 5.012147e-001,
	5.012777e-001, 5.013679e-001, 5.014432e-001, 5.016013e-001, 5.017315e-001, 5.020285e-001, 5.026297e-001, 5.042828e-001, 5.078132e-001, 5.167071e-001, 5.311004e-001, 5.525819e-001, 5.855512e-001,
	6.243486e-001, 6.575813e-001, 6.892231e-001, 7.211460e-001, 7.548921e-001, 7.868558e-001, 8.184602e-001, 8.477510e-001, 8.769230e-001, 9.067805e-001, 9.341401e-001, 9.539176e-001, 9.688564e-001,
	9.800169e-001, 9.884149e-001, 9.936835e-001, 9.969391e-001, 9.990411e-001, 1.000304e+000, 1.001043e+000, 1.001406e+000, 1.001560e+000, 1.001626e+000, 1.001659e+000, 1.001669e+000, 1.001674e+000,
	1.001676e+000, 1.001676e+000, 1.001676e+000, 1.001676e+000, 1.001676e+000, 1.001676e+000, 1.001676e+000, 1.001676e+000, 1.001676e+000, 1.001676e+000, 1.001676e+000 };


	if(mass==NULL || E==NULL || Z==NULL)return 0;

	for(int i=0;i<119;i++)
	{
		if(cumpdf[i]>=urnd)
		{

			if(i==0)*mass=1.0; // pour eviter les problèmes avec i=0
			else{
						if(fabs(cumpdf[i]-urnd)<fabs(cumpdf[i-1]-urnd)) // calcul de la masse la plus proche.
						*mass=(double)(i+1);
						else
						*mass=(double)i;
			}
			break;

		}
	}

	*E=234.0/(*mass)*169.1e6;

	*Z=round(*mass/234.0*92.0);

	return 1;
}
