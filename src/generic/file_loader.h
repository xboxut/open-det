


#define UNIV_LOADER_MAX_ARRAY_LEN   32768 // Taille maximum gérée par la routine


/*****************************************************************
double * universal_loader(int *size,const char *path)
------------------------------------------------------------------
utilité:
Fonction pour charger un fichier de données quelque soit son type 
de séparation et de mise en forme.
Il faut cependant respecter la suite de chiffres x, y.

Les char 123456789 . e E sont réservés
-------------------------------------------------------------
la fonction prend en paramètre:
int *size        pointeur vers un entier ou sera stockée la taille
                 du tableau
cont char *path  pointeur vers la chaine de char contenant le nom du 
                 fichier
				 
				 
La fonction retourne un pointeur vers un tableau 1d contenant
une alternance de x et de y.
*****************************************************************/
double * univ_loader_xytox(int *size,const char *path);


/*****************************************************************
double ** universal_loader(int *size,const char *path)
------------------------------------------------------------------
utilité:
Fonction pour charger un fichier de données quelque soit son type 
de séparation et de mise en forme.
Il faut cependant respecter la suite de chiffres x, y.

Les char 123456789 . e E sont réservés
-------------------------------------------------------------
la fonction prend en paramètre:
int *size        pointeur vers un entier ou sera stockée la taille
                 du tableau
cont char *path  pointeur vers la chaine de char contenant le nom du 
                 fichier
				 
				 
La fonction retourne un pointeur vers un tableau 2d contenant
deux ligne: x et y.
*****************************************************************/
double ** univ_loader_xytoxy(int *size,const char *path);