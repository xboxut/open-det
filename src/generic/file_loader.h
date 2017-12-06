


#define UNIV_LOADER_MAX_ARRAY_LEN   32768 // Taille maximum g�r�e par la routine


/*****************************************************************
double * universal_loader(int *size,const char *path)
------------------------------------------------------------------
utilit�:
Fonction pour charger un fichier de donn�es quelque soit son type 
de s�paration et de mise en forme.
Il faut cependant respecter la suite de chiffres x, y.

Les char 123456789 . e E sont r�serv�s
-------------------------------------------------------------
la fonction prend en param�tre:
int *size        pointeur vers un entier ou sera stock�e la taille
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
utilit�:
Fonction pour charger un fichier de donn�es quelque soit son type 
de s�paration et de mise en forme.
Il faut cependant respecter la suite de chiffres x, y.

Les char 123456789 . e E sont r�serv�s
-------------------------------------------------------------
la fonction prend en param�tre:
int *size        pointeur vers un entier ou sera stock�e la taille
                 du tableau
cont char *path  pointeur vers la chaine de char contenant le nom du 
                 fichier
				 
				 
La fonction retourne un pointeur vers un tableau 2d contenant
deux ligne: x et y.
*****************************************************************/
double ** univ_loader_xytoxy(int *size,const char *path);