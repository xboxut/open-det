#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

#include "file_loader.h"

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
double * univ_loader_xytox(int *size,const char *path)
{

    char *stringloaded;//pointeur qui pointera sur le fichier une fois chargé en memoire
    int  nb_cnt=0;//compteur du nombre de points de mesures présent dans le fichier
    char buffer[32];//buffer pour isoler chacun des nombres du fichier en vue de le convertir en double
    char *ptr_fich;//pointeur de recopie placé sur le fichier
    char *ptr_buff;//pointeur de recopie placé sur le buffer
    double *tab;//pointeur sur tableau qui contiendra l'ensemble des points
    double *return_tab;//pointeur sur le tableau que l'on retournera

    bool copy=false;
    int i;//variable utilisée dans les for

    /*****************************************************
    Chargement du fichier en mémoire
    *****************************************************/

    FILE *temp=fopen(path,"rb");// ouverture du fichier en mode binaire
    if(temp==NULL)return NULL;
	fseek(temp,0,SEEK_END);// recherche de la fin du fichier

    int charlen=ftell(temp); // récupération du nombre d'octet entre le
                             //début et la fin du fichier

    fseek(temp,0,SEEK_SET);// On se replace au début du fichier

	stringloaded=malloc(sizeof(char)*(charlen+1));
    // construction d'une chaine de char qui contiendra l'ensemble du fichier

	fread(stringloaded,charlen,1,temp);
    //chargement du fichier dans la chaine nouvellement crée

    fclose(temp);

    stringloaded[charlen]='\0';


    //printf(stringloaded);
    //Placement de la fin de la chaine à la place de EOF.

    /********************************************************
    traitement du fichier
    ********************************************************/

     tab=malloc(sizeof(double)*2*UNIV_LOADER_MAX_ARRAY_LEN);//création du tableau de double


     ptr_fich=stringloaded;
     ptr_buff=buffer;

     do// on boucle tant que tout le fichier n'est pas traité
     {
          //Si le charactere pointé par le pointeur est un nombre
          // ou une "puissance"
          if(*ptr_fich=='#')
          {
              ptr_fich++;

              while(*ptr_fich!='#')
              {
                  ptr_fich++;
              }
          }
          if((*ptr_fich>='0' && *ptr_fich<='9')|| *ptr_fich=='.'
                || *ptr_fich=='e' || *ptr_fich=='E' || *ptr_fich=='-'|| *ptr_fich=='+')
          {
             copy=true; // on signale que la copie de nombre est en marche
             *ptr_buff=*ptr_fich; //on copie chaque chiffre dans le buffer
             ptr_buff++;//on incremente le pointeur du buff se déplacer dans la chaine

          }
          else //sinon
          {
             if(copy==true)//si l'on a plus rien a copier mais que l'on vient de le faire
             {
                  copy=false;// on arrete la copie
                  *ptr_buff='\0';//on cloture la chaine contenue dans buffer

                  tab[nb_cnt]=atof(buffer);// et on la transforme en double
                //  printf("%s %e\n",buffer,tab[nb_cnt]);
                  //printf("\n%d\n",nb_cnt);
                  nb_cnt++;//on incrémente le nombre de nombre déjà convertis
                  ptr_buff=buffer;// on replace le pointeur au début du buffer pour réécrire dessus
             }
          }

          ptr_fich++;
          //printf("%c\n",*ptr_fich);
     }while(*(ptr_fich-1)!='\0');


     //On dispose maintenant d'un tableau en 1d contenant les x et les y
     // On va trier les x et les y en utilisant un nouveau tableau que
     // l'on va ensuite pouvoir retourner.


     return_tab=malloc(sizeof(double)*nb_cnt);

     for( i=0;i<nb_cnt;i++)//trie des x et des y et stockage de ces dernier dans un tableau
     {
     return_tab[i]=tab[i];
     }

     //avant de terminer la fonction, on spécifie à l'utilisateur le nb de points dans
     // le tableau
   if(size!=NULL)
    *size=nb_cnt;

    free(tab);
    return return_tab;
}



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
double ** univ_loader_xytoxy(int *size,const char *path)
{

    char *stringloaded;//pointeur qui pointera sur le fichier une fois chargé en memoire
    int  nb_cnt=0;//compteur du nombre de points de mesures présent dans le fichier
    char buffer[32];//buffer pour isoler chacun des nombres du fichier en vue de le convertir en double
    char *ptr_fich;//pointeur de recopie placé sur le fichier
    char *ptr_buff;//pointeur de recopie placé sur le buffer
    double *tab;//pointeur sur tableau qui contiendra l'ensemble des points
    double **return_tab;//pointeur sur le tableau que l'on retournera

    bool copy=false;
    int i;//variable utilisée dans les for

    /*****************************************************
    Chargement du fichier en mémoire
    *****************************************************/

    FILE *temp=fopen(path,"rb");// ouverture du fichier en mode binaire

	fseek(temp,0,SEEK_END);// recherche de la fin du fichier

    int charlen=ftell(temp); // récupération du nombre d'octet entre le
                             //début et la fin du fichier

    fseek(temp,0,SEEK_SET);// On se replace au début du fichier

	stringloaded=malloc(sizeof(char)*(charlen+1));
    // construction d'une chaine de char qui contiendra l'ensemble du fichier

	fread(stringloaded,charlen,1,temp);
    //chargement du fichier dans la chaine nouvellement crée

    fclose(temp);

    stringloaded[charlen]='\0';


    //printf(stringloaded);
    //Placement de la fin de la chaine à la place de EOF.

    /********************************************************
    traitement du fichier
    ********************************************************/

     tab=malloc(sizeof(double)*2*UNIV_LOADER_MAX_ARRAY_LEN);//création du tableau de double


     ptr_fich=stringloaded;
     ptr_buff=buffer;

     do// on boucle tant que tout le fichier n'est pas traité
     {
         if(*ptr_fich=='#')
          {
              ptr_fich++;

              while(*ptr_fich!='#')
              {
                  ptr_fich++;
              }
          }
          //Si le charactere pointé par le pointeur est un nombre
          // ou une "puissance"
          if((*ptr_fich>='0' && *ptr_fich<='9')|| *ptr_fich=='.'
                || *ptr_fich=='e' || *ptr_fich=='E' || *ptr_fich=='-' || *ptr_fich=='+')
          {
             copy=true; // on signale que la copie de nombre est en marche
             *ptr_buff=*ptr_fich; //on copie chaque chiffre dans le buffer
             ptr_buff++;//on incremente le pointeur du buff se déplacer dans la chaine

          }
          else //sinon
          {
             if(copy==true)//si l'on a plus rien a copier mais que l'on vient de le faire
             {
                  copy=false;// on arrete la copie
                  *ptr_buff='\0';//on cloture la chaine contenue dans buffer
                  tab[nb_cnt]=atof(buffer);// et on la transforme en double
                  //printf(buffer);
                  //printf("\n%d\n",nb_cnt);
                  nb_cnt++;//on incrémente le nombre de nombre déjà convertis
                  ptr_buff=buffer;// on replace le pointeur au début du buffer pour réécrire dessus
             }
          }

          ptr_fich++;
          //printf("%c\n",*ptr_fich);
     }while(*(ptr_fich-1)!='\0');


     //On dispose maintenant d'un tableau en 1d contenant les x et les y
     // On va trier les x et les y en utilisant un nouveau tableau que
     // l'on va ensuite pouvoir retourner.

     return_tab=malloc(sizeof(double*)*2);
     for( i=0;i<2;i++)
     return_tab[i]=malloc(sizeof(double)*nb_cnt/2);

     for( i=0;i<nb_cnt/2;i++)//trie des x et des y et stockage de ces dernier dans un tableau
     {
     return_tab[0][i]=tab[i*2];
     return_tab[1][i]=tab[i*2+1];
     }

     //avant de terminer la fonction, on spécifie à l'utilisateur le nb de points dans
     // le tableau
   if(size!=NULL)
    *size=nb_cnt/2;

    free(tab);
    return return_tab;
}
