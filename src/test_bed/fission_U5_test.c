
#include <stdlib.h>

#include <stdio.h>
#include <stdint.h>

#include "./../fission/fiss_U5.c"
#include "./../generic/rnd_num_gen.c"





int main()
{

    int mass[119];

    for(int i=0;i<119;i++)mass[i]=0.0;

    double masse,E,Z;
    float rnd=0;
    for(int i=0;i<10000;i++)
    {
    rnd=xorshift128_uniform();
    fiss_U5_get_rnd_ff((double)rnd,&masse, &E, &Z);

    mass[(int)round(masse)-1]++;
   // printf("%e %e, i %d%\n",masse,rnd);
    }

    FILE *fich=fopen("yield_U5.txt","w");

    for(int i=0;i<119;i++)
        fprintf(fich,"%e %d \n",(double)i,mass[i]);

    fclose(fich);



return 0;
}

