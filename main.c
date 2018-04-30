//DUTOIT Valentin - 44551400
//MECA2660_PROJECT 2018

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"functions.h"

int main (int argc, char *argv[])
{
    int M = 2;
    int N = 3;
    double h = 1.0;
    double** a = calloc(M-1, sizeof(double));
        for (int k = 0; k<M-1; k++)
        {
            a[k] = calloc(N,sizeof(double));
        }
    double** P = calloc(M+2, sizeof( double *));
        for (int k = 0; k<M+2; k++)
        {
            P[k] = calloc(N+2,sizeof(double));
        }
    double** u = calloc(M+1, sizeof( double *));
        for (int k = 0; k<M+1; k++)
        {
            u[k] = calloc(N+2,sizeof(double));
        }
    double** v = calloc(M+2, sizeof( double *));
        for (int k = 0; k<M+2; k++)
        {
            v[k] = calloc(N+1,sizeof(double));
        }

//    for (int i=0; i<M+1; i++)
//    {
//            for (int j=0; j<N+2; j++)
//            {
//                u[i][j] = 0;
//            }
//    }

//    for (int i=0; i<M+2; i++)
//    {
//            for (int j=0; j<N+1; j++)
//            {
//                v[i][j] = 2*i+j;
//            }
//    }

//    for( int i = 0 ; i<M;i++)
//    {
//        for (int j = 0; j <N;j++)
//        {
//            printf("%f \n",u[i][j]);
//        }
//    }

    u[0][1] = 4;
    u[1][1] = 1;
    u[2][1] = 0;

    AdvectiveX_fun(u,v,a,h,N,M);

    printf("a00 = %f",a[0][0]);

    free(a);
    free(u);
    free(v);
}
