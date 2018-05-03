//DUTOIT Valentin - 44551400
//MECA2660_PROJECT_2018
#include<stdio.h>
#include<math.h>
#include"functions.h"

//TRIER FONCTIONS DANS PLUSIEURS DOSSIERS + FAIRE HEADER

void dPdx_fun(double** P, double** sol, double h, int M, int N)
{
    //P est M+2xN+2, sol est de dim M-1xN (positionne sur les u) avec M et N le nombre de points dans le domaine frontiere comprise
    for (int j = 0; j<N; j++)
    {
        for (int i = 0; i<M-1; i++)
        {
            sol[i][j] = (P[i+2][j+1]-P[i+1][j+1])/h;
        }
    }
}

void dPdy_fun(double** P, double** sol, double h, int M, int N)
{
    //P est M+2xN+2, sol est de dim MxN-1 (positionne sur les v) avec M et N le nombre de points dans le domaine frontiere comprise
    for (int i = 0; i<M-1; i++)
    {
        for (int j = 0; j<N; j++)
        {
            sol[i][j] = (P[i+1][j+2]-P[i+1][j+1])/h;
        }
    }
}

void d2udx2_fun(double** u, double** sol, double h, int M, int N)
{
    //u est M+1xN+2, sol est de dim M-1xN avec M et N le nombre de points dans le domaine frontiere comprise
    for (int j = 0; j<N; j++)
    {
        for (int i = 0; i<M-1; i++)
        {
            sol[i][j] = (u[i+2][j+1]-2*u[i+1][j+1]+u[i][j+1])/pow(h,2);
        }
    }
}

void d2udy2_fun(double** u, double** sol, double h, int M, int N)
{
    //u est M+1xN+2, sol est de dim M-1xN avec M et N le nombre de points dans le domaine frontiere comprise
    for (int i = 0; i<M-1; i++)
    {
        for (int j = 0; j<N; j++)
        {
            sol[i][j] = (u[i+1][j+2]-2*u[i+1][j+1]+u[i+1][j])/pow(h,2);
        }
    }
}

void d2vdx2_fun(double** v, double** sol, double h, int M, int N) //IDEM QUE POUR U A SUPPRIMER
{
    //v est M+2xN+1, sol est de dim MxN-1 avec M et N le nombre de points dans le domaine frontiere comprise
    for (int j = 0; j<N-1; j++)
    {
        for (int i = 0; i<M; i++)
        {
            sol[i][j] = (v[i+2][j+1]-2*v[i+1][j+1]+v[i][j+1])/pow(h,2);
        }
    }
}

void d2vdy2_fun(double** v, double** sol, double h, int M, int N) //IDEM QUE POUR U A SUPPRIMER
{
    //v est M+2xN+1, sol est de dim MxN-1 avec M et N le nombre de points dans le domaine frontiere comprise
    for (int i = 0; i<M; i++)
    {
        for (int j = 0; j<N-1; j++)
        {
            sol[i][j] = (v[i+1][j+2]-2*v[i+1][j+1]+v[i+1][j])/pow(h,2);
        }
    }
}

void d2Tdx2_fun(double** T, double** sol, double h, int M, int N)
{   //T est M+2xN+2, sol est MxN AUSSI VALABLE POUR P ET PHI
    for (int i = 0; i<M; i++)
    {
        for (int j = 0; j<N; j++)
            sol[i][j] = (T[i+2][j+1]-2*T[i+1][j+1]+T[i][j+1])/pow(h,2);
    }
}

void d2Tdy2_fun(double** T, double** sol, double h, int M, int N)
{   //T est M+2xN+2, sol est MxN AUSSI VALABLE POUR P ET PHI
    for (int i = 0; i<M; i++)
    {
        for (int j = 0; j<N; j++)
            sol[i][j] = (T[i+1][j+2]-2*T[i+1][j+1]+T[i+1][j])/pow(h,2);
    }
}

void AdvectiveX_fun(double** u, double** v, double** a, double h, int M, int N)
{
    //u est M+1xN+2 , v est M+2xN+1, a est M-1xN avec M et N le nombre de points dans le domaine frontiere comprise

    double** Uij = calloc(M, sizeof( double *));
    for (int k = 0; k<M; k++)
    {
        Uij[k] = calloc(N,sizeof(double));
    }
    double** Vij = calloc(M-1, sizeof( double *));
    for (int k = 0; k<M-1; k++)
    {
        Vij[k] = calloc(N+1,sizeof(double));
    }

    for (int i=0; i<M; i++)
    {
        for (int j=0; j<N; j++)
        {
            Uij[i][j] = 0.5 * (u[i+1][j+1] + u[i][j+1]);//U[0][0] est place en (1,1)
            printf("U[%d][%d] = %f \n",i,j,Uij[i][j]);
        }
    }

    for (int i=0; i<M-1; i++)
    {
        for (int j=0; j<N+1; j++)
        {
            Vij[i][j] = 0.5 * (v[i+2][j] + v[i+1][j]); // V[0][0] est place en (1.5,0.5)

        }
    }

    for (int i=0; i<M-1; i++)
    {
        for (int j=0; j<N; j++)
        {
            //a[0][0] est place en (1.5,1)
            a[i][j] = 0.5 *(Uij[i+1][j] * (u[i+2][j+1]-u[i+1][j+1])/h + Uij[i][j] * (u[i+1][j+1]-u[i][j+1])/h) + 0.5 *(Vij[i][j+1] * (u[i+1][j+2]-u[i+1][j+1])/h + Vij[i][j] * (u[i+1][j+1]-u[i+1][j])/h);
        }
    }

    free(Uij);
    free(Vij);

}

void AdvectiveY_fun(double** u, double** v, double** b, double h, int M, int N)
{
    //u est M+1xN+2 , v est M+2xN+1, b est MxN-1 avec M et N le nombre de points dans le domaine frontiere comprise

    double** Uij = calloc(M+1, sizeof( double *));
    for (int k = 0; k<M+1; k++)
    {
        Uij[k] = calloc(N-1,sizeof(double));
    }
    double** Vij = calloc(M, sizeof( double *));
    for (int k = 0; k<M; k++)
    {
        Vij[k] = calloc(N,sizeof(double));
    }

    for (int i=0; i<M; i++)
    {
        for (int j=0; j<N; j++)
        {
            Uij[i][j] = 0.5 * (u[i][j+2] + u[i][j+1]);//U[0][0] est place en (0.5,1.5)
//                printf("U[%d][%d] = %f \n",i,j,Uij[i][j]);
        }
    }

    for (int i=0; i<M; i++)
    {
        for (int j=0; j<N+1; j++)
        {
            Vij[i][j] = 0.5 * (v[i+1][j+1] + v[i+1][j]); // V[0][0] est place en (1,1)

        }
    }

    for (int i=0; i<M; i++)
    {
        for (int j=0; j<N-1; j++)
        {
            //b[0][0] est place en (1,1.5)
            b[i][j] = 0.5 *(Vij[i][j+1] * (v[i+1][j+2]-v[i+1][j+1])/h + Vij[i][j] * (v[i+1][j+1]-v[i+1][j])/h)
                      + 0.5 *(Uij[i+1][j] * (v[i+2][j+1]-v[i+1][j+1])/h + Uij[i][j] * (v[i+1][j+1]-v[i][j+1])/h);
        }
    }

    free(Uij);
    free(Vij);

}

AdvectiveT_fun(double** T, double** u, double** v, double** H, double h, int M, int N)
{
    //T est M+2xN+2, u est M+1xN+2, v est M+2xN+1, H est MxN

    for (int i=0; i<M; i++)
    {
        for (int j=0; j<N; j++ )
        {
            H[i][j] = 0.5*(u[i][j+1]*(T[i+1][j+1]-T[i][j+1])/h + u[i+1][j+1]*(T[i+2][j+1]-T[i+1][j+1])/h)
                      + 0.5*(v[i+1][j]*(T[i+1][j+1]-T[i+1][j])/h + v[i+1][j+1]*(T[i+1][j+2]-T[i+1][j+1])/h);
        }
    }
}

void AB2X_fun(double** a_old, double** a_now, double** sol, int M, int N)
{
    //a_old, a_now et sol sont M-1xN avec M et N le nombre de points dans le domaine frontiere comprise

    for (int i = 0; i<M-1; i++)
    {
        for(int j = 0; j<N; j++)
        {
            sol[i][j] = 1.5*a_now[i][j] - 0.5*a_old[i][j];
        }
    }
}

void AB2Y_fun(double** b_old, double** b_now, double** sol, int M, int N)
{
    //b_old, b_now et sol sont MxN-1 avec M et N le nombre de points dans le domaine frontiere comprise

    for (int i = 0; i<M; i++)
    {
        for(int j = 0; j<N-1; j++)
        {
            sol[i][j] = 1.5*b_now[i][j] - 0.5*b_old[i][j];
        }
    }
}

void AB2T_fun(double** H_old, double** H_now, double** sol, int M, int N)
{
    //H_old, H_now et sol sont MxN

    for (int i = 0; i<M; i++)
    {
        for(int j = 0; j<N; j++)
        {
            sol[i][j] = 1.5*H_now[i][j] - 0.5*H_old[i][j];
        }
    }
}


void dudx_fun(double** u, double** sol, double h, int M, int N)
{
    //derivee centree sur les point i,j comme le champ P
    //u est M+1xN+2, sol est MxN

    for (int i = 0; i<M; i++)
    {
        for (int j = 0; j<N; j++)
        {
            sol[i][j] = (u[i+1][j+1]-u[i][j+1])/h;
        }
    }
}

void dvdy_fun(double** v, double** sol, double h, int M, int N)
{
    //derivee centree sur les point i,j comme le champ P
    //v est M+2xN+1, sol est MxN

    for (int i = 0; i<M; i++)
    {
        for (int j = 0; j<N; j++)
        {
            sol[i][j] = (v[i+1][j+1]-v[i+1][j])/h;
        }
    }
}

void dvdx_fun(double** v, double** sol, double h, int M, int N)
{
    //derivee centree sur les points sans symbole
    //v est M+2xN+1, sol est M-1xN-1

    for (int i = 0; i<M-1; i++)
    {
        for (int j = 0; j<N-1; j++)
        {
            sol[i][j] = (v[i+2][j+1]-v[i+1][j+1])/h;
        }
    }
}

void dudy_fun(double** u, double** sol, double h, int M, int N)
{
    //derivee centree sur les points sans symbole
    //v est M+2xN+1, sol est M-1xN-1

    for (int i = 0; i<M-1; i++)
    {
        for (int j = 0; j<N-1; j++)
        {
            sol[i][j] = (u[i+1][j+2]-u[i+1][j+1])/h;
        }
    }
}

void getNorm(double** u, double** v, double** sol, double h, int M, int N)
{
    double** Uij = calloc(M+1, sizeof( double *));
    for (int k = 0; k<M; k++)
    {
        Uij[k] = calloc(N-1,sizeof(double));
    }
}
