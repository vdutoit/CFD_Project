//DUTOIT Valentin - 44551400
//MECA2660_PROJECT_2018
#include<stdio.h>
#include<math.h>
#include"thomas.h"
#include"functions.h"

void ustar_Solve(double** u, double** a_old, double** a_now, double** P, double** sol, double h, double dt, double nu, int M, int N) //AJOUTER TEMPERATURE
{
    //u est M+1xN+2, a est M-1xN, P est M+2xN+2, sol est M+1xN+2
    double** adv  = calloc(M-1, sizeof( double *));
    for (int k = 0; k<M-1; k++)
    {
        adv[k] = calloc(N,sizeof(double));
    }

    double** dPdx  = calloc(M-1, sizeof( double *));
    for (int k = 0; k<M-1; k++)
    {
        dPdx[k] = calloc(N,sizeof(double));
    }

    double** d2udx2  = calloc(M-1, sizeof( double *));
    for (int k = 0; k<M-1; k++)
    {
        d2udx2[k] = calloc(N,sizeof(double));
    }

    double** d2udy2  = calloc(M-1, sizeof( double *));
    for (int k = 0; k<M-1; k++)
    {
        d2udy2[k] = calloc(N,sizeof(double));
    }

    AB2X_fun(a_old,a_now,adv,M,N);
    dPdx_fun(P,dPdx,h,N,M);
    d2udx2_fun(u,d2udx2,h,N,M);
    d2udy2_fun(u,d2udy2,h,N,M);

    for (int i = 0; i<M-1; i++)
    {
        for (int j = 0; j<N; j++)
        {
            sol[i+1][j+1] = u[i+1][j+1] + dt * (-1*adv[i][j] - dPdx[i][j] + nu * (d2udx2[i][j] + d2udy2[i][j])); // + terme temperature
        }
    }

    //imposer conditions aux limites pour ustar -> 0 sur les paroies laterales (rien changer), u0 en y=0 et u[i][N+1] = u[i][N] sauf en i=0 et i=M

    for (int i = 1; i<M; i++)
    {
        sol[i][N+1] = sol[i][N]; //free surface BC
        sol[i][0] = -0.2*(sol[i][3]-5*sol[i][2]+15*sol[i][1]); //no slip condition at bottom
    }

    free(adv);
    free(dPdx);
    free(d2udx2);
    free(d2udy2);
}

void vstar_Solve(double** v, double** b_old, double** b_now, double** P, double** sol, double h, double dt, double nu, int M, int N) //AJOUTER TEMPERATURE
{
    //v est M+2xN+1, b est MxN-1, P est M+2xN+2, sol est M+2xN+1
    double** adv  = calloc(M, sizeof( double *));
    for (int k = 0; k<M; k++)
    {
        adv[k] = calloc(N-1,sizeof(double));
    }

    double** dPdy  = calloc(M, sizeof( double *));
    for (int k = 0; k<M; k++)
    {
        dPdy[k] = calloc(N-1,sizeof(double));
    }

    double** d2vdx2  = calloc(M, sizeof( double *));
    for (int k = 0; k<M; k++)
    {
        d2vdx2[k] = calloc(N-1,sizeof(double));
    }

    double** d2vdy2  = calloc(M, sizeof( double *));
    for (int k = 0; k<M; k++)
    {
        d2vdy2[k] = calloc(N-1,sizeof(double));
    }

    AB2Y_fun(b_old,b_now,adv,M,N);
    dPdy_fun(P,dPdy,h,N,M);
    d2vdx2_fun(v,d2vdx2,h,N,M);
    d2vdy2_fun(v,d2vdy2,h,N,M);

    for (int i = 0; i<M; i++)
    {
        for (int j = 0; j<N-1; j++)
        {
            sol[i+1][j+1] = v[i+1][j+1] + dt * (-1*adv[i][j] - dPdy[i][j] + nu * (d2vdx2[i][j] + d2vdy2[i][j])); // + terme temperature
        }
    }

    //imposer conditions aux limites pour ustar -> 0 sur les paroies laterales (rien changer), u0 en y=0 et u[i][N+1] = u[i][N] sauf en i=0 et i=M

    for (int j = 1; j<N; j++)
    {
        sol[0][j] = -0.2*(sol[3][j]-5*sol[2][j]+15*sol[1][j]); //no slip condition at left wall
        sol[M+1][j] = -0.2*(sol[M+1-3][j]-5*sol[M+1-2][j]+15*sol[M+1-1][j]); //no slip condition at right wall
    }

    free(adv);
    free(dPdy);
    free(d2vdx2);
    free(d2vdy2);
}

void T_solve(double** T,  double** u, double** v, double** sol, double h, double dt, double alpha, int M, int N)
{
    double** d2Tdx2  = calloc(M, sizeof( double *));
    for (int k = 0; k<M; k++)
    {
        d2Tdx2[k] = calloc(N,sizeof(double));
    }

    double** d2Tdy2  = calloc(M, sizeof( double *));
    for (int k = 0; k<M; k++)
    {
        d2Tdy2[k] = calloc(N,sizeof(double));
    }
    d2Tdx2_fun(T, d2Tdx2, h, N, M)
    d2Tdy2_fun(T, d2Tdy2, h, N, M)

    for (int i = 0; i<M; i++)
    {
        for (int j = 0; j<N; j++)
        {
            sol[i+1][j+1] = T[i+1][j+1] + dt * ( - + alpha * (d2Tdx2[i][j] + d2Tdy2[i][j])); // + terme temperature
        }
    }
}
