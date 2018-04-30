////DUTOIT Valentin - 44551400
////MECA2660_PROJECT_2018
#include<stdio.h>
#include<math.h>
#include"thomas.h"
#include"functions.h"

void ustar_Solve(double** u, double** v, double** adv, double** P, double** sol, double h, double dt, double nu, int M, int N)
{
    //u est M+1xN+2, a est M-1xN, P est M+2xN+2, sol est M+1xN+2, adv est M-1xN

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

    dPdx_fun(P,dPdx,h,N,M);
    d2udx2_fun(u,d2udx2,h,N,M);
    d2udy2_fun(u,d2udy2,h,N,M);

    for (int i = 0; i<M-1; i++)
    {
        for (int j = 0; j<N; j++)
        {
            sol[i+1][j+1] = u[i+1][j+1] + dt * (-1*adv[i][j] - dPdx[i][j]  + nu * (d2udx2[i][j] + d2udy2[i][j]));
        }
    }

    //imposer conditions aux limites pour ustar -> 0 sur les paroies laterales (rien changer), u0 en y=0 et u[i][N+1] = u[i][N] sauf en i=0 et i=M

    for (int i = 1; i<M; i++)
    {
        sol[i][N+1] = sol[i][N]; //free surface BC
        sol[i][0] = -0.2*(sol[i][3]-5*sol[i][2]+15*sol[i][1]); //no slip condition at bottom
    }

    free(dPdx);
    free(d2udx2);
    free(d2udy2);
}

void vstar_Solve(double** u, double** v, double** adv, double** P, double** T, double** sol, double h, double dt, double T0, double nu, double beta, int M, int N) //AJOUTER TEMPERATURE
{
    //v est M+2xN+1, b est MxN-1, P est M+2xN+2, sol est M+2xN+1, adv est MxN-1

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

    dPdy_fun(P,dPdy,h,N,M);
    d2vdx2_fun(v,d2vdx2,h,N,M);
    d2vdy2_fun(v,d2vdy2,h,N,M);

    for (int i = 0; i<M; i++)
    {
        for (int j = 0; j<N-1; j++)
        {
            sol[i+1][j+1] = v[i+1][j+1] + dt * (-1*adv[i][j] - dPdy[i][j] + 9.81 * beta*(T[i+1][j+1] - T0) + nu * (d2vdx2[i][j] + d2vdy2[i][j]));
        }
    }

    //imposer conditions aux limites pour ustar -> 0 sur les paroies laterales (rien changer), u0 en y=0 et u[i][N+1] = u[i][N] sauf en i=0 et i=M

    for (int j = 1; j<N; j++)
    {
        sol[0][j] = -0.2*(sol[3][j]-5*sol[2][j]+15*sol[1][j]); //no slip condition at left wall
        sol[M+1][j] = -0.2*(sol[M+1-3][j]-5*sol[M+1-2][j]+15*sol[M+1-1][j]); //no slip condition at right wall
    }

    free(dPdy);
    free(d2vdx2);
    free(d2vdy2);
}

void T_solve(double** T, double** H, double** sol, double h, double dt, double q_w, double T_inf, double h_barre, double k, double alpha, int M, int N)
{
    //T est M+2xN+2, H est MxN, deja assemble avec AB2 ou euler, sol est M+2xN+2

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

    double* dTdy_e  = calloc(M, sizeof(double));

    d2Tdx2_fun(T, d2Tdx2, h, N, M);
    d2Tdy2_fun(T, d2Tdy2, h, N, M);

    for (int i = 0; i<M; i++)
    {
        for (int j = 0; j<N; j++)
        {
            sol[i+1][j+1] = T[i+1][j+1] + dt * ( -1 * H[i][j] + alpha * (d2Tdx2[i][j] + d2Tdy2[i][j]));
        }
    }

    for (int i = 0; i<N; i++)
    {
        sol[0][i+1] = sol[1][i+1];
        sol[M+1][i+1] = sol[M][i+1];
    }

    double q_e = 0;
    for (int i=1;1<M+1;i++)
    {
        q_e = h_barre*((T[i][N+1]-T[i][N])/2 - T_inf);
        dTdy_e[i-1] = -q_e/k;
    }

    double dTdy_w = -q_w/k;

    for (int i = 0; i<M; i++)
    {
        sol[i+1][0] = sol[i+1][1] - h*dTdy_w;
        sol[i+1][N+1] = sol[i+1][N] + h*dTdy_e[i];
    }

    free(d2Tdy2);
    free(d2Tdx2);

}
