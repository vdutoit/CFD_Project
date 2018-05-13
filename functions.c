//DUTOIT Valentin - 44551400
//MECA2660_PROJECT_2018
#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>
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
    for (int i = 0; i<M; i++)
    {
        for (int j = 0; j<N-1; j++)
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
//    printf("checkpoint7\n");
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
        }
    }
//    printf("checkpoint8\n");
//    printf("M,N: %d, %d\n",M,N);
//    int iter = 1;
    for (int i=0; i<M-1; i++)
    {
        for (int j=0; j<N+1; j++)
        {
//            printf("check\n");
//            printf("v[i+2][j]: %f\n", v[i+2][j]);
            Vij[i][j] = 0.5 * (v[i+2][j] + v[i+1][j]); // V[0][0] est place en (1.5,0.5)
//            printf("iter: %d\n", iter);
//            iter++;

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
    for (int k = 0; k<M; k++)
    {
        free(Uij[k]);
    }
    for (int k = 0; k<M-1; k++)
    {
        free(Vij[k]);
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

    for (int i=0; i<M+1; i++)
    {
        for (int j=0; j<N-1; j++)
        {
            Uij[i][j] = 0.5 * (u[i][j+2] + u[i][j+1]);//U[0][0] est place en (0.5,1.5)
        }
    }

    for (int i=0; i<M; i++)
    {
        for (int j=0; j<N; j++)
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
    for (int k = 0; k<M+1; k++)
    {
        free(Uij[k]);
    }
    for (int k = 0; k<M; k++)
    {
        free(Vij[k]);
    }
    free(Uij);
    free(Vij);

}

void AdvectiveT_fun(double** T, double** u, double** v, double** H, double h, int M, int N)
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

double Vortex(double** u, double** v, double** vortex, double Gr, double h, int M, int N, double H)
{
    double dvdx;
    double dudy;
    double Re_hw;
    double Re_hw_max = 0;
    for (int i = 0; i<M-1; i++)
    {
        for (int j = 0; j<N-1; j++)
        {
            dvdx = (v[i+2][j+1]-v[i+1][j+1])/h;
            dudy = (u[i+1][j+2]-u[i+1][j+1])/h;
            vortex[i][j] = dvdx - dudy;
            Re_hw = sqrt(Gr)*fabs(dvdx-dudy)*pow(h,2)/pow(H,2);
            if (Re_hw > Re_hw_max)
            {
                Re_hw_max = Re_hw;
            }
        }
    }
    return Re_hw_max;
}

double Reynolds(double** u, double** v, double** norm, double Gr, double h, int M, int N, double H)
{
    double u_avg;
    double v_avg;
    double Re_h;
    double Re_h_max = 0;
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            u_avg = (u[i][j+1]+u[i+1][j+1])/2;
            v_avg = (v[i+1][j]+v[i+1][j+1])/2;
            Re_h = sqrt(Gr)*(fabs(u_avg)+fabs(v_avg))*h/H;
            if (Re_h > Re_h_max)
            {
                Re_h_max = Re_h;
            }
            norm[i][j] = sqrt(u_avg*u_avg + v_avg*v_avg);
        }
    }
    return Re_h_max;
}
double AverageT(double** T, double h, int M, int N)
{
    double T_avg = 0;
    for (int i = 0; i<M; i++)
    {
        for (int j = 0; j<N; j++)
        {
            T_avg += T[i+1][j+1];
        }
    }
    T_avg = T_avg/(M*N);
    return T_avg;
}

void T_RMS(double** T, double T_rms, double T_avg, double h, int M, int N)
{
    for (int i = 0; i<M; i++)
    {
        for (int j = 0; j<N; j++)
        {
            T_rms += pow(T[i+1][j+1]-T_avg,2);
        }
    }
    T_rms = sqrt(T_rms/(M*N));
}
