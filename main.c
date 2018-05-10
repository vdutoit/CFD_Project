//DUTOIT Valentin - 44551400
//MECA2660_PROJECT 2018

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include"functions.h"

int main (int argc, char *argv[])
{

    //Partie initialisation
    double H = 1.0;
    double L = 2*H/3.0;
    double h = H/(3.0*1);
    int M = (int) L/h;
    int N = (int) H/h;
    double tol = 1E-3;
    double dt = 0.01;
    double t_tot = 2;
    double nt = t_tot/dt;
    double Gr = 2.0 * 1e10;
    double Pr = 2.0;
    double nu = 1.0; //A SUPPRIMER

    FILE *temperature = fopen("temperature.txt", "w");
    FILE *velocity = fopen("velocity.txt", "w");
    FILE *vorticity = fopen("vorticity.txt", "w");

    if (temperature == NULL || velocity == NULL || vorticity == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

    double** Re_h = calloc(M, sizeof( double *));
    double** norm = calloc(M, sizeof( double *));
    for (int i = 0; i < M; i++)
    {
        Re_h[i] = calloc(N, sizeof( double));
        norm[i] = calloc(N, sizeof( double));
    }
    double** Re_hw = calloc(M-1, sizeof( double *));
    double** vortex = calloc(M-1, sizeof( double *));
    for (int i = 0; i<M-1; i++)
    {
            vortex[i] = calloc(N-1, sizeof( double));
            Re_hw[i] = calloc(N-1, sizeof( double));
    }

    // Pour u et v, pas de remplissage supplémentaire nécessaire, initiallement à 0;
    double** u = calloc(M+1, sizeof( double *));
    double** u_old = calloc(M+1, sizeof( double *));
    double** ustar = calloc(M+1, sizeof( double *));
    for (int k = 0; k<M+1; k++)
    {
        u[k] = calloc(N+2,sizeof(double));
        u_old[k] = calloc(N+2,sizeof(double));
        ustar[k] = calloc(N+2,sizeof(double));
    }

    double** P = calloc(M+2, sizeof( double *));
    double** T = calloc(M+2, sizeof( double *));
    double** T_old = calloc(M+2, sizeof( double *));
    double** v = calloc(M+2, sizeof( double *));
    double** v_old = calloc(M+2, sizeof( double *));
    double** vstar = calloc(M+2, sizeof( double *));
    double** phi = calloc(M+2, sizeof( double *));
    for (int k = 0; k<M+2; k++)
    {
        phi[k] = calloc(N+2,sizeof(double));
        v[k] = calloc(N+1,sizeof(double));
        v_old[k] = calloc(N+1,sizeof(double));
        vstar[k] = calloc(N+1,sizeof(double));
        P[k] = calloc(N+2,sizeof(double));
        T[k] = calloc(N+2,sizeof(double));
        T_old[k] = calloc(N+2,sizeof(double));
    }

    int firstStep = 1;

    double** u_buffer;
    double** v_buffer;
    double** T_buffer;

    for (int i = 0; i < nt; i++)
    {
        printf("time = %f \n", i*dt);
        printf("checkpoint ustarsolve \n");
        ustar_Solve(u, v, u_old, v_old, P, ustar, h, dt, Gr, M, N, firstStep); //Il y aura 2 if a la place d'un vu que les ifs sont incorporés dans la fonction
        printf("checkpoint vstarsolve \n");
        vstar_Solve(u, v, u_old, v_old, P, T, vstar, h, dt, Gr, M, N, firstStep);   //Ça vaut peut etre la peine de juste mettre un if dans la main.
        printf("checkpoint Tsolve \n");
        T_solve(T, T_old, u_old, v_old, u, v, h, dt, Pr, Gr, M, N, firstStep);

        u_buffer = u_old;
        v_buffer = v_old;
        T_buffer = T_old;

        u_old = u; //Attention opérations incorrectes!!
        v_old = v;
        T_old = T;

        u = u_buffer;
        v = v_buffer;
        T = T_buffer;

        SOR(phi, ustar, vstar, tol, 1.97, h, dt, M, N);
        u_Solve(ustar, phi, u, dt, h, M, N);
        v_Solve(vstar, phi, v, dt, h, M, N);
        P_solve(P, phi, M, N);

        for(int n = 0; n < N+2; n++)
        {
            for(int m = 0; m < M+2; m++)
            {
                phi[m][n] = 0;
            }
        }

        Vortex(u, v, Re_hw, vortex, nu, h, M, N);
        Reynolds(u, v, Re_h, norm, nu, h, M, N);
        for (int j = 0; j < N+2; j++)
        {
            for (int k = 0; k < M+2; k++)
            {
                fprintf(temperature,"%f ", T[k][j]);
            }
            fprintf(temperature,"\n");
        }
        for (int i = 0; i<M-1; i++)
        {
            for (int j = 0; j<N-1; j++)
            {
                fprintf(vorticity,"%f ", vortex[i][j]);
            }
            fprintf(vorticity,"\n");
        }
        for (int i = 0; i < M; i++)
        {
            for (int j = 0; j < N; j++)
            {
                fprintf(velocity,"%f ", norm[i][j]);
            }
            fprintf(velocity,"\n");
        }

        firstStep = 0;
    }

    fclose(temperature);
    for (int k = 0; k<M+1; k++)
    {
        free(u[k]);
        free(u_old[k]);
        free(ustar[k]);
    }

    for (int k = 0; k<M+2; k++)
    {
        free(v[k]);
        free(v_old[k]);
        free(vstar[k]);
        free(P[k]);
        free(T[k]);
        free(phi[k]);
    }
    free(u);
    free(u_old);
    free(ustar);
    free(v);
    free(v_old);
    free(vstar);
    free(P);
    free(T);
    free(phi);
    free(T_old);


}
