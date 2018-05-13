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
    int M = 64*2;
    int N = 3*M/2;
    double H = 1.0;
    //double L = 2*H/3.0;
    double h = H/N;
    double tol = 1E-3;
    double dt = 0.1;
    double Gr = 2.0 * 1e10;
    double Pr = 2.0;

    FILE *temperature = fopen("temperature.txt", "w");
    FILE *velocity = fopen("velocity.txt", "w");
    FILE *vorticity = fopen("vorticity.txt", "w");

    if (temperature == NULL || velocity == NULL || vorticity == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

    double** norm = calloc(M, sizeof( double *));
    for (int i = 0; i < M; i++)
    {
        norm[i] = calloc(N, sizeof( double));
    }

    double** vortex = calloc(M-1, sizeof( double *));
    for (int i = 0; i<M-1; i++)
    {
            vortex[i] = calloc(N-1, sizeof( double));
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

    double T_avg = 0;

    double Re_hw;
    double Re_h;
    int timeCounter = 0;
    double r = H*dt/(sqrt(Gr)*h);

    while(T_avg < 3e-3 && Re_h < 25 && Re_hw < 40)
    {

//        printf("checkpoint ustarsolve \n");
        ustar_Solve(u, v, u_old, v_old, P, ustar, h, dt, Gr, M, N, firstStep); //Il y aura 2 if a la place d'un vu que les ifs sont incorporés dans la fonction
//        printf("checkpoint vstarsolve \n");
        vstar_Solve(u, v, u_old, v_old, P, T, vstar, h, dt, Gr, M, N, firstStep);   //Ça vaut peut etre la peine de juste mettre un if dans la main.
//        printf("checkpoint Tsolve \n");
        T_solve(T, T_old, u_old, v_old, u, v, h, dt, Pr, Gr, M, N, firstStep);

        T_avg = AverageT(T, h, M, N);
        printf("T_avg: %f\n",T_avg);
        printf("T matrix:\n");
        for (int j = 1; j < N+1; j++)
        {
            for (int k = 1; k < M+1; k++)
            {
                fprintf(temperature,"%f ", T[k][j]);
                //printf("%f ", T[k][j]);
            }
            fprintf(temperature,"\n");
            //printf("\n");
        }

        u_buffer = u_old;
        v_buffer = v_old;
        T_buffer = T_old;

        u_old = u; //Attention opérations incorrectes!!
        v_old = v;
        T_old = T;

        u = u_buffer;
        v = v_buffer;
        T = T_buffer;
//        printf("checkpoint SOR\n");
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

        Re_hw = Vortex(u, v, vortex, Gr, h, M, N, H);
        Re_h = Reynolds(u, v, norm, Gr, h, M, N, H);

        printf("time = %f \n", (timeCounter+1)*dt);
        printf("Re_hw (<40): %f, Re_h (<25): %f\n",Re_hw,Re_h);
        printf("CFL : %f\n",r*Re_h);
        timeCounter++;

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
    fclose(vorticity);
    fclose(velocity);
    for (int i = 0; i<M-1; i++)
    {
        free(vortex[i]);
    }
    for (int i = 0; i < M; i++)
    {
        free(norm[i]);
    }
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
        free(T_old[k]);
    }
    free(vortex);
    free(norm);
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
