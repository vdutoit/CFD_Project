//DUTOIT Valentin - 44551400
//MECA2660_PROJECT 2018

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"functions.h"

int main (int argc, char *argv[])
{

    //Partie initialisation

    double h = 0.5;
    double L = 1;
    double H = 1.5*L;
    int M = (int) L/h;
    int N = (int) H/h;
    double rho_0 = 1E3;   //densité de l'eau
    double p_ref = 1E5; //[Pa]
    double p_init = 101325;// Un peu au pif
    double g = 9.81;       //[m/s^2]
    double T_0 = 298.15;   //[K]
    double T_inf = 273.15; //[K]
    double dt = 0.1;
    double t_tot = 10;
    double nt = t_tot/dt;
    double k = 1;
    double l_0 = (1E-3)*H;
    double h_barre = k/l_0;
    double deltaT = (T_0 - T_inf)/(5E3);
    double nu = 1E-6;
    double alpha = nu/2;
    double U = sqrt(2E10)*nu/H;
    double beta = pow(U,2)/(g*deltaT*H);
    double q_w = k*deltaT/H;
    double tol = 1E-3;

    FILE *temperature = fopen("temperature.txt", "w");

    if (temperature == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

    // On considère initiallement une matrice petit p identique partout, ici représentée par un scalaire
    for(int i = 0; i<M+2; i++)
    {
        for(int j = 0; j<N+2; j++)
        {
            P[i][j] = ((p_init - p_ref) + rho_0*g*(H - j*h))/rho_0;
            T[i][j] = T_0;
        }
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
    double** phi = calloc(M+2, sizeof( double *));
    double** v = calloc(M+2, sizeof( double *));
    double** v_old = calloc(M+2, sizeof( double *));
    double** vstar = calloc(M+2, sizeof( double *));
    for (int k = 0; k<M+2; k++)
    {
        v[k] = calloc(N+1,sizeof(double));
        v_old[k] = calloc(N+1,sizeof(double));
        vstar[k] = calloc(N+1,sizeof(double));
        P[k] = calloc(N+2,sizeof(double));
        T[k] = calloc(N+2,sizeof(double));
        phi[k] = calloc(N+2,sizeof(double));
    }

    int firstStep = 1;

    double** u_buffer;
    double** v_buffer;

    for (int i = 0; i < nt; i++)
    {
        for (int j = 0; j < N+2; j++)
        {
            for (int k = 0; k < M+2; k++)
            {
                fprintf(temperature,"%f ", (T[k][j]-T_0)/T_inf);
            }
        }

        ustar_Solve(u, v, u_old, v_old, P, ustar, h, dt, nu, M, N, firstStep); //Il y aura 2 if a la place d'un vu que les ifs sont incorporés dans la fonction
        vstar_Solve(u, v, u_old, v_old, P, vstar, h, dt, nu, M, N, firstStep); //Ça vaut peut etre la peine de juste mettre un if dans la main.

        T_solve(T, u_old, v_old, u, v, h, dt, q_w, T_inf, h_barre, k, alpha, M, N, firstStep)

        u_buffer = u_old;
        v_buffer = v_old;

        u_old = u; //Attention opérations incorrectes!!
        v_old = v;

        u = u_buffer;
        v = v_buffer;

        SOR(phi, ustar, vstar, tol, alpha, H, U, L, h, dt, M, N);

        u_Solve(ustar, phi, u, dt, h, M, N);
        v_Solve(vstar, phi, v, dt, h, M, N);

        P_solve(P, phi, M, N);

        firstStep = 0;


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
    }

}
