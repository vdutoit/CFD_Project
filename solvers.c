////DUTOIT Valentin - 44551400
////MECA2660_PROJECT_2018
#include<stdio.h>
#include<math.h>
#include"thomas.h"
#include"functions.h"

void ustar_Solve(double** u, double** v, double** u_old, double** v_old, double** P, double** sol, double h, double dt, double nu, int M, int N, int firstStep)
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

    double** H_now  = calloc(M-1, sizeof( double *));
    for (int k = 0; k<M-1; k++)
    {
        H_now[k] = calloc(N,sizeof(double));
    }

    dPdx_fun(P,dPdx,h,M,N);
    d2udx2_fun(u,d2udx2,h,M,N);
    d2udy2_fun(u,d2udy2,h,M,N);
    AdvectiveX_fun(u,v,H_now,h,M,N);
    if (firstStep == 0)
    {
        double** H_old  = calloc(M-1, sizeof( double *));
        for (int k = 0; k<M-1; k++)
        {
            H_old[k] = calloc(N,sizeof(double));
        }
        AdvectiveX_fun(u_old,v_old,H_old,h,M,N);
        for (int i = 0; i<M-1; i++)
        {
            for (int j = 0; j<N; j++)
            {
                sol[i+1][j+1] = u[i+1][j+1] + dt * (-0.5*(3*H_now[i][j]-H_old[i][j]) - dPdx[i][j]  + nu * (d2udx2[i][j] + d2udy2[i][j]));
            }
        }
        free(H_old);
    }
    else
    {
        for (int i = 0; i<M-1; i++)
        {
            for (int j = 0; j<N; j++)
            {
                sol[i+1][j+1] = u[i+1][j+1] + dt * (-1*H_now[i][j] - dPdx[i][j]  + nu * (d2udx2[i][j] + d2udy2[i][j]));
            }
        }
    }

    //imposer conditions aux limites pour ustar -> 0 sur les paroies laterales (rien changer), u0 en y=0 et u[i][N+1] = u[i][N] sauf en i=0 et i=M

    for (int i = 1; i<M; i++)
    {
        sol[i][N+1] = sol[i][N]; //free surface BC
        sol[i][0] = -0.2*(sol[i][3]-5*sol[i][2]+15*sol[i][1]); //no slip condition at bottom
    }

    for (int k = 0; k<M-1; k++)
    {
         free(dPdx[k]);
         free(d2udx2[k]);
         free(d2udy2[k]);
         free(H_now[k]);
    }
    free(dPdx);
    free(d2udx2);
    free(d2udy2);
    free(H_now);
}

void vstar_Solve(double** u, double** v, double** u_old, double** v_old, double** P, double** T, double** sol, double h, double dt, double T0, double nu, double beta, int M, int N, int firstStep) //AJOUTER TEMPERATURE
{
    //v est M+2xN+1, b est MxN-1, P est M+2xN+2, sol est M+2xN+1, adv est MxN-1

    double** dPdy  = calloc(M, sizeof( double *));
    double** d2vdx2  = calloc(M, sizeof( double *));
    double** d2vdy2  = calloc(M, sizeof( double *));
    double** H_now  = calloc(M, sizeof( double *));
    for (int k = 0; k<M; k++)

    {
        dPdy[k] = calloc(N-1,sizeof(double));
        d2vdx2[k] = calloc(N-1,sizeof(double));
        d2vdy2[k] = calloc(N-1,sizeof(double));
        H_now[k] = calloc(N-1,sizeof(double));
    }

    dPdy_fun(P,dPdy,h,M,N);
    d2vdx2_fun(v,d2vdx2,h,M,N);
    d2vdy2_fun(v,d2vdy2,h,M,N);
    AdvectiveY_fun(u,v,H_now,h,M,N);

    if (firstStep == 0)
    {
        double** H_old  = calloc(M, sizeof( double *));
        for (int k = 0; k<M; k++)
        {
            H_old[k] = calloc(N-1,sizeof(double));
        }
        AdvectiveY_fun(u_old,v_old,H_old,h,M,N);
        for (int i = 0; i<M; i++)
    {
        for (int j = 0; j<N-1; j++)
        {
            sol[i+1][j+1] = v[i+1][j+1] + dt * (-0.5*(3*H_now[i][j]-H_old[i][j]) - dPdy[i][j] + 9.81 * beta*(T[i+1][j+1] - T0) + nu * (d2vdx2[i][j] + d2vdy2[i][j]));
        }
    }
        free(H_old);
    }
    else
    {
        for (int i = 0; i<M; i++)
        {
            for (int j = 0; j<N-1; j++)
            {
                sol[i+1][j+1] = v[i+1][j+1] + dt * (-1*H_now[i][j] - dPdy[i][j] + 9.81 * beta*(T[i+1][j+1] - T0) + nu * (d2vdx2[i][j] + d2vdy2[i][j]));
                //printf("%f\n", dPdy[i][j]);
            }
        }
    }

    //imposer conditions aux limites pour ustar -> 0 sur les paroies laterales (rien changer), u0 en y=0 et u[i][N+1] = u[i][N] sauf en i=0 et i=M

    for (int j = 1; j<N; j++)
    {
        sol[0][j] = -0.2*(sol[3][j]-5*sol[2][j]+15*sol[1][j]); //no slip condition at left wall
        sol[M+1][j] = -0.2*(sol[M+1-3][j]-5*sol[M+1-2][j]+15*sol[M+1-1][j]); //no slip condition at right wall
    }

    if (firstStep == 0)
    {
        double** H_old  = calloc(M, sizeof( double *));
        for (int k = 0; k<M; k++)
        {
            H_old[k] = calloc(N-1,sizeof(double));
        }
        AdvectiveY_fun(u_old,v_old,H_old,h,M,N);
        for (int i = 0; i<M; i++)
    {
        for (int j = 0; j<N-1; j++)
        {
            sol[i+1][j+1] = v[i+1][j+1] + dt * (-0.5*(3*H_now[i][j]-H_old[i][j]) - dPdy[i][j] + 9.81 * beta*(T[i+1][j+1] - T0) + nu * (d2vdx2[i][j] + d2vdy2[i][j]));
        }
    }
        free(H_old);
    }
    else
    {
        for (int i = 0; i<M; i++)
        {
            for (int j = 0; j<N-1; j++)
            {
                sol[i+1][j+1] = v[i+1][j+1] + dt * (-1*H_now[i][j] - dPdy[i][j] + 9.81 * beta*(T[i+1][j+1] - T0) + nu * (d2vdx2[i][j] + d2vdy2[i][j]));
            }
        }
    }
    free(dPdy);
    free(d2vdx2);
    free(d2vdy2);
    free(H_now);
}

void T_solve(double** T, double** u_old, double** v_old, double** u_now, double** v_now, double h, double dt, double q_w, double T_inf, double h_barre, double k, double alpha, int M, int N, int firstStep)
{
    //T est M+2xN+2, H est MxN, deja assemble avec AB2 ou euler, sol est M+2xN+2
    // on a vraiment besoin d un vecteur sol ? pourquoi ne pas directement overwrite T
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

    double** H_now  = calloc(M, sizeof( double *));
    for (int k = 0; k<M; k++)
    {
        H_now[k] = calloc(N,sizeof(double));
    }

    double* dTdy_e  = calloc(M, sizeof(double));

    d2Tdx2_fun(T, d2Tdx2, h, M, N);
    d2Tdy2_fun(T, d2Tdy2, h, M, N);
    AdvectiveT_fun(T,u_now,v_now,H_now,h,M,N);
    if (firstStep == 0)
    {
        double** H_old  = calloc(M, sizeof( double *));
        for (int k = 0; k<M; k++)
        {
            H_old[k] = calloc(N,sizeof(double));
        }
        AdvectiveT_fun(T,u_old,v_old,H_old,h,M,N);
        for (int i = 0; i<M; i++)
        {
            for (int j = 0; j<N; j++)
            {
                T[i+1][j+1] = T[i+1][j+1] + dt * ( -0.5 * (3*H_now[i][j]-H_old[i][j])  + alpha * (d2Tdx2[i][j] + d2Tdy2[i][j]));
            }
        }

        for (int k = 0; k<M; k++)
        {
            free(H_old[k]);
        }
        free(H_old);
    }
    else
    {
        for (int i = 0; i<M; i++)
        {
            for (int j = 0; j<N; j++)
            {
                T[i+1][j+1] = T[i+1][j+1] + dt * ( -1*H_now[i][j]  + alpha * (d2Tdx2[i][j] + d2Tdy2[i][j]));
            }
        }
    }

    //calculer dTdy_e avec les T au temps n ou n+1 ?
    double q_e = 0;
    for (int i=1; i<M+1; i++)
    {
        q_e = h_barre*((T[i][N+1]-T[i][N])/2 - T_inf);
        dTdy_e[i-1] = -q_e/k;
    }

    for (int i = 0; i<N; i++)
    {
        T[0][i+1] = T[1][i+1];
        T[M+1][i+1] = T[M][i+1];
    }

    double dTdy_w = -q_w/k;

    for (int i = 0; i<M; i++)
    {
        T[i+1][0] = T[i+1][1] - h*dTdy_w;
        T[i+1][N+1] = T[i+1][N] + h*dTdy_e[i];
    }

    for (int k = 0; k<M; k++)
        {
            free(H_now[k]);
            free(d2Tdx2[k]);
            free(d2Tdy2[k]);
        }
    free(H_now);
    free(d2Tdy2);
    free(d2Tdx2);
    free(dTdy_e);
}

void SOR(double** phi, double** ustar, double** vstar, double tol, double alpha, double H, double U, double L, double h, double dt, int M, int N)
{
    double** phiStar  = calloc(M, sizeof( double *));
    double** dustardx  = calloc(M, sizeof( double *));
    double** dvstardy  = calloc(M, sizeof( double *));
    double** d2phidx2  = calloc(M, sizeof( double *));
    double** d2phidy2  = calloc(M, sizeof( double *));
    for (int k = 0; k<M; k++)
    {
        phiStar[k] = calloc(N,sizeof(double));
        dustardx[k] = calloc(N,sizeof(double));
        dvstardy[k] = calloc(N,sizeof(double));
        d2phidx2[k] = calloc(N,sizeof(double));
        d2phidy2[k] = calloc(N,sizeof(double));
    }

    double R = 0;
    double sumR = 0;
    double error = 1;
    int iter = 0;
    dudx_fun(ustar,dustardx,h,M,N); //Erreur avec vstar!!!!@@@@@ attention!
    dvdy_fun(vstar,dvstardy,h,M,N);
    int leftWall;
    int bottomWall;

    while (error > tol && iter < 100)
    {
        sumR = 0;
        bottomWall = 1;
        for (int j=0; j<N; j++)
        {
            leftWall = 1;
            for (int i=0; i<M; i++)
            {
                phiStar[i][j] = 0.25* (-1* pow(h,2)*(dustardx[i][j]+dvstardy[i][j])/dt + phi[i+2][j+1] + phi[i][j+1] + phi[i+1][j+2] + phi[i+1][j]);
                //printf("dvstardy = %f\n",phiStar[i][j]);
                phi[i+1][j+1] = (alpha*phiStar[i][j] + (1-alpha) * phi[i+1][j+1]); //remplacer phistar direct ?
                leftWall = 0;

            }
            //cond limite paroies lat
            phi[0][j+1] = phi[1][j+1];
            phi[M+1][j+1] = phi[M][j+1];
            bottomWall = 0;
        }

        //cond limite paroies horizontales
        for (int i=0; i<M; i++)
        {
            phi[i+1][0] = phi[i+1][1];
            phi[i+1][N+1] = phi[i+1][N];
        }
        for (int j=0; j<N; j++)
        {
            for (int i=0; i<M; i++)
            {
                d2Tdx2_fun(phi,d2phidx2,h,M,N);
                d2Tdy2_fun(phi,d2phidy2,h,M,N);
                R = (d2phidx2[i][j]+d2phidy2[i][j]) - 1/dt*(dustardx[i][j]+dvstardy[i][j]);
                printf("%f ",phi[i+1][j+1]-phiStar[i][j]);
                sumR += R*R;
            }
            printf("\n");
        }
        error = (dt*H/U)*sqrt(1/(L*H)*sumR*h*h);
        printf("SOR global error = %f\n",error);
        printf("iter = %d\n",iter);
        //printf("SOR yoobal error = %f\n",alpha);
        iter++;


    }
}

void u_Solve(double** ustar, double** phi, double** sol, double dt, double h,int M, int N)
{
    double** dphidx  = calloc(M-1, sizeof( double *));
    for (int k = 0; k<M-1; k++)
    {
        dphidx[k] = calloc(N,sizeof(double));
    }

    dPdx_fun(phi,dphidx,h,M,N);

    for (int i = 0; i<M-1; i++)
    {
        for (int j = 0; j<N; j++)
        {
            sol[i+1][j+1] = ustar[i+1][j+1] - dt * dphidx[i][j];
        }
        sol[i+1][0] = -0.2*(sol[i+1][3]-5*sol[i+1][2]+15*sol[i+1][1]);
        sol[i+1][N+1] = sol[i+1][N];
    }

    for (int k = 0; k<M-1; k++)
    {
        free(dphidx[k]);
    }
    free(dphidx);
}

void v_Solve(double** vstar, double** phi, double** sol, double dt, double h,int M, int N)
{
    double** dphidy  = calloc(M, sizeof( double *));
    for (int k = 0; k<M; k++)
    {
        dphidy[k] = calloc(N-1,sizeof(double));
    }

    dPdy_fun(phi,dphidy,h,M,N);

    for (int j = 0; j<N-1; j++)
    {
        for (int i = 0; i<M; i++)
        {
            sol[i+1][j+1] = vstar[i+1][j+1] - dt * dphidy[i][j];
        }
        sol[0][j+1] = -0.2*(sol[3][j+1]-5*sol[2][j+1]+15*sol[1][j+1]);;
        sol[M+1][j+1] = -0.2*(sol[M+1-3][j+1]-5*sol[M+1-2][j+1]+15*sol[M+1-1][j+1]);;
    }

    for (int k = 0; k<M; k++)
    {
        free(dphidy[k]);
    }
    free(dphidy);
}

void P_solve(double** P, double** phi, int M, int N)
{
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            P[i+1][j+1] = P[i+1][j+1] + phi[i+1][j+1];
        }
        P[i+1][0] = P[i+1][1];
        P[i+1][N+1] = P[i+1][N];
    }
    for (int j = 0; j < N; j++)
    {
        P[0][j+1] = P[1][j+1];
        P[M+1][j+1] = P[M][j+1];
    }
}
