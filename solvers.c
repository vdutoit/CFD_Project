////DUTOIT Valentin - 44551400
////MECA2660_PROJECT_2018
#include<stdio.h>
#include<math.h>
#include"thomas.h"
#include"functions.h"

void ustar_Solve(double** u, double** v, double** u_old, double** v_old, double** P, double** sol, double h, double dt, double Gr, int M, int N, int firstStep)
{
    //u est M+1xN+2, a est M-1xN, P est M+2xN+2, sol est M+1xN+2, adv est M-1xN

    double** dPdx  = calloc(M-1, sizeof( double *));
    double** d2udx2  = calloc(M-1, sizeof( double *));
    double** d2udy2  = calloc(M-1, sizeof( double *));
    double** H_now  = calloc(M-1, sizeof( double *));
    for (int k = 0; k<M-1; k++)
    {
        dPdx[k] = calloc(N,sizeof(double));
        d2udx2[k] = calloc(N,sizeof(double));
        d2udy2[k] = calloc(N,sizeof(double));
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
                sol[i+1][j+1] = u[i+1][j+1] + dt * (-0.5*(3*H_now[i][j]-H_old[i][j]) - dPdx[i][j]  + 1/sqrt(Gr) * (d2udx2[i][j] + d2udy2[i][j]));
            }
        }
        for (int k = 0; k<M-1; k++)
        {
            free(H_old[k]);
        }
        free(H_old);
    }
    else
    {

        for (int i = 0; i<M-1; i++)
        {
            for (int j = 0; j<N; j++)
            {
                sol[i+1][j+1] = u[i+1][j+1] + dt * (-1*H_now[i][j] - dPdx[i][j]  + 1/sqrt(Gr) * (d2udx2[i][j] + d2udy2[i][j]));
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

void vstar_Solve(double** u, double** v, double** u_old, double** v_old, double** P, double** T, double** sol, double h, double dt, double Gr, int M, int N, int firstStep)
{
    //v est M+2xN+1, b est MxN-1, P est M+2xN+2, sol est M+2xN+1, adv est MxN-1
//    printf("checkpoint M = %d \n",M);
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
//    printf("checkpoint vstar dPdy\n");
    dPdy_fun(P,dPdy,h,M,N);
//    printf("checkpoint vstar d2vdx2\n");
    d2vdx2_fun(v,d2vdx2,h,M,N);
//    printf("checkpoint vstar d2vdy2\n");
    d2vdy2_fun(v,d2vdy2,h,M,N);
//    printf("checkpoint vstar AdvectiveY\n");
    AdvectiveY_fun(u,v,H_now,h,M,N);
//    printf("checkpoint vstar 2");

    if (firstStep == 0)
    {
//        printf("checkpoint vstar 2\n");
        double** H_old  = calloc(M, sizeof( double *));
        for (int k = 0; k<M; k++)
        {
            H_old[k] = calloc(N-1,sizeof(double));
        }
//        printf("checkpoint vstar 3\n");
        AdvectiveY_fun(u_old,v_old,H_old,h,M,N);
//        printf("checkpoint vstar 4\n");
        for (int i = 0; i<M; i++)
        {
            for (int j = 0; j<N-1; j++)
            {
                sol[i+1][j+1] = v[i+1][j+1] + dt * (-0.5*(3*H_now[i][j]-H_old[i][j]) - dPdy[i][j] + (T[i+1][j+2]+T[i+1][j+1])*0.5 + 1/sqrt(Gr) * (d2vdx2[i][j] + d2vdy2[i][j]));
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
        for (int m = 0; m<M; m++)
        {
            for (int n = 0; n<N-1; n++)
            {
//                printf("%d",n);//BIZARRE
                sol[m+1][n+1] = v[m+1][n+1] + dt * (-1*H_now[m][n] - dPdy[m][n] + (T[m+1][n+2]+T[m+1][n+1])*0.5 + 1/sqrt(Gr) * (d2vdx2[m][n] + d2vdy2[m][n]));
            }
        }
//        printf("checkpoint vstar elsefin\n");
    }

    //imposer conditions aux limites pour ustar -> 0 sur les paroies laterales (rien changer), u0 en y=0 et u[i][N+1] = u[i][N] sauf en i=0 et i=M

    for (int j = 1; j<N; j++)
    {
//        printf("checkpoint vstar 3\n");
        sol[0][j] = -0.2*(sol[3][j]-5*sol[2][j]+15*sol[1][j]); //no slip condition at left wall
//        printf("checkpoint vstar 32\n");
        sol[M+1][j] = -0.2*(sol[M+1-3][j]-5*sol[M+1-2][j]+15*sol[M+1-1][j]); //no slip condition at right wall
    }
//    printf("checkpoint vstar 4\n");
    for (int k = 0; k<M; k++)
    {
        free(dPdy[k]);
        free(d2vdx2[k]);
        free(d2vdy2[k]);
        free(H_now[k]);
    }
    free(dPdy);
    free(d2vdx2);
    free(d2vdy2);
    free(H_now);
}

void T_solve(double** T, double** T_old, double** u_old, double** v_old, double** u_now, double** v_now, double h, double dt, double Pr, double Gr, int M, int N, int firstStep)
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
        AdvectiveT_fun(T_old,u_old,v_old,H_old,h,M,N);

        for (int i = 0; i<M; i++)
        {
            for (int j = 0; j<N; j++)
            {
                T[i+1][j+1] = T[i+1][j+1] + dt * (-0.5 * (3*H_now[i][j]-H_old[i][j])  + 1/(Pr*sqrt(Gr)) * (d2Tdx2[i][j] + d2Tdy2[i][j]));
//                printf("pr = %f ",(d2Tdy2[i][j]));
            }
//            printf("\n");
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
                T[i+1][j+1] = T[i+1][j+1] + dt * (-1*H_now[i][j] + 1/(Pr*sqrt(Gr)) * (d2Tdx2[i][j] + d2Tdy2[i][j]));
            }
        }
    }

    for (int i = 0; i<M; i++)
    {
        T[i+1][N+1] = ((1/h - 1e3/2.0)*T[i+1][N] - 5)/(1/h + 1e3/2.0);
        T[i+1][0] = T[i+1][1] + h;
    }
    for (int j = 0; j<N; j++)
    {
        T[0][j+1] = T[1][j+1];
        T[M+1][j+1] = T[M][j+1];
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


void SOR(double** phi, double** ustar, double** vstar, double tol, double alpha, double h, double dt, int M, int N)
{
//    printf("checkpoint calloc\n");
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
//    printf("checkpoint div\n");
    dudx_fun(ustar,dustardx,h,M,N); //Erreur avec vstar!!!!@@@@@ attention!
    dvdy_fun(vstar,dvstardy,h,M,N);
    int leftWall;
    int bottomWall;

    while (error > tol)
    {
        sumR = 0;
        bottomWall = 1;
//        printf("checkpoint algo\n");
        for (int j=0; j<N; j++)
        {
            leftWall = 1;
            for (int i=0; i<M; i++)
            {
                phiStar[i][j] = 0.25 * (-1* pow(h,2)*(dustardx[i][j]+dvstardy[i][j])/dt + phi[i+2][j+1] + (1-leftWall)*phi[i][j+1] + phi[i+1][j+2] + (1-bottomWall)*phi[i+1][j]);
                phi[i+1][j+1] = (alpha*phiStar[i][j] + (1-alpha) * phi[i+1][j+1])/(1-0.25*alpha*(bottomWall+leftWall)); //remplacer phistar direct ?
                leftWall = 0;
            }
            phi[0][j+1] = phi[1][j+1];
            phi[M+1][j+1] = phi[M][j+1];
            bottomWall = 0;
        }
//        printf("checkpoint algofin\n");
        //cond limite paroies horizontales
//        printf("checkpoint CL\n");
        for (int i=0; i<M; i++)
        {
            phi[i+1][0] = phi[i+1][1];
            phi[i+1][N+1] = phi[i+1][N];
        }

        d2Tdx2_fun(phi,d2phidx2,h,M,N);
        d2Tdy2_fun(phi,d2phidy2,h,M,N);
        for (int j=0; j<N; j++)
        {
            for (int i=0; i<M; i++)
            {
                R = (d2phidx2[i][j]+d2phidy2[i][j]) - 1/dt*(dustardx[i][j]+dvstardy[i][j]);
                sumR += R*R;
            }
        }
        error = (dt)*sqrt(1.5*sumR*h*h);
        printf("SOR global error = %f, iter: %d\n",error,iter);
        iter++;
    }

    for (int k = 0; k<M; k++)
    {
        free(phiStar[k]);
        free(dustardx[k]);
        free(dvstardy[k]);
        free(d2phidx2[k]);
        free(d2phidy2[k]);
    }
    free(phiStar);
    free(dustardx);
    free(dvstardy);
    free(d2phidx2);
    free(d2phidy2);
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
