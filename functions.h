//Finite differences
void dPdx_fun(double** P, double** sol, double h, int M, int N);
void dPdy_fun(double** P, double** sol, double h, int M, int N);
void d2udx2_fun(double** u, double** sol, double h, int M, int N);
void d2udy2_fun(double** u, double** sol, double h, int M, int N);
void d2vdx2_fun(double** v, double** sol, double h, int M, int N); //IDEM QUE POUR U A SUPPRIMER
void d2vdy2_fun(double** v, double** sol, double h, int M, int N); //IDEM QUE POUR U A SUPPRIMER
void AdvectiveX_fun(double** u, double** v, double** a, double h, int M, int N);
void AdvectiveY_fun(double** u, double** v, double** b, double h, int M, int N);
void AB2X_fun(double** a_old, double** a_now, double** sol, int M, int N);
void AB2Y_fun(double** b_old, double** b_now, double** sol, int M, int N);
void AB2T_fun(double** H_old, double** H_now, double** sol, int M, int N);
void dudx_fun(double** u, double** sol, double h, int M, int N);
void dvdy_fun(double** v, double** sol, double h, int M, int N);
void d2Tdx2_fun(double** T, double** sol, double h, int N, int M);
void d2Tdy2_fun(double** T, double** sol, double h, int N, int M);
void AdvectiveT_fun(double** T, double** u, double** v, double** H, double h, int M, int N);
void dudy_fun(double** u, double** sol, double h, int M, int N);
void dvdx_fun(double** v, double** sol, double h, int M, int N);
double Vortex(double** u, double** v, double** vortex, double Gr, double h, int M, int N, double H);
double Reynolds(double** u, double** v, double** norm, double Gr, double h, int M, int N, double H);
double AverageT(double** T, double h, int M, int N);
void T_RMS(double** T, double T_rms, double T_avg, double h, int M, int N);

//Solvers
void ustar_Solve(double** u, double** v, double** u_old, double** v_old, double** P, double** sol, double h, double dt, double Gr, int M, int N, int firstStep);
void vstar_Solve(double** u, double** v, double** u_old, double** v_old, double** P, double** T, double** sol, double h, double dt, double Gr, int M, int N, int firstStep);
void SOR(double** phi, double** ustar, double** vstar, double tol, double alpha, double h, double dt, int M, int N);
void u_Solve(double** ustar, double** phi, double** sol, double dt, double h, int M, int N);
void v_Solve(double** vstar, double** phi, double** sol, double dt, double h, int M, int N);
void P_solve(double** P, double** phi, int M, int N);
void T_solve(double** T, double** T_old, double** u_old, double** v_old, double** u_now, double** v_now, double h, double dt, double Pr, double Gr, int M, int N, int firstStep);
