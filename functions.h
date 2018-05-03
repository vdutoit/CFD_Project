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
void getNorm(double** u, double** v, double** sol, double h, int M, int N);

//Solvers
void ustar_Solve(double** u, double** v, double** u_old, double** v_old, double** P, double** sol, double h, double dt, double nu, int M, int N, int firstStep);
void vstar_Solve(double** u, double** v, double** u_old, double** v_old, double** P, double** T, double** sol, double h, double dt, double T0, double nu, double beta, int M, int N, int firstStep);
void SOR(double** phi, double** ustar, double** vstar, double tol, double alpha, double H, double U, double L, double h, double dt, int M, int N);
void u_Solve(double** ustar, double** phi, double** sol, double dt, double h, int M, int N);
void v_Solve(double** vstar, double** phi, double** sol, double dt, double h, int M, int N);
void P_solve(double** P, double** phi, int M, int N);
void T_solve(double** T, double** u_old, double** v_old, double** u_now, double** v_now, double h, double dt, double q_w, double T_inf, double h_barre, double k, double alpha, int M, int N, int firstStep);
