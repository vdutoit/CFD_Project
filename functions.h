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
void Div_star_fun(double** ustar, double** vstar, double** sol, double h, int M, int N); //vraiment utile ? je pense pas autant faire direct dudx + dvdy
void d2Tdx2_fun(double** T, double** sol, double h, int N, int M);
void d2Tdy2_fun(double** T, double** sol, double h, int N, int M);
void AdvectiveT_fun(double** T, double** u, double** v, double** H, double h, int M, int N);

//Solvers
void ustar_Solve(double** u, double** v, double** adv, double** P, double** sol, double h, double dt, double nu, int M, int N);
void vstar_Solve(double** u, double** v, double** adv, double** P, double** T, double** sol, double h, double dt, double T0, double nu, double beta, int M, int N); //AJOUTER TEMPERATURE
void T_solve(double** T, double** H, double h, double dt, double q_w, double T_inf, double h_barre, double k, double alpha, int M, int N);
void SOR(double** phi, double** ustar, double** vstar, double tol, double alpha, double h, double dt, int M, int N);
