#ifndef SOLVER_H_
#define SOLVER_H_

void coefficient_matrix (double E[], double F[], double D[], double B[], double H[], double f[], double d[], double b[], double e[], double c[], int K, int nx, double ALPHA);
double SIP_Solver_x (double E[], double F[], double D[], double B[], double H[], double Q[], double f[], double d[], double b[], double e[], double c[], double res[], double vec_Y[], double vec_delta[], int K, int nx, int inner, double variable[], double **original, double *coeff_p);
double SIP_Solver_y (double E[], double F[], double D[], double B[], double H[], double Q[], double f[], double d[], double b[], double e[], double c[], double res[], double vec_Y[], double vec_delta[], int K, int nx, int inner, double variable[], double **original, double *coeff_p);
double SIP_Solver_Centre (double E[], double F[], double D[], double B[], double H[], double Q[], double f[], double d[], double b[], double e[], double c[], double res[], double vec_Y[], double vec_delta[], int K, int nx, double variable[], double **original);
double BiCGSTAB_Pre(double *p, double *rtilde, double *res_old, double *s, double *v, double *t, double *Ep, double *Fp, double *Hp, double *Dp, double *Bp, double *variable, double *res, double *vec_Y_t, double *vec_Y_s, double *vec_delta_p, double *vec_delta_s, double *b, double *c, double *d, double *e, double *f, int K, int nx, double **original);
double BiCGSTAB_Pre_y(double *p, double *rtilde, double *res_old, double *s, double *v, double *t, double *Ep, double *Fp, double *Hp, double *Dp, double *Bp, double *variable, double *res, double *vec_Y_t, double *vec_Y_s, double *vec_delta_p, double *vec_delta_s, double *b, double *c, double *d, double *e, double *f, int K, int nx, int inner, double **original, double *coeff_p);
double BiCGSTAB_Pre_x(double *p, double *rtilde, double *res_old, double *s, double *v, double *t, double *Ep, double *Fp, double *Hp, double *Dp, double *Bp, double *variable, double *res, double *vec_Y_t, double *vec_Y_s, double *vec_delta_p, double *vec_delta_s, double *b, double *c, double *d, double *e, double *f, int K, int nx, int inner, double **original, double *coeff_p);
void F_B(double *vec_Y, double *res, double *vec_delta, double *b, double *c, double *d, double *e, double *f, int nx, int K);
void L(double *vec_Y, double *res, double *b, double *c, double *d, double *e, double *f, int nx, int K);
double VC_x(int i, int j, double **uc_ap, double **ud_ap, double **alpha_c, double **alpha_d, double **B_kx, int is_it_continuous);
double VC_y(int i, int j, double **vc_ap, double **vd_ap, double **alpha_c, double **alpha_d, double **B_ky, int is_it_continuous);
double PCC_x(int i, int j, double **uc_ap, double **ud_ap, double alpha_c, double alpha_d, double **B_kx, int is_it_continuous);
double PCC_y(int i, int j, double **vc_ap, double **vd_ap, double alpha_c, double alpha_d, double **B_ky, int is_it_continuous);
double psi_pos_x(double **alpha_d, int i, int j, int is_it_boundary, int east_or_west);
double psi_neg_x(double **alpha_d, int i, int j, int is_it_boundary, int east_or_west);
double psi_pos_y(double **alpha_d, int i, int j, int is_it_boundary, int north_or_south, double *alphad_inlet);
double psi_neg_y(double **alpha_d, int i, int j, int is_it_boundary, int north_or_south);
double psi_y_pos_x(double **alpha_d, int i, int j, int is_it_boundary, int east_or_west);
double psi_y_neg_x(double **alpha_d, int i, int j, int is_it_boundary, int east_or_west);
double psi_y_pos_y(double **alpha_d, int i, int j, int is_it_boundary, int north_or_south, double *alphad_inlet);
double psi_y_neg_y(double **alpha_d, int i, int j, int is_it_boundary, int north_or_south);
double psi_x_pos_x(double **alpha_d, int i, int j, int is_it_boundary, int east_or_west);
double psi_x_neg_x(double **alpha_d, int i, int j, int is_it_boundary, int east_or_west);
double psi_x_pos_y(double **alpha_d, int i, int j, int is_it_boundary, int north_or_south, double *alphad_inlet);
double psi_x_neg_y(double **alpha_d, int i, int j, int is_it_boundary, int north_or_south);
double min_mod(double r);
double van_leer(double r);
void ILU0(double E[], double F[], double D[], double B[], double H[], double f[], double d[], double b[], double e[], double c[], int K, int nx);
void MATRIX_VECTOR_MULTIPLY_Y(double *Ep, double *Dp, double *Fp, double *Hp, double *Bp, double *vec_delta_p, double *v, int nx);
void MATRIX_VECTOR_MULTIPLY_X(double *Ep, double *Dp, double *Fp, double *Hp, double *Bp, double *vec_delta_p, double *v, int nx);
void MATRIX_VECTOR_MULTIPLY_Center(double *Ep, double *Dp, double *Fp, double *Hp, double *Bp, double *vec_delta_p, double *v, int nx);
void residual_center(double *Q, double *E, double *F, double *H, double *D, double *B, double *res, double *variable, int nx);
void residual_X(double *Q, double *E, double *F, double *H, double *D, double *B, double *res, double *variable, double **original, int nx);
void residual_Y(double *Q, double *E, double *F, double *H, double *D, double *B, double *res, double *variable, double **original, int nx);


#endif