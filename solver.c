#include "functions.h"
#include <math.h>

#define DX2 2e-2 * 2e-2
#define UV_COR 0.7
#define NX 25
#define NY 75
#define min(A, B) ( (A) < (B) ? (A) : (B) )

double min_mod(double r)
{
    if(r > 0)
        return min(r,1);
    else
        return 0;
}


double BiCGSTAB_Pre(double *p, double *rtilde, double *res_old, double *s, double *v, double *t, double *Ep, double *Fp, double *Hp, double *Dp, double *Bp, double *variable, double *res, double *vec_Y_t, double *vec_Y_s, double *vec_delta_p, double *vec_delta_s, double *b, double *c, double *d, double *e, double *f, int K, int nx, double **original)
{
    int i,j,m;

    F_B(vec_Y_t, p, vec_delta_p, b, c, d, e, f, nx, K);

    //// Calculating alpha

    double numer_p = 0;

    for (m=0;m<K;m++)
    {
        numer_p += res_old[m]*rtilde[m];
    }

    MATRIX_VECTOR_MULTIPLY_Center(Ep, Dp, Fp, Hp, Bp, vec_delta_p, v, nx);


    double denom_p = 0;

    for (m=0;m<K;m++)
    {
        denom_p += v[m]*rtilde[m];
    }

    double alpha = numer_p/denom_p;


    //// Calculating G vector

    for (m=0;m<K;m++)
    {
        s[m] = res_old[m] - alpha*v[m];
    }


    F_B(vec_Y_s, s, vec_delta_s, b, c, d, e, f, nx, K);


    MATRIX_VECTOR_MULTIPLY_Center(Ep, Dp, Fp, Hp, Bp, vec_delta_s, t, nx);


    L(vec_Y_t, t, b, c, d, e, f, nx, K);

    L(vec_Y_s, s, b, c, d, e, f, nx, K);


    double numer = 0;

    for (i=0;i<K;i++)
    {
        numer += vec_Y_t[i]*vec_Y_s[i];
    }

    double denom = 0;

    for (i=0;i<K;i++)
    {
        denom += vec_Y_t[i]*vec_Y_t[i];
    }

    double omega = numer/denom;

    //// UpdatinGp Solution

    for (i=0;i<K;i++)
    {
        variable[i] += alpha*vec_delta_p[i] + omega*vec_delta_s[i];
    }


    for (i=0;i<K;i++)
    {
        res[i] = s[i] - omega*t[i];
    }


    double res_sum = 0;

    for (m=0;m<K;m++)
    {
        res_sum += res[m]*res[m];
    }


    double R2 = sqrt(res_sum);


    //// ComputinGp beta

    double betanumer = 0;

    for (m=0;m<K;m++)
    {
        betanumer += res[m]*rtilde[m];
    }

    double beta = (alpha/omega)*(betanumer/numer_p);

    //// UpdatinGp ConjuGpate direction vector

    for (m=0;m<K;m++)
    {
        p[m] = res[m] + beta*( p[m] - omega*v[m] );
    }


    for (m=0;m<K;m++)
    {
        res_old[m] = res[m];
    }

    for (i=0;i<NX;i++)
    {
        for (j=0;j<NY;j++)
        {
            m = j*nx + i;

            original[i][j] = variable[m];
        }
    }

    return R2;
}


double BiCGSTAB_Pre_x(double *p, double *rtilde, double *res_old, double *s, double *v, double *t, double *Ep, double *Fp, double *Hp, double *Dp, double *Bp, double *variable, double *res, double *vec_Y_t, double *vec_Y_s, double *vec_delta_p, double *vec_delta_s, double *b, double *c, double *d, double *e, double *f, int K, int nx, int inner, double **original, double *coeff_p)
{
    int i,j,m;

    F_B(vec_Y_t, p, vec_delta_p, b, c, d, e, f, nx, K);

    //// Calculating alpha

    double numer_p = 0;

    for (m=0;m<K;m++)
    {
        numer_p += res_old[m]*rtilde[m];
    }

    MATRIX_VECTOR_MULTIPLY_X(Ep, Dp, Fp, Hp, Bp, vec_delta_p, v, nx);


    double denom_p = 0;

    for (m=0;m<K;m++)
    {
        denom_p += v[m]*rtilde[m];
    }

    double alpha = numer_p/denom_p;


    //// Calculating G vector

    for (m=0;m<K;m++)
    {
        s[m] = res_old[m] - alpha*v[m];
    }


    F_B(vec_Y_s, s, vec_delta_s, b, c, d, e, f, nx, K);


    MATRIX_VECTOR_MULTIPLY_X(Ep, Dp, Fp, Hp, Bp, vec_delta_s, t, nx);

    L(vec_Y_t, t, b, c, d, e, f, nx, K);

    L(vec_Y_s, s, b, c, d, e, f, nx, K);


    double numer = 0;

    for (i=0;i<K;i++)
    {
        numer += vec_Y_t[i]*vec_Y_s[i];
    }

    double denom = 0;

    for (i=0;i<K;i++)
    {
        denom += vec_Y_t[i]*vec_Y_t[i];
    }

    double omega = numer/denom;

    //// UpdatinGp Solution

    for (i=0;i<K;i++)
    {
        variable[i] += alpha*vec_delta_p[i] + omega*vec_delta_s[i];
    }


    for (i=0;i<K;i++)
    {
        res[i] = s[i] - omega*t[i];
    }


    double res_sum = 0;

    for (m=0;m<K;m++)
    {
        res_sum += res[m]*res[m];
    }


    double R2 = sqrt(res_sum);


    //// ComputinGp beta

    double betanumer = 0;

    for (m=0;m<K;m++)
    {
        betanumer += res[m]*rtilde[m];
    }

    double beta = (alpha/omega)*(betanumer/numer_p);

    //// UpdatinGp ConjuGpate direction vector

    for (m=0;m<K;m++)
    {
        p[m] = res[m] + beta*( p[m] - omega*v[m] );
    }


    for (m=0;m<K;m++)
    {
        res_old[m] = res[m];
    }

    for (i=1;i<NX;i++)
    {
        for (j=0;j<NY;j++)
        {
            m = j*nx + i-1;

            original[i][j] = variable[m];

        }
    }

    double nor = 0;
    
    for (i=1;i<NX;i++)
    {
        for(j=0;j<NY;j++)
        {
            m = j*nx + i-1;

            nor += fabs(coeff_p[m]*original[i][j]);
        }
    }
    
    R2 = R2/nor;

    return R2;

}


double BiCGSTAB_Pre_y(double *p, double *rtilde, double *res_old, double *s, double *v, double *t, double *Ep, double *Fp, double *Hp, double *Dp, double *Bp, double *variable, double *res, double *vec_Y_t, double *vec_Y_s, double *vec_delta_p, double *vec_delta_s, double *b, double *c, double *d, double *e, double *f, int K, int nx, int inner, double **original, double *coeff_p)
{
    int i,j,m;

    F_B(vec_Y_t, p, vec_delta_p, b, c, d, e, f, nx, K);

    //// Calculating alpha

    double numer_p = 0;

    for (m=0;m<K;m++)
    {
        numer_p += res_old[m]*rtilde[m];
    }

    MATRIX_VECTOR_MULTIPLY_Y(Ep, Dp, Fp, Hp, Bp, vec_delta_p, v, nx);


    double denom_p = 0;

    for (m=0;m<K;m++)
    {
        denom_p += v[m]*rtilde[m];
    }

    double alpha = numer_p/denom_p;


    //// Calculating G vector

    for (m=0;m<K;m++)
    {
        s[m] = res_old[m] - alpha*v[m];
    }


    F_B(vec_Y_s, s, vec_delta_s, b, c, d, e, f, nx, K);

    MATRIX_VECTOR_MULTIPLY_Y(Ep, Dp, Fp, Hp, Bp, vec_delta_s, t, nx);


    L(vec_Y_t, t, b, c, d, e, f, nx, K);

    L(vec_Y_s, s, b, c, d, e, f, nx, K);


    double numer = 0;

    for (i=0;i<K;i++)
    {
        numer += vec_Y_t[i]*vec_Y_s[i];
    }

    double denom = 0;

    for (i=0;i<K;i++)
    {
        denom += vec_Y_t[i]*vec_Y_t[i];
    }

    double omega = numer/denom;

    //// UpdatinGp Solution

    for (i=0;i<K;i++)
    {
        variable[i] += alpha*vec_delta_p[i] + omega*vec_delta_s[i];
    }


    for (i=0;i<K;i++)
    {
        res[i] = s[i] - omega*t[i];
    }


    double res_sum = 0;

    for (m=0;m<K;m++)
    {
        res_sum += res[m]*res[m];
    }


    double R2 = sqrt(res_sum);


    //// ComputinGp beta

    double betanumer = 0;

    for (m=0;m<K;m++)
    {
        betanumer += res[m]*rtilde[m];
    }

    double beta = (alpha/omega)*(betanumer/numer_p);

    //// UpdatinGp ConjuGpate direction vector

    for (m=0;m<K;m++)
    {
        p[m] = res[m] + beta*( p[m] - omega*v[m] );
    }


    for (m=0;m<K;m++)
    {
        res_old[m] = res[m];
    }


    for (i=0;i<NX;i++)
    {
        for (j=1;j<NY;j++)
        {
            m = (j-1)*nx + i;

            original[i][j] = variable[m];

        }
    }

    double nor = 0;
    
    for (i=0;i<NX;i++)
    {
        for(j=1;j<NY;j++)
        {
            m = (j-1)*nx + i;

            nor += fabs(coeff_p[m]*original[i][j]);
        }
    }
    
    R2 = R2/nor;

    return R2;

}



void F_B(double *vec_Y, double *res, double *vec_delta, double *b, double *c, double *d, double *e, double *f, int nx, int K)
{

    int i;
    
    //// Forward-Sub

    vec_Y[0] = res[0];

    for (i=1;i<K;i++)
    {
        if (i>=nx)
        {
            vec_Y[i] = ( res[i] - c[i]*vec_Y[i-1] - b[i]*vec_Y[i-nx] );
        }
        else
        {
            vec_Y[i] = ( res[i] - c[i]*vec_Y[i-1] );
        }
    }


    //// Backward-Sub

    vec_delta[K-1] = vec_Y[K-1]/d[K-1];

    for (i=K-2;i>=0;i--)
    {
        if (i<=K-1-nx)
        {    
            vec_delta[i] = ( vec_Y[i] - e[i]*vec_delta[i+1] - f[i]*vec_delta[i+nx] )/d[i];
        }
        else
        {   
            vec_delta[i] = ( vec_Y[i] - e[i]*vec_delta[i+1] )/d[i];
        }

    }

}

void L(double *vec_Y, double *res, double *b, double *c, double *d, double *e, double *f, int nx, int K)
{
    int i;
    
    //// Forward Sub

    vec_Y[0] = res[0];

    for (i=1;i<K;i++)
    {
        if (i>=nx)
        {
            vec_Y[i] = ( res[i] - c[i]*vec_Y[i-1] - b[i]*vec_Y[i-nx] );
        }
        else
        {
            vec_Y[i] = ( res[i] - c[i]*vec_Y[i-1] );
        }
    }
}

double SIP_Solver_x (double E[], double F[], double D[], double B[], double H[], double Q[], double f[], double d[], double b[], double e[], double c[], double res[], double vec_Y[], double vec_delta[], int K, int nx, int inner, double variable[], double **original, double *coeff_p)
{
    int i,j,m;

    residual_X(Q, E, F, H, D, B, res, variable, original, nx);

    double res_sum = 0;

    for (m=0;m<K;m++)
    {
        res_sum += res[m]*res[m];
    }

    double R2 = sqrt(res_sum);

    //// Forward Sub


    vec_Y[0] = res[0]/d[0];

    for (i=1;i<K;i++)
    {
        if (i>=nx)
        {
            vec_Y[i] = ( res[i] - c[i]*vec_Y[i-1] - b[i]*vec_Y[i-nx] )/d[i];
        }
        else
        {
            vec_Y[i] = ( res[i] - c[i]*vec_Y[i-1] )/d[i];
        }
    }

    //// Backward sub

    vec_delta[K-1] = vec_Y[K-1];

    for (i=K-2;i>=0;i--)
    {
        if (i<=K-1-nx)
        {    
            vec_delta[i] = ( vec_Y[i] - e[i]*vec_delta[i+1] - f[i]*vec_delta[i+nx] );
        }
        else
        {
            vec_delta[i] = ( vec_Y[i] - e[i]*vec_delta[i+1] );
        }

    }

    for (i=0;i<K;i++)
    {
        variable[i] += vec_delta[i];
    }

    for (i=1;i<NX;i++)
    {
        for (j=0;j<NY;j++)
        {
            m = j*nx + i-1;

            original[i][j] = variable[m];

        }
    }

    double nor = 0;
    
    for (i=1;i<NX;i++)
    {
        for(j=0;j<NY;j++)
        {
            m = j*nx + i-1;

            nor += fabs(coeff_p[m]*original[i][j]);
        }
    }
    
    R2 = R2/nor;

    return R2;
}


double SIP_Solver_y (double E[], double F[], double D[], double B[], double H[], double Q[], double f[], double d[], double b[], double e[], double c[], double res[], double vec_Y[], double vec_delta[], int K, int nx, int inner, double variable[], double **original, double *coeff_p)
{
    int i,j,m;

    residual_Y(Q, E, F, H, D, B, res, variable, original, nx);


    double res_sum = 0;

    for (m=0;m<K;m++)
    {
        res_sum += res[m]*res[m];
    }

    double R2 = sqrt(res_sum);

    //// Forward Sub

    vec_Y[0] = res[0]/d[0];

    for (i=1;i<K;i++)
    {
        if (i>=nx)
        {
            vec_Y[i] = ( res[i] - c[i]*vec_Y[i-1] - b[i]*vec_Y[i-nx] )/d[i];
        }
        else
        {
            vec_Y[i] = ( res[i] - c[i]*vec_Y[i-1] )/d[i];
        }
    }

    //// Backward sub

    vec_delta[K-1] = vec_Y[K-1];

    for (i=K-2;i>=0;i--)
    {
        if (i<=K-1-nx)
        {    
            vec_delta[i] = ( vec_Y[i] - e[i]*vec_delta[i+1] - f[i]*vec_delta[i+nx] );
        }
        else
        {
            vec_delta[i] = ( vec_Y[i] - e[i]*vec_delta[i+1] );
        }

    }

    for (i=0;i<K;i++)
    {
        variable[i] += vec_delta[i];
    }

    for (i=0;i<NX;i++)
    {
        for (j=1;j<NY;j++)
        {
            m = (j-1)*nx + i;

            original[i][j] = variable[m];

        }
    }

    double nor = 0;
    
    for (i=0;i<NX;i++)
    {
        for(j=1;j<NY;j++)
        {
            m = (j-1)*nx + i;

            nor += fabs(coeff_p[m]*original[i][j]);
        }
    }
    
    R2 = R2/nor;

    return R2;
}


double SIP_Solver_Centre (double E[], double F[], double D[], double B[], double H[], double Q[], double f[], double d[], double b[], double e[], double c[], double res[], double vec_Y[], double vec_delta[], int K, int nx, double variable[], double **original)
{
    int i,j,m;

    residual_center(Q, E, F, H, D, B, res, variable, nx);

    double res_sum = 0;

    for (m=0;m<K;m++)
    {
        res_sum += res[m]*res[m];
    }

    double R2 = sqrt(res_sum);

    //// Forward Sub

    vec_Y[0] = res[0]/d[0];

    for (i=1;i<K;i++)
    {
        if (i>=nx)
        {
            vec_Y[i] = ( res[i] - c[i]*vec_Y[i-1] - b[i]*vec_Y[i-nx] )/d[i];
        }
        else
        {
            vec_Y[i] = ( res[i] - c[i]*vec_Y[i-1] )/d[i];
        }
    }

    //// Bckward sub

    vec_delta[K-1] = vec_Y[K-1];

    for (i=K-2;i>=0;i--)
    {
        if (i<=K-1-nx)
        {    
            vec_delta[i] = ( vec_Y[i] - e[i]*vec_delta[i+1] - f[i]*vec_delta[i+nx] );
        }
        else
        {   
            vec_delta[i] = ( vec_Y[i] - e[i]*vec_delta[i+1] );
        }

    }


    for (i=0;i<K;i++)
    {
        variable[i] += vec_delta[i];
    }

    for (i=0;i<NX;i++)
    {
        for (j=0;j<NY;j++)
        {
            m = j*nx + i;

            original[i][j] = variable[m];
        }
    }

    return R2;

}

void coefficient_matrix (double E[], double F[], double D[], double B[], double H[], double f[], double d[], double b[], double e[], double c[], int K, int nx, double ALPHA)
{
    //// Filling out L and U matrix

    d[0] = E[0];
    e[0] = F[0]/d[0];
    f[0] = H[0]/d[0];

    for (int l=1;l<K;l++)
    {
        if (l>=nx)
        {
            b[l] = B[l]/( 1 + ALPHA*e[l-nx] );
        }
        
        c[l] = D[l]/( 1 + ALPHA*f[l-1] );

        if (l>=nx)
        {
            d[l] = E[l] + ALPHA*( b[l]*e[l-nx] + c[l]*f[l-1] ) - b[l]*f[l-nx] - c[l]*e[l-1]; 
        }
        else
        {
            d[l] = E[l] + ALPHA*( c[l]*f[l-1] ) - c[l]*e[l-1];
        }

        if (l>=nx)
        {
            e[l] = ( F[l] - ALPHA*b[l]*e[l-nx] )/d[l];
        }
        else
        {
            e[l] = F[l]/d[l];
        }

        f[l] = ( H[l] - ALPHA*c[l]*f[l-1] )/d[l];
    }


}


double VC_x(int i, int j, double **uc_ap, double **ud_ap, double **alpha_c, double **alpha_d, double **B_kx, int is_it_continuous)
{
    double alphac_centre = (alpha_c[i-1][j] + alpha_c[i][j])/2;
    double alphad_centre = (alpha_d[i-1][j] + alpha_d[i][j])/2;


    if (is_it_continuous)
    {
        return UV_COR*( alphac_centre*DX2*( ud_ap[i][j] + B_kx[i][j] ) + alphad_centre*DX2*B_kx[i][j] )/( uc_ap[i][j]*ud_ap[i][j] + B_kx[i][j]*( uc_ap[i][j] + ud_ap[i][j] ) );
    }
    else
    {
        return UV_COR*( alphad_centre*DX2*( uc_ap[i][j] + B_kx[i][j] ) + alphac_centre*DX2*B_kx[i][j] )/( uc_ap[i][j]*ud_ap[i][j] + B_kx[i][j]*( uc_ap[i][j] + ud_ap[i][j] ) );
    }
}

double VC_y(int i, int j, double **vc_ap, double **vd_ap, double **alpha_c, double **alpha_d, double **B_ky, int is_it_continuous)
{
    double alphac_centre = (alpha_c[i][j-1] + alpha_c[i][j])/2;
    double alphad_centre = (alpha_d[i][j-1] + alpha_d[i][j])/2;


    if (is_it_continuous)
    {
        return UV_COR*( alphac_centre*DX2*( vd_ap[i][j] + B_ky[i][j] ) + alphad_centre*DX2*B_ky[i][j] )/( vc_ap[i][j]*vd_ap[i][j] + B_ky[i][j]*( vc_ap[i][j] + vd_ap[i][j] ) );
    }
    else
    {
        return UV_COR*( alphad_centre*DX2*( vc_ap[i][j] + B_ky[i][j] ) + alphac_centre*DX2*B_ky[i][j] )/( vc_ap[i][j]*vd_ap[i][j] + B_ky[i][j]*( vc_ap[i][j] + vd_ap[i][j] ) );
    }
}

double PCC_x(int i, int j, double **uc_ap, double **ud_ap, double alpha_c, double alpha_d, double **B_kx, int is_it_continuous)
{
    
    if (is_it_continuous)
    {
        return UV_COR*( alpha_c*DX2*( ud_ap[i][j] + B_kx[i][j] ) + alpha_d*DX2*B_kx[i][j] )/( uc_ap[i][j]*ud_ap[i][j] + B_kx[i][j]*( uc_ap[i][j] + ud_ap[i][j] ) );
    }
    else
    {
        return UV_COR*( alpha_d*DX2*( uc_ap[i][j] + B_kx[i][j] ) + alpha_c*DX2*B_kx[i][j] )/( uc_ap[i][j]*ud_ap[i][j] + B_kx[i][j]*( uc_ap[i][j] + ud_ap[i][j] ) );
    }

}

double PCC_y(int i, int j, double **vc_ap, double **vd_ap, double alpha_c, double alpha_d, double **B_ky, int is_it_continuous)
{

    if (is_it_continuous)
    {
        return UV_COR*( alpha_c*DX2*( vd_ap[i][j] + B_ky[i][j] ) + alpha_d*DX2*B_ky[i][j] )/( vc_ap[i][j]*vd_ap[i][j] + B_ky[i][j]*( vc_ap[i][j] + vd_ap[i][j] ) );
    }
    else
    {
        return UV_COR*( alpha_d*DX2*( vc_ap[i][j] + B_ky[i][j] ) + alpha_c*DX2*B_ky[i][j] )/( vc_ap[i][j]*vd_ap[i][j] + B_ky[i][j]*( vc_ap[i][j] + vd_ap[i][j] ) );
    }

}


double psi_pos_x(double **alpha_d, int i, int j, int is_it_boundary, int east_or_west)
{
    double r_pos;

    if( (east_or_west == 0 && is_it_boundary == 1) || i == 0 || (alpha_d[i+1][j] - alpha_d[i][j]) == 0)
    {
        if((alpha_d[i+1][j] - alpha_d[i][j]) == 0)
            return 1;
        else
            return 0;
    }
    else
    {
        r_pos = ( alpha_d[i][j] - alpha_d[i-1][j] )/( alpha_d[i+1][j] - alpha_d[i][j] );
        return min_mod(r_pos);
    }
}

double psi_neg_x(double **alpha_d, int i, int j, int is_it_boundary, int east_or_west)
{
    double r_neg;

    if( ( east_or_west == 1 && is_it_boundary == 1 ) || i == NX-1 || ( alpha_d[i][j] - alpha_d[i-1][j] ) == 0)
    {
        if(( alpha_d[i][j] - alpha_d[i-1][j] ) == 0 )
            return 1;
        else
            return 0;
    }
    else
    {
        r_neg = ( alpha_d[i+1][j] - alpha_d[i][j] )/( alpha_d[i][j] - alpha_d[i-1][j] );
        return min_mod(r_neg);
    }
}

double psi_pos_y(double **alpha_d, int i, int j, int is_it_boundary, int north_or_south, double *alphad_inlet)
{
    double r_pos;

    if( (north_or_south == 0 && is_it_boundary == 1) || j == 0 || ( alpha_d[i][j+1] - alpha_d[i][j] ) == 0)
    {
        if(( alpha_d[i][j+1] - alpha_d[i][j] ) == 0 )
            return 1;
        else
        {
            r_pos = 2*( alpha_d[i][j] - alphad_inlet[i] )/( alpha_d[i][j+1] - alpha_d[i][j] );
            return min_mod(r_pos);

        }

    }
    else
    {
        r_pos = ( alpha_d[i][j] - alpha_d[i][j-1] )/( alpha_d[i][j+1] - alpha_d[i][j] );
        return min_mod(r_pos);
    }
}

double psi_neg_y(double **alpha_d, int i, int j, int is_it_boundary, int north_or_south)
{
    double r_neg;

    if((north_or_south == 1 && is_it_boundary == 1) || j == NY-1 || ( alpha_d[i][j] - alpha_d[i][j-1] ) == 0)
    {
        if(( alpha_d[i][j] - alpha_d[i][j-1] ) == 0 )
            return 1;
        else 
            return 0;
    }
    else
    {
        r_neg = ( alpha_d[i][j+1] - alpha_d[i][j] )/( alpha_d[i][j] - alpha_d[i][j-1] );
        return min_mod(r_neg);
    }
}


double psi_y_pos_x(double **alpha_d, int i, int j, int is_it_boundary, int east_or_west)             //// function for y-momentum equation for finding flux limiter in x-direction
{
    double r_pos;

    if((east_or_west == 0 && is_it_boundary == 1) || i == 0 || (alpha_d[i+1][j] - alpha_d[i][j]) == 0)
    {
        if ((alpha_d[i+1][j] - alpha_d[i][j]) == 0)
            return 1;
        else
        {
            r_pos = 2*( alpha_d[i][j] )/( alpha_d[i+1][j] - alpha_d[i][j] );
            return min_mod(r_pos);
        }
        return 0;
    }   
    else
    {
        r_pos = ( alpha_d[i][j] - alpha_d[i-1][j] )/( alpha_d[i+1][j] - alpha_d[i][j] );
        return min_mod(r_pos);
    }
}


double psi_y_neg_x(double **alpha_d, int i, int j, int is_it_boundary, int east_or_west)         //// function for y-momentum equation for finding flux limiter in x-direction
{
    double r_neg;

    if((east_or_west == 1 && is_it_boundary == 1) || i == NX-1 || ( alpha_d[i][j] - alpha_d[i-1][j] ) == 0)
    {
        if(( alpha_d[i][j] - alpha_d[i-1][j] ) == 0 )
            return 1;
        else
        {
            r_neg = 2*( - alpha_d[i][j] )/( alpha_d[i][j] - alpha_d[i-1][j] );
            return min_mod(r_neg);
        }
        return 0;
    }
    else
    {
        r_neg = ( alpha_d[i+1][j] - alpha_d[i][j] )/( alpha_d[i][j] - alpha_d[i-1][j] );
        return min_mod(r_neg);
    }
}

double psi_x_pos_x(double **alpha_d, int i, int j, int is_it_boundary, int east_or_west)            //// function for x-momentum equation for finding flux limiter in x-direction
{
    double r_pos;

    if(i == 0 || (alpha_d[i+1][j] - alpha_d[i][j]) == 0 )
    {
        if((alpha_d[i+1][j] - alpha_d[i][j]) == 0 )
            return 1;
        else
        {
            r_pos = ( alpha_d[i][j] + alpha_d[i+1][j] )/( alpha_d[i+1][j] - alpha_d[i][j] );
            return min_mod(r_pos);
        } 
        // return 1;
    }
    else
    {
        r_pos = ( alpha_d[i][j] - alpha_d[i-1][j] )/( alpha_d[i+1][j] - alpha_d[i][j] );
        return min_mod(r_pos);
    }
}


double psi_x_neg_x(double **alpha_d, int i, int j, int is_it_boundary, int east_or_west)                //// function for x-momentum equation for finding flux limiter in x-direction
{
    double r_neg;

    if(i == NX || ( alpha_d[i][j] - alpha_d[i-1][j] ) == 0)
    {
        if(( alpha_d[i][j] - alpha_d[i-1][j] ) == 0 )
            return 1;
        else
        {
            r_neg = ( - alpha_d[i-1][j] - alpha_d[i][j] )/( alpha_d[i][j] - alpha_d[i-1][j] );
            return min_mod(r_neg);
        } 

    } 
    else
    {
        r_neg = ( alpha_d[i+1][j] - alpha_d[i][j] )/( alpha_d[i][j] - alpha_d[i-1][j] );
        return min_mod(r_neg);
    }
}


double psi_x_pos_y(double **alpha_d, int i, int j, int is_it_boundary, int north_or_south, double *alphad_inlet)            //// function for x-momentum equation for finding flux limiter in y-direction
{
    double r_pos;

    if((north_or_south == 0 && is_it_boundary == 1) || j == 0 || ( alpha_d[i][j+1] - alpha_d[i][j] ) == 0)
    {
        if(( alpha_d[i][j+1] - alpha_d[i][j] ) == 0 )
            return 1;
        else
        {
            r_pos = 2*( alpha_d[i][j] - alphad_inlet[i] )/( alpha_d[i][j+1] - alpha_d[i][j] );
            return min_mod(r_pos);
        }
        // return 0;
    }
    else
    {
        r_pos = ( alpha_d[i][j] - alpha_d[i][j-1] )/( alpha_d[i][j+1] - alpha_d[i][j] );
        return min_mod(r_pos);
    }
}

double psi_x_neg_y(double **alpha_d, int i, int j, int is_it_boundary, int north_or_south)                  //// function for x-momentum equation for finding flux limiter in y-direction
{
    double r_neg;

    if((north_or_south == 1 && is_it_boundary == 1) || j == NY-1 || ( alpha_d[i][j] - alpha_d[i][j-1] ) == 0)
    {
        if(( alpha_d[i][j] - alpha_d[i][j-1] ) == 0 )
            return 1;
        else 
            return 0;
    }
    else
    {
        r_neg = ( alpha_d[i][j+1] - alpha_d[i][j] )/( alpha_d[i][j] - alpha_d[i][j-1] );
        return min_mod(r_neg);
    }
}


double psi_y_pos_y(double **alpha_d, int i, int j, int is_it_boundary, int north_or_south, double *alphad_inlet)                //// function for y-momentum equation for finding flux limiter in y-direction
{
    double r_pos;

    if(j == 0 || ( alpha_d[i][j+1] - alpha_d[i][j] ) == 0)
    {
        if(( alpha_d[i][j+1] - alpha_d[i][j] ) == 0 )
            return 1;
        else
        {
            r_pos = ( alpha_d[i][j] + alpha_d[i][j+1] )/( alpha_d[i][j+1] - alpha_d[i][j] );
            return min_mod(r_pos);
        }

    }
    else
    {
        r_pos = ( alpha_d[i][j] - alpha_d[i][j-1] )/( alpha_d[i][j+1] - alpha_d[i][j] );
        return min_mod(r_pos);
    }
}

double psi_y_neg_y(double **alpha_d, int i, int j, int is_it_boundary, int north_or_south)                  //// function for y-momentum equation for finding flux limiter in y-direction
{
    double r_neg;

    if(j == NY || ( alpha_d[i][j] - alpha_d[i][j-1] ) == 0)
    {
        if(( alpha_d[i][j] - alpha_d[i][j-1] ) == 0 )
            return 1;
        else
        {
            r_neg = ( - alpha_d[i][j-1] - alpha_d[i][j] )/( alpha_d[i][j] - alpha_d[i][j-1] );
            return min_mod(r_neg);
        }
        // return 1;
    }
    else
    {
        r_neg = ( alpha_d[i][j+1] - alpha_d[i][j] )/( alpha_d[i][j] - alpha_d[i][j-1] );
        return min_mod(r_neg);
    }
}

void ILU0(double E[], double F[], double D[], double B[], double H[], double f[], double d[], double b[], double e[], double c[], int K, int nx)
{
    //// Filling out L and U matrix

    e[0] = F[0];
    f[0] = H[0];
    d[0] = E[0];

    for(int l=1;l<K;l++)
    {
        e[l] = F[l];
        f[l] = H[l];
        c[l] = D[l]/d[l-1]; 

        if(l>=nx)
        {
            b[l] = B[l]/d[l-nx];
        }

        if(l>=nx)
        {
            d[l] = E[l] - (b[l]*f[l-nx] + c[l]*e[l-1]);
        }
        else
        {
            d[l] = E[l] - (c[l]*e[l-1]);
        }

    }
}

void MATRIX_VECTOR_MULTIPLY_Center(double *Ep, double *Dp, double *Fp, double *Hp, double *Bp, double *vec_delta_p, double *v, int nx)
{
    int i,j,m;

    for (i=1;i<NX-1;i++)
    {
        for (j=1;j<NY-1;j++)
        {
            m = j*nx + i;

            v[m] = Ep[m]*vec_delta_p[m] + Fp[m]*vec_delta_p[m+1] + Hp[m]*vec_delta_p[m+nx] + Dp[m]*vec_delta_p[m-1] + Bp[m]*vec_delta_p[m-nx];
        
        }
    }


    // Residuals in Bottom Cells

    j=0;

    for (i=1;i<NX-1;i++)
    {
        m = j*nx + i;

        v[m] = Ep[m]*vec_delta_p[m] + Fp[m]*vec_delta_p[m+1] + Hp[m]*vec_delta_p[m+nx] + Dp[m]*vec_delta_p[m-1];

    }
    
    // Residuals in Top cells

    j=NY-1;

    for (i=1;i<NX-1;i++)
    {
        m = j*nx + i;

        v[m] = Ep[m]*vec_delta_p[m] + Fp[m]*vec_delta_p[m+1] + Dp[m]*vec_delta_p[m-1] + Bp[m]*vec_delta_p[m-nx];

    }

    // Residual in Left cells

    i=0;

    for (j=1;j<NY-1;j++)
    {
        m = j*nx + i;

        v[m] = Ep[m]*vec_delta_p[m] + Fp[m]*vec_delta_p[m+1] + Hp[m]*vec_delta_p[m+nx] + Bp[m]*vec_delta_p[m-nx];

    }

    // Residual in Outlet Cells

    i=NX-1;

    for (j=1;j<NY-1;j++)
    {
        m = j*nx + i;

        v[m] = Ep[m]*vec_delta_p[m] + Dp[m]*vec_delta_p[m-1] + Hp[m]*vec_delta_p[m+nx] + Bp[m]*vec_delta_p[m-nx];

    }

    // Residual in Bottom-Left

    i=0;
    j=0;

    m = j*nx + i;

    v[m] = Ep[m]*vec_delta_p[m] + Fp[m]*vec_delta_p[m+1] + Hp[m]*vec_delta_p[m+nx];

    // Residual in Top-Left

    j=NY-1;

    m = j*nx + i;

    v[m] = Ep[m]*vec_delta_p[m] + Fp[m]*vec_delta_p[m+1] + Bp[m]*vec_delta_p[m-nx];

    // Residual in Bottom-Rivec_delta_pht

    i=NX-1;
    j=0;

    m = j*nx + i;

    v[m] = Ep[m]*vec_delta_p[m] + Dp[m]*vec_delta_p[m-1] + Hp[m]*vec_delta_p[m+nx];

    // Residual in Top-Rivec_delta_pht

    j=NY-1;

    m = j*nx + i;

    v[m] = Ep[m]*vec_delta_p[m] + Dp[m]*vec_delta_p[m-1] + Bp[m]*vec_delta_p[m-nx];


}


void MATRIX_VECTOR_MULTIPLY_X(double *Ep, double *Dp, double *Fp, double *Hp, double *Bp, double *vec_delta_p, double *v, int nx)
{
    int i,j,m;

    for (i=2;i<NX-1;i++)
    {
        for (j=1;j<NY-1;j++)
        {
            m = j*nx + i-1;

            v[m] = Ep[m]*vec_delta_p[m] + Fp[m]*vec_delta_p[m+1] + Hp[m]*vec_delta_p[m+nx] + Dp[m]*vec_delta_p[m-1] + Bp[m]*vec_delta_p[m-nx];
        
        }
    }


    // Residuals in Bottom Cells

    j=0;

    for (i=2;i<NX-1;i++)
    {
        m = j*nx + i-1;

        v[m] = Ep[m]*vec_delta_p[m] + Fp[m]*vec_delta_p[m+1] + Hp[m]*vec_delta_p[m+nx] + Dp[m]*vec_delta_p[m-1];

    }
    
    // Residuals in Top cells

    j=NY-1;

    for (i=2;i<NX-1;i++)
    {
        m = j*nx + i-1;

        v[m] = Ep[m]*vec_delta_p[m] + Fp[m]*vec_delta_p[m+1] + Dp[m]*vec_delta_p[m-1] + Bp[m]*vec_delta_p[m-nx];

    }

    // Residual in Left cells

    i=1;

    for (j=1;j<NY-1;j++)
    {
        m = j*nx + i-1;

        v[m] = Ep[m]*vec_delta_p[m] + Fp[m]*vec_delta_p[m+1] + Hp[m]*vec_delta_p[m+nx] + Bp[m]*vec_delta_p[m-nx];

    }

    // Residual in Outlet Cells

    i=NX-1;

    for (j=1;j<NY-1;j++)
    {
        m = j*nx + i-1;

        v[m] = Ep[m]*vec_delta_p[m] + Dp[m]*vec_delta_p[m-1] + Hp[m]*vec_delta_p[m+nx] + Bp[m]*vec_delta_p[m-nx];

    }

    // Residual in Bottom-Left

    i=1;
    j=0;

    m = j*nx + i-1;

    v[m] = Ep[m]*vec_delta_p[m] + Fp[m]*vec_delta_p[m+1] + Hp[m]*vec_delta_p[m+nx];

    // Residual in Top-Left

    j=NY-1;

    m = j*nx + i-1;

    v[m] = Ep[m]*vec_delta_p[m] + Fp[m]*vec_delta_p[m+1] + Bp[m]*vec_delta_p[m-nx];

    // Residual in Bottom-Rivec_delta_pht

    i=NX-1;
    j=0;

    m = j*nx + i-1;

    v[m] = Ep[m]*vec_delta_p[m] + Dp[m]*vec_delta_p[m-1] + Hp[m]*vec_delta_p[m+nx];

    // Residual in Top-Rivec_delta_pht

    j=NY-1;

    m = j*nx + i-1;

    v[m] = Ep[m]*vec_delta_p[m] + Dp[m]*vec_delta_p[m-1] + Bp[m]*vec_delta_p[m-nx];


}


void MATRIX_VECTOR_MULTIPLY_Y(double *Ep, double *Dp, double *Fp, double *Hp, double *Bp, double *vec_delta_p, double *v, int nx)
{
    int i,j,m;

    for (i=1;i<NX-1;i++)
    {
        for (j=2;j<NY-1;j++)
        {
            m = (j-1)*nx + i;

            v[m] = Ep[m]*vec_delta_p[m] + Fp[m]*vec_delta_p[m+1] + Hp[m]*vec_delta_p[m+nx] + Dp[m]*vec_delta_p[m-1] + Bp[m]*vec_delta_p[m-nx];
        
        }
    }


    // Residuals in Bottom Cells

    j=1;

    for (i=1;i<NX-1;i++)
    {
        m = (j-1)*nx + i;

        v[m] = Ep[m]*vec_delta_p[m] + Fp[m]*vec_delta_p[m+1] + Hp[m]*vec_delta_p[m+nx] + Dp[m]*vec_delta_p[m-1];

    }
    
    // Residuals in Top cells

    j=NY-1;

    for (i=1;i<NX-1;i++)
    {
        m = (j-1)*nx + i;

        v[m] = Ep[m]*vec_delta_p[m] + Fp[m]*vec_delta_p[m+1] + Dp[m]*vec_delta_p[m-1] + Bp[m]*vec_delta_p[m-nx];

    }

    // Residual in Left cells

    i=0;

    for (j=2;j<NY-1;j++)
    {
        m = (j-1)*nx + i;

        v[m] = Ep[m]*vec_delta_p[m] + Fp[m]*vec_delta_p[m+1] + Hp[m]*vec_delta_p[m+nx] + Bp[m]*vec_delta_p[m-nx];

    }

    // Residual in Outlet Cells

    i=NX-1;

    for (j=2;j<NY-1;j++)
    {
        m = (j-1)*nx + i;

        v[m] = Ep[m]*vec_delta_p[m] + Dp[m]*vec_delta_p[m-1] + Hp[m]*vec_delta_p[m+nx] + Bp[m]*vec_delta_p[m-nx];

    }

    // Residual in Bottom-Left

    i=0;
    j=1;

    m = (j-1)*nx + i;

    v[m] = Ep[m]*vec_delta_p[m] + Fp[m]*vec_delta_p[m+1] + Hp[m]*vec_delta_p[m+nx];

    // Residual in Top-Left

    j=NY-1;

    m = (j-1)*nx + i;

    v[m] = Ep[m]*vec_delta_p[m] + Fp[m]*vec_delta_p[m+1] + Bp[m]*vec_delta_p[m-nx];

    // Residual in Bottom-Rivec_delta_pht

    i=NX-1;
    j=1;

    m = (j-1)*nx + i;

    v[m] = Ep[m]*vec_delta_p[m] + Dp[m]*vec_delta_p[m-1] + Hp[m]*vec_delta_p[m+nx];

    // Residual in Top-Rivec_delta_pht

    j=NY-1;

    m = (j-1)*nx + i;

    v[m] = Ep[m]*vec_delta_p[m] + Dp[m]*vec_delta_p[m-1] + Bp[m]*vec_delta_p[m-nx];

}

void residual_center(double *Q, double *E, double *F, double *H, double *D, double *B, double *res, double *variable, int nx)
{
    int i,j,m;

    // Residuals in Interior Cells

    for (i=1;i<NX-1;i++)
    {
        for (j=1;j<NY-1;j++)
        {
            m = j*nx + i;

            res[m] = Q[m] - E[m]*variable[m] - F[m]*variable[m+1] - H[m]*variable[m+nx] - D[m]*variable[m-1] - B[m]*variable[m-nx];
        
        }
    }


    // Residuals in Bottom Cells

    j=0;

    for (i=1;i<NX-1;i++)
    {
        m = j*nx + i;

        res[m] = Q[m] - E[m]*variable[m] - F[m]*variable[m+1] - H[m]*variable[m+nx] - D[m]*variable[m-1];

    }
    
    // Residuals in Toa cells

    j=NY-1;

    for (i=1;i<NX-1;i++)
    {
        m = j*nx + i;

        res[m] = Q[m] - E[m]*variable[m] - F[m]*variable[m+1] - D[m]*variable[m-1] - B[m]*variable[m-nx];

    }

    // Residual in Left cells

    i=0;

    for (j=1;j<NY-1;j++)
    {
        m = j*nx + i;

        res[m] = Q[m] - E[m]*variable[m] - F[m]*variable[m+1] - H[m]*variable[m+nx] - B[m]*variable[m-nx];

    }

    // Residual in Outlet Cells

    i=NX-1;

    for (j=1;j<NY-1;j++)
    {
        m = j*nx + i;

        res[m] = Q[m] - E[m]*variable[m] - D[m]*variable[m-1] - H[m]*variable[m+nx] - B[m]*variable[m-nx];

    }

    // Residual in Bottom-Left

    i=0;
    j=0;

    m = j*nx + i;

    res[m] = Q[m] - E[m]*variable[m] - F[m]*variable[m+1] - H[m]*variable[m+nx];

    // Residual in Toa-Left

    j=NY-1;

    m = j*nx + i;

    res[m] = Q[m] - E[m]*variable[m] - F[m]*variable[m+1] - B[m]*variable[m-nx];

    // Residual in Bottom-Right

    i=NX-1;
    j=0;

    m = j*nx + i;

    res[m] = Q[m] - E[m]*variable[m] - D[m]*variable[m-1] - H[m]*variable[m+nx];

    // Residual in Toa-Right

    j=NY-1;

    m = j*nx + i;

    res[m] = Q[m] - E[m]*variable[m] - D[m]*variable[m-1] - B[m]*variable[m-nx];

}

void residual_X(double *Q, double *E, double *F, double *H, double *D, double *B, double *res, double *variable, double **original, int nx)
{
    int i,j,m;

        // Residual in Interior Cells

        for (i=2;i<NX-1;i++)
        {
            for (j=1;j<NY-1;j++)
            {
                m = j*nx + i-1;

                res[m] = Q[m] - E[m]*variable[m] - F[m]*variable[m+1] - H[m]*variable[m+nx] - D[m]*variable[m-1] - B[m]*variable[m-nx];
            
            }
        }

    // Residual in Left Boundary cells

    i=1;

    for (j=1;j<NY-1;j++)
    {
        m = j*nx + i-1;

        res[m] = Q[m] - E[m]*variable[m] - F[m]*variable[m+1] - H[m]*variable[m+nx] - B[m]*variable[m-nx] - D[m]*original[i-1][j];
    
    }
    // Residual in Right boundary Cells

    i=NX-1;

    for (j=1;j<NY-1;j++)
    {
        m = j*nx + i-1;

        res[m] = Q[m] - E[m]*variable[m] - H[m]*variable[m+nx] - D[m]*variable[m-1] - B[m]*variable[m-nx] - F[m]*original[i+1][j];
    
    }

    // Residual in Bottom Boundary Cells

    j=0;

    for (i=2;i<NX-1;i++)
    {
        m = j*nx + i-1;

        res[m] = Q[m] - E[m]*variable[m] - F[m]*variable[m+1] - H[m]*variable[m+nx] - D[m]*variable[m-1];

    }

    // Residual in Top Boundary Cells

    j=NY-1;

    for (i=2;i<NX-1;i++)
    {
        m = j*nx + i-1;

        res[m] = Q[m] - E[m]*variable[m] - F[m]*variable[m+1] - D[m]*variable[m-1] - B[m]*variable[m-nx];

    }

    // Residual in Left-Bottom cell

    i=1;
    j=0;

    m = j*nx + i-1;

    res[m] = Q[m] - E[m]*variable[m] - F[m]*variable[m+1] - H[m]*variable[m+nx] - D[m]*original[i-1][j];

    // Residual in Top-Left cell

    j=NY-1;

    m = j*nx + i-1;

    res[m] = Q[m] - E[m]*variable[m] - F[m]*variable[m+1] - B[m]*variable[m-nx] - D[m]*original[i-1][j];

    // Residual in Bottom-Right cell

    i=NX-1;
    j=0;

    m = j*nx + i-1;

    res[m] = Q[m] - E[m]*variable[m] - H[m]*variable[m+nx] - D[m]*variable[m-1] - F[m]*original[i+1][j];

    // Residual in Top-Right Cell

    j=NY-1;

    m = j*nx + i-1;

    res[m] = Q[m] - E[m]*variable[m] - D[m]*variable[m-1] - B[m]*variable[m-nx] - F[m]*original[i+1][j];

}


void residual_Y(double *Q, double *E, double *F, double *H, double *D, double *B, double *res, double *variable, double **original, int nx)
{
    int i,j,m;

        // Residuals in Interior Cells

        for (i=1;i<NX-1;i++)
        {
            for (j=2;j<NY-1;j++)
            {
                m = (j-1)*NX + i;

                res[m] = Q[m] - E[m]*variable[m] - F[m]*variable[m+1] - H[m]*variable[m+nx] - D[m]*variable[m-1] - B[m]*variable[m-nx];
            
            }
        }


        // Residuals in Bottom Cells

        j=1;

        for (i=1;i<NX-1;i++)
        {
            m = (j-1)*NX + i;

            res[m] = Q[m] - E[m]*variable[m] - F[m]*variable[m+1] - H[m]*variable[m+nx] - D[m]*variable[m-1] - B[m]*original[i][j-1];

        }
    
        // Residuals in Top cells

        j=NY-1;

        for (i=1;i<NX-1;i++)
        {
            m = (j-1)*NX + i;

            res[m] = Q[m] - E[m]*variable[m] - F[m]*variable[m+1] - D[m]*variable[m-1] - B[m]*variable[m-nx] - H[m]*original[i][j+1];

        }

        // Residual in Left cells

        i=0;

        for (j=2;j<NY-1;j++)
        {
            m = (j-1)*NX + i;

            res[m] = Q[m] - E[m]*variable[m] - F[m]*variable[m+1] - H[m]*variable[m+nx] - B[m]*variable[m-nx];

        }

        // Residual in Outlet Cells

        i=NX-1;

        for (j=2;j<NY-1;j++)
        {
            m = (j-1)*NX + i;

            res[m] = Q[m] - E[m]*variable[m] - D[m]*variable[m-1] - H[m]*variable[m+nx] - B[m]*variable[m-nx];

        }

        // Residual in Bottom-Left

        i=0;
        j=1;

        m = (j-1)*NX + i;

        res[m] = Q[m] - E[m]*variable[m] - F[m]*variable[m+1] - H[m]*variable[m+nx] - B[m]*original[i][j-1];

        // Residual in Top-Left

        j=NY-1;

        m = (j-1)*NX + i;

        res[m] = Q[m] - E[m]*variable[m] - F[m]*variable[m+1] - B[m]*variable[m-nx] - H[m]*original[i][j+1];

        // Residual in Bottom-Right

        i=NX-1;
        j=1;

        m = (j-1)*NX + i;

        res[m] = Q[m] - E[m]*variable[m] - D[m]*variable[m-1] - H[m]*variable[m+nx] - B[m]*original[i][j-1];

        // Residual in Top-Right

        j=NY-1;

        m = (j-1)*NX + i;

        res[m] = Q[m] - E[m]*variable[m] - D[m]*variable[m-1] - B[m]*variable[m-nx] - H[m]*original[i][j+1];

}
