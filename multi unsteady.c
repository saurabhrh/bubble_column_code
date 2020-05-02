#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <time.h>
#include "functions.h"

#define NX 25
#define NY 75
#define DX 2e-2
#define DY 2e-2
#define dt 5e-2
#define OUTER 20
#define ITERATION_X 5
#define ITERATION_Y 5
#define ITERATION_P 50 
#define ITERATION_KEPS 5
#define ITERATION_ALPHA 5
#define TIMESTEPS 500
#define SAVING_INTERVAL 600
#define C_MU 0.09
#define C1 1.44
#define C2 1.92
#define K_COR 0.3
#define EPS_COR 0.3
#define ALPHA_COR 0.3
#define RHO_C 998.2
#define RHO_D 1.225
#define UV_COR 0.7
#define P_COR 0.3
#define KAPPA 0.41
#define EE 9.8
#define SIGMAK 1.0
#define SIGMAEPS 1.3
#define V_INLET 0.2
#define TOL_X_C 1e-6
#define TOL_Y_C 1e-6
#define TOL_X_D 1e-6
#define TOL_Y_D 1e-6
#define TOL_KEPS 1e-6
#define TOL_ALPHA 1e-16
#define TOLERANCE 1e-5
#define RELATIVE_TOL 0.0001
#define DP 3e-3
#define g -9.81
#define SURF_TENS 0.07
#define mu_dispersed 18.03e-6
#define mu 1e-3
#define ALPHA_XY 0.92
#define ALPHA_P 0.92
#define ALPHA_KEPS 0.92
#define I 0.05
#define max(A, B) ( (A) > (B) ? (A) : (B) )

int main()
{
    int i,j,nx=NX-1,ny=NY-1,inner,outer,Kx=nx*NY,Ky=NX*ny,Kp=NX*NY,m,time,numsaves=TIMESTEPS/SAVING_INTERVAL,count=0,z_e,z_w,z_n,z_s;
    FILE *resuc;
    FILE *resvc;
    FILE *resud;
    FILE *resvd;
    FILE *resP;
    FILE *resalpha;
    FILE *resK;
    FILE *resEPS;
    FILE *uc[numsaves];
    FILE *vc[numsaves];
    FILE *ud[numsaves];
    FILE *vd[numsaves];
    FILE *Pr[numsaves];
    FILE *Kfile[numsaves];
    FILE *EpsFile[numsaves];
    FILE *mutC[numsaves];
    FILE *mutD[numsaves];
    FILE *alphac[numsaves];
    FILE *alphad[numsaves];
    char statement_xc[30] = "Residual_x_c at iteration";
    char statement_yc[30] = "Residual_y_c at iteration";
    char statement_xd[30] = "Residual_x_d at iteration";
    char statement_yd[30] = "Residual_y_d at iteration";
    char statement_P[30] = "Residual_P at iteration";
    char statement_k[30] = "Residual_k at iteration";
    char statement_eps[30] = "Residual_Eps at iteration";
    char statement_alpha[30] = "Residual_alpha at iteration";
    double DX2=DX*DX,DY2=DY*DY,ar=DY/DX,ra=1/ar,U_x,V_y,U_y,V_x,top,left,right,bottom,AR=ar/SIGMAK,RA=ra/SIGMAK,R=ra/SIGMAEPS,r=ar/SIGMAEPS,alphacn_centre,alphadn_centre,xp=DX/2,up,vp,tau_y,mu_w,mu_s,mu_e,mu_n,v_plus,x_plus,kp,fE,fW,fN,fS,mut_e,mut_n,mut_s,mut_w,nor,residual_y_c,residual_y_d,residual_P,uE,uW,vN,vS,residual_eps,residual_k,residual_x_c,residual_x_d,residual_alpha,alc_e,alc_w,alc_n,alc_s,E,Eo,Re_b,C_D,alpe_c,alpn_c,alps_c,alpw_c,alpe_d,alpw_d,alpn_d,alps_d,ald_e,ald_w,ald_n,ald_s,U_mc,U_md,U_diffmag,V_md,V_mc,mdot_c,mdot_d,a_if,R_1,R2eps,R2p,R2k,R2a,alphad_centre,alphac_centre,Edash,C_D_Churn,C_D_Distorted,C_D_Newton,U_diffmagx,U_diffmagy,kappa = KAPPA*xp,DT=1/dt,cpu_time;
    clock_t start,end;

    double ** u_c = calloc((NX+1),sizeof(double *));
    double ** u_d = calloc((NX+1),sizeof(double *));
    double ** B_kx = calloc((NX+1),sizeof(double *));
    double ** DIV_xc = calloc((NX+1),sizeof(double *));
    double ** DIV_xd = calloc((NX+1),sizeof(double *));
    double ** un_d = calloc((NX+1),sizeof(double *));
    double ** un_c = calloc((NX+1),sizeof(double *));
    double ** F_added_x = calloc((NX+1),sizeof(double *));
    double ** mat_der_x = calloc((NX+1),sizeof(double *));
    for (i=0;i<NX+1;i++)
    {
        u_c[i] = calloc(NY,sizeof(double));
        u_d[i] = calloc(NY,sizeof(double));
        B_kx[i] = calloc(NY,sizeof(double));
        DIV_xc[i] = calloc(NY,sizeof(double));
        DIV_xd[i] = calloc(NY,sizeof(double));
        un_d[i] = calloc(NY,sizeof(double));
        un_c[i] = calloc(NY,sizeof(double));
        F_added_x[i] = calloc(NY,sizeof(double));
        mat_der_x[i] = calloc(NY,sizeof(double));
    }
      
    double ** v_c = calloc(NX,sizeof(double *));
    double ** v_d = calloc(NX,sizeof(double *));
    double ** B_ky = calloc(NX,sizeof(double *));
    double ** DIV_yc = calloc(NX,sizeof(double *));
    double ** DIV_yd = calloc(NX,sizeof(double *));
    double ** vn_d = calloc(NX,sizeof(double *));
    double ** vn_c = calloc(NX,sizeof(double *));
    double ** F_added_y = calloc(NX,sizeof(double *));
    double ** mat_der_y = calloc(NX,sizeof(double *));
    for (i=0;i<NX;i++)
    {
        v_c[i] = calloc((NY+1),sizeof(double));
        v_d[i] = calloc((NY+1),sizeof(double));
        B_ky[i] = calloc((NY+1),sizeof(double));
        DIV_yc[i] = calloc((NY+1),sizeof(double));
        DIV_yd[i] = calloc((NY+1),sizeof(double));
        vn_d[i] = calloc((NY+1),sizeof(double));
        vn_c[i] = calloc((NY+1),sizeof(double));
        F_added_y[i] = calloc((NY+1),sizeof(double));
        mat_der_y[i] = calloc((NY+1),sizeof(double));
    }

    double ** alpha_d = calloc(NX,sizeof(double *));
    double ** alpha_c = calloc(NX,sizeof(double *));
    double ** P = calloc(NX,sizeof(double *));
    double ** mut_C = calloc(NX,sizeof(double *));
    double ** k = calloc(NX,sizeof(double *));
    double ** epsilon = calloc(NX,sizeof(double *));
    double ** mu_total_C = calloc(NX,sizeof(double *));
    double ** mu_total_D = calloc(NX,sizeof(double *));
    double ** uslip_mag = calloc(NX,sizeof(double *));
    double ** DIV = calloc(NX,sizeof(double *));
    double ** source = calloc(NX,sizeof(double *));
    double ** production = calloc(NX,sizeof(double *));
    double ** kn = calloc(NX,sizeof(double *));
    double ** epsn = calloc(NX,sizeof(double *));
    double ** alphan_d = calloc(NX,sizeof(double *));
    double ** alphan_c = calloc(NX,sizeof(double *));
    for (i=0;i<NX;i++)
    {
        alpha_d[i] = calloc(NY,sizeof(double));
        alpha_c[i] = calloc(NY,sizeof(double));
        P[i] = calloc(NY,sizeof(double));
        mut_C[i] = calloc(NY,sizeof(double));
        k[i] = calloc(NY,sizeof(double));
        epsilon[i] = calloc(NY,sizeof(double));
        mu_total_C[i] = calloc(NY,sizeof(double));
        mu_total_D[i] = calloc(NY,sizeof(double));
        uslip_mag[i] = calloc(NY,sizeof(double));
        DIV[i] = calloc(NY,sizeof(double));
        source[i] = calloc(NY,sizeof(double));
        production[i] = calloc(NY,sizeof(double));
        kn[i] = calloc(NY,sizeof(double));
        epsn[i] = calloc(NY,sizeof(double));
        alphan_d[i] = calloc(NY,sizeof(double));
        alphan_c[i] = calloc(NY,sizeof(double));
    }

    double * k_inlet = calloc(NX,sizeof(double));
    double * v_inlet = calloc(NX,sizeof(double));
    double * eps_inlet = calloc(NX,sizeof(double));
    double * mut_c_inlet = calloc(NX,sizeof(double));
    double * mu_total_inletc = calloc(NX,sizeof(double));
    double * mu_total_inletd = calloc(NX,sizeof(double));
    double * muxc_inlet = calloc(NX,sizeof(double));
    double * muxd_inlet = calloc(NX,sizeof(double));
    double * alphad_inlet = calloc(NX,sizeof(double));
    double * alphac_inlet = calloc(NX,sizeof(double));
    double * u_cinlet = calloc(NX+1,sizeof(double));
    double * u_dinlet = calloc(NX+1,sizeof(double));
    double * Uc = calloc(Kx,sizeof(double));
    double * Vc = calloc(Ky,sizeof(double));
    double * Ud = calloc(Kx,sizeof(double));
    double * Vd = calloc(Ky,sizeof(double));
    double * p = calloc(Kp,sizeof(double));
    double * K = calloc(Kp,sizeof(double));
    double * Eps = calloc(Kp,sizeof(double));
    double * a = calloc(Kp,sizeof(double));

    double * Exd = calloc(Kx,sizeof(double));
    double * Fxd = calloc(Kx,sizeof(double));
    double * Dxd = calloc(Kx,sizeof(double));
    double * Bxd = calloc(Kx,sizeof(double));
    double * Hxd = calloc(Kx,sizeof(double));
    double * Qxd = calloc(Kx,sizeof(double));
    double * Qxdmod = calloc(Kx,sizeof(double));
    double * fxd = calloc(Kx,sizeof(double));
    double * dxd = calloc(Kx,sizeof(double));
    double * exd = calloc(Kx,sizeof(double));
    double * bxd = calloc(Kx,sizeof(double));
    double * cxd = calloc(Kx,sizeof(double));
    double * resxd = calloc(Kx,sizeof(double));
    double * vec_Yxd = calloc(Kx,sizeof(double));
    double * vec_deltaxd = calloc(Kx,sizeof(double));

    double * Exc = calloc(Kx,sizeof(double));
    double * Fxc = calloc(Kx,sizeof(double));
    double * Dxc = calloc(Kx,sizeof(double));
    double * Bxc = calloc(Kx,sizeof(double));
    double * Hxc = calloc(Kx,sizeof(double));
    double * Qxc = calloc(Kx,sizeof(double));
    double * Qxcmod = calloc(Kx,sizeof(double));
    double * fxc = calloc(Kx,sizeof(double));
    double * dxc = calloc(Kx,sizeof(double));
    double * exc = calloc(Kx,sizeof(double));
    double * bxc = calloc(Kx,sizeof(double));
    double * cxc = calloc(Kx,sizeof(double));
    double * resxc = calloc(Kx,sizeof(double));
    double * vec_Yxc = calloc(Kx,sizeof(double));
    double * vec_deltaxc = calloc(Kx,sizeof(double));

    double * Eyc = calloc(Ky,sizeof(double));
    double * Fyc = calloc(Ky,sizeof(double));
    double * Dyc = calloc(Ky,sizeof(double));
    double * Byc = calloc(Ky,sizeof(double));
    double * Hyc = calloc(Ky,sizeof(double));
    double * Qyc = calloc(Ky,sizeof(double));
    double * Qycmod = calloc(Ky,sizeof(double));
    double * fyc = calloc(Ky,sizeof(double));
    double * dyc = calloc(Ky,sizeof(double));
    double * eyc = calloc(Ky,sizeof(double));
    double * byc = calloc(Ky,sizeof(double));
    double * cyc = calloc(Ky,sizeof(double));
    double * resyc = calloc(Ky,sizeof(double));
    double * vec_Yyc = calloc(Ky,sizeof(double));
    double * vec_deltayc = calloc(Ky,sizeof(double));

    double * Eyd = calloc(Ky,sizeof(double));
    double * Fyd = calloc(Ky,sizeof(double));
    double * Dyd = calloc(Ky,sizeof(double));
    double * Byd = calloc(Ky,sizeof(double));
    double * Hyd = calloc(Ky,sizeof(double));
    double * Qyd = calloc(Ky,sizeof(double));
    double * Qydmod = calloc(Ky,sizeof(double));
    double * fyd = calloc(Ky,sizeof(double));
    double * dyd = calloc(Ky,sizeof(double));
    double * eyd = calloc(Ky,sizeof(double));
    double * byd = calloc(Ky,sizeof(double));
    double * cyd = calloc(Ky,sizeof(double));
    double * resyd = calloc(Ky,sizeof(double));
    double * vec_Yyd = calloc(Ky,sizeof(double));
    double * vec_deltayd = calloc(Ky,sizeof(double));

    double * Ep = calloc(Kp,sizeof(double));
    double * Fp = calloc(Kp,sizeof(double));
    double * Dp = calloc(Kp,sizeof(double));
    double * Bp = calloc(Kp,sizeof(double));
    double * Hp = calloc(Kp,sizeof(double));
    double * Qp = calloc(Kp,sizeof(double));
    double * fp = calloc(Kp,sizeof(double));
    double * dp = calloc(Kp,sizeof(double));
    double * ep = calloc(Kp,sizeof(double));
    double * bp = calloc(Kp,sizeof(double));
    double * cp = calloc(Kp,sizeof(double));
    double * resp = calloc(Kp,sizeof(double));
    double * vec_Yp = calloc(Kp,sizeof(double));
    double * vec_deltap = calloc(Kp,sizeof(double));

    double * Ek = calloc(Kp,sizeof(double));
    double * Fk = calloc(Kp,sizeof(double));
    double * Dk = calloc(Kp,sizeof(double));
    double * Bk = calloc(Kp,sizeof(double));
    double * Hk = calloc(Kp,sizeof(double));
    double * Qk = calloc(Kp,sizeof(double));
    double * fk = calloc(Kp,sizeof(double));
    double * dk = calloc(Kp,sizeof(double));
    double * ek = calloc(Kp,sizeof(double));
    double * bk = calloc(Kp,sizeof(double));
    double * ck = calloc(Kp,sizeof(double));
    double * resk = calloc(Kp,sizeof(double));
    double * vec_Yk = calloc(Kp,sizeof(double));
    double * vec_deltak = calloc(Kp,sizeof(double));

    double * Ee = calloc(Kp,sizeof(double));
    double * Fe = calloc(Kp,sizeof(double));
    double * De = calloc(Kp,sizeof(double));
    double * Be = calloc(Kp,sizeof(double));
    double * He = calloc(Kp,sizeof(double));
    double * Qe = calloc(Kp,sizeof(double));
    double * fe = calloc(Kp,sizeof(double));
    double * de = calloc(Kp,sizeof(double));
    double * ee = calloc(Kp,sizeof(double));
    double * be = calloc(Kp,sizeof(double));
    double * ce = calloc(Kp,sizeof(double));
    double * reseps = calloc(Kp,sizeof(double));
    double * vec_Yeps = calloc(Kp,sizeof(double));
    double * vec_deltaeps = calloc(Kp,sizeof(double));

    double * Ea = calloc(Kp,sizeof(double));
    double * Fa = calloc(Kp,sizeof(double));
    double * Da = calloc(Kp,sizeof(double));
    double * Ba = calloc(Kp,sizeof(double));
    double * Ha = calloc(Kp,sizeof(double));
    double * Qa = calloc(Kp,sizeof(double));
    double * fa = calloc(Kp,sizeof(double));
    double * da = calloc(Kp,sizeof(double));
    double * ea = calloc(Kp,sizeof(double));
    double * ba = calloc(Kp,sizeof(double));
    double * ca = calloc(Kp,sizeof(double));
    double * resa = calloc(Kp,sizeof(double));
    double * vec_Ya = calloc(Kp,sizeof(double));
    double * vec_deltaa = calloc(Kp,sizeof(double));

    // double dt = 1e-2;
    // double DT = 1/dt;

//// Initialisation and Boundary Conditions

    for (i=0;i<NX;i++)
    {
        for (j=0;j<NY;j++)
        {
            k[i][j] = 1.5*V_INLET*V_INLET*I*I;
            epsilon[i][j] = pow(C_MU,0.75)*pow(k[i][j],1.5)/(0.07*DP);
            mut_C[i][j] = RHO_C*C_MU*pow(k[i][j],2)/epsilon[i][j];
            mu_total_C[i][j] = mu + mut_C[i][j];
            mu_total_D[i][j] = mu_total_C[i][j]*RHO_D/RHO_C;
            alpha_c[i][j] = 1 - 1e-10;
            alpha_d[i][j] = 1e-10;
            // P[i][j] = 1e5;
        }
    }    

    for (i=1;i<NX;i++)
    {
        for (j=0;j<NY;j++)
        {
            u_c[i][j] = 0;
            u_d[i][j] = 0;
        }
    }

    for (i=0;i<NX;i++)
    {
        for (j=1;j<NY;j++)
        {
            v_c[i][j] = 0;
            v_d[i][j] = 1e-100;
        }
    }

    // Inlet boundary conditions

    for (i=12;i<24;i++)
    {
        v_d[i][0] = V_INLET;
        v_inlet[i] = v_d[i][0];
        k_inlet[i] = 1.5*V_INLET*V_INLET*I*I;
        eps_inlet[i] = pow(C_MU,0.75)*pow(k_inlet[i],1.5)/(0.07*DP);
        mut_c_inlet[i] = RHO_C*C_MU*pow(k_inlet[i],2)/eps_inlet[i];
        mu_total_inletc[i] = mu + mut_c_inlet[i];
        mu_total_inletd[i] = mu_total_inletc[i]*RHO_D/RHO_C;
        alphad_inlet[i] = 0.0417;
    }

    for (i=1;i<NX;i++)
    {
        muxc_inlet[i] = ( mu_total_inletc[i] + mu_total_inletc[i-1] )/2;
        muxd_inlet[i] = ( mu_total_inletd[i] + mu_total_inletd[i-1] )/2;
    }

    for (i=0;i<NX;i++)
    {
        alphac_inlet[i] = 1 - alphad_inlet[i];
    }

    /// Inlet boundary conditions end


    for (i=0;i<NX+1;i++)
    {
        for (j=0;j<NY;j++)
        {
            un_d[i][j] = u_d[i][j];
            un_c[i][j] = u_c[i][j];
        }
    }
    
    for (i=0;i<NX;i++)
    {
        for (j=0;j<NY+1;j++)
        {
            vn_d[i][j] = v_d[i][j];
            vn_c[i][j] = v_c[i][j];
        }
    }
        
    for (i=0;i<NX;i++)
    {
        for (j=0;j<NY;j++)
        {
            alphan_d[i][j] = alpha_d[i][j];
            alphan_c[i][j] = alpha_c[i][j];
            kn[i][j] = k[i][j];
            epsn[i][j] = epsilon[i][j];
        }
    }
    

for (time=1;time<=TIMESTEPS;time++)
{
    // if (time>2180)
    // {
    //     dt = 1e-5;
    //     DT = 1/dt;
    // }
    // start = clock();

    // Outer Iteration Loop

    for (outer=1;outer<=OUTER;outer++)
    {
        double ** uc_ap = calloc(NX,sizeof(double *));
        double ** uc_an = calloc(NX,sizeof(double *));
        double ** uc_as = calloc(NX,sizeof(double *));
        double ** uc_ae = calloc(NX,sizeof(double *));
        double ** uc_aw = calloc(NX,sizeof(double *));
        for (i=0;i<NX;i++)
        {
            uc_ap[i] = calloc(NY,sizeof(double));
            uc_an[i] = calloc(NY,sizeof(double));
            uc_as[i] = calloc(NY,sizeof(double));
            uc_ae[i] = calloc(NY,sizeof(double));
            uc_aw[i] = calloc(NY,sizeof(double));
        }
                
        // x-momentum links
        

        for (i=1;i<NX;i++)
        {
            for (j=1;j<NY-1;j++)
            {
            	V_mc = (v_c[i][j] + v_c[i][j+1] + v_c[i-1][j] + v_c[i-1][j+1])/4;
            	V_md = (v_d[i][j] + v_d[i][j+1] + v_d[i-1][j] + v_d[i-1][j+1])/4;

            	mat_der_x[i][j] = u_c[i][j]*( u_c[i+1][j] - u_c[i-1][j] )/(2*DX) +  V_mc*( u_c[i][j+1] - u_c[i][j-1] )/(2*DX) - u_d[i][j]*( u_d[i+1][j] - u_d[i-1][j] )/(2*DX) - V_md*( u_d[i][j+1] - u_d[i][j-1] )/(2*DX);
            }   
        }

        // Bottom Boundary
        
        j = 0;
        
        for (i=1;i<NX;i++)
        {
        	V_mc = (v_c[i][j] + v_c[i][j+1] + v_c[i-1][j] + v_c[i-1][j+1])/4;
        	V_md = (v_d[i][j] + v_d[i][j+1] + v_d[i-1][j] + v_d[i-1][j+1])/4;

        	mat_der_x[i][j] = u_c[i][j]*( u_c[i+1][j] - u_c[i-1][j] )/(2*DX) +  V_mc*( (u_c[i][j+1] + u_c[i][j])/2 )/(2*DX) - u_d[i][j]*( u_d[i+1][j] + u_d[i-1][j] )/(2*DX) - V_md*( (u_d[i][j+1] + u_d[i][j])/2 )/(2*DX);

        }
                
        // Outlet
        
        j = NY-1;
        
        for (i=1;i<NX;i++)
        {
        	V_mc = (v_c[i][j] + v_c[i][j+1] + v_c[i-1][j] + v_c[i-1][j+1])/4;
        	V_md = (v_d[i][j] + v_d[i][j+1] + v_d[i-1][j] + v_d[i-1][j+1])/4;

        	mat_der_x[i][j] = u_c[i][j]*( u_c[i+1][j] - u_c[i-1][j] )/(2*DX) +  V_mc*( u_c[i][j] - (u_c[i][j] + u_c[i][j-1])/2 )/(2*DX) - u_d[i][j]*( u_d[i+1][j] - u_d[i-1][j] )/(2*DX) - V_md*( u_d[i][j] - (u_d[i][j] + u_d[i][j-1])/2 )/(2*DX);
        }


        for (i=1;i<NX;i++)
        {
            for (j=0;j<NY;j++)
            {
                alphad_centre = ( alpha_d[i-1][j] + alpha_d[i][j] )/2;
                alphac_centre = ( alpha_c[i-1][j] + alpha_c[i][j] )/2;

                U_mc = u_c[i][j];
                V_mc = (v_c[i][j] + v_c[i][j+1] + v_c[i-1][j] + v_c[i-1][j+1])/4;

                U_md = u_d[i][j];
                V_md = (v_d[i][j] + v_d[i][j+1] + v_d[i-1][j] + v_d[i-1][j+1])/4;

                U_diffmagx = sqrt( pow((U_md - U_mc),2) + pow((V_md - V_mc),2) );

                Re_b = RHO_C*U_diffmagx*DP/mu;

                if (Re_b<1000)
                {
                    C_D = 24*( 1 + 0.15*pow(Re_b,0.687) )/Re_b;
                }
                else
                {
                    C_D = 0.44;
                }

                a_if = 6*alphad_centre/DP;

                B_kx[i][j] = C_D*a_if*RHO_C*U_diffmagx*DX2/8;

                DIV_xc[i][j] = ( u_c[i+1][j] - u_c[i-1][j] )/(2*DX) + ( v_c[i][j+1] + v_c[i-1][j+1] - v_c[i-1][j] - v_c[i][j] )/(2*DY);

                DIV_xd[i][j] = ( u_d[i+1][j] - u_d[i-1][j] )/(2*DX) + ( v_d[i][j+1] + v_d[i-1][j+1] - v_d[i-1][j] - v_d[i][j] )/(2*DY);

                F_added_x[i][j] = 0.5*RHO_C*alphad_centre*( ((u_c[i][j] - u_d[i][j]) - (un_c[i][j] - un_d[i][j]))*DT + mat_der_x[i][j])*DX2;// + u_c[i][j]*( u_c[i+1][j] - u_c[i-1][j] )/(2*DX) +  V_mc*( u_c[i][j+1] - u_c[i][j-1] )/(2*DX) - u_d[i][j]*( u_d[i+1][j] - u_d[i-1][j] )/(2*DX) +  V_md*( u_d[i][j+1] - u_d[i][j-1] )/(2*DX) )*DX2;

            }
        }

        /// Interior Cells
        
        for (i=1;i<NX;i++)
        {
            for (j=1;j<NY-1;j++)
            {
                uE = DY*RHO_C*(u_c[i][j] + u_c[i+1][j])/2;
                uW = DY*RHO_C*(u_c[i][j] + u_c[i-1][j])/2;
                vS = DX*RHO_C*(v_c[i][j] + v_c[i-1][j])/2;
                vN = DX*RHO_C*(v_c[i][j+1] + v_c[i-1][j+1])/2;

                alc_w = alpha_c[i-1][j];
                alc_e = alpha_c[i][j];
                alc_n = ( alpha_c[i-1][j] + alpha_c[i-1][j+1] + alpha_c[i][j] + alpha_c[i][j+1] )/4;
                alc_s = ( alpha_c[i-1][j] + alpha_c[i-1][j-1] + alpha_c[i][j] + alpha_c[i][j-1] )/4;

                mu_w = alc_w*mu_total_C[i-1][j];
                mu_e = alc_e*mu_total_C[i][j];
                mu_n = alc_n*( mu_total_C[i-1][j+1] + mu_total_C[i][j+1] + mu_total_C[i-1][j] + mu_total_C[i][j] )/4;
                mu_s = alc_s*( mu_total_C[i-1][j-1] + mu_total_C[i][j-1] + mu_total_C[i-1][j] + mu_total_C[i][j] )/4;

                m = j*nx + i-1;
                alphac_centre = ( alpha_c[i-1][j] + alpha_c[i][j] )/2;
                alphacn_centre = ( alphan_c[i-1][j] + alphan_c[i][j] )/2;
                
                fE = uE*alc_e;
                fW = uW*alc_w;
                fN = vN*alc_n;
                fS = vS*alc_s;

                z_e = fE > 0;
                z_w = fW > 0;
                z_n = fN > 0;
                z_s = fS > 0;

                uc_ae[i][j] = -(max(-fE,0) + mu_e*ar);
                uc_aw[i][j] = -(max(fW,0) + mu_w*ar);
                uc_an[i][j] = -(max(-fN,0) + mu_n*ra);
                uc_as[i][j] = -(max(fS,0) + mu_s*ra);
                uc_ap[i][j] = -(uc_ae[i][j] + uc_aw[i][j] + uc_an[i][j] + uc_as[i][j]) + alphac_centre*RHO_C*DX2*DT + fE - fW + fN - fS;
                Qxc[m] = alphac_centre*( P[i-1][j] - P[i][j] )*DY + alphacn_centre*RHO_C*un_c[i][j]*DX2*DT - ( 0.67*( alpha_c[i][j]*mu_total_C[i][j]*DIV_xc[i][j] - alpha_c[i-1][j]*mu_total_C[i-1][j]*DIV_xc[i-1][j] ) )*DX - F_added_x[i][j] + 0.5*fE*( u_c[i+1][j] - u_c[i][j] )*( (1 - z_e)*psi_x_neg_x(u_c,i+1,j,0,0) - z_e*psi_x_pos_x(u_c,i,j,0,0) ) + 0.5*fW*( u_c[i][j] - u_c[i-1][j] )*( z_w*psi_x_pos_x(u_c,i-1,j,0,0) - (1 - z_w)*psi_x_neg_x(u_c,i,j,0,0) ) + 0.5*fN*( u_c[i][j+1] - u_c[i][j] )*( (1 - z_n)*psi_x_neg_y(u_c,i,j+1,0,0) - z_n*psi_x_pos_y(u_c,i,j,0,0,u_cinlet) ) + 0.5*fS*( u_c[i][j] - u_c[i][j-1] )*( z_s*psi_x_pos_y(u_c,i,j-1,0,0,u_cinlet) - (1 - z_s)*psi_x_neg_y(u_c,i,j,0,0) );// - ( 0.67*RHO_C*alphac_centre*( k[i][j] - k[i-1][j] ) )*DX;

            }   
        }
        
        // Bottom Boundary
        
        j = 0;
        
        for (i=1;i<NX;i++)
        {
            uE = DY*RHO_C*(u_c[i][j] + u_c[i+1][j])/2;
            uW = DY*RHO_C*(u_c[i][j] + u_c[i-1][j])/2;
            vS = DX*RHO_C*(v_c[i][j] + v_c[i-1][j])/2;
            vN = DX*RHO_C*(v_c[i][j+1] + v_c[i-1][j+1])/2;
            
            alc_w = alpha_c[i-1][j];
            alc_e = alpha_c[i][j];
            alc_n = ( alpha_c[i-1][j] + alpha_c[i-1][j+1] + alpha_c[i][j] + alpha_c[i][j+1] )/4;
            alc_s = ( alphac_inlet[i] + alphac_inlet[i-1] )/2;
            
            mu_w = alc_w*mu_total_C[i-1][j];
            mu_e = alc_e*mu_total_C[i][j];
            mu_n = alc_n*( mu_total_C[i-1][j+1] + mu_total_C[i][j+1] + mu_total_C[i-1][j] + mu_total_C[i][j] )/4;
            mu_s = alc_s*muxc_inlet[i];

            m = j*nx + i-1;
            alphac_centre = ( alpha_c[i-1][j] + alpha_c[i][j] )/2;
            alphacn_centre = ( alphan_c[i-1][j] + alphan_c[i][j] )/2;
        
            fE = uE*alc_e;
            fW = uW*alc_w;
            fN = vN*alc_n;
            fS = vS*alc_s;

            z_e = fE > 0;
            z_w = fW > 0;
            z_n = fN > 0;
            z_s = fS > 0;

            uc_ap[i][j] = max(-fE,0) + max(fW,0) + max(-fN,0) + mu_e*ar + mu_w*ar + mu_n*ra + 3*mu_s*ra + alphac_centre*RHO_C*DX2*DT + fE - fW + fN - fS;
            uc_ae[i][j] = -(max(-fE,0) + mu_e*ar);
            uc_aw[i][j] = -(max(fW,0) + mu_w*ar);
            uc_as[i][j] = 0;
            uc_an[i][j] = -(max(-fN,0) + mu_n*ra + mu_s*ra/3);
            Qxc[m] = alphac_centre*( P[i-1][j] - P[i][j] )*DY + alphacn_centre*RHO_C*un_c[i][j]*DX2*DT - ( 0.67*( alpha_c[i][j]*mu_total_C[i][j]*DIV_xc[i][j] - alpha_c[i-1][j]*mu_total_C[i-1][j]*DIV_xc[i-1][j] ) )*DX - F_added_x[i][j] + 0.5*fE*( u_c[i+1][j] - u_c[i][j] )*( (1 - z_e)*psi_x_neg_x(u_c,i+1,j,0,0) - z_e*psi_x_pos_x(u_c,i,j,0,0) ) + 0.5*fW*( u_c[i][j] - u_c[i-1][j] )*( z_w*psi_x_pos_x(u_c,i-1,j,0,0) - (1 - z_w)*psi_x_neg_x(u_c,i,j,0,0) ) + 0.5*fN*( u_c[i][j+1] - u_c[i][j] )*( (1 - z_n)*psi_x_neg_y(u_c,i,j+1,1,0) - z_n*psi_x_pos_y(u_c,i,j,1,0,u_cinlet) );// - ( 0.67*RHO_C*alphac_centre*( k[i][j] - k[i-1][j] ) )*DX;
        }
                
        // Outlet
        
        j = NY-1;
        
        for (i=1;i<NX;i++)
        {
            uE = DY*RHO_C*(u_c[i][j] + u_c[i+1][j])/2;
            uW = DY*RHO_C*(u_c[i][j] + u_c[i-1][j])/2;
            vS = DX*RHO_C*(v_c[i][j] + v_c[i-1][j])/2;
            vN = DX*RHO_C*(v_c[i][j+1] + v_c[i-1][j+1])/2;
            
            alc_w = alpha_c[i-1][j];
            alc_e = alpha_c[i][j];
            alc_s = ( alpha_c[i-1][j] + alpha_c[i-1][j-1] + alpha_c[i][j] + alpha_c[i][j-1] )/4;
            alc_n = ( alc_w + alc_e )/2;
            
            mu_w = alc_w*mu_total_C[i-1][j];
            mu_e = alc_e*mu_total_C[i][j];
            mu_s = alc_s*( mu_total_C[i-1][j-1] + mu_total_C[i][j-1] + mu_total_C[i-1][j] + mu_total_C[i][j] )/4;
            mu_n = alc_n*( mu_total_C[i][j] + mu_total_C[i-1][j] )/2;

            m = j*nx + i-1;
            alphac_centre = ( alpha_c[i-1][j] + alpha_c[i][j] )/2;
            alphacn_centre = ( alphan_c[i-1][j] + alphan_c[i][j] )/2;
            
            fE = uE*alc_e;
            fW = uW*alc_w;
            fN = vN*alc_n;
            fS = vS*alc_s;

            z_e = fE > 0;
            z_w = fW > 0;
            z_n = fN > 0;
            z_s = fS > 0;

            uc_ae[i][j] = -(max(-fE,0) + mu_e*ar);
            uc_aw[i][j] = -(max(fW,0) + mu_w*ar);
            uc_an[i][j] = 0;
            uc_as[i][j] = -(max(fS,0) + mu_s*ra);
            uc_ap[i][j] = -(uc_ae[i][j] + uc_aw[i][j] + uc_an[i][j] + uc_as[i][j]) + alphac_centre*RHO_C*DX2*DT + fE - fW + fN - fS;
            Qxc[m] = alphac_centre*( P[i-1][j] - P[i][j] )*DY + alphacn_centre*RHO_C*un_c[i][j]*DX2*DT - ( 0.67*( alpha_c[i][j]*mu_total_C[i][j]*DIV_xc[i][j] - alpha_c[i-1][j]*mu_total_C[i-1][j]*DIV_xc[i-1][j] ) )*DX - F_added_x[i][j] + 0.5*fE*( u_c[i+1][j] - u_c[i][j] )*( (1 - z_e)*psi_x_neg_x(u_c,i+1,j,0,0) - z_e*psi_x_pos_x(u_c,i,j,0,0) ) + 0.5*fW*( u_c[i][j] - u_c[i-1][j] )*( z_w*psi_x_pos_x(u_c,i-1,j,0,0) - (1 - z_w)*psi_x_neg_x(u_c,i,j,0,0) ) + 0.5*fS*( u_c[i][j] - u_c[i][j-1] )*( z_s*psi_x_pos_y(u_c,i,j-1,1,1,u_cinlet) - (1 - z_s)*psi_x_neg_y(u_c,i,j,1,1) );// - ( 0.67*RHO_C*alphac_centre*( k[i][j] - k[i-1][j] ) )*DX;
        }
        
        // y-momentum links
        
        double ** vc_ap = calloc(NX,sizeof(double *));
        double ** vc_an = calloc(NX,sizeof(double *));
        double ** vc_as = calloc(NX,sizeof(double *));
        double ** vc_ae = calloc(NX,sizeof(double *));
        double ** vc_aw = calloc(NX,sizeof(double *));
        for (i=0;i<NX;i++)
        {
            vc_ap[i] = calloc(NY,sizeof(double));
            vc_an[i] = calloc(NY,sizeof(double));
            vc_as[i] = calloc(NY,sizeof(double));
            vc_ae[i] = calloc(NY,sizeof(double));
            vc_aw[i] = calloc(NY,sizeof(double));
        }


        /// Interior Cells
        
        for (i=1;i<NX-1;i++)
        {
            for (j=1;j<NY;j++)
            {
            	U_mc = (u_c[i][j] + u_c[i+1][j] + u_c[i+1][j-1] + u_c[i][j-1])/4;
            	U_md = (u_d[i][j] + u_d[i+1][j] + u_d[i+1][j-1] + u_d[i][j-1])/4;

            	mat_der_y[i][j] = U_mc*( v_c[i+1][j] - v_c[i-1][j] )/(2*DY) + v_c[i][j]*( v_c[i][j+1] - v_c[i][j-1] )/(2*DY) - U_md*( v_d[i+1][j] - v_d[i-1][j] )/(2*DY) + v_d[i][j]*( v_d[i][j+1] - v_d[i][j-1] )/(2*DY);

            }
        }
        
        /// Left Wall
        
        i = 0;
        
        for (j=1;j<NY;j++)
        {
        	U_mc = (u_c[i][j] + u_c[i+1][j] + u_c[i+1][j-1] + u_c[i][j-1])/4;
        	U_md = (u_d[i][j] + u_d[i+1][j] + u_d[i+1][j-1] + u_d[i][j-1])/4;

        	mat_der_y[i][j] = U_mc*( (v_c[i+1][j] + v_c[i][j])/2 )/(2*DY) + v_c[i][j]*( v_c[i][j+1] - v_c[i][j-1] )/(2*DY) - U_md*( (v_d[i+1][j] + v_d[i][j])/2 )/(2*DY) - v_d[i][j]*( v_d[i][j+1] - v_d[i][j-1] )/(2*DY);
        }
                
        /// Right wall
        
        i = NX-1;
        
        for (j=1;j<NY;j++)
        {

        	U_mc = (u_c[i][j] + u_c[i+1][j] + u_c[i+1][j-1] + u_c[i][j-1])/4;
        	U_md = (u_d[i][j] + u_d[i+1][j] + u_d[i+1][j-1] + u_d[i][j-1])/4;

        	mat_der_y[i][j] = U_mc*( - (v_c[i-1][j] + v_c[i][j])/2 )/(2*DY) + v_c[i][j]*( v_c[i][j+1] - v_c[i][j-1] )/(2*DY) - U_md*( - (v_d[i-1][j] + v_d[i][j])/2 )/(2*DY) - v_d[i][j]*( v_d[i][j+1] - v_d[i][j-1] )/(2*DY);
        }

        
        for (i=0;i<NX;i++)
        {
            for (j=1;j<NY;j++)
            {
                alphad_centre = ( alpha_d[i][j-1] + alpha_d[i][j] )/2;
                alphac_centre = ( alpha_c[i][j-1] + alpha_c[i][j] )/2;
                U_mc = (u_c[i][j] + u_c[i+1][j] + u_c[i+1][j-1] + u_c[i][j-1])/4;
                V_mc = v_c[i][j];

                U_md = (u_d[i][j] + u_d[i+1][j] + u_d[i+1][j-1] + u_d[i][j-1])/4;
                V_md = v_d[i][j];

                U_diffmagy = sqrt( pow((U_md - U_mc),2) + pow((V_md - V_mc),2) );

                Re_b = RHO_C*U_diffmagy*DP/mu;

                if (Re_b<1000)
                {
                    C_D = 24*( 1 + 0.15*pow(Re_b,0.687) )/Re_b;
                }
                else
                {
                    C_D = 0.44;
                }

                a_if = 6*alphad_centre/DP;

                B_ky[i][j] = C_D*a_if*RHO_C*U_diffmagy*DX2/8;

                DIV_yc[i][j] = ( u_c[i+1][j-1] + u_c[i+1][j] - u_c[i][j-1] - u_c[i][j] )/(2*DX) + ( v_c[i][j+1] - v_c[i][j-1] )/(2*DY);

                DIV_yd[i][j] = ( u_d[i+1][j-1] + u_d[i+1][j] - u_d[i][j-1] - u_d[i][j] )/(2*DX) + ( v_d[i][j+1] - v_d[i][j-1] )/(2*DY);

                F_added_y[i][j] = 0.5*RHO_C*alphad_centre*( ( (v_c[i][j] - v_d[i][j]) - (vn_c[i][j] - vn_d[i][j]) )*DT + mat_der_y[i][j] )*DX2;// + U_mc*( v_c[i+1][j] - v_c[i-1][j] )/(2*DY) + v_c[i][j]*( v_c[i][j+1] - v_c[i][j-1] )/(2*DY) - U_md*( v_d[i+1][j] - v_d[i-1][j] )/(2*DY) + v_d[i][j]*( v_d[i][j+1] - v_d[i][j-1] )/(2*DY) )*DX2;
            }
        }

        /// Interior Cells
        
        for (i=1;i<NX-1;i++)
        {
            for (j=1;j<NY;j++)
            {
                vN = DX*RHO_C*(v_c[i][j] + v_c[i][j+1])/2;
                vS = DX*RHO_C*(v_c[i][j] + v_c[i][j-1])/2;
                uE = DY*RHO_C*(u_c[i+1][j-1] + u_c[i+1][j])/2;
                uW = DY*RHO_C*(u_c[i][j-1] + u_c[i][j])/2;

                alc_n = alpha_c[i][j];
                alc_s = alpha_c[i][j-1];
                alc_e = ( alpha_c[i][j] + alpha_c[i+1][j] + alpha_c[i][j-1] + alpha_c[i+1][j-1] )/4;
                alc_w = ( alpha_c[i][j] + alpha_c[i-1][j] + alpha_c[i][j-1] + alpha_c[i-1][j-1] )/4;
                
                mu_n = alc_n*mu_total_C[i][j];
                mu_s = alc_s*mu_total_C[i][j-1];
                mu_e = alc_e*( mu_total_C[i][j-1] + mu_total_C[i][j] + mu_total_C[i+1][j-1] + mu_total_C[i+1][j] )/4;
                mu_w = alc_w*( mu_total_C[i][j-1] + mu_total_C[i][j] + mu_total_C[i-1][j-1] + mu_total_C[i-1][j] )/4;

                m = (j-1)*NX + i;
                alphac_centre = ( alpha_c[i][j-1] + alpha_c[i][j] )/2;
                alphacn_centre = ( alphan_c[i][j-1] + alphan_c[i][j] )/2;
                
                fE = uE*alc_e;
                fW = uW*alc_w;
                fN = vN*alc_n;
                fS = vS*alc_s;

                z_e = fE > 0;
                z_w = fW > 0;
                z_n = fN > 0;
                z_s = fS > 0;

                vc_ae[i][j] = -(max(-fE,0) + mu_e*ar);
                vc_aw[i][j] = -(max(fW,0) + mu_w*ar);
                vc_an[i][j] = -(max(-fN,0) + mu_n*ra);
                vc_as[i][j] = -(max(fS,0) + mu_s*ra);
                vc_ap[i][j] = -(vc_ae[i][j] + vc_aw[i][j] + vc_an[i][j] + vc_as[i][j]) + alphac_centre*RHO_C*DX2*DT + fE - fW + fN - fS;
                Qyc[m] = alphac_centre*( P[i][j-1] - P[i][j] )*DX + alphac_centre*RHO_C*g*DX2 + alphacn_centre*RHO_C*vn_c[i][j]*DX2*DT - 0.67*( alpha_c[i][j]*mu_total_C[i][j]*DIV_yc[i][j] - alpha_c[i][j-1]*mu_total_C[i][j-1]*DIV_yc[i][j-1] )*DX - F_added_y[i][j] + 0.5*fE*( v_c[i+1][j] - v_c[i][j] )*( (1 - z_e)*psi_y_neg_x(v_c,i+1,j,0,0) - z_e*psi_y_pos_x(v_c,i,j,0,0) ) + 0.5*fW*( v_c[i][j] - v_c[i-1][j] )*( z_w*psi_y_pos_x(v_c,i-1,j,0,0) - (1 - z_w)*psi_y_neg_x(v_c,i,j,0,0) ) + 0.5*fN*( v_c[i][j+1] - v_c[i][j] )*( (1 - z_n)*psi_y_neg_y(v_c,i,j+1,0,0) - z_n*psi_y_pos_y(v_c,i,j,0,0,u_cinlet) ) + 0.5*fS*( v_c[i][j] - v_c[i][j-1] )*( z_s*psi_y_pos_y(v_c,i,j-1,0,0,u_cinlet) - (1 - z_s)*psi_y_neg_y(v_c,i,j,0,0) );// - ( 0.67*RHO_C*alphac_centre*( k[i][j] - k[i][j-1] ) )*DX;
            }
        }
        
        /// Left Wall
        
        i = 0;
        
        for (j=1;j<NY;j++)
        {
            vN = DX*RHO_C*(v_c[i][j] + v_c[i][j+1])/2;
            vS = DX*RHO_C*(v_c[i][j] + v_c[i][j-1])/2;
            uE = DY*RHO_C*(u_c[i+1][j-1] + u_c[i+1][j])/2;
            uW = DY*RHO_C*(u_c[i][j-1] + u_c[i][j])/2;

            alc_n = alpha_c[i][j];
            alc_s = alpha_c[i][j-1];
            alc_e = ( alpha_c[i][j] + alpha_c[i+1][j] + alpha_c[i][j-1] + alpha_c[i+1][j-1] )/4;
            alc_w = ( alc_s + alc_n )/2;

            mu_n = alc_n*mu_total_C[i][j];
            mu_s = alc_s*mu_total_C[i][j-1];
            mu_e = alc_e*( mu_total_C[i][j-1] + mu_total_C[i][j] + mu_total_C[i+1][j-1] + mu_total_C[i+1][j] )/4;
            mu_w = alc_w*mu;

            m = (j-1)*NX + i;
            alphac_centre = ( alpha_c[i][j-1] + alpha_c[i][j] )/2;
            alphacn_centre = ( alphan_c[i][j-1] + alphan_c[i][j] )/2;

            kp = (k[i][j] + k[i][j-1])/2;
            x_plus = RHO_C*xp*pow(C_MU,0.25)*pow(kp,0.5)/mu;
            v_plus = log(EE*x_plus)/KAPPA;
            
            fE = uE*alc_e;
            fW = uW*alc_w;
            fN = vN*alc_n;
            fS = vS*alc_s;

            z_e = fE > 0;
            z_w = fW > 0;
            z_n = fN > 0;
            z_s = fS > 0;

            vc_ap[i][j] = max(-fN,0) + max(fS,0) + max(-fE,0) + mu_n*ra + mu_s*ra + 3*mu_w*ar + mu_e*ar + RHO_C*(pow(C_MU,0.25)*pow(kp,0.5))*DX/v_plus + alphac_centre*RHO_C*DX2*DT + fE - fW + fN - fS;
            vc_ae[i][j] = -(max(-fE,0) + mu_e*ar + mu_w*ar/3);
            vc_aw[i][j] = 0;
            vc_an[i][j] = -(max(-fN,0) + mu_n*ra);
            vc_as[i][j] = -(max(fS,0) + mu_s*ra);
            Qyc[m] = alphac_centre*( P[i][j-1] - P[i][j] )*DX + alphac_centre*RHO_C*g*DX2 + alphacn_centre*RHO_C*vn_c[i][j]*DX2*DT - 0.67*( alpha_c[i][j]*mu_total_C[i][j]*DIV_yc[i][j] - alpha_c[i][j-1]*mu_total_C[i][j-1]*DIV_yc[i][j-1] )*DX - F_added_y[i][j] + 0.5*fE*( v_c[i+1][j] - v_c[i][j] )*( (1 - z_e)*psi_y_neg_x(v_c,i+1,j,1,0) - z_e*psi_y_pos_x(v_c,i,j,1,0) ) + 0.5*fN*( v_c[i][j+1] - v_c[i][j] )*( (1 - z_n)*psi_y_neg_y(v_c,i,j+1,0,0) - z_n*psi_y_pos_y(v_c,i,j,0,0,u_cinlet) ) + 0.5*fS*( v_c[i][j] - v_c[i][j-1] )*( z_s*psi_y_pos_y(v_c,i,j-1,0,0,u_cinlet) - (1 - z_s)*psi_y_neg_y(v_c,i,j,0,0) );// - ( 0.67*RHO_C*alphac_centre*( k[i][j] - k[i][j-1] ) )*DX;
        }
                
        /// Right wall
        
        i = NX-1;
        
        for (j=1;j<NY;j++)
        {
            vN = DX*RHO_C*(v_c[i][j] + v_c[i][j+1])/2;
            vS = DX*RHO_C*(v_c[i][j] + v_c[i][j-1])/2;
            uE = DY*RHO_C*(u_c[i+1][j-1] + u_c[i+1][j])/2;
            uW = DY*RHO_C*(u_c[i][j-1] + u_c[i][j])/2;

            alc_n = alpha_c[i][j];
            alc_s = alpha_c[i][j-1];
            alc_w = ( alpha_c[i][j] + alpha_c[i-1][j] + alpha_c[i][j-1] + alpha_c[i-1][j-1])/4;
            alc_e = (alc_s + alc_n)/2;

            mu_n = alc_n*mu_total_C[i][j];
            mu_s = alc_s*mu_total_C[i][j-1];
            mu_w = alc_w*( mu_total_C[i][j-1] + mu_total_C[i][j] + mu_total_C[i-1][j-1] + mu_total_C[i-1][j] )/4;
            mu_e = alc_e*mu;

            m = (j-1)*NX + i;
            alphac_centre = ( alpha_c[i][j-1] + alpha_c[i][j] )/2;
            alphacn_centre = ( alphan_c[i][j-1] + alphan_c[i][j] )/2;
            
            kp = (k[i][j] + k[i][j-1])/2;
            x_plus = RHO_C*xp*pow(C_MU,0.25)*pow(kp,0.5)/mu;
            v_plus = log(EE*x_plus)/KAPPA;
            // printf("%12.5e\n", x_plus);

            fE = uE*alc_e;
            fW = uW*alc_w;
            fN = vN*alc_n;
            fS = vS*alc_s;

            z_e = fE > 0;
            z_w = fW > 0;
            z_n = fN > 0;
            z_s = fS > 0;

            vc_ap[i][j] = max(fW,0) + max(-fN,0) + max(fS,0) + mu_n*ra + mu_s*ra + mu_w*ar + 3*mu_e*ar + RHO_C*(pow(C_MU,0.25)*pow(kp,0.5))*DX/v_plus + alphac_centre*RHO_C*DX2*DT + fE - fW + fN - fS;
            vc_ae[i][j] = 0;
            vc_aw[i][j] = -(max(fW,0) + mu_w*ar + mu_e*ar/3);
            vc_an[i][j] = -(max(-fN,0) + mu_n*ra);
            vc_as[i][j] = -(max(fS,0) + mu_s*ra);
            Qyc[m] = alphac_centre*( P[i][j-1] - P[i][j] )*DX + alphac_centre*RHO_C*g*DX2 + alphacn_centre*RHO_C*vn_c[i][j]*DX2*DT - 0.67*( alpha_c[i][j]*mu_total_C[i][j]*DIV_yc[i][j] - alpha_c[i][j-1]*mu_total_C[i][j-1]*DIV_yc[i][j-1] )*DX - F_added_y[i][j] + 0.5*fW*( v_c[i][j] - v_c[i-1][j] )*( z_w*psi_y_pos_x(v_c,i-1,j,1,1) - (1 - z_w)*psi_y_neg_x(v_c,i,j,1,1) ) + 0.5*fN*( v_c[i][j+1] - v_c[i][j] )*( (1 - z_n)*psi_y_neg_y(v_c,i,j+1,0,0) - z_n*psi_y_pos_y(v_c,i,j,0,0,u_cinlet) ) + 0.5*fS*( v_c[i][j] - v_c[i][j-1] )*( z_s*psi_y_pos_y(v_c,i,j-1,0,0,u_cinlet) - (1 - z_s)*psi_y_neg_y(v_c,i,j,0,0) );// - ( 0.67*RHO_C*alphac_centre*( k[i][j] - k[i][j-1] ) )*DX;

        }
        
        //// x-momentum co-efficients for dispersed phase

        double ** ud_ap = calloc(NX,sizeof(double *));
        double ** ud_an = calloc(NX,sizeof(double *));
        double ** ud_as = calloc(NX,sizeof(double *));
        double ** ud_ae = calloc(NX,sizeof(double *));
        double ** ud_aw = calloc(NX,sizeof(double *));
        for (i=0;i<NX;i++)
        {
            ud_ap[i] = calloc(NY,sizeof(double));
            ud_an[i] = calloc(NY,sizeof(double));
            ud_as[i] = calloc(NY,sizeof(double));
            ud_ae[i] = calloc(NY,sizeof(double));
            ud_aw[i] = calloc(NY,sizeof(double));
        }
        
        /// Interior Cells

        
        for (i=1;i<NX;i++)
        {
            for (j=1;j<NY-1;j++)
            {
                uE = DY*RHO_D*(u_d[i][j] + u_d[i+1][j])/2;
                uW = DY*RHO_D*(u_d[i][j] + u_d[i-1][j])/2;
                vS = DX*RHO_D*(v_d[i][j] + v_d[i-1][j])/2;
                vN = DX*RHO_D*(v_d[i][j+1] + v_d[i-1][j+1])/2;

                ald_w = alpha_d[i-1][j];
                ald_e = alpha_d[i][j];
                ald_n = ( alpha_d[i-1][j] + alpha_d[i-1][j+1] + alpha_d[i][j] + alpha_d[i][j+1] )/4;
                ald_s = ( alpha_d[i-1][j] + alpha_d[i-1][j-1] + alpha_d[i][j] + alpha_d[i][j-1] )/4;
                
                mu_w = ald_w*mu_total_D[i-1][j];
                mu_e = ald_e*mu_total_D[i][j];
                mu_n = ald_n*( mu_total_D[i-1][j+1] + mu_total_D[i][j+1] + mu_total_D[i-1][j] + mu_total_D[i][j] )/4;
                mu_s = ald_s*( mu_total_D[i-1][j-1] + mu_total_D[i][j-1] + mu_total_D[i-1][j] + mu_total_D[i][j] )/4;

                m = j*nx + i-1;
                alphad_centre = ( alpha_d[i-1][j] + alpha_d[i][j] )/2;
                alphadn_centre = ( alphan_d[i-1][j] + alphan_d[i][j] )/2;

                fE = uE*ald_e;
                fW = uW*ald_w;
                fN = vN*ald_n;
                fS = vS*ald_s;

                z_e = fE > 0;
                z_w = fW > 0;
                z_n = fN > 0;
                z_s = fS > 0;

                ud_ae[i][j] = -(max(-fE,0) + mu_e*ar);
                ud_aw[i][j] = -(max(fW,0) + mu_w*ar);
                ud_an[i][j] = -(max(-fN,0) + mu_n*ra);
                ud_as[i][j] = -(max(fS,0) + mu_s*ra);
                ud_ap[i][j] = -(ud_ae[i][j] + ud_aw[i][j] + ud_an[i][j] + ud_as[i][j]) + alphad_centre*RHO_D*DX2*DT + fE - fW + fN - fS;
                Qxd[m] = alphad_centre*( P[i-1][j] - P[i][j] )*DY + alphadn_centre*RHO_D*un_d[i][j]*DX2*DT - 0.67*( alpha_d[i][j]*mu_total_D[i][j]*DIV_xd[i][j] - alpha_d[i-1][j]*mu_total_D[i-1][j]*DIV_xd[i-1][j] )*DX + F_added_x[i][j] + 0.5*fE*( u_d[i+1][j] - u_d[i][j] )*( (1 - z_e)*psi_x_neg_x(u_d,i+1,j,0,0) - z_e*psi_x_pos_x(u_d,i,j,0,0) ) + 0.5*fW*( u_d[i][j] - u_d[i-1][j] )*( z_w*psi_x_pos_x(u_d,i-1,j,0,0) - (1 - z_w)*psi_x_neg_x(u_d,i,j,0,0) ) + 0.5*fN*( u_d[i][j+1] - u_d[i][j] )*( (1 - z_n)*psi_x_neg_y(u_d,i,j+1,0,0) - z_n*psi_x_pos_y(u_d,i,j,0,0,u_dinlet) ) + 0.5*fS*( u_d[i][j] - u_d[i][j-1] )*( z_s*psi_x_pos_y(u_d,i,j-1,0,0,u_dinlet) - (1 - z_s)*psi_x_neg_y(u_d,i,j,0,0) );

            }
            
        }

        // Bottom Boundary
        
        j = 0;
        
        for (i=1;i<NX;i++)
        {
            uE = DY*RHO_D*(u_d[i][j] + u_d[i+1][j])/2;
            uW = DY*RHO_D*(u_d[i][j] + u_d[i-1][j])/2;
            vS = DX*RHO_D*(v_d[i][j] + v_d[i-1][j])/2;
            vN = DX*RHO_D*(v_d[i][j+1] + v_d[i-1][j+1])/2;
            
            ald_w = alpha_d[i-1][j];
            ald_e = alpha_d[i][j];
            ald_n = ( alpha_d[i-1][j] + alpha_d[i-1][j+1] + alpha_d[i][j] + alpha_d[i][j+1] )/4;
            ald_s = ( alphad_inlet[i] + alphad_inlet[i-1] )/2;

            mu_w = ald_w*mu_total_D[i-1][j];
            mu_e = ald_e*mu_total_D[i][j];
            mu_n = ald_n*( mu_total_D[i-1][j+1] + mu_total_D[i][j+1] + mu_total_D[i-1][j] + mu_total_D[i][j] )/4;
            mu_s = ald_s*muxd_inlet[i];

            m = j*nx + i-1;
            alphad_centre = ( alpha_d[i-1][j] + alpha_d[i][j] )/2;
            alphadn_centre = ( alphan_d[i-1][j] + alphan_d[i][j] )/2;
            
            fE = uE*ald_e;
            fW = uW*ald_w;
            fN = vN*ald_n;
            fS = vS*ald_s;

            z_e = fE > 0;
            z_w = fW > 0;
            z_n = fN > 0;
            z_s = fS > 0;

            ud_ap[i][j] = max(-fE,0) + max(fW,0) + max(-fN,0) + mu_e*ar + mu_w*ar + mu_n*ra + 3*mu_s*ra + alphad_centre*RHO_D*DX2*DT + fE - fW + fN - fS;
            ud_ae[i][j] = -(max(-fE,0) + mu_e*ar);
            ud_aw[i][j] = -(max(fW,0) + mu_w*ar);
            ud_as[i][j] = 0;
            ud_an[i][j] = -(max(-fN,0) + mu_n*ra + mu_s*ra/3);
            Qxd[m] = alphad_centre*( P[i-1][j] - P[i][j] )*DY + alphadn_centre*RHO_D*un_d[i][j]*DX2*DT - 0.67*( alpha_d[i][j]*mu_total_D[i][j]*DIV_xd[i][j] - alpha_d[i-1][j]*mu_total_D[i-1][j]*DIV_xd[i-1][j] )*DX + F_added_x[i][j] + 0.5*fE*( u_d[i+1][j] - u_d[i][j] )*( (1 - z_e)*psi_x_neg_x(u_d,i+1,j,0,0) - z_e*psi_x_pos_x(u_d,i,j,0,0) ) + 0.5*fW*( u_d[i][j] - u_d[i-1][j] )*( z_w*psi_x_pos_x(u_d,i-1,j,0,0) - (1 - z_w)*psi_x_neg_x(u_d,i,j,0,0) ) + 0.5*fN*( u_d[i][j+1] - u_d[i][j] )*( (1 - z_n)*psi_x_neg_y(u_d,i,j+1,1,0) - z_n*psi_x_pos_y(u_d,i,j,1,0,u_dinlet) );
        }
                
        // Outlet
        
        j = NY-1;
        
        for (i=1;i<NX;i++)
        {
            uE = DY*RHO_D*(u_d[i][j] + u_d[i+1][j])/2;
            uW = DY*RHO_D*(u_d[i][j] + u_d[i-1][j])/2;
            vS = DX*RHO_D*(v_d[i][j] + v_d[i-1][j])/2;
            vN = DX*RHO_D*(v_d[i][j+1] + v_d[i-1][j+1])/2;
            
            ald_w = alpha_d[i-1][j];
            ald_e = alpha_d[i][j];
            ald_s = ( alpha_d[i-1][j] + alpha_d[i-1][j-1] + alpha_d[i][j] + alpha_d[i][j-1] )/4;
            ald_n = ( ald_w + ald_e )/2;

            mu_w = ald_w*mu_total_D[i-1][j];
            mu_e = ald_e*mu_total_D[i][j];
            mu_s = ald_s*( mu_total_D[i-1][j-1] + mu_total_D[i][j-1] + mu_total_D[i-1][j] + mu_total_D[i][j] )/4;
            mu_n = ald_n*( mu_total_D[i][j] + mu_total_D[i-1][j] )/2;

            m = j*nx + i-1;
            alphad_centre = ( alpha_d[i-1][j] + alpha_d[i][j] )/2;
            alphadn_centre = ( alphan_d[i-1][j] + alphan_d[i][j] )/2;
            
            fE = uE*ald_e;
            fW = uW*ald_w;
            fN = vN*ald_n;
            fS = vS*ald_s;

            z_e = fE > 0;
            z_w = fW > 0;
            z_n = fN > 0;
            z_s = fS > 0;

            ud_ae[i][j] = -(max(-fE,0) + mu_e*ar);
            ud_aw[i][j] = -(max(fW,0) + mu_w*ar);
            ud_an[i][j] = 0;
            ud_as[i][j] = -(max(fS,0) + mu_s*ra);
            ud_ap[i][j] = -(ud_ae[i][j] + ud_aw[i][j] + ud_an[i][j] + ud_as[i][j]) + alphad_centre*RHO_D*DX2*DT + fE - fW + fN - fS;
            Qxd[m] = alphad_centre*( P[i-1][j] - P[i][j] )*DY + alphadn_centre*RHO_D*un_d[i][j]*DX2*DT - 0.67*( alpha_d[i][j]*mu_total_D[i][j]*DIV_xd[i][j] - alpha_d[i-1][j]*mu_total_D[i-1][j]*DIV_xd[i-1][j] )*DX + F_added_x[i][j] + 0.5*fE*( u_d[i+1][j] - u_d[i][j] )*( (1 - z_e)*psi_x_neg_x(u_d,i+1,j,0,0) - z_e*psi_x_pos_x(u_d,i,j,0,0) ) + 0.5*fW*( u_d[i][j] - u_d[i-1][j] )*( z_w*psi_x_pos_x(u_d,i-1,j,0,0) - (1 - z_w)*psi_x_neg_x(u_d,i,j,0,0) ) + 0.5*fS*( u_d[i][j] - u_d[i][j-1] )*( z_s*psi_x_pos_y(u_d,i,j-1,1,1,u_dinlet) - (1 - z_s)*psi_x_neg_y(u_d,i,j,1,1) );
        }

        //// y-momentum co-efficients for dispersed phase

        double ** vd_ap = calloc(NX,sizeof(double *));
        double ** vd_an = calloc(NX,sizeof(double *));
        double ** vd_as = calloc(NX,sizeof(double *));
        double ** vd_ae = calloc(NX,sizeof(double *));
        double ** vd_aw = calloc(NX,sizeof(double *));
        for (i=0;i<NX;i++)
        {
            vd_ap[i] = calloc(NY,sizeof(double));
            vd_an[i] = calloc(NY,sizeof(double));
            vd_as[i] = calloc(NY,sizeof(double));
            vd_ae[i] = calloc(NY,sizeof(double));
            vd_aw[i] = calloc(NY,sizeof(double));
        }
        /// Interior Cells
        
        for (i=1;i<NX-1;i++)
        {
            for (j=1;j<NY;j++)
            {
                vN = DX*RHO_D*(v_d[i][j] + v_d[i][j+1])/2;
                vS = DX*RHO_D*(v_d[i][j] + v_d[i][j-1])/2;
                uE = DY*RHO_D*(u_d[i+1][j-1] + u_d[i+1][j])/2;
                uW = DY*RHO_D*(u_d[i][j-1] + u_d[i][j])/2;

                ald_n = alpha_d[i][j];
                ald_s = alpha_d[i][j-1];
                ald_e = ( alpha_d[i][j] + alpha_d[i+1][j] + alpha_d[i][j-1] + alpha_d[i+1][j-1])/4;
                ald_w = ( alpha_d[i][j] + alpha_d[i-1][j] + alpha_d[i][j-1] + alpha_d[i-1][j-1])/4;
                
                mu_n = ald_n*mu_total_D[i][j];
                mu_s = ald_s*mu_total_D[i][j-1];
                mu_e = ald_e*( mu_total_D[i][j-1] + mu_total_D[i][j] + mu_total_D[i+1][j-1] + mu_total_D[i+1][j] )/4;
                mu_w = ald_w*( mu_total_D[i][j-1] + mu_total_D[i][j] + mu_total_D[i-1][j-1] + mu_total_D[i-1][j] )/4;

                m = (j-1)*NX + i;
                alphad_centre = ( alpha_d[i][j-1] + alpha_d[i][j] )/2;
                alphadn_centre = ( alphan_d[i][j-1] + alphan_d[i][j] )/2;

                fE = uE*ald_e;
                fW = uW*ald_w;
                fN = vN*ald_n;
                fS = vS*ald_s;

                z_e = fE > 0;
                z_w = fW > 0;
                z_n = fN > 0;
                z_s = fS > 0;

                vd_ae[i][j] = -(max(-fE,0) + mu_e*ar);
                vd_aw[i][j] = -(max(fW,0) + mu_w*ar);
                vd_an[i][j] = -(max(-fN,0) + mu_n*ra);
                vd_as[i][j] = -(max(fS,0) + mu_s*ra);
                vd_ap[i][j] = -(vd_ae[i][j] + vd_aw[i][j] + vd_an[i][j] + vd_as[i][j]) + alphad_centre*RHO_D*DX2*DT + fE - fW + fN - fS;
                Qyd[m] = alphad_centre*( P[i][j-1] - P[i][j] )*DX + alphad_centre*RHO_D*g*DX2 + alphadn_centre*RHO_D*vn_d[i][j]*DX2*DT - 0.67*( alpha_d[i][j]*mu_total_D[i][j]*DIV_yd[i][j] - alpha_d[i][j-1]*mu_total_D[i][j-1]*DIV_yd[i][j-1] )*DX + F_added_y[i][j] + 0.5*fE*( v_d[i+1][j] - v_d[i][j] )*( (1 - z_e)*psi_y_neg_x(v_d,i+1,j,0,0) - z_e*psi_y_pos_x(v_d,i,j,0,0) ) + 0.5*fW*( v_d[i][j] - v_d[i-1][j] )*( z_w*psi_y_pos_x(v_d,i-1,j,0,0) - (1 - z_w)*psi_y_neg_x(v_d,i,j,0,0) ) + 0.5*fN*( v_d[i][j+1] - v_d[i][j] )*( (1 - z_n)*psi_y_neg_y(v_d,i,j+1,0,0) - z_n*psi_y_pos_y(v_d,i,j,0,0,v_inlet) ) + 0.5*fS*( v_d[i][j] - v_d[i][j-1] )*( z_s*psi_y_pos_y(v_d,i,j-1,0,0,v_inlet) - (1 - z_s)*psi_y_neg_y(v_d,i,j,0,0) );

            }
        }
        
        /// Left wall
        
        i = 0;
        
        for (j=1;j<NY;j++)
        {
            vN = DX*RHO_D*(v_d[i][j] + v_d[i][j+1])/2;
            vS = DX*RHO_D*(v_d[i][j] + v_d[i][j-1])/2;
            uE = DY*RHO_D*(u_d[i+1][j-1] + u_d[i+1][j])/2;
            uW = DY*RHO_D*(u_d[i][j-1] + u_d[i][j])/2;

            ald_n = alpha_d[i][j];
            ald_s = alpha_d[i][j-1];
            ald_e = ( alpha_d[i][j] + alpha_d[i+1][j] + alpha_d[i][j-1] + alpha_d[i+1][j-1])/4;
            ald_w = ( ald_s + ald_n )/2;

            mu_n = ald_n*mu_total_D[i][j];
            mu_s = ald_s*mu_total_D[i][j-1];
            mu_e = ald_e*( mu_total_D[i][j-1] + mu_total_D[i][j] + mu_total_D[i+1][j-1] + mu_total_D[i+1][j] )/4;
            mu_w = ald_w*mu_dispersed;

            m = (j-1)*NX + i;
            alphad_centre = ( alpha_d[i][j-1] + alpha_d[i][j] )/2;
            alphadn_centre = ( alphan_d[i][j-1] + alphan_d[i][j] )/2;
            
            fE = uE*ald_e;
            fW = uW*ald_w;
            fN = vN*ald_n;
            fS = vS*ald_s;

            z_e = fE > 0;
            z_w = fW > 0;
            z_n = fN > 0;
            z_s = fS > 0;

            vd_ap[i][j] = max(-fN,0) + max(fS,0) + max(-fE,0) + mu_n*ra + mu_s*ra + 3*mu_w*ar + mu_e*ar + alphad_centre*RHO_D*DX2*DT + fE - fW + fN - fS;
            vd_ae[i][j] = -(max(-fE,0) + mu_e*ar + mu_w*ar/3);
            vd_aw[i][j] = 0;
            vd_an[i][j] = -(max(-fN,0) + mu_n*ra);
            vd_as[i][j] = -(max(fS,0) + mu_s*ra);
            Qyd[m] = alphad_centre*( P[i][j-1] - P[i][j] )*DX + alphad_centre*RHO_D*g*DX2 + alphadn_centre*RHO_D*vn_d[i][j]*DX2*DT - 0.67*( alpha_d[i][j]*mu_total_D[i][j]*DIV_yd[i][j] - alpha_d[i][j-1]*mu_total_D[i][j-1]*DIV_yd[i][j-1] )*DX + F_added_y[i][j] + 0.5*fE*( v_d[i+1][j] - v_d[i][j] )*( (1 - z_e)*psi_y_neg_x(v_d,i+1,j,1,0) - z_e*psi_y_pos_x(v_d,i,j,1,0) ) + 0.5*fN*( v_d[i][j+1] - v_d[i][j] )*( (1 - z_n)*psi_y_neg_y(v_d,i,j+1,0,0) - z_n*psi_y_pos_y(v_d,i,j,0,0,v_inlet) ) + 0.5*fS*( v_d[i][j] - v_d[i][j-1] )*( z_s*psi_y_pos_y(v_d,i,j-1,0,0,v_inlet) - (1 - z_s)*psi_y_neg_y(v_d,i,j,0,0) );
        }
                
        /// Right wall
        
        i = NX-1;
        
        for (j=1;j<NY;j++)
        {
            vN = DX*RHO_D*(v_d[i][j] + v_d[i][j+1])/2;
            vS = DX*RHO_D*(v_d[i][j] + v_d[i][j-1])/2;
            uE = DY*RHO_D*(u_d[i+1][j-1] + u_d[i+1][j])/2;
            uW = DY*RHO_D*(u_d[i][j-1] + u_d[i][j])/2;

            ald_n = alpha_d[i][j];
            ald_s = alpha_d[i][j-1];
            ald_w = ( alpha_d[i][j] + alpha_d[i-1][j] + alpha_d[i][j-1] + alpha_d[i-1][j-1])/4;
            ald_e = ( ald_n + ald_s )/2;
            
            mu_n = ald_n*mu_total_D[i][j];
            mu_s = ald_s*mu_total_D[i][j-1];
            mu_w = ald_w*( mu_total_D[i][j-1] + mu_total_D[i][j] + mu_total_D[i-1][j-1] + mu_total_D[i-1][j] )/4;
            mu_e = ald_e*mu_dispersed;
            
            m = (j-1)*NX + i;
            alphad_centre = ( alpha_d[i][j-1] + alpha_d[i][j] )/2;
            alphadn_centre = ( alphan_d[i][j-1] + alphan_d[i][j] )/2;
            
            fE = uE*ald_e;
            fW = uW*ald_w;
            fN = vN*ald_n;
            fS = vS*ald_s;

            z_e = fE > 0;
            z_w = fW > 0;
            z_n = fN > 0;
            z_s = fS > 0;

            vd_ap[i][j] = max(fW,0) + max(-fN,0) + max(fS,0) + mu_n*ra + mu_s*ra + mu_w*ar + 3*mu_e*ar + alphad_centre*RHO_D*DX2*DT + fE - fW + fN - fS;
            vd_ae[i][j] = 0;
            vd_aw[i][j] = -(max(fW,0) + mu_w*ar + mu_e*ar/3);
            vd_an[i][j] = -(max(-fN,0) + mu_n*ra);
            vd_as[i][j] = -(max(fS,0) + mu_s*ra);
            Qyd[m] = alphad_centre*( P[i][j-1] - P[i][j] )*DX + alphad_centre*RHO_D*g*DX2 + alphadn_centre*RHO_D*vn_d[i][j]*DX2*DT - 0.67*( alpha_d[i][j]*mu_total_D[i][j]*DIV_yd[i][j] - alpha_d[i][j-1]*mu_total_D[i][j-1]*DIV_yd[i][j-1] )*DX + F_added_y[i][j] + 0.5*fW*( v_d[i][j] - v_d[i-1][j] )*( z_w*psi_y_pos_x(v_d,i-1,j,1,1) - (1 - z_w)*psi_y_neg_x(v_d,i,j,1,1) ) + 0.5*fN*( v_d[i][j+1] - v_d[i][j] )*( (1 - z_n)*psi_y_neg_y(v_d,i,j+1,0,0) - z_n*psi_y_pos_y(v_d,i,j,0,0,v_inlet) ) + 0.5*fS*( v_d[i][j] - v_d[i][j-1] )*( z_s*psi_y_pos_y(v_d,i,j-1,0,0,v_inlet) - (1 - z_s)*psi_y_neg_y(v_d,i,j,0,0) );
        }

                //// Stone's ILU solver for x-momentum

                for (i=1;i<NX;i++)
                {
                    for (j=0;j<NY;j++)
                    {
                        m = j*nx + i-1;

                        Exc[m] = ( uc_ap[i][j]*ud_ap[i][j] + B_kx[i][j]*( uc_ap[i][j] + ud_ap[i][j] ) )/UV_COR;
                        Fxc[m] = uc_ae[i][j]*( ud_ap[i][j] + B_kx[i][j] );
                        Dxc[m] = uc_aw[i][j]*( ud_ap[i][j] + B_kx[i][j] );
                        Bxc[m] = uc_as[i][j]*( ud_ap[i][j] + B_kx[i][j] );
                        Hxc[m] = uc_an[i][j]*( ud_ap[i][j] + B_kx[i][j] );
                        Uc[m] = u_c[i][j];

                        Exd[m] = ( ud_ap[i][j]*uc_ap[i][j] + B_kx[i][j]*( uc_ap[i][j] + ud_ap[i][j] ) )/UV_COR;
                        Fxd[m] = ud_ae[i][j]*( uc_ap[i][j] + B_kx[i][j] );
                        Dxd[m] = ud_aw[i][j]*( uc_ap[i][j] + B_kx[i][j] );
                        Bxd[m] = ud_as[i][j]*( uc_ap[i][j] + B_kx[i][j] );
                        Hxd[m] = ud_an[i][j]*( uc_ap[i][j] + B_kx[i][j] );
                        Ud[m] = u_d[i][j];
                    }
                }

                for (i=1;i<NX;i++)
                {
                    for (j=1;j<NY-1;j++)
                    {
                        m = j*nx + i-1;

                        Qxcmod[m] = Qxc[m]*( ud_ap[i][j] + B_kx[i][j] ) + B_kx[i][j]*Qxd[m] - B_kx[i][j]*( ud_ae[i][j]*u_d[i+1][j] + ud_aw[i][j]*u_d[i-1][j] + ud_an[i][j]*u_d[i][j+1] + ud_as[i][j]*u_d[i][j-1] ) + (1-UV_COR)*Exc[m]*u_c[i][j];

                        Qxdmod[m] = Qxd[m]*( uc_ap[i][j] + B_kx[i][j] ) + B_kx[i][j]*Qxc[m] - B_kx[i][j]*( uc_ae[i][j]*u_c[i+1][j] + uc_aw[i][j]*u_c[i-1][j] + uc_an[i][j]*u_c[i][j+1] + uc_as[i][j]*u_c[i][j-1] ) + (1-UV_COR)*Exd[m]*u_d[i][j]; 
                    }
                }

                j = 0;

                for (i=1;i<NX;i++)
                {
                    
                    m = j*nx + i-1;

                    Qxcmod[m] = Qxc[m]*( ud_ap[i][j] + B_kx[i][j] ) + B_kx[i][j]*Qxd[m] - B_kx[i][j]*( ud_ae[i][j]*u_d[i+1][j] + ud_aw[i][j]*u_d[i-1][j] + ud_an[i][j]*u_d[i][j+1] ) + (1-UV_COR)*Exc[m]*u_c[i][j];

                    Qxdmod[m] = Qxd[m]*( uc_ap[i][j] + B_kx[i][j] ) + B_kx[i][j]*Qxc[m] - B_kx[i][j]*( uc_ae[i][j]*u_c[i+1][j] + uc_aw[i][j]*u_c[i-1][j] + uc_an[i][j]*u_c[i][j+1] ) + (1-UV_COR)*Exd[m]*u_d[i][j];                     
                }

                j = NY-1;

                for (i=1;i<NX;i++)
                {
                    m = j*nx + i-1;

                    Qxcmod[m] = Qxc[m]*( ud_ap[i][j] + B_kx[i][j] ) + B_kx[i][j]*Qxd[m] - B_kx[i][j]*( ud_ae[i][j]*u_d[i+1][j] + ud_aw[i][j]*u_d[i-1][j] + ud_as[i][j]*u_d[i][j-1] ) + (1-UV_COR)*Exc[m]*u_c[i][j];

                    Qxdmod[m] = Qxd[m]*( uc_ap[i][j] + B_kx[i][j] ) + B_kx[i][j]*Qxc[m] - B_kx[i][j]*( uc_ae[i][j]*u_c[i+1][j] + uc_aw[i][j]*u_c[i-1][j] + uc_as[i][j]*u_c[i][j-1] ) + (1-UV_COR)*Exd[m]*u_d[i][j]; 
                }

                coefficient_matrix(Exc, Fxc, Dxc, Bxc, Hxc, fxc, dxc, bxc, exc, cxc, Kx, nx, ALPHA_XY);

                //// Inner iteration loop begins

                for (inner=1;inner<=ITERATION_X;inner++)
                {
                    // start = clock();

                    residual_x_c = SIP_Solver_x (Exc, Fxc, Dxc, Bxc, Hxc, Qxcmod, fxc, dxc, bxc, exc, cxc, resxc, vec_Yxc, vec_deltaxc, Kx, nx, inner, Uc, u_c, Exc);

                    // end = clock();

                    // if(inner == 1)
                        // printf("%12.5e\n", ((double)(end-start))/CLOCKS_PER_SEC);

                    if (residual_x_c<TOL_X_C)
                    {
                        printf("residual_x_c at iteration %d, outer, %d, time, %d, = %12.5e\n",inner,outer,time,residual_x_c);
                        printf("Solution converged!\n");
                        // resuc=fopen("res_xc.txt","a");
                        // fprintf(resuc,"%s %d,outer, %d, time, %d, = % 12.5e\n",statement_xc,inner,outer,time,residual_x_c);
                        // fclose(resuc);
                        break;
                    }

                
                    else if (inner == ITERATION_X)
                    {
                        // resuc=fopen("res_xc.txt","a");
                        // fprintf(resuc,"%s , %d, outer, %d, time, %d, = % 12.5e\n",statement_xc,inner,outer,time,residual_x_c);
                        // fclose(resuc);
                        printf("residual_x_c at iteration %d at outer %d at time %d = % 12.5e\n",inner,outer,time,residual_x_c);
                    }   
                
                }


        for (i=0;i<NX;i++)
        {
            free(uc_aw[i]);
            free(uc_as[i]);
            free(uc_an[i]);
            free(uc_ae[i]);
        }
        free(uc_aw);
        free(uc_as);
        free(uc_an);
        free(uc_ae);

        ///// Stone's ILU Solver for y-momentum

            for (i=0;i<NX;i++)
            {
                for (j=1;j<NY;j++)
                {
                    m = (j-1)*NX + i;

                    Eyc[m] = ( vc_ap[i][j]*vd_ap[i][j] + B_ky[i][j]*( vc_ap[i][j] + vd_ap[i][j] ) )/UV_COR;
                    Fyc[m] = vc_ae[i][j]*( vd_ap[i][j] + B_ky[i][j] );
                    Dyc[m] = vc_aw[i][j]*( vd_ap[i][j] + B_ky[i][j] );
                    Byc[m] = vc_as[i][j]*( vd_ap[i][j] + B_ky[i][j] );
                    Hyc[m] = vc_an[i][j]*( vd_ap[i][j] + B_ky[i][j] );
                    Vc[m] = v_c[i][j];

                    Eyd[m] = ( vd_ap[i][j]*vc_ap[i][j] + B_ky[i][j]*( vc_ap[i][j] + vd_ap[i][j] ) )/UV_COR;
                    Fyd[m] = vd_ae[i][j]*( vc_ap[i][j] + B_ky[i][j] );
                    Dyd[m] = vd_aw[i][j]*( vc_ap[i][j] + B_ky[i][j] );
                    Byd[m] = vd_as[i][j]*( vc_ap[i][j] + B_ky[i][j] );
                    Hyd[m] = vd_an[i][j]*( vc_ap[i][j] + B_ky[i][j] );
                    Vd[m] = v_d[i][j];
                }
            }


            for (i=1;i<NX-1;i++)
            {
                for (j=1;j<NY;j++)
                {
                    m = (j-1)*NX + i;

                    Qycmod[m] = Qyc[m]*( vd_ap[i][j] + B_ky[i][j] ) + B_ky[i][j]*Qyd[m] - B_ky[i][j]*( vd_ae[i][j]*v_d[i+1][j] + vd_aw[i][j]*v_d[i-1][j] + vd_an[i][j]*v_d[i][j+1] + vd_as[i][j]*v_d[i][j-1] ) + (1-UV_COR)*Eyc[m]*v_c[i][j];

                    Qydmod[m] = Qyd[m]*( vc_ap[i][j] + B_ky[i][j] ) + B_ky[i][j]*Qyc[m] - B_ky[i][j]*( vc_ae[i][j]*v_c[i+1][j] + vc_aw[i][j]*v_c[i-1][j] + vc_an[i][j]*v_c[i][j+1] + vc_as[i][j]*v_c[i][j-1] ) + (1-UV_COR)*Eyd[m]*v_d[i][j];
                }
            }

            i = 0;

            for (j=1;j<NY;j++)
            {
                m = (j-1)*NX + i;

                Qycmod[m] = Qyc[m]*( vd_ap[i][j] + B_ky[i][j] ) + B_ky[i][j]*Qyd[m] - B_ky[i][j]*( vd_ae[i][j]*v_d[i+1][j] + vd_an[i][j]*v_d[i][j+1] + vd_as[i][j]*v_d[i][j-1] ) + (1-UV_COR)*Eyc[m]*v_c[i][j];

                Qydmod[m] = Qyd[m]*( vc_ap[i][j] + B_ky[i][j] ) + B_ky[i][j]*Qyc[m] - B_ky[i][j]*( vc_ae[i][j]*v_c[i+1][j] + vc_an[i][j]*v_c[i][j+1] + vc_as[i][j]*v_c[i][j-1] ) + (1-UV_COR)*Eyd[m]*v_d[i][j];
            }

            i = NX-1;

            for (j=1;j<NY;j++)
            {
                m = (j-1)*NX + i;

                Qycmod[m] = Qyc[m]*( vd_ap[i][j] + B_ky[i][j] ) + B_ky[i][j]*Qyd[m] - B_ky[i][j]*( vd_aw[i][j]*v_d[i-1][j] + vd_an[i][j]*v_d[i][j+1] + vd_as[i][j]*v_d[i][j-1] ) + (1-UV_COR)*Eyc[m]*v_c[i][j];

                Qydmod[m] = Qyd[m]*( vc_ap[i][j] + B_ky[i][j] ) + B_ky[i][j]*Qyc[m] - B_ky[i][j]*( vc_aw[i][j]*v_c[i-1][j] + vc_an[i][j]*v_c[i][j+1] + vc_as[i][j]*v_c[i][j-1] ) + (1-UV_COR)*Eyd[m]*v_d[i][j];
            }


            coefficient_matrix(Eyc, Fyc, Dyc, Byc, Hyc, fyc, dyc, byc, eyc, cyc, Ky, NX, ALPHA_XY);

            //// Inner loop begins

            // start  = clock();

            for (inner=1;inner<=ITERATION_Y;inner++)
            {

                residual_y_c = SIP_Solver_y (Eyc, Fyc, Dyc, Byc, Hyc, Qycmod, fyc, dyc, byc, eyc, cyc, resyc, vec_Yyc, vec_deltayc, Ky, NX, inner, Vc, v_c, Eyc);

                if (residual_y_c<TOL_Y_C)
                {
                    printf("residual_y_c at iteration %d, outer, %d, time, %d, = %12.5e\n",inner,outer,time,residual_y_c);
                    printf("Solution converged!\n");
                    // resvc=fopen("res_yc.txt","a");
                    // fprintf(resvc,"%s, %d, time, %d, = % 12.5e\n",statement_yc,outer,time,residual_y_c);
                    // fclose(resvc);
                    break;
                }

            
                if (inner == ITERATION_Y)
                {
                    // resvc=fopen("res_yc.txt","a");
                    // fprintf(resvc,"%s , %d, outer, %d, time, %d, = % 12.5e\n",statement_yc,inner,outer,time,residual_y_c);
                    // fclose(resvc);
                    printf("residual_y_c at iteration %d at outer %d at time %d = % 12.5e\n",inner,outer,time,residual_y_c);
                }   
                
            }

            // end = clock();

            // printf("%lf\n", ((double)(end-start))/CLOCKS_PER_SEC);

        for (i=0;i<NX;i++)
        {
            free(vc_aw[i]);
            free(vc_ae[i]);
            free(vc_as[i]);
            free(vc_an[i]);
        }
        free(vc_aw);
        free(vc_ae);
        free(vc_as);
        free(vc_an);

        //// x-momentum solver for dispersed phase

        coefficient_matrix(Exd, Fxd, Dxd, Bxd, Hxd, fxd, dxd, bxd, exd, cxd, Kx, nx, ALPHA_XY);

        //// Inner iteration loop begins

        for (inner=1;inner<=ITERATION_X;inner++)
        {

            residual_x_d = SIP_Solver_x (Exd, Fxd, Dxd, Bxd, Hxd, Qxdmod, fxd, dxd, bxd, exd, cxd, resxd, vec_Yxd, vec_deltaxd, Kx, nx, inner, Ud, u_d, Exd);

            if (residual_x_d<TOL_X_D)
            {
                printf("residual_x_d at iteration %d, outer, %d, time, %d, = %12.5e\n",inner,outer,time,residual_x_d);
                printf("Solution converged!\n");
                // resud=fopen("res_xd.txt","a");
                // fprintf(resud,"%s %d,outer, %d, time, %d, = % 12.5e\n",statement_xd,inner,outer,time,residual_x_d);
                // fclose(resud);
                break;
            }

        
            else if (inner == ITERATION_X)
            {
                // resud=fopen("res_xd.txt","a");
                // fprintf(resud,"%s , %d, outer, %d, time, %d, = % 12.5e\n",statement_xd,inner,outer,time,residual_x_d);
                // fclose(resud);
                printf("residual_x_d at iteration %d at outer %d at time %d = % 12.5e\n",inner,outer,time,residual_x_d);
            }   
            
        }


        for (i=0;i<NX;i++)
        {
            free(ud_aw[i]);
            free(ud_as[i]);
            free(ud_an[i]);
            free(ud_ae[i]);
        }
        free(ud_aw);
        free(ud_as);
        free(ud_an);
        free(ud_ae);

        
        //// y-momentum solver for dispersed phase

            coefficient_matrix(Eyd, Fyd, Dyd, Byd, Hyd, fyd, dyd, byd, eyd, cyd, Ky, NX, ALPHA_XY);

            //// Inner loop begins

            for (inner=1;inner<=ITERATION_Y;inner++)
            {

                residual_y_d = SIP_Solver_y (Eyd, Fyd, Dyd, Byd, Hyd, Qydmod, fyd, dyd, byd, eyd, cyd, resyd, vec_Yyd, vec_deltayd, Ky, NX, inner, Vd, v_d, Eyd);

                if (residual_y_d<TOL_Y_D)
                {
                    printf("residual_y_d at iteration %d, outer, %d, time, %d, = %12.5e\n",inner,outer,time,residual_y_d);
                    printf("Solution converged!\n");
                    // resvd=fopen("res_yd.txt","a");
                    // fprintf(resvd,"%s, %d, time, %d, = % 12.5e\n",statement_yd,outer,time,residual_y_d);
                    // fclose(resvd);
                    break;
                }

            
                if (inner == ITERATION_Y)
                {
                    // resvd=fopen("res_yd.txt","a");
                    // fprintf(resvd,"%s , %d, outer, %d, time, %d, = % 12.5e\n",statement_yd,inner,outer,time,residual_y_d);
                    // fclose(resvd);
                    printf("residual_y_d at iteration %d at outer %d at time %d = % 12.5e\n",inner,outer,time,residual_y_d);
                }   
                
            }

        for (i=0;i<NX;i++)
        {
            free(vd_aw[i]);
            free(vd_ae[i]);
            free(vd_as[i]);
            free(vd_an[i]);
        }
        free(vd_aw);
        free(vd_ae);
        free(vd_as);
        free(vd_an);


        double ** pp = calloc(NX,sizeof(double *));
        double ** pc_ap = calloc(NX,sizeof(double *));
        double ** pc_an = calloc(NX,sizeof(double *));
        double ** pc_as = calloc(NX,sizeof(double *));
        double ** pc_ae = calloc(NX,sizeof(double *));
        double ** pc_aw = calloc(NX,sizeof(double *));
        double ** mdot_c = calloc(NX,sizeof(double *));
        double ** pd_ap = calloc(NX,sizeof(double *));
        double ** pd_an = calloc(NX,sizeof(double *));
        double ** pd_as = calloc(NX,sizeof(double *));
        double ** pd_ae = calloc(NX,sizeof(double *));
        double ** pd_aw = calloc(NX,sizeof(double *));
        double ** mdot_d = calloc(NX,sizeof(double *));
        for (i=0;i<NX;i++)
        {
            pp[i] = calloc(NY,sizeof(double));
            pc_ap[i] = calloc(NY,sizeof(double));
            pc_an[i] = calloc(NY,sizeof(double));
            pc_as[i] = calloc(NY,sizeof(double));
            pc_ae[i] = calloc(NY,sizeof(double));
            pc_aw[i] = calloc(NY,sizeof(double));
            mdot_c[i] = calloc(NY,sizeof(double));
            pd_ap[i] = calloc(NY,sizeof(double));
            pd_an[i] = calloc(NY,sizeof(double));
            pd_as[i] = calloc(NY,sizeof(double));
            pd_ae[i] = calloc(NY,sizeof(double));
            pd_aw[i] = calloc(NY,sizeof(double));
            mdot_d[i] = calloc(NY,sizeof(double));
        }

        //// Interior Cells

        for (i=1;i<NX-1;i++)
        {
            for (j=1;j<NY-1;j++)
            {
                alpe_c = (alpha_c[i][j] + alpha_c[i+1][j])/2;
                alpw_c = (alpha_c[i][j] + alpha_c[i-1][j])/2;
                alpn_c = (alpha_c[i][j] + alpha_c[i][j+1])/2;
                alps_c = (alpha_c[i][j] + alpha_c[i][j-1])/2;

                alpe_d = (alpha_d[i][j] + alpha_d[i+1][j])/2;
                alpw_d = (alpha_d[i][j] + alpha_d[i-1][j])/2;
                alpn_d = (alpha_d[i][j] + alpha_d[i][j+1])/2;
                alps_d = (alpha_d[i][j] + alpha_d[i][j-1])/2;

                pc_ae[i][j] = -alpe_c*PCC_x(i+1,j,uc_ap,ud_ap,alpe_c,alpe_d,B_kx,1);
                pc_aw[i][j] = -alpw_c*(PCC_x(i,j,uc_ap,ud_ap,alpw_c,alpw_d,B_kx,1));
                pc_an[i][j] = -alpn_c*(PCC_y(i,j+1,vc_ap,vd_ap,alpn_c,alpn_d,B_ky,1));
                pc_as[i][j] = -alps_c*(PCC_y(i,j,vc_ap,vd_ap,alps_c,alps_d,B_ky,1));
                pc_ap[i][j] = -(pc_ae[i][j] + pc_aw[i][j] + pc_an[i][j] + pc_as[i][j]);
                mdot_c[i][j] = (u_c[i+1][j]*alpe_c - u_c[i][j]*alpw_c)*DY + (v_c[i][j+1]*alpn_c - v_c[i][j]*alps_c)*DX;
            }
        }

        //// Left Boundary

        i=0;

        for (j=1;j<NY-1;j++)
        {
            alpe_c = (alpha_c[i][j] + alpha_c[i+1][j])/2;
            alpn_c = (alpha_c[i][j] + alpha_c[i][j+1])/2;
            alps_c = (alpha_c[i][j] + alpha_c[i][j-1])/2;

            alpe_d = (alpha_d[i][j] + alpha_d[i+1][j])/2;
            alpn_d = (alpha_d[i][j] + alpha_d[i][j+1])/2;
            alps_d = (alpha_d[i][j] + alpha_d[i][j-1])/2;

            pc_ae[i][j] = -alpe_c*(PCC_x(i+1,j,uc_ap,ud_ap,alpe_c,alpe_d,B_kx,1));
            pc_aw[i][j] = 0;
            pc_an[i][j] = -alpn_c*(PCC_y(i,j+1,vc_ap,vd_ap,alpn_c,alpn_d,B_ky,1));
            pc_as[i][j] = -alps_c*(PCC_y(i,j,vc_ap,vd_ap,alps_c,alps_d,B_ky,1));
            pc_ap[i][j] = -(pc_ae[i][j] + pc_aw[i][j] + pc_an[i][j] + pc_as[i][j]);
            mdot_c[i][j] = (u_c[i+1][j]*alpe_c)*DY + (v_c[i][j+1]*alpn_c - v_c[i][j]*alps_c)*DX;
        }

        //// Right Boundary

        i=NX-1;

        for (j=1;j<NY-1;j++)
        {
            alpw_c = (alpha_c[i][j] + alpha_c[i-1][j])/2;
            alpn_c = (alpha_c[i][j] + alpha_c[i][j+1])/2;
            alps_c = (alpha_c[i][j] + alpha_c[i][j-1])/2;

            alpw_d = (alpha_d[i][j] + alpha_d[i-1][j])/2;
            alpn_d = (alpha_d[i][j] + alpha_d[i][j+1])/2;
            alps_d = (alpha_d[i][j] + alpha_d[i][j-1])/2;

            pc_ae[i][j] = 0;
            pc_aw[i][j] = -alpw_c*(PCC_x(i,j,uc_ap,ud_ap,alpw_c,alpw_d,B_kx,1));
            pc_an[i][j] = -alpn_c*(PCC_y(i,j+1,vc_ap,vd_ap,alpn_c,alpn_d,B_ky,1));
            pc_as[i][j] = -alps_c*(PCC_y(i,j,vc_ap,vd_ap,alps_c,alps_d,B_ky,1));
            pc_ap[i][j] = -(pc_ae[i][j] + pc_aw[i][j] + pc_an[i][j] + pc_as[i][j]);
            mdot_c[i][j] = ( - u_c[i][j]*alpw_c)*DY + (v_c[i][j+1]*alpn_c - v_c[i][j]*alps_c)*DX;
        }

        //// Bottom Bpundary

        j=0;

        for (i=1;i<NX-1;i++)
        {
            alpe_c = (alpha_c[i][j] + alpha_c[i+1][j])/2;
            alpw_c = (alpha_c[i][j] + alpha_c[i-1][j])/2;
            alpn_c = (alpha_c[i][j] + alpha_c[i][j+1])/2;
            alps_c = alphac_inlet[i];

            alpe_d = (alpha_d[i][j] + alpha_d[i+1][j])/2;
            alpw_d = (alpha_d[i][j] + alpha_d[i-1][j])/2;
            alpn_d = (alpha_d[i][j] + alpha_d[i][j+1])/2;
            alps_d = alphad_inlet[i];

            pc_ae[i][j] = -alpe_c*(PCC_x(i+1,j,uc_ap,ud_ap,alpe_c,alpe_d,B_kx,1));
            pc_aw[i][j] = -alpw_c*(PCC_x(i,j,uc_ap,ud_ap,alpw_c,alpw_d,B_kx,1));
            pc_an[i][j] = -alpn_c*(PCC_y(i,j+1,vc_ap,vd_ap,alpn_c,alpn_d,B_ky,1));
            pc_as[i][j] = 0;
            pc_ap[i][j] = -(pc_ae[i][j] + pc_aw[i][j] + pc_an[i][j] + pc_as[i][j]);
            mdot_c[i][j] = (u_c[i+1][j]*alpe_c - u_c[i][j]*alpw_c)*DY + (v_c[i][j+1]*alpn_c - v_c[i][j]*alps_c)*DX;
        }

        //// Top Boundary

        j=NY-1;

        for (i=1;i<NX-1;i++)
        {
            alpe_c = (alpha_c[i][j] + alpha_c[i+1][j])/2;
            alpw_c = (alpha_c[i][j] + alpha_c[i-1][j])/2;
            alpn_c = alpha_c[i][j];
            alps_c = (alpha_c[i][j] + alpha_c[i][j-1])/2;

            alpe_d = (alpha_d[i][j] + alpha_d[i+1][j])/2;
            alpw_d = (alpha_d[i][j] + alpha_d[i-1][j])/2;
            alpn_d = alpha_d[i][j];
            alps_d = (alpha_d[i][j] + alpha_d[i][j-1])/2;

            pc_ae[i][j] = -alpe_c*(PCC_x(i+1,j,uc_ap,ud_ap,alpe_c,alpe_d,B_kx,1));
            pc_aw[i][j] = -alpw_c*(PCC_x(i,j,uc_ap,ud_ap,alpw_c,alpw_d,B_kx,1));
            pc_an[i][j] = 0;
            pc_as[i][j] = -alps_c*(PCC_y(i,j,vc_ap,vd_ap,alps_c,alps_d,B_ky,1));
            pc_ap[i][j] = -(pc_ae[i][j] + pc_aw[i][j] + pc_an[i][j] + pc_as[i][j]);
            mdot_c[i][j] = (u_c[i+1][j]*alpe_c - u_c[i][j]*alpw_c)*DY + (v_c[i][j+1]*alpn_c - v_c[i][j]*alps_c)*DX;
        }

        //// Bottom-Left

        i=0;
        j=0;

        alpe_c = (alpha_c[i][j] + alpha_c[i+1][j])/2;
        alpn_c = (alpha_c[i][j] + alpha_c[i][j+1])/2;
        alps_c = alphac_inlet[i];

        alpe_d = (alpha_d[i][j] + alpha_d[i+1][j])/2;
        alpn_d = (alpha_d[i][j] + alpha_d[i][j+1])/2;
        alps_d = alphad_inlet[i];

        pc_ae[i][j] = -alpe_c*(PCC_x(i+1,j,uc_ap,ud_ap,alpe_c,alpe_d,B_kx,1));
        pc_aw[i][j] = 0;
        pc_an[i][j] = -alpn_c*(PCC_y(i,j+1,vc_ap,vd_ap,alpn_c,alpn_d,B_ky,1));
        pc_as[i][j] = 0;
        pc_ap[i][j] = -(pc_ae[i][j] + pc_aw[i][j] + pc_an[i][j] + pc_as[i][j]);
        mdot_c[i][j] = (u_c[i+1][j]*alpe_c)*DY + (v_c[i][j+1]*alpn_c - v_c[i][j]*alps_c)*DX;


        //// Top-Left

        j=NY-1;

        alpe_c = (alpha_c[i][j] + alpha_c[i+1][j])/2;
        alpn_c = alpha_c[i][j];
        alps_c = (alpha_c[i][j] + alpha_c[i][j-1])/2;

        alpe_d = (alpha_d[i][j] + alpha_d[i+1][j])/2;
        alpn_d = alpha_d[i][j];
        alps_d = (alpha_d[i][j] + alpha_d[i][j-1])/2;

        pc_ae[i][j] = -alpe_c*(PCC_x(i+1,j,uc_ap,ud_ap,alpe_c,alpe_d,B_kx,1));
        pc_aw[i][j] = 0;
        pc_an[i][j] = 0;
        pc_as[i][j] = -alps_c*(PCC_y(i,j,vc_ap,vd_ap,alps_c,alps_d,B_ky,1));
        pc_ap[i][j] = -(pc_ae[i][j] + pc_aw[i][j] + pc_an[i][j] + pc_as[i][j]);
        mdot_c[i][j] = (u_c[i+1][j]*alpe_c)*DY + (v_c[i][j+1]*alpn_c - v_c[i][j]*alps_c)*DX;

        //// Bottom-Right

        i=NX-1;
        j=0;

        alpw_c = (alpha_c[i][j] + alpha_c[i-1][j])/2;
        alpn_c = (alpha_c[i][j] + alpha_c[i][j+1])/2;
        alps_c = alphac_inlet[i];

        alpw_d = (alpha_d[i][j] + alpha_d[i-1][j])/2;
        alpn_d = (alpha_d[i][j] + alpha_d[i][j+1])/2;
        alps_d = alphad_inlet[i];

        pc_ae[i][j] = 0;
        pc_aw[i][j] = -alpw_c*(PCC_x(i,j,uc_ap,ud_ap,alpw_c,alpw_d,B_kx,1));
        pc_an[i][j] = -alpn_c*(PCC_y(i,j+1,vc_ap,vd_ap,alpn_c,alpn_d,B_ky,1));
        pc_as[i][j] = 0;
        pc_ap[i][j] = -(pc_ae[i][j] + pc_aw[i][j] + pc_an[i][j] + pc_as[i][j]);
        mdot_c[i][j] = ( - u_c[i][j]*alpw_c)*DY + (v_c[i][j+1]*alpn_c - v_c[i][j]*alps_c)*DX;

        //// Top-Right

        j=NY-1;

        alpw_c = (alpha_c[i][j] + alpha_c[i-1][j])/2;
        alpn_c = alpha_c[i][j];
        alps_c = (alpha_c[i][j] + alpha_c[i][j-1])/2;

        alpw_d = (alpha_d[i][j] + alpha_d[i-1][j])/2;
        alpn_d = alpha_d[i][j];
        alps_d = (alpha_d[i][j] + alpha_d[i][j-1])/2;

        pc_ae[i][j] = 0;
        pc_aw[i][j] = -alpw_c*(PCC_x(i,j,uc_ap,ud_ap,alpw_c,alpw_d,B_kx,1));
        pc_an[i][j] = 0;
        pc_as[i][j] = -alps_c*(PCC_y(i,j,vc_ap,vd_ap,alps_c,alps_d,B_ky,1));
        pc_ap[i][j] = -(pc_ae[i][j] + pc_aw[i][j] + pc_an[i][j] + pc_as[i][j]);
        mdot_c[i][j] = ( - u_c[i][j]*alpw_c)*DY + (v_c[i][j+1]*alpn_c - v_c[i][j]*alps_c)*DX;


        //// Interior Cells

        for (i=1;i<NX-1;i++)
        {
            for (j=1;j<NY-1;j++)
            {
                alpe_d = (alpha_d[i][j] + alpha_d[i+1][j])/2;
                alpw_d = (alpha_d[i][j] + alpha_d[i-1][j])/2;
                alpn_d = (alpha_d[i][j] + alpha_d[i][j+1])/2;
                alps_d = (alpha_d[i][j] + alpha_d[i][j-1])/2;

                alpe_c = (alpha_c[i][j] + alpha_c[i+1][j])/2;
                alpw_c = (alpha_c[i][j] + alpha_c[i-1][j])/2;
                alpn_c = (alpha_c[i][j] + alpha_c[i][j+1])/2;
                alps_c = (alpha_c[i][j] + alpha_c[i][j-1])/2;

                pd_ae[i][j] = -alpe_d*PCC_x(i+1,j,uc_ap,ud_ap,alpe_c,alpe_d,B_kx,0);
                pd_aw[i][j] = -alpw_d*(PCC_x(i,j,uc_ap,ud_ap,alpw_c,alpw_d,B_kx,0));
                pd_an[i][j] = -alpn_d*(PCC_y(i,j+1,vc_ap,vd_ap,alpn_c,alpn_d,B_ky,0));
                pd_as[i][j] = -alps_d*(PCC_y(i,j,vc_ap,vd_ap,alps_c,alps_d,B_ky,0));
                pd_ap[i][j] = -(pd_ae[i][j] + pd_aw[i][j] + pd_an[i][j] + pd_as[i][j]);
                mdot_d[i][j] = (u_d[i+1][j]*alpe_d - u_d[i][j]*alpw_d)*DY + (v_d[i][j+1]*alpn_d - v_d[i][j]*alps_d)*DX;
            }
        }

        //// Left Boundary

        i=0;

        for (j=1;j<NY-1;j++)
        {
            alpe_d = (alpha_d[i][j] + alpha_d[i+1][j])/2;
            alpn_d = (alpha_d[i][j] + alpha_d[i][j+1])/2;
            alps_d = (alpha_d[i][j] + alpha_d[i][j-1])/2;

            alpe_c = (alpha_c[i][j] + alpha_c[i+1][j])/2;
            alpn_c = (alpha_c[i][j] + alpha_c[i][j+1])/2;
            alps_c = (alpha_c[i][j] + alpha_c[i][j-1])/2;

            pd_ae[i][j] = -alpe_d*(PCC_x(i+1,j,uc_ap,ud_ap,alpe_c,alpe_d,B_kx,0));
            pd_aw[i][j] = 0;
            pd_an[i][j] = -alpn_d*(PCC_y(i,j+1,vc_ap,vd_ap,alpn_c,alpn_d,B_ky,0));
            pd_as[i][j] = -alps_d*(PCC_y(i,j,vc_ap,vd_ap,alps_c,alps_d,B_ky,0));
            pd_ap[i][j] = -(pd_ae[i][j] + pd_aw[i][j] + pd_an[i][j] + pd_as[i][j]);
            mdot_d[i][j] = (u_d[i+1][j]*alpe_d)*DY + (v_d[i][j+1]*alpn_d - v_d[i][j]*alps_d)*DX;
        }

        //// Right Boundary

        i=NX-1;

        for (j=1;j<NY-1;j++)
        {
            alpw_d = (alpha_d[i][j] + alpha_d[i-1][j])/2;
            alpn_d = (alpha_d[i][j] + alpha_d[i][j+1])/2;
            alps_d = (alpha_d[i][j] + alpha_d[i][j-1])/2;

            alpw_c = (alpha_c[i][j] + alpha_c[i-1][j])/2;
            alpn_c = (alpha_c[i][j] + alpha_c[i][j+1])/2;
            alps_c = (alpha_c[i][j] + alpha_c[i][j-1])/2;

            pd_ae[i][j] = 0;
            pd_aw[i][j] = -alpw_d*(PCC_x(i,j,uc_ap,ud_ap,alpw_c,alpw_d,B_kx,0));
            pd_an[i][j] = -alpn_d*(PCC_y(i,j+1,vc_ap,vd_ap,alpn_c,alpn_d,B_ky,0));
            pd_as[i][j] = -alps_d*(PCC_y(i,j,vc_ap,vd_ap,alps_c,alps_d,B_ky,0));
            pd_ap[i][j] = -(pd_ae[i][j] + pd_aw[i][j] + pd_an[i][j] + pd_as[i][j]);
            mdot_d[i][j] = ( - u_d[i][j]*alpw_d)*DY + (v_d[i][j+1]*alpn_d - v_d[i][j]*alps_d)*DX;
        }

        //// Bottom Bpundary

        j=0;

        for (i=1;i<NX-1;i++)
        {
            alpe_d = (alpha_d[i][j] + alpha_d[i+1][j])/2;
            alpw_d = (alpha_d[i][j] + alpha_d[i-1][j])/2;
            alpn_d = (alpha_d[i][j] + alpha_d[i][j+1])/2;
            alps_d = alphad_inlet[i];

            alpe_c = (alpha_c[i][j] + alpha_c[i+1][j])/2;
            alpw_c = (alpha_c[i][j] + alpha_c[i-1][j])/2;
            alpn_c = (alpha_c[i][j] + alpha_c[i][j+1])/2;
            alps_c = alphac_inlet[i];

            pd_ae[i][j] = -alpe_d*(PCC_x(i+1,j,uc_ap,ud_ap,alpe_c,alpe_d,B_kx,0));
            pd_aw[i][j] = -alpw_d*(PCC_x(i,j,uc_ap,ud_ap,alpw_c,alpw_d,B_kx,0));
            pd_an[i][j] = -alpn_d*(PCC_y(i,j+1,vc_ap,vd_ap,alpn_c,alpn_d,B_ky,0));
            pd_as[i][j] = 0;
            pd_ap[i][j] = -(pd_ae[i][j] + pd_aw[i][j] + pd_an[i][j] + pd_as[i][j]);
            mdot_d[i][j] = (u_d[i+1][j]*alpe_d - u_d[i][j]*alpw_d)*DY + (v_d[i][j+1]*alpn_d - v_d[i][j]*alps_d)*DX;
        }

        //// Top Boundary

        j=NY-1;

        for (i=1;i<NX-1;i++)
        {
            alpe_d = (alpha_d[i][j] + alpha_d[i+1][j])/2;
            alpw_d = (alpha_d[i][j] + alpha_d[i-1][j])/2;
            alpn_d = alpha_d[i][j];
            alps_d = (alpha_d[i][j] + alpha_d[i][j-1])/2;

            alpe_c = (alpha_c[i][j] + alpha_c[i+1][j])/2;
            alpw_c = (alpha_c[i][j] + alpha_c[i-1][j])/2;
            alpn_c = alpha_c[i][j];
            alps_c = (alpha_c[i][j] + alpha_c[i][j-1])/2;

            pd_ae[i][j] = -alpe_d*(PCC_x(i+1,j,uc_ap,ud_ap,alpe_c,alpe_d,B_kx,0));
            pd_aw[i][j] = -alpw_d*(PCC_x(i,j,uc_ap,ud_ap,alpw_c,alpw_d,B_kx,0));
            pd_an[i][j] = 0;
            pd_as[i][j] = -alps_d*(PCC_y(i,j,vc_ap,vd_ap,alps_c,alps_d,B_ky,0));
            pd_ap[i][j] = -(pd_ae[i][j] + pd_aw[i][j] + pd_an[i][j] + pd_as[i][j]);
            mdot_d[i][j] = (u_d[i+1][j]*alpe_d - u_d[i][j]*alpw_d)*DY + (v_d[i][j+1]*alpn_d - v_d[i][j]*alps_d)*DX + alpha_d[i][j]*DX*v_d[i][j];
        }

        //// Bottom-Left

        i=0;
        j=0;

        alpe_d = (alpha_d[i][j] + alpha_d[i+1][j])/2;
        alpn_d = (alpha_d[i][j] + alpha_d[i][j+1])/2;
        alps_d = alphad_inlet[i];

        alpe_c = (alpha_c[i][j] + alpha_c[i+1][j])/2;
        alpn_c = (alpha_c[i][j] + alpha_c[i][j+1])/2;
        alps_c = alphac_inlet[i];

        pd_ae[i][j] = -alpe_d*(PCC_x(i+1,j,uc_ap,ud_ap,alpe_c,alpe_d,B_kx,0));
        pd_aw[i][j] = 0;
        pd_an[i][j] = -alpn_d*(PCC_y(i,j+1,vc_ap,vd_ap,alpn_c,alpn_d,B_ky,0));
        pd_as[i][j] = 0;
        pd_ap[i][j] = -(pd_ae[i][j] + pd_aw[i][j] + pd_an[i][j] + pd_as[i][j]);
        mdot_d[i][j] = (u_d[i+1][j]*alpe_d)*DY + (v_d[i][j+1]*alpn_d - v_d[i][j]*alps_d)*DX;

        //// Top-Left

        j=NY-1;

        alpe_d = (alpha_d[i][j] + alpha_d[i+1][j])/2;
        alpn_d = alpha_d[i][j];
        alps_d = (alpha_d[i][j] + alpha_d[i][j-1])/2;

        alpe_c = (alpha_c[i][j] + alpha_c[i+1][j])/2;
        alpn_c = alpha_c[i][j];
        alps_c = (alpha_c[i][j] + alpha_c[i][j-1])/2;

        pd_ae[i][j] = -alpe_d*(PCC_x(i+1,j,uc_ap,ud_ap,alpe_c,alpe_d,B_kx,0));
        pd_aw[i][j] = 0;
        pd_an[i][j] = 0;
        pd_as[i][j] = -alps_d*(PCC_y(i,j,vc_ap,vd_ap,alps_c,alps_d,B_ky,0));
        pd_ap[i][j] = -(pd_ae[i][j] + pd_aw[i][j] + pd_an[i][j] + pd_as[i][j]);
        mdot_d[i][j] = (u_d[i+1][j]*alpe_d)*DY + (v_d[i][j+1]*alpn_d - v_d[i][j]*alps_d)*DX + alpha_d[i][j]*DX*v_d[i][j];

        //// Bottom-Right

        i=NX-1;
        j=0;

        alpw_d = (alpha_d[i][j] + alpha_d[i-1][j])/2;
        alpn_d = (alpha_d[i][j] + alpha_d[i][j+1])/2;
        alps_d = alphad_inlet[i];

        alpw_c = (alpha_c[i][j] + alpha_c[i-1][j])/2;
        alpn_c = (alpha_c[i][j] + alpha_c[i][j+1])/2;
        alps_c = alphac_inlet[i];

        pd_ae[i][j] = 0;
        pd_aw[i][j] = -alpw_d*(PCC_x(i,j,uc_ap,ud_ap,alpw_c,alpw_d,B_kx,0));
        pd_an[i][j] = -alpn_d*(PCC_y(i,j+1,vc_ap,vd_ap,alpn_c,alpn_d,B_ky,0));
        pd_as[i][j] = 0;
        pd_ap[i][j] = -(pd_ae[i][j] + pd_aw[i][j] + pd_an[i][j] + pd_as[i][j]);
        mdot_d[i][j] = (- u_d[i][j]*alpw_d)*DY + (v_d[i][j+1]*alpn_d - v_d[i][j]*alps_d)*DX;

        //// Top-Right

        j=NY-1;

        alpw_d = (alpha_d[i][j] + alpha_d[i-1][j])/2;
        alpn_d = alpha_d[i][j];
        alps_d = (alpha_d[i][j] + alpha_d[i][j-1])/2;

        alpw_c = (alpha_c[i][j] + alpha_c[i-1][j])/2;
        alpn_c = alpha_c[i][j];
        alps_c = (alpha_c[i][j] + alpha_c[i][j-1])/2;

        pd_ae[i][j] = 0;
        pd_aw[i][j] = -alpw_d*(PCC_x(i,j,uc_ap,ud_ap,alpw_c,alpw_d,B_kx,0));
        pd_an[i][j] = 0;
        pd_as[i][j] = -alps_d*(PCC_y(i,j,vc_ap,vd_ap,alps_c,alps_d,B_ky,0));
        pd_ap[i][j] = -(pd_ae[i][j] + pd_aw[i][j] + pd_an[i][j] + pd_as[i][j]);
        mdot_d[i][j] = ( - u_d[i][j]*alpw_d)*DY + (v_d[i][j+1]*alpn_d - v_d[i][j]*alps_d)*DX + alpha_d[i][j]*DX*v_d[i][j];


        ///// Co-efficients end
     
        //// Pressure-Correction equation solver

        for (i=0;i<NX;i++)
        {
            for (j=0;j<NY;j++)
            {
                m = j*NX + i;

                Ep[m] = (pc_ap[i][j] + pd_ap[i][j]);
                Fp[m] = (pc_ae[i][j] + pd_ae[i][j]);
                Dp[m] = (pc_aw[i][j] + pd_aw[i][j]);
                Bp[m] = (pc_as[i][j] + pd_as[i][j]);
                Hp[m] = (pc_an[i][j] + pd_an[i][j]);
                Qp[m] = - ( mdot_c[i][j] + mdot_d[i][j] );// + ( alpha_c[i][j] - alphan_c[i][j] )*DX2*DT + ( alpha_d[i][j] - alphan_d[i][j] )*DX2*DT );
                p[m] = pp[i][j];
            }
        }

        coefficient_matrix(Ep, Fp, Dp, Bp, Hp, fp, dp, bp, ep, cp, Kp, NX, ALPHA_P);

        //// Inner loop begins

        start = clock();

        for (inner=1;inner<=ITERATION_P;inner++)
        {

            R2p = SIP_Solver_Centre (Ep, Fp, Dp, Bp, Hp, Qp, fp, dp, bp, ep, cp, resp, vec_Yp, vec_deltap, Kp, NX, p, pp);

            if ((outer==1)&&(inner==1))
            {
                R_1 = R2p;
            }

            residual_P = R2p/R_1;

            // if(time == 1000)
            // {
            //     resP=fopen("res_p.txt","a");
            //     fprintf(resP,"%s ,%d,%d, time, %d, =, % 12.5e\n",statement_P,outer,inner,time,residual_P);
            //     fclose(resP);
            // }




                if (residual_P<1e-3)
                {
                    printf("residual_p at iteration %d, outer, %d, time, %d, = %12.5e\n",inner,outer,time,residual_P);
                    printf("Solution converged!\n");
                    // resP=fopen("res_p.txt","a");
                    // fprintf(resP,"%s ,%d, time, %d, = % 12.5e\n",statement_P,outer,time,residual_P);
                    // fclose(resP);
                    break;
                }

            
            if (inner == ITERATION_P)
            {
                // resP=fopen("res_p.txt","a");
                // fprintf(resP,"%s , %d, outer, %d, time, %d, = % 12.5e\n",statement_P,inner,outer,time,residual_P);
                // fclose(resP);
                printf("residual_p at iteration %d at outer %d at time %d = % 12.5e\n",inner,outer,time,residual_P);
            }   

        }

        // end = clock();

        // up = ((double)(end-start))/CLOCKS_PER_SEC;
        // // printf("%lf\n", ((double)(end-start))/CLOCKS_PER_SEC);

        // if(time == 1)
        // {
        //     resP=fopen("res_p.txt","a");
        //     fprintf(resP,"%d, =, %lf\n",outer,up);
        //     fclose(resP);
        // }

        
        //// Correct u_c-velocity
        
        for (i=1;i<NX;i++)
        {
            for (j=0;j<NY;j++)
            {
                u_c[i][j] += ( ( pp[i-1][j] - pp[i][j] )/DX )*VC_x(i,j,uc_ap,ud_ap,alpha_c,alpha_d,B_kx,1);
                u_d[i][j] += ( ( pp[i-1][j] - pp[i][j] )/DX )*VC_x(i,j,uc_ap,ud_ap,alpha_c,alpha_d,B_kx,0);
            }
        }
        
        //// Correct v_c-velocity
        
        for (i=0;i<NX;i++)
        {
            for (j=1;j<NY;j++)
            {
                v_c[i][j] += ( ( pp[i][j-1] - pp[i][j] )/DY )*VC_y(i,j,vc_ap,vd_ap,alpha_c,alpha_d,B_ky,1);
                v_d[i][j] += ( ( pp[i][j-1] - pp[i][j] )/DY )*VC_y(i,j,vc_ap,vd_ap,alpha_c,alpha_d,B_ky,0);
            }
        }

        //// Correct Pressure
        
        for (i=0;i<NX;i++)
        {
            for (j=0;j<NY;j++)
            {
                P[i][j] += P_COR*pp[i][j];
            }
        }
        
        for (i=0;i<NX;i++)
        {
            free(uc_ap[i]);
            free(ud_ap[i]);
            free(vc_ap[i]);
            free(vd_ap[i]);
            free(pp[i]);
            free(pc_aw[i]);
            free(pc_ae[i]);
            free(pc_as[i]);
            free(pc_an[i]);
            free(pc_ap[i]);
            free(mdot_c[i]);
            free(pd_aw[i]);
            free(pd_ae[i]);
            free(pd_as[i]);
            free(pd_an[i]);
            free(pd_ap[i]);
            free(mdot_d[i]);
        }
        free(uc_ap);
        free(ud_ap);
        free(vc_ap);
        free(vd_ap);
        free(pp);
        free(pc_aw);
        free(pc_ae);
        free(pc_as);
        free(pc_an);
        free(pc_ap);
        free(mdot_c);
        free(pd_aw);
        free(pd_ae);
        free(pd_as);
        free(pd_an);
        free(pd_ap);
        free(mdot_d);

        // Phasic-Continuity equation for dispersed phase

        double ** alphad_p = calloc(NX,sizeof(double *));
        double ** alphad_e = calloc(NX,sizeof(double *));
        double ** alphad_w = calloc(NX,sizeof(double *));
        double ** alphad_n = calloc(NX,sizeof(double *));
        double ** alphad_s = calloc(NX,sizeof(double *));
        for (i=0;i<NX;i++)
        {
            alphad_p[i] = calloc(NY,sizeof(double));
            alphad_e[i] = calloc(NY,sizeof(double));
            alphad_w[i] = calloc(NY,sizeof(double));
            alphad_n[i] = calloc(NY,sizeof(double));
            alphad_s[i] = calloc(NY,sizeof(double));
        }
        // Coefficients for Phasic Continuity

        //// Interior Cells

        for (i=1;i<NX-1;i++)
        {
            for (j=1;j<NY-1;j++)
            {
                uE = u_d[i+1][j]*DY;
                uW = u_d[i][j]*DY;
                vN = v_d[i][j+1]*DX;
                vS = v_d[i][j]*DX;

                z_e = uE > 0;
                z_w = uW > 0;
                z_n = vN > 0;
                z_s = vS > 0;

                m = j*NX + i;
                
                alphad_e[i][j] = -( max(-uE,0) );
                alphad_w[i][j] = -( max(uW,0) );
                alphad_n[i][j] = -( max(-vN,0) );
                alphad_s[i][j] = -( max(vS,0) );
                alphad_p[i][j] = -(alphad_e[i][j] + alphad_w[i][j] + alphad_n[i][j] + alphad_s[i][j]) + DX2*DT + uE - uW + vN - vS;
                Qa[m] = (1-ALPHA_COR)*(alphad_p[i][j]/ALPHA_COR)*alpha_d[i][j] + alphan_d[i][j]*DX2*DT + 0.5*uE*( alpha_d[i+1][j] - alpha_d[i][j] )*( (1 - z_e)*psi_neg_x(alpha_d,i+1,j,0,0) - z_e*psi_pos_x(alpha_d,i,j,0,0) ) + 0.5*uW*( alpha_d[i][j] - alpha_d[i-1][j] )*( z_w*psi_pos_x(alpha_d,i-1,j,0,0) - (1 - z_w)*psi_neg_x(alpha_d,i,j,0,0) ) + 0.5*vN*( alpha_d[i][j+1] - alpha_d[i][j] )*( (1 - z_n)*psi_neg_y(alpha_d,i,j+1,0,0) - z_n*psi_pos_y(alpha_d,i,j,0,0,alphad_inlet) ) + 0.5*vS*( alpha_d[i][j] - alpha_d[i][j-1] )*( z_s*psi_pos_y(alpha_d,i,j-1,0,0,alphad_inlet) - (1 - z_s)*psi_neg_y(alpha_d,i,j,0,0) );            
            }
        }

        // Bottom Boundary
        
        j=0;
        
        for (i=1;i<NX-1;i++)
        {
            uE = u_d[i+1][j]*DY;
            uW = u_d[i][j]*DY;
            vN = v_d[i][j+1]*DX;
            vS = v_d[i][j]*DX;

            z_e = uE > 0;
            z_w = uW > 0;
            z_n = vN > 0;
            z_s = vS > 0;

            m = j*NX + i;

            alphad_e[i][j] = -( max(-uE,0) );
            alphad_w[i][j] = -( max(uW,0) );
            alphad_n[i][j] = -( max(-vN,0) );
            alphad_s[i][j] = 0;
            alphad_p[i][j] = -(alphad_e[i][j] + alphad_w[i][j] + alphad_n[i][j] + alphad_s[i][j]) + DX2*DT + uE - uW + vN - vS;
            Qa[m] = alphad_inlet[i]*vS + (1-ALPHA_COR)*(alphad_p[i][j]/ALPHA_COR)*alpha_d[i][j] + alphan_d[i][j]*DX2*DT + 0.5*uE*( alpha_d[i+1][j] - alpha_d[i][j] )*( (1 - z_e)*psi_neg_x(alpha_d,i+1,j,0,0) - z_e*psi_pos_x(alpha_d,i,j,0,0) ) + 0.5*uW*( alpha_d[i][j] - alpha_d[i-1][j] )*( z_w*psi_pos_x(alpha_d,i-1,j,0,0) - (1 - z_w)*psi_neg_x(alpha_d,i,j,0,0) ) + 0.5*vN*( alpha_d[i][j+1] - alpha_d[i][j] )*( (1 - z_n)*psi_neg_y(alpha_d,i,j+1,1,0) - z_n*psi_pos_y(alpha_d,i,j,1,0,alphad_inlet) );    
        }

        // Outlet
        
        j=NY-1;
        
        for (i=1;i<NX-1;i++)
        {
            uE = u_d[i+1][j]*DY;
            uW = u_d[i][j]*DY;
            vN = v_d[i][j+1]*DX;
            vS = v_d[i][j]*DX;

            z_e = uE > 0;
            z_w = uW > 0;
            z_n = vN > 0;
            z_s = vS > 0;

            m = j*NX + i;

            alphad_e[i][j] = -( max(-uE,0) );
            alphad_w[i][j] = -( max(uW,0) );
            alphad_n[i][j] = 0;
            alphad_s[i][j] = -( max(vS,0) );
            alphad_p[i][j] = -(alphad_e[i][j] + alphad_w[i][j] + alphad_n[i][j] + alphad_s[i][j]) + DX2*DT + uE - uW + vN - vS + DX*v_d[i][j];
            Qa[m] = (1-ALPHA_COR)*(alphad_p[i][j]/ALPHA_COR)*alpha_d[i][j] + alphan_d[i][j]*DX2*DT + 0.5*uE*( alpha_d[i+1][j] - alpha_d[i][j] )*( (1 - z_e)*psi_neg_x(alpha_d,i+1,j,0,0) - z_e*psi_pos_x(alpha_d,i,j,0,0) ) + 0.5*uW*( alpha_d[i][j] - alpha_d[i-1][j] )*( z_w*psi_pos_x(alpha_d,i-1,j,0,0) - (1 - z_w)*psi_neg_x(alpha_d,i,j,0,0) ) + 0.5*vS*( alpha_d[i][j] - alpha_d[i][j-1] )*( z_s*psi_pos_y(alpha_d,i,j-1,1,1,alphad_inlet) - (1 - z_s)*psi_neg_y(alpha_d,i,j,1,1) );    
        }

        // Right wall
        
        i=NX-1;
        
        for (j=1;j<NY-1;j++)
        {
            uE = u_d[i+1][j]*DY;
            uW = u_d[i][j]*DY;
            vN = v_d[i][j+1]*DX;
            vS = v_d[i][j]*DX;

            z_e = uE > 0;
            z_w = uW > 0;
            z_n = vN > 0;
            z_s = vS > 0;

            m = j*NX + i;

            alphad_e[i][j] = 0;
            alphad_w[i][j] = -( max(uW,0) );
            alphad_n[i][j] = -( max(-vN,0) );
            alphad_s[i][j] = -( max(vS,0) );
            alphad_p[i][j] = -(alphad_e[i][j] + alphad_w[i][j] + alphad_n[i][j] + alphad_s[i][j]) + DX2*DT + uE - uW + vN - vS;
            Qa[m] = (1-ALPHA_COR)*(alphad_p[i][j]/ALPHA_COR)*alpha_d[i][j] + alphan_d[i][j]*DX2*DT + 0.5*uW*( alpha_d[i][j] - alpha_d[i-1][j] )*( z_w*psi_pos_x(alpha_d,i-1,j,1,1) - (1 - z_w)*psi_neg_x(alpha_d,i,j,1,1) ) + 0.5*vN*( alpha_d[i][j+1] - alpha_d[i][j] )*( (1 - z_n)*psi_neg_y(alpha_d,i,j+1,0,0) - z_n*psi_pos_y(alpha_d,i,j,0,0,alphad_inlet) ) + 0.5*vS*( alpha_d[i][j] - alpha_d[i][j-1] )*( z_s*psi_pos_y(alpha_d,i,j-1,0,0,alphad_inlet) - (1 - z_s)*psi_neg_y(alpha_d,i,j,0,0) );            
        }
        
        // Left wall
        
        i=0;
        
        for (j=1;j<NY-1;j++)
        {
            uE = u_d[i+1][j]*DY;
            uW = u_d[i][j]*DY;
            vN = v_d[i][j+1]*DX;
            vS = v_d[i][j]*DX;

            z_e = uE > 0;
            z_w = uW > 0;
            z_n = vN > 0;
            z_s = vS > 0;

            m = j*NX + i;

            alphad_e[i][j] = -( max(-uE,0) );
            alphad_w[i][j] = 0;
            alphad_n[i][j] = -( max(-vN,0) );
            alphad_s[i][j] = -( max(vS,0) );
            alphad_p[i][j] = -(alphad_e[i][j] + alphad_w[i][j] + alphad_n[i][j] + alphad_s[i][j]) + DX2*DT + uE - uW + vN - vS;
            Qa[m] = (1-ALPHA_COR)*(alphad_p[i][j]/ALPHA_COR)*alpha_d[i][j] + alphan_d[i][j]*DX2*DT + 0.5*uE*( alpha_d[i+1][j] - alpha_d[i][j] )*( (1 - z_e)*psi_neg_x(alpha_d,i+1,j,1,0) - z_e*psi_pos_x(alpha_d,i,j,1,0) ) + 0.5*vN*( alpha_d[i][j+1] - alpha_d[i][j] )*( (1 - z_n)*psi_neg_y(alpha_d,i,j+1,0,0) - z_n*psi_pos_y(alpha_d,i,j,0,0,alphad_inlet) ) + 0.5*vS*( alpha_d[i][j] - alpha_d[i][j-1] )*( z_s*psi_pos_y(alpha_d,i,j-1,0,0,alphad_inlet) - (1 - z_s)*psi_neg_y(alpha_d,i,j,0,0) );            
        }


        // Top Left
        
        j=NY-1;
        
        uE = u_d[i+1][j]*DY;
        uW = u_d[i][j]*DY;
        vN = v_d[i][j+1]*DX;
        vS = v_d[i][j]*DX;

        z_e = uE > 0;
        z_w = uW > 0;
        z_n = vN > 0;
        z_s = vS > 0;

        m = j*NX + i;

        alphad_e[i][j] = -( max(-uE,0) );
        alphad_w[i][j] = 0;
        alphad_n[i][j] = 0;
        alphad_s[i][j] = -( max(vS,0) );
        alphad_p[i][j] = -(alphad_e[i][j] + alphad_w[i][j] + alphad_n[i][j] + alphad_s[i][j]) + DX2*DT + uE - uW + vN - vS + DX*v_d[i][j];
        Qa[m] = (1-ALPHA_COR)*(alphad_p[i][j]/ALPHA_COR)*alpha_d[i][j] + alphan_d[i][j]*DX2*DT + 0.5*uE*( alpha_d[i+1][j] - alpha_d[i][j] )*( (1 - z_e)*psi_neg_x(alpha_d,i+1,j,1,0) - z_e*psi_pos_x(alpha_d,i,j,1,0) ) + 0.5*vS*( alpha_d[i][j] - alpha_d[i][j-1] )*( z_s*psi_pos_y(alpha_d,i,j-1,1,1,alphad_inlet) - (1 - z_s)*psi_neg_y(alpha_d,i,j,1,1) );    
        
        
        // Bottom Left
        
        j=0;

        uE = u_d[i+1][j]*DY;
        uW = u_d[i][j]*DY;
        vN = v_d[i][j+1]*DX;
        vS = v_d[i][j]*DX;

        z_e = uE > 0;
        z_w = uW > 0;
        z_n = vN > 0;
        z_s = vS > 0;

        m = j*NX + i;

        alphad_e[i][j] = -( max(-uE,0) );
        alphad_w[i][j] = 0;
        alphad_n[i][j] = -( max(-vN,0) );
        alphad_s[i][j] = 0;
        alphad_p[i][j] = -(alphad_e[i][j] + alphad_w[i][j] + alphad_n[i][j] + alphad_s[i][j]) + DX2*DT + uE - uW + vN - vS;
        Qa[m] = alphad_inlet[i]*vS + (1-ALPHA_COR)*(alphad_p[i][j]/ALPHA_COR)*alpha_d[i][j] + alphan_d[i][j]*DX2*DT + 0.5*uE*( alpha_d[i+1][j] - alpha_d[i][j] )*( (1 - z_e)*psi_neg_x(alpha_d,i+1,j,1,0) - z_e*psi_pos_x(alpha_d,i,j,1,0) ) + 0.5*vN*( alpha_d[i][j+1] - alpha_d[i][j] )*( (1 - z_n)*psi_neg_y(alpha_d,i,j+1,1,0) - z_n*psi_pos_y(alpha_d,i,j,1,0,alphad_inlet) );    


        // Top Right
        
        i=NX-1;
        j=NY-1;
        
        m = j*NX + i;

        uE = u_d[i+1][j]*DY;
        uW = u_d[i][j]*DY;
        vN = v_d[i][j+1]*DX;
        vS = v_d[i][j]*DX;

        z_e = uE > 0;
        z_w = uW > 0;
        z_n = vN > 0;
        z_s = vS > 0;

        alphad_e[i][j] = 0;
        alphad_w[i][j] = -( max(uW,0) );
        alphad_n[i][j] = 0;
        alphad_s[i][j] = -( max(vS,0) );
        alphad_p[i][j] = -(alphad_e[i][j] + alphad_w[i][j] + alphad_n[i][j] + alphad_s[i][j]) + DX2*DT + uE - uW + vN - vS + DX*v_d[i][j];
        Qa[m] = (1-ALPHA_COR)*(alphad_p[i][j]/ALPHA_COR)*alpha_d[i][j] + alphan_d[i][j]*DX2*DT + 0.5*uW*( alpha_d[i][j] - alpha_d[i-1][j] )*( z_w*psi_pos_x(alpha_d,i-1,j,1,1) - (1 - z_w)*psi_neg_x(alpha_d,i,j,1,1) ) + 0.5*vS*( alpha_d[i][j] - alpha_d[i][j-1] )*( z_s*psi_pos_y(alpha_d,i,j-1,1,1,alphad_inlet) - (1 - z_s)*psi_neg_y(alpha_d,i,j,1,1) );    
        
        // Bottom Right
        
        j=0;
        
        uE = u_d[i+1][j]*DY;
        uW = u_d[i][j]*DY;
        vN = v_d[i][j+1]*DX;
        vS = v_d[i][j]*DX;

        z_e = uE > 0;
        z_w = uW > 0;
        z_n = vN > 0;
        z_s = vS > 0;

        m = j*NX + i;

        alphad_e[i][j] = 0;
        alphad_w[i][j] = -( max(uW,0) );
        alphad_n[i][j] = -( max(-vN,0) );
        alphad_s[i][j] = 0;
        alphad_p[i][j] = -(alphad_e[i][j] + alphad_w[i][j] + alphad_n[i][j] + alphad_s[i][j]) + DX2*DT + uE - uW + vN - vS;
        Qa[m] = alphad_inlet[i]*vS + (1-ALPHA_COR)*(alphad_p[i][j]/ALPHA_COR)*alpha_d[i][j] + alphan_d[i][j]*DX2*DT + 0.5*uW*( alpha_d[i][j] - alpha_d[i-1][j] )*( z_w*psi_pos_x(alpha_d,i-1,j,1,1) - (1 - z_w)*psi_neg_x(alpha_d,i,j,1,1) ) + 0.5*vN*( alpha_d[i][j+1] - alpha_d[i][j] )*( (1 - z_n)*psi_neg_y(alpha_d,i,j+1,1,0) - z_n*psi_pos_y(alpha_d,i,j,1,0,alphad_inlet) );    
        

        ////// End of alphad-coefficients

        ///// Solver for Phasic Continuity


        for (i=0;i<NX;i++)
        {
            for (j=0;j<NY;j++)
            {
                m = j*NX + i;

                Ea[m] = alphad_p[i][j]/ALPHA_COR;
                Fa[m] = alphad_e[i][j];
                Da[m] = alphad_w[i][j];
                Ba[m] = alphad_s[i][j];
                Ha[m] = alphad_n[i][j];
                a[m] = alpha_d[i][j];
            }
        }

        coefficient_matrix(Ea, Fa, Da, Ba, Ha, fa, da, ba, ea, ca, Kp, NX, ALPHA_P);

        //// Inner loop begins

        for (inner=1;inner<=ITERATION_ALPHA;inner++)
        {

            residual_alpha = SIP_Solver_Centre (Ea, Fa, Da, Ba, Ha, Qa, fa, da, ba, ea, ca, resa, vec_Ya, vec_deltaa, Kp, NX, a, alpha_d);

            nor = 0;
            
            for (i=0;i<NX;i++)
            {
                for(j=0;j<NY;j++)
                {
                    nor += fabs(alphad_p[i][j]*alpha_d[i][j]);
                }
            }
                
            residual_alpha = residual_alpha/nor;

            if (residual_alpha<TOL_ALPHA)
            {
                printf("residual_alpha at iteration %d, outer, %d, time, %d, = %12.5e\n",inner,outer,time,residual_alpha);
                printf("Solution converged!\n");
                // resalpha=fopen("res_alpha.txt","a");
                // fprintf(resalpha,"%s ,%d, time, %d, = % 12.5e\n",statement_alpha,outer,time,residual_alpha);
                // fclose(resalpha);
                break;
            }

            
            if (inner == ITERATION_ALPHA)
            {
                // resalpha=fopen("res_alpha.txt","a");
                // fprintf(resalpha,"%s , %d, outer, %d, time, %d, = % 12.5e\n",statement_alpha,inner,outer,time,residual_alpha);
                // fclose(resalpha);
                printf("residual_alpha at iteration %d at outer %d at time %d = % 12.5e\n",inner,outer,time,residual_alpha);
            }   
            
        }


        for (i=0;i<NX;i++)
        {
            free(alphad_p[i]);
            free(alphad_e[i]);
            free(alphad_w[i]);
            free(alphad_n[i]);
            free(alphad_s[i]);
        }
        free(alphad_p);
        free(alphad_e);
        free(alphad_w);
        free(alphad_n);
        free(alphad_s);

        for (i=0;i<NX;i++)
        {
            for (j=0;j<NY;j++)
            {
                alpha_c[i][j] = 1 - alpha_d[i][j];
                if (alpha_d[i][j] < 0) 
                {
                    printf("%d \t%d\n",i,j );
                    // printf("Negative Volume Fraction Detected\n");
                    // exit(1);
                }
            }
        }
        
        //// Turbulence Modeling
        
        /// k-epsilon Turbulence Model
        
        /// k-equation coefficients
        

        for (i=0;i<NX;i++)
        {
            for (j=0;j<NY;j++)
            {
                U_mc = (u_c[i][j] + u_c[i+1][j])/2;
                V_mc = (v_c[i][j] + v_c[i][j+1])/2;

                U_md = (u_d[i][j] + u_d[i+1][j])/2;
                V_md = (v_d[i][j] + v_d[i][j+1])/2;

                uslip_mag[i][j] = sqrt( pow((U_md - U_mc),2) + pow((V_md - V_mc),2) );

                DIV[i][j] = ( u_c[i+1][j] - u_c[i][j] )/DX + ( v_c[i][j+1] - v_c[i][j] )/DY;
            }
        }

        //// Source term calculation

        for (i=1;i<NX-1;i++)
        {
            for (j=1;j<NY-1;j++)
            {
                U_x = (u_c[i+1][j] - u_c[i][j])/DX;
                V_y = (v_c[i][j+1] - v_c[i][j])/DY;
                top = (u_c[i][j+1] + u_c[i][j] + u_c[i+1][j] + u_c[i+1][j+1])/4;
                bottom = (u_c[i][j-1] + u_c[i][j] + u_c[i+1][j-1] + u_c[i+1][j])/4;
                left = (v_c[i-1][j+1] + v_c[i][j+1] + v_c[i-1][j] + v_c[i][j])/4;
                right = (v_c[i][j+1] + v_c[i+1][j+1] + v_c[i][j] + v_c[i+1][j])/4;
                U_y = (top-bottom)/DY;
                V_x = (right-left)/DX;
                source[i][j] =  pow((U_y + V_x),2)/2 + pow(U_x,2) + pow(V_y,2);
                production[i][j] = 2*mut_C[i][j]*source[i][j];// - 0.67*DIV[i][j]*mut_C[i][j]*DIV[i][j];// - 0.67*DIV[i][j]*RHO_C*k[i][j];
            }
        }

        j=0;

        for (i=1;i<NX-1;i++)
        {
            U_x = (u_c[i+1][j] - u_c[i][j])/DX;
            V_y = (v_c[i][j+1] - v_c[i][j])/DY;
            top = (u_c[i][j+1] + u_c[i][j] + u_c[i+1][j] + u_c[i+1][j+1])/4;
            left = (v_c[i-1][j+1] + v_c[i][j+1] + v_c[i-1][j] + v_c[i][j])/4;
            right = (v_c[i][j+1] + v_c[i+1][j+1] + v_c[i][j] + v_c[i+1][j])/4;
            bottom = 0;//(u_c[i][j] + u_c[i+1][j])/4;
            U_y = (top-bottom)/DY;
            V_x = (right-left)/DX;
            source[i][j] =  pow((U_y + V_x),2)/2 + pow(U_x,2) + pow(V_y,2);
            production[i][j] = 2*mut_C[i][j]*source[i][j];// - 0.67*DIV[i][j]*mut_C[i][j]*DIV[i][j];// - 0.67*DIV[i][j]*RHO_C*k[i][j];
        }

        j=NY-1;

        for (i=1;i<NX-1;i++)
        {
            U_x = (u_c[i+1][j] - u_c[i][j])/DX;
            V_y = (v_c[i][j+1] - v_c[i][j])/DY;
            left = (v_c[i-1][j+1] + v_c[i][j+1] + v_c[i-1][j] + v_c[i][j])/4;
            right = (v_c[i][j+1] + v_c[i+1][j+1] + v_c[i][j] + v_c[i+1][j])/4;
            top = (u_c[i][j] + u_c[i+1][j])/2;
            bottom = (u_c[i][j-1] + u_c[i][j] + u_c[i+1][j-1] + u_c[i+1][j])/4;
            U_y = (top-bottom)/DY;
            V_x = (right-left)/DX;
            source[i][j] =  pow((U_y + V_x),2)/2 + pow(U_x,2) + pow(V_y,2);
            production[i][j] = 2*mut_C[i][j]*source[i][j];// - 0.67*DIV[i][j]*mut_C[i][j]*DIV[i][j];// - 0.67*DIV[i][j]*RHO_C*k[i][j];
        }

            
        double ** k_p = calloc(NX,sizeof(double *));
        double ** k_e = calloc(NX,sizeof(double *));
        double ** k_w = calloc(NX,sizeof(double *));
        double ** k_n = calloc(NX,sizeof(double *));
        double ** k_s = calloc(NX,sizeof(double *));
        for (i=0;i<NX;i++)
        {
            k_p[i] = calloc(NY,sizeof(double));
            k_e[i] = calloc(NY,sizeof(double));
            k_w[i] = calloc(NY,sizeof(double));
            k_n[i] = calloc(NY,sizeof(double));
            k_s[i] = calloc(NY,sizeof(double));
        }
        // Interior cells
        
        for (i=1;i<NX-1;i++)
        {
            for (j=1;j<NY-1;j++)
            {
                alpe_c = (alpha_c[i][j] + alpha_c[i+1][j])/2;
                alpw_c = (alpha_c[i][j] + alpha_c[i-1][j])/2;
                alpn_c = (alpha_c[i][j] + alpha_c[i][j+1])/2;
                alps_c = (alpha_c[i][j] + alpha_c[i][j-1])/2;

                mut_e = alpe_c*(mut_C[i][j] + mut_C[i+1][j])/2;
                mut_w = alpw_c*(mut_C[i][j] + mut_C[i-1][j])/2;
                mut_n = alpn_c*(mut_C[i][j] + mut_C[i][j+1])/2;
                mut_s = alps_c*(mut_C[i][j] + mut_C[i][j-1])/2;

                uE = RHO_C*alpe_c*u_c[i+1][j]*DY;
                uW = RHO_C*alpw_c*u_c[i][j]*DY;
                vN = RHO_C*alpn_c*v_c[i][j+1]*DX;
                vS = RHO_C*alps_c*v_c[i][j]*DX;

                z_e = uE > 0;
                z_w = uW > 0;
                z_n = vN > 0;
                z_s = vS > 0;
                
                k_e[i][j] = -(max(-uE,0) + mut_e*AR);
                k_w[i][j] = -(max(uW,0) + mut_w*AR);
                k_n[i][j] = -(max(-vN,0) + mut_n*RA);
                k_s[i][j] = -(max(vS,0) + mut_s*RA);
                k_p[i][j] = -(k_e[i][j] + k_w[i][j] + k_n[i][j] + k_s[i][j]) + alpha_c[i][j]*RHO_C*DX2*DT + uE - uW + vN - vS;
                m = j*NX + i;

                Qk[m] = alpha_c[i][j]*production[i][j]*DX2 - alpha_c[i][j]*RHO_C*epsilon[i][j]*DX2 + (1-K_COR)*(k_p[i][j]/K_COR)*k[i][j] + alphan_c[i][j]*RHO_C*kn[i][j]*DX2*DT + 0.5*uE*( k[i+1][j] - k[i][j] )*( (1 - z_e)*psi_neg_x(k,i+1,j,0,0) - z_e*psi_pos_x(k,i,j,0,0) ) + 0.5*uW*( k[i][j] - k[i-1][j] )*( z_w*psi_pos_x(k,i-1,j,0,0) - (1 - z_w)*psi_neg_x(k,i,j,0,0) ) + 0.5*vN*( k[i][j+1] - k[i][j] )*( (1 - z_n)*psi_neg_y(k,i,j+1,0,0) - z_n*psi_pos_y(k,i,j,0,0,k_inlet) ) + 0.5*vS*( k[i][j] - k[i][j-1] )*( z_s*psi_pos_y(k,i,j-1,0,0,k_inlet) - (1 - z_s)*psi_neg_y(k,i,j,0,0) );
        
            }
        }
        
        // Left Boundary
        
        i=0;
        
        for (j=1;j<NY-1;j++)
        {
            alpe_c = (alpha_c[i][j] + alpha_c[i+1][j])/2;
            alpw_c = alpha_c[i][j];
            alpn_c = (alpha_c[i][j] + alpha_c[i][j+1])/2;
            alps_c = (alpha_c[i][j] + alpha_c[i][j-1])/2;

            mut_e = alpe_c*(mut_C[i][j] + mut_C[i+1][j])/2;
            mut_n = alpn_c*(mut_C[i][j] + mut_C[i][j+1])/2;
            mut_s = alps_c*(mut_C[i][j] + mut_C[i][j-1])/2;
            x_plus = RHO_C*xp*pow(C_MU,0.25)*pow(k[i][j],0.5)/mu;
            v_plus = log(EE*x_plus)/KAPPA;
            vp = ( v_c[i][j] + v_c[i][j+1] )/2;
            tau_y = RHO_C*pow(C_MU,0.25)*pow(k[i][j],0.5)*vp/v_plus;
            // printf("%12.5e\n",x_plus );

            uE = RHO_C*alpe_c*u_c[i+1][j]*DY;
            uW = RHO_C*alpw_c*u_c[i][j]*DY;
            vN = RHO_C*alpn_c*v_c[i][j+1]*DX;
            vS = RHO_C*alps_c*v_c[i][j]*DX;

            z_e = uE > 0;
            z_w = uW > 0;
            z_n = vN > 0;
            z_s = vS > 0;

            k_e[i][j] = -(max(-uE,0) + mut_e*AR);
            k_w[i][j] = 0;
            k_n[i][j] = -(max(-vN,0) + mut_n*RA);
            k_s[i][j] = -(max(vS,0) + mut_s*RA);
            k_p[i][j] = -(k_e[i][j] + k_w[i][j] + k_n[i][j] + k_s[i][j]) + RHO_C*pow(C_MU,0.75)*pow(k[i][j],0.5)*DX2*v_plus/xp + alpha_c[i][j]*RHO_C*DX2*DT + uE - uW + vN - vS;
            m = j*NX + i;
            
            Qk[m] = tau_y*vp*DX2/xp + (1-K_COR)*(k_p[i][j]/K_COR)*k[i][j] + alphan_c[i][j]*RHO_C*kn[i][j]*DX2*DT + 0.5*uE*( k[i+1][j] - k[i][j] )*( (1 - z_e)*psi_neg_x(k,i+1,j,1,0) - z_e*psi_pos_x(k,i,j,1,0) ) + 0.5*vN*( k[i][j+1] - k[i][j] )*( (1 - z_n)*psi_neg_y(k,i,j+1,0,0) - z_n*psi_pos_y(k,i,j,0,0,k_inlet) ) + 0.5*vS*( k[i][j] - k[i][j-1] )*( z_s*psi_pos_y(k,i,j-1,0,0,k_inlet) - (1 - z_s)*psi_neg_y(k,i,j,0,0) );

        }
        
        // Outlet
        
        j=NY-1;
        
        for (i=1;i<NX-1;i++)
        {            
            alpe_c = (alpha_c[i][j] + alpha_c[i+1][j])/2;
            alpw_c = (alpha_c[i][j] + alpha_c[i-1][j])/2;
            alpn_c = alpha_c[i][j];
            alps_c = (alpha_c[i][j] + alpha_c[i][j-1])/2;

            mut_e = alpe_c*(mut_C[i][j] + mut_C[i+1][j])/2;
            mut_w = alpw_c*(mut_C[i][j] + mut_C[i-1][j])/2;
            mut_s = alps_c*(mut_C[i][j] + mut_C[i][j-1])/2;

            uE = RHO_C*alpe_c*u_c[i+1][j]*DY;
            uW = RHO_C*alpw_c*u_c[i][j]*DY;
            vN = RHO_C*alpn_c*v_c[i][j+1]*DX;
            vS = RHO_C*alps_c*v_c[i][j]*DX;

            z_e = uE > 0;
            z_w = uW > 0;
            z_n = vN > 0;
            z_s = vS > 0;

            k_e[i][j] = -(max(-uE,0) + mut_e*AR);
            k_w[i][j] = -(max(uW,0) + mut_w*AR);
            k_n[i][j] = 0;
            k_s[i][j] = -(max(vS,0) + mut_s*RA);
            k_p[i][j] = -(k_e[i][j] + k_w[i][j] + k_n[i][j] + k_s[i][j]) + alpha_c[i][j]*RHO_C*DX2*DT + uE - uW + vN - vS;
            m = j*NX + i;

            Qk[m] = alpha_c[i][j]*production[i][j]*DX2 - alpha_c[i][j]*RHO_C*epsilon[i][j]*DX2 + (1-K_COR)*(k_p[i][j]/K_COR)*k[i][j] + alphan_c[i][j]*RHO_C*kn[i][j]*DX2*DT + 0.5*uE*( k[i+1][j] - k[i][j] )*( (1 - z_e)*psi_neg_x(k,i+1,j,0,0) - z_e*psi_pos_x(k,i,j,0,0) ) + 0.5*uW*( k[i][j] - k[i-1][j] )*( z_w*psi_pos_x(k,i-1,j,0,0) - (1 - z_w)*psi_neg_x(k,i,j,0,0) ) + 0.5*vS*( k[i][j] - k[i][j-1] )*( z_s*psi_pos_y(k,i,j-1,1,1,k_inlet) - (1 - z_s)*psi_neg_y(k,i,j,1,1) );
        }
        
        // Right Boundary
        
        i=NX-1;
        
        for (j=1;j<NY-1;j++)
        {            
            alpe_c = alpha_c[i][j];
            alpw_c = (alpha_c[i][j] + alpha_c[i-1][j])/2;
            alpn_c = (alpha_c[i][j] + alpha_c[i][j+1])/2;
            alps_c = (alpha_c[i][j] + alpha_c[i][j-1])/2;

            mut_w = alpw_c*(mut_C[i][j] + mut_C[i-1][j])/2;
            mut_n = alpn_c*(mut_C[i][j] + mut_C[i][j+1])/2;
            mut_s = alps_c*(mut_C[i][j] + mut_C[i][j-1])/2;
            x_plus = RHO_C*xp*pow(C_MU,0.25)*pow(k[i][j],0.5)/mu;
            v_plus = log(EE*x_plus)/KAPPA;
            vp = ( v_c[i][j] + v_c[i][j+1] )/2;
            tau_y = RHO_C*pow(C_MU,0.25)*pow(k[i][j],0.5)*vp/v_plus;

            uE = RHO_C*alpe_c*u_c[i+1][j]*DY;
            uW = RHO_C*alpw_c*u_c[i][j]*DY;
            vN = RHO_C*alpn_c*v_c[i][j+1]*DX;
            vS = RHO_C*alps_c*v_c[i][j]*DX;

            z_e = uE > 0;
            z_w = uW > 0;
            z_n = vN > 0;
            z_s = vS > 0;

            k_e[i][j] = 0;
            k_w[i][j] = -(max(uW,0) + mut_w*AR);
            k_n[i][j] = -(max(-vN,0) + mut_n*RA);
            k_s[i][j] = -(max(vS,0) + mut_s*RA);
            k_p[i][j] = -(k_e[i][j] + k_w[i][j] + k_n[i][j] + k_s[i][j]) + RHO_C*pow(C_MU,0.75)*pow(k[i][j],0.5)*DX2*v_plus/xp + alpha_c[i][j]*RHO_C*DX2*DT + uE - uW + vN - vS;
            m = j*NX + i;
            
            Qk[m] = tau_y*vp*DX2/xp + (1-K_COR)*(k_p[i][j]/K_COR)*k[i][j] + alphan_c[i][j]*RHO_C*kn[i][j]*DX2*DT + 0.5*uW*( k[i][j] - k[i-1][j] )*( z_w*psi_pos_x(k,i-1,j,1,1) - (1 - z_w)*psi_neg_x(k,i,j,1,1) ) + 0.5*vN*( k[i][j+1] - k[i][j] )*( (1 - z_n)*psi_neg_y(k,i,j+1,0,0) - z_n*psi_pos_y(k,i,j,0,0,k_inlet) ) + 0.5*vS*( k[i][j] - k[i][j-1] )*( z_s*psi_pos_y(k,i,j-1,0,0,k_inlet) - (1 - z_s)*psi_neg_y(k,i,j,0,0) );

        }
        
        // Inlet
        
        j=0;
        
        for (i=1;i<NX-1;i++)
        {
            alpe_c = (alpha_c[i][j] + alpha_c[i+1][j])/2;
            alpw_c = (alpha_c[i][j] + alpha_c[i-1][j])/2;
            alpn_c = (alpha_c[i][j] + alpha_c[i][j+1])/2;
            alps_c = alphac_inlet[i];
            
            mut_e = alpe_c*(mut_C[i][j] + mut_C[i+1][j])/2;
            mut_w = alpw_c*(mut_C[i][j] + mut_C[i-1][j])/2;
            mut_n = alpn_c*(mut_C[i][j] + mut_C[i][j+1])/2;
            
            uE = RHO_C*alpe_c*u_c[i+1][j]*DY;
            uW = RHO_C*alpw_c*u_c[i][j]*DY;
            vN = RHO_C*alpn_c*v_c[i][j+1]*DX;
            vS = RHO_C*alps_c*v_c[i][j]*DX;

            z_e = uE > 0;
            z_w = uW > 0;
            z_n = vN > 0;
            z_s = vS > 0;

            k_e[i][j] = -(max(-uE,0) + mut_e*AR);
            k_w[i][j] = -(max(uW,0) + mut_w*AR);
            k_n[i][j] = -(max(-vN,0) + mut_n*RA + mut_c_inlet[i]*alps_c*RA/3);
            k_s[i][j] = 0;
            k_p[i][j] = max(-uE,0) + max(uW,0) + max(-vN,0) + mut_n*RA + mut_e*AR + mut_w*AR + 3*mut_c_inlet[i]*alps_c*RA + alpha_c[i][j]*RHO_C*DX2*DT + uE - uW + vN - vS;
            m = j*NX + i;

            Qk[m] = alps_c*8*k_inlet[i]*mut_c_inlet[i]*RA/3 + alpha_c[i][j]*production[i][j]*DX2 - alpha_c[i][j]*RHO_C*epsilon[i][j]*DX2 + (1-K_COR)*(k_p[i][j]/K_COR)*k[i][j] + alphan_c[i][j]*RHO_C*kn[i][j]*DX2*DT + 0.5*uE*( k[i+1][j] - k[i][j] )*( (1 - z_e)*psi_neg_x(k,i+1,j,0,0) - z_e*psi_pos_x(k,i,j,0,0) ) + 0.5*uW*( k[i][j] - k[i-1][j] )*( z_w*psi_pos_x(k,i-1,j,0,0) - (1 - z_w)*psi_neg_x(k,i,j,0,0) ) + 0.5*vN*( k[i][j+1] - k[i][j] )*( (1 - z_n)*psi_neg_y(k,i,j+1,1,0) - z_n*psi_pos_y(k,i,j,1,0,k_inlet) );//alpha_c[i][j]*8*k_inlet[i]*mut_c_inlet[i]*RA/3 +
        }


        // Top Left

        i = 0;
        j = NY-1;

        alpe_c = (alpha_c[i][j] + alpha_c[i+1][j])/2;
        alpw_c = alpha_c[i][j];
        alpn_c = alpha_c[i][j];
        alps_c = (alpha_c[i][j] + alpha_c[i][j-1])/2;

        mut_e = alpe_c*(mut_C[i][j] + mut_C[i+1][j])/2;
        mut_s = alps_c*(mut_C[i][j] + mut_C[i][j-1])/2;
        x_plus = RHO_C*xp*pow(C_MU,0.25)*pow(k[i][j],0.5)/mu;
        v_plus = log(EE*x_plus)/KAPPA;
        vp = ( v_c[i][j] + v_c[i][j+1] )/2;
        tau_y = RHO_C*pow(C_MU,0.25)*pow(k[i][j],0.5)*vp/v_plus;
        // printf("%12.5e\n", x_plus);

        uE = RHO_C*alpe_c*u_c[i+1][j]*DY;
        uW = RHO_C*alpw_c*u_c[i][j]*DY;
        vN = RHO_C*alpn_c*v_c[i][j+1]*DX;
        vS = RHO_C*alps_c*v_c[i][j]*DX;

        z_e = uE > 0;
        z_w = uW > 0;
        z_n = vN > 0;
        z_s = vS > 0;
        
        k_e[i][j] = -(max(-uE,0) + mut_e*AR);
        k_w[i][j] = 0;
        k_n[i][j] = 0;
        k_s[i][j] = -(max(vS,0) + mut_s*RA);
        k_p[i][j] = -(k_e[i][j] + k_w[i][j] + k_n[i][j] + k_s[i][j]) + RHO_C*pow(C_MU,0.75)*pow(k[i][j],0.5)*DX2*v_plus/xp + alpha_c[i][j]*RHO_C*DX2*DT + uE - uW + vN - vS;
        m = j*NX + i;
        
        Qk[m] = tau_y*vp*DX2/xp + (1-K_COR)*(k_p[i][j]/K_COR)*k[i][j] + alphan_c[i][j]*RHO_C*kn[i][j]*DX2*DT + 0.5*uE*( k[i+1][j] - k[i][j] )*( (1 - z_e)*psi_neg_x(k,i+1,j,1,0) - z_e*psi_pos_x(k,i,j,1,0) ) + 0.5*vS*( k[i][j] - k[i][j-1] )*( z_s*psi_pos_y(k,i,j-1,1,1,k_inlet) - (1 - z_s)*psi_neg_y(k,i,j,1,1) );

        
        
        // Bottom Left
        
        j=0;

        alpe_c = (alpha_c[i][j] + alpha_c[i+1][j])/2;
        alpw_c = alpha_c[i][j];
        alpn_c = (alpha_c[i][j] + alpha_c[i][j+1])/2;
        alps_c = alphac_inlet[i];

        mut_e = alpe_c*(mut_C[i][j] + mut_C[i+1][j])/2;
        mut_n = alpn_c*(mut_C[i][j] + mut_C[i][j+1])/2;
        x_plus = RHO_C*xp*pow(C_MU,0.25)*pow(k[i][j],0.5)/mu;
        v_plus = log(EE*x_plus)/KAPPA;
        vp = ( v_c[i][j] + v_c[i][j+1] )/2;
        tau_y = RHO_C*pow(C_MU,0.25)*pow(k[i][j],0.5)*vp/v_plus;

        uE = RHO_C*alpe_c*u_c[i+1][j]*DY;
        uW = RHO_C*alpw_c*u_c[i][j]*DY;
        vN = RHO_C*alpn_c*v_c[i][j+1]*DX;
        vS = RHO_C*alps_c*v_c[i][j]*DX;

        z_e = uE > 0;
        z_w = uW > 0;
        z_n = vN > 0;
        z_s = vS > 0;

        k_e[i][j] = -(max(-uE,0) + mut_e*AR);
        k_w[i][j] = 0;
        k_n[i][j] = -(max(-vN,0) + mut_n*RA);
        k_s[i][j] = 0;
        k_p[i][j] = -(k_e[i][j] + k_w[i][j] + k_n[i][j] + k_s[i][j]) + RHO_C*pow(C_MU,0.75)*pow(k[i][j],0.5)*DX2*v_plus/xp + alpha_c[i][j]*RHO_C*DX2*DT + uE - uW + vN - vS;
        m = j*NX + i;
        
        Qk[m] = tau_y*vp*DX2/xp + (1-K_COR)*(k_p[i][j]/K_COR)*k[i][j] + alphan_c[i][j]*RHO_C*kn[i][j]*DX2*DT + 0.5*uE*( k[i+1][j] - k[i][j] )*( (1 - z_e)*psi_neg_x(k,i+1,j,1,0) - z_e*psi_pos_x(k,i,j,1,0) ) + 0.5*vN*( k[i][j+1] - k[i][j] )*( (1 - z_n)*psi_neg_y(k,i,j+1,1,0) - z_n*psi_pos_y(k,i,j,1,0,k_inlet) );



        // Top Right
        
        i=NX-1;
        j=NY-1;

        alpe_c = alpha_c[i][j];
        alpw_c = (alpha_c[i][j] + alpha_c[i-1][j])/2;
        alpn_c = alpha_c[i][j];
        alps_c = (alpha_c[i][j] + alpha_c[i][j-1])/2;

        mut_w = alpw_c*(mut_C[i][j] + mut_C[i-1][j])/2;
        mut_s = alps_c*(mut_C[i][j] + mut_C[i][j-1])/2;
        x_plus = RHO_C*xp*pow(C_MU,0.25)*pow(k[i][j],0.5)/mu;
        v_plus = log(EE*x_plus)/KAPPA;
        vp = ( v_c[i][j] + v_c[i][j+1] )/2;
        tau_y = RHO_C*pow(C_MU,0.25)*pow(k[i][j],0.5)*vp/v_plus;
        
        uE = RHO_C*alpe_c*u_c[i+1][j]*DY;
        uW = RHO_C*alpw_c*u_c[i][j]*DY;
        vN = RHO_C*alpn_c*v_c[i][j+1]*DX;
        vS = RHO_C*alps_c*v_c[i][j]*DX;

        z_e = uE > 0;
        z_w = uW > 0;
        z_n = vN > 0;
        z_s = vS > 0;

        k_e[i][j] = 0;
        k_w[i][j] = -(max(uW,0) + mut_w*AR);
        k_n[i][j] = 0;
        k_s[i][j] = -(max(vS,0) + mut_s*RA);
        k_p[i][j] = -(k_e[i][j] + k_w[i][j] + k_n[i][j] + k_s[i][j]) + RHO_C*pow(C_MU,0.75)*pow(k[i][j],0.5)*DX2*v_plus/xp + alpha_c[i][j]*RHO_C*DX2*DT + uE - uW + vN - vS;
        m = j*NX + i;
        
        Qk[m] = tau_y*vp*DX2/xp + (1-K_COR)*(k_p[i][j]/K_COR)*k[i][j] + alphan_c[i][j]*RHO_C*kn[i][j]*DX2*DT + 0.5*uW*( k[i][j] - k[i-1][j] )*( z_w*psi_pos_x(k,i-1,j,1,1) - (1 - z_w)*psi_neg_x(k,i,j,1,1) ) + 0.5*vS*( k[i][j] - k[i][j-1] )*( z_s*psi_pos_y(k,i,j-1,1,1,k_inlet) - (1 - z_s)*psi_neg_y(k,i,j,1,1) );

        
        // Bottom Right
        
        j=0;

        alpe_c = alpha_c[i][j];
        alpw_c = (alpha_c[i][j] + alpha_c[i-1][j])/2;
        alpn_c = (alpha_c[i][j] + alpha_c[i][j+1])/2;
        alps_c = alphac_inlet[i];

        mut_w = alpw_c*(mut_C[i][j] + mut_C[i-1][j])/2;
        mut_n = alpn_c*(mut_C[i][j] + mut_C[i][j+1])/2;
        x_plus = RHO_C*xp*pow(C_MU,0.25)*pow(k[i][j],0.5)/mu;
        v_plus = log(EE*x_plus)/KAPPA;
        vp = ( v_c[i][j] + v_c[i][j+1] )/2;
        tau_y = RHO_C*pow(C_MU,0.25)*pow(k[i][j],0.5)*vp/v_plus;

        uE = RHO_C*alpe_c*u_c[i+1][j]*DY;
        uW = RHO_C*alpw_c*u_c[i][j]*DY;
        vN = RHO_C*alpn_c*v_c[i][j+1]*DX;
        vS = RHO_C*alps_c*v_c[i][j]*DX;

        z_e = uE > 0;
        z_w = uW > 0;
        z_n = vN > 0;
        z_s = vS > 0;

        k_e[i][j] = 0;
        k_w[i][j] = -(max(uW,0) + mut_w*AR);
        k_n[i][j] = -(max(-vN,0) + mut_n*RA);
        k_s[i][j] = 0;
        k_p[i][j] = -(k_e[i][j] + k_w[i][j] + k_n[i][j] + k_s[i][j]) + RHO_C*pow(C_MU,0.75)*pow(k[i][j],0.5)*DX2*v_plus/xp + alpha_c[i][j]*RHO_C*DX2*DT + uE - uW + vN - vS;
        m = j*NX + i;
        
        Qk[m] = tau_y*vp*DX2/xp + (1-K_COR)*(k_p[i][j]/K_COR)*k[i][j] + alphan_c[i][j]*RHO_C*kn[i][j]*DX2*DT + 0.5*uW*( k[i][j] - k[i-1][j] )*( z_w*psi_pos_x(k,i-1,j,1,1) - (1 - z_w)*psi_neg_x(k,i,j,1,1) ) + 0.5*vN*( k[i][j+1] - k[i][j] )*( (1 - z_n)*psi_neg_y(k,i,j+1,1,0) - z_n*psi_pos_y(k,i,j,1,0,k_inlet) );

            
        
        ////// End of K-coefficients
        
        //// epsilon-equation coefficients
        
        double ** eps_p = calloc(NX,sizeof(double *));
        double ** eps_e = calloc(NX,sizeof(double *));
        double ** eps_w = calloc(NX,sizeof(double *));
        double ** eps_n = calloc(NX,sizeof(double *));
        double ** eps_s = calloc(NX,sizeof(double *));
        for (i=0;i<NX;i++)
        {
            eps_p[i] = calloc(NY,sizeof(double));
            eps_e[i] = calloc(NY,sizeof(double));
            eps_w[i] = calloc(NY,sizeof(double));
            eps_n[i] = calloc(NY,sizeof(double));
            eps_s[i] = calloc(NY,sizeof(double));
        }
        
        /// Interior Cells

        for (i=1;i<NX-1;i++)
        {
            for (j=1;j<NY-1;j++)
            {
                alpe_c = (alpha_c[i][j] + alpha_c[i+1][j])/2;
                alpw_c = (alpha_c[i][j] + alpha_c[i-1][j])/2;
                alpn_c = (alpha_c[i][j] + alpha_c[i][j+1])/2;
                alps_c = (alpha_c[i][j] + alpha_c[i][j-1])/2;

                mut_e = alpe_c*(mut_C[i][j] + mut_C[i+1][j])/2;
                mut_w = alpw_c*(mut_C[i][j] + mut_C[i-1][j])/2;
                mut_n = alpn_c*(mut_C[i][j] + mut_C[i][j+1])/2;
                mut_s = alps_c*(mut_C[i][j] + mut_C[i][j-1])/2;

                uE = RHO_C*alpe_c*u_c[i+1][j]*DY;
                uW = RHO_C*alpw_c*u_c[i][j]*DY;
                vN = RHO_C*alpn_c*v_c[i][j+1]*DX;
                vS = RHO_C*alps_c*v_c[i][j]*DX;

                z_e = uE > 0;
                z_w = uW > 0;
                z_n = vN > 0;
                z_s = vS > 0;

                eps_e[i][j] = -(max(-uE,0) + mut_e*r);
                eps_w[i][j] = -(max(uW,0) + mut_w*r);
                eps_n[i][j] = -(max(-vN,0) + mut_n*R);
                eps_s[i][j] = -(max(vS,0) + mut_s*R);
                eps_p[i][j] = -(eps_e[i][j] + eps_w[i][j] + eps_n[i][j] + eps_s[i][j]) + (C2*epsilon[i][j]*DX2*RHO_C*alpha_c[i][j])/k[i][j] + alpha_c[i][j]*RHO_C*DX2*DT + uE - uW + vN - vS;
                m = j*NX + i;

                Qe[m] = ( C1*alpha_c[i][j]*production[i][j]*epsilon[i][j]*DX2 )/k[i][j] + (1-EPS_COR)*(eps_p[i][j]/EPS_COR)*epsilon[i][j] + alphan_c[i][j]*RHO_C*epsn[i][j]*DX2*DT + 0.5*uE*( epsilon[i+1][j] - epsilon[i][j] )*( (1 - z_e)*psi_neg_x(epsilon,i+1,j,0,0) - z_e*psi_pos_x(epsilon,i,j,0,0) ) + 0.5*uW*( epsilon[i][j] - epsilon[i-1][j] )*( z_w*psi_pos_x(epsilon,i-1,j,0,0) - (1 - z_w)*psi_neg_x(epsilon,i,j,0,0) ) + 0.5*vN*( epsilon[i][j+1] - epsilon[i][j] )*( (1 - z_n)*psi_neg_y(epsilon,i,j+1,0,0) - z_n*psi_pos_y(epsilon,i,j,0,0,eps_inlet) ) + 0.5*vS*( epsilon[i][j] - epsilon[i][j-1] )*( z_s*psi_pos_y(epsilon,i,j-1,0,0,eps_inlet) - (1 - z_s)*psi_neg_y(epsilon,i,j,0,0) );
            }
        }

        // Left Boundary

        i=0;

        for (j=1;j<NY-1;j++)
        {
            alpe_c = (alpha_c[i][j] + alpha_c[i+1][j])/2;
            alpw_c = alpha_c[i][j];
            alpn_c = (alpha_c[i][j] + alpha_c[i][j+1])/2;
            alps_c = (alpha_c[i][j] + alpha_c[i][j-1])/2;

            mut_e = alpe_c*(mut_C[i][j] + mut_C[i+1][j])/2;
            mut_n = alpn_c*(mut_C[i][j] + mut_C[i][j+1])/2;
            mut_s = alps_c*(mut_C[i][j] + mut_C[i][j-1])/2;
            
            uE = RHO_C*alpe_c*u_c[i+1][j]*DY;
            uW = RHO_C*alpw_c*u_c[i][j]*DY;
            vN = RHO_C*alpn_c*v_c[i][j+1]*DX;
            vS = RHO_C*alps_c*v_c[i][j]*DX;

            z_e = uE > 0;
            z_w = uW > 0;
            z_n = vN > 0;
            z_s = vS > 0;

            eps_e[i][j] = -(max(-uE,0) + mut_e*r);
            eps_w[i][j] = 0;
            eps_n[i][j] = -(max(-vN,0) + mut_n*R);
            eps_s[i][j] = -(max(vS,0) + mut_s*R);
            eps_p[i][j] = -(eps_e[i][j] + eps_w[i][j] + eps_n[i][j] + eps_s[i][j]) + pow(10,30) + alpha_c[i][j]*RHO_C*DX2*DT + uE - uW + vN - vS;
            m = j*NX + i;
            
            Qe[m] = ( pow(C_MU,0.75)*pow(k[i][j],1.5)/kappa )*pow(10,30) + (1-EPS_COR)*(eps_p[i][j]/EPS_COR)*epsilon[i][j] + alphan_c[i][j]*RHO_C*epsn[i][j]*DX2*DT + 0.5*uE*( epsilon[i+1][j] - epsilon[i][j] )*( (1 - z_e)*psi_neg_x(epsilon,i+1,j,1,0) - z_e*psi_pos_x(epsilon,i,j,1,0) ) + 0.5*vN*( epsilon[i][j+1] - epsilon[i][j] )*( (1 - z_n)*psi_neg_y(epsilon,i,j+1,0,0) - z_n*psi_pos_y(epsilon,i,j,0,0,eps_inlet) ) + 0.5*vS*( epsilon[i][j] - epsilon[i][j-1] )*( z_s*psi_pos_y(epsilon,i,j-1,0,0,eps_inlet) - (1 - z_s)*psi_neg_y(epsilon,i,j,0,0) );

        }

        // Bottom-Left

        j=0;
                                    
            alpe_c = (alpha_c[i][j] + alpha_c[i+1][j])/2;
            alpw_c = alpha_c[i][j];
            alpn_c = (alpha_c[i][j] + alpha_c[i][j+1])/2;
            alps_c = alphac_inlet[i];
                        
            mut_e = alpe_c*(mut_C[i][j] + mut_C[i+1][j])/2;
            mut_n = alpn_c*(mut_C[i][j] + mut_C[i][j+1])/2;
            
            uE = RHO_C*alpe_c*u_c[i+1][j]*DY;
            uW = RHO_C*alpw_c*u_c[i][j]*DY;
            vN = RHO_C*alpn_c*v_c[i][j+1]*DX;
            vS = RHO_C*alps_c*v_c[i][j]*DX;

            z_e = uE > 0;
            z_w = uW > 0;
            z_n = vN > 0;
            z_s = vS > 0;
            
            eps_e[i][j] = -(max(-uE,0) + mut_e*r);
            eps_w[i][j] = 0;
            eps_n[i][j] = -(max(-vN,0) + mut_n*R);
            eps_s[i][j] = 0;
            eps_p[i][j] = -(eps_e[i][j] + eps_w[i][j] + eps_n[i][j] + eps_s[i][j]) + pow(10,30) + alpha_c[i][j]*RHO_C*DX2*DT + uE - uW + vN - vS;
            m = j*NX + i;
            
            Qe[m] = ( pow(C_MU,0.75)*pow(k[i][j],1.5)/kappa )*pow(10,30) + (1-EPS_COR)*(eps_p[i][j]/EPS_COR)*epsilon[i][j] + alphan_c[i][j]*RHO_C*epsn[i][j]*DX2*DT + 0.5*uE*( epsilon[i+1][j] - epsilon[i][j] )*( (1 - z_e)*psi_neg_x(epsilon,i+1,j,1,0) - z_e*psi_pos_x(epsilon,i,j,1,0) ) + 0.5*vN*( epsilon[i][j+1] - epsilon[i][j] )*( (1 - z_n)*psi_neg_y(epsilon,i,j+1,1,0) - z_n*psi_pos_y(epsilon,i,j,1,0,eps_inlet) );


        // Top-Left    

        j=NY-1;

            alpe_c = (alpha_c[i][j] + alpha_c[i+1][j])/2;
            alpw_c = alpha_c[i][j];
            alpn_c = alpha_c[i][j];
            alps_c = (alpha_c[i][j] + alpha_c[i][j-1])/2;

            mut_e = alpe_c*(mut_C[i][j] + mut_C[i+1][j])/2;
            mut_s = alps_c*(mut_C[i][j] + mut_C[i][j-1])/2;
            
            uE = RHO_C*alpe_c*u_c[i+1][j]*DY;
            uW = RHO_C*alpw_c*u_c[i][j]*DY;
            vN = RHO_C*alpn_c*v_c[i][j+1]*DX;
            vS = RHO_C*alps_c*v_c[i][j]*DX;

            z_e = uE > 0;
            z_w = uW > 0;
            z_n = vN > 0;
            z_s = vS > 0;
            
            eps_e[i][j] = -(max(-uE,0) + mut_e*r);
            eps_w[i][j] = 0;
            eps_n[i][j] = 0;
            eps_s[i][j] = -(max(vS,0) + mut_s*R);
            eps_p[i][j] = -(eps_e[i][j] + eps_w[i][j] + eps_n[i][j] + eps_s[i][j]) + pow(10,30) + alpha_c[i][j]*RHO_C*DX2*DT + uE - uW + vN - vS;
            m = j*NX + i;
            
            Qe[m] = ( pow(C_MU,0.75)*pow(k[i][j],1.5)/kappa )*pow(10,30) + (1-EPS_COR)*(eps_p[i][j]/EPS_COR)*epsilon[i][j] + alphan_c[i][j]*RHO_C*epsn[i][j]*DX2*DT + 0.5*uE*( epsilon[i+1][j] - epsilon[i][j] )*( (1 - z_e)*psi_neg_x(epsilon,i+1,j,1,0) - z_e*psi_pos_x(epsilon,i,j,1,0) ) + 0.5*vS*( epsilon[i][j] - epsilon[i][j-1] )*( z_s*psi_pos_y(epsilon,i,j-1,1,1,eps_inlet) - (1 - z_s)*psi_neg_y(epsilon,i,j,1,1) );


        // Inlet

        j=0;

            for (i=1;i<NX-1;i++)
            {
                alpe_c = (alpha_c[i][j] + alpha_c[i+1][j])/2;
                alpw_c = (alpha_c[i][j] + alpha_c[i-1][j])/2;
                alpn_c = (alpha_c[i][j] + alpha_c[i][j+1])/2;
                alps_c = alphac_inlet[i];

                mut_e = alpe_c*(mut_C[i][j] + mut_C[i+1][j])/2;
                mut_w = alpw_c*(mut_C[i][j] + mut_C[i-1][j])/2;
                mut_n = alpn_c*(mut_C[i][j] + mut_C[i][j+1])/2;
                
                uE = RHO_C*alpe_c*u_c[i+1][j]*DY;
                uW = RHO_C*alpw_c*u_c[i][j]*DY;
                vN = RHO_C*alpn_c*v_c[i][j+1]*DX;
                vS = RHO_C*alps_c*v_c[i][j]*DX;

                z_e = uE > 0;
                z_w = uW > 0;
                z_n = vN > 0;
                z_s = vS > 0;
                
                eps_p[i][j] = max(-uE,0) + max(uW,0) + max(-vN,0) + mut_e*r + mut_w*r + mut_n*R + 3*mut_c_inlet[i]*alps_c*R + (C2*epsilon[i][j]*DX2*RHO_C*alpha_c[i][j])/k[i][j] + alpha_c[i][j]*RHO_C*DX2*DT + uE - uW + vN - vS;
                eps_e[i][j] = -(max(-uE,0) + mut_e*r);
                eps_w[i][j] = -(max(uW,0) + mut_w*r);
                eps_n[i][j] = -(max(-vN,0) + mut_n*R + mut_c_inlet[i]*alps_c*R/3);
                eps_s[i][j] = 0;
                m = j*NX + i;

                Qe[m] = alps_c*8*eps_inlet[i]*mut_c_inlet[i]*R/3 + ( C1*alpha_c[i][j]*production[i][j]*epsilon[i][j]*DX2 )/k[i][j] + (1-EPS_COR)*(eps_p[i][j]/EPS_COR)*epsilon[i][j] + alphan_c[i][j]*RHO_C*epsn[i][j]*DX2*DT + 0.5*uE*( epsilon[i+1][j] - epsilon[i][j] )*( (1 - z_e)*psi_neg_x(epsilon,i+1,j,0,0) - z_e*psi_pos_x(epsilon,i,j,0,0) ) + 0.5*uW*( epsilon[i][j] - epsilon[i-1][j] )*( z_w*psi_pos_x(epsilon,i-1,j,0,0) - (1 - z_w)*psi_neg_x(epsilon,i,j,0,0) ) + 0.5*vN*( epsilon[i][j+1] - epsilon[i][j] )*( (1 - z_n)*psi_neg_y(epsilon,i,j+1,1,0) - z_n*psi_pos_y(epsilon,i,j,1,0,eps_inlet) );//8*mut_c_inlet[i]*eps_inlet[i]*R*alpha_c[i][j]/3 +

            }

            // Outlet

            j=NY-1;

            for (i=1;i<NX-1;i++)
            {
                alpe_c = (alpha_c[i][j] + alpha_c[i+1][j])/2;
                alpw_c = (alpha_c[i][j] + alpha_c[i-1][j])/2;
                alpn_c = alpha_c[i][j];
                alps_c = (alpha_c[i][j] + alpha_c[i][j-1])/2;

                mut_e = alpe_c*(mut_C[i][j] + mut_C[i+1][j])/2;
                mut_w = alpw_c*(mut_C[i][j] + mut_C[i-1][j])/2;
                mut_s = alps_c*(mut_C[i][j] + mut_C[i][j-1])/2;
                
                uE = RHO_C*alpe_c*u_c[i+1][j]*DY;
                uW = RHO_C*alpw_c*u_c[i][j]*DY;
                vN = RHO_C*alpn_c*v_c[i][j+1]*DX;
                vS = RHO_C*alps_c*v_c[i][j]*DX;

                z_e = uE > 0;
                z_w = uW > 0;
                z_n = vN > 0;
                z_s = vS > 0;
                
                eps_e[i][j] = -(max(-uE,0) + mut_e*r);
                eps_w[i][j] = -(max(uW,0) + mut_w*r);
                eps_n[i][j] = 0;
                eps_s[i][j] = -(max(vS,0) + mut_s*R);
                eps_p[i][j] = -(eps_e[i][j] + eps_w[i][j] + eps_n[i][j] + eps_s[i][j]) + (C2*epsilon[i][j]*DX2*RHO_C*alpha_c[i][j])/k[i][j] + alpha_c[i][j]*RHO_C*DX2*DT + uE - uW + vN - vS;
                m = j*NX + i;

                Qe[m] = ( C1*alpha_c[i][j]*production[i][j]*epsilon[i][j]*DX2 )/k[i][j] + (1-EPS_COR)*(eps_p[i][j]/EPS_COR)*epsilon[i][j] + alphan_c[i][j]*RHO_C*epsn[i][j]*DX2*DT + 0.5*uE*( epsilon[i+1][j] - epsilon[i][j] )*( (1 - z_e)*psi_neg_x(epsilon,i+1,j,0,0) - z_e*psi_pos_x(epsilon,i,j,0,0) ) + 0.5*uW*( epsilon[i][j] - epsilon[i-1][j] )*( z_w*psi_pos_x(epsilon,i-1,j,0,0) - (1 - z_w)*psi_neg_x(epsilon,i,j,0,0) ) + 0.5*vS*( epsilon[i][j] - epsilon[i][j-1] )*( z_s*psi_pos_y(epsilon,i,j-1,1,1,eps_inlet) - (1 - z_s)*psi_neg_y(epsilon,i,j,1,1) );
            }
        
        // Right Boundary
        
        i=NX-1;
        
        for (j=1;j<NY-1;j++)
        {
            alpe_c = alpha_c[i][j];
            alpw_c = (alpha_c[i][j] + alpha_c[i-1][j])/2;
            alpn_c = (alpha_c[i][j] + alpha_c[i][j+1])/2;
            alps_c = (alpha_c[i][j] + alpha_c[i][j-1])/2;

            mut_w = alpw_c*(mut_C[i][j] + mut_C[i-1][j])/2;
            mut_n = alpn_c*(mut_C[i][j] + mut_C[i][j+1])/2;
            mut_s = alps_c*(mut_C[i][j] + mut_C[i][j-1])/2;
            
            uE = RHO_C*alpe_c*u_c[i+1][j]*DY;
            uW = RHO_C*alpw_c*u_c[i][j]*DY;
            vN = RHO_C*alpn_c*v_c[i][j+1]*DX;
            vS = RHO_C*alps_c*v_c[i][j]*DX;

            z_e = uE > 0;
            z_w = uW > 0;
            z_n = vN > 0;
            z_s = vS > 0;
            
            eps_e[i][j] = 0;
            eps_w[i][j] = -(max(uW,0) + mut_w*r);
            eps_n[i][j] = -(max(-vN,0) + mut_n*R);
            eps_s[i][j] = -(max(vS,0) + mut_s*R);
            eps_p[i][j] = -(eps_e[i][j] + eps_w[i][j] + eps_n[i][j] + eps_s[i][j]) + pow(10,30) + alpha_c[i][j]*RHO_C*DX2*DT + uE - uW + vN - vS;
            m = j*NX + i;
            
            Qe[m] = ( pow(C_MU,0.75)*pow(k[i][j],1.5)/kappa )*pow(10,30) + (1-EPS_COR)*(eps_p[i][j]/EPS_COR)*epsilon[i][j] + alphan_c[i][j]*RHO_C*epsn[i][j]*DX2*DT + 0.5*uW*( epsilon[i][j] - epsilon[i-1][j] )*( z_w*psi_pos_x(epsilon,i-1,j,1,1) - (1 - z_w)*psi_neg_x(epsilon,i,j,1,1) ) + 0.5*vN*( epsilon[i][j+1] - epsilon[i][j] )*( (1 - z_n)*psi_neg_y(epsilon,i,j+1,0,0) - z_n*psi_pos_y(epsilon,i,j,0,0,eps_inlet) ) + 0.5*vS*( epsilon[i][j] - epsilon[i][j-1] )*( z_s*psi_pos_y(epsilon,i,j-1,0,0,eps_inlet) - (1 - z_s)*psi_neg_y(epsilon,i,j,0,0) );

        }
        
        // Top-Right
        
        i=NX-1;
        j=NY-1;

        alpe_c = alpha_c[i][j];
        alpw_c = (alpha_c[i][j] + alpha_c[i-1][j])/2;
        alpn_c = alpha_c[i][j];
        alps_c = (alpha_c[i][j] + alpha_c[i][j-1])/2;        

        mut_w = alpw_c*(mut_C[i][j] + mut_C[i-1][j])/2;
        mut_s = alps_c*(mut_C[i][j] + mut_C[i][j-1])/2; 
        
        uE = RHO_C*alpe_c*u_c[i+1][j]*DY;
        uW = RHO_C*alpw_c*u_c[i][j]*DY;
        vN = RHO_C*alpn_c*v_c[i][j+1]*DX;
        vS = RHO_C*alps_c*v_c[i][j]*DX;

        z_e = uE > 0;
        z_w = uW > 0;
        z_n = vN > 0;
        z_s = vS > 0;

        eps_e[i][j] = 0;
        eps_w[i][j] = -(max(uW,0) + mut_w*r);
        eps_n[i][j] = 0;
        eps_s[i][j] = -(max(vS,0) + mut_s*R);
        eps_p[i][j] = -(eps_e[i][j] + eps_w[i][j] + eps_n[i][j] + eps_s[i][j]) + pow(10,30) + alpha_c[i][j]*RHO_C*DX2*DT + uE - uW + vN - vS;
        m = j*NX + i;
        
        Qe[m] = ( pow(C_MU,0.75)*pow(k[i][j],1.5)/kappa )*pow(10,30) + (1-EPS_COR)*(eps_p[i][j]/EPS_COR)*epsilon[i][j] + alphan_c[i][j]*RHO_C*epsn[i][j]*DX2*DT + 0.5*uW*( epsilon[i][j] - epsilon[i-1][j] )*( z_w*psi_pos_x(epsilon,i-1,j,1,1) - (1 - z_w)*psi_neg_x(epsilon,i,j,1,1) ) + 0.5*vS*( epsilon[i][j] - epsilon[i][j-1] )*( z_s*psi_pos_y(epsilon,i,j-1,1,1,eps_inlet) - (1 - z_s)*psi_neg_y(epsilon,i,j,1,1) );

        
        // Bottom-Right
        
        j=0;

        alpe_c = alpha_c[i][j];
        alpw_c = (alpha_c[i][j] + alpha_c[i-1][j])/2;
        alpn_c = (alpha_c[i][j] + alpha_c[i][j+1])/2;
        alps_c = alphac_inlet[i];

        mut_w = alpw_c*(mut_C[i][j] + mut_C[i-1][j])/2;
        mut_n = alpn_c*(mut_C[i][j] + mut_C[i][j+1])/2;

        uE = RHO_C*alpe_c*u_c[i+1][j]*DY;
        uW = RHO_C*alpw_c*u_c[i][j]*DY;
        vN = RHO_C*alpn_c*v_c[i][j+1]*DX;
        vS = RHO_C*alps_c*v_c[i][j]*DX;

        z_e = uE > 0;
        z_w = uW > 0;
        z_n = vN > 0;
        z_s = vS > 0;

        eps_e[i][j] = 0;
        eps_w[i][j] = -(max(uW,0) + mut_w*r);
        eps_n[i][j] = -(max(-vN,0) + mut_n*R);
        eps_s[i][j] = 0;
        eps_p[i][j] = -(eps_e[i][j] + eps_w[i][j] + eps_n[i][j] + eps_s[i][j]) + pow(10,30) + alpha_c[i][j]*RHO_C*DX2*DT + uE - uW + vN - vS;
        m = j*NX + i;
        
        Qe[m] = ( pow(C_MU,0.75)*pow(k[i][j],1.5)/kappa )*pow(10,30) + (1-EPS_COR)*(eps_p[i][j]/EPS_COR)*epsilon[i][j] + alphan_c[i][j]*RHO_C*epsn[i][j]*DX2*DT + 0.5*uW*( epsilon[i][j] - epsilon[i-1][j] )*( z_w*psi_pos_x(epsilon,i-1,j,1,1) - (1 - z_w)*psi_neg_x(epsilon,i,j,1,1) ) + 0.5*vN*( epsilon[i][j+1] - epsilon[i][j] )*( (1 - z_n)*psi_neg_y(epsilon,i,j+1,1,0) - z_n*psi_pos_y(epsilon,i,j,1,0,eps_inlet) );

        ////// End of epsilon coefficients


        ///// Stone's ILU Method for k equation

        for (i=0;i<NX;i++)
        {
            for (j=0;j<NY;j++)
            {
                m = j*NX + i;

                Ek[m] = k_p[i][j]/K_COR;
                Fk[m] = k_e[i][j];
                Dk[m] = k_w[i][j];
                Bk[m] = k_s[i][j];
                Hk[m] = k_n[i][j];
                K[m] = k[i][j];
            }
        }

        coefficient_matrix(Ek, Fk, Dk, Bk, Hk, fk, dk, bk, ek, ck, Kp, NX, ALPHA_P);

        //// Inner loop begins

        for (inner=1;inner<=ITERATION_KEPS;inner++)
        {
            R2k = SIP_Solver_Centre (Ek, Fk, Dk, Bk, Hk, Qk, fk, dk, bk, ek, ck, resk, vec_Yk, vec_deltak, Kp, NX, K, k);

            nor = 0;
            
            for (i=0;i<NX;i++)
            {
                for(j=0;j<NY;j++)
                {
                    nor += fabs(k_p[i][j]*k[i][j]);
                }
            }
            
            residual_k = R2k/nor;

            if (residual_k<TOL_KEPS)
            {
                printf("residual_k at iteration %d, outer, %d, time, %d, = %12.5e\n",inner,outer,time,residual_k);
                printf("Solution converged!\n");
                // resK=fopen("res_k.txt","a");
                // fprintf(resK,"%s ,%d, time, %d, = % 12.5e\n",statement_k,outer,time,residual_k);
                // fclose(resK);
                break;
            }

        
            if (inner == ITERATION_KEPS)
            {
                // resK=fopen("res_k.txt","a");
                // fprintf(resK,"%s , %d, outer, %d, time, %d, = % 12.5e\n",statement_k,inner,outer,time,residual_k);
                // fclose(resK);
                printf("residual_k at iteration %d at outer %d at time %d = % 12.5e\n",inner,outer,time,residual_k);
            }   
            
        }


        //// End of k-equation solver
        
        for (i=0;i<NX;i++)
        {
            free(k_w[i]);
            free(k_s[i]);
            free(k_n[i]);
            free(k_e[i]);
            free(k_p[i]);
        }
        free(k_w);
        free(k_s);
        free(k_n);
        free(k_e);
        free(k_p);

        ///// Stone's ILU Method for Epsilon equation solver

        for (i=0;i<NX;i++)
        {
            for (j=0;j<NY;j++)
            {
                m = j*NX + i;

                Ee[m] = eps_p[i][j]/EPS_COR;
                Fe[m] = eps_e[i][j];
                De[m] = eps_w[i][j];
                Be[m] = eps_s[i][j];
                He[m] = eps_n[i][j];
                Eps[m] = epsilon[i][j];
            }
        }

        coefficient_matrix(Ee, Fe, De, Be, He, fe, de, be, ee, ce, Kp, NX, ALPHA_P);

        //// Inner loop begins

        for (inner=1;inner<=ITERATION_KEPS;inner++)
        {

            R2eps = SIP_Solver_Centre (Ee, Fe, De, Be, He, Qe, fe, de, be, ee, ce, reseps, vec_Yeps, vec_deltaeps, Kp, NX, Eps, epsilon);

            nor = 0;
            
            for (i=0;i<NX;i++)
            {
                for(j=0;j<NY;j++)
                {
                    nor += fabs(eps_p[i][j]*epsilon[i][j]);
                }
            }
            
            residual_eps = R2eps/nor;

            if (residual_eps<TOL_KEPS)
            {
                printf("residual_eps at iteration %d, outer, %d, time, %d, = %12.5e\n",inner,outer,time,residual_eps);
                printf("Solution converged!\n");
                // resEPS=fopen("res_eps.txt","a");
                // fprintf(resEPS,"%s ,%d,time, %d = % 12.5e\n",statement_eps,outer,time,residual_eps);
                // fclose(resEPS);
                break;
            }

        
            if (inner == ITERATION_KEPS)
            {
                // resEPS=fopen("res_eps.txt","a");
                // fprintf(resEPS,"%s , %d, outer, %d,time, %d, = % 12.5e\n",statement_eps,inner,outer,time,residual_eps);
                // fclose(resEPS);
                printf("residual_eps at iteration %d at outer %d at time %d = % 12.5e\n",inner,outer,time,residual_eps);
            }   
            
        }

        
        for (i=0;i<NX;i++)
        {
            free(eps_w[i]);
            free(eps_s[i]);
            free(eps_n[i]);
            free(eps_e[i]);
            free(eps_p[i]);
        }
        free(eps_w);
        free(eps_s);
        free(eps_n);
        free(eps_e);
        free(eps_p);        

    // Outer iteration loop end
        
        for (i=0;i<NX;i++)
        {
            for (j=0;j<NY;j++)
            {
                mut_C[i][j] = RHO_C*C_MU*pow(k[i][j],2)/epsilon[i][j] + 0.6*RHO_C*alpha_d[i][j]*DP*uslip_mag[i][j];
                mu_total_C[i][j] = mu + mut_C[i][j];
                mu_total_D[i][j] = mu_total_C[i][j]*RHO_D/RHO_C;
            }
        }

        if ((residual_x_c<TOL_X_C) && (residual_y_c<TOL_Y_C) && (residual_x_d<TOL_X_D) && (residual_y_d<TOL_Y_D) && (residual_P<RELATIVE_TOL) && (residual_k<TOL_KEPS) && (residual_eps<TOL_KEPS) && (residual_alpha<TOL_ALPHA))
        {
            printf("Solution converged after %d outer iterations\n",outer);
            break;
        }
            //// Outer iteration loop end
    }

    for (i=0;i<NX+1;i++)
    {
        for (j=0;j<NY;j++)
        {
            un_d[i][j] = u_d[i][j];
            un_c[i][j] = u_c[i][j];
        }
    }
    
    for (i=0;i<NX;i++)
    {
        for (j=0;j<NY+1;j++)
        {
            vn_d[i][j] = v_d[i][j];
            vn_c[i][j] = v_c[i][j];
        }
    }
        
    for (i=0;i<NX;i++)
    {
        for (j=0;j<NY;j++)
        {
            alphan_d[i][j] = alpha_d[i][j];
            alphan_c[i][j] = alpha_c[i][j];
            kn[i][j] = k[i][j];
            epsn[i][j] = epsilon[i][j];
        }
    }

    if ( (residual_y_d > 1e5) || (residual_x_d > 1e5) || (residual_y_c > 1e5) || (residual_x_c > 1e5) || (residual_P > 1e5) )
    {
        printf("Divergence Detected\nAborting Solution\n");
        exit(1);
    }

if (time % SAVING_INTERVAL == 0)
{
 
 char filename[100];


 //// Exporting u-velocity

    sprintf(filename,"uc%d.txt",time/SAVING_INTERVAL);
    uc[count]=fopen(filename,"a");

    for (i=0;i<NX+1;i++)
    {
        for (j=0;j<NY;j++)
        {
            fprintf(uc[count],"% 12.5e\t",u_c[i][j]);
        }
        fprintf(uc[count],"\n");
    }
    fclose(uc[count]);

    //// Exporting v-velocity
    
    sprintf(filename,"vc%d.txt",time/SAVING_INTERVAL);
    vc[count]=fopen(filename,"a");

    for (i=0;i<NX;i++)
    {
        for (j=0;j<NY+1;j++)
        {
            fprintf(vc[count],"% 12.5e\t",v_c[i][j]);
        }
        fprintf(vc[count],"\n");
    }
    fclose(vc[count]);
    
    //// Exporting u-velocity

    sprintf(filename,"ud%d.txt",time/SAVING_INTERVAL);
    ud[count]=fopen(filename,"a");

    for (i=0;i<NX+1;i++)
    {
        for (j=0;j<NY;j++)
        {
            fprintf(ud[count],"% 12.5e\t",u_d[i][j]);
        }
        fprintf(ud[count],"\n");
    }
    fclose(ud[count]);
    
    //// Exporting v-velocity

    sprintf(filename,"vd%d.txt",time/SAVING_INTERVAL);
    vd[count]=fopen(filename,"a");
    
    for (i=0;i<NX;i++)
    {
        for (j=0;j<NY+1;j++)
        {
            fprintf(vd[count],"% 12.5e\t",v_d[i][j]);
        }
        fprintf(vd[count],"\n");
    }
    fclose(vd[count]);

    
    //// Exporting Pressure

    sprintf(filename,"p%d.txt",time/SAVING_INTERVAL);
    Pr[count]=fopen(filename,"a");
    
    for (i=0;i<NX;i++)
    {
        for (j=0;j<NY;j++)
        {
            fprintf(Pr[count],"% 12.5e\t",P[i][j]);
        }
        fprintf(Pr[count],"\n");
    }
    fclose(Pr[count]);


    //// Exporting k

    sprintf(filename,"k%d.txt",time/SAVING_INTERVAL);
    Kfile[count]=fopen(filename,"a");

    for (i=0;i<NX;i++)
    {
        for (j=0;j<NY;j++)
        {
            fprintf(Kfile[count],"% 12.5e\t",k[i][j]);
        }
        fprintf(Kfile[count],"\n");
    }
    fclose(Kfile[count]);

    //// Exporting Epsilon

    sprintf(filename,"eps%d.txt",time/SAVING_INTERVAL);
    EpsFile[count]=fopen(filename,"a");

    for (i=0;i<NX;i++)
    {
        for (j=0;j<NY;j++)
        {
            fprintf(EpsFile[count],"% 12.5e\t",epsilon[i][j]);
        }
        fprintf(EpsFile[count],"\n");
    }
    fclose(EpsFile[count]);

    //// Exporting Turbulent Viscosity

    sprintf(filename,"mut%d.txt",time/SAVING_INTERVAL);
    mutC[count]=fopen(filename,"a");

    for (i=0;i<NX;i++)
    {
        for (j=0;j<NY;j++)
        {
            fprintf(mutC[count],"% 12.5e\t",mut_C[i][j]);
        }
        fprintf(mutC[count],"\n");
    }
    fclose(mutC[count]);

    //// Exporting Volume Fraction -- Disperse Phase

    sprintf(filename,"alphad%d.txt",time/SAVING_INTERVAL);
    alphad[count]=fopen(filename,"a");

    for (i=0;i<NX;i++)
    {
        for (j=0;j<NY;j++)
        {
            fprintf(alphad[count],"% 12.5e\t",alpha_d[i][j]);
        }
        fprintf(alphad[count],"\n");
    }
    fclose(alphad[count]);

    //// Exporting Volume Fraction -- Continuous Phase

    sprintf(filename,"alphac%d.txt",time/SAVING_INTERVAL);
    alphac[count]=fopen(filename,"a");

    for (i=0;i<NX;i++)
    {
        for (j=0;j<NY;j++)
        {
            fprintf(alphac[count],"% 12.5e\t",alpha_c[i][j]);
        }
        fprintf(alphac[count],"\n");
    }
    fclose(alphac[count]);


    count++;
 }

//// Time loop end
}

// end = clock();

// printf("%lf\n", ((double)(end-start))/CLOCKS_PER_SEC);
//// main function ends
    return 0;
}







