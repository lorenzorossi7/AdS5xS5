//============================================================================
// in polar coordinates t,x,y for x in [0,1] and y in [0,1]
// using r=x/(1-x) compactification and \chi=\pi y rescaling 
//
// dynamical coordinates are t,x,y
//
// coordinate labels are (1:t),(2:x),(3:y),(4,5,6,:theta1,theta2,theta3),
// (7,8,9,10:phi2,phi3,phi4,phi5) 
// AdS_5 sector (t,x,theta1,theta2,theta3)
// S^5 sector (y,phi2,phi3,phi4,phi5)
//
// source functions are full ten-dimensional source functions
// H_a = H3_a + dA/2 gA^{-1} g_A,a + dB/2 gB^{-1} g_B,a
// where dA=3 for the 3-sphere, dB=4 for the remainder of the 5-sphere
//
// application interface functions for AdS5xS5
//=============================================================================

#include <stdlib.h>
#include <stdio.h>
#include <pamr.h>
#include <amrd.h>
#include <math.h>
#include <m_util_r8.h>
#include <bbhutil.h>
#include <mpi.h>
#include "AdS5xS5.h"

//=============================================================================
// set in fparam for now
//=============================================================================
real AdS_L;

//=============================================================================
// Carsten's "constraint-damping" parameters
//=============================================================================
real kappa_cd,rho_cd;

//=============================================================================
// id (and other) parameters
//=============================================================================

// gaussians
real phi1_amp_1,phi1_r0_1,phi1_delta_1,phi1_x0_1[2],phi1_width_1[2],phi1_amp2_1,phi1_r02_1,phi1_delta2_1,phi1_x02_1[2],phi1_width2_1[2];
real phi1_amp_2,phi1_r0_2,phi1_delta_2,phi1_x0_2[2],phi1_width_2[2],phi1_amp2_2,phi1_r02_2,phi1_delta2_2,phi1_x02_2[2],phi1_width2_2[2];
real phi1_amp_3,phi1_r0_3,phi1_delta_3,phi1_x0_3[2],phi1_width_3[2],phi1_amp2_3,phi1_r02_3,phi1_delta2_3,phi1_x02_3[2],phi1_width2_3[2];

// if nonzero, initialize with exact BH
real ief_bh_r0;

// excision parameters
real ex_rbuf[MAX_BHS];

// "internal" excision parameters, set by AH finder 
real ex_r[MAX_BHS][2],ex_xc[MAX_BHS][2];

int background,skip_constraints;

// new parameters in rtfile
int interptype,i_shift;

//gauge parameters
int gauge_t;
int gauge_i;
real rho1_t,rho2_t,xi1_t,xi2_t;
real rho1_i,rho2_i,xi1_i,xi2_i;
real rhoa,rhob;
real a1,a2;
real b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12;
real c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13;

int cp_version; 

//=============================================================================
// some convenient, "local" global variables
//=============================================================================

real *phi1,*phi1_n,*phi1_np1,*phi1_nm1; // MGH, AMRH n/np1/nm1
real *phi1_t,*phi1_t_n;

real *gb_tt,*gb_tt_n,*gb_tt_np1,*gb_tt_nm1;
real *gb_tx,*gb_tx_n,*gb_tx_np1,*gb_tx_nm1;
real *gb_xx,*gb_xx_n,*gb_xx_np1,*gb_xx_nm1;
real *psi,*psi_n,*psi_np1,*psi_nm1;
real *omega,*omega_n,*omega_np1,*omega_nm1;
real *gb_ty,*gb_ty_n,*gb_ty_np1,*gb_ty_nm1;
real *gb_xy,*gb_xy_n,*gb_xy_np1,*gb_xy_nm1;
real *gb_yy,*gb_yy_n,*gb_yy_np1,*gb_yy_nm1;
real *gb_tt_t,*gb_tt_t_n;
real *gb_tx_t,*gb_tx_t_n;
real *gb_xx_t,*gb_xx_t_n;
real *psi_t,*psi_t_n;
real *gb_ty_t,*gb_ty_t_n;
real *gb_xy_t,*gb_xy_t_n;
real *gb_yy_t,*gb_yy_t_n;
real *omega_t,*omega_t_n;
real *gb_res;
real *cl_res;

real *Hb_t,*Hb_t_n,*Hb_t_np1,*Hb_t_nm1;
real *Hb_x,*Hb_x_n,*Hb_x_np1,*Hb_x_nm1;
real *Hb_y,*Hb_y_n,*Hb_y_np1,*Hb_y_nm1;
real *Hb_t_t,*Hb_t_t_n;
real *Hb_x_t,*Hb_x_t_n;
real *Hb_y_t,*Hb_y_t_n;
real *hb_t_res,*hb_i_res;

real *Hb_t_0,*Hb_x_0,*Hb_y_0;

real *zetab,*zetab_res,*zetab_lop,*zetab_rhs;

real *w1,*mg_w1;
real *w2,*mg_w2;
real *w3,*mg_w3;
real *w4,*mg_w4;

real *mask,*mask_mg,*chr,*chr_mg;

real *kg_ires,*iresphi1;
real *efe_all_ires,*iresall;

real *g_norms;

real *x,*y;
int shape[2],ghost_width[4],Nx,Ny,phys_bdy[4],size,g_rank;
real base_bbox[4],bbox[4],dx,dy,dt,dx_Lc,dy_Lc;
int g_L;

int phi1_gfn,phi1_n_gfn,phi1_np1_gfn,phi1_nm1_gfn; 
int phi1_t_gfn,phi1_t_n_gfn;

int gb_tt_gfn,gb_tt_n_gfn,gb_tt_np1_gfn,gb_tt_nm1_gfn;
int gb_tx_gfn,gb_tx_n_gfn,gb_tx_np1_gfn,gb_tx_nm1_gfn;
int gb_xx_gfn,gb_xx_n_gfn,gb_xx_np1_gfn,gb_xx_nm1_gfn;
int psi_gfn,psi_n_gfn,psi_np1_gfn,psi_nm1_gfn;
int omega_gfn,omega_n_gfn,omega_np1_gfn,omega_nm1_gfn;
int gb_ty_gfn,gb_ty_n_gfn,gb_ty_np1_gfn,gb_ty_nm1_gfn;
int gb_xy_gfn,gb_xy_n_gfn,gb_xy_np1_gfn,gb_xy_nm1_gfn;
int gb_yy_gfn,gb_yy_n_gfn,gb_yy_np1_gfn,gb_yy_nm1_gfn;
int gb_tt_t_gfn,gb_tt_t_n_gfn;
int gb_tx_t_gfn,gb_tx_t_n_gfn;
int gb_xx_t_gfn,gb_xx_t_n_gfn;
int psi_t_gfn,psi_t_n_gfn;
int gb_ty_t_gfn,gb_ty_t_n_gfn;
int gb_xy_t_gfn,gb_xy_t_n_gfn;
int gb_yy_t_gfn,gb_yy_t_n_gfn;
int omega_t_gfn,omega_t_n_gfn;
int gb_res_gfn;
int cl_res_gfn;

int Hb_t_gfn,Hb_t_n_gfn,Hb_t_np1_gfn,Hb_t_nm1_gfn;
int Hb_x_gfn,Hb_x_n_gfn,Hb_x_np1_gfn,Hb_x_nm1_gfn;
int Hb_y_gfn,Hb_y_n_gfn,Hb_y_np1_gfn,Hb_y_nm1_gfn;
int Hb_t_t_gfn,Hb_t_t_n_gfn;
int Hb_x_t_gfn,Hb_x_t_n_gfn;
int Hb_y_t_gfn,Hb_y_t_n_gfn;
int hb_t_res_gfn,hb_i_res_gfn;

int Hb_t_0_gfn,Hb_x_0_gfn,Hb_y_0_gfn;	

int zetab_gfn,zetab_res_gfn,zetab_lop_gfn,zetab_rhs_gfn;

int w1_gfn,mg_w1_gfn;
int w2_gfn,mg_w2_gfn;
int w3_gfn,mg_w3_gfn;
int w4_gfn,mg_w4_gfn;

int mask_gfn,mask_mg_gfn,chr_gfn,chr_mg_gfn;

int kg_ires_gfn,iresphi1_gfn;
int efe_all_ires_gfn,iresall_gfn;

//=============================================================================
// call after variables have been defined
//=============================================================================
void set_gfns(void)
{
    if ((phi1_gfn     = PAMR_get_gfn("phi1",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((phi1_nm1_gfn = PAMR_get_gfn("phi1",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((phi1_n_gfn   = PAMR_get_gfn("phi1",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((phi1_np1_gfn = PAMR_get_gfn("phi1",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((phi1_t_gfn   = PAMR_get_gfn("phi1_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((phi1_t_n_gfn = PAMR_get_gfn("phi1_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

    if ((gb_tt_gfn     = PAMR_get_gfn("gb_tt",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tt_nm1_gfn = PAMR_get_gfn("gb_tt",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tt_n_gfn   = PAMR_get_gfn("gb_tt",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tt_np1_gfn = PAMR_get_gfn("gb_tt",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tx_gfn     = PAMR_get_gfn("gb_tx",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tx_nm1_gfn = PAMR_get_gfn("gb_tx",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tx_n_gfn   = PAMR_get_gfn("gb_tx",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tx_np1_gfn = PAMR_get_gfn("gb_tx",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xx_gfn     = PAMR_get_gfn("gb_xx",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xx_nm1_gfn = PAMR_get_gfn("gb_xx",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xx_n_gfn   = PAMR_get_gfn("gb_xx",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xx_np1_gfn = PAMR_get_gfn("gb_xx",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((psi_gfn       = PAMR_get_gfn("psi",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((psi_nm1_gfn   = PAMR_get_gfn("psi",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((psi_n_gfn     = PAMR_get_gfn("psi",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((psi_np1_gfn   = PAMR_get_gfn("psi",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((omega_gfn     = PAMR_get_gfn("omega",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((omega_nm1_gfn = PAMR_get_gfn("omega",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((omega_n_gfn   = PAMR_get_gfn("omega",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((omega_np1_gfn = PAMR_get_gfn("omega",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_ty_gfn     = PAMR_get_gfn("gb_ty",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_ty_nm1_gfn = PAMR_get_gfn("gb_ty",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_ty_n_gfn   = PAMR_get_gfn("gb_ty",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_ty_np1_gfn = PAMR_get_gfn("gb_ty",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xy_gfn     = PAMR_get_gfn("gb_xy",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xy_nm1_gfn = PAMR_get_gfn("gb_xy",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xy_n_gfn   = PAMR_get_gfn("gb_xy",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xy_np1_gfn = PAMR_get_gfn("gb_xy",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_yy_gfn     = PAMR_get_gfn("gb_yy",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_yy_nm1_gfn = PAMR_get_gfn("gb_yy",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_yy_n_gfn   = PAMR_get_gfn("gb_yy",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_yy_np1_gfn = PAMR_get_gfn("gb_yy",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tt_t_gfn   = PAMR_get_gfn("gb_tt_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tt_t_n_gfn = PAMR_get_gfn("gb_tt_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tx_t_gfn   = PAMR_get_gfn("gb_tx_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_tx_t_n_gfn = PAMR_get_gfn("gb_tx_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xx_t_gfn   = PAMR_get_gfn("gb_xx_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xx_t_n_gfn = PAMR_get_gfn("gb_xx_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((psi_t_gfn     = PAMR_get_gfn("psi_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((psi_t_n_gfn   = PAMR_get_gfn("psi_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_ty_t_gfn   = PAMR_get_gfn("gb_ty_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_ty_t_n_gfn = PAMR_get_gfn("gb_ty_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xy_t_gfn   = PAMR_get_gfn("gb_xy_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_xy_t_n_gfn = PAMR_get_gfn("gb_xy_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_yy_t_gfn   = PAMR_get_gfn("gb_yy_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_yy_t_n_gfn = PAMR_get_gfn("gb_yy_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((omega_t_gfn   = PAMR_get_gfn("omega_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((omega_t_n_gfn = PAMR_get_gfn("omega_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gb_res_gfn    = PAMR_get_gfn("gb_res",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((cl_res_gfn    = PAMR_get_gfn("cl_res",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((Hb_t_gfn      = PAMR_get_gfn("Hb_t",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_t_nm1_gfn  = PAMR_get_gfn("Hb_t",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_t_n_gfn    = PAMR_get_gfn("Hb_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_t_np1_gfn  = PAMR_get_gfn("Hb_t",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_x_gfn      = PAMR_get_gfn("Hb_x",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_x_nm1_gfn  = PAMR_get_gfn("Hb_x",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_x_n_gfn    = PAMR_get_gfn("Hb_x",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_x_np1_gfn  = PAMR_get_gfn("Hb_x",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_y_gfn      = PAMR_get_gfn("Hb_y",PAMR_MGH, 0))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_y_nm1_gfn  = PAMR_get_gfn("Hb_y",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_y_n_gfn    = PAMR_get_gfn("Hb_y",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_y_np1_gfn  = PAMR_get_gfn("Hb_y",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_t_t_gfn  = PAMR_get_gfn("Hb_t_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_t_t_n_gfn  = PAMR_get_gfn("Hb_t_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_x_t_gfn  = PAMR_get_gfn("Hb_x_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_x_t_n_gfn  = PAMR_get_gfn("Hb_x_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_y_t_gfn  = PAMR_get_gfn("Hb_y_t",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_y_t_n_gfn  = PAMR_get_gfn("Hb_y_t",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((hb_t_res_gfn  = PAMR_get_gfn("hb_t_res",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((hb_i_res_gfn  = PAMR_get_gfn("hb_i_res",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((Hb_t_0_gfn  = PAMR_get_gfn("Hb_t_0",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_x_0_gfn  = PAMR_get_gfn("Hb_x_0",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((Hb_y_0_gfn  = PAMR_get_gfn("Hb_y_0",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((zetab_gfn     = PAMR_get_gfn("zetab",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((zetab_res_gfn = PAMR_get_gfn("zetab_res",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((zetab_lop_gfn = PAMR_get_gfn("zetab_lop",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((zetab_rhs_gfn = PAMR_get_gfn("zetab_rhs",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);

    if ((w1_gfn   = PAMR_get_gfn("w1",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((w2_gfn   = PAMR_get_gfn("w2",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((w3_gfn   = PAMR_get_gfn("w3",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((w4_gfn   = PAMR_get_gfn("w4",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((mg_w1_gfn   = PAMR_get_gfn("mg_w1",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((mg_w2_gfn   = PAMR_get_gfn("mg_w2",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((mg_w3_gfn   = PAMR_get_gfn("mg_w3",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((mg_w4_gfn   = PAMR_get_gfn("mg_w4",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);

    if ((kg_ires_gfn= PAMR_get_gfn("kg_ires",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((iresphi1_gfn  = PAMR_get_gfn("iresphi1",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((efe_all_ires_gfn= PAMR_get_gfn("efe_all_ires",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((iresall_gfn  = PAMR_get_gfn("iresall",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((mask_mg_gfn = PAMR_get_gfn("cmask",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((mask_gfn    = PAMR_get_gfn("cmask",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((chr_gfn     = PAMR_get_gfn("chr",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((chr_mg_gfn  = PAMR_get_gfn("chr",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);

    g_norms=AMRD_get_global_norms();
}

//=============================================================================
// call with valid iter to set up globals:
//=============================================================================
void ldptr_bbox(void)
{
   real dx0[2];
   static int first=1;

   if (first) 
   {
      first=0; 
      set_gfns();
      PAMR_get_global_bbox(base_bbox);
      if (PAMR_get_max_lev(PAMR_AMRH)>1) PAMR_get_dxdt(2,dx0,&dt); else PAMR_get_dxdt(1,dx0,&dt);
      dx_Lc=dx0[0];
      dy_Lc=dx0[1];
   }

   PAMR_get_g_rank(&g_rank);
   PAMR_get_g_shape(shape);
   PAMR_get_g_bbox(bbox);
   PAMR_get_g_ghost_width(ghost_width);
   PAMR_get_g_level(&g_L);
   PAMR_get_dxdt(g_L,dx0,&dt);
   dx=dx0[0];
   dy=dx0[1];

   if ((bbox[0]-base_bbox[0])<dx/2) phys_bdy[0]=1; else phys_bdy[0]=0;
   if ((base_bbox[1]-bbox[1])<dx/2) phys_bdy[1]=1; else phys_bdy[1]=0;
   if ((bbox[2]-base_bbox[2])<dy/2) phys_bdy[2]=1; else phys_bdy[2]=0;
   if ((base_bbox[3]-bbox[3])<dy/2) phys_bdy[3]=1; else phys_bdy[3]=0;

   Nx=shape[0];
   Ny=shape[1];

   size=Nx*Ny;
}

void ldptr(void)
{
   real *x0[2],*gfs[PAMR_MAX_GFNS];

   ldptr_bbox();

   PAMR_get_g_x(x0);

   x=x0[0];
   y=x0[1];

   PAMR_get_g_gfs(gfs);

   phi1     = gfs[phi1_gfn-1];
   phi1_n   = gfs[phi1_n_gfn-1];
   phi1_np1 = gfs[phi1_np1_gfn-1];
   phi1_nm1 = gfs[phi1_nm1_gfn-1];
   phi1_t   = gfs[phi1_t_gfn-1];
   phi1_t_n = gfs[phi1_t_n_gfn-1];

   gb_tt     = gfs[gb_tt_gfn-1];
   gb_tt_n   = gfs[gb_tt_n_gfn-1];
   gb_tt_np1 = gfs[gb_tt_np1_gfn-1];
   gb_tt_nm1 = gfs[gb_tt_nm1_gfn-1];
   gb_tx     = gfs[gb_tx_gfn-1];
   gb_tx_n   = gfs[gb_tx_n_gfn-1];
   gb_tx_np1 = gfs[gb_tx_np1_gfn-1];
   gb_tx_nm1 = gfs[gb_tx_nm1_gfn-1];
   gb_xx     = gfs[gb_xx_gfn-1];
   gb_xx_n   = gfs[gb_xx_n_gfn-1];
   gb_xx_np1 = gfs[gb_xx_np1_gfn-1];
   gb_xx_nm1 = gfs[gb_xx_nm1_gfn-1];
   psi       = gfs[psi_gfn-1];
   psi_n     = gfs[psi_n_gfn-1];
   psi_np1   = gfs[psi_np1_gfn-1];
   psi_nm1   = gfs[psi_nm1_gfn-1];
   omega     = gfs[omega_gfn-1];
   omega_n   = gfs[omega_n_gfn-1];
   omega_np1 = gfs[omega_np1_gfn-1];
   omega_nm1 = gfs[omega_nm1_gfn-1];
   gb_ty     = gfs[gb_ty_gfn-1];
   gb_ty_n   = gfs[gb_ty_n_gfn-1];
   gb_ty_np1 = gfs[gb_ty_np1_gfn-1];
   gb_ty_nm1 = gfs[gb_ty_nm1_gfn-1];
   gb_xy     = gfs[gb_xy_gfn-1];
   gb_xy_n   = gfs[gb_xy_n_gfn-1];
   gb_xy_np1 = gfs[gb_xy_np1_gfn-1];
   gb_xy_nm1 = gfs[gb_xy_nm1_gfn-1];
   gb_yy     = gfs[gb_yy_gfn-1];
   gb_yy_n   = gfs[gb_yy_n_gfn-1];
   gb_yy_np1 = gfs[gb_yy_np1_gfn-1];
   gb_yy_nm1 = gfs[gb_yy_nm1_gfn-1];
   gb_tt_t   = gfs[gb_tt_t_gfn-1];
   gb_tt_t_n = gfs[gb_tt_t_n_gfn-1];
   gb_tx_t   = gfs[gb_tx_t_gfn-1];
   gb_tx_t_n = gfs[gb_tx_t_n_gfn-1];
   gb_xx_t   = gfs[gb_xx_t_gfn-1];
   gb_xx_t_n = gfs[gb_xx_t_n_gfn-1];
   psi_t     = gfs[psi_t_gfn-1];
   psi_t_n   = gfs[psi_t_n_gfn-1];
   gb_ty_t   = gfs[gb_ty_t_gfn-1];
   gb_ty_t_n = gfs[gb_ty_t_n_gfn-1];
   gb_xy_t   = gfs[gb_xy_t_gfn-1];
   gb_xy_t_n = gfs[gb_xy_t_n_gfn-1];
   gb_yy_t   = gfs[gb_yy_t_gfn-1];
   gb_yy_t_n = gfs[gb_yy_t_n_gfn-1];
   omega_t   = gfs[omega_t_gfn-1];
   omega_t_n = gfs[omega_t_n_gfn-1];
   gb_res    = gfs[gb_res_gfn-1];
   cl_res    = gfs[cl_res_gfn-1];

   Hb_t      = gfs[Hb_t_gfn-1];
   Hb_t_n    = gfs[Hb_t_n_gfn-1];
   Hb_t_nm1  = gfs[Hb_t_nm1_gfn-1];
   Hb_t_np1  = gfs[Hb_t_np1_gfn-1];
   Hb_x      = gfs[Hb_x_gfn-1];
   Hb_x_n    = gfs[Hb_x_n_gfn-1];
   Hb_x_nm1  = gfs[Hb_x_nm1_gfn-1];
   Hb_x_np1  = gfs[Hb_x_np1_gfn-1];
   Hb_y      = gfs[Hb_y_gfn-1];
   Hb_y_n    = gfs[Hb_y_n_gfn-1];
   Hb_y_nm1  = gfs[Hb_y_nm1_gfn-1];
   Hb_y_np1  = gfs[Hb_y_np1_gfn-1];
   Hb_t_t  = gfs[Hb_t_t_gfn-1];
   Hb_t_t_n  = gfs[Hb_t_t_n_gfn-1];
   Hb_x_t  = gfs[Hb_x_t_gfn-1];
   Hb_x_t_n  = gfs[Hb_x_t_n_gfn-1];
   Hb_y_t  = gfs[Hb_y_t_gfn-1];
   Hb_y_t_n  = gfs[Hb_y_t_n_gfn-1];
   hb_t_res  = gfs[hb_t_res_gfn-1];
   hb_i_res  = gfs[hb_i_res_gfn-1];

   Hb_t_0  = gfs[Hb_t_0_gfn-1];
   Hb_x_0  = gfs[Hb_x_0_gfn-1];
   Hb_y_0  = gfs[Hb_y_0_gfn-1];

   zetab     = gfs[zetab_gfn-1];
   zetab_lop = gfs[zetab_lop_gfn-1];
   zetab_res = gfs[zetab_res_gfn-1];
   zetab_rhs = gfs[zetab_rhs_gfn-1];

   w1 = gfs[w1_gfn-1];
   w2 = gfs[w2_gfn-1];
   w3 = gfs[w3_gfn-1];
   w4 = gfs[w4_gfn-1];

   mg_w1    =gfs[mg_w1_gfn-1]; 
   mg_w2    =gfs[mg_w2_gfn-1]; 
   mg_w3    =gfs[mg_w3_gfn-1]; 
   mg_w4    =gfs[mg_w4_gfn-1]; 

   kg_ires = gfs[kg_ires_gfn-1];
   iresphi1 = gfs[iresphi1_gfn-1];
   efe_all_ires = gfs[efe_all_ires_gfn-1];
   iresall = gfs[iresall_gfn-1];

   mask    = gfs[mask_gfn-1];
   mask_mg = gfs[mask_mg_gfn-1];
   chr = gfs[chr_gfn-1]; 
   chr_mg = gfs[chr_mg_gfn-1]; 
}

//=============================================================================
// PAMR_get_dxdt() only works with AMR hierarchy levels ... here we use
// lambda for dt, but this only works if rhosp=rhotm
//=============================================================================
void ldptr_mg(void)
{
   real lambda;

   ldptr();

   dx=x[1]-x[0]; dy=y[1]-y[0];
   PAMR_get_lambda(&lambda);
   dt=lambda*dx;
}

//=============================================================================
// utility routines
//=============================================================================
void const_f(real *f, real c)
{
   int i;

   for (i=0; i<Nx*Ny; i++) f[i]=c;
}

void zero_f(real *f)
{
   const_f(f,0);
}

void zero_f_ex(real *f, real *chr)
{
   int i;

   for (i=0; i<Nx*Ny; i++) if (chr[i]==AMRD_ex) f[i]=0;
}

real norm_l2(real *f, real *cmask, real *chr)
{
   int i;
   real norm=0;
   int sum=0;

   for (i=0; i<Nx*Ny; i++) 
      if (cmask[i]==AMRD_CMASK_ON && (chr[i]!=AMRD_ex)) { sum++; norm+=f[i]*f[i]; }

   if (!sum) sum=1;
   return (sqrt(norm/sum));
}

//=============================================================================
// Routines required by amrd:
//=============================================================================

//=============================================================================
// Returns 0 to use default mechanism, or is expected to calculate
// the correct initial hierarchy and return 1:
//=============================================================================
int AdS5xS5_id(void)
{
   return 0;
}

//=============================================================================
// Sets custom parameters, variables, etc. Split up into two segments,
// one called before the pamr context is initialized and standard
// parameters are read, and the other afterwards
//=============================================================================
void AdS5xS5_var_pre_init(char *pfile)
{
   AMRD_echo_params=1;
   AMRD_int_param(pfile,"echo_params",&AMRD_echo_params,1);

   cp_version=ADS5DP_CP_VERSION;
   AMRD_int_param(pfile,"cp_version",&cp_version,1);
   
   AdS_L=1.0; AMRD_real_param(pfile,"AdS_L",&AdS_L,1);

   interptype=0; AMRD_int_param(pfile,"interptype",&interptype,1);
   i_shift=0; AMRD_int_param(pfile,"i_shift",&i_shift,1);

   return;
}

void AdS5xS5_var_post_init(char *pfile)
{
   if (my_rank==0)
   {
      printf("===================================================================\n");
      printf("Reading AdS5xS5 parameters:\n\n");
   }

   phi1_amp_1=phi1_r0_1=phi1_x0_1[0]=phi1_x0_1[1]=phi1_amp2_1=phi1_r02_1=phi1_x02_1[0]=phi1_x02_1[1]=0;
   phi1_amp_2=phi1_r0_2=phi1_x0_2[0]=phi1_x0_2[1]=phi1_amp2_2=phi1_r02_2=phi1_x02_2[0]=phi1_x02_2[1]=0;
   phi1_amp_3=phi1_r0_3=phi1_x0_3[0]=phi1_x0_3[1]=phi1_amp2_3=phi1_r02_3=phi1_x02_3[0]=phi1_x02_3[1]=0;
   phi1_width_1[0]=phi1_width_1[1]=phi1_width2_1[0]=phi1_width2_1[1]=1;
   phi1_width_2[0]=phi1_width_2[1]=phi1_width2_2[0]=phi1_width2_2[1]=1;
   phi1_width_3[0]=phi1_width_3[1]=phi1_width2_3[0]=phi1_width2_3[1]=1;

   AMRD_real_param(pfile,"phi1_amp_1",&phi1_amp_1,1);
   AMRD_real_param(pfile,"phi1_r0_1",&phi1_r0_1,1);
   AMRD_real_param(pfile,"phi1_delta_1",&phi1_delta_1,1);
   AMRD_real_param(pfile,"phi1_x0_1",phi1_x0_1,AMRD_dim);
   AMRD_real_param(pfile,"phi1_width_1",phi1_width_1,AMRD_dim);
   AMRD_real_param(pfile,"phi1_amp2_1",&phi1_amp2_1,1);
   AMRD_real_param(pfile,"phi1_r02_1",&phi1_r02_1,1);
   AMRD_real_param(pfile,"phi1_delta2_1",&phi1_delta2_1,1);
   AMRD_real_param(pfile,"phi1_x02_1",phi1_x02_1,AMRD_dim);
   AMRD_real_param(pfile,"phi1_width2_1",phi1_width2_1,AMRD_dim);

   AMRD_real_param(pfile,"phi1_amp_2",&phi1_amp_2,1);
   AMRD_real_param(pfile,"phi1_r0_2",&phi1_r0_2,1);
   AMRD_real_param(pfile,"phi1_delta_2",&phi1_delta_2,1);
   AMRD_real_param(pfile,"phi1_x0_2",phi1_x0_2,AMRD_dim);
   AMRD_real_param(pfile,"phi1_width_2",phi1_width_2,AMRD_dim);
   AMRD_real_param(pfile,"phi1_amp2_2",&phi1_amp2_2,1);
   AMRD_real_param(pfile,"phi1_r02_2",&phi1_r02_2,1);
   AMRD_real_param(pfile,"phi1_delta2_2",&phi1_delta2_2,1);
   AMRD_real_param(pfile,"phi1_x02_2",phi1_x02_2,AMRD_dim);
   AMRD_real_param(pfile,"phi1_width2_2",phi1_width2_2,AMRD_dim);

   AMRD_real_param(pfile,"phi1_amp_3",&phi1_amp_3,1);
   AMRD_real_param(pfile,"phi1_r0_3",&phi1_r0_3,1);
   AMRD_real_param(pfile,"phi1_delta_3",&phi1_delta_3,1);
   AMRD_real_param(pfile,"phi1_x0_3",phi1_x0_3,AMRD_dim);
   AMRD_real_param(pfile,"phi1_width_3",phi1_width_3,AMRD_dim);
   AMRD_real_param(pfile,"phi1_amp2_3",&phi1_amp2_3,1);
   AMRD_real_param(pfile,"phi1_r02_3",&phi1_r02_3,1);
   AMRD_real_param(pfile,"phi1_delta2_3",&phi1_delta2_3,1);
   AMRD_real_param(pfile,"phi1_x02_3",phi1_x02_3,AMRD_dim);
   AMRD_real_param(pfile,"phi1_width2_3",phi1_width2_3,AMRD_dim);

   kappa_cd=0; AMRD_real_param(pfile,"kappa_cd",&kappa_cd,1);
   rho_cd=0; AMRD_real_param(pfile,"rho_cd",&rho_cd,1);

   background=0; AMRD_int_param(pfile,"background",&background,1);
   skip_constraints=0; AMRD_int_param(pfile,"skip_constraints",&skip_constraints,1);

   gauge_t=0; AMRD_int_param(pfile,"gauge_t",&gauge_t,1);
   gauge_i=0; AMRD_int_param(pfile,"gauge_i",&gauge_i,1);
   rho1_t=1; AMRD_real_param(pfile,"rho1_t",&rho1_t,1);
   rho2_t=1; AMRD_real_param(pfile,"rho2_t",&rho2_t,1);
   xi1_t=1; AMRD_real_param(pfile,"xi1_t",&xi1_t,1);
   xi2_t=1; AMRD_real_param(pfile,"xi2_t",&xi2_t,1);
   rho1_i=1; AMRD_real_param(pfile,"rho1_i",&rho1_i,1);
   rho2_i=1; AMRD_real_param(pfile,"rho2_i",&rho2_i,1);
   xi1_i=1; AMRD_real_param(pfile,"xi1_i",&xi1_i,1);
   xi2_i=1; AMRD_real_param(pfile,"xi2_i",&xi2_i,1);
   a1=1; AMRD_real_param(pfile,"a1",&a1,1);
   a2=0; AMRD_real_param(pfile,"a2",&a2,1);
   b1=1; AMRD_real_param(pfile,"b1",&b1,1);
   b2=1; AMRD_real_param(pfile,"b2",&b2,1);
   b3=1; AMRD_real_param(pfile,"b3",&b3,1);
   b4=1; AMRD_real_param(pfile,"b4",&b4,1);
   b5=1; AMRD_real_param(pfile,"b5",&b5,1);
   b6=1; AMRD_real_param(pfile,"b6",&b6,1);
   b7=1; AMRD_real_param(pfile,"b7",&b7,1);
   b8=1; AMRD_real_param(pfile,"b8",&b8,1);
   b9=1; AMRD_real_param(pfile,"b9",&b9,1);
   b10=1; AMRD_real_param(pfile,"b10",&b10,1);
   b11=1; AMRD_real_param(pfile,"b11",&b11,1);
   b12=1; AMRD_real_param(pfile,"b12",&b12,1);
   c1=0; AMRD_real_param(pfile,"c1",&c1,1);
   c2=0; AMRD_real_param(pfile,"c2",&c2,1);
   c3=0; AMRD_real_param(pfile,"c3",&c3,1);
   c4=0; AMRD_real_param(pfile,"c4",&c4,1);
   c5=0; AMRD_real_param(pfile,"c5",&c5,1);
   c6=0; AMRD_real_param(pfile,"c6",&c6,1);
   c7=0; AMRD_real_param(pfile,"c7",&c7,1);
   c8=0; AMRD_real_param(pfile,"c8",&c8,1);
   c9=0; AMRD_real_param(pfile,"c9",&c9,1);
   c10=0; AMRD_real_param(pfile,"c10",&c10,1);
   c11=0; AMRD_real_param(pfile,"c11",&c11,1);
   c12=0; AMRD_real_param(pfile,"c12",&c12,1);
   c13=0; AMRD_real_param(pfile,"c13",&c13,1);

   rhoa=1; AMRD_real_param(pfile,"rhoa",&rhoa,1);
   rhob=1; AMRD_real_param(pfile,"rhob",&rhob,1);

   // set fraction, 1-ex_rbuf, of AH radius to be excised
   int j;
   char buf[64];
   for (j=0; j<MAX_BHS; j++)
   {
      if (j==0) { if (!AMRD_cp_restart) ex_rbuf[j]=0; sprintf(buf,"ex_rbuf"); }
      else { if (!AMRD_cp_restart) ex_rbuf[j]=ex_rbuf[0]; sprintf(buf,"ex_rbuf_%i",j+1); }
      AMRD_real_param(pfile,buf,&ex_rbuf[j],1);
      if (ex_rbuf[j]<0 || ex_rbuf[j]>1 ) printf("WARNING ... ex_rbuf[%i]=%lf is outside of standard bbox\n",j,ex_rbuf[j]);
   }

   ief_bh_r0=0; AMRD_real_param(pfile,"ief_bh_r0",&ief_bh_r0,1);

   // (the following is for a single analytic BH j=0)
   // xh is radius in uncompcatified coordicates, xh_c in compactified,
   // ief_bh_r0 is BH radius parameter, ex_r[c_BH][0] is excision radius
   if (ief_bh_r0!=0)
   {
     real xh,xh_c;
     xh=AdS_L*sqrt(sqrt(1+4*ief_bh_r0*ief_bh_r0/AdS_L/AdS_L)/2-0.5);
     xh_c=xh/(1+xh);
     ex_r[0][0]=xh_c*(1-ex_rbuf[0]);
     if (my_rank==0) printf("\nBH initial data\n xh/L=%lf\n"
                            "Initial BH radius=%lf, (%lf in compactified (code) coords)\n"
                            "Initial excision radius=%lf\n\n",xh/AdS_L,xh,xh_c,ex_r[0][0]);
   }

   if (AMRD_do_ex==0) AMRD_stop("require excision to be on","");

   PAMR_excision_on("chr",&AdS5xS5_fill_ex_mask,AMRD_ex,1);

   if (my_rank==0) printf("===================================================================\n");
   return;
}

//=============================================================================
// Sets all variables to their 'zero' values:
//=============================================================================
void AdS5xS5_AMRH_var_clear(void)
{
   ldptr();

   zero_f(phi1_n); 

   return;
}

//=============================================================================
// Initial data for free fields: (at tn=2) ... following vars also used in 
// t0_cnst_data
//=============================================================================
void AdS5xS5_free_data(void)
{
   int i;

   ldptr();

   AdS5xS5_AMRH_var_clear(); // constrained variables are set post-MG

   zero_f(phi1_t_n); // sets initial time derivatives for ID
   zero_f(gb_tt_t_n);
   zero_f(gb_tx_t_n);
   zero_f(gb_ty_t_n);
   zero_f(gb_xx_t_n);
   zero_f(gb_xy_t_n);
   zero_f(gb_yy_t_n);
   zero_f(psi_t_n);
   zero_f(omega_t_n);
   zero_f(Hb_t_t_n);
   zero_f(Hb_x_t_n);
   zero_f(Hb_y_t_n);

   gauss2d_(phi1_n,
            &phi1_amp_1,&phi1_r0_1,&phi1_delta_1,&phi1_x0_1[0],&phi1_x0_1[0],&phi1_width_1[0],&phi1_width_1[0],
            &phi1_amp2_1,&phi1_r02_1,&phi1_delta2_1,&phi1_x02_1[0],&phi1_x02_1[1],&phi1_width2_1[0],&phi1_width2_1[1],
            &AdS_L,x,y,&Nx,&Ny);

   gauss2d_(w1,
            &phi1_amp_2,&phi1_r0_2,&phi1_delta_2,&phi1_x0_2[0],&phi1_x0_2[1],&phi1_width_2[0],&phi1_width_2[1],
            &phi1_amp2_2,&phi1_r02_2,&phi1_delta2_2,&phi1_x02_2[0],&phi1_x02_2[1],&phi1_width2_2[0],&phi1_width2_2[1],
            &AdS_L,x,y,&Nx,&Ny);

   gauss2d_(w2,
            &phi1_amp_3,&phi1_r0_3,&phi1_delta_3,&phi1_x0_3[0],&phi1_x0_3[1],&phi1_width_3[0],&phi1_width_3[1],
            &phi1_amp2_3,&phi1_r02_3,&phi1_delta2_3,&phi1_x02_3[0],&phi1_x02_3[1],&phi1_width2_3[0],&phi1_width2_3[1],
            &AdS_L,x,y,&Nx,&Ny);

   for (i=0; i<size; i++) phi1_n[i]+=w1[i]+w2[i]; 

   return;
}  

//=============================================================================
// Initialize any "elliptic_vars_t0" post construction of MGH, but before
// the start of vcycling.
//=============================================================================
void AdS5xS5_elliptic_vars_t0_init(void)
{
   // initializes dt, dx
   ldptr_mg();

   // initializes zetab conformal factor variable to zero
   const_f(zetab,0);
}

//=============================================================================
// Initial constraint data --- called after each MG iteration.
//
// Here we also initialize past time level information if
// AMRD_id_pl_method==3
//
// NOTE: np1,n,nm1 variables are allocated only at the top level of the MG hierarchy, 
//       so do an if(f_nm1){...}, for example, to make sure we are at the top level
//=============================================================================
void AdS5xS5_t0_cnst_data(void)
{
   int i;

   ldptr_mg();

   // initialize gbars
   if (ief_bh_r0==0 && skip_constraints==0)
   {
     init_ghb_(zetab,phi1,gb_tt,gb_tx,gb_ty,gb_xx,gb_xy,gb_yy,psi,omega,
               &rhoa,&rhob,&AdS_L,phys_bdy,chr_mg,&AMRD_ex,x,y,&Nx,&Ny);
   }
   else
   {
     init_schw_(gb_tt,gb_tx,gb_ty,gb_xx,gb_xy,gb_yy,psi,omega,&ief_bh_r0,
                &AdS_L,phys_bdy,chr_mg,&AMRD_ex,x,y,&Nx,&Ny);
   }   

   // initialize nm1,np1 time levels and hbars
   if (AMRD_id_pl_method==3 && phi1_nm1)
   {

     init_hb_(gb_tt_np1,gb_tt_n,gb_tt_nm1,
              gb_tx_np1,gb_tx_n,gb_tx_nm1,
              gb_ty_np1,gb_ty_n,gb_ty_nm1,
              gb_xx_np1,gb_xx_n,gb_xx_nm1,
              gb_xy_np1,gb_xy_n,gb_xy_nm1,
              gb_yy_np1,gb_yy_n,gb_yy_nm1,
              psi_np1,psi_n,psi_nm1,
              omega_np1,omega_n,omega_nm1,
              Hb_t_n,Hb_x_n,Hb_y_n,
              &AdS_L,phys_bdy,x,y,&dt,chr,&AMRD_ex,&Nx,&Ny);

     init_nm1_(gb_tt_np1,gb_tt_n,gb_tt_nm1,gb_tt_t_n,
               gb_tx_np1,gb_tx_n,gb_tx_nm1,gb_tx_t_n,
               gb_ty_np1,gb_ty_n,gb_ty_nm1,gb_ty_t_n,
               gb_xx_np1,gb_xx_n,gb_xx_nm1,gb_xx_t_n,
               gb_xy_np1,gb_xy_n,gb_xy_nm1,gb_xy_t_n,
               gb_yy_np1,gb_yy_n,gb_yy_nm1,gb_yy_t_n,
               psi_np1,psi_n,psi_nm1,psi_t_n,
               omega_np1,omega_n,omega_nm1,omega_t_n,
               Hb_t_np1,Hb_t_n,Hb_t_nm1,Hb_t_t_n,
               Hb_x_np1,Hb_x_n,Hb_x_nm1,Hb_x_t_n,
               Hb_y_np1,Hb_y_n,Hb_y_nm1,Hb_y_t_n,
               phi1_np1,phi1_n,phi1_nm1,phi1_t_n,
               &AdS_L,phys_bdy,x,y,&dt,chr,&AMRD_ex,&Nx,&Ny);
//     // straight copies 
//     for (i=0; i<size; i++)
//     {
//       gb_tt_np1[i]=gb_tt_nm1[i]=gb_tt[i];
//       gb_tx_np1[i]=gb_tx_nm1[i]=gb_tx[i];
//       gb_ty_np1[i]=gb_ty_nm1[i]=gb_ty[i];
//       gb_xx_np1[i]=gb_xx_nm1[i]=gb_xx[i];
//       gb_xy_np1[i]=gb_xy_nm1[i]=gb_xy[i];
//       gb_yy_np1[i]=gb_yy_nm1[i]=gb_yy[i];
//       psi_np1[i]=psi_nm1[i]=psi[i];
//       omega_np1[i]=omega_nm1[i]=omega[i];
//       Hb_t_np1[i]=Hb_t_nm1[i]=Hb_t_n[i];
//       Hb_x_np1[i]=Hb_x_nm1[i]=Hb_x_n[i];
//       Hb_y_np1[i]=Hb_y_nm1[i]=Hb_y_n[i];
//       phi1_np1[i]=phi1_nm1[i]=phi1[i];
//     }

     // store initial source functions
     for (i=0; i<size; i++)
     {
       Hb_t_0[i]=Hb_t[i];
       Hb_x_0[i]=Hb_x[i];
       Hb_y_0[i]=Hb_y[i];
     }

   }

   return;
}

//=============================================================================
// Calculations prior to saving info to disk.
//
// NOTE: at this point, the time sequence is: n,nm1,np1 (unless t=0)
//=============================================================================
void AdS5xS5_pre_io_calc(void)
{
   ldptr();

   int i,j;
   real ct;

   dx=x[1]-x[0]; dy=y[1]-y[0];
   ct=PAMR_get_time(g_L);

   // compute independent residuals of the AdS5D system
   if (ct!=0)
   {
     //(NOTE: for t>0, have cycled time sequence np1,n,nm1 to time sequence n,nm1,np1,
     // so here, time level n is the most advanced time level)
     ires_(efe_all_ires,kg_ires,
        gb_tt_n,gb_tt_nm1,gb_tt_np1,
        gb_tx_n,gb_tx_nm1,gb_tx_np1,
        gb_ty_n,gb_ty_nm1,gb_ty_np1,
        gb_xx_n,gb_xx_nm1,gb_xx_np1,
        gb_xy_n,gb_xy_nm1,gb_xy_np1,
        gb_yy_n,gb_yy_nm1,gb_yy_np1,
        psi_n,psi_nm1,psi_np1,
        omega_n,omega_nm1,omega_np1,
        phi1_n,phi1_nm1,phi1_np1,
        x,y,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,phys_bdy,ghost_width);
   }
   else
   {
     //(NOTE: for t=0, have *not* cycled time sequence, so still np1,n,nm1,
     // so here, time level np1 is the most advanced time level)
     ires_(efe_all_ires,kg_ires,
        gb_tt_np1,gb_tt_n,gb_tt_nm1,
        gb_tx_np1,gb_tx_n,gb_tx_nm1,
        gb_ty_np1,gb_ty_n,gb_ty_nm1,
        gb_xx_np1,gb_xx_n,gb_xx_nm1,
        gb_xy_np1,gb_xy_n,gb_xy_nm1,
        gb_yy_np1,gb_yy_n,gb_yy_nm1,
        psi_np1,psi_n,psi_nm1,
        omega_np1,omega_n,omega_nm1,
        phi1_np1,phi1_n,phi1_nm1,
        x,y,&dt,chr,&AdS_L,&AMRD_ex,&Nx,&Ny,phys_bdy,ghost_width);
   }

   // fill in ires arrays with independent residuals
   for (i=0; i<Nx; i++)
   {
      for (j=0; j<Ny; j++)
      {
         if (chr[i+j*Nx]==AMRD_ex) 
         {
            iresphi1[i+j*Nx]=0;  
            iresall[i+j*Nx]=0;
         }
         else
         {
            iresphi1[i+j*Nx]=kg_ires[i+j*Nx];
            iresall[i+j*Nx]=efe_all_ires[i+j*Nx];
         }
      }
   }

   return;
}

// to be callable from fortran
void check_nan_(real *x, int *is_nan)
{
   if (isnan(*x)) *is_nan=1; else *is_nan=0;
}

//=============================================================================
// Returns some norm of the residual for the evolution variables ... just
// use the value from the most recent evolution step
//=============================================================================
#define LIN_ZERO_BND 1
int lin_zero_bnd_all=1;
real AdS5xS5_evo_residual(void)
{
   real l2norm=0,l2norm_gb,l2norm_hb_t,l2norm_hb_i;
   int is_nan;

   ldptr();

   if (LIN_ZERO_BND) 
   {
      lin_zero_bnd_res_(gb_res,phys_bdy,&lin_zero_bnd_all,&Nx,&Ny);
      lin_zero_bnd_res_(hb_t_res,phys_bdy,&lin_zero_bnd_all,&Nx,&Ny);
      lin_zero_bnd_res_(hb_i_res,phys_bdy,&lin_zero_bnd_all,&Nx,&Ny);
   }

   l2norm_gb  =norm_l2(gb_res,mask,chr);
   l2norm_hb_t=norm_l2(hb_t_res,mask,chr);
   l2norm_hb_i=norm_l2(hb_i_res,mask,chr);

   l2norm=l2norm_gb+l2norm_hb_t+l2norm_hb_i;

   check_nan_(&l2norm,&is_nan);

   if (is_nan)
   {
      printf("\nl2norm_gb=%lf, l2norm_hb_t=%lf, l2norm_hb_i=%lf\n",
              l2norm_gb,l2norm_hb_t,l2norm_hb_i);
      printf("Nx=%i,Ny=%i,L=%i\n",Nx,Ny,g_L);
      AMRD_stop("l2norm is nan ... stopping","");
      l2norm=0;
   }

   return l2norm;
}

//=============================================================================
// Performs 1 iteration of the evolution equations 
//
// NOTE: at this point, the time sequence is: np1,n,nm1 
//=============================================================================
void AdS5xS5_evolve(int iter)
{
   int ltrace=0;
   real ct;

   ldptr();

   ct=PAMR_get_time(g_L);

   if (ltrace) printf("AdS5xS5_evolve: iter=%i , time=%lf, lev=%i, rank=%i\n",iter,ct,g_L,my_rank);

   hb_t_evo_(hb_t_res,
             gb_tt_np1,gb_tt_n,gb_tt_nm1,
             gb_tx_np1,gb_tx_n,gb_tx_nm1,
             gb_ty_np1,gb_ty_n,gb_ty_nm1,
             gb_xx_np1,gb_xx_n,gb_xx_nm1,
             gb_xy_np1,gb_xy_n,gb_xy_nm1,
             gb_yy_np1,gb_yy_n,gb_yy_nm1,
             psi_np1,psi_n,psi_nm1,
             omega_np1,omega_n,omega_nm1,
             Hb_t_np1,Hb_t_n,Hb_t_nm1,
             Hb_x_np1,Hb_x_n,Hb_x_nm1,
             Hb_y_np1,Hb_y_n,Hb_y_nm1,
             phi1_np1,phi1_n,phi1_nm1,
             &AdS_L,x,y,&dt,chr,&AMRD_ex,
             phys_bdy,ghost_width,&Nx,&Ny,
             Hb_t_0,Hb_x_0,Hb_y_0,
             &gauge_t,&ct,&rho1_t,&rho2_t,&xi1_t,&xi2_t,
             &a1,&a2);

   hb_i_evo_(hb_i_res,
             gb_tt_np1,gb_tt_n,gb_tt_nm1,
             gb_tx_np1,gb_tx_n,gb_tx_nm1,
             gb_ty_np1,gb_ty_n,gb_ty_nm1,
             gb_xx_np1,gb_xx_n,gb_xx_nm1,
             gb_xy_np1,gb_xy_n,gb_xy_nm1,
             gb_yy_np1,gb_yy_n,gb_yy_nm1,
             psi_np1,psi_n,psi_nm1,
             omega_np1,omega_n,omega_nm1,
             Hb_t_np1,Hb_t_n,Hb_t_nm1,
             Hb_x_np1,Hb_x_n,Hb_x_nm1,
             Hb_y_np1,Hb_y_n,Hb_y_nm1,
             phi1_np1,phi1_n,phi1_nm1,
             &AdS_L,x,y,&dt,chr,&AMRD_ex,
             phys_bdy,ghost_width,&Nx,&Ny,
             Hb_t_0,Hb_x_0,Hb_y_0,
             &gauge_i,&ct,&rho1_i,&rho2_i,&xi1_i,&xi2_i,
             &b1,&b2,&b3,&b4,&b5,&b6,&b7,&b8,&b9,&b10,&b11,&b12,
             &c1,&c2,&c3,&c4,&c5,&c6,&c7,&c8,&c9,&c10,&c11,&c12,&c13);

   g_evo_opt_(gb_res,cl_res,
              gb_tt_np1,gb_tt_n,gb_tt_nm1,
              gb_tx_np1,gb_tx_n,gb_tx_nm1,
              gb_ty_np1,gb_ty_n,gb_ty_nm1,
              gb_xx_np1,gb_xx_n,gb_xx_nm1,
              gb_xy_np1,gb_xy_n,gb_xy_nm1,
              gb_yy_np1,gb_yy_n,gb_yy_nm1,
              psi_np1,psi_n,psi_nm1,
              omega_np1,omega_n,omega_nm1,
              Hb_t_np1,Hb_t_n,Hb_t_nm1,
              Hb_x_np1,Hb_x_n,Hb_x_nm1,
              Hb_y_np1,Hb_y_n,Hb_y_nm1,
              phi1_np1,phi1_n,phi1_nm1,
              &AdS_L,x,y,&dt,chr,&AMRD_ex,
              phys_bdy,ghost_width,&Nx,&Ny,
              &background,&kappa_cd,&rho_cd);

   return;
}

//=============================================================================
// sets excision mask (NO ITERATOR, SO DON'T LOAD POINTERS!!!)
// 
// in polar code, only excised regions are those inside horizons
//=============================================================================
void AdS5xS5_fill_ex_mask(real *mask, int dim, int *shape, real *bbox, real excised)
{
   int i,j;
   int l,n,a,b;
   int axiregpts;
   int rho_sp[PAMR_MAX_LEVS],rho_tm[PAMR_MAX_LEVS];
   int Lf,rho_sp0,rho_tm0;
   real x,y,dx,dy;

   PAMR_get_rho(rho_sp,rho_tm,PAMR_MAX_LEVS);
   Lf=PAMR_get_max_lev(PAMR_AMRH);
   rho_sp0=rho_sp[Lf]; rho_tm0=rho_tm[Lf];

   dx=(bbox[1]-bbox[0])/(shape[0]-1);
   dy=(bbox[3]-bbox[2])/(shape[1]-1);

   for (i=0; i<shape[0]; i++)
   {
      for (j=0; j<shape[1]; j++)
      {
         x=bbox[0]+i*dx;
         y=bbox[1]+j*dy;
         mask[i+j*shape[0]]=excised-1;
         for (l=0; l<MAX_BHS; l++)
         {
            if (ief_bh_r0!=0) // only activates for analytic BH initial data (nonzero ief_bh_r0)
            {
              if (x/ex_r[l][0]<1) mask[i+j*shape[0]]=excised; 
            }
         }
      }
   }

}

//=============================================================================
//=============================================================================
void AdS5xS5_fill_bh_bboxes(real *bbox, int *num, int max_num)
{
   if (max_num<MAX_BHS) AMRD_stop("AdS5xS5_fill_bh_bboxes: error max_num too small\n","");

   *num=0;
}

//=============================================================================
// The following routine searches for AH's, manages excision, 
// and if t==0, determines past time level information using evolve
//
// NOTE: at this point, the time sequence is: n,nm1,np1 (unless t=0)
//=============================================================================
int AH_count[MAX_BHS],found_AH[MAX_BHS],freq0[MAX_BHS];
int pre_tstep_global_first=1,valid;

void AdS5xS5_pre_tstep(int L)
{
   char name[256];
   int AH[MAX_BHS];
   int AH_shape[1],got_an_AH,do_reinit_ex,do_repop;
   real M,J,c_equat,c_polar;
   real AH_bbox[2],AH_min_resid0;

   real ct;
   real new_rbuf;
   int n,i,Lf,Lc;

   ct=PAMR_get_time(L);

   Lf=PAMR_get_max_lev(PAMR_AMRH);
   Lc=PAMR_get_min_lev(PAMR_AMRH);  //if (PAMR_get_max_lev(PAMR_AMRH)>1) Lc=2; else Lc=1;

   return;
}


//=============================================================================
// The following routine prints diagnostic quantities
//
// NOTE: at this point, the time sequence is: n,nm1,np1 (unless t=0)
//=============================================================================
void AdS5xS5_post_tstep(int L)
{
   int itrace=1,valid;
   double lE4,gE4,liphi4,giphi4,cE4;
   char out1_name[256];
   static FILE *out1;
   static int local_first = 1;

   real ct;
   int n,i,Lf,Lc;

   ct=PAMR_get_time(L);

   Lf=PAMR_get_max_lev(PAMR_AMRH);  
   Lc=PAMR_get_min_lev(PAMR_AMRH);  //if (PAMR_get_max_lev(PAMR_AMRH)>1) Lc=2; else Lc=1;

   return;
}


#define ACTION_RELAX 1
#define ACTION_LOP 3
#define ACTION_RESIDUAL 2
//=============================================================================
// Returns some norm of the residual for the MG variables, *AND* 
// stores the point-wise residual in "f_res" for each MG variable "f" (for
// computing new RHS's)
//=============================================================================
real AdS5xS5_MG_residual(void)
{
   int action=ACTION_RESIDUAL;
   real norm;

   ldptr_mg();

   if (skip_constraints)
   {
      zero_f(zetab_res);
      return 0;
   }

   // solves for zetab conformal factor at t=0; residual
   mg_sup_(&action,zetab,zetab_rhs,zetab_lop,zetab_res,
           phi1,&AdS_L,mask_mg,phys_bdy,chr_mg,&AMRD_ex,x,y,&norm,&Nx,&Ny);

   return norm;
}


//=============================================================================
// Performs 1 relaxation sweep of the MG equations, and returns an estimate
// of the norm of the residual.
//=============================================================================
real AdS5xS5_MG_relax(void)
{
   int action=ACTION_RELAX;
   real norm;

   ldptr_mg();

   if (skip_constraints)
   {
      const_f(zetab,1);
      return 0;
   } 

   // solves for zetab conformal factor at t=0; residual
   mg_sup_(&action,zetab,zetab_rhs,zetab_lop,zetab_res,
           phi1,&AdS_L,mask_mg,phys_bdy,chr_mg,&AMRD_ex,x,y,&norm,&Nx,&Ny);

   return norm;
}

//=============================================================================
// Computes the differential operator L acting on the current grid,
// in a region specified by the grid function "mask". Stores the result
// in "f_lop" for each MG variable "f"
//=============================================================================
void AdS5xS5_L_op(void)
{
   int action=ACTION_LOP;
   real norm;

   ldptr_mg();

   if (skip_constraints)
   {
      zero_f(zetab_lop);
      return;
   }

   // solves for zetab conformal factor at t=0; residual
   mg_sup_(&action,zetab,zetab_rhs,zetab_lop,zetab_res,
           phi1,&AdS_L,mask_mg,phys_bdy,chr_mg,&AMRD_ex,x,y,&norm,&Nx,&Ny);

   return;
}

//=============================================================================
//=============================================================================
void AdS5xS5_scale_tre(void)
{
}

//=============================================================================
// post-regrid initialization of constant functions
//=============================================================================
void AdS5xS5_post_regrid(void)
{
   int i;

   if (!background || skip_constraints) return;

   ldptr();
}

//=============================================================================
//check-pointing
//=============================================================================
#define CP_DATA_SIZE 50000
void AdS5xS5_copy_block(char *r, char **q, int n, int dir, int *tot_size)
{
   char *r0;
   int n0;

   if (n==0) return;

   if ((*tot_size+n) > (CP_DATA_SIZE))
      AMRD_stop("AdS5xS5_copy_block: error ... CP_DATA_SIZE too small\n","");
   *tot_size+=n;

   n0=n;
   r0=r;

   if (dir==AMRD_CP_SAVE) while(n0--) *(*q)++=*r0++;
   else while(n0--) *r0++=*(*q)++;

   return;
}

void AdS5xS5_cp(int dir, char *data)
{
   int size=0;
   
   if (dir==AMRD_CP_SAVE)
   {
       cp_version=ADS5DP_CP_VERSION;
       AdS5xS5_copy_block((char *)&cp_version,&data,sizeof(int),dir,&size);
   }
}

//=============================================================================
int main(int argc, char **argv)
{
   amrd_set_app_user_cp_hook(AdS5xS5_cp,CP_DATA_SIZE);
   amrd_set_app_pre_tstep_hook(AdS5xS5_pre_tstep);
   amrd_set_elliptic_vars_t0_init(AdS5xS5_elliptic_vars_t0_init);
   amrd(argc,argv,&AdS5xS5_id,&AdS5xS5_var_pre_init,
        &AdS5xS5_var_post_init, &AdS5xS5_AMRH_var_clear,
        &AdS5xS5_free_data, &AdS5xS5_t0_cnst_data,
        &AdS5xS5_evo_residual, &AdS5xS5_MG_residual,
        &AdS5xS5_evolve, &AdS5xS5_MG_relax, &AdS5xS5_L_op, 
        &AdS5xS5_pre_io_calc, &AdS5xS5_scale_tre, 
        &AdS5xS5_post_regrid, &AdS5xS5_post_tstep,
        &AdS5xS5_fill_ex_mask, &AdS5xS5_fill_bh_bboxes);
}

