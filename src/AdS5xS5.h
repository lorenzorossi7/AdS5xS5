#ifndef _ADS5DP_H
#define _ADS5DP_H

/*============================================================================= */
/* global variables and prototypes for AdS5xS5                                    */
/*============================================================================= */

extern real AdS_L;

extern int pre_tstep_global_first;

/*============================================================================= */
/* id parameters                                                                */
/*============================================================================= */

#define MAX_BHS 4

extern real phi1_amp_1,phi1_r0_1,phi1_delta_1,phi1_x0_1[1],phi1_width_1[1];
extern real phi1_amp_2,phi1_r0_2,phi1_delta_2,phi1_x0_2[1],phi1_width_2[1];

extern real ief_bh_r0;

extern int background,skip_constraints;

extern int interptype, i_shift;

extern real rhoa,rhob;

#define ADS5DP_CP_VERSION 1
extern int cp_version;

extern real *phi1,*phi1_n,*phi1_np1,*phi1_nm1;
extern real *phi1_t,*phi1_t_n;
extern real *kg_res;

extern real *gb_tt,*gb_tt_n,*gb_tt_np1,*gb_tt_nm1;
extern real *gb_tx,*gb_tx_n,*gb_tx_np1,*gb_tx_nm1;
extern real *gb_xx,*gb_xx_n,*gb_xx_np1,*gb_xx_nm1;
extern real *psi,*psi_n,*psi_np1,*psi_nm1;
extern real *gb_ty,*gb_ty_n,*gb_ty_np1,*gb_ty_nm1;
extern real *gb_xy,*gb_xy_n,*gb_xy_np1,*gb_xy_nm1;
extern real *gb_yy,*gb_yy_n,*gb_yy_np1,*gb_yy_nm1;
extern real *omega,*omega_n,*omega_np1,*omega_nm1;
extern real *fb_t,*fb_t_n,*fb_t_np1,*fb_t_nm1;
extern real *fb_x,*fb_x_n,*fb_x_np1,*fb_x_nm1;
extern real *fb_y,*fb_y_n,*fb_y_np1,*fb_y_nm1;
extern real *gb_tt_t,*gb_tt_t_n;
extern real *gb_tx_t,*gb_tx_t_n;
extern real *gb_xx_t,*gb_xx_t_n;
extern real *psi_t,*psi_t_n;
extern real *gb_ty_t,*gb_ty_t_n;
extern real *gb_xy_t,*gb_xy_t_n;
extern real *gb_yy_t,*gb_yy_t_n;
extern real *omega_t,*omega_t_n;
extern real *fb_t_t,*fb_t_t_n;
extern real *fb_x_t,*fb_x_t_n;
extern real *fb_y_t,*fb_y_t_n;
extern real *gb_res;

extern real *Hb_t,*Hb_t_n,*Hb_t_np1,*Hb_t_nm1;
extern real *Hb_x,*Hb_x_n,*Hb_x_np1,*Hb_x_nm1;
extern real *Hb_t_t,*Hb_t_t_n;
extern real *Hb_x_t,*Hb_x_t_n;

extern real *zetab,*zetab_res,*zetab_lop,*zetab_rhs;

extern real *w1,*mg_w1;
extern real *w2,*mg_w2;
extern real *w3,*mg_w3;
extern real *w4,*mg_w4;

extern real *mask,*mask_mg,*chr,*chr_mg,AMRD_ex;
extern real *efe_all_ires,*kg_ires;
extern real *iresphi1,*iresall;

extern real *x;
extern int shape[1],ghost_width[2],Nx,phys_bdy[2],size,g_rank;
extern real base_bbox[2],bbox[2],dx,dt,dx_Lc;
extern int g_L;

extern int phi1_gfn,phi1_n_gfn,phi1_np1_gfn,phi1_nm1_gfn; 
extern int kg_res_gfn;

extern int gb_tt_gfn,gb_tt_n_gfn,gb_tt_np1_gfn,gb_tt_nm1_gfn;
extern int gb_tx_gfn,gb_tx_n_gfn,gb_tx_np1_gfn,gb_tx_nm1_gfn;
extern int gb_xx_gfn,gb_xx_n_gfn,gb_xx_np1_gfn,gb_xx_nm1_gfn;
extern int psi_gfn,psi_n_gfn,psi_np1_gfn,psi_nm1_gfn;
extern int gb_ty_gfn,gb_ty_n_gfn,gb_ty_np1_gfn,gb_ty_nm1_gfn;
extern int gb_xy_gfn,gb_xy_n_gfn,gb_xy_np1_gfn,gb_xy_nm1_gfn;
extern int gb_yy_gfn,gb_yy_n_gfn,gb_yy_np1_gfn,gb_yy_nm1_gfn;
extern int omega_gfn,omega_n_gfn,omega_np1_gfn,omega_nm1_gfn;
extern int fb_t_gfn,fb_t_n_gfn,fb_t_np1_gfn,fb_t_nm1_gfn;
extern int fb_x_gfn,fb_x_n_gfn,fb_x_np1_gfn,fb_x_nm1_gfn;
extern int fb_y_gfn,fb_y_n_gfn,fb_y_np1_gfn,fb_y_nm1_gfn;
extern int gb_tt_t_gfn,gb_tt_t_n_gfn;
extern int gb_tx_t_gfn,gb_tx_t_n_gfn;
extern int gb_xx_t_gfn,gb_xx_t_n_gfn;
extern int psi_t_gfn,psi_t_n_gfn;
extern int gb_ty_t_gfn,gb_ty_t_n_gfn;
extern int gb_xy_t_gfn,gb_xy_t_n_gfn;
extern int gb_yy_t_gfn,gb_yy_t_n_gfn;
extern int omega_t_gfn,omega_t_n_gfn;
extern int fb_t_t_gfn,fb_t_t_n_gfn;
extern int fb_x_t_gfn,fb_x_t_n_gfn;
extern int fb_y_t_gfn,fb_y_t_n_gfn;
extern int gb_res_gfn;

extern int Hb_t_gfn,Hb_t_n_gfn,Hb_t_np1_gfn,Hb_t_nm1_gfn;
extern int Hb_x_gfn,Hb_x_n_gfn,Hb_x_np1_gfn,Hb_x_nm1_gfn;
extern int Hb_t_t_gfn,Hb_t_t_n_gfn;
extern int Hb_x_t_gfn,Hb_x_t_n_gfn;

extern int zetab_gfn,zetab_res_gfn,zetab_lop_gfn,zetab_rhs_gfn;

extern int w1_gfn,mg_w1_gfn;
extern int w2_gfn,mg_w2_gfn;
extern int w3_gfn,mg_w3_gfn;
extern int w4_gfn,mg_w4_gfn;

extern int mask_gfn,mask_mg_gfn,chr_gfn,chr_mg_gfn;
extern int efe_all_ires_gfn,kg_ires_gfn;
extern int iresphi1_gfn,iresall_gfn;

#define MAX_N 2049

void set_gfns(void);
void ldptr_bbox(void);
void ldptr(void);
void ldptr_mg(void);
void const_f(real *f, real c);
void zero_f(real *f);
real norm_l2(real *f, real *cmask, real *chr);
void calc_gbu(void);
void calc_gbu_nm1(void);

/* prototypes for the various fortran functions we use: */
void g_evo_opt_(real *gb_res, real *kg_res, 
                real *gb_tt_np1, real *gb_tt_n, real *gb_tt_nm1,
                real *gb_tx_np1, real *gb_tx_n, real *gb_tx_nm1,
                real *gb_xx_np1, real *gb_xx_n, real *gb_xx_nm1,
                real *psi_np1, real *psi_n, real *psi_nm1,
                real *Hb_t_np1, real *Hb_t_n, real *Hb_t_nm1,
                real *Hb_x_np1, real *Hb_x_n, real *Hb_x_nm1,
                real *phi1_np1, real *phi1_n, real *phi1_nm1,
                real *AdS_L, real *x, real *dt,real *chr, real *ex, 
                int *phys_bdy, int *ghost_width, int *Nx,
                int *background, real *kappa_cd, real *rho_cd); 

void gauss2d_(real *f, 
              real *amp, real *r0, real *delta, real *xu0, real *ax,  
              real *amp2, real *r02, real *delta2, real *xu02, real *ax2,  
              real *AdS_L, real *x, int *Nx);

void AdS5xS5_fill_ex_mask(real *mask, int dim, int *shape, real *bbox, real excised);

void axi_reg_phi_(real *phi1, real *chr, real *ex, real *AdS_L, real *x, int *Nx);

void init_nm1_(real *gb_tt_np1, real *gb_tt_n, real *gb_tt_nm1, real *gb_tt_t_n,
               real *gb_tx_np1, real *gb_tx_n, real *gb_tx_nm1, real *gb_tx_t_n,
               real *gb_xx_np1, real *gb_xx_n, real *gb_xx_nm1, real *gb_xx_t_n,
               real *Hb_t_np1, real *Hb_t_n, real *Hb_t_nm1, real *Hb_t_t_n,
               real *Hb_x_np1, real *Hb_x_n, real *Hb_x_nm1, real *Hb_x_t_n,
               real *psi_np1, real *psi_n, real *psi_nm1, real *psi_t_n,
               real *phi1_np1, real *phi1_n, real *phi1_nm1, real *phi1_t_n,
               real *AdS_L, int *phys_bdy, real *x, real *dt, 
               real *chr, real *ex, int *Nx); 

void lin_zero_bnd_res_(real *f, int *phys_bdy, int *all, int *Nx);

void ires_(real *efe_all_ires,real *kg_ires, 
           real *gb_tt_np1, real *gb_tt_n, real *gb_tt_nm1,
           real *gb_tx_np1, real *gb_tx_n, real *gb_tx_nm1,
           real *gb_xx_np1, real *gb_xx_n, real *gb_xx_nm1,
           real *psi_np1, real *psi_n, real *psi_nm1,
           real *phi1_np1, real *phi1_n, real *phi1_nm1,
           real *x, real *dt, real *chr,
           real *AdS_L, real *ex, int *Nx, int *phys_bdy, int *ghost_width);

void mg_sup_(int *action, real *zetab, real *zetab_rhs, real *zetab_lop,real *zetab_res, 
             real *phi1, real *AdS_L, real *cmask, int *phys_bdy, real *chr, real *ex, real *x, real *norm, int *Nx);

void init_hb_(real *gb_tt_np1, real *gb_tt_n, real *gb_tt_nm1,
              real *gb_tx_np1, real *gb_tx_n, real *gb_tx_nm1, 
              real *gb_xx_np1, real *gb_xx_n, real *gb_xx_nm1, 
              real *psi_np1, real *psi_n, real *psi_nm1, 
              real *Hb_t_n, real *Hb_x_n, 
              real *AdS_L, int *phys_bdy, real *x, real *dt, real *chr, real *ex, int *Nx);

void hb_t_evo_(real *res, 
               real *gb_tt_np1, real *gb_tt_n, real *gb_tt_nm1,
               real *gb_tx_np1, real *gb_tx_n, real *gb_tx_nm1,
               real *gb_xx_np1, real *gb_xx_n, real *gb_xx_nm1,
               real *psi_np1, real *psi_n, real *psi_nm1,
               real *Hb_t_np1, real *Hb_t_n, real *Hb_t_nm1,
               real *Hb_x_np1, real *Hb_x_n, real *Hb_x_nm1,
               real *phi1_np1, real *phi1_n, real *phi1_nm1,
               real *AdS_L, real *x, real *dt, real *chr, real *ex,
               int *phys_bdy, int *ghost_width, int *Nx,
               real *Hb_t_0, real *Hb_x_0,
               int *gauge, real *t_n, real *rho1, real *rho2, real *xi1, real *xi2);

void hb_i_evo_(real *res, 
               real *gb_tt_np1, real *gb_tt_n, real *gb_tt_nm1,
               real *gb_tx_np1, real *gb_tx_n, real *gb_tx_nm1,
               real *gb_xx_np1, real *gb_xx_n, real *gb_xx_nm1,
               real *psi_np1, real *psi_n, real *psi_nm1,
               real *Hb_t_np1, real *Hb_t_n, real *Hb_t_nm1,
               real *Hb_x_np1, real *Hb_x_n, real *Hb_x_nm1,
               real *phi1_np1, real *phi1_n, real *phi1_nm1,
               real *AdS_L, real *x, real *dt, real *chr, real *ex,
               int *phys_bdy, int *ghost_width, int *Nx,
               real *Hb_t_0, real *Hb_x_0,
               int *gauge, real *t_n, real *rho1, real *rho2, real *xi1, real *xi2);

void init_ghb_(real *zetab, real *phi1, real *gb_tt, real *gb_tx, real *gb_xx, real *psi, real* rhoa, real* rhob,
               real *AdS_L, int *phys_bdy, real *chr, real *ex, real *x, int *Nx);

void init_schw(real *gb_tt, real *gb_tx, real *gb_xx, real *psi, real *ief_bh_r0, 
               real *AdS_L, int *phys_bdy, real *chr, real *ex, real *x, int *Nx);

#endif
