c----------------------------------------------------------------------
c in polar coordinates t,x for x in [0,1] 
c
c evolution routine for phi1 (not yet gb), computing
c the residual at the time just prior to updated
c
c L below is AdS_L (the length scale)
c----------------------------------------------------------------------
        subroutine g_evo_opt(gb_res,cl_res,
     &                       gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                       gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                       gb_ty_np1,gb_ty_n,gb_ty_nm1,
     &                       gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                       gb_xy_np1,gb_xy_n,gb_xy_nm1,
     &                       gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &                       psi_np1,psi_n,psi_nm1,
     &                       omega_np1,omega_n,omega_nm1,
     &                       Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &                       Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &                       Hb_y_np1,Hb_y_n,Hb_y_nm1,
     &                       phi1_np1,phi1_n,phi1_nm1,
     &                       L,x,y,dt,chr,ex,
     &                       phys_bdy,ghost_width,Nx,Ny,
     &                       background,kappa_cd,rho_cd)
        implicit none
        integer Nx,Ny
        integer phys_bdy(4),ghost_width(4)
        integer background
        real*8 kappa_cd,rho_cd
        real*8 gb_res(Nx,Ny),cl_res(Nx,Ny)
        real*8 gb_tt_np1(Nx,Ny),gb_tt_n(Nx,Ny),gb_tt_nm1(Nx,Ny)
        real*8 gb_tx_np1(Nx,Ny),gb_tx_n(Nx,Ny),gb_tx_nm1(Nx,Ny)
        real*8 gb_ty_np1(Nx,Ny),gb_ty_n(Nx,Ny),gb_ty_nm1(Nx,Ny)
        real*8 gb_xx_np1(Nx,Ny),gb_xx_n(Nx,Ny),gb_xx_nm1(Nx,Ny)
        real*8 gb_xy_np1(Nx,Ny),gb_xy_n(Nx,Ny),gb_xy_nm1(Nx,Ny)
        real*8 gb_yy_np1(Nx,Ny),gb_yy_n(Nx,Ny),gb_yy_nm1(Nx,Ny)
        real*8 psi_np1(Nx,Ny),psi_n(Nx,Ny),psi_nm1(Nx,Ny)
        real*8 omega_np1(Nx,Ny),omega_n(Nx,Ny),omega_nm1(Nx,Ny)
        real*8 Hb_t_np1(Nx,Ny),Hb_t_n(Nx,Ny),Hb_t_nm1(Nx,Ny)
        real*8 Hb_x_np1(Nx,Ny),Hb_x_n(Nx,Ny),Hb_x_nm1(Nx,Ny)
        real*8 Hb_y_np1(Nx,Ny),Hb_y_n(Nx,Ny),Hb_y_nm1(Nx,Ny)
        real*8 phi1_np1(Nx,Ny),phi1_n(Nx,Ny),phi1_nm1(Nx,Ny)
        real*8 L
        real*8 x(Nx),y(Ny),dt,chr(Nx,Ny),ex

        integer a,b,c,d,e

        integer rb,i,j
        real*8 PI
        parameter (PI=3.141592653589793d0)
        real*8 dx,dy

        integer max_ghost_width

        real*8 boxx_u(3),boxx_l(3)

        real*8 x0,y0

        real*8 phi1_res,phi1_J
        real*8 phi10_res,phi10_J

        integer is,ie,js,je,ks,ke,is_a_nan

        logical ltrace,is_nan,dump,first_nan
        parameter (ltrace=.false.)
        data first_nan/.true./

        !--------------------------------------------------------------
        ! using g0_ab = gads_ab + h0_ab
        !       g0^ab = gads^ab + h0^ab
        !       H0_a  = Hads_a + A_a
        ! where h0_ab, h0^ab are not inverses of each other
        !
        ! g0_ll(a,b)           = g0_ab           !use for diagnostics  !
        ! g0_uu(a,b)           = g0^ab           !eg: compare g0_ll    !
        ! g0_ll_x(a,b,c)       = g0_ab_,c        !    with gads_ll+h_ll!  
        ! g0_ll_xx(a,b,c,d)    = g0_ab_,cd
        ! g0_uu_x(a,b,c)       = g0^ab_,c
        !                      = -g0^ad g0^be g0_de_,c 
        !
        ! gads_ll(a,b)         = gads_ab
        ! gads_uu(a,b)         = gads^ab
        ! gads_ll_x(a,b,c)     = gads_ab_,c
        ! gads_ll_xx(a,b,c,d)  = gads_ab_,cd
        ! gads_uu_x(a,b,c)     = gads^ab_,c
        !                      = -gads^ad gads^be gads_de_,c 
        !
        ! h0_ll(a,b)           = (1-x0**2)*gb_ab   
        ! h0_uu(a,b)           = inverse(gads+(1-x0**2)*gb)^ab - gads^ab
        ! h0_ll_x(a,b,c)       = ((1-x0**2)*gb_ab)_,c
        ! h0_ll_xx(a,b,c,d)    = ((1-x0**2)*gb_ab)_,cd
        ! h0_uu_x(a,b,c)       = g^ab_,c - gads^ab_,c
        !
        ! cuuuu(a,b,c,d) = gads_uu(a,b)*gads_uu(c,d)
        !                   +h0_uu(a,b)*h0_uu(c,d) 
        !                   +gads_uu(a,b)*h0_uu(c,d) 
        !                   +h0_uu(a,b)*gads_uu(c,d) 
        !
        ! dlll(a,b,c)    = g0_ll_x(b,c,a)-g0_ll_x(a,b,c)+g0_ll_x(c,a,b)
        !
        ! given z_t=(1-x0^2)^3, z_x=(1-x0^2)^2, z_y=(1-x0^2)^3
        ! A_l(a)     = (1-x0**2)**2*Hb_a   
        ! A_l_x(a,b) = ((1-x0**2)**2*Hb_a)_,b
        ! 
        ! phi10_x(a)  = phi1_,a     
        ! phi4rx(a) = phi4r_,a      
        ! phi4ix(a) = phi4i_,a
        !
        ! grad_phi1_sq = g^cd*phi1_,c*phi1_,d
        ! grad_phi4_sq = g^cd*phi4r_,c*phi4r_,d + g^cd*phi4i_,c*phi4i_,d
        !
        ! set_ab = 2*phi1_,a*phi1_,b - g_ab*grad_phi1_sq
        !          + 2*phi4r_,a*phi4r_,b + 2*phi4i_,a*phi4i_,b
        !          - g_ab*grad_phi4_sq 
        !          - g_ab*( mass4*(phi4r**2+phi4i**2)
        !                  -lambda3*(phi4r**2+phi4i**2)**(1.5d0)
        !                  +lambda4*(phi4r**2+phi4i**2)**2 )
        !
        ! efe(a,b) = residual ... hardcoded expressions (see below)
        !
        ! t,x,y=1,2,3
        ! 
        ! NOTE: g0_ll_xx,gads_ll_xx,h0_ll_xx,efe,efe_J
        !       do *NOT* symmetric components filled in
        !
        !--------------------------------------------------------------
        real*8 efe(3,3),efe_J(3,3)
        real*8 afe,afe_J
        real*8 bfe,bfe_J
        real*8 ffe(3),ffe_J(3)
        real*8 term1(3,3),term2(3,3),term3(3,3),term4(3,3)
        real*8 term5(3,3),term6(3,3),term7(3,3),term8(3,3)
        real*8 term9(3,3)
 
        real*8 ndotc,n_l(3),n_u(3),c_l(3)
        real*8 cd_ll(3,3),cd_J_ll(3,3)

        real*8 grad_phi1_sq
        
        real*8 grad_phi4_sq

        real*8 g0u_tt_ads0,g0u_xx_ads0,g0u_xy_ads0,g0u_yy_ads0

        real*8 H0_t_ads0,H0_x_ads0,H0_y_ads0

        real*8 dgb_J,ddgb_J,ddgb_J_tx,ddgb_J_ty
        real*8 dphi1_J,ddphi1_J,ddphi1_J_tx
        real*8 dc_J

        real*8 dimA,dimB

        !--------------------------------------------------------------
        ! variables for tensor manipulations 
        !(indices are t,x,y,theta,phi)
        !--------------------------------------------------------------
        real*8 g0_ll(3,3),g0_uu(3,3)
        real*8 g0_ll_x(3,3,3),g0_uu_x(3,3,3),g0_ll_xx(3,3,3,3)
        real*8 gads_ll(3,3),gads_uu(3,3)
        real*8 gads_ll_x(3,3,3),gads_uu_x(3,3,3),gads_ll_xx(3,3,3,3)
        real*8 gA,gA_x(3),gA_xx(3,3)
        real*8 gB,gB_x(3),gB_xx(3,3)
        real*8 gAads,gAads_x(3),gAads_xx(3,3)
        real*8 gBads,gBads_x(3),gBads_xx(3,3)
        real*8 h0_ll(3,3),h0_uu(3,3)
        real*8 h0_ll_x(3,3,3),h0_uu_x(3,3,3),h0_ll_xx(3,3,3,3)
        real*8 hA,hAu,hA_x(3),hA_xx(3,3)
        real*8 hB,hBu,hB_x(3),hB_xx(3,3)
        real*8 gamma_ull(3,3,3),gamma_ull_x(3,3,3,3)
        real*8 riemann_ulll(3,3,3,3)
        real*8 ricci_ll(3,3),ricci_lu(3,3),ricci
        real*8 H0_l(3),H0_l_x(3,3)
        real*8 Hads_l(3),Hads_l_x(3,3),A_l(3),A_l_x(3,3)
        real*8 phi10_x(3),phi10_xx(3,3)

        real*8 f1_l(3),f1_l_x(3,3),f2_ll(3,3),f2_ll_x(3,3,3)
        real*8 f1ads_l(3),f1ads_l_x(3,3),f2ads_ll(3,3),f2ads_ll_x(3,3,3)
        real*8 h1_l(3),h1_l_x(3,3),h2_ll(3,3),h2_ll_x(3,3,3)
        real*8 sqrtdetg,sqrtdetg_x(3)
        real*8 sqrtdetgads,sqrtdetgads_x(3)
        real*8 sqrtdeth,sqrtdeth_x(3)

        real*8 gammagg(3,3,3),gammahh(3,3,3)
        real*8 gammagh(3,3,3),gammahg(3,3,3)
        real*8 cuuuu(3,3,3,3),dlll(3,3,3)

        !--------------------------------------------------------------
        ! initialize fixed-size variables 
        !--------------------------------------------------------------
        data a,b,c,d,e/0,0,0,0,0/
        data rb,i/0,0/
        data max_ghost_width/0/
        data is,ie,js,je,ks,ke,is_a_nan/0,0,0,0,0,0,0/

        data g0u_tt_ads0,g0u_xx_ads0/0.0,0.0/
        data g0u_xy_ads0,g0u_yy_ads0/0.0,0.0/

        data H0_t_ads0,H0_x_ads0,H0_y_ads0/0.0,0.0,0.0/

        data dgb_J,ddgb_J,ddgb_J_tx,ddgb_J_ty/0.0,0.0,0.0,0.0/
        data dphi1_J,ddphi1_J/0.0,0.0/
        data dc_J/0.0/

        data term1,term2/9*0.0,9*0.0/
        data term3,term4/9*0.0,9*0.0/
        data term5,term6/9*0.0,9*0.0/
        data term7,term8/9*0.0,9*0.0/
        data term9/9*0.0/

        data efe,efe_J/9*0.0,9*0.0/
        data afe,afe_J/0.0,0.0/
        data bfe,bfe_J/0.0,0.0/
        data ffe,ffe_J/3*0.0,3*0.0/
        data cd_ll,cd_J_ll/9*0.0,9*0.0/

        data x0/0.0/
        data y0/0.0/

        data dx/0.0/
        data dy/0.0/

        data boxx_u,boxx_l/3*0.0,3*0.0/

        data phi1_res,phi1_J/0.0,0.0/
        data phi10_res,phi10_J/0.0,0.0/

        data n_l,n_u,c_l/3*0.0,3*0.0,3*0.0/

        data g0_ll,g0_uu/9*0.0,9*0.0/
        data gads_ll,gads_uu/9*0.0,9*0.0/
        data h0_ll,h0_uu/9*0.0,9*0.0/
        data gamma_ull/27*0.0/
        data gamma_ull_x/81*0.0/

        data g0_ll_x,g0_uu_x/27*0.0,27*0.0/
        data gads_ll_x,gads_uu_x/27*0.0,27*0.0/
        data h0_ll_x,h0_uu_x/27*0.0,27*0.0/

        data g0_ll_xx/81*0.0/
        data gads_ll_xx/81*0.0/
        data h0_ll_xx/81*0.0/

        data ricci/0.0/
        data ricci_ll,ricci_lu/9*0.0,9*0.0/
        data riemann_ulll/81*0.0/

        data f1_l,f1_l_x/3*0.0,9*0.0/
        data f2_ll,f2_ll_x/9*0.0,27*0.0/
        data f1ads_l,f1ads_l_x/3*0.0,9*0.0/
        data f2ads_ll,f2ads_ll_x/9*0.0,27*0.0/
        data h1_l,h1_l_x/3*0.0,9*0.0/
        data h2_ll,h2_ll_x/9*0.0,27*0.0/
        data sqrtdetg,sqrtdetg_x/0.0,3*0.0/
        data sqrtdetgads,sqrtdetgads_x/0.0,3*0.0/
        data sqrtdeth,sqrtdeth_x/0.0,3*0.0/

        data H0_l,Hads_l,A_l/3*0.0,3*0.0,3*0.0/
        data H0_l_x,Hads_l_x,A_l_x/9*0.0,9*0.0,9*0.0/

        data phi10_x/3*0.0/
        data phi10_xx/9*0.0/

        data gA,gB/0.0,0.0/
        data gA_x,gB_x/3*0.0,3*0.0/
        data gA_xx,gB_xx/9*0.0,9*0.0/
        data gAads,gBads/0.0,0.0/
        data gAads_x,gBads_x/3*0.0,3*0.0/
        data gAads_xx,gBads_xx/9*0.0,9*0.0/
        data hA,hB/0.0,0.0/
        data hAu,hBu/0.0,0.0/
        data hA_x,hB_x/3*0.0,3*0.0/
        data hA_xx,hB_xx/9*0.0,9*0.0/

        data gammagg,gammahh/27*0.0,27*0.0/
        data gammagh,gammahg/27*0.0,27*0.0/
        data dlll/27*0.0/
        data cuuuu/81*0.0/

        !--------------------------------------------------------------
        if (ltrace) write(*,*) 'gb_psi_evo ... N=',Nx

        dx=x(2)-x(1)
        dy=y(2)-y(1)

        ! set dimensions of S3 and S4 subspaces
        dimA=3d0 
        dimB=4d0

        ! set index bounds for main loop
        is=2
        ie=Nx-1
        js=2
        je=Ny-1

        ! adjust index bounds to compensate for ghost_width
        if (ghost_width(1).gt.0) is=is+ghost_width(1)-1
        if (ghost_width(2).gt.0) ie=ie-(ghost_width(2)-1)
        if (ghost_width(3).gt.0) js=js+ghost_width(3)-1
        if (ghost_width(4).gt.0) je=je-(ghost_width(4)-1)

        !(MAIN LOOP) red-black loop through spacetime pts x(i)
        do rb=0,1

          do j=js,je
            do i=is+mod(j+rb,2),ie,2
              x0=x(i)
              y0=y(j)

              dump=.false.

              if (ltrace) write(*,*) 'i,j:',i,j

              !(REGION) interior points; evolve 
              if (chr(i,j).ne.ex) then

                ! computes tensors at point i,j
                call tensor_init(
     &                  gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                  gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                  gb_ty_np1,gb_ty_n,gb_ty_nm1,
     &                  gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                  gb_xy_np1,gb_xy_n,gb_xy_nm1,
     &                  gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &                  psi_np1,psi_n,psi_nm1,
     &                  omega_np1,omega_n,omega_nm1,
     &                  Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &                  Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &                  Hb_y_np1,Hb_y_n,Hb_y_nm1,
     &                  phi1_np1,phi1_n,phi1_nm1,
     &                  g0_ll,g0_uu,g0_ll_x,g0_uu_x,g0_ll_xx,
     &                  gads_ll,gads_uu,gads_ll_x,gads_uu_x,gads_ll_xx,
     &                  gA,gB,gA_x,gB_x,gA_xx,gB_xx,
     &                  gAads,gBads,gAads_x,gBads_x,gAads_xx,gBads_xx,
     &                  sqrtdetg,sqrtdetg_x,
     &                  sqrtdetgads,sqrtdetgads_x,
     &                  h0_ll,h0_uu,h0_ll_x,h0_uu_x,h0_ll_xx,
     &                  hA,hB,hAu,hBu,hA_x,hB_x,hA_xx,hB_xx,
     &                  sqrtdeth,sqrtdeth_x,
     &                  H0_l,H0_l_x,Hads_l,Hads_l_x,A_l,A_l_x,
     &                  gamma_ull,gamma_ull_x,
     &                  riemann_ulll,ricci_ll,ricci_lu,ricci,
     &                  f1_l,f1_l_x,f2_ll,f2_ll_x,
     &                  f1ads_l,f1ads_l_x,f2ads_ll,f2ads_ll_x,
     &                  h1_l,h1_l_x,h2_ll,h2_ll_x,
     &                  phi10_x,phi10_xx,
     &                  gammagg,gammahh,gammagh,gammahg,dlll,cuuuu,
     &                  x,y,dt,chr,L,ex,Nx,Ny,i,j)

                do c=1,3
                  boxx_u(c)=-( gamma_ull(c,1,1)*g0_uu(1,1)+
     &                         gamma_ull(c,2,2)*g0_uu(2,2)+
     &                         gamma_ull(c,3,3)*g0_uu(3,3)+
     &                      2*(gamma_ull(c,1,2)*g0_uu(1,2)+
     &                         gamma_ull(c,1,3)*g0_uu(1,3)+
     &                         gamma_ull(c,2,3)*g0_uu(2,3)) )
                end do

                do a=1,3
                  boxx_l(a)=boxx_u(1)*g0_ll(a,1)+
     &                      boxx_u(2)*g0_ll(a,2)+
     &                      boxx_u(3)*g0_ll(a,3)
                end do
    
                do a=1,3
                  c_l(a)=(Hads_l(a)+A_l(a))-boxx_l(a)
                end do
 
                n_l(1)=-1/sqrt(-g0_uu(1,1))
                do a=1,3
                  n_u(a)=n_l(1)*g0_uu(a,1)+
     &                   n_l(2)*g0_uu(a,2)+
     &                   n_l(3)*g0_uu(a,3)
                end do

                ndotc  =n_l(1)*c_l(1)*g0_uu(1,1)+
     &                  n_l(2)*c_l(2)*g0_uu(2,2)+
     &                  n_l(3)*c_l(3)*g0_uu(3,3)+
     &               2*(n_l(1)*c_l(2)*g0_uu(1,2)+
     &                  n_l(1)*c_l(3)*g0_uu(1,3)+
     &                  n_l(2)*c_l(3)*g0_uu(2,3))

                grad_phi1_sq=phi10_x(1)*phi10_x(1)*g0_uu(1,1)+
     &                       phi10_x(2)*phi10_x(2)*g0_uu(2,2)+
     &                       phi10_x(3)*phi10_x(3)*g0_uu(3,3)+
     &                    2*(phi10_x(1)*phi10_x(2)*g0_uu(1,2)+
     &                       phi10_x(1)*phi10_x(3)*g0_uu(1,3)+
     &                       phi10_x(2)*phi10_x(3)*g0_uu(2,3))

                !---------------------------------------------------------------- 
                ! Included pure AdS terms in the EFEs,
                ! efe_ab =   term1_ab + term2_ab + term3_ab + term4_ab 
                !          + term5_ab + term6_ab + term7_ab + term8_ab 
                !          + term9_ab
                ! for
                ! 
                ! term1_ab = -1/2 g^cd g_ab,cd  
                ! term2_ab = -1/2 g^cd,a g_bc,d
                ! term3_ab = -1/2 g^cd,b g_ac,d
                ! term4_ab = -1/2 H_a,b
                ! term5_ab = -1/2 H_b,a
                ! term6_ab = H_c G^c_ab
                ! term7_ab = -G^c_db G^d_ca
                ! term8_ab = -dimA*(gA_,a gA_,b)/(4 gA^2)
                !            -dimB*(gB_,a gB_,b)/(4 gB^2)
                ! term9_ab = -(2*phi1,a phi1,b - g_ab g^cd phi1,c phi1,d)/(4 gB^4)
                !
                ! where G   = guu(g_ll_x-g_ll_x+g_ll_x)
                !
                ! and the indices are
                ! a,b = (t,x,y)
                !
                !----------------------------------------------------------------
                do a=1,3
                  do b=a,3

                    term1(a,b)=-0.5d0*(                             
     &                            g0_uu(1,1)*g0_ll_xx(a,b,1,1)+
     &                            g0_uu(2,2)*g0_ll_xx(a,b,2,2)+
     &                            g0_uu(3,3)*g0_ll_xx(a,b,3,3)+
     &                         2*(g0_uu(1,2)*g0_ll_xx(a,b,1,2)+
     &                            g0_uu(1,3)*g0_ll_xx(a,b,1,3)+
     &                            g0_uu(2,3)*g0_ll_xx(a,b,2,3))
     &                                )
     &
                    term2(a,b)=-0.5d0*(                             
     &                            g0_uu_x(1,1,a)* g0_ll_x(b,1,1) +
     &                            g0_uu_x(1,2,a)*(g0_ll_x(b,1,2) +
     &                                            g0_ll_x(b,2,1))+
     &                            g0_uu_x(1,3,a)*(g0_ll_x(b,1,3) +
     &                                            g0_ll_x(b,3,1))+
     &                            g0_uu_x(2,2,a)* g0_ll_x(b,2,2) +
     &                            g0_uu_x(2,3,a)*(g0_ll_x(b,2,3) +
     &                                            g0_ll_x(b,3,2))+
     &                            g0_uu_x(3,3,a)* g0_ll_x(b,3,3) 
     &                              )
     &
                    term3(a,b)=-0.5d0*(                            
     &                            g0_uu_x(1,1,b)* g0_ll_x(a,1,1) +
     &                            g0_uu_x(1,2,b)*(g0_ll_x(a,1,2) +
     &                                            g0_ll_x(a,2,1))+
     &                            g0_uu_x(1,3,b)*(g0_ll_x(a,1,3) +
     &                                            g0_ll_x(a,3,1))+
     &                            g0_uu_x(2,2,b)* g0_ll_x(a,2,2) +
     &                            g0_uu_x(2,3,b)*(g0_ll_x(a,2,3) +
     &                                            g0_ll_x(a,3,2))+
     &                            g0_uu_x(3,3,b)* g0_ll_x(a,3,3) 
     &                              )
     &
                    term4(a,b)=-0.5d0*H0_l_x(a,b)
     &
                    term5(a,b)=-0.5d0*H0_l_x(b,a)
     &
                    term6(a,b)=     (           
     &                            H0_l(1)*gamma_ull(1,a,b)+
     &                            H0_l(2)*gamma_ull(2,a,b)+
     &                            H0_l(3)*gamma_ull(3,a,b)
     &                              ) 
     &                             
                    term7(a,b)=    -(
     &                            gamma_ull(1,1,b)*gamma_ull(1,1,a)+
     &                            gamma_ull(1,2,b)*gamma_ull(2,1,a)+
     &                            gamma_ull(1,3,b)*gamma_ull(3,1,a)+
     &                            gamma_ull(2,1,b)*gamma_ull(1,2,a)+
     &                            gamma_ull(2,2,b)*gamma_ull(2,2,a)+
     &                            gamma_ull(2,3,b)*gamma_ull(3,2,a)+
     &                            gamma_ull(3,1,b)*gamma_ull(1,3,a)+
     &                            gamma_ull(3,2,b)*gamma_ull(2,3,a)+
     &                            gamma_ull(3,3,b)*gamma_ull(3,3,a)
     &                              )
     &
                    term8(a,b)=-dimA*gA_x(a)*gA_x(b)/(4*gA**2)
     &                         -dimB*gB_x(a)*gB_x(b)/(4*gB**2)
     &
                    term9(a,b)=-(2*phi10_x(a)*phi10_x(b)
     &                          -g0_ll(a,b)*grad_phi1_sq
     &                          )/4/gB**4
     &
                    efe(a,b)=term1(a,b)+term2(a,b)+term3(a,b)+term4(a,b)
     &                      +term5(a,b)+term6(a,b)+term7(a,b)+term8(a,b)
     &                      +term9(a,b)

                  end do
                end do

                !--------------------------------------------------------------------------
                ! afe      = -g^ab gA_,ab/2 - H_a gA_,b/2 + (gA_,a gA_,b)/(2 gA)
                !            +(dimA-1)
                !            +gA*(g^ab phi1,a phi1,b)/4
                !
                ! where G   = guu(g_ll_x-g_ll_x+g_ll_x)
                !--------------------------------------------------------------------------

                afe=-g0_uu(1,1)*gA_xx(1,1)/2
     &              -g0_uu(1,2)*gA_xx(1,2)/2
     &              -g0_uu(1,3)*gA_xx(1,3)/2
     &              -g0_uu(2,1)*gA_xx(2,1)/2
     &              -g0_uu(2,2)*gA_xx(2,2)/2
     &              -g0_uu(2,3)*gA_xx(2,3)/2
     &              -g0_uu(3,1)*gA_xx(3,1)/2
     &              -g0_uu(3,2)*gA_xx(3,2)/2
     &              -g0_uu(3,3)*gA_xx(3,3)/2
     &              -g0_uu(1,1)*H0_l(1)*gA_x(1)/2
     &              -g0_uu(1,2)*H0_l(1)*gA_x(2)/2
     &              -g0_uu(1,3)*H0_l(1)*gA_x(3)/2
     &              -g0_uu(2,1)*H0_l(2)*gA_x(1)/2
     &              -g0_uu(2,2)*H0_l(2)*gA_x(2)/2
     &              -g0_uu(2,3)*H0_l(2)*gA_x(3)/2
     &              -g0_uu(3,1)*H0_l(3)*gA_x(1)/2
     &              -g0_uu(3,2)*H0_l(3)*gA_x(2)/2
     &              -g0_uu(3,3)*H0_l(3)*gA_x(3)/2
     &              +g0_uu(1,1)*(gA_x(1)*gA_x(1))/(2*gA)
     &              +g0_uu(1,2)*(gA_x(1)*gA_x(2))/(2*gA)
     &              +g0_uu(1,3)*(gA_x(1)*gA_x(3))/(2*gA)
     &              +g0_uu(2,1)*(gA_x(2)*gA_x(1))/(2*gA)
     &              +g0_uu(2,2)*(gA_x(2)*gA_x(2))/(2*gA)
     &              +g0_uu(2,3)*(gA_x(2)*gA_x(3))/(2*gA)
     &              +g0_uu(3,1)*(gA_x(3)*gA_x(1))/(2*gA)
     &              +g0_uu(3,2)*(gA_x(3)*gA_x(2))/(2*gA)
     &              +g0_uu(3,3)*(gA_x(3)*gA_x(3))/(2*gA)
     &
     &              +(dimA-1d0)
     &
     &              +gA*grad_phi1_sq/4/gB**4

                !--------------------------------------------------------------------------
                ! bfe      = -g^ab gB_,ab/2 - H_a gB_,b/2 + (gB_,a gB_,b)/(2 gB)
                !            +(dimB-1)
                !            -gB*(g^ab phi1,a phi1,b)/4
                !
                ! where G   = guu(g_ll_x-g_ll_x+g_ll_x)
                !--------------------------------------------------------------------------

                bfe=-g0_uu(1,1)*gB_xx(1,1)/2
     &              -g0_uu(1,2)*gB_xx(1,2)/2
     &              -g0_uu(1,3)*gB_xx(1,3)/2
     &              -g0_uu(2,1)*gB_xx(2,1)/2
     &              -g0_uu(2,2)*gB_xx(2,2)/2
     &              -g0_uu(2,3)*gB_xx(2,3)/2
     &              -g0_uu(3,1)*gB_xx(3,1)/2
     &              -g0_uu(3,2)*gB_xx(3,2)/2
     &              -g0_uu(3,3)*gB_xx(3,3)/2
     &              -g0_uu(1,1)*H0_l(1)*gB_x(1)/2
     &              -g0_uu(1,2)*H0_l(1)*gB_x(2)/2
     &              -g0_uu(1,3)*H0_l(1)*gB_x(3)/2
     &              -g0_uu(2,1)*H0_l(2)*gB_x(1)/2
     &              -g0_uu(2,2)*H0_l(2)*gB_x(2)/2
     &              -g0_uu(2,3)*H0_l(2)*gB_x(3)/2
     &              -g0_uu(3,1)*H0_l(3)*gB_x(1)/2
     &              -g0_uu(3,2)*H0_l(3)*gB_x(2)/2
     &              -g0_uu(3,3)*H0_l(3)*gB_x(3)/2
     &              +g0_uu(1,1)*(gB_x(1)*gB_x(1))/(2*gB)
     &              +g0_uu(1,2)*(gB_x(1)*gB_x(2))/(2*gB)
     &              +g0_uu(1,3)*(gB_x(1)*gB_x(3))/(2*gB)
     &              +g0_uu(2,1)*(gB_x(2)*gB_x(1))/(2*gB)
     &              +g0_uu(2,2)*(gB_x(2)*gB_x(2))/(2*gB)
     &              +g0_uu(2,3)*(gB_x(2)*gB_x(3))/(2*gB)
     &              +g0_uu(3,1)*(gB_x(3)*gB_x(1))/(2*gB)
     &              +g0_uu(3,2)*(gB_x(3)*gB_x(2))/(2*gB)
     &              +g0_uu(3,3)*(gB_x(3)*gB_x(3))/(2*gB)
     &
     &              +(dimB-1d0)
     &
     &              -grad_phi1_sq/4/gB**3

                !---------------------------------------------------------------- 
                ! ffe_t = -sqrtdetg g^tt f_t,t 
                !         -sqrtdetg g^tx (f_t,x + 2 f_t gB_,x/gB - 2 f_x gB_,t/gB)
                !         -sqrtdetg g^ty (f_t,y + 2 f_t gB_,y/gB - 2 f_y gB_,t/gB)
                !         -sqrtdetg g^ta_,t f_a 
                !         -sqrtdetg_,t g^ta f_a 
                !         -f_ty,x + f_tx,y 
                !         +1.5d0 (f_xy gA_,t/gA - f_ty gA_,x/gA + f_tx gA_,y/gA)
                ! ffe_x = -f_x,t + f_t,x + 2 f_t gB_,x/gB - 2 f_x gB_,t/gB 
                ! ffe_y = -f_y,t + f_t,y + 2 f_t gB_,y/gB - 2 f_y gB_,t/gB 
                !
                !---------------------------------------------------------------- 
                ffe(1)=-sqrtdetg*g0_uu(1,1)*f1_l_x(1,1)
     &                 -sqrtdetg*g0_uu(1,2)*(f1_l_x(1,2)
     &                                    +2*f1_l(1)*gB_x(2)/gB
     &                                    -2*f1_l(2)*gB_x(1)/gB)
     &                 -sqrtdetg*g0_uu(1,3)*(f1_l_x(1,3)
     &                                    +2*f1_l(1)*gB_x(3)/gB
     &                                    -2*f1_l(3)*gB_x(1)/gB)
     &                 -sqrtdetg           *(g0_uu_x(1,1,1)*f1_l(1) 
     &                                      +g0_uu_x(1,2,1)*f1_l(2)
     &                                      +g0_uu_x(1,3,1)*f1_l(3))
     &                 -sqrtdetg_x(1)      *(g0_uu(1,1)*f1_l(1)
     &                                      +g0_uu(1,2)*f1_l(2)
     &                                      +g0_uu(1,3)*f1_l(3))
     &                 -f2_ll_x(1,3,2)+f2_ll_x(1,2,3)
     &                 +1.5d0              *(f2_ll(2,3)*gA_x(1)/gA
     &                                      -f2_ll(1,3)*gA_x(2)/gA 
     &                                      +f2_ll(1,2)*gA_x(3)/gA)

                ffe(2)=-f1_l_x(2,1)+f1_l_x(1,2)
     &                 +2*f1_l(1)*gB_x(2)/gB-2*f1_l(2)*gB_x(1)/gB 

                ffe(3)=-f1_l_x(3,1)+f1_l_x(1,3)
     &                 +2*f1_l(1)*gB_x(3)/gB-2*f1_l(3)*gB_x(1)/gB 

                !--------------------------------------------------------------------------
                ! phi1_res = g^ab phi1,ab - g^ab gamma^c_ab phi1,c
                !           +1.5d0/gA*g^ab gA,a phi1,b - 2.0d0/gB*g^ab gB,a phi1,b
                !--------------------------------------------------------------------------
                phi1_res= phi10_xx(1,1)*g0_uu(1,1)+
     &                    phi10_xx(2,2)*g0_uu(2,2)+
     &                    phi10_xx(3,3)*g0_uu(3,3)+
     &                 2*(phi10_xx(1,2)*g0_uu(1,2)+
     &                    phi10_xx(1,3)*g0_uu(1,3)+
     &                    phi10_xx(2,3)*g0_uu(2,3))
     &                  +
     &                     g0_uu(1,1)*H0_l(1)*phi10_x(1)+
     &                     g0_uu(1,2)*H0_l(1)*phi10_x(2)+
     &                     g0_uu(1,3)*H0_l(1)*phi10_x(3)+
     &                     g0_uu(2,1)*H0_l(2)*phi10_x(1)+
     &                     g0_uu(2,2)*H0_l(2)*phi10_x(2)+
     &                     g0_uu(2,3)*H0_l(2)*phi10_x(3)+
     &                     g0_uu(3,1)*H0_l(3)*phi10_x(1)+
     &                     g0_uu(3,2)*H0_l(3)*phi10_x(2)+
     &                     g0_uu(3,3)*H0_l(3)*phi10_x(3)
     &                  +
     &                     (3d0-dimA)/2/gA*(
     &                      gA_x(1)*phi10_x(1)*g0_uu(1,1)+
     &                      gA_x(2)*phi10_x(2)*g0_uu(2,2)+
     &                      gA_x(3)*phi10_x(3)*g0_uu(3,3)+
     &                      gA_x(1)*phi10_x(2)*g0_uu(1,2)+
     &                      gA_x(1)*phi10_x(3)*g0_uu(1,3)+
     &                      gA_x(2)*phi10_x(3)*g0_uu(2,3)+
     &                      gA_x(2)*phi10_x(1)*g0_uu(1,2)+
     &                      gA_x(3)*phi10_x(1)*g0_uu(1,3)+
     &                      gA_x(3)*phi10_x(2)*g0_uu(2,3)
     &                                               )
     &                  -
     &                     (4d0+dimB)/2/gB*(
     &                      gB_x(1)*phi10_x(1)*g0_uu(1,1)+
     &                      gB_x(2)*phi10_x(2)*g0_uu(2,2)+
     &                      gB_x(3)*phi10_x(3)*g0_uu(3,3)+
     &                      gB_x(1)*phi10_x(2)*g0_uu(1,2)+
     &                      gB_x(1)*phi10_x(3)*g0_uu(1,3)+
     &                      gB_x(2)*phi10_x(3)*g0_uu(2,3)+
     &                      gB_x(2)*phi10_x(1)*g0_uu(1,2)+
     &                      gB_x(3)*phi10_x(1)*g0_uu(1,3)+
     &                      gB_x(3)*phi10_x(2)*g0_uu(2,3)
     &                                               )

                !---------------------------------------------------------------- 
                ! computes diag. Jacobian of g_np1->L.g_np1 transformation
                ! by differentiating L.g wrt. g(a,b)_ij_np1 diag. entries
                ! 
                ! ddgb_J_tx differs from ddgb_J due to forward/backward stencils
                ! at excision surfaces that affect the cross-derivative tx
                ! (this is the only contribution, since the diag Jacobian is diff wrt. g_ij_np1)
                !---------------------------------------------------------------- 
                dgb_J=1d0/2d0/dt
                ddgb_J=1d0/dt/dt

                if (i.eq.1.or.(chr(i-1,j).eq.ex)) then
                   if (i.le.(Nx-3)
     &                 .and.((chr(i+1,j).ne.ex
     &                 .and.chr(i+2,j).ne.ex
     &                 .and.chr(i+3,j).ne.ex))) then
                      ddgb_J_tx=-1d0/dt/dx
                   else if (i.le.(Nx-2)
     &                      .and.((chr(i+1,j).ne.ex
     &                      .and.chr(i+2,j).ne.ex))) then
                      ddgb_J_tx=-3d0/4d0/dt/dx
                   else if (i.le.(Nx-1).and.chr(i+1,j).ne.ex) then
                      ddgb_J_tx=-1d0/2d0/dt/dx
                   else
                      write(*,*) 'g_evo_opt: error in chr stencil (A)'
                      write(*,*) '    i,Nx,dx=',i,Nx,dx
                      write(*,*) '    (first error only)'
                      ddgb_J_tx=0
                   end if
                else if (i.eq.Nx.or.(chr(i+1,j).eq.ex)) then
                   if (i.ge.4
     &                 .and.((chr(i-1,j).ne.ex
     &                 .and.chr(i-2,j).ne.ex
     &                 .and.chr(i-3,j).ne.ex))) then
                      ddgb_J_tx=1d0/dt/dx
                   else if (i.ge.3
     &                      .and.((chr(i-1,j).ne.ex
     &                      .and.chr(i-2,j).ne.ex))) then
                      ddgb_J_tx=3d0/4d0/dt/dx
                   else if (i.ge.2.and.chr(i-1,j).ne.ex) then
                      ddgb_J_tx=1d0/2d0/dt/dx
                   else
                      write(*,*) 'g_evo_opt: error in chr stencil (B)'
                      write(*,*) '    i,Nx,dx=',i,Nx,dx
                      write(*,*) '    (first error only)'
                      ddgb_J_tx=0
                   end if
                else
                   if ((chr(i+1,j).ne.ex.and.chr(i-1,j).ne.ex)) then
                      ddgb_J_tx=0
                   else
                      write(*,*) 'g_evo_opt: error in chr stencil (C)'
                      write(*,*) '    i,Nx,dx=',i,Nx,dx
                      write(*,*) '    (first error only)'
                      ddgb_J_tx=0
                   end if
                end if

                if ((j.eq.1).or.(chr(i,j-1).eq.ex)) then
                   if (j.le.(Ny-3)
     &                 .and.((chr(i,j+1).ne.ex
     &                 .and.chr(i,j+2).ne.ex
     &                 .and.chr(i,j+3).ne.ex))) then
                      ddgb_J_ty=-1d0/dt/dy
                   else if (j.le.(Ny-2).and.((chr(i,j+1).ne.ex
     &                      .and.chr(i,j+2).ne.ex))) then
                      ddgb_J_ty=-3d0/4d0/dt/dy
                   else if (j.le.(Ny-1).and.chr(i,j+1).ne.ex) then
                      ddgb_J_ty=-1d0/2d0/dt/dy
                   else
                      write(*,*) 'g_evo_opt: error in chr stencil (D)'
                      write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
                      write(*,*) '    (first error only)'
                      ddgb_J_ty=0
                   end if
                else if ((j.eq.Ny).or.(chr(i,j+1).eq.ex)) then
                   if (j.ge.4
     &                 .and.((chr(i,j-1).ne.ex
     &                 .and.chr(i,j-2).ne.ex
     &                 .and.chr(i,j-3).ne.ex))) then
                      ddgb_J_ty=1d0/dt/dy
                   else if (j.ge.3.and.((chr(i,j-1).ne.ex
     &                      .and.chr(i,j-2).ne.ex))) then
                      ddgb_J_ty=3d0/4d0/dt/dy
                   else if (j.ge.2.and.chr(i,j-1).ne.ex) then
                      ddgb_J_ty=1d0/2d0/dt/dy
                   else
                      write(*,*) 'g_evo_opt: error in chr stencil (E)'
                      write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
                      write(*,*) '    (first error only)'
                      ddgb_J_ty=0
                   end if
                else
                   if ((chr(i,j+1).ne.ex.and.chr(i,j-1).ne.ex)) then
                      ddgb_J_ty=0
                   else
                      write(*,*) 'g_evo_opt: error in chr stencil (F)'
                      write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
                      write(*,*) '    (first error only)'
                      ddgb_J_ty=0
                   end if
                end if

                efe_J(1,1)=    -0.5d0*(
     &                            g0_uu(1,1)*ddgb_J
     &                                      *(1-x0**2)
     &                            -4*x0*g0_uu(1,2)*dgb_J
     &                            +2*g0_uu(1,2)*ddgb_J_tx
     &                                      *(1-x0**2)
     &                            +2*g0_uu(1,3)*ddgb_J_ty
     &                                      *(1-x0**2)
     &                                )
     &                      
     &                         -0.5d0*(
     &                            -dgb_J*(1-x0**2)*
     &                            (g0_uu(1,1)*g0_uu(1,1)*g0_ll_x(1,1,1)+
     &                             g0_uu(1,1)*g0_uu(2,1)*g0_ll_x(1,1,2)+
     &                             g0_uu(1,1)*g0_uu(3,1)*g0_ll_x(1,1,3)+
     &                             g0_uu(2,1)*g0_uu(1,1)*g0_ll_x(1,2,1)+
     &                             g0_uu(2,1)*g0_uu(2,1)*g0_ll_x(1,2,2)+
     &                             g0_uu(2,1)*g0_uu(3,1)*g0_ll_x(1,2,3)+
     &                             g0_uu(3,1)*g0_uu(1,1)*g0_ll_x(1,3,1)+
     &                             g0_uu(3,1)*g0_uu(2,1)*g0_ll_x(1,3,2)+
     &                             g0_uu(3,1)*g0_uu(3,1)*g0_ll_x(1,3,3))
     &                            +dgb_J*(1-x0**2)*
     &                            (g0_uu_x(1,1,1))
     &                                )*2
     &
     &                         +      (
     &                            0.5d0*dgb_J*(1-x0**2)*
     &                            ((Hads_l(1)+A_l(1))*g0_uu(1,1)+
     &                             (Hads_l(2)+A_l(2))*g0_uu(2,1)+
     &                             (Hads_l(3)+A_l(3))*g0_uu(3,1))
     &                                )  
     &                      
     &                         -      (
     &                            0.25d0*dgb_J*(1-x0**2)*
     &                            (cuuuu(1,1,1,1)*dlll(1,1,1)+
     &                             cuuuu(1,1,1,2)*dlll(1,1,2)+
     &                             cuuuu(1,1,1,3)*dlll(1,1,3)+
     &                             cuuuu(2,1,1,1)*dlll(2,1,1)+
     &                             cuuuu(2,1,1,2)*dlll(2,1,2)+
     &                             cuuuu(2,1,1,3)*dlll(2,1,3)+
     &                             cuuuu(3,1,1,1)*dlll(3,1,1)+
     &                             cuuuu(3,1,1,2)*dlll(3,1,2)+
     &                             cuuuu(3,1,1,3)*dlll(3,1,3)
     &                           +
     &                             cuuuu(1,1,1,1)*dlll(1,1,1)+
     &                             cuuuu(1,1,2,1)*dlll(2,1,1)+
     &                             cuuuu(1,1,3,1)*dlll(3,1,1)+
     &                             cuuuu(1,2,1,1)*dlll(1,1,2)+
     &                             cuuuu(1,2,2,1)*dlll(2,1,2)+
     &                             cuuuu(1,2,3,1)*dlll(3,1,2)+
     &                             cuuuu(1,3,1,1)*dlll(1,1,3)+
     &                             cuuuu(1,3,2,1)*dlll(2,1,3)+
     &                             cuuuu(1,3,3,1)*dlll(3,1,3))
     &                                )

                efe_J(1,2)=  -0.5d0*(
     &                          g0_uu(1,1)*ddgb_J
     &                                    *(1-x0**2)
     &                          -4*x0*g0_uu(1,2)*dgb_J
     &                          +2*g0_uu(1,2)*ddgb_J_tx
     &                                    *(1-x0**2)
     &                          +2*g0_uu(1,3)*ddgb_J_ty
     &                                    *(1-x0**2)
     &                              )
     &                       
     &                       -0.5d0*(
     &                          -dgb_J*(1-x0**2)*
     &                          (g0_uu(1,1)*g0_uu(1,2)*g0_ll_x(2,1,1)+
     &                           g0_uu(1,1)*g0_uu(2,2)*g0_ll_x(2,1,2)+
     &                           g0_uu(1,1)*g0_uu(3,2)*g0_ll_x(2,1,3)+
     &                           g0_uu(2,1)*g0_uu(1,2)*g0_ll_x(2,2,1)+
     &                           g0_uu(2,1)*g0_uu(2,2)*g0_ll_x(2,2,2)+
     &                           g0_uu(2,1)*g0_uu(3,2)*g0_ll_x(2,2,3)+
     &                           g0_uu(3,1)*g0_uu(1,2)*g0_ll_x(2,3,1)+
     &                           g0_uu(3,1)*g0_uu(2,2)*g0_ll_x(2,3,2)+
     &                           g0_uu(3,1)*g0_uu(3,2)*g0_ll_x(2,3,3))
     &                          +dgb_J*(1-x0**2)*
     &                          (g0_uu_x(1,1,1))
     &                              )
     &                    
     &                       -0.5d0*(
     &                          dgb_J*(1-x0**2)*
     &                          (g0_uu_x(2,1,2))
     &                              )
     &
     &                       -      (
     &                          0.5d0*dgb_J*(1-x0**2)*
     &                          (cuuuu(1,1,1,2)*dlll(1,2,1)+
     &                           cuuuu(1,1,2,2)*dlll(2,2,1)+
     &                           cuuuu(1,1,3,2)*dlll(3,2,1)+
     &                           cuuuu(1,2,1,2)*dlll(1,2,2)+
     &                           cuuuu(1,2,2,2)*dlll(2,2,2)+
     &                           cuuuu(1,2,3,2)*dlll(3,2,2)+
     &                           cuuuu(1,3,1,2)*dlll(1,2,3)+
     &                           cuuuu(1,3,2,2)*dlll(2,2,3)+
     &                           cuuuu(1,3,3,2)*dlll(3,2,3))
     &                              )

                efe_J(1,3)=  -0.5d0*(
     &                          g0_uu(1,1)*ddgb_J
     &                                    *(1-x0**2)**2
     &                          -4*(2)*x0*g0_uu(1,2)*dgb_J
     &                                    *(1-x0**2)
     &                          +2*g0_uu(1,2)*ddgb_J_tx
     &                                    *(1-x0**2)**2
     &                          +2*g0_uu(1,3)*ddgb_J_ty
     &                                    *(1-x0**2)**2
     &                              )
     &
     &                       -0.5d0*(
     &                          -dgb_J*(1-x0**2)**2*
     &                          (g0_uu(1,1)*g0_uu(1,3)*g0_ll_x(3,1,1)+
     &                           g0_uu(1,1)*g0_uu(2,3)*g0_ll_x(3,1,2)+
     &                           g0_uu(1,1)*g0_uu(3,3)*g0_ll_x(3,1,3)+
     &                           g0_uu(2,1)*g0_uu(1,3)*g0_ll_x(3,2,1)+
     &                           g0_uu(2,1)*g0_uu(2,3)*g0_ll_x(3,2,2)+
     &                           g0_uu(2,1)*g0_uu(3,3)*g0_ll_x(3,2,3)+
     &                           g0_uu(3,1)*g0_uu(1,3)*g0_ll_x(3,3,1)+
     &                           g0_uu(3,1)*g0_uu(2,3)*g0_ll_x(3,3,2)+
     &                           g0_uu(3,1)*g0_uu(3,3)*g0_ll_x(3,3,3))
     &                          +dgb_J*(1-x0**2)**2*
     &                          (g0_uu_x(1,1,1))
     &                              )
     &
     &                       -0.5d0*(
     &                          dgb_J*(1-x0**2)**2*
     &                          (g0_uu_x(3,1,3))
     &                              )
     &
     &                       -      (
     &                          0.5d0*dgb_J*(1-x0**2)**2*
     &                          (cuuuu(1,1,1,3)*dlll(1,3,1)+
     &                           cuuuu(1,1,2,3)*dlll(2,3,1)+
     &                           cuuuu(1,1,3,3)*dlll(3,3,1)+
     &                           cuuuu(1,2,1,3)*dlll(1,3,2)+
     &                           cuuuu(1,2,2,3)*dlll(2,3,2)+
     &                           cuuuu(1,2,3,3)*dlll(3,3,2)+
     &                           cuuuu(1,3,1,3)*dlll(1,3,3)+
     &                           cuuuu(1,3,2,3)*dlll(2,3,3)+
     &                           cuuuu(1,3,3,3)*dlll(3,3,3))
     &                              )

                efe_J(2,2)=-0.5d0*(
     &                        g0_uu(1,1)*ddgb_J
     &                                  *(1-x0**2)
     &                        -4*x0*g0_uu(1,2)*dgb_J
     &                        +2*g0_uu(1,2)*ddgb_J_tx
     &                                  *(1-x0**2)
     &                        +2*g0_uu(1,2)*ddgb_J_ty
     &                                  *(1-x0**2)
     &                            )
     &                  
     &                     -0.5d0*(
     &                        +dgb_J*(1-x0**2)*
     &                        (g0_uu_x(2,1,2))
     &                            )
     &                  
     &                     -0.5d0*(
     &                        +dgb_J*(1-x0**2)*
     &                        (g0_uu_x(2,1,2))
     &                            )
     &
     &                     +      (
     &                        -0.5d0*dgb_J*(1-x0**2)*
     &                        ((Hads_l(1)+A_l(1))*g0_uu(1,1)+
     &                         (Hads_l(2)+A_l(2))*g0_uu(2,1)+
     &                         (Hads_l(3)+A_l(3))*g0_uu(3,1))
     &                            )  
     &                  
     &                     -      (
     &                        0.25d0*dgb_J*(1-x0**2)*
     &                        (cuuuu(1,2,1,1)*dlll(1,2,1)+
     &                         cuuuu(1,2,1,2)*dlll(1,2,2)+
     &                         cuuuu(1,2,1,3)*dlll(1,2,3)+
     &                         cuuuu(2,2,1,1)*dlll(2,2,1)+
     &                         cuuuu(2,2,1,2)*dlll(2,2,2)+
     &                         cuuuu(2,2,1,3)*dlll(2,2,3)+
     &                         cuuuu(3,2,1,1)*dlll(3,2,1)+
     &                         cuuuu(3,2,1,2)*dlll(3,2,2)+
     &                         cuuuu(3,2,1,3)*dlll(3,2,3)+
     &                         cuuuu(1,1,2,1)*dlll(1,2,1)-
     &                         cuuuu(1,1,2,2)*dlll(1,2,2)-
     &                         cuuuu(1,1,2,3)*dlll(1,2,3)-
     &                         cuuuu(2,1,2,1)*dlll(2,2,1)-
     &                         cuuuu(2,1,2,2)*dlll(2,2,2)-
     &                         cuuuu(2,1,2,3)*dlll(2,2,3)-
     &                         cuuuu(3,1,2,1)*dlll(3,2,1)-
     &                         cuuuu(3,1,2,2)*dlll(3,2,2)-
     &                         cuuuu(3,1,2,3)*dlll(3,2,3)
     &                       +
     &                         cuuuu(1,1,1,2)*dlll(1,2,1)+
     &                         cuuuu(1,1,2,2)*dlll(2,2,1)+
     &                         cuuuu(1,1,3,2)*dlll(3,2,1)+
     &                         cuuuu(1,2,1,2)*dlll(1,2,2)+
     &                         cuuuu(1,2,2,2)*dlll(2,2,2)+
     &                         cuuuu(1,2,3,2)*dlll(3,2,2)+
     &                         cuuuu(1,3,1,2)*dlll(1,2,3)+
     &                         cuuuu(1,3,2,2)*dlll(2,2,3)+
     &                         cuuuu(1,3,3,2)*dlll(3,2,3)+
     &                         cuuuu(2,1,1,1)*dlll(1,2,1)-
     &                         cuuuu(2,1,2,1)*dlll(2,2,1)-
     &                         cuuuu(2,1,3,1)*dlll(3,2,1)-
     &                         cuuuu(2,2,1,1)*dlll(1,2,2)-
     &                         cuuuu(2,2,2,1)*dlll(2,2,2)-
     &                         cuuuu(2,2,3,1)*dlll(3,2,2)-
     &                         cuuuu(2,3,1,1)*dlll(1,2,3)-
     &                         cuuuu(2,3,2,1)*dlll(2,2,3)-
     &                         cuuuu(2,3,3,1)*dlll(3,2,3))
     &                            )

                  efe_J(2,3)=-0.5d0*(
     &                          g0_uu(1,1)*ddgb_J
     &                                    *(1-x0**2)**2
     &                          -4*(2)*x0*g0_uu(1,2)*dgb_J
     &                                    *(1-x0**2)
     &                          +2*g0_uu(1,2)*ddgb_J_tx
     &                                    *(1-x0**2)**2
     &                          +2*g0_uu(1,3)*ddgb_J_ty
     &                                    *(1-x0**2)**2
     &                              )
     &
     &                       -0.5d0*(
     &                          +dgb_J*(1-x0**2)**2*
     &                          (g0_uu_x(2,1,2))
     &                              )
     &
     &                       -0.5d0*(
     &                          +dgb_J*(1-x0**2)**2*
     &                          (g0_uu_x(3,1,3))
     &                              )
     &
     &                       +      (
     &                          -0.5d0*dgb_J*(1-x0**2)**2*
     &                          ((Hads_l(1)+A_l(1))*g0_uu(1,1)+
     &                           (Hads_l(2)+A_l(2))*g0_uu(2,1)+
     &                           (Hads_l(3)+A_l(3))*g0_uu(3,1))
     &                              )
     &
     &                       -      (
     &                          0.25d0*dgb_J*(1-x0**2)**2*
     &                          (cuuuu(1,2,1,1)*dlll(1,2,1)+
     &                           cuuuu(1,2,1,2)*dlll(1,2,2)+
     &                           cuuuu(1,2,1,3)*dlll(1,2,3)+
     &                           cuuuu(2,2,1,1)*dlll(2,2,1)+
     &                           cuuuu(2,2,1,2)*dlll(2,2,2)+
     &                           cuuuu(2,2,1,3)*dlll(2,2,3)+
     &                           cuuuu(3,2,1,1)*dlll(3,2,1)+
     &                           cuuuu(3,2,1,2)*dlll(3,2,2)+
     &                           cuuuu(3,2,1,3)*dlll(3,2,3)+
     &                           cuuuu(1,1,2,1)*dlll(1,2,1)-
     &                           cuuuu(1,1,2,2)*dlll(1,2,2)-
     &                           cuuuu(1,1,2,3)*dlll(1,2,3)-
     &                           cuuuu(2,1,2,1)*dlll(2,2,1)-
     &                           cuuuu(2,1,2,2)*dlll(2,2,2)-
     &                           cuuuu(2,1,2,3)*dlll(2,2,3)-
     &                           cuuuu(3,1,2,1)*dlll(3,2,1)-
     &                           cuuuu(3,1,2,2)*dlll(3,2,2)-
     &                           cuuuu(3,1,2,3)*dlll(3,2,3)
     &                         +
     &                           cuuuu(1,1,1,3)*dlll(1,3,1)+
     &                           cuuuu(1,1,2,3)*dlll(2,3,1)+
     &                           cuuuu(1,1,3,3)*dlll(3,3,1)+
     &                           cuuuu(1,2,1,3)*dlll(1,3,2)+
     &                           cuuuu(1,2,2,3)*dlll(2,3,2)+
     &                           cuuuu(1,2,3,3)*dlll(3,3,2)+
     &                           cuuuu(1,3,1,3)*dlll(1,3,3)+
     &                           cuuuu(1,3,2,3)*dlll(2,3,3)+
     &                           cuuuu(1,3,3,3)*dlll(3,3,3)+
     &                           cuuuu(3,1,1,1)*dlll(1,3,1)-
     &                           cuuuu(3,1,2,1)*dlll(2,3,1)-
     &                           cuuuu(3,1,3,1)*dlll(3,3,1)-
     &                           cuuuu(3,2,1,1)*dlll(1,3,2)-
     &                           cuuuu(3,2,2,1)*dlll(2,3,2)-
     &                           cuuuu(3,2,3,1)*dlll(3,3,2)-
     &                           cuuuu(3,3,1,1)*dlll(1,3,3)-
     &                           cuuuu(3,3,2,1)*dlll(2,3,3)-
     &                           cuuuu(3,3,3,1)*dlll(3,3,3))
     &                              )

                efe_J(3,3)=    -0.5d0*(
     &                            g0_uu(1,1)*ddgb_J
     &                                          *(1-x0**2)**3
     &                            +2*(3*(-2*x0))*g0_uu(1,2)*dgb_J
     &                                          *(1-x0**2)**2
     &                            +2*g0_uu(1,2)*ddgb_J_tx
     &                                          *(1-x0**2)**3
     &                            +2*g0_uu(1,3)*ddgb_J_ty
     &                                          *(1-x0**2)**3
     &                                )
     &                      
     &                         -0.5d0*(
     &                            +dgb_J
     &                                        *(g0_uu_x(3,1,3))
     &                                        *(1-x0**2)**3
     &                                )
     &                      
     &                         -0.5d0*(
     &                            +dgb_J
     &                                        *(g0_uu_x(3,1,3))
     &                                        *(1-x0**2)**3
     &                                )
     &
     &                         +      (
     &                            -0.5d0*dgb_J
     &                                        *(1-x0**2)**3*
     &                            ((Hads_l(1)+A_l(1))*g0_uu(1,1)+
     &                             (Hads_l(2)+A_l(2))*g0_uu(2,1)+
     &                             (Hads_l(3)+A_l(3))*g0_uu(3,1))
     &                                )  
     &                      
     &                         -      (
     &                            0.25d0*dgb_J
     &                                        *(1-x0**2)**3*
     &                            (cuuuu(1,3,1,1)*dlll(1,3,1)+
     &                             cuuuu(1,3,1,2)*dlll(1,3,2)+
     &                             cuuuu(1,3,1,3)*dlll(1,3,3)+
     &                             cuuuu(2,3,1,1)*dlll(2,3,1)+
     &                             cuuuu(2,3,1,2)*dlll(2,3,2)+
     &                             cuuuu(2,3,1,3)*dlll(2,3,3)+
     &                             cuuuu(3,3,1,1)*dlll(3,3,1)+
     &                             cuuuu(3,3,1,2)*dlll(3,3,2)+
     &                             cuuuu(3,3,1,3)*dlll(3,3,3)+
     &                             cuuuu(1,1,3,1)*dlll(1,3,1)-
     &                             cuuuu(1,1,3,2)*dlll(1,3,2)-
     &                             cuuuu(1,1,3,3)*dlll(1,3,3)-
     &                             cuuuu(2,1,3,1)*dlll(2,3,1)-
     &                             cuuuu(2,1,3,2)*dlll(2,3,2)-
     &                             cuuuu(2,1,3,3)*dlll(2,3,3)-
     &                             cuuuu(3,1,3,1)*dlll(3,3,1)-
     &                             cuuuu(3,1,3,2)*dlll(3,3,2)-
     &                             cuuuu(3,1,3,3)*dlll(3,3,3)
     &                           +
     &                             cuuuu(1,1,1,3)*dlll(1,3,1)+
     &                             cuuuu(1,1,2,3)*dlll(2,3,1)+
     &                             cuuuu(1,1,3,3)*dlll(3,3,1)+
     &                             cuuuu(1,2,1,3)*dlll(1,3,2)+
     &                             cuuuu(1,2,2,3)*dlll(2,3,2)+
     &                             cuuuu(1,2,3,3)*dlll(3,3,2)+
     &                             cuuuu(1,3,1,3)*dlll(1,3,3)+
     &                             cuuuu(1,3,2,3)*dlll(2,3,3)+
     &                             cuuuu(1,3,3,3)*dlll(3,3,3)+
     &                             cuuuu(3,1,1,1)*dlll(1,3,1)-
     &                             cuuuu(3,1,2,1)*dlll(2,3,1)-
     &                             cuuuu(3,1,3,1)*dlll(3,3,1)-
     &                             cuuuu(3,2,1,1)*dlll(1,3,2)-
     &                             cuuuu(3,2,2,1)*dlll(2,3,2)-
     &                             cuuuu(3,2,3,1)*dlll(3,3,2)-
     &                             cuuuu(3,3,1,1)*dlll(1,3,3)-
     &                             cuuuu(3,3,2,1)*dlll(2,3,3)-
     &                             cuuuu(3,3,3,1)*dlll(3,3,3))
     &                                )

                !---------------------------------------------------------------- 
                ! computes diag. Jacobian of psi_np1->L.psi_np1 transformation
                ! by differentiating L.psi wrt. psi_ij_np1 diag. entries
                ! and remember: gA=gA_ads0+psi*(1-x0**2)*x0**2
                !---------------------------------------------------------------- 

                afe_J=(
     &                -g0_uu(1,1)*ddgb_J/2
     &                +g0_uu(1,1)*gamma_ull(1,1,1)*dgb_J/2
     &                +g0_uu(1,2)*gamma_ull(1,2,1)*dgb_J/2
     &                +g0_uu(1,3)*gamma_ull(1,3,1)*dgb_J/2
     &                +g0_uu(2,1)*gamma_ull(1,1,2)*dgb_J/2
     &                +g0_uu(2,2)*gamma_ull(1,2,2)*dgb_J/2
     &                +g0_uu(2,3)*gamma_ull(1,3,2)*dgb_J/2
     &                +g0_uu(3,1)*gamma_ull(1,1,3)*dgb_J/2
     &                +g0_uu(3,2)*gamma_ull(1,2,3)*dgb_J/2
     &                +g0_uu(3,3)*gamma_ull(1,3,3)*dgb_J/2
     &                -g0_uu(1,1)*(2*dgb_J*gA_x(1))/(4*gA)*(dimA-2d0)
     &                -g0_uu(1,2)*(dgb_J*gA_x(2))/(4*gA)*(dimA-2d0)
     &                -g0_uu(1,3)*(dgb_J*gA_x(3))/(4*gA)*(dimA-2d0)
     &                -g0_uu(2,1)*(gA_x(2)*dgb_J)/(4*gA)*(dimA-2d0)
     &                -g0_uu(3,1)*(gA_x(3)*dgb_J)/(4*gA)*(dimA-2d0)
     &                -g0_uu(1,1)*(dgb_J*gB_x(1))/(4*gB)*dimB
     &                -g0_uu(1,2)*(dgb_J*gB_x(2))/(4*gB)*dimB
     &                -g0_uu(1,3)*(dgb_J*gB_x(3))/(4*gB)*dimB
     &                )*(1-x0**2)*x0**2

                !---------------------------------------------------------------- 
                ! computes diag. Jacobian of gB_np1->L.gB_np1 transformation
                ! by differentiating L.gB wrt. gB_ij_np1 diag. entries
                ! and remember: gB=gB_ads0+omega*(1-x0**2)**3
                !---------------------------------------------------------------- 

                bfe_J=(
     &                -g0_uu(1,1)*ddgb_J/2
     &                +g0_uu(1,1)*gamma_ull(1,1,1)*dgb_J/2
     &                +g0_uu(1,2)*gamma_ull(1,2,1)*dgb_J/2
     &                +g0_uu(1,3)*gamma_ull(1,3,1)*dgb_J/2
     &                +g0_uu(2,1)*gamma_ull(1,1,2)*dgb_J/2
     &                +g0_uu(2,2)*gamma_ull(1,2,2)*dgb_J/2
     &                +g0_uu(2,3)*gamma_ull(1,3,2)*dgb_J/2
     &                +g0_uu(3,1)*gamma_ull(1,1,3)*dgb_J/2
     &                +g0_uu(3,2)*gamma_ull(1,2,3)*dgb_J/2
     &                +g0_uu(3,3)*gamma_ull(1,3,3)*dgb_J/2
     &                -g0_uu(1,1)*2*dgb_J*gB_x(1)/(4*gB)*(dimB-2d0)
     &                -g0_uu(1,2)*dgb_J*gB_x(2)/(4*gB)*(dimB-2d0)
     &                -g0_uu(1,3)*dgb_J*gB_x(3)/(4*gB)*(dimB-2d0)
     &                -g0_uu(2,1)*gB_x(2)*dgb_J/(4*gB)*(dimB-2d0)
     &                -g0_uu(3,1)*gB_x(3)*dgb_J/(4*gB)*(dimB-2d0)
     &                -g0_uu(1,1)*dgb_J*gA_x(1)/(4*gA)*dimA 
     &                -g0_uu(1,2)*dgb_J*gA_x(2)/(4*gA)*dimA
     &                -g0_uu(1,3)*dgb_J*gA_x(3)/(4*gA)*dimA 
     &                )*(1-x0**2)**3

                !----------------------------------------------------------------
                ! computes diag. Jacobian of phi1_np1->L.phi1_np1 transformation
                ! by differentiating L.phi1=box.phi1-dV/dphi1 wrt. phi1_np1
                ! and remember: phi10=phi1*(1-x0**2)**3*sin(PI*y0)**5 
                ! 
                ! ddphi1_J_tx differs from ddphi_J due to forward/backward stencils
                ! at excision surfaces that affect the cross-derivative tx 
                ! (these are the only contributions, since the diag Jacobian is diff wrt. phi1_np1)
                !----------------------------------------------------------------
                dphi1_J=1d0/2d0/dt
                ddphi1_J=1d0/dt/dt

                if (i.eq.1.or.(chr(i-1,j).eq.ex)) then
                   if (i.le.(Nx-3)
     &                 .and.((chr(i+1,j).ne.ex
     &                 .and.chr(i+2,j).ne.ex
     &                 .and.chr(i+3,j).ne.ex))) then
                      ddphi1_J_tx=-1d0/dt/dx
                   else if (i.le.(Nx-2)
     &                      .and.((chr(i+1,j).ne.ex
     &                      .and.chr(i+2,j).ne.ex))) then
                      ddphi1_J_tx=-3d0/4d0/dt/dx
                   else if (i.le.(Nx-1).and.chr(i+1,j).ne.ex) then
                      ddphi1_J_tx=-1d0/2d0/dt/dx
                   else
                      write(*,*) 'g_evo_opt: error in chr stencil (A)'
                      write(*,*) '    i,Nx,dx=',i,Nx,dx
                      write(*,*) '    (first error only)'
                      ddphi1_J_tx=0
                   end if
                else if (i.eq.Nx.or.(chr(i+1,j).eq.ex)) then
                   if (i.ge.4
     &                 .and.((chr(i-1,j).ne.ex
     &                 .and.chr(i-2,j).ne.ex
     &                 .and.chr(i-3,j).ne.ex))) then
                      ddphi1_J_tx=1d0/dt/dx
                   else if (i.ge.3
     &                      .and.((chr(i-1,j).ne.ex
     &                      .and.chr(i-2,j).ne.ex))) then
                      ddphi1_J_tx=3d0/4d0/dt/dx
                   else if (i.ge.2.and.chr(i-1,j).ne.ex) then
                      ddphi1_J_tx=1d0/2d0/dt/dx
                   else
                      write(*,*) 'g_evo_opt: error in chr stencil (B)'
                      write(*,*) '    i,Nx,dx=',i,Nx,dx
                      write(*,*) '    (first error only)'
                      ddphi1_J_tx=0
                   end if
                else
                   if ((chr(i+1,j).ne.ex.and.chr(i-1,j).ne.ex)) then
                      ddphi1_J_tx=0
                   else
                      write(*,*) 'g_evo_opt: error in chr stencil (C)'
                      write(*,*) '    i,Nx,dx=',i,Nx,dx
                      write(*,*) '    (first error only)'
                      ddphi1_J_tx=0
                   end if
                end if

                phi1_J=    (
     &                 g0_uu(1,1)*ddphi1_J
     &                           *(1-x0**2)**3*sin(PI*y0)**5
     &                 -4*(3)*x0*g0_uu(1,2)*dphi1_J
     &                           *(1-x0**2)**3*sin(PI*y0)**5
     &                 +2*g0_uu(1,2)*ddphi1_J_tx
     &                           *(1-x0**2)**3*sin(PI*y0)**5
     &                     )
     &                +
     &                  dphi1_J*(1-x0**2)**3*sin(PI*y0)**5
     &                  *(
     &                    g0_uu_x(1,1,1)+
     &                    g0_uu_x(2,1,2)+
     &                    g0_uu_x(3,1,3)
     &                   )
     &                +
     &                  dphi1_J*(1-x0**2)**3*sin(PI*y0)**5
     &                  *(
     &                    g0_uu(1,1)*gamma_ull(1,1,1)+
     &                    g0_uu(1,1)*gamma_ull(2,2,1)+
     &                    g0_uu(1,1)*gamma_ull(3,3,1)+
     &                    g0_uu(1,2)*gamma_ull(1,1,2)+
     &                    g0_uu(1,2)*gamma_ull(2,2,2)+
     &                    g0_uu(1,2)*gamma_ull(3,3,2)+
     &                    g0_uu(1,3)*gamma_ull(1,1,3)+
     &                    g0_uu(1,3)*gamma_ull(2,2,3)+
     &                    g0_uu(1,3)*gamma_ull(3,3,3)
     &                   )
     &                +
     &                   dphi1_J*(1-x0**2)**3*sin(PI*y0)**5
     &                          *1.5d0/gA*(
     &                    gA_x(1)*g0_uu(1,1)+
     &                    gA_x(2)*g0_uu(1,2)+
     &                    gA_x(3)*g0_uu(1,3)
     &                                             )
     &                -
     &                   dphi1_J*(1-x0**2)**3*sin(PI*y0)**5
     &                          *2.0d0/gB*(
     &                    gB_x(1)*g0_uu(1,1)+
     &                    gB_x(2)*g0_uu(1,2)+
     &                    gB_x(3)*g0_uu(1,3)
     &                                             )

                ! constraint damping terms added to efe,efe_J
                do a=1,3
                  do b=1,3
                    cd_ll(a,b)=-kappa_cd*
     &                  ( n_l(a)*c_l(b)+n_l(b)*c_l(a)
     &                    -(1+rho_cd)*g0_ll(a,b)*ndotc )
                  end do
                end do

                ! differentiating cd(a,b) wrt. g(a,b)_ij_np1 diag. entries
                dc_J=1d0/2d0/dt

                cd_J_ll(1,1)=-kappa_cd*
     &              ( n_l(1)*g0_uu(1,1)*dc_J*(1-x0**2)
     &                -(1+rho_cd)*g0_ll(1,1)*n_u(1)*0.5d0*
     &                            g0_uu(1,1)*dc_J*(1-x0**2) )
                cd_J_ll(1,2)=-kappa_cd*
     &              ( n_l(1)*g0_uu(1,1)*dc_J*(1-x0**2)
     &                -(1+rho_cd)*g0_ll(1,2)*n_u(2)*g0_uu(1,1)*
     &                            dc_J*(1-x0**2) )
                cd_J_ll(1,3)=-kappa_cd*
     &              ( n_l(1)*g0_uu(1,1)*dc_J*(1-x0**2)**2
     &                -(1+rho_cd)*g0_ll(1,3)*n_u(3)*g0_uu(1,1)*
     &                            dc_J*(1-x0**2)**2 )
                cd_J_ll(2,2)=-kappa_cd*
     &              ( 2*n_l(2)*g0_uu(1,2)*dc_J*(1-x0**2)
     &                -(1+rho_cd)*g0_ll(2,2)*
     &                (-n_u(1)*0.5d0*g0_uu(2,2)*dc_J*(1-x0**2)
     &                 +n_u(2)*g0_uu(1,2)*dc_J*(1-x0**2)) )
                cd_J_ll(2,3)=-kappa_cd*
     &              ( n_l(2)*g0_uu(1,2)*dc_J*(1-x0**2)**2
     &               +n_l(3)*g0_uu(1,3)*dc_J*(1-x0**2)**2
     &                -(1+rho_cd)*g0_ll(2,3)*
     &                (-n_u(1)*0.5d0*g0_uu(2,3)*dc_J*(1-x0**2)**2
     &                 +n_u(2)*g0_uu(1,3)*dc_J*(1-x0**2)**2
     &                 +n_u(3)*g0_uu(1,2)*dc_J*(1-x0**2)**2) ) 
                cd_J_ll(3,3)=-kappa_cd*
     &              ( 2*n_l(3)*g0_uu(1,3)*dc_J*(1-x0**2)**3
     &                -(1+rho_cd)*g0_ll(3,3)*
     &                (-n_u(1)*0.5d0*g0_uu(3,3)*dc_J*(1-x0**2)**3
     &                 +n_u(3)*g0_uu(1,3)*dc_J*(1-x0**2)**3) )
 
                if (kappa_cd.ne.0) then
                  efe(1,1)=efe(1,1)+cd_ll(1,1)
                  efe(1,2)=efe(1,2)+cd_ll(1,2)
                  efe(1,3)=efe(1,3)+cd_ll(1,3)
                  efe(2,2)=efe(2,2)+cd_ll(2,2)
                  efe(2,3)=efe(2,3)+cd_ll(2,3)
                  efe(3,3)=efe(3,3)+cd_ll(3,3)
                  efe_J(1,1)=efe_J(1,1)+cd_J_ll(1,1)
                  efe_J(1,2)=efe_J(1,2)+cd_J_ll(1,2)
                  efe_J(1,3)=efe_J(1,3)+cd_J_ll(1,3)
                  efe_J(2,2)=efe_J(2,2)+cd_J_ll(2,2)
                  efe_J(2,3)=efe_J(2,3)+cd_J_ll(2,3)
                  efe_J(3,3)=efe_J(3,3)+cd_J_ll(3,3)
                end if

                ! update gbars 
                if (is_nan(efe(1,1)).or.is_nan(efe_J(1,1)).or.
     &            efe_J(1,1).eq.0) then
                  dump=.true.
                else
                  gb_tt_np1(i,j)=gb_tt_np1(i,j)-efe(1,1)/efe_J(1,1)
                end if

                if (is_nan(efe(1,2)).or.is_nan(efe_J(1,2)).or.
     &            efe_J(1,2).eq.0) then
                  dump=.true.
                else
                  gb_tx_np1(i,j)=gb_tx_np1(i,j)-efe(1,2)/efe_J(1,2)
                end if

                if (is_nan(efe(1,3)).or.is_nan(efe_J(1,3)).or.
     &            efe_J(1,3).eq.0) then
                  dump=.true.
                else
                  gb_ty_np1(i,j)=gb_ty_np1(i,j)-efe(1,3)/efe_J(1,3)
                end if

                if (is_nan(efe(2,2)).or.is_nan(efe_J(2,2)).or.
     &            efe_J(2,2).eq.0) then
                  dump=.true.
                else
                  gb_xx_np1(i,j)=gb_xx_np1(i,j)-efe(2,2)/efe_J(2,2)
                end if

                if (is_nan(efe(2,3)).or.is_nan(efe_J(2,3)).or.
     &            efe_J(2,3).eq.0) then
                  dump=.true.
                else
                  gb_xy_np1(i,j)=gb_xy_np1(i,j)-efe(2,3)/efe_J(2,3)
                end if

                if (is_nan(efe(3,3)).or.is_nan(efe_J(3,3)).or.
     &            efe_J(3,3).eq.0) then
                  dump=.true.
                else
                  gb_yy_np1(i,j)=gb_yy_np1(i,j)-efe(3,3)/efe_J(3,3)
                end if

                ! update psi
                if (is_nan(afe).or.is_nan(afe_J).or.
     &            afe_J.eq.0) then
                  dump=.true.
                else
                  psi_np1(i,j)=psi_np1(i,j)-afe/afe_J
                end if

                ! update omega
                if (is_nan(bfe).or.is_nan(bfe_J).or.
     &            bfe_J.eq.0) then
                  dump=.true.
                else
                  omega_np1(i,j)=omega_np1(i,j)-bfe/bfe_J
                end if

                ! update phi1 
                if (is_nan(phi1_res).or.is_nan(phi1_J)) then
                  dump=.true.
                else
                  phi1_np1(i,j)=phi1_np1(i,j)-phi1_res/phi1_J 
                end if

                gb_res(i,j) =
     &            max(abs(efe(1,1)/efe_J(1,1)),
     &                abs(efe(1,2)/efe_J(1,2)),
     &                abs(efe(1,3)/efe_J(1,3)),
     &                abs(efe(2,2)/efe_J(2,2)),
     &                abs(efe(2,3)/efe_J(2,3)),
     &                abs(efe(3,3)/efe_J(3,3)),
     &                abs(afe/afe_J),
     &                abs(bfe/bfe_J))
                cl_res(i,j)=max(abs(c_l(1)),abs(c_l(2)),abs(c_l(3)))

                ! check for NaNs
                if (dump.and.first_nan.or.ltrace) then
                  first_nan=.false.
                  write(*,*)
                  write(*,*) 'g_evo_opt: Nan/zero at i,j,Nx,Ny=',
     &                                              i,j,Nx,Ny
                  write(*,*) 'x,y=',x(i),y(j)
                  write(*,*) 'dt,dx,dy=',dt,dx,dy
                  write(*,*) 'x0,y0',x0,y0

                  write(*,*) ' at tn:'
                  write(*,*) ' gb_tt np1,n,nm1:',gb_tt_np1(i,j),
     &                   gb_tt_n(i,j),gb_tt_nm1(i,j)
                  write(*,*) ' gb_tx np1,n,nm1:',gb_tx_np1(i,j),
     &                   gb_tx_n(i,j),gb_tx_nm1(i,j)
                  write(*,*) ' gb_xx np1,n,nm1:',gb_xx_np1(i,j),
     &                   gb_xx_n(i,j),gb_xx_nm1(i,j)
                  write(*,*) ' psi np1,n,nm1:',psi_np1(i,j),
     &                   psi_n(i,j),psi_nm1(i,j)
                  write(*,*) ' g0_11 :',g0_ll(1,1)
                  write(*,*) ' g0_12 :',g0_ll(1,2)
                  write(*,*) ' g0_22 :',g0_ll(2,2)
                  write(*,*) ' g0_33:',g0_ll(3,3)
                  write(*,*) ' g0u_11 :',g0_uu(1,1)
                  write(*,*) ' g0u_12 :',g0_uu(1,2)
                  write(*,*) ' g0u_22 :',g0_uu(2,2)
                  write(*,*) ' g0u_33:',g0_uu(3,3)
                  write(*,*) ' phi np1,n,nm1:',phi1_np1(i,j),
     &                     phi1_n(i,j),phi1_nm1(i,j)
                  write(*,*) ' phi1_t/x:',phi10_x(1),phi10_x(2)
                  write(*,*) ' phi1_ti..:',phi10_xx(1,1),phi10_xx(1,2)

                  write(*,*) ' res J:'
                  write(*,*) ' tt:',efe(1,1),efe_J(1,1)
                  write(*,*) ' tx:',efe(1,2),efe_J(1,2)
                  write(*,*) ' xx:',efe(2,2),efe_J(2,2)
                  write(*,*) ' psi:',efe(3,3),efe_J(3,3)
                  write(*,*) ' phi1:',phi1_res,phi1_J
                end if

              ! (REGION) non-interior points; set to zero prior to applying bcs 
              else 
                gb_tt_np1(i,j) = 0
                gb_tx_np1(i,j) = 0
                gb_ty_np1(i,j) = 0
                gb_xx_np1(i,j) = 0
                gb_xy_np1(i,j) = 0
                gb_yy_np1(i,j) = 0
                psi_np1(i,j)   = 0
                omega_np1(i,j) = 0
                phi1_np1(i,j)  = 0 

              endif ! (near start of main loop)

            end do
          end do

        end do

        ! (REGION) x=0; impose regularity conditions 
        call axi_reg_phi(phi1_np1,chr,ex,L,x,y,Nx,Ny)
        call axi_reg_g(gb_tt_np1,gb_tx_np1,gb_ty_np1,gb_xx_np1,
     &                 gb_xy_np1,gb_yy_np1,psi_np1,omega_np1,
     &                 chr,ex,L,x,y,Nx,Ny)
        call axi_reg_Hb(Hb_t_np1,Hb_x_np1,Hb_y_np1,chr,ex,L,x,y,Nx,Ny)

        return
        end

c----------------------------------------------------------------------
        logical function is_nan(x)
        implicit none
        real*8 x

        integer is_a_nan

        call check_nan(x,is_a_nan)

        if (is_a_nan.eq.0) then
           is_nan=.false.
        else
           is_nan=.true.
        end if

        return
        end
c----------------------------------------------------------------------
