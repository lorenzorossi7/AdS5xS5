c----------------------------------------------------------------------
c this routine initializes the past time level to O(h^3),
c given [gb,phi] at t=t0, and d[[gb,phi]]/dt
c No gauge evolution considered for H yet.
c
c note that at interior boundaries the spatial derivatives
c are calculated using backwards/forwards stencils. (can't 
c inject from the parent level, as they are out of sync)
c----------------------------------------------------------------------
        subroutine init_nm1(gb_tt_np1,gb_tt_n,gb_tt_nm1,gb_tt_t_n,
     &                      gb_tx_np1,gb_tx_n,gb_tx_nm1,gb_tx_t_n,
     &                      gb_ty_np1,gb_ty_n,gb_ty_nm1,gb_ty_t_n,
     &                      gb_xx_np1,gb_xx_n,gb_xx_nm1,gb_xx_t_n,
     &                      gb_xy_np1,gb_xy_n,gb_xy_nm1,gb_xy_t_n,
     &                      gb_yy_np1,gb_yy_n,gb_yy_nm1,gb_yy_t_n,
     &                      psi_np1,psi_n,psi_nm1,psi_t_n,
     &                      omega_np1,omega_n,omega_nm1,omega_t_n,
     &                      fb_t_np1,fb_t_n,fb_t_nm1,fb_t_t_n,
     &                      fb_x_np1,fb_x_n,fb_x_nm1,fb_x_t_n,
     &                      fb_y_np1,fb_y_n,fb_y_nm1,fb_y_t_n,
     &                      Hb_t_np1,Hb_t_n,Hb_t_nm1,Hb_t_t_n,
     &                      Hb_x_np1,Hb_x_n,Hb_x_nm1,Hb_x_t_n,
     &                      Hb_y_np1,Hb_y_n,Hb_y_nm1,Hb_y_t_n,
     &                      phi1_np1,phi1_n,phi1_nm1,phi1_t_n,
     &                      L,phys_bdy,x,y,dt,chr,ex,Nx,Ny)
        implicit none
        integer Nx,Ny
        integer phys_bdy(4)
        real*8 dt,ex,L
        real*8 chr(Nx,Ny)
        real*8 gb_tt_np1(Nx,Ny),gb_tt_n(Nx,Ny),gb_tt_nm1(Nx,Ny)
        real*8 gb_tt_t_n(Nx,Ny)
        real*8 gb_tx_np1(Nx,Ny),gb_tx_n(Nx,Ny),gb_tx_nm1(Nx,Ny)
        real*8 gb_tx_t_n(Nx,Ny)
        real*8 gb_ty_np1(Nx,Ny),gb_ty_n(Nx,Ny),gb_ty_nm1(Nx,Ny)
        real*8 gb_ty_t_n(Nx,Ny)
        real*8 gb_xx_np1(Nx,Ny),gb_xx_n(Nx,Ny),gb_xx_nm1(Nx,Ny)
        real*8 gb_xx_t_n(Nx,Ny)
        real*8 gb_xy_np1(Nx,Ny),gb_xy_n(Nx,Ny),gb_xy_nm1(Nx,Ny)
        real*8 gb_xy_t_n(Nx,Ny)
        real*8 gb_yy_np1(Nx,Ny),gb_yy_n(Nx,Ny),gb_yy_nm1(Nx,Ny)
        real*8 gb_yy_t_n(Nx,Ny)
        real*8 psi_np1(Nx,Ny),psi_n(Nx,Ny),psi_nm1(Nx,Ny)
        real*8 psi_t_n(Nx,Ny)
        real*8 omega_np1(Nx,Ny),omega_n(Nx,Ny),omega_nm1(Nx,Ny)
        real*8 omega_t_n(Nx,Ny)
        real*8 fb_t_np1(Nx,Ny),fb_t_n(Nx,Ny),fb_t_nm1(Nx,Ny)
        real*8 fb_t_t_n(Nx,Ny)
        real*8 fb_x_np1(Nx,Ny),fb_x_n(Nx,Ny),fb_x_nm1(Nx,Ny)
        real*8 fb_x_t_n(Nx,Ny)
        real*8 fb_y_np1(Nx,Ny),fb_y_n(Nx,Ny),fb_y_nm1(Nx,Ny)
        real*8 fb_y_t_n(Nx,Ny)
        real*8 Hb_t_np1(Nx,Ny),Hb_t_n(Nx,Ny),Hb_t_nm1(Nx,Ny)
        real*8 Hb_t_t_n(Nx,Ny)
        real*8 Hb_x_np1(Nx,Ny),Hb_x_n(Nx,Ny),Hb_x_nm1(Nx,Ny)
        real*8 Hb_x_t_n(Nx,Ny)
        real*8 Hb_y_np1(Nx,Ny),Hb_y_n(Nx,Ny),Hb_y_nm1(Nx,Ny)
        real*8 Hb_y_t_n(Nx,Ny)
        real*8 phi1_np1(Nx,Ny),phi1_n(Nx,Ny),phi1_nm1(Nx,Ny)
        real*8 phi1_t_n(Nx,Ny)

        real*8 x(Nx),y(Ny)

        logical is_nan

        !--------------------------------------------------------------
        ! the following are first and second time derivatives of *n*
        ! level variables, and as these are the only derivatives we
        ! use we drop any _n identifier
        !--------------------------------------------------------------
        real*8 gb_tt_t,gb_tt_tt
        real*8 gb_tx_t,gb_tx_tt
        real*8 gb_ty_t,gb_ty_tt
        real*8 gb_xx_t,gb_xx_tt
        real*8 gb_xy_t,gb_xy_tt
        real*8 gb_yy_t,gb_yy_tt
        real*8 psi_t,psi_tt
        real*8 omega_t,omega_tt
        real*8 fb_t_t,fb_x_t,fb_y_t
        real*8 Hb_t_t,Hb_x_t,Hb_y_t
        real*8 phi1_t,phi1_tt

        real*8 h0_ll_tt(3,3)
        real*8 gA_tt,gB_tt
        real*8 phi10_tt

        real*8 x0,y0

        integer i,is,ie
        integer j,js,je
        integer a,b,c,d

        real*8 dV1_dphi10 ! NEED TO HAVE THIS AS AN INPUT ARGUMENT

        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 dx,dy

        logical ltrace
        parameter (ltrace=.false.)

        real*8 dimA,dimB

        real*8 grad_phi1_sq

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
        data i,is,ie/0,0,0/
        data j,js,je/0,0,0/

        data x0/0.0/

        data dx/0.0/

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

        !---------------------------------------------------------------

        dx=x(2)-x(1)
        dy=y(2)-y(1)

        ! set dimensions of S3 and S4 subspaces
        dimA=3d0 
        dimB=4d0

        do i=1,Nx
          do j=1,Ny
            phi1_nm1(i,j)=phi1_n(i,j)
          end do
        end do

        is=1
        ie=Nx
        js=1
        je=Ny
        if (phys_bdy(1).eq.1) is=2
        if (phys_bdy(2).eq.1) ie=Nx-1
        if (phys_bdy(3).eq.1) js=2
        if (phys_bdy(4).eq.1) je=Ny-1

        do i=is,ie
          do j=js,je
            x0=x(i)
            y0=y(j)

            if (chr(i,j).ne.ex) then

              !-----------------------------------------------------------
              ! some other initializion, which needs to be done before
              ! temporal derivatives are calculated
              !-----------------------------------------------------------

              ! computes tensors at point i,j 
              call tensor_init(
     &                gb_tt_n,gb_tt_n,gb_tt_n,
     &                gb_tx_n,gb_tx_n,gb_tx_n,
     &                gb_ty_n,gb_ty_n,gb_ty_n,
     &                gb_xx_n,gb_xx_n,gb_xx_n,
     &                gb_xy_n,gb_xy_n,gb_xy_n,
     &                gb_yy_n,gb_yy_n,gb_yy_n,
     &                psi_n,psi_n,psi_n,
     &                omega_n,omega_n,omega_n,
     &                fb_t_np1,fb_t_n,fb_t_nm1,
     &                fb_x_np1,fb_x_n,fb_x_nm1,
     &                fb_y_np1,fb_y_n,fb_y_nm1,
     &                Hb_t_n,Hb_t_n,Hb_t_n,
     &                Hb_x_n,Hb_x_n,Hb_x_n,
     &                Hb_y_n,Hb_y_n,Hb_y_n,
     &                phi1_n,phi1_n,phi1_n,
     &                g0_ll,g0_uu,g0_ll_x,g0_uu_x,g0_ll_xx,
     &                gads_ll,gads_uu,gads_ll_x,gads_uu_x,gads_ll_xx,
     &                gA,gB,gA_x,gB_x,gA_xx,gB_xx,
     &                gAads,gBads,gAads_x,gBads_x,gAads_xx,gBads_xx,
     &                sqrtdetg,sqrtdetg_x,
     &                sqrtdetgads,sqrtdetgads_x,
     &                h0_ll,h0_uu,h0_ll_x,h0_uu_x,h0_ll_xx,
     &                hA,hB,hAu,hBu,hA_x,hB_x,hA_xx,hB_xx,
     &                sqrtdeth,sqrtdeth_x,
     &                H0_l,H0_l_x,Hads_l,Hads_l_x,A_l,A_l_x,
     &                gamma_ull,gamma_ull_x,
     &                riemann_ulll,ricci_ll,ricci_lu,ricci,
     &                f1_l,f1_l_x,f2_ll,f2_ll_x,
     &                f1ads_l,f1ads_l_x,f2ads_ll,f2ads_ll_x,
     &                h1_l,h1_l_x,h2_ll,h2_ll_x,
     &                phi10_x,phi10_xx,
     &                gammagg,gammahh,gammagh,gammahg,dlll,cuuuu,
     &                x,y,dt,chr,L,ex,Nx,Ny,i,j)

              ! initial first time derivatives; gb_ii_t_n,Hb_i_t_n,phi1_t_n were set in AdS5xS5_free_data()

              ! need this in h0_ll_tt,psi_tt,omega_tt,phi10_tt calculations
              phi10_x(1)    =phi1_t_n(i,j)*(1-x0**2)**3   
              h0_ll_x(1,1,1)=gb_tt_t_n(i,j)*(1-x0**2)  
              h0_ll_x(1,2,1)=gb_tx_t_n(i,j)*(1-x0**2)
              h0_ll_x(1,3,1)=gb_ty_t_n(i,j)*(1-x0**2)**2
              h0_ll_x(2,2,1)=gb_xx_t_n(i,j)*(1-x0**2)  
              h0_ll_x(2,3,1)=gb_xy_t_n(i,j)*(1-x0**2)**2
              h0_ll_x(3,3,1)=gb_yy_t_n(i,j)*(1-x0**2)**3
              A_l_x(1,1)    =Hb_t_t_n(i,j)*(1-x0**2)**2
              A_l_x(2,1)    =Hb_x_t_n(i,j)*(1-x0**2)**2
              A_l_x(3,1)    =Hb_y_t_n(i,j)*(1-x0**2)**2
              f1_l_x(1,1)   =fb_t_t_n(i,j)*(1-x0**2)**3
              f1_l_x(2,1)   =fb_x_t_n(i,j)*(1-x0**2)**3
              f1_l_x(3,1)   =fb_y_t_n(i,j)*(1-x0**2)**3
              gA_x(1)       =psi_t_n(i,j)*(1-x0**2)*x0**2
              gB_x(1)       =omega_t_n(i,j)*(1-x0**2)**3

              ! need this in gb_ii_nm1/np1,Hb_i_nm1/np1,phi1_nm1/np1 updates
              phi1_t =phi1_t_n(i,j)                
              gb_tt_t=gb_tt_t_n(i,j) 
              gb_tx_t=gb_tx_t_n(i,j) 
              gb_ty_t=gb_ty_t_n(i,j) 
              gb_xx_t=gb_xx_t_n(i,j) 
              gb_xy_t=gb_xy_t_n(i,j) 
              gb_yy_t=gb_yy_t_n(i,j) 
              psi_t=psi_t_n(i,j) 
              omega_t=omega_t_n(i,j) 
              fb_t_t =fb_t_t_n(i,j)
              fb_x_t =fb_x_t_n(i,j)
              fb_y_t =fb_y_t_n(i,j)
              Hb_t_t =Hb_t_t_n(i,j)
              Hb_x_t =Hb_x_t_n(i,j)
              Hb_y_t =Hb_y_t_n(i,j)

              grad_phi1_sq=phi10_x(1)*phi10_x(1)*g0_uu(1,1)+
     &                     phi10_x(2)*phi10_x(2)*g0_uu(2,2)+
     &                     phi10_x(3)*phi10_x(3)*g0_uu(3,3)+
     &                  2*(phi10_x(1)*phi10_x(2)*g0_uu(1,2)+
     &                     phi10_x(1)*phi10_x(3)*g0_uu(1,3)+
     &                     phi10_x(2)*phi10_x(3)*g0_uu(2,3))

              ! 0 = efe_ab
              do a=1,3
                do b=a,3
                  h0_ll_tt(a,b)=2/g0_uu(1,1)*
     &                ( 
     &                           -0.5d0*(                             
     &                              h0_uu(2,2)*h0_ll_xx(a,b,2,2)+
     &                              h0_uu(3,3)*h0_ll_xx(a,b,3,3)+
     &                           2*(h0_uu(1,2)*h0_ll_xx(a,b,1,2)+
     &                              h0_uu(1,3)*h0_ll_xx(a,b,1,3)+
     &                              h0_uu(2,3)*h0_ll_xx(a,b,2,3))
     &                           +
     &                              gads_uu(2,2)*h0_ll_xx(a,b,2,2)+
     &                              gads_uu(3,3)*h0_ll_xx(a,b,3,3)+
     &                           2*(gads_uu(1,2)*h0_ll_xx(a,b,1,2)+
     &                              gads_uu(1,3)*h0_ll_xx(a,b,1,3)+
     &                              gads_uu(2,3)*h0_ll_xx(a,b,2,3))
     &                           +
     &                              h0_uu(1,1)*gads_ll_xx(a,b,1,1)+
     &                              h0_uu(2,2)*gads_ll_xx(a,b,2,2)+
     &                              h0_uu(3,3)*gads_ll_xx(a,b,3,3)+
     &                           2*(h0_uu(1,2)*gads_ll_xx(a,b,1,2)+
     &                              h0_uu(1,3)*gads_ll_xx(a,b,1,3)+
     &                              h0_uu(2,3)*gads_ll_xx(a,b,2,3))
     &                           +
     &                              gads_uu(1,1)*gads_ll_xx(a,b,1,1)+
     &                              gads_uu(2,2)*gads_ll_xx(a,b,2,2)+
     &                              gads_uu(3,3)*gads_ll_xx(a,b,3,3)+
     &                           2*(gads_uu(1,2)*gads_ll_xx(a,b,1,2)+
     &                              gads_uu(1,3)*gads_ll_xx(a,b,1,3)+
     &                              gads_uu(2,3)*gads_ll_xx(a,b,2,3))
     &                                  )
     &
     &                           -0.5d0*(
     &                              g0_uu_x(1,1,a)* g0_ll_x(b,1,1) +
     &                              g0_uu_x(1,2,a)*(g0_ll_x(b,1,2) +
     &                                              g0_ll_x(b,2,1))+
     &                              g0_uu_x(1,3,a)*(g0_ll_x(b,1,3) +
     &                                              g0_ll_x(b,3,1))+
     &                              g0_uu_x(2,2,a)* g0_ll_x(b,2,2) +
     &                              g0_uu_x(2,3,a)*(g0_ll_x(b,2,3) +
     &                                              g0_ll_x(b,3,2))+
     &                              g0_uu_x(3,3,a)* g0_ll_x(b,3,3) 
     &                                )
     &
     &                           -0.5d0*(
     &                              g0_uu_x(1,1,b)* g0_ll_x(a,1,1) +
     &                              g0_uu_x(1,2,b)*(g0_ll_x(a,1,2) +
     &                                              g0_ll_x(a,2,1))+
     &                              g0_uu_x(1,3,b)*(g0_ll_x(a,1,3) +
     &                                              g0_ll_x(a,3,1))+
     &                              g0_uu_x(2,2,b)* g0_ll_x(a,2,2) +
     &                              g0_uu_x(2,3,b)*(g0_ll_x(a,2,3) +
     &                                              g0_ll_x(a,3,2))+
     &                              g0_uu_x(3,3,b)* g0_ll_x(a,3,3) 
     &                                )
     &
     &                           -0.5d0*H0_l_x(a,b)
     &
     &                           -0.5d0*H0_l_x(b,a)
     &
     &                               +(
     &                              H0_l(1)*gamma_ull(1,a,b)+
     &                              H0_l(2)*gamma_ull(2,a,b)+
     &                              H0_l(3)*gamma_ull(3,a,b)
     &                                )
     &
     &                               -(
     &                              gamma_ull(1,1,b)*gamma_ull(1,1,a)+
     &                              gamma_ull(1,2,b)*gamma_ull(2,1,a)+
     &                              gamma_ull(1,3,b)*gamma_ull(3,1,a)+
     &                              gamma_ull(2,1,b)*gamma_ull(1,2,a)+
     &                              gamma_ull(2,2,b)*gamma_ull(2,2,a)+
     &                              gamma_ull(2,3,b)*gamma_ull(3,2,a)+
     &                              gamma_ull(3,1,b)*gamma_ull(1,3,a)+
     &                              gamma_ull(3,2,b)*gamma_ull(2,3,a)+
     &                              gamma_ull(3,3,b)*gamma_ull(3,3,a)
     &                                )
     &
     &                           -dimA*( (gA_xx(a,b) ! OLD TERMS
     &                                   -gamma_ull(1,b,a)*gA_x(1)
     &                                   -gamma_ull(2,b,a)*gA_x(2)
     &                                   -gamma_ull(3,b,a)*gA_x(3)
     &                                   )/(2*gA)
     &                                 - (gA_x(a)*gA_x(b))/(4*gA**2) )
     &                           -dimB*( (gB_xx(a,b)
     &                                   -gamma_ull(1,b,a)*gB_x(1)
     &                                   -gamma_ull(2,b,a)*gB_x(2)
     &                                   -gamma_ull(3,b,a)*gB_x(3)
     &                                   )/(2*gB)
     &                                 - (gB_x(a)*gB_x(b))/(4*gB**2) )
     &
     &                           -(f1_l(a)*f1_l(b)
     &                            +g0_uu(1,1)*f2_ll(a,1)*f2_ll(b,1)
     &                            +g0_uu(1,2)*f2_ll(a,1)*f2_ll(b,2)
     &                            +g0_uu(1,3)*f2_ll(a,1)*f2_ll(b,3)
     &                            +g0_uu(2,1)*f2_ll(a,2)*f2_ll(b,1)
     &                            +g0_uu(2,2)*f2_ll(a,2)*f2_ll(b,2)
     &                            +g0_uu(2,3)*f2_ll(a,2)*f2_ll(b,3)
     &                            +g0_uu(3,1)*f2_ll(a,3)*f2_ll(b,1)
     &                            +g0_uu(3,2)*f2_ll(a,3)*f2_ll(b,2)
     &                            +g0_uu(3,3)*f2_ll(a,3)*f2_ll(b,3)
     &                            )/4
     & 
!     &                           -dimA*(gA_x(a)*gA_x(b))/(4*gA**2) !NEW TERMS
!     &                           -dimB*(gB_x(a)*gB_x(b))/(4*gB**2)
!     &
!     &                           -(2*phi10_x(a)*phi10_x(b)
!     &                           -g0_ll(a,b)*grad_phi1_sq
!     &                            )/4/gB**4
     &                )
                end do
              end do          

              ! 0 = afe
              gA_tt=2/g0_uu(1,1)*
     &                (
     &         -g0_uu(1,2)*gA_xx(1,2)/2
     &         -g0_uu(1,3)*gA_xx(1,3)/2
     &         -g0_uu(2,1)*gA_xx(2,1)/2
     &         -g0_uu(2,2)*gA_xx(2,2)/2
     &         -g0_uu(2,3)*gA_xx(2,3)/2
     &         -g0_uu(3,1)*gA_xx(3,1)/2
     &         -g0_uu(3,2)*gA_xx(3,2)/2
     &         -g0_uu(3,3)*gA_xx(3,3)/2
     &         -g0_uu(1,1)*H0_l(1)*gA_x(1)/2
     &         -g0_uu(1,2)*H0_l(1)*gA_x(2)/2
     &         -g0_uu(1,3)*H0_l(1)*gA_x(3)/2
     &         -g0_uu(2,1)*H0_l(2)*gA_x(1)/2
     &         -g0_uu(2,2)*H0_l(2)*gA_x(2)/2
     &         -g0_uu(2,3)*H0_l(2)*gA_x(3)/2
     &         -g0_uu(3,1)*H0_l(3)*gA_x(1)/2
     &         -g0_uu(3,2)*H0_l(3)*gA_x(2)/2
     &         -g0_uu(3,3)*H0_l(3)*gA_x(3)/2
     &         +g0_uu(1,1)*(gA_x(1)*gA_x(1))/(2*gA)
     &         +g0_uu(1,2)*(gA_x(1)*gA_x(2))/(2*gA)
     &         +g0_uu(1,3)*(gA_x(1)*gA_x(3))/(2*gA)
     &         +g0_uu(2,1)*(gA_x(2)*gA_x(1))/(2*gA)
     &         +g0_uu(2,2)*(gA_x(2)*gA_x(2))/(2*gA)
     &         +g0_uu(2,3)*(gA_x(2)*gA_x(3))/(2*gA)
     &         +g0_uu(3,1)*(gA_x(3)*gA_x(1))/(2*gA)
     &         +g0_uu(3,2)*(gA_x(3)*gA_x(2))/(2*gA)
     &         +g0_uu(3,3)*(gA_x(3)*gA_x(3))/(2*gA)
     &
     &         +(dimA-1d0)
     &
     &         -g0_uu(1,1)*(gA_x(1)*gA_x(1))/(4*gA)*(dimA) !OLD TERMS
     &         -g0_uu(1,2)*(gA_x(1)*gA_x(2))/(4*gA)*(dimA)
     &         -g0_uu(1,3)*(gA_x(1)*gA_x(3))/(4*gA)*(dimA)
     &         -g0_uu(2,1)*(gA_x(2)*gA_x(1))/(4*gA)*(dimA)
     &         -g0_uu(2,2)*(gA_x(2)*gA_x(2))/(4*gA)*(dimA)
     &         -g0_uu(2,3)*(gA_x(2)*gA_x(3))/(4*gA)*(dimA)
     &         -g0_uu(3,1)*(gA_x(3)*gA_x(1))/(4*gA)*(dimA)
     &         -g0_uu(3,2)*(gA_x(3)*gA_x(2))/(4*gA)*(dimA)
     &         -g0_uu(3,3)*(gA_x(3)*gA_x(3))/(4*gA)*(dimA)
     &         -g0_uu(1,1)*(gA_x(1)*gB_x(1))/(4*gB)*dimB
     &         -g0_uu(1,2)*(gA_x(1)*gB_x(2))/(4*gB)*dimB
     &         -g0_uu(1,3)*(gA_x(1)*gB_x(3))/(4*gB)*dimB
     &         -g0_uu(2,1)*(gA_x(2)*gB_x(1))/(4*gB)*dimB
     &         -g0_uu(2,2)*(gA_x(2)*gB_x(2))/(4*gB)*dimB
     &         -g0_uu(2,3)*(gA_x(2)*gB_x(3))/(4*gB)*dimB
     &         -g0_uu(3,1)*(gA_x(3)*gB_x(1))/(4*gB)*dimB
     &         -g0_uu(3,2)*(gA_x(3)*gB_x(2))/(4*gB)*dimB
     &         -g0_uu(3,3)*(gA_x(3)*gB_x(3))/(4*gB)*dimB
     &         -gA*g0_uu(1,1)*g0_uu(1,1)*f2_ll(1,1)*f2_ll(1,1)/8
     &         -gA*g0_uu(1,2)*g0_uu(1,1)*f2_ll(1,1)*f2_ll(2,1)/8
     &         -gA*g0_uu(1,3)*g0_uu(1,1)*f2_ll(1,1)*f2_ll(3,1)/8
     &         -gA*g0_uu(2,1)*g0_uu(1,1)*f2_ll(2,1)*f2_ll(1,1)/8
     &         -gA*g0_uu(2,2)*g0_uu(1,1)*f2_ll(2,1)*f2_ll(2,1)/8
     &         -gA*g0_uu(2,3)*g0_uu(1,1)*f2_ll(2,1)*f2_ll(3,1)/8
     &         -gA*g0_uu(3,1)*g0_uu(1,1)*f2_ll(3,1)*f2_ll(1,1)/8
     &         -gA*g0_uu(3,2)*g0_uu(1,1)*f2_ll(3,1)*f2_ll(2,1)/8
     &         -gA*g0_uu(3,3)*g0_uu(1,1)*f2_ll(3,1)*f2_ll(3,1)/8
     &         -gA*g0_uu(1,1)*g0_uu(1,2)*f2_ll(1,1)*f2_ll(1,2)/8
     &         -gA*g0_uu(1,2)*g0_uu(1,2)*f2_ll(1,1)*f2_ll(2,2)/8
     &         -gA*g0_uu(1,3)*g0_uu(1,2)*f2_ll(1,1)*f2_ll(3,2)/8
     &         -gA*g0_uu(2,1)*g0_uu(1,2)*f2_ll(2,1)*f2_ll(1,2)/8
     &         -gA*g0_uu(2,2)*g0_uu(1,2)*f2_ll(2,1)*f2_ll(2,2)/8
     &         -gA*g0_uu(2,3)*g0_uu(1,2)*f2_ll(2,1)*f2_ll(3,2)/8
     &         -gA*g0_uu(3,1)*g0_uu(1,2)*f2_ll(3,1)*f2_ll(1,2)/8
     &         -gA*g0_uu(3,2)*g0_uu(1,2)*f2_ll(3,1)*f2_ll(2,2)/8
     &         -gA*g0_uu(3,3)*g0_uu(1,2)*f2_ll(3,1)*f2_ll(3,2)/8
     &         -gA*g0_uu(1,1)*g0_uu(1,3)*f2_ll(1,1)*f2_ll(1,3)/8
     &         -gA*g0_uu(1,2)*g0_uu(1,3)*f2_ll(1,1)*f2_ll(2,3)/8
     &         -gA*g0_uu(1,3)*g0_uu(1,3)*f2_ll(1,1)*f2_ll(3,3)/8
     &         -gA*g0_uu(2,1)*g0_uu(1,3)*f2_ll(2,1)*f2_ll(1,3)/8
     &         -gA*g0_uu(2,2)*g0_uu(1,3)*f2_ll(2,1)*f2_ll(2,3)/8
     &         -gA*g0_uu(2,3)*g0_uu(1,3)*f2_ll(2,1)*f2_ll(3,3)/8
     &         -gA*g0_uu(3,1)*g0_uu(1,3)*f2_ll(3,1)*f2_ll(1,3)/8
     &         -gA*g0_uu(3,2)*g0_uu(1,3)*f2_ll(3,1)*f2_ll(2,3)/8
     &         -gA*g0_uu(3,3)*g0_uu(1,3)*f2_ll(3,1)*f2_ll(3,3)/8
     &         -gA*g0_uu(1,1)*g0_uu(2,1)*f2_ll(1,2)*f2_ll(1,1)/8
     &         -gA*g0_uu(1,2)*g0_uu(2,1)*f2_ll(1,2)*f2_ll(2,1)/8
     &         -gA*g0_uu(1,3)*g0_uu(2,1)*f2_ll(1,2)*f2_ll(3,1)/8
     &         -gA*g0_uu(2,1)*g0_uu(2,1)*f2_ll(2,2)*f2_ll(1,1)/8
     &         -gA*g0_uu(2,2)*g0_uu(2,1)*f2_ll(2,2)*f2_ll(2,1)/8
     &         -gA*g0_uu(2,3)*g0_uu(2,1)*f2_ll(2,2)*f2_ll(3,1)/8
     &         -gA*g0_uu(3,1)*g0_uu(2,1)*f2_ll(3,2)*f2_ll(1,1)/8
     &         -gA*g0_uu(3,2)*g0_uu(2,1)*f2_ll(3,2)*f2_ll(2,1)/8
     &         -gA*g0_uu(3,3)*g0_uu(2,1)*f2_ll(3,2)*f2_ll(3,1)/8
     &         -gA*g0_uu(1,1)*g0_uu(2,2)*f2_ll(1,2)*f2_ll(1,2)/8
     &         -gA*g0_uu(1,2)*g0_uu(2,2)*f2_ll(1,2)*f2_ll(2,2)/8
     &         -gA*g0_uu(1,3)*g0_uu(2,2)*f2_ll(1,2)*f2_ll(3,2)/8
     &         -gA*g0_uu(2,1)*g0_uu(2,2)*f2_ll(2,2)*f2_ll(1,2)/8
     &         -gA*g0_uu(2,2)*g0_uu(2,2)*f2_ll(2,2)*f2_ll(2,2)/8
     &         -gA*g0_uu(2,3)*g0_uu(2,2)*f2_ll(2,2)*f2_ll(3,2)/8
     &         -gA*g0_uu(3,1)*g0_uu(2,2)*f2_ll(3,2)*f2_ll(1,2)/8
     &         -gA*g0_uu(3,2)*g0_uu(2,2)*f2_ll(3,2)*f2_ll(2,2)/8
     &         -gA*g0_uu(3,3)*g0_uu(2,2)*f2_ll(3,2)*f2_ll(3,2)/8
     &         -gA*g0_uu(1,1)*g0_uu(2,3)*f2_ll(1,2)*f2_ll(1,3)/8
     &         -gA*g0_uu(1,2)*g0_uu(2,3)*f2_ll(1,2)*f2_ll(2,3)/8
     &         -gA*g0_uu(1,3)*g0_uu(2,3)*f2_ll(1,2)*f2_ll(3,3)/8
     &         -gA*g0_uu(2,1)*g0_uu(2,3)*f2_ll(2,2)*f2_ll(1,3)/8
     &         -gA*g0_uu(2,2)*g0_uu(2,3)*f2_ll(2,2)*f2_ll(2,3)/8
     &         -gA*g0_uu(2,3)*g0_uu(2,3)*f2_ll(2,2)*f2_ll(3,3)/8
     &         -gA*g0_uu(3,1)*g0_uu(2,3)*f2_ll(3,2)*f2_ll(1,3)/8
     &         -gA*g0_uu(3,2)*g0_uu(2,3)*f2_ll(3,2)*f2_ll(2,3)/8
     &         -gA*g0_uu(3,3)*g0_uu(2,3)*f2_ll(3,2)*f2_ll(3,3)/8
     &         -gA*g0_uu(1,1)*g0_uu(3,1)*f2_ll(1,3)*f2_ll(1,1)/8
     &         -gA*g0_uu(1,2)*g0_uu(3,1)*f2_ll(1,3)*f2_ll(2,1)/8
     &         -gA*g0_uu(1,3)*g0_uu(3,1)*f2_ll(1,3)*f2_ll(3,1)/8
     &         -gA*g0_uu(2,1)*g0_uu(3,1)*f2_ll(2,3)*f2_ll(1,1)/8
     &         -gA*g0_uu(2,2)*g0_uu(3,1)*f2_ll(2,3)*f2_ll(2,1)/8
     &         -gA*g0_uu(2,3)*g0_uu(3,1)*f2_ll(2,3)*f2_ll(3,1)/8
     &         -gA*g0_uu(3,1)*g0_uu(3,1)*f2_ll(3,3)*f2_ll(1,1)/8
     &         -gA*g0_uu(3,2)*g0_uu(3,1)*f2_ll(3,3)*f2_ll(2,1)/8
     &         -gA*g0_uu(3,3)*g0_uu(3,1)*f2_ll(3,3)*f2_ll(3,1)/8
     &         -gA*g0_uu(1,1)*g0_uu(3,2)*f2_ll(1,3)*f2_ll(1,2)/8
     &         -gA*g0_uu(1,2)*g0_uu(3,2)*f2_ll(1,3)*f2_ll(2,2)/8
     &         -gA*g0_uu(1,3)*g0_uu(3,2)*f2_ll(1,3)*f2_ll(3,2)/8
     &         -gA*g0_uu(2,1)*g0_uu(3,2)*f2_ll(2,3)*f2_ll(1,2)/8
     &         -gA*g0_uu(2,2)*g0_uu(3,2)*f2_ll(2,3)*f2_ll(2,2)/8
     &         -gA*g0_uu(2,3)*g0_uu(3,2)*f2_ll(2,3)*f2_ll(3,2)/8
     &         -gA*g0_uu(3,1)*g0_uu(3,2)*f2_ll(3,3)*f2_ll(1,2)/8
     &         -gA*g0_uu(3,2)*g0_uu(3,2)*f2_ll(3,3)*f2_ll(2,2)/8
     &         -gA*g0_uu(3,3)*g0_uu(3,2)*f2_ll(3,3)*f2_ll(3,2)/8
     &         -gA*g0_uu(1,1)*g0_uu(3,3)*f2_ll(1,3)*f2_ll(1,3)/8
     &         -gA*g0_uu(1,2)*g0_uu(3,3)*f2_ll(1,3)*f2_ll(2,3)/8
     &         -gA*g0_uu(1,3)*g0_uu(3,3)*f2_ll(1,3)*f2_ll(3,3)/8
     &         -gA*g0_uu(2,1)*g0_uu(3,3)*f2_ll(2,3)*f2_ll(1,3)/8
     &         -gA*g0_uu(2,2)*g0_uu(3,3)*f2_ll(2,3)*f2_ll(2,3)/8
     &         -gA*g0_uu(2,3)*g0_uu(3,3)*f2_ll(2,3)*f2_ll(3,3)/8
     &         -gA*g0_uu(3,1)*g0_uu(3,3)*f2_ll(3,3)*f2_ll(1,3)/8
     &         -gA*g0_uu(3,2)*g0_uu(3,3)*f2_ll(3,3)*f2_ll(2,3)/8
     &         -gA*g0_uu(3,3)*g0_uu(3,3)*f2_ll(3,3)*f2_ll(3,3)/8
     & 
!     &         +gA*grad_phi1_sq/4/gB**4 ! NEW TERMS
     &                )

!uncommenting Analytic removal terms below would Analytically remove pure ads5xs5 terms but 
!neglecting the h^{-1}*g contribution, and so implicitly assumes that h^{-1}=0
!i.e. that g0_uu(3,3)==gads_uu(3,3) for the initial data being used,
!this is Analytically true for bh initial data for which g0_uu(3,3)=1/PI^2=gads_uu(3,3) 
!but the g0_uu(3,3) computed via tensor_init() has imperfect cancellations which 

              ! 0 = bfe
              gB_tt=2/g0_uu(1,1)*
     &                (
     &         -g0_uu(1,2)*gB_xx(1,2)/2
     &         -g0_uu(1,3)*gB_xx(1,3)/2
     &         -g0_uu(2,1)*gB_xx(2,1)/2
     &         -g0_uu(2,2)*gB_xx(2,2)/2
     &         -g0_uu(2,3)*gB_xx(2,3)/2
     &         -g0_uu(3,1)*gB_xx(3,1)/2
     &         -g0_uu(3,2)*gB_xx(3,2)/2
     &         -g0_uu(3,3)*gB_xx(3,3)/2
     &         -g0_uu(1,1)*H0_l(1)*gB_x(1)/2
     &         -g0_uu(1,2)*H0_l(1)*gB_x(2)/2
     &         -g0_uu(1,3)*H0_l(1)*gB_x(3)/2
     &         -g0_uu(2,1)*H0_l(2)*gB_x(1)/2
     &         -g0_uu(2,2)*H0_l(2)*gB_x(2)/2
     &         -g0_uu(2,3)*H0_l(2)*gB_x(3)/2
     &         -g0_uu(3,1)*H0_l(3)*gB_x(1)/2
     &         -g0_uu(3,2)*H0_l(3)*gB_x(2)/2
     &         -g0_uu(3,3)*H0_l(3)*gB_x(3)/2
     &         +g0_uu(1,1)*(gB_x(1)*gB_x(1))/(2*gB)
     &         +g0_uu(1,2)*(gB_x(1)*gB_x(2))/(2*gB)
     &         +g0_uu(1,3)*(gB_x(1)*gB_x(3))/(2*gB)
     &         +g0_uu(2,1)*(gB_x(2)*gB_x(1))/(2*gB)
     &         +g0_uu(2,2)*(gB_x(2)*gB_x(2))/(2*gB)
     &         +g0_uu(2,3)*(gB_x(2)*gB_x(3))/(2*gB)
     &         +g0_uu(3,1)*(gB_x(3)*gB_x(1))/(2*gB)
     &         +g0_uu(3,2)*(gB_x(3)*gB_x(2))/(2*gB)
     &         +g0_uu(3,3)*(gB_x(3)*gB_x(3))/(2*gB)
     &
     &         +(dimB-1d0) !analytic removal
     &
     &         -g0_uu(1,1)*(gB_x(1)*gB_x(1))/(4*gB)*(dimB) !OLD TERMS
     &         -g0_uu(1,2)*(gB_x(1)*gB_x(2))/(4*gB)*(dimB)
     &         -g0_uu(1,3)*(gB_x(1)*gB_x(3))/(4*gB)*(dimB)
     &         -g0_uu(2,1)*(gB_x(2)*gB_x(1))/(4*gB)*(dimB)
     &         -g0_uu(2,2)*(gB_x(2)*gB_x(2))/(4*gB)*(dimB)
     &         -g0_uu(2,3)*(gB_x(2)*gB_x(3))/(4*gB)*(dimB)
     &         -g0_uu(3,1)*(gB_x(3)*gB_x(1))/(4*gB)*(dimB)
     &         -g0_uu(3,2)*(gB_x(3)*gB_x(2))/(4*gB)*(dimB)
     &         -g0_uu(3,3)*(gB_x(3)*gB_x(3))/(4*gB)*(dimB)
     &         -g0_uu(1,1)*(gA_x(1)*gB_x(1))/(4*gA)*dimA
     &         -g0_uu(1,2)*(gA_x(1)*gB_x(2))/(4*gA)*dimA
     &         -g0_uu(1,3)*(gA_x(1)*gB_x(3))/(4*gA)*dimA
     &         -g0_uu(2,1)*(gA_x(2)*gB_x(1))/(4*gA)*dimA
     &         -g0_uu(2,2)*(gA_x(2)*gB_x(2))/(4*gA)*dimA
     &         -g0_uu(2,3)*(gA_x(2)*gB_x(3))/(4*gA)*dimA
     &         -g0_uu(3,1)*(gA_x(3)*gB_x(1))/(4*gA)*dimA
     &         -g0_uu(3,2)*(gA_x(3)*gB_x(2))/(4*gA)*dimA
     &         -g0_uu(3,3)*(gA_x(3)*gB_x(3))/(4*gA)*dimA
     &         -gB*g0_uu(1,1)*f1_l(1)*f1_l(1)/4
     &         -gB*g0_uu(1,2)*f1_l(1)*f1_l(2)/4
     &         -gB*g0_uu(1,3)*f1_l(1)*f1_l(3)/4
     &         -gB*g0_uu(2,1)*f1_l(2)*f1_l(1)/4
     &         -gB*g0_uu(2,2)*f1_l(2)*f1_l(2)/4
     &         -gB*g0_uu(2,3)*f1_l(2)*f1_l(3)/4
     &         -gB*g0_uu(3,1)*f1_l(3)*f1_l(1)/4
     &         -gB*g0_uu(3,2)*f1_l(3)*f1_l(2)/4
     &         -gB*g0_uu(3,3)*f1_l(3)*f1_l(3)/4
     & 
!     &         -grad_phi1_sq/4/gB**3 !NEW TERMS
     &                )

              ! 0 = g^ab phi1,ab - g^ab gamma^c_ab phi1,c 
              phi10_tt=-1/g0_uu(1,1)
     &                *(    
     &                      phi10_xx(2,2)*g0_uu(2,2)+
     &                      phi10_xx(3,3)*g0_uu(3,3)+
     &                   2*(phi10_xx(1,2)*g0_uu(1,2)+
     &                      phi10_xx(1,3)*g0_uu(1,3)+
     &                      phi10_xx(2,3)*g0_uu(2,3))
     &                  -
     &                      phi10_x(1)*( gamma_ull(1,1,1)*g0_uu(1,1)+
     &                                   gamma_ull(1,2,2)*g0_uu(2,2)+
     &                                   gamma_ull(1,3,3)*g0_uu(3,3)+
     &                                2*(gamma_ull(1,1,2)*g0_uu(1,2)+
     &                                   gamma_ull(1,1,3)*g0_uu(1,3)+
     &                                   gamma_ull(1,2,3)*g0_uu(2,3)) )
     &                  -
     &                      phi10_x(2)*( gamma_ull(2,1,1)*g0_uu(1,1)+
     &                                   gamma_ull(2,2,2)*g0_uu(2,2)+
     &                                   gamma_ull(2,3,3)*g0_uu(3,3)+
     &                                2*(gamma_ull(2,1,2)*g0_uu(1,2)+
     &                                   gamma_ull(2,1,3)*g0_uu(1,3)+
     &                                   gamma_ull(2,2,3)*g0_uu(2,3)) )
     &                  -
     &                      phi10_x(3)*( gamma_ull(3,1,1)*g0_uu(1,1)+
     &                                   gamma_ull(3,2,2)*g0_uu(2,2)+
     &                                   gamma_ull(3,3,3)*g0_uu(3,3)+
     &                                2*(gamma_ull(3,1,2)*g0_uu(1,2)+
     &                                   gamma_ull(3,1,3)*g0_uu(1,3)+
     &                                   gamma_ull(3,2,3)*g0_uu(2,3)) )
     &                  +
     &                     1.5d0/gA*(
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
     &                     2.0d0/gB*(
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
     &                                                                 )

              if (is_nan(h0_ll_tt(1,1)).or.is_nan(h0_ll_tt(1,2)) 
     &        .or.is_nan(h0_ll_tt(1,3)).or.is_nan(h0_ll_tt(2,2))
     &        .or.is_nan(h0_ll_tt(2,3)).or.is_nan(h0_ll_tt(3,3))
     &        .or.is_nan(gA_tt).or.is_nan(gB_tt)) then
                write(*,*) 'h0_ll_tt(1,1)=',h0_ll_tt(1,1)
                write(*,*) 'h0_ll_tt(1,2)=',h0_ll_tt(1,2)
                write(*,*) 'h0_ll_tt(1,3)=',h0_ll_tt(1,3)
                write(*,*) 'h0_ll_tt(2,2)=',h0_ll_tt(2,2)
                write(*,*) 'h0_ll_tt(2,3)=',h0_ll_tt(2,3)
                write(*,*) 'h0_ll_tt(3,3)=',h0_ll_tt(3,3)
                write(*,*) 'gA_tt=',gA_tt
                write(*,*) 'gB_tt=',gB_tt
                stop
              end if

              ! initial second time derivatives
              gb_tt_tt=h0_ll_tt(1,1)/(1-x0**2) 
              gb_tx_tt=h0_ll_tt(1,2)/(1-x0**2)
              gb_ty_tt=h0_ll_tt(1,3)/(1-x0**2)**2
              gb_xx_tt=h0_ll_tt(2,2)/(1-x0**2) 
              gb_xy_tt=h0_ll_tt(2,3)/(1-x0**2)**2
              gb_yy_tt=h0_ll_tt(3,3)/(1-x0**2)**3
              psi_tt=gA_tt/(1-x0**2)/x0**2
              omega_tt=gB_tt/(1-x0**2)**3/sin(PI*y0/L)**2
              phi1_tt =phi10_tt/(1-x0**2)**3

              ! initialize past time level by O(h^3) expansion
              gb_tt_nm1(i,j)=gb_tt_n(i,j) - gb_tt_t*dt
     &                     + gb_tt_tt*dt**2/2
              gb_tx_nm1(i,j)=gb_tx_n(i,j) - gb_tx_t*dt
     &                     + gb_tx_tt*dt**2/2
              gb_ty_nm1(i,j)=gb_ty_n(i,j) - gb_ty_t*dt
     &                     + gb_ty_tt*dt**2/2
              gb_xx_nm1(i,j)=gb_xx_n(i,j) - gb_xx_t*dt
     &                     + gb_xx_tt*dt**2/2
              gb_xy_nm1(i,j)=gb_xy_n(i,j) - gb_xy_t*dt
     &                     + gb_xy_tt*dt**2/2
              gb_yy_nm1(i,j)=gb_yy_n(i,j) - gb_yy_t*dt
     &                     + gb_yy_tt*dt**2/2
              psi_nm1(i,j)  =psi_n(i,j) - psi_t*dt  
     &                     + psi_tt*dt**2/2
              omega_nm1(i,j)=omega_n(i,j) - omega_t*dt  
     &                     + omega_tt*dt**2/2
              fb_t_nm1(i,j) =fb_t_n(i,j) - fb_t_t*dt
              fb_x_nm1(i,j) =fb_x_n(i,j) - fb_x_t*dt
              fb_y_nm1(i,j) =fb_y_n(i,j) - fb_y_t*dt
              Hb_t_nm1(i,j) =Hb_t_n(i,j) - Hb_t_t*dt
              Hb_x_nm1(i,j) =Hb_x_n(i,j) - Hb_x_t*dt
              Hb_y_nm1(i,j) =Hb_y_n(i,j) - Hb_y_t*dt
              phi1_nm1(i,j) =phi1_n(i,j) - phi1_t*dt
     &                     + phi1_tt*dt**2/2
        
              ! initialize future time level by O(h^3) expansion
              gb_tt_np1(i,j)=gb_tt_n(i,j) + gb_tt_t*dt
     &                     + gb_tt_tt*dt**2/2
              gb_tx_np1(i,j)=gb_tx_n(i,j) + gb_tx_t*dt
     &                     + gb_tx_tt*dt**2/2
              gb_ty_np1(i,j)=gb_ty_n(i,j) + gb_ty_t*dt
     &                     + gb_ty_tt*dt**2/2
              gb_xx_np1(i,j)=gb_xx_n(i,j) + gb_xx_t*dt
     &                     + gb_xx_tt*dt**2/2
              gb_xy_np1(i,j)=gb_xy_n(i,j) + gb_xy_t*dt
     &                     + gb_xy_tt*dt**2/2
              gb_yy_np1(i,j)=gb_yy_n(i,j) + gb_yy_t*dt
     &                     + gb_yy_tt*dt**2/2
              psi_np1(i,j)  =psi_n(i,j) + psi_t*dt  
     &                     + psi_tt*dt**2/2
              omega_np1(i,j)=omega_n(i,j) + omega_t*dt  
     &                     + omega_tt*dt**2/2
              Hb_t_np1(i,j) =Hb_t_n(i,j) + Hb_t_t*dt
              Hb_x_np1(i,j) =Hb_x_n(i,j) + Hb_x_t*dt     
              Hb_y_np1(i,j) =Hb_y_n(i,j) + Hb_y_t*dt     
              phi1_np1(i,j) =phi1_n(i,j) + phi1_t*dt
     &                     + phi1_tt*dt**2/2
  
            end if

          end do
        end do

        call axi_reg_phi(phi1_nm1,chr,ex,L,x,y,Nx,Ny)
        call axi_reg_phi(phi1_np1,chr,ex,L,x,y,Nx,Ny)
        call axi_reg_g(gb_tt_nm1,gb_tx_nm1,gb_ty_nm1,gb_xx_nm1,
     &                 gb_xy_nm1,gb_yy_nm1,psi_nm1,omega_nm1,
     &                 chr,ex,L,x,y,Nx,Ny)
        call axi_reg_g(gb_tt_np1,gb_tx_np1,gb_ty_np1,gb_xx_np1,
     &                 gb_xy_np1,gb_yy_np1,psi_np1,omega_np1,
     &                 chr,ex,L,x,y,Nx,Ny)
        call axi_reg_Hb(Hb_t_nm1,Hb_x_nm1,Hb_y_nm1,chr,ex,L,x,y,Nx,Ny)
        call axi_reg_Hb(Hb_t_np1,Hb_x_np1,Hb_y_np1,chr,ex,L,x,y,Nx,Ny)
        call axi_reg_fb(fb_t_nm1,fb_x_nm1,fb_y_nm1,chr,ex,L,x,y,Nx,Ny)
        call axi_reg_fb(fb_t_np1,fb_x_np1,fb_y_np1,chr,ex,L,x,y,Nx,Ny)

        return
        end
