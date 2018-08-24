c----------------------------------------------------------------------
c in polar coordinates t==t, x==rho
c
c this routine calculates Hb, given gb, d(gb)dt.
c----------------------------------------------------------------------
        subroutine init_hb(gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                     gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                     gb_ty_np1,gb_ty_n,gb_ty_nm1,
     &                     gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                     gb_xy_np1,gb_xy_n,gb_xy_nm1,
     &                     gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &                     psi_np1,psi_n,psi_nm1,
     &                     omega_np1,omega_n,omega_nm1,
     &                     Hb_t_n,Hb_x_n,Hb_y_n,
     &                     L,phys_bdy,x,y,dt,chr,ex,Nx,Ny)
        implicit none
        integer Nx,Ny
        integer phys_bdy(4)
        real*8 dt,ex,L
        real*8 chr(Nx,Ny)
        real*8 Hb_t_n(Nx,Ny),Hb_x_n(Nx,Ny),Hb_y_n(Nx,Ny)
        real*8 gb_tt_np1(Nx,Ny),gb_tt_n(Nx,Ny),gb_tt_nm1(Nx,Ny)
        real*8 gb_tx_np1(Nx,Ny),gb_tx_n(Nx,Ny),gb_tx_nm1(Nx,Ny)
        real*8 gb_ty_np1(Nx,Ny),gb_ty_n(Nx,Ny),gb_ty_nm1(Nx,Ny)
        real*8 gb_xx_np1(Nx,Ny),gb_xx_n(Nx,Ny),gb_xx_nm1(Nx,Ny)
        real*8 gb_xy_np1(Nx,Ny),gb_xy_n(Nx,Ny),gb_xy_nm1(Nx,Ny)
        real*8 gb_yy_np1(Nx,Ny),gb_yy_n(Nx,Ny),gb_yy_nm1(Nx,Ny)
        real*8 psi_np1(Nx,Ny),psi_n(Nx,Ny),psi_nm1(Nx,Ny)
        real*8 omega_np1(Nx,Ny),omega_n(Nx,Ny),omega_nm1(Nx,Ny)
        real*8 x(Nx),y(Ny)

        integer i,j,is,ie,js,je,a,b,c,d

        real*8 x0,y0

        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 dx,dy

        logical ltrace,extrap_int_boundaries
        parameter (ltrace=.false.,extrap_int_boundaries=.true.)

        real*8 boxx_u(3),boxx_l(3)

        real*8 zeros(Nx,Ny)

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
        real*8 Hads_l(3),Hads_l_x(3,3),A_l(3),A_l_x(3,3)
        real*8 phi10_x(3),phi10_xx(3,3)

        real*8 f1_l(3),f1_l_x(3,3),f2_ll(3,3),f2_ll_x(3,3,3)
        real*8 f1ads_l(3),f1ads_l_x(3,3),f2ads_ll(3,3),f2ads_ll_x(3,3,3)
        real*8 h1_l(3),h1_l_x(3,3),h2_ll(3,3),h2_ll_x(3,3,3)
        real*8 sA,sB,tA,tB
        real*8 sqrtdetg,sqrtdetg_x(3)
        real*8 sqrtdetgads,sqrtdetgads_x(3)
        real*8 sqrtdeth,sqrtdeth_x(3)

        !--------------------------------------------------------------
        ! initialize fixed-size variables 
        !--------------------------------------------------------------
        data i,j,is,ie,js,je,a,b,c,d/0,0,0,0,0,0,0,0,0,0/

        data boxx_u,boxx_l/3*0.0,3*0.0/

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
        data sA,sB,tA,tB/0.0,0.0,0.0,0.0/
        data sqrtdetg,sqrtdetg_x/0.0,3*0.0/
        data sqrtdetgads,sqrtdetgads_x/0.0,3*0.0/
        data sqrtdeth,sqrtdeth_x/0.0,3*0.0/

        data A_l,Hads_l/3*0.0,3*0.0/
        data A_l_x,Hads_l_x/9*0.0,9*0.0/

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

        !---------------------------------------------------------------

        dx=x(2)-x(1)
        dy=y(2)-y(1)

        do i=1,Nx
          do j=1,Ny
            Hb_t_n(i,j)=0
            Hb_x_n(i,j)=0
            Hb_y_n(i,j)=0
            zeros(i,j)=0
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
            if (chr(i,j).ne.ex) then
              x0=x(i)
              y0=y(j)

              !-----------------------------------------------------------
              ! some other initializion
              !-----------------------------------------------------------

              ! computes tensors at point i,j 
              call tensor_init(
     &                gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                gb_ty_np1,gb_ty_n,gb_ty_nm1,
     &                gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                gb_xy_np1,gb_xy_n,gb_xy_nm1,
     &                gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &                psi_np1,psi_n,psi_nm1,
     &                omega_np1,omega_n,omega_nm1,
     &                zeros,zeros,zeros,
     &                zeros,zeros,zeros,
     &                zeros,zeros,zeros,
     &                zeros,zeros,zeros,
     &                zeros,zeros,zeros,
     &                zeros,zeros,zeros,
     &                zeros,zeros,zeros,
     &                g0_ll,g0_uu,g0_ll_x,g0_uu_x,g0_ll_xx,
     &                gads_ll,gads_uu,gads_ll_x,gads_uu_x,gads_ll_xx,
     &                gA,gB,gA_x,gB_x,gA_xx,gB_xx,
     &                gAads,gBads,gAads_x,gBads_x,gAads_xx,gBads_xx,
     &                sqrtdetg,sqrtdetg_x,
     &                sqrtdetgads,sqrtdetgads_x,
     &                h0_ll,h0_uu,h0_ll_x,h0_uu_x,h0_ll_xx,
     &                hA,hB,hAu,hBu,hA_x,hB_x,hA_xx,hB_xx,
     &                sqrtdeth,sqrtdeth_x,
     &                A_l,A_l_x,Hads_l,Hads_l_x,
     &                gamma_ull,gamma_ull_x,
     &                riemann_ulll,ricci_ll,ricci_lu,ricci,
     &                f1_l,f1_l_x,f2_ll,f2_ll_x,
     &                f1ads_l,f1ads_l_x,f2ads_ll,f2ads_ll_x,
     &                h1_l,h1_l_x,h2_ll,h2_ll_x,
     &                sA,sB,tA,tB,
     &                phi10_x,phi10_xx,
     &                x,y,dt,chr,L,ex,Nx,Ny,i,j)

              ! calculate boxx^c at point i 
              ! (boxx^c = -g^ab gamma^c_ab)
              do c=1,3
                boxx_u(c)=-( gamma_ull(c,1,1)*g0_uu(1,1)+
     &                       gamma_ull(c,2,2)*g0_uu(2,2)+
     &                       gamma_ull(c,3,3)*g0_uu(3,3)+
     &                    2*(gamma_ull(c,1,2)*g0_uu(1,2)+
     &                       gamma_ull(c,1,3)*g0_uu(1,3)+
     &                       gamma_ull(c,2,3)*g0_uu(2,3)) )
              end do

              ! calculate boxx_a at point i 
              ! (boxx_a = g_ab boxx^b)
              do a=1,3
                boxx_l(a)=boxx_u(1)*g0_ll(a,1)+
     &                    boxx_u(2)*g0_ll(a,2)+
     &                    boxx_u(3)*g0_ll(a,3)
              end do

              ! here have \box{x}_a = Hads_a + (1-x0**2)**2*H0_a
              !           \box{x}_m = Hads_m + (1-x0**2)**3*H0_m
              ! for x^a=(t,x,\Omega_3), x^m=(y,\Omega_4)
              Hb_t_n(i,j)=(boxx_l(1)-Hads_l(1))/(1-x0**2)**2
              Hb_x_n(i,j)=(boxx_l(2)-Hads_l(2))/(1-x0**2)**2
              Hb_y_n(i,j)=(boxx_l(3)-Hads_l(3))/(1-x0**2)**3

            end if
          end do
        end do

        if (extrap_int_boundaries) then
           if (phys_bdy(1).eq.0) then
              do j=1,Ny
                    Hb_t_n(1,j)=4*Hb_t_n(2,j) - 6*Hb_t_n(3,j) +
     &                            4*Hb_t_n(4,j) -   Hb_t_n(5,j)
                    Hb_x_n(1,j)=4*Hb_x_n(2,j) - 6*Hb_x_n(3,j) +
     &                            4*Hb_x_n(4,j) -   Hb_x_n(5,j)
                    Hb_y_n(1,j)=4*Hb_y_n(2,j) - 6*Hb_y_n(3,j) +
     &                            4*Hb_y_n(4,j) -   Hb_y_n(5,j)
              end do
           end if
           if (phys_bdy(2).eq.0) then
              do j=1,Ny
                    Hb_t_n(Nx,j)=4*Hb_t_n(Nx-1,j)-6*Hb_t_n(Nx-2,j)
     &                            +4*Hb_t_n(Nx-3,j)-  Hb_t_n(Nx-4,j)
                    Hb_x_n(Nx,j)=4*Hb_x_n(Nx-1,j)-6*Hb_x_n(Nx-2,j)
     &                            +4*Hb_x_n(Nx-3,j)-  Hb_x_n(Nx-4,j)
                    Hb_y_n(Nx,j)=4*Hb_y_n(Nx-1,j)-6*Hb_y_n(Nx-2,j)
     &                            +4*Hb_y_n(Nx-3,j)-  Hb_y_n(Nx-4,j)
              end do
           end if
           if (phys_bdy(3).eq.0) then
              do i=1,Nx
                    Hb_t_n(i,1)=4*Hb_t_n(i,2) - 6*Hb_t_n(i,3) +
     &                            4*Hb_t_n(i,4) -   Hb_t_n(i,5)
                    Hb_x_n(i,1)=4*Hb_x_n(i,2) - 6*Hb_x_n(i,3) +
     &                            4*Hb_x_n(i,4) -   Hb_x_n(i,5)
                    Hb_y_n(i,1)=4*Hb_y_n(i,2) - 6*Hb_y_n(i,3) +
     &                            4*Hb_y_n(i,4) -   Hb_y_n(i,5)
              end do
           end if
           if (phys_bdy(4).eq.0) then
              do i=1,Nx
                    Hb_t_n(i,Ny)=4*Hb_t_n(i,Ny-1)-6*Hb_t_n(i,Ny-2)
     &                            +4*Hb_t_n(i,Ny-3)-  Hb_t_n(i,Ny-4)
                    Hb_x_n(i,Ny)=4*Hb_x_n(i,Ny-1)-6*Hb_x_n(i,Ny-2)
     &                            +4*Hb_x_n(i,Ny-3)-  Hb_x_n(i,Ny-4)
                    Hb_y_n(i,Ny)=4*Hb_y_n(i,Ny-1)-6*Hb_y_n(i,Ny-2)
     &                            +4*Hb_y_n(i,Ny-3)-  Hb_y_n(i,Ny-4)
              end do
           end if
        end if

        call axi_reg_Hb(Hb_t_n,Hb_x_n,Hb_y_n,chr,ex,L,x,y,Nx,Ny)

        return
        end
