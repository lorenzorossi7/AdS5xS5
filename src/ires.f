c----------------------------------------------------------------------
c in polar coordinates t==t, x==rho
c
c routine for computing independent residuals of the AdS5D system
c
c NOTE: at the moment, iresphi1 is only for pure Poincare AdS runs
c (fixed using correct g=gads+(1-x^2)gb factorization)
c----------------------------------------------------------------------
        subroutine ires(efe_all_ires,kg_ires,
     &                  gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                  gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                  gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                  gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &                  psi_np1,psi_n,psi_nm1,
     &                  omega_np1,omega_n,omega_nm1,
     &                  phi1_np1,phi1_n,phi1_nm1,
     &                  x,dt,chr,L,ex,Nx,phys_bdy,ghost_width)
        implicit none
        integer Nx
        integer i
        integer phys_bdy(2),ghost_width(2)
        real*8 efe_all_ires(Nx),kg_ires(Nx)
        real*8 chr(Nx),ex
        real*8 x(Nx),dt,L
        real*8 gb_tt_np1(Nx),gb_tt_n(Nx),gb_tt_nm1(Nx)
        real*8 gb_tx_np1(Nx),gb_tx_n(Nx),gb_tx_nm1(Nx)
        real*8 gb_xx_np1(Nx),gb_xx_n(Nx),gb_xx_nm1(Nx)
        real*8 gb_yy_np1(Nx),gb_yy_n(Nx),gb_yy_nm1(Nx)
        real*8 psi_np1(Nx),psi_n(Nx),psi_nm1(Nx)
        real*8 omega_np1(Nx),omega_n(Nx),omega_nm1(Nx)
        real*8 phi1_np1(Nx),phi1_n(Nx),phi1_nm1(Nx)

        integer dimA,dimB
        integer is,ie
        integer a,b,c,d,e,f,g,h

        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 dx

        real*8 boxx_u(3)
        real*8 sqrtdetmg0

        real*8 zeros(Nx)

        real*8 efe_ires(3,3)

        real*8 fterm(3,3)

        !--------------------------------------------------------------
        ! the following are first and second time derivatives of *n*
        ! level variables, and as these are the only derivatives we
        ! use we drop any _n identifier
        !--------------------------------------------------------------
        real*8 phi1_t,phi1_x,phi1_y,phi1_z
        real*8 phi1_tt,phi1_tx,phi1_ty,phi1_tz
        real*8 phi1_xx,phi1_xy,phi1_xz
        real*8 phi1_yy,phi1_yz
        real*8 phi1_zz

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
        real*8 gamma_ull(3,3,3),gamma_ull_x(3,3,3,3)
        real*8 riemann_ulll(3,3,3,3)
        real*8 ricci_ll(3,3),ricci_lu(3,3),ricci
        real*8 einstein_ll(3,3),set_ll(3,3)
        real*8 f1_l(3),f2_ll(3,3)
        real*8 Hads_l(3),Hads_l_x(3,3),A_l(3),A_l_x(3,3)
        real*8 phi10_x(3),phi10_xx(3,3)

        !--------------------------------------------------------------
        ! initialize fixed-size variables 
        !--------------------------------------------------------------
        data is,ie/0,0/
        data a,b,c,d,e,f,g,h/0,0,0,0,0,0,0,0/

        data dx/0.0/

        data fterm/9*0.0/

        data phi1_t,phi1_x,phi1_y/0.0,0.0,0.0/
        data phi1_tt,phi1_tx,phi1_ty,phi1_tz/0.0,0.0,0.0,0.0/
        data phi1_xx,phi1_xy,phi1_xz/0.0,0.0,0.0/
        data phi1_yy,phi1_yz/0.0,0.0/
        data phi1_zz/0.0/

        data boxx_u/3*0.0/

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
        data einstein_ll,set_ll/9*0.0,9*0.0/
        data f1_l,f2_ll/3*0.0,9*0.0/
        data riemann_ulll/81*0.0/

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

!----------------------------------------------------------------------
        
        dx=(x(2)-x(1))

        ! set dimensions of S3 and S4 subspaces
        dimA=3
        dimB=4

        ! initialize 
        do i=1,Nx
          zeros(i)=0
        end do

        ! set index bounds for main loop
        is=2
        ie=Nx-1

        ! adjust index bounds to compensate for ghost_width
        if (ghost_width(1).gt.0) is=is+ghost_width(1)-2
        if (ghost_width(2).gt.0) ie=ie-(ghost_width(2)-2)

        ! (MAIN LOOP) loop through spacetime points x(i),y(j),z(k)
        do i=is,ie

          if (chr(i).ne.ex) then

            ! computes tensors at point i
            call tensor_init(
     &              gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &              gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &              gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &              gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &              psi_np1,psi_n,psi_nm1,
     &              omega_np1,omega_n,omega_nm1,
     &              zeros,zeros,zeros,
     &              zeros,zeros,zeros,
     &              phi1_np1,phi1_n,phi1_nm1,
     &              g0_ll,g0_uu,g0_ll_x,g0_uu_x,g0_ll_xx,
     &              gads_ll,gads_uu,gads_ll_x,gads_uu_x,gads_ll_xx,
     &              gA,gB,gA_x,gB_x,gA_xx,gB_xx,
     &              gAads,gBads,gAads_x,gBads_x,gAads_xx,gBads_xx,
     &              h0_ll,h0_uu,h0_ll_x,h0_uu_x,h0_ll_xx,
     &              A_l,A_l_x,Hads_l,Hads_l_x,
     &              gamma_ull,gamma_ull_x,
     &              riemann_ulll,ricci_ll,ricci_lu,ricci,
     &              einstein_ll,set_ll,f1_l,f2_ll,
     &              phi10_x,phi10_xx,
     &              x,dt,chr,L,ex,Nx,i)

              ! calculates efe_ires functions at point i,j
              !(efe_ires_ab=G_ab+fterm_ab)
              do a=1,3
                do b=a,3
                  fterm(a,b)=-(f1_l(a)*f1_l(b)
     &                        +g0_uu(1,1)*f2_ll(a,1)*f2_ll(b,1)
     &                        +g0_uu(1,2)*f2_ll(a,1)*f2_ll(b,2)
     &                        +g0_uu(1,3)*f2_ll(a,1)*f2_ll(b,3)
     &                        +g0_uu(2,1)*f2_ll(a,2)*f2_ll(b,1)
     &                        +g0_uu(2,2)*f2_ll(a,2)*f2_ll(b,2)
     &                        +g0_uu(2,3)*f2_ll(a,2)*f2_ll(b,3)
     &                        +g0_uu(3,1)*f2_ll(a,3)*f2_ll(b,1)
     &                        +g0_uu(3,2)*f2_ll(a,3)*f2_ll(b,2)
     &                        +g0_uu(3,3)*f2_ll(a,3)*f2_ll(b,3)
     &                        )/4

                  efe_ires(a,b)=einstein_ll(a,b)+fterm(a,b)
                end do
              end do

              ! calculate efe_all_ires function at point i,j
              efe_all_ires(i)=
     &        max(abs(efe_ires(1,1)),abs(efe_ires(1,2)),
     &            abs(efe_ires(2,2)),abs(efe_ires(3,3)))

          end if

        end do

        return
        end

