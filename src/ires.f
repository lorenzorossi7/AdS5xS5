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
     &                  psi_np1,psi_n,psi_nm1,
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
        real*8 psi_np1(Nx),psi_n(Nx),psi_nm1(Nx)
        real*8 phi1_np1(Nx),phi1_n(Nx),phi1_nm1(Nx)

        real*8 lambda5

        integer is,ie
        integer a,b,c,d,e,f,g,h

        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 dx

        real*8 boxx_u(5)
        real*8 sqrtdetmg0

        real*8 zeros(Nx)

        real*8 efe_ires(5,5)

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
        real*8 g0_ll(5,5),g0_uu(5,5)
        real*8 g0_ll_x(5,5,5),g0_uu_x(5,5,5),g0_ll_xx(5,5,5,5)
        real*8 gads_ll(5,5),gads_uu(5,5)
        real*8 gads_ll_x(5,5,5),gads_uu_x(5,5,5),gads_ll_xx(5,5,5,5)
        real*8 h0_ll(5,5),h0_uu(5,5)
        real*8 h0_ll_x(5,5,5),h0_uu_x(5,5,5),h0_ll_xx(5,5,5,5)
        real*8 gamma_ull(5,5,5),gamma_ull_x(5,5,5,5)
        real*8 riemann_ulll(5,5,5,5)
        real*8 ricci_ll(5,5),ricci_lu(5,5),ricci
        real*8 einstein_ll(5,5),set_ll(5,5)
        real*8 Hads_l(5),A_l(5),A_l_x(5,5)
        real*8 phi10_x(5),phi10_xx(5,5)
        real*8 ff_ll(10,10)

        !--------------------------------------------------------------
        ! initialize fixed-size variables 
        !--------------------------------------------------------------
        data is,ie/0,0/
        data a,b,c,d,e,f,g,h/0,0,0,0,0,0,0,0/

        data dx/0.0/

        data phi1_t,phi1_x,phi1_y/0.0,0.0,0.0/
        data phi1_tt,phi1_tx,phi1_ty,phi1_tz/0.0,0.0,0.0,0.0/
        data phi1_xx,phi1_xy,phi1_xz/0.0,0.0,0.0/
        data phi1_yy,phi1_yz/0.0,0.0/
        data phi1_zz/0.0/

        data boxx_u/5*0.0/

        data g0_ll,g0_uu/25*0.0,25*0.0/
        data gads_ll,gads_uu/25*0.0,25*0.0/
        data h0_ll,h0_uu/25*0.0,25*0.0/
        data gamma_ull/125*0.0/
        data gamma_ull_x/625*0.0/

        data g0_ll_x,g0_uu_x/125*0.0,125*0.0/
        data gads_ll_x,gads_uu_x/125*0.0,125*0.0/
        data h0_ll_x,h0_uu_x/125*0.0,125*0.0/

        data g0_ll_xx/625*0.0/
        data gads_ll_xx/625*0.0/
        data h0_ll_xx/625*0.0/

        data ricci/0.0/
        data ricci_ll,ricci_lu/25*0.0,25*0.0/
        data einstein_ll,set_ll/25*0.0,25*0.0/
        data riemann_ulll/625*0.0/

        data A_l,Hads_l/5*0.0,5*0.0/
        data A_l_x/25*0.0/

        data phi10_x/5*0.0/
        data phi10_xx/25*0.0/

        data ff_ll/100*0.0/

!----------------------------------------------------------------------
        
        dx=(x(2)-x(1))

        ! AdS5D cosmological constant
        !(lambda5=-(n-1)(n-2)/L^2) for n=5 dimensional AdS)
        lambda5=-6/L/L

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
     &              psi_np1,psi_n,psi_nm1,
     &              zeros,zeros,zeros,
     &              zeros,zeros,zeros,
     &              phi1_np1,phi1_n,phi1_nm1,
     &              g0_ll,g0_uu,g0_ll_x,g0_uu_x,g0_ll_xx,
     &              gads_ll,gads_uu,gads_ll_x,gads_uu_x,gads_ll_xx,
     &              h0_ll,h0_uu,h0_ll_x,h0_uu_x,h0_ll_xx,
     &              A_l,A_l_x,Hads_l,
     &              gamma_ull,gamma_ull_x,
     &              riemann_ulll,ricci_ll,ricci_lu,ricci,
     &              einstein_ll,set_ll,
     &              phi10_x,phi10_xx,
!     &              ff_ll,
     &              x,dt,chr,L,ex,Nx,i)

              ! calculates efe_ires functions at point i,j
              !(efe_ires_ab=G_ab+lambda5*g_ab-8*PI*T_ab)
              do a=1,5
                do b=a,5
                  efe_ires(a,b)=einstein_ll(a,b)+lambda5*g0_ll(a,b)
     &                                          -8*PI*set_ll(a,b)
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

