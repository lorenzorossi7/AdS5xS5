c----------------------------------------------------------------------
c in polar coordinates t==t, x==rho
c
c this routine calculates Hb, given gb, d(gb)dt.
c----------------------------------------------------------------------
        subroutine init_hb(gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                     gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                     gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                     psi_np1,psi_n,psi_nm1,
     &                     Hb_t_n,Hb_x_n,
     &                     L,phys_bdy,x,dt,chr,ex,Nx)
        implicit none
        integer Nx
        integer phys_bdy(2)
        real*8 dt,ex,L
        real*8 chr(Nx)
        real*8 Hb_t_n(Nx),Hb_x_n(Nx)
        real*8 gb_tt_np1(Nx),gb_tt_n(Nx),gb_tt_nm1(Nx)
        real*8 gb_tx_np1(Nx),gb_tx_n(Nx),gb_tx_nm1(Nx)
        real*8 gb_xx_np1(Nx),gb_xx_n(Nx),gb_xx_nm1(Nx)
        real*8 psi_np1(Nx),psi_n(Nx),psi_nm1(Nx)
        real*8 x(Nx)

        integer i,is,ie,a,b,c,d

        real*8 x0

        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 dx

        logical ltrace,extrap_int_boundaries
        parameter (ltrace=.false.,extrap_int_boundaries=.true.)

        real*8 boxx_u(5),boxx_l(5)

        real*8 zeros(Nx)

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
        data i,is,ie,a,b,c,d/0,0,0,0,0,0,0/

        data boxx_u,boxx_l/5*0.0,5*0.0/

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

        !---------------------------------------------------------------

        dx=x(2)-x(1)

        do i=1,Nx
          Hb_t_n(i)=0
          Hb_x_n(i)=0
          zeros(i)=0
        end do

        is=1
        ie=Nx
        if (phys_bdy(1).eq.1) is=2
        if (phys_bdy(2).eq.1) ie=Nx-1 

        do i=is,ie
          if (chr(i).ne.ex) then
            x0=x(i)

            !-----------------------------------------------------------
            ! some other initializion
            !-----------------------------------------------------------

            ! computes tensors at point i 
            call tensor_init(
     &              gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &              gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &              gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &              psi_np1,psi_n,psi_nm1,
     &              zeros,zeros,zeros,
     &              zeros,zeros,zeros,
     &              zeros,zeros,zeros,
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

            ! calculate boxx^c at point i 
            ! (boxx^c = -g^ab gamma^c_ab)
            do c=1,5
              boxx_u(c)=-( gamma_ull(c,1,1)*g0_uu(1,1)+
     &                     gamma_ull(c,2,2)*g0_uu(2,2)+
     &                     gamma_ull(c,3,3)*g0_uu(3,3)+
     &                     gamma_ull(c,4,4)*g0_uu(4,4)+
     &                     gamma_ull(c,5,5)*g0_uu(5,5)+
     &                  2*(gamma_ull(c,1,2)*g0_uu(1,2)+
     &                     gamma_ull(c,1,3)*g0_uu(1,3)+
     &                     gamma_ull(c,1,4)*g0_uu(1,4)+
     &                     gamma_ull(c,1,5)*g0_uu(1,5)+
     &                     gamma_ull(c,2,3)*g0_uu(2,3)+
     &                     gamma_ull(c,2,4)*g0_uu(2,4)+
     &                     gamma_ull(c,2,5)*g0_uu(2,5)+
     &                     gamma_ull(c,3,4)*g0_uu(3,4)+
     &                     gamma_ull(c,3,5)*g0_uu(3,5)+
     &                     gamma_ull(c,4,5)*g0_uu(4,5)) )
            end do

            ! calculate boxx_a at point i 
            ! (boxx_a = g_ab boxx^b)
            do a=1,5
              boxx_l(a)=boxx_u(1)*g0_ll(a,1)+
     &                  boxx_u(2)*g0_ll(a,2)+
     &                  boxx_u(3)*g0_ll(a,3)+
     &                  boxx_u(4)*g0_ll(a,4)+
     &                  boxx_u(5)*g0_ll(a,5)
            end do

            ! here have \box{x}_a = Hads_a + e_a*H0_a, where
            ! e_t=(1-x0^2)^3, e_x=(1-x0^2)^2 
            Hb_t_n(i)=(boxx_l(1)-Hads_l(1))/(1-x0**2)**3
            Hb_x_n(i)=(boxx_l(2)-Hads_l(2))/(1-x0**2)**2

          end if
        end do

        if (extrap_int_boundaries) then
          if (phys_bdy(1).eq.0) then
            Hb_t_n(1)=4*Hb_t_n(2)-6*Hb_t_n(3)
     &               +4*Hb_t_n(4)-Hb_t_n(5)
            Hb_x_n(1)=4*Hb_x_n(2)-6*Hb_x_n(3)
     &               +4*Hb_x_n(4)-Hb_x_n(5)
          end if
          if (phys_bdy(2).eq.0) then
            Hb_t_n(Nx)=4*Hb_t_n(Nx-1)-6*Hb_t_n(Nx-2)
     &                +4*Hb_t_n(Nx-3)-Hb_t_n(Nx-4)
            Hb_x_n(Nx)=4*Hb_x_n(Nx-1)-6*Hb_x_n(Nx-2) 
     &                +4*Hb_x_n(Nx-3)-Hb_x_n(Nx-4)
          end if
        end if

        call axi_reg_Hb(Hb_t_n,Hb_x_n,chr,ex,L,x,Nx)

        return
        end
