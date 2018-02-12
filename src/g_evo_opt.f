c----------------------------------------------------------------------
c in polar coordinates t,x for x in [0,1] 
c
c evolution routine for phi1 (not yet gb), computing
c the residual at the time just prior to updated
c
c L below is AdS_L (the length scale)
c----------------------------------------------------------------------
        subroutine g_evo_opt(gb_res,kg_res,
     &                       gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                       gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                       gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                       psi_np1,psi_n,psi_nm1,
     &                       Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &                       Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &                       phi1_np1,phi1_n,phi1_nm1,
     &                       L,x,dt,chr,ex,
     &                       phys_bdy,ghost_width,Nx,
     &                       background,kappa_cd,rho_cd)
        implicit none
        integer Nx
        integer phys_bdy(2),ghost_width(2)
        integer background
        real*8 kappa_cd,rho_cd
        real*8 gb_res(Nx),kg_res(Nx)
        real*8 gb_tt_np1(Nx),gb_tt_n(Nx),gb_tt_nm1(Nx)
        real*8 gb_tx_np1(Nx),gb_tx_n(Nx),gb_tx_nm1(Nx)
        real*8 gb_xx_np1(Nx),gb_xx_n(Nx),gb_xx_nm1(Nx)
        real*8 psi_np1(Nx),psi_n(Nx),psi_nm1(Nx)
        real*8 Hb_t_np1(Nx),Hb_t_n(Nx),Hb_t_nm1(Nx)
        real*8 Hb_x_np1(Nx),Hb_x_n(Nx),Hb_x_nm1(Nx)
        real*8 phi1_np1(Nx),phi1_n(Nx),phi1_nm1(Nx)
        real*8 L
        real*8 x(Nx),dt,chr(Nx),ex

        integer a,b,c,d,e

        integer rb,i
        real*8 PI
        parameter (PI=3.141592653589793d0)
        real*8 dx

        integer max_ghost_width

        real*8 boxx_u(5),boxx_l(5)

        !--------------------------------------------------------------
        ! the following are first and second time derivatives of *n*
        ! level variables, and as these are the only derivatives we
        ! use we drop any _n identifier
        !--------------------------------------------------------------
        real*8 phi1_t,phi1_x
        real*8 phi1_tt,phi1_tx
        real*8 phi1_xx

        real*8 x0

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
        ! given z_xx=(1-x0**2), z_xm=(1-x0**2)**2, z_mn=(1-x0**2)
        ! h0_ll(a,b)           = z_ab*gb_ab   
        ! h0_uu(a,b)           = inverse(gads+z_ab*gb)^ab - gads^ab
        ! h0_ll_x(a,b,c)       = (z_ab*gb_ab)_,c
        ! h0_ll_xx(a,b,c,d)    = (z_ab*gb_ab)_,cd
        ! h0_uu_x(a,b,c)       = g^ab_,c - gads^ab_,c
        !
        ! gammagg(a,b,c) = 0.5d0*gads_uu(a,d) 
        !                  *(gads_ll_x(c,d,b)-gads_ll_x(b,c,d)+gads_ll_x(d,b,c))
        ! gammahh(a,b,c) = 0.5d0*h0_uu(a,d) 
        !                  *(h0_ll_x(c,d,b)-h0_ll_x(b,c,d)+h0_ll_x(d,b,c))
        ! gammagh(a,b,c) = 0.5d0*gads_uu(a,d) 
        !                  *(h0_ll_x(c,d,b)-h0_ll_x(b,c,d)+h0_ll_x(d,b,c))
        ! gammahg(a,b,c) = 0.5d0*h0_uu(a,d) 
        !                  *(gads_ll_x(c,d,b)-gads_ll_x(b,c,d)+gads_ll_x(d,b,c))
        !
        ! cuuuu(a,b,c,d) = gads_uu(a,b)*gads_uu(c,d)
        !                   +h0_uu(a,b)*h0_uu(c,d) 
        !                   +gads_uu(a,b)*h0_uu(c,d) 
        !                   +h0_uu(a,b)*gads_uu(c,d) 
        !
        ! dlll(a,b,c)    = g0_ll_x(b,c,a)-g0_ll_x(a,b,c)+g0_ll_x(c,a,b)
        !
        ! given z_t=(1-x0^2)^3, z_x=(1-x0^2)^2, z_y=(1-x0^2)^3
        ! A_l(a)     = z_a*Hb_a   
        ! A_l_x(a,b) = (z_a*Hb_a)_,b
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
        ! tr_set= g^cd*set_cd 
        !
        ! efe(a,b) = residual ... hardcoded expressions (see below)
        !
        ! t,x,y=1,2,3
        ! 
        ! NOTE: g0_ll_xx,gads_ll_xx,h0_ll_xx,efe,efe_J
        !       do *NOT* symmetric components filled in
        !
        !--------------------------------------------------------------
        real*8 efe(5,5),efe_J(5,5)
        real*8 term1(5,5),term2(5,5),term3(5,5),term4(5,5)
        real*8 term5(5,5),term6(5,5),term7(5,5),term8(5,5)
        real*8 gammagg(5,5,5),gammahh(5,5,5)
        real*8 gammagh(5,5,5),gammahg(5,5,5) 
        real*8 cuuuu(5,5,5,5),dlll(5,5,5)
 
        real*8 alpha,ndotc,n_l(5),n_u(5),c_l(5),c_J_l(5)
        real*8 cd_ll(5,5),cd_J_ll(5,5)

        real*8 tr_set,grad_phi1_sq
        
        real*8 grad_phi4_sq

        real*8 g0u_tt_ads0,g0u_xx_ads0,g0u_xy_ads0,g0u_yy_ads0

        real*8 H0_t_ads0,H0_x_ads0,H0_y_ads0

        real*8 dgb_J,ddgb_J,ddgb_J_tx
        real*8 dphi1_J,ddphi1_J,ddphi1_J_tx

        real*8 lambda5

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

        !--------------------------------------------------------------
        ! initialize fixed-size variables 
        !--------------------------------------------------------------
        data a,b,c,d,e/0,0,0,0,0/
        data rb,i/0,0/
        data max_ghost_width/0/
        data is,ie,js,je,ks,ke,is_a_nan/0,0,0,0,0,0,0/

        data lambda5/0.0/

        data g0u_tt_ads0,g0u_xx_ads0/0.0,0.0/
        data g0u_xy_ads0,g0u_yy_ads0/0.0,0.0/

        data H0_t_ads0,H0_x_ads0,H0_y_ads0/0.0,0.0,0.0/

        data dgb_J,ddgb_J,ddgb_J_tx/0.0,0.0,0.0/
        data dphi1_J,ddphi1_J/0.0,0.0/

        data dlll/125*0.0/
        data cuuuu/625*0.0/

        data term1,term2/25*0.0,25*0.0/
        data term3,term4/25*0.0,25*0.0/
        data term5,term6/25*0.0,25*0.0/
        data term7,term8/25*0.0,25*0.0/

        data efe,efe_J/25*0.0,25*0.0/
        data cd_ll,cd_J_ll/25*0.0,25*0.0/

        data phi1_t,phi1_x/0.0,0.0/
        data phi1_tt,phi1_tx/0.0,0.0/
        data phi1_xx/0.0/

        data x0/0.0/

        data dx/0.0/

        data boxx_u,boxx_l/5*0.0,5*0.0/

        data phi1_res,phi1_J/0.0,0.0/
        data phi10_res,phi10_J/0.0,0.0/

        data n_l,n_u,c_l,c_J_l/5*0.0,5*0.0,5*0.0,5*0.0/

        data gammagg,gammahh/125*0.0,125*0.0/
        data gammagh,gammahg/125*0.0,125*0.0/

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

        !--------------------------------------------------------------
        if (ltrace) write(*,*) 'gb_psi_evo ... N=',Nx

        dx=x(2)-x(1)

        ! AdS5D cosmological constant
        !(lambda5=-(n-1)(n-2)/L^2) for n=5 dimensional AdS)
        lambda5=-6/L/L

        ! initialize output variables
        do i=1,Nx
          kg_res(i)=0
          if (chr(i).eq.ex) then
            phi1_np1(i) = 0
          end if
        end do

        ! set index bounds for main loop
        is=2
        ie=Nx-1

        ! adjust index bounds to compensate for ghost_width
        if (ghost_width(1).gt.0) is=is+ghost_width(1)-1
        if (ghost_width(2).gt.0) ie=ie-(ghost_width(2)-1)

        !(MAIN LOOP) red-black loop through spacetime pts x(i)
        do rb=0,1

          do i=is+mod(rb,2),ie,2
            x0=x(i)

            dump=.false.

            if (ltrace) write(*,*) 'i:',i

            !(REGION) interior points; evolve 
            if (chr(i).ne.ex) then

              ! computes tensors at point i
              call tensor_init(
     &                gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                psi_np1,psi_n,psi_nm1,
     &                Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &                Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &                phi1_np1,phi1_n,phi1_nm1,
     &                g0_ll,g0_uu,g0_ll_x,g0_uu_x,g0_ll_xx,
     &                gads_ll,gads_uu,gads_ll_x,gads_uu_x,gads_ll_xx,
     &                h0_ll,h0_uu,h0_ll_x,h0_uu_x,h0_ll_xx,
     &                A_l,A_l_x,Hads_l,
     &                gamma_ull,gamma_ull_x,
     &                riemann_ulll,ricci_ll,ricci_lu,ricci,
     &                einstein_ll,set_ll,
     &                phi10_x,phi10_xx,
     &                x,dt,chr,L,ex,Nx,i)

              ! computes auxiliary objects at point i
              do a=1,5
                do b=1,5
                  do c=1,5
                    dlll(a,b,c)=
     &                  g0_ll_x(b,c,a)-g0_ll_x(a,b,c)+g0_ll_x(c,a,b)
                    gammagg(a,b,c)=0
                    gammahh(a,b,c)=0
                    gammagh(a,b,c)=0
                    gammahg(a,b,c)=0
                    do d=1,5
                      gammagg(a,b,c)=gammagg(a,b,c)
     &                               +0.5d0*gads_uu(a,d)
     &                                    *(gads_ll_x(c,d,b)
     &                                     -gads_ll_x(b,c,d)
     &                                     +gads_ll_x(d,b,c))
                      gammahh(a,b,c)=gammahh(a,b,c)
     &                               +0.5d0*h0_uu(a,d)
     &                                    *(h0_ll_x(c,d,b)
     &                                     -h0_ll_x(b,c,d)
     &                                     +h0_ll_x(d,b,c))
                      gammagh(a,b,c)=gammagh(a,b,c)
     &                               +0.5d0*gads_uu(a,d)
     &                                    *(h0_ll_x(c,d,b)
     &                                     -h0_ll_x(b,c,d)
     &                                     +h0_ll_x(d,b,c))
                      gammahg(a,b,c)=gammahg(a,b,c)
     &                               +0.5d0*h0_uu(a,d)
     &                                    *(gads_ll_x(c,d,b)
     &                                     -gads_ll_x(b,c,d)
     &                                     +gads_ll_x(d,b,c))
                      cuuuu(a,b,c,d)=gads_uu(a,b)*gads_uu(c,d)+
     &                               h0_uu(a,b)*h0_uu(c,d)+
     &                               gads_uu(a,b)*h0_uu(c,d)+
     &                               h0_uu(a,b)*gads_uu(c,d)
                    end do
                  end do
                end do
              end do

              do c=1,5
                boxx_u(c)=-( gamma_ull(c,1,1)*g0_uu(1,1)+
     &                       gamma_ull(c,2,2)*g0_uu(2,2)+
     &                       gamma_ull(c,3,3)*g0_uu(3,3)+
     &                       gamma_ull(c,4,4)*g0_uu(4,4)+
     &                       gamma_ull(c,5,5)*g0_uu(5,5)+
     &                    2*(gamma_ull(c,1,2)*g0_uu(1,2)+
     &                       gamma_ull(c,1,3)*g0_uu(1,3)+
     &                       gamma_ull(c,1,4)*g0_uu(1,4)+
     &                       gamma_ull(c,1,5)*g0_uu(1,5)+
     &                       gamma_ull(c,2,3)*g0_uu(2,3)+
     &                       gamma_ull(c,2,4)*g0_uu(2,4)+
     &                       gamma_ull(c,2,5)*g0_uu(2,5)+
     &                       gamma_ull(c,3,4)*g0_uu(3,4)+
     &                       gamma_ull(c,3,5)*g0_uu(3,5)+
     &                       gamma_ull(c,4,5)*g0_uu(4,5)) )
              end do

              do a=1,5
                boxx_l(a)=boxx_u(1)*g0_ll(a,1)+
     &                    boxx_u(2)*g0_ll(a,2)+
     &                    boxx_u(3)*g0_ll(a,3)+
     &                    boxx_u(4)*g0_ll(a,4)+
     &                    boxx_u(5)*g0_ll(a,5)
              end do
    
              do a=1,5
                c_l(a)=(Hads_l(a)+A_l(a))-boxx_l(a)
              end do
 
              alpha=1/sqrt(-g0_uu(1,1))
              n_l(1)=-alpha
              do a=1,5
                n_u(a)=n_l(1)*g0_uu(a,1)+
     &                 n_l(2)*g0_uu(a,2)+
     &                 n_l(3)*g0_uu(a,3)+
     &                 n_l(4)*g0_uu(a,4)+
     &                 n_l(5)*g0_uu(a,5)
              end do

              ndotc  =n_l(1)*c_l(1)*g0_uu(1,1)+
     &                n_l(2)*c_l(2)*g0_uu(2,2)+
     &                n_l(3)*c_l(3)*g0_uu(3,3)+
     &                n_l(4)*c_l(4)*g0_uu(4,4)+
     &                n_l(5)*c_l(5)*g0_uu(5,5)+
     &             2*(n_l(1)*c_l(2)*g0_uu(1,2)+
     &                n_l(1)*c_l(3)*g0_uu(1,3)+
     &                n_l(1)*c_l(4)*g0_uu(1,4)+
     &                n_l(1)*c_l(5)*g0_uu(1,5)+
     &                n_l(2)*c_l(3)*g0_uu(2,3)+
     &                n_l(2)*c_l(4)*g0_uu(2,4)+
     &                n_l(2)*c_l(5)*g0_uu(2,5)+
     &                n_l(3)*c_l(4)*g0_uu(3,4)+
     &                n_l(3)*c_l(5)*g0_uu(3,5)+
     &                n_l(4)*c_l(5)*g0_uu(4,5))

              grad_phi1_sq=phi10_x(1)*phi10_x(1)*g0_uu(1,1)+
     &                     phi10_x(2)*phi10_x(2)*g0_uu(2,2)+
     &                     phi10_x(3)*phi10_x(3)*g0_uu(3,3)+
     &                     phi10_x(4)*phi10_x(4)*g0_uu(4,4)+
     &                     phi10_x(5)*phi10_x(5)*g0_uu(5,5)+
     &                  2*(phi10_x(1)*phi10_x(2)*g0_uu(1,2)+
     &                     phi10_x(1)*phi10_x(3)*g0_uu(1,3)+
     &                     phi10_x(1)*phi10_x(4)*g0_uu(1,4)+
     &                     phi10_x(1)*phi10_x(5)*g0_uu(1,5)+
     &                     phi10_x(2)*phi10_x(3)*g0_uu(2,3)+
     &                     phi10_x(2)*phi10_x(4)*g0_uu(2,4)+
     &                     phi10_x(2)*phi10_x(5)*g0_uu(2,5)+
     &                     phi10_x(3)*phi10_x(4)*g0_uu(3,4)+
     &                     phi10_x(3)*phi10_x(5)*g0_uu(3,5)+
     &                     phi10_x(4)*phi10_x(5)*g0_uu(4,5))
         
              do a=1,5
                do b=1,5
                  set_ll(a,b)=phi10_x(a)*phi10_x(b)
     &                       -g0_ll(a,b)*(grad_phi1_sq/2)
                end do
              end do

              tr_set =set_ll(1,1)*g0_uu(1,1)+
     &                set_ll(2,2)*g0_uu(2,2)+
     &                set_ll(3,3)*g0_uu(3,3)+
     &                set_ll(4,4)*g0_uu(4,4)+
     &                set_ll(5,5)*g0_uu(5,5)+
     &             2*(set_ll(1,2)*g0_uu(1,2)+
     &                set_ll(1,3)*g0_uu(1,3)+
     &                set_ll(1,4)*g0_uu(1,4)+
     &                set_ll(1,5)*g0_uu(1,5)+
     &                set_ll(2,3)*g0_uu(2,3)+
     &                set_ll(2,4)*g0_uu(2,4)+
     &                set_ll(2,5)*g0_uu(2,5)+
     &                set_ll(3,4)*g0_uu(3,4)+
     &                set_ll(3,5)*g0_uu(3,5)+
     &                set_ll(4,5)*g0_uu(4,5))

              !---------------------------------------------------------------- 
              ! Analytically remove the pure AdS terms in the EFEs,
              ! denoted ->, so 
              ! efe_ab =   term1_ab + term2_ab + term3_ab + term4_ab 
              !          + term5_ab + term6_ab + term7_ab + term8_ab 
              !          - 8*PI*(se_ab-1/3*tr(se)*g_ab)
              ! for
              ! 
              ! term1_ab = -1/2 g^cd g_ab,cd  
              !          ->-1/2 (h^cd h_ab,cd + gads^cd h_ab,cd + h^cd
              !          gads_ab,cd) 
              ! term2_ab = -1/2 g^cd,a g_bc,d
              !          ->-1/2 (h^cd,a h_bc,d + gads^cd,a h_bc,d +
              !          h^cd,a gads_bc,d)
              ! term3_ab = -1/2 g^cd,b g_ac,d
              !          ->-1/2 (h^cd,b h_ac,d + gads^cd,b h_ac,d +
              !          h^cd,b gads_ac,d)
              ! term4_ab = -1/2 H_a,b
              !          ->-1/2 Hb_a,b
              ! term5_ab = -1/2 H_b,a
              !          ->-1/2 Hb_b,a
              ! term6_ab = H_c G^c_ab
              !          ->Hb_c Ghh^c_ab + Hb_c Ggh^c_ab + Hb_c
              !          Ghg^c_ab  
              ! term7_ab = -G^c_db G^d_ca
              !          ->-(  Ghh^c_db Ghh^d_ca + Ggg^c_db Ghh^d_ca
              !              + Ggg^c_db Ggh^d_ca + Ggg^c_db Ghg^d_ca
              !              + Ghh^c_db Ggg^d_ca + Ghh^c_db Ggh^d_ca
              !              + Ghh^c_db Ghg^d_ca + Ggh^c_db Ggg^d_ca
              !              + Ggh^c_db Ghh^d_ca + Ggh^c_db Ggh^d_ca
              !              + Ggh^c_db Ghg^d_ca + Ghg^c_db Ggg^d_ca
              !              + Ghg^c_db Ghh^d_ca + Ghg^c_db Ggh^d_ca
              !              + Ghg^c_db Ghg^d_ca  )
              ! term8_ab = -2/3 lambda5 g_ab
              !          ->-2/3 lambda5 h_ab
              !
              ! where G   = guu(g_ll_x-g_ll_x+g_ll_x)
              ! where Ggg = gadsuu(gads_ll_x-gads_ll_x+gads_ll_x)
              ! where Ghh = huu(h_ll_x-h_ll_x+h_ll_x)
              ! where Ggh = gadsuu(h_ll_x-h_ll_x+h_ll_x)
              ! where Ghg = huu(gads_ll_x-gads_ll_x+gads_ll_x)
              !
              !----------------------------------------------------------------
              do a=1,3
                do b=a,3
                  term1(a,b)=-0.5d0*(                             
     &                          h0_uu(1,1)*h0_ll_xx(a,b,1,1)+
     &                          h0_uu(2,2)*h0_ll_xx(a,b,2,2)+
     &                          h0_uu(3,3)*h0_ll_xx(a,b,3,3)+
     &                          h0_uu(4,4)*h0_ll_xx(a,b,4,4)+
     &                          h0_uu(5,5)*h0_ll_xx(a,b,5,5)+
     &                       2*(h0_uu(1,2)*h0_ll_xx(a,b,1,2)+
     &                          h0_uu(1,3)*h0_ll_xx(a,b,1,3)+
     &                          h0_uu(1,4)*h0_ll_xx(a,b,1,4)+
     &                          h0_uu(1,5)*h0_ll_xx(a,b,1,5)+
     &                          h0_uu(2,3)*h0_ll_xx(a,b,2,3)+
     &                          h0_uu(2,4)*h0_ll_xx(a,b,2,4)+
     &                          h0_uu(2,5)*h0_ll_xx(a,b,2,5)+
     &                          h0_uu(3,4)*h0_ll_xx(a,b,3,4)+
     &                          h0_uu(3,5)*h0_ll_xx(a,b,3,5)+
     &                          h0_uu(4,5)*h0_ll_xx(a,b,4,5))
     &                       +
     &                          gads_uu(1,1)*h0_ll_xx(a,b,1,1)+
     &                          gads_uu(2,2)*h0_ll_xx(a,b,2,2)+
     &                          gads_uu(3,3)*h0_ll_xx(a,b,3,3)+
     &                          gads_uu(4,4)*h0_ll_xx(a,b,4,4)+
     &                          gads_uu(5,5)*h0_ll_xx(a,b,5,5)+
     &                       2*(gads_uu(1,2)*h0_ll_xx(a,b,1,2)+
     &                          gads_uu(1,3)*h0_ll_xx(a,b,1,3)+
     &                          gads_uu(1,4)*h0_ll_xx(a,b,1,4)+
     &                          gads_uu(1,5)*h0_ll_xx(a,b,1,5)+
     &                          gads_uu(2,3)*h0_ll_xx(a,b,2,3)+
     &                          gads_uu(2,4)*h0_ll_xx(a,b,2,4)+
     &                          gads_uu(2,5)*h0_ll_xx(a,b,2,5)+
     &                          gads_uu(3,4)*h0_ll_xx(a,b,3,4)+
     &                          gads_uu(3,5)*h0_ll_xx(a,b,3,5)+
     &                          gads_uu(4,5)*h0_ll_xx(a,b,4,5))
     &                       +
     &                          h0_uu(1,1)*gads_ll_xx(a,b,1,1)+
     &                          h0_uu(2,2)*gads_ll_xx(a,b,2,2)+
     &                          h0_uu(3,3)*gads_ll_xx(a,b,3,3)+
     &                          h0_uu(4,4)*gads_ll_xx(a,b,4,4)+
     &                          h0_uu(5,5)*gads_ll_xx(a,b,5,5)+
     &                       2*(h0_uu(1,2)*gads_ll_xx(a,b,1,2)+
     &                          h0_uu(1,3)*gads_ll_xx(a,b,1,3)+
     &                          h0_uu(1,4)*gads_ll_xx(a,b,1,4)+
     &                          h0_uu(1,5)*gads_ll_xx(a,b,1,5)+
     &                          h0_uu(2,3)*gads_ll_xx(a,b,2,3)+
     &                          h0_uu(2,4)*gads_ll_xx(a,b,2,4)+
     &                          h0_uu(2,5)*gads_ll_xx(a,b,2,5)+
     &                          h0_uu(3,4)*gads_ll_xx(a,b,3,4)+
     &                          h0_uu(3,5)*gads_ll_xx(a,b,3,5)+
     &                          h0_uu(4,5)*gads_ll_xx(a,b,4,5))
     &                              )
     &
                  term2(a,b)=-0.5d0*(                             
     &                          h0_uu_x(1,1,a)* h0_ll_x(b,1,1) +
     &                          h0_uu_x(1,2,a)*(h0_ll_x(b,1,2) +
     &                                          h0_ll_x(b,2,1))+
     &                          h0_uu_x(1,3,a)*(h0_ll_x(b,1,3) +
     &                                          h0_ll_x(b,3,1))+
     &                          h0_uu_x(1,4,a)*(h0_ll_x(b,1,4) +
     &                                          h0_ll_x(b,4,1))+
     &                          h0_uu_x(1,5,a)*(h0_ll_x(b,1,5) +
     &                                          h0_ll_x(b,5,1))+
     &                          h0_uu_x(2,2,a)* h0_ll_x(b,2,2) +
     &                          h0_uu_x(2,3,a)*(h0_ll_x(b,2,3) +
     &                                          h0_ll_x(b,3,2))+
     &                          h0_uu_x(2,4,a)*(h0_ll_x(b,2,4) +
     &                                          h0_ll_x(b,4,2))+
     &                          h0_uu_x(2,5,a)*(h0_ll_x(b,2,5) +
     &                                          h0_ll_x(b,5,2))+
     &                          h0_uu_x(3,3,a)* h0_ll_x(b,3,3) +
     &                          h0_uu_x(3,4,a)*(h0_ll_x(b,3,4) +
     &                                          h0_ll_x(b,4,3))+
     &                          h0_uu_x(3,5,a)*(h0_ll_x(b,3,5) +
     &                                          h0_ll_x(b,5,3))+
     &                          h0_uu_x(4,4,a)* h0_ll_x(b,4,4) +
     &                          h0_uu_x(4,5,a)*(h0_ll_x(b,4,5) +
     &                                          h0_ll_x(b,5,4))+
     &                          h0_uu_x(5,5,a)* h0_ll_x(b,5,5)
     &                       +
     &                          gads_uu_x(1,1,a)* h0_ll_x(b,1,1) +
     &                          gads_uu_x(1,2,a)*(h0_ll_x(b,1,2) +
     &                                            h0_ll_x(b,2,1))+
     &                          gads_uu_x(1,3,a)*(h0_ll_x(b,1,3) +
     &                                            h0_ll_x(b,3,1))+
     &                          gads_uu_x(1,4,a)*(h0_ll_x(b,1,4) +
     &                                            h0_ll_x(b,4,1))+
     &                          gads_uu_x(1,5,a)*(h0_ll_x(b,1,5) +
     &                                            h0_ll_x(b,5,1))+
     &                          gads_uu_x(2,2,a)* h0_ll_x(b,2,2) +
     &                          gads_uu_x(2,3,a)*(h0_ll_x(b,2,3) +
     &                                            h0_ll_x(b,3,2))+
     &                          gads_uu_x(2,4,a)*(h0_ll_x(b,2,4) +
     &                                            h0_ll_x(b,4,2))+
     &                          gads_uu_x(2,5,a)*(h0_ll_x(b,2,5) +
     &                                            h0_ll_x(b,5,2))+
     &                          gads_uu_x(3,3,a)* h0_ll_x(b,3,3) +
     &                          gads_uu_x(3,4,a)*(h0_ll_x(b,3,4) +
     &                                            h0_ll_x(b,4,3))+
     &                          gads_uu_x(3,5,a)*(h0_ll_x(b,3,5) +
     &                                            h0_ll_x(b,5,3))+
     &                          gads_uu_x(4,4,a)* h0_ll_x(b,4,4) +
     &                          gads_uu_x(4,5,a)*(h0_ll_x(b,4,5) +
     &                                            h0_ll_x(b,5,4))+
     &                          gads_uu_x(5,5,a)* h0_ll_x(b,5,5)
     &                       +
     &                          h0_uu_x(1,1,a)* gads_ll_x(b,1,1) +
     &                          h0_uu_x(1,2,a)*(gads_ll_x(b,1,2) + 
     &                                          gads_ll_x(b,2,1))+ 
     &                          h0_uu_x(1,3,a)*(gads_ll_x(b,1,3) + 
     &                                          gads_ll_x(b,3,1))+ 
     &                          h0_uu_x(1,4,a)*(gads_ll_x(b,1,4) +
     &                                          gads_ll_x(b,4,1))+
     &                          h0_uu_x(1,5,a)*(gads_ll_x(b,1,5) +
     &                                          gads_ll_x(b,5,1))+
     &                          h0_uu_x(2,2,a)* gads_ll_x(b,2,2) +
     &                          h0_uu_x(2,3,a)*(gads_ll_x(b,2,3) +
     &                                          gads_ll_x(b,3,2))+
     &                          h0_uu_x(2,4,a)*(gads_ll_x(b,2,4) +
     &                                          gads_ll_x(b,4,2))+
     &                          h0_uu_x(2,5,a)*(gads_ll_x(b,2,5) +
     &                                          gads_ll_x(b,5,2))+
     &                          h0_uu_x(3,3,a)* gads_ll_x(b,3,3) +
     &                          h0_uu_x(3,4,a)*(gads_ll_x(b,3,4) +
     &                                          gads_ll_x(b,4,3))+
     &                          h0_uu_x(3,5,a)*(gads_ll_x(b,3,5) +
     &                                          gads_ll_x(b,5,3))+
     &                          h0_uu_x(4,4,a)* gads_ll_x(b,4,4) +
     &                          h0_uu_x(4,5,a)*(gads_ll_x(b,4,5) +
     &                                          gads_ll_x(b,5,4))+
     &                          h0_uu_x(5,5,a)* gads_ll_x(b,5,5)
     &                            )
     &
                  term3(a,b)=-0.5d0*(                            
     &                          h0_uu_x(1,1,b)* h0_ll_x(a,1,1) +
     &                          h0_uu_x(1,2,b)*(h0_ll_x(a,1,2) +
     &                                          h0_ll_x(a,2,1))+
     &                          h0_uu_x(1,3,b)*(h0_ll_x(a,1,3) +
     &                                          h0_ll_x(a,3,1))+
     &                          h0_uu_x(1,4,b)*(h0_ll_x(a,1,4) +
     &                                          h0_ll_x(a,4,1))+
     &                          h0_uu_x(1,5,b)*(h0_ll_x(a,1,5) +
     &                                          h0_ll_x(a,5,1))+
     &                          h0_uu_x(2,2,b)* h0_ll_x(a,2,2) +
     &                          h0_uu_x(2,3,b)*(h0_ll_x(a,2,3) +
     &                                          h0_ll_x(a,3,2))+
     &                          h0_uu_x(2,4,b)*(h0_ll_x(a,2,4) +
     &                                          h0_ll_x(a,4,2))+
     &                          h0_uu_x(2,5,b)*(h0_ll_x(a,2,5) +
     &                                          h0_ll_x(a,5,2))+
     &                          h0_uu_x(3,3,b)* h0_ll_x(a,3,3) +
     &                          h0_uu_x(3,4,b)*(h0_ll_x(a,3,4) +
     &                                          h0_ll_x(a,4,3))+
     &                          h0_uu_x(3,5,b)*(h0_ll_x(a,3,5) +
     &                                          h0_ll_x(a,5,3))+
     &                          h0_uu_x(4,4,b)* h0_ll_x(a,4,4) +
     &                          h0_uu_x(4,5,b)*(h0_ll_x(a,4,5) +
     &                                          h0_ll_x(a,5,4))+
     &                          h0_uu_x(5,5,b)* h0_ll_x(a,5,5)
     &                       +
     &                          gads_uu_x(1,1,b)* h0_ll_x(a,1,1) +
     &                          gads_uu_x(1,2,b)*(h0_ll_x(a,1,2) +
     &                                            h0_ll_x(a,2,1))+
     &                          gads_uu_x(1,3,b)*(h0_ll_x(a,1,3) +
     &                                            h0_ll_x(a,3,1))+
     &                          gads_uu_x(1,4,b)*(h0_ll_x(a,1,4) +
     &                                            h0_ll_x(a,4,1))+
     &                          gads_uu_x(1,5,b)*(h0_ll_x(a,1,5) +
     &                                            h0_ll_x(a,5,1))+
     &                          gads_uu_x(2,2,b)* h0_ll_x(a,2,2) +
     &                          gads_uu_x(2,3,b)*(h0_ll_x(a,2,3) +
     &                                            h0_ll_x(a,3,2))+
     &                          gads_uu_x(2,4,b)*(h0_ll_x(a,2,4) +
     &                                            h0_ll_x(a,4,2))+
     &                          gads_uu_x(2,5,b)*(h0_ll_x(a,2,5) +
     &                                            h0_ll_x(a,5,2))+
     &                          gads_uu_x(3,3,b)* h0_ll_x(a,3,3) +
     &                          gads_uu_x(3,4,b)*(h0_ll_x(a,3,4) +
     &                                            h0_ll_x(a,4,3))+
     &                          gads_uu_x(3,5,b)*(h0_ll_x(a,3,5) +
     &                                            h0_ll_x(a,5,3))+
     &                          gads_uu_x(4,4,b)* h0_ll_x(a,4,4) +
     &                          gads_uu_x(4,5,b)*(h0_ll_x(a,4,5) +
     &                                            h0_ll_x(a,5,4))+
     &                          gads_uu_x(5,5,b)* h0_ll_x(a,5,5)
     &                       +
     &                          h0_uu_x(1,1,b)* gads_ll_x(a,1,1) +
     &                          h0_uu_x(1,2,b)*(gads_ll_x(a,1,2) +  
     &                                          gads_ll_x(a,2,1))+  
     &                          h0_uu_x(1,3,b)*(gads_ll_x(a,1,3) +  
     &                                          gads_ll_x(a,3,1))+   
     &                          h0_uu_x(1,4,b)*(gads_ll_x(a,1,4) +
     &                                          gads_ll_x(a,4,1))+
     &                          h0_uu_x(1,5,b)*(gads_ll_x(a,1,5) +
     &                                          gads_ll_x(a,5,1))+
     &                          h0_uu_x(2,2,b)* gads_ll_x(a,2,2) +
     &                          h0_uu_x(2,3,b)*(gads_ll_x(a,2,3) +
     &                                          gads_ll_x(a,3,2))+
     &                          h0_uu_x(2,4,b)*(gads_ll_x(a,2,4) +
     &                                          gads_ll_x(a,4,2))+
     &                          h0_uu_x(2,5,b)*(gads_ll_x(a,2,5) +
     &                                          gads_ll_x(a,5,2))+
     &                          h0_uu_x(3,3,b)* gads_ll_x(a,3,3) +
     &                          h0_uu_x(3,4,b)*(gads_ll_x(a,3,4) +
     &                                          gads_ll_x(a,4,3))+
     &                          h0_uu_x(3,5,b)*(gads_ll_x(a,3,5) +
     &                                          gads_ll_x(a,5,3))+
     &                          h0_uu_x(4,4,b)* gads_ll_x(a,4,4) +
     &                          h0_uu_x(4,5,b)*(gads_ll_x(a,4,5) +
     &                                          gads_ll_x(a,5,4))+
     &                          h0_uu_x(5,5,b)* gads_ll_x(a,5,5)
     &                            )
     &
                  term4(a,b)=-0.5d0*A_l_x(a,b)                  
     &
                  term5(a,b)=-0.5d0*A_l_x(b,a)           
     &
                  term6(a,b)=     (           
     &                          Hads_l(1)*gammahh(1,a,b)+      
     &                          Hads_l(2)*gammahh(2,a,b)+
     &                          Hads_l(3)*gammahh(3,a,b)+
     &                          Hads_l(4)*gammahh(4,a,b)+
     &                          Hads_l(5)*gammahh(5,a,b)
     &                       +
     &                          Hads_l(1)*gammagh(1,a,b)+
     &                          Hads_l(2)*gammagh(2,a,b)+
     &                          Hads_l(3)*gammagh(3,a,b)+
     &                          Hads_l(4)*gammagh(4,a,b)+
     &                          Hads_l(5)*gammagh(5,a,b)
     &                       +
     &                          Hads_l(1)*gammahg(1,a,b)+
     &                          Hads_l(2)*gammahg(2,a,b)+
     &                          Hads_l(3)*gammahg(3,a,b)+
     &                          Hads_l(4)*gammahg(4,a,b)+
     &                          Hads_l(5)*gammahg(5,a,b)
     &                       +                  
     &                          A_l(1)*gammagg(1,a,b)  +
     &                          A_l(2)*gammagg(2,a,b)  +
     &                          A_l(3)*gammagg(3,a,b)  +
     &                          A_l(4)*gammagg(4,a,b)  +
     &                          A_l(5)*gammagg(5,a,b)
     &                       +
     &                          A_l(1)*gammahh(1,a,b)  +
     &                          A_l(2)*gammahh(2,a,b)  +
     &                          A_l(3)*gammahh(3,a,b)  +
     &                          A_l(4)*gammahh(4,a,b)  +
     &                          A_l(5)*gammahh(5,a,b)   
     &                       +
     &                          A_l(1)*gammagh(1,a,b)  +
     &                          A_l(2)*gammagh(2,a,b)  +
     &                          A_l(3)*gammagh(3,a,b)  +
     &                          A_l(4)*gammagh(4,a,b)  +
     &                          A_l(5)*gammagh(5,a,b)
     &                       +
     &                          A_l(1)*gammahg(1,a,b)  +
     &                          A_l(2)*gammahg(2,a,b)  +
     &                          A_l(3)*gammahg(3,a,b)  +
     &                          A_l(4)*gammahg(4,a,b)  +
     &                          A_l(5)*gammahg(5,a,b)
     &                            ) 
     &                           
                  term7(a,b)=    -(
     &                          gammahh(1,1,b)*gammahh(1,1,a)+
     &                          gammahh(1,2,b)*gammahh(2,1,a)+
     &                          gammahh(1,3,b)*gammahh(3,1,a)+
     &                          gammahh(1,4,b)*gammahh(4,1,a)+
     &                          gammahh(1,5,b)*gammahh(5,1,a)+
     &                          gammahh(2,1,b)*gammahh(1,2,a)+
     &                          gammahh(2,2,b)*gammahh(2,2,a)+
     &                          gammahh(2,3,b)*gammahh(3,2,a)+
     &                          gammahh(2,4,b)*gammahh(4,2,a)+
     &                          gammahh(2,5,b)*gammahh(5,2,a)+
     &                          gammahh(3,1,b)*gammahh(1,3,a)+
     &                          gammahh(3,2,b)*gammahh(2,3,a)+
     &                          gammahh(3,3,b)*gammahh(3,3,a)+
     &                          gammahh(3,4,b)*gammahh(4,3,a)+
     &                          gammahh(3,5,b)*gammahh(5,3,a)+
     &                          gammahh(4,1,b)*gammahh(1,4,a)+
     &                          gammahh(4,2,b)*gammahh(2,4,a)+
     &                          gammahh(4,3,b)*gammahh(3,4,a)+
     &                          gammahh(4,4,b)*gammahh(4,4,a)+
     &                          gammahh(4,5,b)*gammahh(5,4,a)+
     &                          gammahh(5,1,b)*gammahh(1,5,a)+
     &                          gammahh(5,2,b)*gammahh(2,5,a)+
     &                          gammahh(5,3,b)*gammahh(3,5,a)+
     &                          gammahh(5,4,b)*gammahh(4,5,a)+
     &                          gammahh(5,5,b)*gammahh(5,5,a)
     &                       +
     &                          gammagg(1,1,b)*gammahh(1,1,a)+
     &                          gammagg(1,2,b)*gammahh(2,1,a)+
     &                          gammagg(1,3,b)*gammahh(3,1,a)+
     &                          gammagg(1,4,b)*gammahh(4,1,a)+
     &                          gammagg(1,5,b)*gammahh(5,1,a)+
     &                          gammagg(2,1,b)*gammahh(1,2,a)+
     &                          gammagg(2,2,b)*gammahh(2,2,a)+
     &                          gammagg(2,3,b)*gammahh(3,2,a)+
     &                          gammagg(2,4,b)*gammahh(4,2,a)+
     &                          gammagg(2,5,b)*gammahh(5,2,a)+
     &                          gammagg(3,1,b)*gammahh(1,3,a)+
     &                          gammagg(3,2,b)*gammahh(2,3,a)+
     &                          gammagg(3,3,b)*gammahh(3,3,a)+
     &                          gammagg(3,4,b)*gammahh(4,3,a)+
     &                          gammagg(3,5,b)*gammahh(5,3,a)+
     &                          gammagg(4,1,b)*gammahh(1,4,a)+
     &                          gammagg(4,2,b)*gammahh(2,4,a)+
     &                          gammagg(4,3,b)*gammahh(3,4,a)+
     &                          gammagg(4,4,b)*gammahh(4,4,a)+
     &                          gammagg(4,5,b)*gammahh(5,4,a)+
     &                          gammagg(5,1,b)*gammahh(1,5,a)+
     &                          gammagg(5,2,b)*gammahh(2,5,a)+
     &                          gammagg(5,3,b)*gammahh(3,5,a)+
     &                          gammagg(5,4,b)*gammahh(4,5,a)+
     &                          gammagg(5,5,b)*gammahh(5,5,a)
     &                       +
     &                          gammagg(1,1,b)*gammagh(1,1,a)+
     &                          gammagg(1,2,b)*gammagh(2,1,a)+
     &                          gammagg(1,3,b)*gammagh(3,1,a)+
     &                          gammagg(1,4,b)*gammagh(4,1,a)+
     &                          gammagg(1,5,b)*gammagh(5,1,a)+
     &                          gammagg(2,1,b)*gammagh(1,2,a)+
     &                          gammagg(2,2,b)*gammagh(2,2,a)+
     &                          gammagg(2,3,b)*gammagh(3,2,a)+
     &                          gammagg(2,4,b)*gammagh(4,2,a)+
     &                          gammagg(2,5,b)*gammagh(5,2,a)+
     &                          gammagg(3,1,b)*gammagh(1,3,a)+
     &                          gammagg(3,2,b)*gammagh(2,3,a)+
     &                          gammagg(3,3,b)*gammagh(3,3,a)+
     &                          gammagg(3,4,b)*gammagh(4,3,a)+
     &                          gammagg(3,5,b)*gammagh(5,3,a)+
     &                          gammagg(4,1,b)*gammagh(1,4,a)+
     &                          gammagg(4,2,b)*gammagh(2,4,a)+
     &                          gammagg(4,3,b)*gammagh(3,4,a)+
     &                          gammagg(4,4,b)*gammagh(4,4,a)+
     &                          gammagg(4,5,b)*gammagh(5,4,a)+
     &                          gammagg(5,1,b)*gammagh(1,5,a)+
     &                          gammagg(5,2,b)*gammagh(2,5,a)+
     &                          gammagg(5,3,b)*gammagh(3,5,a)+
     &                          gammagg(5,4,b)*gammagh(4,5,a)+
     &                          gammagg(5,5,b)*gammagh(5,5,a)
     &                       +
     &                          gammagg(1,1,b)*gammahg(1,1,a)+
     &                          gammagg(1,2,b)*gammahg(2,1,a)+
     &                          gammagg(1,3,b)*gammahg(3,1,a)+
     &                          gammagg(1,4,b)*gammahg(4,1,a)+
     &                          gammagg(1,5,b)*gammahg(5,1,a)+
     &                          gammagg(2,1,b)*gammahg(1,2,a)+
     &                          gammagg(2,2,b)*gammahg(2,2,a)+
     &                          gammagg(2,3,b)*gammahg(3,2,a)+
     &                          gammagg(2,4,b)*gammahg(4,2,a)+
     &                          gammagg(2,5,b)*gammahg(5,2,a)+
     &                          gammagg(3,1,b)*gammahg(1,3,a)+
     &                          gammagg(3,2,b)*gammahg(2,3,a)+
     &                          gammagg(3,3,b)*gammahg(3,3,a)+
     &                          gammagg(3,4,b)*gammahg(4,3,a)+
     &                          gammagg(3,5,b)*gammahg(5,3,a)+
     &                          gammagg(4,1,b)*gammahg(1,4,a)+
     &                          gammagg(4,2,b)*gammahg(2,4,a)+
     &                          gammagg(4,3,b)*gammahg(3,4,a)+
     &                          gammagg(4,4,b)*gammahg(4,4,a)+
     &                          gammagg(4,5,b)*gammahg(5,4,a)+
     &                          gammagg(5,1,b)*gammahg(1,5,a)+
     &                          gammagg(5,2,b)*gammahg(2,5,a)+
     &                          gammagg(5,3,b)*gammahg(3,5,a)+
     &                          gammagg(5,4,b)*gammahg(4,5,a)+
     &                          gammagg(5,5,b)*gammahg(5,5,a)
     &                       +
     &                          gammahh(1,1,b)*gammagg(1,1,a)+
     &                          gammahh(1,2,b)*gammagg(2,1,a)+
     &                          gammahh(1,3,b)*gammagg(3,1,a)+
     &                          gammahh(1,4,b)*gammagg(4,1,a)+
     &                          gammahh(1,5,b)*gammagg(5,1,a)+
     &                          gammahh(2,1,b)*gammagg(1,2,a)+
     &                          gammahh(2,2,b)*gammagg(2,2,a)+
     &                          gammahh(2,3,b)*gammagg(3,2,a)+
     &                          gammahh(2,4,b)*gammagg(4,2,a)+
     &                          gammahh(2,5,b)*gammagg(5,2,a)+
     &                          gammahh(3,1,b)*gammagg(1,3,a)+
     &                          gammahh(3,2,b)*gammagg(2,3,a)+
     &                          gammahh(3,3,b)*gammagg(3,3,a)+
     &                          gammahh(3,4,b)*gammagg(4,3,a)+
     &                          gammahh(3,5,b)*gammagg(5,3,a)+
     &                          gammahh(4,1,b)*gammagg(1,4,a)+
     &                          gammahh(4,2,b)*gammagg(2,4,a)+
     &                          gammahh(4,3,b)*gammagg(3,4,a)+
     &                          gammahh(4,4,b)*gammagg(4,4,a)+
     &                          gammahh(4,5,b)*gammagg(5,4,a)+
     &                          gammahh(5,1,b)*gammagg(1,5,a)+
     &                          gammahh(5,2,b)*gammagg(2,5,a)+
     &                          gammahh(5,3,b)*gammagg(3,5,a)+
     &                          gammahh(5,4,b)*gammagg(4,5,a)+
     &                          gammahh(5,5,b)*gammagg(5,5,a)
     &                       +
     &                          gammahh(1,1,b)*gammagh(1,1,a)+
     &                          gammahh(1,2,b)*gammagh(2,1,a)+
     &                          gammahh(1,3,b)*gammagh(3,1,a)+
     &                          gammahh(1,4,b)*gammagh(4,1,a)+
     &                          gammahh(1,5,b)*gammagh(5,1,a)+
     &                          gammahh(2,1,b)*gammagh(1,2,a)+
     &                          gammahh(2,2,b)*gammagh(2,2,a)+
     &                          gammahh(2,3,b)*gammagh(3,2,a)+
     &                          gammahh(2,4,b)*gammagh(4,2,a)+
     &                          gammahh(2,5,b)*gammagh(5,2,a)+
     &                          gammahh(3,1,b)*gammagh(1,3,a)+
     &                          gammahh(3,2,b)*gammagh(2,3,a)+
     &                          gammahh(3,3,b)*gammagh(3,3,a)+
     &                          gammahh(3,4,b)*gammagh(4,3,a)+
     &                          gammahh(3,5,b)*gammagh(5,3,a)+
     &                          gammahh(4,1,b)*gammagh(1,4,a)+
     &                          gammahh(4,2,b)*gammagh(2,4,a)+
     &                          gammahh(4,3,b)*gammagh(3,4,a)+
     &                          gammahh(4,4,b)*gammagh(4,4,a)+
     &                          gammahh(4,5,b)*gammagh(5,4,a)+
     &                          gammahh(5,1,b)*gammagh(1,5,a)+
     &                          gammahh(5,2,b)*gammagh(2,5,a)+
     &                          gammahh(5,3,b)*gammagh(3,5,a)+
     &                          gammahh(5,4,b)*gammagh(4,5,a)+
     &                          gammahh(5,5,b)*gammagh(5,5,a)
     &                       +
     &                          gammahh(1,1,b)*gammahg(1,1,a)+
     &                          gammahh(1,2,b)*gammahg(2,1,a)+
     &                          gammahh(1,3,b)*gammahg(3,1,a)+
     &                          gammahh(1,4,b)*gammahg(4,1,a)+
     &                          gammahh(1,5,b)*gammahg(5,1,a)+
     &                          gammahh(2,1,b)*gammahg(1,2,a)+
     &                          gammahh(2,2,b)*gammahg(2,2,a)+
     &                          gammahh(2,3,b)*gammahg(3,2,a)+
     &                          gammahh(2,4,b)*gammahg(4,2,a)+
     &                          gammahh(2,5,b)*gammahg(5,2,a)+
     &                          gammahh(3,1,b)*gammahg(1,3,a)+
     &                          gammahh(3,2,b)*gammahg(2,3,a)+
     &                          gammahh(3,3,b)*gammahg(3,3,a)+
     &                          gammahh(3,4,b)*gammahg(4,3,a)+
     &                          gammahh(3,5,b)*gammahg(5,3,a)+
     &                          gammahh(4,1,b)*gammahg(1,4,a)+
     &                          gammahh(4,2,b)*gammahg(2,4,a)+
     &                          gammahh(4,3,b)*gammahg(3,4,a)+
     &                          gammahh(4,4,b)*gammahg(4,4,a)+
     &                          gammahh(4,5,b)*gammahg(5,4,a)+
     &                          gammahh(5,1,b)*gammahg(1,5,a)+
     &                          gammahh(5,2,b)*gammahg(2,5,a)+
     &                          gammahh(5,3,b)*gammahg(3,5,a)+
     &                          gammahh(5,4,b)*gammahg(4,5,a)+
     &                          gammahh(5,5,b)*gammahg(5,5,a)
     &                       +
     &                          gammagh(1,1,b)*gammagg(1,1,a)+
     &                          gammagh(1,2,b)*gammagg(2,1,a)+
     &                          gammagh(1,3,b)*gammagg(3,1,a)+
     &                          gammagh(1,4,b)*gammagg(4,1,a)+
     &                          gammagh(1,5,b)*gammagg(5,1,a)+
     &                          gammagh(2,1,b)*gammagg(1,2,a)+
     &                          gammagh(2,2,b)*gammagg(2,2,a)+
     &                          gammagh(2,3,b)*gammagg(3,2,a)+
     &                          gammagh(2,4,b)*gammagg(4,2,a)+
     &                          gammagh(2,5,b)*gammagg(5,2,a)+
     &                          gammagh(3,1,b)*gammagg(1,3,a)+
     &                          gammagh(3,2,b)*gammagg(2,3,a)+
     &                          gammagh(3,3,b)*gammagg(3,3,a)+
     &                          gammagh(3,4,b)*gammagg(4,3,a)+
     &                          gammagh(3,5,b)*gammagg(5,3,a)+
     &                          gammagh(4,1,b)*gammagg(1,4,a)+
     &                          gammagh(4,2,b)*gammagg(2,4,a)+
     &                          gammagh(4,3,b)*gammagg(3,4,a)+
     &                          gammagh(4,4,b)*gammagg(4,4,a)+
     &                          gammagh(4,5,b)*gammagg(5,4,a)+
     &                          gammagh(5,1,b)*gammagg(1,5,a)+
     &                          gammagh(5,2,b)*gammagg(2,5,a)+
     &                          gammagh(5,3,b)*gammagg(3,5,a)+
     &                          gammagh(5,4,b)*gammagg(4,5,a)+
     &                          gammagh(5,5,b)*gammagg(5,5,a)
     &                       +
     &                          gammagh(1,1,b)*gammahh(1,1,a)+
     &                          gammagh(1,2,b)*gammahh(2,1,a)+
     &                          gammagh(1,3,b)*gammahh(3,1,a)+
     &                          gammagh(1,4,b)*gammahh(4,1,a)+
     &                          gammagh(1,5,b)*gammahh(5,1,a)+
     &                          gammagh(2,1,b)*gammahh(1,2,a)+
     &                          gammagh(2,2,b)*gammahh(2,2,a)+
     &                          gammagh(2,3,b)*gammahh(3,2,a)+
     &                          gammagh(2,4,b)*gammahh(4,2,a)+
     &                          gammagh(2,5,b)*gammahh(5,2,a)+
     &                          gammagh(3,1,b)*gammahh(1,3,a)+
     &                          gammagh(3,2,b)*gammahh(2,3,a)+
     &                          gammagh(3,3,b)*gammahh(3,3,a)+
     &                          gammagh(3,4,b)*gammahh(4,3,a)+
     &                          gammagh(3,5,b)*gammahh(5,3,a)+
     &                          gammagh(4,1,b)*gammahh(1,4,a)+
     &                          gammagh(4,2,b)*gammahh(2,4,a)+
     &                          gammagh(4,3,b)*gammahh(3,4,a)+
     &                          gammagh(4,4,b)*gammahh(4,4,a)+
     &                          gammagh(4,5,b)*gammahh(5,4,a)+
     &                          gammagh(5,1,b)*gammahh(1,5,a)+
     &                          gammagh(5,2,b)*gammahh(2,5,a)+
     &                          gammagh(5,3,b)*gammahh(3,5,a)+
     &                          gammagh(5,4,b)*gammahh(4,5,a)+
     &                          gammagh(5,5,b)*gammahh(5,5,a)
     &                       +
     &                          gammagh(1,1,b)*gammagh(1,1,a)+
     &                          gammagh(1,2,b)*gammagh(2,1,a)+
     &                          gammagh(1,3,b)*gammagh(3,1,a)+
     &                          gammagh(1,4,b)*gammagh(4,1,a)+
     &                          gammagh(1,5,b)*gammagh(5,1,a)+
     &                          gammagh(2,1,b)*gammagh(1,2,a)+
     &                          gammagh(2,2,b)*gammagh(2,2,a)+
     &                          gammagh(2,3,b)*gammagh(3,2,a)+
     &                          gammagh(2,4,b)*gammagh(4,2,a)+
     &                          gammagh(2,5,b)*gammagh(5,2,a)+
     &                          gammagh(3,1,b)*gammagh(1,3,a)+
     &                          gammagh(3,2,b)*gammagh(2,3,a)+
     &                          gammagh(3,3,b)*gammagh(3,3,a)+
     &                          gammagh(3,4,b)*gammagh(4,3,a)+
     &                          gammagh(3,5,b)*gammagh(5,3,a)+
     &                          gammagh(4,1,b)*gammagh(1,4,a)+
     &                          gammagh(4,2,b)*gammagh(2,4,a)+
     &                          gammagh(4,3,b)*gammagh(3,4,a)+
     &                          gammagh(4,4,b)*gammagh(4,4,a)+
     &                          gammagh(4,5,b)*gammagh(5,4,a)+
     &                          gammagh(5,1,b)*gammagh(1,5,a)+
     &                          gammagh(5,2,b)*gammagh(2,5,a)+
     &                          gammagh(5,3,b)*gammagh(3,5,a)+
     &                          gammagh(5,4,b)*gammagh(4,5,a)+
     &                          gammagh(5,5,b)*gammagh(5,5,a)
     &                       +
     &                          gammagh(1,1,b)*gammahg(1,1,a)+
     &                          gammagh(1,2,b)*gammahg(2,1,a)+
     &                          gammagh(1,3,b)*gammahg(3,1,a)+
     &                          gammagh(1,4,b)*gammahg(4,1,a)+
     &                          gammagh(1,5,b)*gammahg(5,1,a)+
     &                          gammagh(2,1,b)*gammahg(1,2,a)+
     &                          gammagh(2,2,b)*gammahg(2,2,a)+
     &                          gammagh(2,3,b)*gammahg(3,2,a)+
     &                          gammagh(2,4,b)*gammahg(4,2,a)+
     &                          gammagh(2,5,b)*gammahg(5,2,a)+
     &                          gammagh(3,1,b)*gammahg(1,3,a)+
     &                          gammagh(3,2,b)*gammahg(2,3,a)+
     &                          gammagh(3,3,b)*gammahg(3,3,a)+
     &                          gammagh(3,4,b)*gammahg(4,3,a)+
     &                          gammagh(3,5,b)*gammahg(5,3,a)+
     &                          gammagh(4,1,b)*gammahg(1,4,a)+
     &                          gammagh(4,2,b)*gammahg(2,4,a)+
     &                          gammagh(4,3,b)*gammahg(3,4,a)+
     &                          gammagh(4,4,b)*gammahg(4,4,a)+
     &                          gammagh(4,5,b)*gammahg(5,4,a)+
     &                          gammagh(5,1,b)*gammahg(1,5,a)+
     &                          gammagh(5,2,b)*gammahg(2,5,a)+
     &                          gammagh(5,3,b)*gammahg(3,5,a)+
     &                          gammagh(5,4,b)*gammahg(4,5,a)+
     &                          gammagh(5,5,b)*gammahg(5,5,a)
     &                       +
     &                          gammahg(1,1,b)*gammagg(1,1,a)+
     &                          gammahg(1,2,b)*gammagg(2,1,a)+
     &                          gammahg(1,3,b)*gammagg(3,1,a)+
     &                          gammahg(1,4,b)*gammagg(4,1,a)+
     &                          gammahg(1,5,b)*gammagg(5,1,a)+
     &                          gammahg(2,1,b)*gammagg(1,2,a)+
     &                          gammahg(2,2,b)*gammagg(2,2,a)+
     &                          gammahg(2,3,b)*gammagg(3,2,a)+
     &                          gammahg(2,4,b)*gammagg(4,2,a)+
     &                          gammahg(2,5,b)*gammagg(5,2,a)+
     &                          gammahg(3,1,b)*gammagg(1,3,a)+
     &                          gammahg(3,2,b)*gammagg(2,3,a)+
     &                          gammahg(3,3,b)*gammagg(3,3,a)+
     &                          gammahg(3,4,b)*gammagg(4,3,a)+
     &                          gammahg(3,5,b)*gammagg(5,3,a)+
     &                          gammahg(4,1,b)*gammagg(1,4,a)+
     &                          gammahg(4,2,b)*gammagg(2,4,a)+
     &                          gammahg(4,3,b)*gammagg(3,4,a)+
     &                          gammahg(4,4,b)*gammagg(4,4,a)+
     &                          gammahg(4,5,b)*gammagg(5,4,a)+
     &                          gammahg(5,1,b)*gammagg(1,5,a)+
     &                          gammahg(5,2,b)*gammagg(2,5,a)+
     &                          gammahg(5,3,b)*gammagg(3,5,a)+
     &                          gammahg(5,4,b)*gammagg(4,5,a)+
     &                          gammahg(5,5,b)*gammagg(5,5,a)
     &                       +
     &                          gammahg(1,1,b)*gammahh(1,1,a)+
     &                          gammahg(1,2,b)*gammahh(2,1,a)+
     &                          gammahg(1,3,b)*gammahh(3,1,a)+
     &                          gammahg(1,4,b)*gammahh(4,1,a)+
     &                          gammahg(1,5,b)*gammahh(5,1,a)+
     &                          gammahg(2,1,b)*gammahh(1,2,a)+
     &                          gammahg(2,2,b)*gammahh(2,2,a)+
     &                          gammahg(2,3,b)*gammahh(3,2,a)+
     &                          gammahg(2,4,b)*gammahh(4,2,a)+
     &                          gammahg(2,5,b)*gammahh(5,2,a)+
     &                          gammahg(3,1,b)*gammahh(1,3,a)+
     &                          gammahg(3,2,b)*gammahh(2,3,a)+
     &                          gammahg(3,3,b)*gammahh(3,3,a)+
     &                          gammahg(3,4,b)*gammahh(4,3,a)+
     &                          gammahg(3,5,b)*gammahh(5,3,a)+
     &                          gammahg(4,1,b)*gammahh(1,4,a)+
     &                          gammahg(4,2,b)*gammahh(2,4,a)+
     &                          gammahg(4,3,b)*gammahh(3,4,a)+
     &                          gammahg(4,4,b)*gammahh(4,4,a)+
     &                          gammahg(4,5,b)*gammahh(5,4,a)+
     &                          gammahg(5,1,b)*gammahh(1,5,a)+
     &                          gammahg(5,2,b)*gammahh(2,5,a)+
     &                          gammahg(5,3,b)*gammahh(3,5,a)+
     &                          gammahg(5,4,b)*gammahh(4,5,a)+
     &                          gammahg(5,5,b)*gammahh(5,5,a)
     &                       +
     &                          gammahg(1,1,b)*gammagh(1,1,a)+
     &                          gammahg(1,2,b)*gammagh(2,1,a)+
     &                          gammahg(1,3,b)*gammagh(3,1,a)+
     &                          gammahg(1,4,b)*gammagh(4,1,a)+
     &                          gammahg(1,5,b)*gammagh(5,1,a)+
     &                          gammahg(2,1,b)*gammagh(1,2,a)+
     &                          gammahg(2,2,b)*gammagh(2,2,a)+
     &                          gammahg(2,3,b)*gammagh(3,2,a)+
     &                          gammahg(2,4,b)*gammagh(4,2,a)+
     &                          gammahg(2,5,b)*gammagh(5,2,a)+
     &                          gammahg(3,1,b)*gammagh(1,3,a)+
     &                          gammahg(3,2,b)*gammagh(2,3,a)+
     &                          gammahg(3,3,b)*gammagh(3,3,a)+
     &                          gammahg(3,4,b)*gammagh(4,3,a)+
     &                          gammahg(3,5,b)*gammagh(5,3,a)+
     &                          gammahg(4,1,b)*gammagh(1,4,a)+
     &                          gammahg(4,2,b)*gammagh(2,4,a)+
     &                          gammahg(4,3,b)*gammagh(3,4,a)+
     &                          gammahg(4,4,b)*gammagh(4,4,a)+
     &                          gammahg(4,5,b)*gammagh(5,4,a)+
     &                          gammahg(5,1,b)*gammagh(1,5,a)+
     &                          gammahg(5,2,b)*gammagh(2,5,a)+
     &                          gammahg(5,3,b)*gammagh(3,5,a)+
     &                          gammahg(5,4,b)*gammagh(4,5,a)+
     &                          gammahg(5,5,b)*gammagh(5,5,a)
     &                       +
     &                          gammahg(1,1,b)*gammahg(1,1,a)+
     &                          gammahg(1,2,b)*gammahg(2,1,a)+
     &                          gammahg(1,3,b)*gammahg(3,1,a)+
     &                          gammahg(1,4,b)*gammahg(4,1,a)+
     &                          gammahg(1,5,b)*gammahg(5,1,a)+
     &                          gammahg(2,1,b)*gammahg(1,2,a)+
     &                          gammahg(2,2,b)*gammahg(2,2,a)+
     &                          gammahg(2,3,b)*gammahg(3,2,a)+
     &                          gammahg(2,4,b)*gammahg(4,2,a)+
     &                          gammahg(2,5,b)*gammahg(5,2,a)+
     &                          gammahg(3,1,b)*gammahg(1,3,a)+
     &                          gammahg(3,2,b)*gammahg(2,3,a)+
     &                          gammahg(3,3,b)*gammahg(3,3,a)+
     &                          gammahg(3,4,b)*gammahg(4,3,a)+
     &                          gammahg(3,5,b)*gammahg(5,3,a)+
     &                          gammahg(4,1,b)*gammahg(1,4,a)+
     &                          gammahg(4,2,b)*gammahg(2,4,a)+
     &                          gammahg(4,3,b)*gammahg(3,4,a)+
     &                          gammahg(4,4,b)*gammahg(4,4,a)+
     &                          gammahg(4,5,b)*gammahg(5,4,a)+
     &                          gammahg(5,1,b)*gammahg(1,5,a)+
     &                          gammahg(5,2,b)*gammahg(2,5,a)+
     &                          gammahg(5,3,b)*gammahg(3,5,a)+
     &                          gammahg(5,4,b)*gammahg(4,5,a)+
     &                          gammahg(5,5,b)*gammahg(5,5,a)
     &                            )
                  term8(a,b)=-2*lambda5*h0_ll(a,b)/3
     &
                  efe(a,b)=term1(a,b)+term2(a,b)+term3(a,b)+term4(a,b)
     &                    +term5(a,b)+term6(a,b)+term7(a,b)+term8(a,b)
     &                    -8*PI*(set_ll(a,b)-tr_set*g0_ll(a,b)/3)

                end do
              end do

              !--------------------------------------------------------------------------
              ! phi1_res = g^ab phi1,ab + g^ab,a phi1,b + g^cb gamma^a_ab phi1,c  
              !         (= g^ab phi1,ab - g^ab gamma^c_ab phi1,c) 
              !--------------------------------------------------------------------------
              phi1_res= phi10_xx(1,1)*g0_uu(1,1)+
     &                  phi10_xx(2,2)*g0_uu(2,2)+
     &                  phi10_xx(3,3)*g0_uu(3,3)+
     &                  phi10_xx(4,4)*g0_uu(4,4)+
     &                  phi10_xx(5,5)*g0_uu(5,5)+
     &               2*(phi10_xx(1,2)*g0_uu(1,2)+
     &                  phi10_xx(1,3)*g0_uu(1,3)+
     &                  phi10_xx(1,4)*g0_uu(1,4)+
     &                  phi10_xx(1,5)*g0_uu(1,5)+
     &                  phi10_xx(2,3)*g0_uu(2,3)+
     &                  phi10_xx(2,4)*g0_uu(2,4)+
     &                  phi10_xx(2,5)*g0_uu(2,5)+
     &                  phi10_xx(3,4)*g0_uu(3,4)+
     &                  phi10_xx(3,5)*g0_uu(3,5)+
     &                  phi10_xx(4,5)*g0_uu(4,5))
     &              +
     &                  phi10_x(1)*g0_uu_x(1,1,1)+
     &                  phi10_x(1)*g0_uu_x(2,1,2)+
     &                  phi10_x(1)*g0_uu_x(3,1,3)+
     &                  phi10_x(1)*g0_uu_x(4,1,4)+
     &                  phi10_x(1)*g0_uu_x(5,1,5)+
     &                  phi10_x(2)*g0_uu_x(1,2,1)+
     &                  phi10_x(2)*g0_uu_x(2,2,2)+
     &                  phi10_x(2)*g0_uu_x(3,2,3)+
     &                  phi10_x(2)*g0_uu_x(4,2,4)+
     &                  phi10_x(2)*g0_uu_x(5,2,5)+
     &                  phi10_x(3)*g0_uu_x(1,3,1)+
     &                  phi10_x(3)*g0_uu_x(2,3,2)+
     &                  phi10_x(3)*g0_uu_x(3,3,3)+
     &                  phi10_x(3)*g0_uu_x(4,3,4)+
     &                  phi10_x(3)*g0_uu_x(5,3,5)+
     &                  phi10_x(4)*g0_uu_x(1,4,1)+
     &                  phi10_x(4)*g0_uu_x(2,4,2)+
     &                  phi10_x(4)*g0_uu_x(3,4,3)+
     &                  phi10_x(4)*g0_uu_x(4,4,4)+
     &                  phi10_x(4)*g0_uu_x(5,4,5)+
     &                  phi10_x(5)*g0_uu_x(1,5,1)+
     &                  phi10_x(5)*g0_uu_x(2,5,2)+
     &                  phi10_x(5)*g0_uu_x(3,5,3)+
     &                  phi10_x(5)*g0_uu_x(4,5,4)+
     &                  phi10_x(5)*g0_uu_x(5,5,5)
     &              +
     &                  phi10_x(1)*g0_uu(1,1)*gamma_ull(1,1,1)+
     &                  phi10_x(1)*g0_uu(1,1)*gamma_ull(2,2,1)+ 
     &                  phi10_x(1)*g0_uu(1,1)*gamma_ull(3,3,1)+
     &                  phi10_x(1)*g0_uu(1,1)*gamma_ull(4,4,1)+
     &                  phi10_x(1)*g0_uu(1,1)*gamma_ull(5,5,1)+
     &                  phi10_x(1)*g0_uu(1,2)*gamma_ull(1,1,2)+
     &                  phi10_x(1)*g0_uu(1,2)*gamma_ull(2,2,2)+
     &                  phi10_x(1)*g0_uu(1,2)*gamma_ull(3,3,2)+
     &                  phi10_x(1)*g0_uu(1,2)*gamma_ull(4,4,2)+
     &                  phi10_x(1)*g0_uu(1,2)*gamma_ull(5,5,2)+
     &                  phi10_x(1)*g0_uu(1,3)*gamma_ull(1,1,3)+
     &                  phi10_x(1)*g0_uu(1,3)*gamma_ull(2,2,3)+
     &                  phi10_x(1)*g0_uu(1,3)*gamma_ull(3,3,3)+
     &                  phi10_x(1)*g0_uu(1,3)*gamma_ull(4,4,3)+
     &                  phi10_x(1)*g0_uu(1,3)*gamma_ull(5,5,3)+
     &                  phi10_x(1)*g0_uu(1,4)*gamma_ull(1,1,4)+
     &                  phi10_x(1)*g0_uu(1,4)*gamma_ull(2,2,4)+
     &                  phi10_x(1)*g0_uu(1,4)*gamma_ull(3,3,4)+
     &                  phi10_x(1)*g0_uu(1,4)*gamma_ull(4,4,4)+
     &                  phi10_x(1)*g0_uu(1,4)*gamma_ull(5,5,4)+
     &                  phi10_x(1)*g0_uu(1,5)*gamma_ull(1,1,5)+
     &                  phi10_x(1)*g0_uu(1,5)*gamma_ull(2,2,5)+
     &                  phi10_x(1)*g0_uu(1,5)*gamma_ull(3,3,5)+
     &                  phi10_x(1)*g0_uu(1,5)*gamma_ull(4,4,5)+
     &                  phi10_x(1)*g0_uu(1,5)*gamma_ull(5,5,5)+
     &                  phi10_x(2)*g0_uu(2,1)*gamma_ull(1,1,1)+
     &                  phi10_x(2)*g0_uu(2,1)*gamma_ull(2,2,1)+ 
     &                  phi10_x(2)*g0_uu(2,1)*gamma_ull(3,3,1)+
     &                  phi10_x(2)*g0_uu(2,1)*gamma_ull(4,4,1)+
     &                  phi10_x(2)*g0_uu(2,1)*gamma_ull(5,5,1)+
     &                  phi10_x(2)*g0_uu(2,2)*gamma_ull(1,1,2)+
     &                  phi10_x(2)*g0_uu(2,2)*gamma_ull(2,2,2)+
     &                  phi10_x(2)*g0_uu(2,2)*gamma_ull(3,3,2)+
     &                  phi10_x(2)*g0_uu(2,2)*gamma_ull(4,4,2)+
     &                  phi10_x(2)*g0_uu(2,2)*gamma_ull(5,5,2)+
     &                  phi10_x(2)*g0_uu(2,3)*gamma_ull(1,1,3)+
     &                  phi10_x(2)*g0_uu(2,3)*gamma_ull(2,2,3)+
     &                  phi10_x(2)*g0_uu(2,3)*gamma_ull(3,3,3)+
     &                  phi10_x(2)*g0_uu(2,3)*gamma_ull(4,4,3)+
     &                  phi10_x(2)*g0_uu(2,3)*gamma_ull(5,5,3)+
     &                  phi10_x(2)*g0_uu(2,4)*gamma_ull(1,1,4)+
     &                  phi10_x(2)*g0_uu(2,4)*gamma_ull(2,2,4)+
     &                  phi10_x(2)*g0_uu(2,4)*gamma_ull(3,3,4)+
     &                  phi10_x(2)*g0_uu(2,4)*gamma_ull(4,4,4)+
     &                  phi10_x(2)*g0_uu(2,4)*gamma_ull(5,5,4)+
     &                  phi10_x(2)*g0_uu(2,5)*gamma_ull(1,1,5)+
     &                  phi10_x(2)*g0_uu(2,5)*gamma_ull(2,2,5)+
     &                  phi10_x(2)*g0_uu(2,5)*gamma_ull(3,3,5)+
     &                  phi10_x(2)*g0_uu(2,5)*gamma_ull(4,4,5)+
     &                  phi10_x(2)*g0_uu(2,5)*gamma_ull(5,5,5)+
     &                  phi10_x(3)*g0_uu(3,1)*gamma_ull(1,1,1)+
     &                  phi10_x(3)*g0_uu(3,1)*gamma_ull(2,2,1)+ 
     &                  phi10_x(3)*g0_uu(3,1)*gamma_ull(3,3,1)+
     &                  phi10_x(3)*g0_uu(3,1)*gamma_ull(4,4,1)+
     &                  phi10_x(3)*g0_uu(3,1)*gamma_ull(5,5,1)+
     &                  phi10_x(3)*g0_uu(3,2)*gamma_ull(1,1,2)+
     &                  phi10_x(3)*g0_uu(3,2)*gamma_ull(2,2,2)+
     &                  phi10_x(3)*g0_uu(3,2)*gamma_ull(3,3,2)+
     &                  phi10_x(3)*g0_uu(3,2)*gamma_ull(4,4,2)+
     &                  phi10_x(3)*g0_uu(3,2)*gamma_ull(5,5,2)+
     &                  phi10_x(3)*g0_uu(3,3)*gamma_ull(1,1,3)+
     &                  phi10_x(3)*g0_uu(3,3)*gamma_ull(2,2,3)+
     &                  phi10_x(3)*g0_uu(3,3)*gamma_ull(3,3,3)+
     &                  phi10_x(3)*g0_uu(3,3)*gamma_ull(4,4,3)+
     &                  phi10_x(3)*g0_uu(3,3)*gamma_ull(5,5,3)+
     &                  phi10_x(3)*g0_uu(3,4)*gamma_ull(1,1,4)+
     &                  phi10_x(3)*g0_uu(3,4)*gamma_ull(2,2,4)+
     &                  phi10_x(3)*g0_uu(3,4)*gamma_ull(3,3,4)+
     &                  phi10_x(3)*g0_uu(3,4)*gamma_ull(4,4,4)+
     &                  phi10_x(3)*g0_uu(3,4)*gamma_ull(5,5,4)+
     &                  phi10_x(3)*g0_uu(3,5)*gamma_ull(1,1,5)+
     &                  phi10_x(3)*g0_uu(3,5)*gamma_ull(2,2,5)+
     &                  phi10_x(3)*g0_uu(3,5)*gamma_ull(3,3,5)+
     &                  phi10_x(3)*g0_uu(3,5)*gamma_ull(4,4,5)+
     &                  phi10_x(3)*g0_uu(3,5)*gamma_ull(5,5,5)+
     &                  phi10_x(4)*g0_uu(4,1)*gamma_ull(1,1,1)+
     &                  phi10_x(4)*g0_uu(4,1)*gamma_ull(2,2,1)+ 
     &                  phi10_x(4)*g0_uu(4,1)*gamma_ull(3,3,1)+
     &                  phi10_x(4)*g0_uu(4,1)*gamma_ull(4,4,1)+
     &                  phi10_x(4)*g0_uu(4,1)*gamma_ull(5,5,1)+
     &                  phi10_x(4)*g0_uu(4,2)*gamma_ull(1,1,2)+
     &                  phi10_x(4)*g0_uu(4,2)*gamma_ull(2,2,2)+
     &                  phi10_x(4)*g0_uu(4,2)*gamma_ull(3,3,2)+
     &                  phi10_x(4)*g0_uu(4,2)*gamma_ull(4,4,2)+
     &                  phi10_x(4)*g0_uu(4,2)*gamma_ull(5,5,2)+
     &                  phi10_x(4)*g0_uu(4,3)*gamma_ull(1,1,3)+
     &                  phi10_x(4)*g0_uu(4,3)*gamma_ull(2,2,3)+
     &                  phi10_x(4)*g0_uu(4,3)*gamma_ull(3,3,3)+
     &                  phi10_x(4)*g0_uu(4,3)*gamma_ull(4,4,3)+
     &                  phi10_x(4)*g0_uu(4,3)*gamma_ull(5,5,3)+
     &                  phi10_x(4)*g0_uu(4,4)*gamma_ull(1,1,4)+
     &                  phi10_x(4)*g0_uu(4,4)*gamma_ull(2,2,4)+
     &                  phi10_x(4)*g0_uu(4,4)*gamma_ull(3,3,4)+
     &                  phi10_x(4)*g0_uu(4,4)*gamma_ull(4,4,4)+
     &                  phi10_x(4)*g0_uu(4,4)*gamma_ull(5,5,4)+
     &                  phi10_x(4)*g0_uu(4,5)*gamma_ull(1,1,5)+
     &                  phi10_x(4)*g0_uu(4,5)*gamma_ull(2,2,5)+
     &                  phi10_x(4)*g0_uu(4,5)*gamma_ull(3,3,5)+
     &                  phi10_x(4)*g0_uu(4,5)*gamma_ull(4,4,5)+
     &                  phi10_x(4)*g0_uu(4,5)*gamma_ull(5,5,5)+
     &                  phi10_x(5)*g0_uu(5,1)*gamma_ull(1,1,1)+
     &                  phi10_x(5)*g0_uu(5,1)*gamma_ull(2,2,1)+ 
     &                  phi10_x(5)*g0_uu(5,1)*gamma_ull(3,3,1)+
     &                  phi10_x(5)*g0_uu(5,1)*gamma_ull(4,4,1)+
     &                  phi10_x(5)*g0_uu(5,1)*gamma_ull(5,5,1)+
     &                  phi10_x(5)*g0_uu(5,2)*gamma_ull(1,1,2)+
     &                  phi10_x(5)*g0_uu(5,2)*gamma_ull(2,2,2)+
     &                  phi10_x(5)*g0_uu(5,2)*gamma_ull(3,3,2)+
     &                  phi10_x(5)*g0_uu(5,2)*gamma_ull(4,4,2)+
     &                  phi10_x(5)*g0_uu(5,2)*gamma_ull(5,5,2)+
     &                  phi10_x(5)*g0_uu(5,3)*gamma_ull(1,1,3)+
     &                  phi10_x(5)*g0_uu(5,3)*gamma_ull(2,2,3)+
     &                  phi10_x(5)*g0_uu(5,3)*gamma_ull(3,3,3)+
     &                  phi10_x(5)*g0_uu(5,3)*gamma_ull(4,4,3)+
     &                  phi10_x(5)*g0_uu(5,3)*gamma_ull(5,5,3)+
     &                  phi10_x(5)*g0_uu(5,4)*gamma_ull(1,1,4)+
     &                  phi10_x(5)*g0_uu(5,4)*gamma_ull(2,2,4)+
     &                  phi10_x(5)*g0_uu(5,4)*gamma_ull(3,3,4)+
     &                  phi10_x(5)*g0_uu(5,4)*gamma_ull(4,4,4)+
     &                  phi10_x(5)*g0_uu(5,4)*gamma_ull(5,5,4)+
     &                  phi10_x(5)*g0_uu(5,5)*gamma_ull(1,1,5)+
     &                  phi10_x(5)*g0_uu(5,5)*gamma_ull(2,2,5)+
     &                  phi10_x(5)*g0_uu(5,5)*gamma_ull(3,3,5)+
     &                  phi10_x(5)*g0_uu(5,5)*gamma_ull(4,4,5)+
     &                  phi10_x(5)*g0_uu(5,5)*gamma_ull(5,5,5)

              !---------------------------------------------------------------- 
              ! computes diag. Jacobian of g_np1->L.g_np1 transformation
              ! by differentiating L.g wrt. g(a,b)_ij_np1 diag. entries
              ! 
              ! ddgb_J_tx differs from ddgb_J due to forward/backward stencils
              ! at excision surfaces that affect the cross-derivative tx
              ! (this is the only contribution, since the diag Jacobian is diff wrt. g_ij_np1)
              !---------------------------------------------------------------- 
              dgb_J=1/2/dt
              ddgb_J=1/dt/dt

              if (i.eq.1.or.(chr(i-1).eq.ex)) then
                 if (i.le.(Nx-3)
     &               .and.((chr(i+1).ne.ex
     &               .and.chr(i+2).ne.ex
     &               .and.chr(i+3).ne.ex))) then
                    ddgb_J_tx=-1/dt/dx
                 else if (i.le.(Nx-2)
     &                    .and.((chr(i+1).ne.ex
     &                    .and.chr(i+2).ne.ex))) then
                    ddgb_J_tx=-3/4/dt/dx
                 else if (i.le.(Nx-1).and.chr(i+1).ne.ex) then
                    ddgb_J_tx=-1/2/dt/dx
                 else
                    write(*,*) 'g_evo_opt: error in chr stencil (A)'
                    write(*,*) '    i,Nx,dx=',i,Nx,dx
                    write(*,*) '    (first error only)'
                    ddgb_J_tx=0
                 end if
              else if (i.eq.Nx.or.(chr(i+1).eq.ex)) then
                 if (i.ge.4
     &               .and.((chr(i-1).ne.ex
     &               .and.chr(i-2).ne.ex
     &               .and.chr(i-3).ne.ex))) then
                    ddgb_J_tx=1/dt/dx
                 else if (i.ge.3
     &                    .and.((chr(i-1).ne.ex
     &                    .and.chr(i-2).ne.ex))) then
                    ddgb_J_tx=3/4/dt/dx
                 else if (i.ge.2.and.chr(i-1).ne.ex) then
                    ddgb_J_tx=1/2/dt/dx
                 else
                    write(*,*) 'g_evo_opt: error in chr stencil (B)'
                    write(*,*) '    i,Nx,dx=',i,Nx,dx
                    write(*,*) '    (first error only)'
                    ddgb_J_tx=0
                 end if
              else
                 if ((chr(i+1).ne.ex.and.chr(i-1).ne.ex)) then
                    ddgb_J_tx=0
                 else
                    write(*,*) 'g_evo_opt: error in chr stencil (C)'
                    write(*,*) '    i,Nx,dx=',i,Nx,dx
                    write(*,*) '    (first error only)'
                    ddgb_J_tx=0
                 end if
              end if

              efe_J(1,1)=    -0.5d0*(
     &                          g0_uu(1,1)*ddgb_J
     &                                    *(1-x0**2)
     &                          -4*x0*g0_uu(1,2)*dgb_J
     &                          +2*g0_uu(1,2)*ddgb_J_tx
     &                                    *(1-x0**2)
     &                              )
     &                    
     &                       -0.5d0*(
     &                          -dgb_J*(1-x0**2)*
     &                          (g0_uu(1,1)*g0_uu(1,1)*g0_ll_x(1,1,1)+
     &                           g0_uu(1,1)*g0_uu(2,1)*g0_ll_x(1,1,2)+
     &                           g0_uu(1,1)*g0_uu(3,1)*g0_ll_x(1,1,3)+
     &                           g0_uu(1,1)*g0_uu(4,1)*g0_ll_x(1,1,4)+
     &                           g0_uu(1,1)*g0_uu(5,1)*g0_ll_x(1,1,5)+
     &                           g0_uu(2,1)*g0_uu(1,1)*g0_ll_x(1,2,1)+
     &                           g0_uu(2,1)*g0_uu(2,1)*g0_ll_x(1,2,2)+
     &                           g0_uu(2,1)*g0_uu(3,1)*g0_ll_x(1,2,3)+
     &                           g0_uu(2,1)*g0_uu(4,1)*g0_ll_x(1,2,4)+
     &                           g0_uu(2,1)*g0_uu(5,1)*g0_ll_x(1,2,5)+
     &                           g0_uu(3,1)*g0_uu(1,1)*g0_ll_x(1,3,1)+
     &                           g0_uu(3,1)*g0_uu(2,1)*g0_ll_x(1,3,2)+
     &                           g0_uu(3,1)*g0_uu(3,1)*g0_ll_x(1,3,3)+
     &                           g0_uu(3,1)*g0_uu(4,1)*g0_ll_x(1,3,4)+
     &                           g0_uu(3,1)*g0_uu(5,1)*g0_ll_x(1,3,5)+
     &                           g0_uu(4,1)*g0_uu(1,1)*g0_ll_x(1,4,1)+
     &                           g0_uu(4,1)*g0_uu(2,1)*g0_ll_x(1,4,2)+
     &                           g0_uu(4,1)*g0_uu(3,1)*g0_ll_x(1,4,3)+
     &                           g0_uu(4,1)*g0_uu(4,1)*g0_ll_x(1,4,4)+
     &                           g0_uu(4,1)*g0_uu(5,1)*g0_ll_x(1,4,5)+
     &                           g0_uu(5,1)*g0_uu(1,1)*g0_ll_x(1,5,1)+
     &                           g0_uu(5,1)*g0_uu(2,1)*g0_ll_x(1,5,2)+
     &                           g0_uu(5,1)*g0_uu(3,1)*g0_ll_x(1,5,3)+
     &                           g0_uu(5,1)*g0_uu(4,1)*g0_ll_x(1,5,4)+
     &                           g0_uu(5,1)*g0_uu(5,1)*g0_ll_x(1,5,5))
     &                          +dgb_J*(1-x0**2)*
     &                          (g0_uu_x(1,1,1))
     &                              )*2
     &
     &                       +      (
     &                          0.5d0*dgb_J*(1-x0**2)*
     &                          ((Hads_l(1)+A_l(1))*g0_uu(1,1)+
     &                           (Hads_l(2)+A_l(2))*g0_uu(2,1)+
     &                           (Hads_l(3)+A_l(3))*g0_uu(3,1)+
     &                           (Hads_l(4)+A_l(4))*g0_uu(4,1)+
     &                           (Hads_l(5)+A_l(5))*g0_uu(5,1))
     &                              )  
     &                    
     &                       -      (
     &                          0.25d0*dgb_J*(1-x0**2)*
     &                          (cuuuu(1,1,1,1)*dlll(1,1,1)+
     &                           cuuuu(1,1,1,2)*dlll(1,1,2)+
     &                           cuuuu(1,1,1,3)*dlll(1,1,3)+
     &                           cuuuu(1,1,1,4)*dlll(1,1,4)+
     &                           cuuuu(1,1,1,5)*dlll(1,1,5)+
     &                           cuuuu(2,1,1,1)*dlll(2,1,1)+
     &                           cuuuu(2,1,1,2)*dlll(2,1,2)+
     &                           cuuuu(2,1,1,3)*dlll(2,1,3)+
     &                           cuuuu(2,1,1,4)*dlll(2,1,4)+
     &                           cuuuu(2,1,1,5)*dlll(2,1,5)+
     &                           cuuuu(3,1,1,1)*dlll(3,1,1)+
     &                           cuuuu(3,1,1,2)*dlll(3,1,2)+
     &                           cuuuu(3,1,1,3)*dlll(3,1,3)+
     &                           cuuuu(3,1,1,4)*dlll(3,1,4)+
     &                           cuuuu(3,1,1,5)*dlll(3,1,5)+
     &                           cuuuu(4,1,1,1)*dlll(4,1,1)+
     &                           cuuuu(4,1,1,2)*dlll(4,1,2)+
     &                           cuuuu(4,1,1,3)*dlll(4,1,3)+
     &                           cuuuu(4,1,1,4)*dlll(4,1,4)+
     &                           cuuuu(4,1,1,5)*dlll(4,1,5)+
     &                           cuuuu(5,1,1,1)*dlll(5,1,1)+
     &                           cuuuu(5,1,1,2)*dlll(5,1,2)+
     &                           cuuuu(5,1,1,3)*dlll(5,1,3)+
     &                           cuuuu(5,1,1,4)*dlll(5,1,4)+
     &                           cuuuu(5,1,1,5)*dlll(5,1,5)
     &                         +
     &                           cuuuu(1,1,1,1)*dlll(1,1,1)+
     &                           cuuuu(1,1,2,1)*dlll(2,1,1)+
     &                           cuuuu(1,1,3,1)*dlll(3,1,1)+
     &                           cuuuu(1,1,4,1)*dlll(4,1,1)+
     &                           cuuuu(1,1,5,1)*dlll(5,1,1)+
     &                           cuuuu(1,2,1,1)*dlll(1,1,2)+
     &                           cuuuu(1,2,2,1)*dlll(2,1,2)+
     &                           cuuuu(1,2,3,1)*dlll(3,1,2)+
     &                           cuuuu(1,2,4,1)*dlll(4,1,2)+
     &                           cuuuu(1,2,5,1)*dlll(5,1,2)+
     &                           cuuuu(1,3,1,1)*dlll(1,1,3)+
     &                           cuuuu(1,3,2,1)*dlll(2,1,3)+
     &                           cuuuu(1,3,3,1)*dlll(3,1,3)+
     &                           cuuuu(1,3,4,1)*dlll(4,1,3)+
     &                           cuuuu(1,3,5,1)*dlll(5,1,3)+
     &                           cuuuu(1,4,1,1)*dlll(1,1,4)+
     &                           cuuuu(1,4,2,1)*dlll(2,1,4)+
     &                           cuuuu(1,4,3,1)*dlll(3,1,4)+
     &                           cuuuu(1,4,4,1)*dlll(4,1,4)+
     &                           cuuuu(1,4,5,1)*dlll(5,1,4)+
     &                           cuuuu(1,5,1,1)*dlll(1,1,5)+
     &                           cuuuu(1,5,2,1)*dlll(2,1,5)+
     &                           cuuuu(1,5,3,1)*dlll(3,1,5)+
     &                           cuuuu(1,5,4,1)*dlll(4,1,5)+
     &                           cuuuu(1,5,5,1)*dlll(5,1,5))
     &                              )

              efe_J(1,2)=  -0.5d0*(
     &                        g0_uu(1,1)*ddgb_J
     &                                  *(1-x0**2)**2
     &                        -4*(2)*x0*g0_uu(1,2)*dgb_J
     &                                  *(1-x0**2)
     &                        +2*g0_uu(1,2)*ddgb_J_tx
     &                                  *(1-x0**2)**2
     &                            )
     &                     
     &                     -0.5d0*(
     &                        -dgb_J*(1-x0**2)**2*
     &                        (g0_uu(1,1)*g0_uu(1,2)*g0_ll_x(2,1,1)+
     &                         g0_uu(1,1)*g0_uu(2,2)*g0_ll_x(2,1,2)+
     &                         g0_uu(1,1)*g0_uu(3,2)*g0_ll_x(2,1,3)+
     &                         g0_uu(1,1)*g0_uu(4,2)*g0_ll_x(2,1,4)+
     &                         g0_uu(1,1)*g0_uu(5,2)*g0_ll_x(2,1,5)+
     &                         g0_uu(2,1)*g0_uu(1,2)*g0_ll_x(2,2,1)+
     &                         g0_uu(2,1)*g0_uu(2,2)*g0_ll_x(2,2,2)+
     &                         g0_uu(2,1)*g0_uu(3,2)*g0_ll_x(2,2,3)+
     &                         g0_uu(2,1)*g0_uu(4,2)*g0_ll_x(2,2,4)+
     &                         g0_uu(2,1)*g0_uu(5,2)*g0_ll_x(2,2,5)+
     &                         g0_uu(3,1)*g0_uu(1,2)*g0_ll_x(2,3,1)+
     &                         g0_uu(3,1)*g0_uu(2,2)*g0_ll_x(2,3,2)+
     &                         g0_uu(3,1)*g0_uu(3,2)*g0_ll_x(2,3,3)+
     &                         g0_uu(3,1)*g0_uu(4,2)*g0_ll_x(2,3,4)+
     &                         g0_uu(3,1)*g0_uu(5,2)*g0_ll_x(2,3,5)+
     &                         g0_uu(4,1)*g0_uu(1,2)*g0_ll_x(2,4,1)+
     &                         g0_uu(4,1)*g0_uu(2,2)*g0_ll_x(2,4,2)+
     &                         g0_uu(4,1)*g0_uu(3,2)*g0_ll_x(2,4,3)+
     &                         g0_uu(4,1)*g0_uu(4,2)*g0_ll_x(2,4,4)+
     &                         g0_uu(4,1)*g0_uu(5,2)*g0_ll_x(2,4,5)+
     &                         g0_uu(5,1)*g0_uu(1,2)*g0_ll_x(2,5,1)+
     &                         g0_uu(5,1)*g0_uu(2,2)*g0_ll_x(2,5,2)+
     &                         g0_uu(5,1)*g0_uu(3,2)*g0_ll_x(2,5,3)+
     &                         g0_uu(5,1)*g0_uu(4,2)*g0_ll_x(2,5,4)+
     &                         g0_uu(5,1)*g0_uu(5,2)*g0_ll_x(2,5,5))
     &                        +dgb_J*(1-x0**2)**2*
     &                        (g0_uu_x(1,1,1))
     &                            )
     &                  
     &                     -0.5d0*(
     &                        dgb_J*(1-x0**2)**2*
     &                        (g0_uu_x(2,1,2))
     &                            )
     &
     &                     -      (
     &                        0.5d0*dgb_J*(1-x0**2)**2*
     &                        (cuuuu(1,1,1,2)*dlll(1,2,1)+
     &                         cuuuu(1,1,2,2)*dlll(2,2,1)+
     &                         cuuuu(1,1,3,2)*dlll(3,2,1)+
     &                         cuuuu(1,1,4,2)*dlll(4,2,1)+
     &                         cuuuu(1,1,5,2)*dlll(5,2,1)+
     &                         cuuuu(1,2,1,2)*dlll(1,2,2)+
     &                         cuuuu(1,2,2,2)*dlll(2,2,2)+
     &                         cuuuu(1,2,3,2)*dlll(3,2,2)+
     &                         cuuuu(1,2,4,2)*dlll(4,2,2)+
     &                         cuuuu(1,2,5,2)*dlll(5,2,2)+
     &                         cuuuu(1,3,1,2)*dlll(1,2,3)+
     &                         cuuuu(1,3,2,2)*dlll(2,2,3)+
     &                         cuuuu(1,3,3,2)*dlll(3,2,3)+
     &                         cuuuu(1,3,4,2)*dlll(4,2,3)+
     &                         cuuuu(1,3,5,2)*dlll(5,2,3)+
     &                         cuuuu(1,4,1,2)*dlll(1,2,4)+
     &                         cuuuu(1,4,2,2)*dlll(2,2,4)+
     &                         cuuuu(1,4,3,2)*dlll(3,2,4)+
     &                         cuuuu(1,4,4,2)*dlll(4,2,4)+
     &                         cuuuu(1,4,5,2)*dlll(5,2,4)+
     &                         cuuuu(1,5,1,2)*dlll(1,2,5)+
     &                         cuuuu(1,5,2,2)*dlll(2,2,5)+
     &                         cuuuu(1,5,3,2)*dlll(3,2,5)+
     &                         cuuuu(1,5,4,2)*dlll(4,2,5)+
     &                         cuuuu(1,5,5,2)*dlll(5,2,5))
     &                            )

              efe_J(2,2)=-0.5d0*(
     &                      g0_uu(1,1)*ddgb_J
     &                                *(1-x0**2)
     &                      -4*x0*g0_uu(1,2)*dgb_J
     &                      +2*g0_uu(1,2)*ddgb_J_tx
     &                                *(1-x0**2)
     &                          )
     &                
     &                   -0.5d0*(
     &                      +dgb_J*(1-x0**2)*
     &                      (g0_uu_x(2,1,2))
     &                          )
     &                
     &                   -0.5d0*(
     &                      +dgb_J*(1-x0**2)*
     &                      (g0_uu_x(2,1,2))
     &                          )
     &
     &                   +      (
     &                      -0.5d0*dgb_J*(1-x0**2)*
     &                      ((Hads_l(1)+A_l(1))*g0_uu(1,1)+
     &                       (Hads_l(2)+A_l(2))*g0_uu(2,1)+
     &                       (Hads_l(3)+A_l(3))*g0_uu(3,1)+
     &                       (Hads_l(4)+A_l(4))*g0_uu(4,1)+
     &                       (Hads_l(5)+A_l(5))*g0_uu(5,1))
     &                          )  
     &                
     &                   -      (
     &                      0.25d0*dgb_J*(1-x0**2)*
     &                      (cuuuu(1,2,1,1)*dlll(1,2,1)+
     &                       cuuuu(1,2,1,2)*dlll(1,2,2)+
     &                       cuuuu(1,2,1,3)*dlll(1,2,3)+
     &                       cuuuu(1,2,1,4)*dlll(1,2,4)+
     &                       cuuuu(1,2,1,5)*dlll(1,2,5)+
     &                       cuuuu(2,2,1,1)*dlll(2,2,1)+
     &                       cuuuu(2,2,1,2)*dlll(2,2,2)+
     &                       cuuuu(2,2,1,3)*dlll(2,2,3)+
     &                       cuuuu(2,2,1,4)*dlll(2,2,4)+
     &                       cuuuu(2,2,1,5)*dlll(2,2,5)+
     &                       cuuuu(3,2,1,1)*dlll(3,2,1)+
     &                       cuuuu(3,2,1,2)*dlll(3,2,2)+
     &                       cuuuu(3,2,1,3)*dlll(3,2,3)+
     &                       cuuuu(3,2,1,4)*dlll(3,2,4)+
     &                       cuuuu(3,2,1,5)*dlll(3,2,5)+
     &                       cuuuu(4,2,1,1)*dlll(4,2,1)+
     &                       cuuuu(4,2,1,2)*dlll(4,2,2)+
     &                       cuuuu(4,2,1,3)*dlll(4,2,3)+
     &                       cuuuu(4,2,1,4)*dlll(4,2,4)+
     &                       cuuuu(4,2,1,5)*dlll(4,2,5)+
     &                       cuuuu(5,2,1,1)*dlll(5,2,1)+
     &                       cuuuu(5,2,1,2)*dlll(5,2,2)+
     &                       cuuuu(5,2,1,3)*dlll(5,2,3)+
     &                       cuuuu(5,2,1,4)*dlll(5,2,4)+
     &                       cuuuu(5,2,1,5)*dlll(5,2,5)-
     &                       cuuuu(1,1,2,1)*dlll(1,2,1)-
     &                       cuuuu(1,1,2,2)*dlll(1,2,2)-
     &                       cuuuu(1,1,2,3)*dlll(1,2,3)-
     &                       cuuuu(1,1,2,4)*dlll(1,2,4)-
     &                       cuuuu(1,1,2,5)*dlll(1,2,5)-
     &                       cuuuu(2,1,2,1)*dlll(2,2,1)-
     &                       cuuuu(2,1,2,2)*dlll(2,2,2)-
     &                       cuuuu(2,1,2,3)*dlll(2,2,3)-
     &                       cuuuu(2,1,2,4)*dlll(2,2,4)-
     &                       cuuuu(2,1,2,5)*dlll(2,2,5)-
     &                       cuuuu(3,1,2,1)*dlll(3,2,1)-
     &                       cuuuu(3,1,2,2)*dlll(3,2,2)-
     &                       cuuuu(3,1,2,3)*dlll(3,2,3)-
     &                       cuuuu(3,1,2,4)*dlll(3,2,4)-
     &                       cuuuu(3,1,2,5)*dlll(3,2,5)-
     &                       cuuuu(4,1,2,1)*dlll(4,2,1)-
     &                       cuuuu(4,1,2,2)*dlll(4,2,2)-
     &                       cuuuu(4,1,2,3)*dlll(4,2,3)-
     &                       cuuuu(4,1,2,4)*dlll(4,2,4)-
     &                       cuuuu(4,1,2,5)*dlll(4,2,5)-
     &                       cuuuu(5,1,2,1)*dlll(5,2,1)-
     &                       cuuuu(5,1,2,2)*dlll(5,2,2)-
     &                       cuuuu(5,1,2,3)*dlll(5,2,3)-
     &                       cuuuu(5,1,2,4)*dlll(5,2,4)-
     &                       cuuuu(5,1,2,5)*dlll(5,2,5)
     &                     +
     &                       cuuuu(1,1,1,2)*dlll(1,2,1)+
     &                       cuuuu(1,1,2,2)*dlll(2,2,1)+
     &                       cuuuu(1,1,3,2)*dlll(3,2,1)+
     &                       cuuuu(1,1,4,2)*dlll(4,2,1)+
     &                       cuuuu(1,1,5,2)*dlll(5,2,1)+
     &                       cuuuu(1,2,1,2)*dlll(1,2,2)+
     &                       cuuuu(1,2,2,2)*dlll(2,2,2)+
     &                       cuuuu(1,2,3,2)*dlll(3,2,2)+
     &                       cuuuu(1,2,4,2)*dlll(4,2,2)+
     &                       cuuuu(1,2,5,2)*dlll(5,2,2)+
     &                       cuuuu(1,3,1,2)*dlll(1,2,3)+
     &                       cuuuu(1,3,2,2)*dlll(2,2,3)+
     &                       cuuuu(1,3,3,2)*dlll(3,2,3)+
     &                       cuuuu(1,3,4,2)*dlll(4,2,3)+
     &                       cuuuu(1,3,5,2)*dlll(5,2,3)+
     &                       cuuuu(1,4,1,2)*dlll(1,2,4)+
     &                       cuuuu(1,4,2,2)*dlll(2,2,4)+
     &                       cuuuu(1,4,3,2)*dlll(3,2,4)+
     &                       cuuuu(1,4,4,2)*dlll(4,2,4)+
     &                       cuuuu(1,4,5,2)*dlll(5,2,4)+
     &                       cuuuu(1,5,1,2)*dlll(1,2,5)+
     &                       cuuuu(1,5,2,2)*dlll(2,2,5)+
     &                       cuuuu(1,5,3,2)*dlll(3,2,5)+
     &                       cuuuu(1,5,4,2)*dlll(4,2,5)+
     &                       cuuuu(1,5,5,2)*dlll(5,2,5)-
     &                       cuuuu(2,1,1,1)*dlll(1,2,1)-
     &                       cuuuu(2,1,2,1)*dlll(2,2,1)-
     &                       cuuuu(2,1,3,1)*dlll(3,2,1)-
     &                       cuuuu(2,1,4,1)*dlll(4,2,1)-
     &                       cuuuu(2,1,5,1)*dlll(5,2,1)-
     &                       cuuuu(2,2,1,1)*dlll(1,2,2)-
     &                       cuuuu(2,2,2,1)*dlll(2,2,2)-
     &                       cuuuu(2,2,3,1)*dlll(3,2,2)-
     &                       cuuuu(2,2,4,1)*dlll(4,2,2)-
     &                       cuuuu(2,2,5,1)*dlll(5,2,2)-
     &                       cuuuu(2,3,1,1)*dlll(1,2,3)-
     &                       cuuuu(2,3,2,1)*dlll(2,2,3)-
     &                       cuuuu(2,3,3,1)*dlll(3,2,3)-
     &                       cuuuu(2,3,4,1)*dlll(4,2,3)-
     &                       cuuuu(2,3,5,1)*dlll(5,2,3)-
     &                       cuuuu(2,4,1,1)*dlll(1,2,4)-
     &                       cuuuu(2,4,2,1)*dlll(2,2,4)-
     &                       cuuuu(2,4,3,1)*dlll(3,2,4)-
     &                       cuuuu(2,4,4,1)*dlll(4,2,4)-
     &                       cuuuu(2,4,5,1)*dlll(5,2,4)-
     &                       cuuuu(2,5,1,1)*dlll(1,2,5)-
     &                       cuuuu(2,5,2,1)*dlll(2,2,5)-
     &                       cuuuu(2,5,3,1)*dlll(3,2,5)-
     &                       cuuuu(2,5,4,1)*dlll(4,2,5)-
     &                       cuuuu(2,5,5,1)*dlll(5,2,5))
     &                          )

              efe_J(3,3)=    -0.5d0*(
     &                          x0**2*g0_uu(1,1)*ddgb_J
     &                                        *(1-x0**2)
     &                          +2*(2*x0-4*x0**3)*g0_uu(1,2)*dgb_J
     &                                        *(1-x0**2)
     &                          +2*x0**2*(1-x0**2)*g0_uu(1,2)*ddgb_J_tx
     &                              )
     &                    
     &                       -0.5d0*(
     &                          +x0**2*dgb_J
     &                                      *(g0_uu_x(3,1,3))
     &                                      *(1-x0**2)
     &                              )
     &                    
     &                       -0.5d0*(
     &                          +x0**2*dgb_J
     &                                      *(g0_uu_x(3,1,3))
     &                                      *(1-x0**2)
     &                              )
     &
     &                       +      (
     &                          -0.5d0*x0**2*dgb_J
     &                                      *(1-x0**2)*
     &                          ((Hads_l(1)+A_l(1))*g0_uu(1,1)+
     &                           (Hads_l(2)+A_l(2))*g0_uu(2,1)+
     &                           (Hads_l(3)+A_l(3))*g0_uu(3,1)+
     &                           (Hads_l(4)+A_l(4))*g0_uu(4,1)+
     &                           (Hads_l(5)+A_l(5))*g0_uu(5,1))
     &                              )  
     &                    
     &                       -      (
     &                          0.25d0*x0**2*dgb_J
     &                                      *(1-x0**2)*
     &                          (cuuuu(1,3,1,1)*dlll(1,3,1)+
     &                           cuuuu(1,3,1,2)*dlll(1,3,2)+
     &                           cuuuu(1,3,1,3)*dlll(1,3,3)+
     &                           cuuuu(1,3,1,4)*dlll(1,3,4)+
     &                           cuuuu(1,3,1,5)*dlll(1,3,5)+
     &                           cuuuu(2,3,1,1)*dlll(2,3,1)+
     &                           cuuuu(2,3,1,2)*dlll(2,3,2)+
     &                           cuuuu(2,3,1,3)*dlll(2,3,3)+
     &                           cuuuu(2,3,1,4)*dlll(2,3,4)+
     &                           cuuuu(2,3,1,5)*dlll(2,3,5)+
     &                           cuuuu(3,3,1,1)*dlll(3,3,1)+
     &                           cuuuu(3,3,1,2)*dlll(3,3,2)+
     &                           cuuuu(3,3,1,3)*dlll(3,3,3)+
     &                           cuuuu(3,3,1,4)*dlll(3,3,4)+
     &                           cuuuu(3,3,1,5)*dlll(3,3,5)+
     &                           cuuuu(4,3,1,1)*dlll(4,3,1)+
     &                           cuuuu(4,3,1,2)*dlll(4,3,2)+
     &                           cuuuu(4,3,1,3)*dlll(4,3,3)+
     &                           cuuuu(4,3,1,4)*dlll(4,3,4)+
     &                           cuuuu(4,3,1,5)*dlll(4,3,5)+
     &                           cuuuu(5,3,1,1)*dlll(5,3,1)+
     &                           cuuuu(5,3,1,2)*dlll(5,3,2)+
     &                           cuuuu(5,3,1,3)*dlll(5,3,3)+
     &                           cuuuu(5,3,1,4)*dlll(5,3,4)+
     &                           cuuuu(5,3,1,5)*dlll(5,3,5)-
     &                           cuuuu(1,1,3,1)*dlll(1,3,1)-
     &                           cuuuu(1,1,3,2)*dlll(1,3,2)-
     &                           cuuuu(1,1,3,3)*dlll(1,3,3)-
     &                           cuuuu(1,1,3,4)*dlll(1,3,4)-
     &                           cuuuu(1,1,3,5)*dlll(1,3,5)-
     &                           cuuuu(2,1,3,1)*dlll(2,3,1)-
     &                           cuuuu(2,1,3,2)*dlll(2,3,2)-
     &                           cuuuu(2,1,3,3)*dlll(2,3,3)-
     &                           cuuuu(2,1,3,4)*dlll(2,3,4)-
     &                           cuuuu(2,1,3,5)*dlll(2,3,5)-
     &                           cuuuu(3,1,3,1)*dlll(3,3,1)-
     &                           cuuuu(3,1,3,2)*dlll(3,3,2)-
     &                           cuuuu(3,1,3,3)*dlll(3,3,3)-
     &                           cuuuu(3,1,3,4)*dlll(3,3,4)-
     &                           cuuuu(3,1,3,5)*dlll(3,3,5)-
     &                           cuuuu(4,1,3,1)*dlll(4,3,1)-
     &                           cuuuu(4,1,3,2)*dlll(4,3,2)-
     &                           cuuuu(4,1,3,3)*dlll(4,3,3)-
     &                           cuuuu(4,1,3,4)*dlll(4,3,4)-
     &                           cuuuu(4,1,3,5)*dlll(4,3,5)-
     &                           cuuuu(5,1,3,1)*dlll(5,3,1)-
     &                           cuuuu(5,1,3,2)*dlll(5,3,2)-
     &                           cuuuu(5,1,3,3)*dlll(5,3,3)-
     &                           cuuuu(5,1,3,4)*dlll(5,3,4)-
     &                           cuuuu(5,1,3,5)*dlll(5,3,5)
     &                         +
     &                           cuuuu(1,1,1,3)*dlll(1,3,1)+
     &                           cuuuu(1,1,2,3)*dlll(2,3,1)+
     &                           cuuuu(1,1,3,3)*dlll(3,3,1)+
     &                           cuuuu(1,1,4,3)*dlll(4,3,1)+
     &                           cuuuu(1,1,5,3)*dlll(5,3,1)+
     &                           cuuuu(1,2,1,3)*dlll(1,3,2)+
     &                           cuuuu(1,2,2,3)*dlll(2,3,2)+
     &                           cuuuu(1,2,3,3)*dlll(3,3,2)+
     &                           cuuuu(1,2,4,3)*dlll(4,3,2)+
     &                           cuuuu(1,2,5,3)*dlll(5,3,2)+
     &                           cuuuu(1,3,1,3)*dlll(1,3,3)+
     &                           cuuuu(1,3,2,3)*dlll(2,3,3)+
     &                           cuuuu(1,3,3,3)*dlll(3,3,3)+
     &                           cuuuu(1,3,4,3)*dlll(4,3,3)+
     &                           cuuuu(1,3,5,3)*dlll(5,3,3)+
     &                           cuuuu(1,4,1,3)*dlll(1,3,4)+
     &                           cuuuu(1,4,2,3)*dlll(2,3,4)+
     &                           cuuuu(1,4,3,3)*dlll(3,3,4)+
     &                           cuuuu(1,4,4,3)*dlll(4,3,4)+
     &                           cuuuu(1,4,5,3)*dlll(5,3,4)+
     &                           cuuuu(1,5,1,3)*dlll(1,3,5)+
     &                           cuuuu(1,5,2,3)*dlll(2,3,5)+
     &                           cuuuu(1,5,3,3)*dlll(3,3,5)+
     &                           cuuuu(1,5,4,3)*dlll(4,3,5)+
     &                           cuuuu(1,5,5,3)*dlll(5,3,5)-
     &                           cuuuu(3,1,1,1)*dlll(1,3,1)-
     &                           cuuuu(3,1,2,1)*dlll(2,3,1)-
     &                           cuuuu(3,1,3,1)*dlll(3,3,1)-
     &                           cuuuu(3,1,4,1)*dlll(4,3,1)-
     &                           cuuuu(3,1,5,1)*dlll(5,3,1)-
     &                           cuuuu(3,2,1,1)*dlll(1,3,2)-
     &                           cuuuu(3,2,2,1)*dlll(2,3,2)-
     &                           cuuuu(3,2,3,1)*dlll(3,3,2)-
     &                           cuuuu(3,2,4,1)*dlll(4,3,2)-
     &                           cuuuu(3,2,5,1)*dlll(5,3,2)-
     &                           cuuuu(3,3,1,1)*dlll(1,3,3)-
     &                           cuuuu(3,3,2,1)*dlll(2,3,3)-
     &                           cuuuu(3,3,3,1)*dlll(3,3,3)-
     &                           cuuuu(3,3,4,1)*dlll(4,3,3)-
     &                           cuuuu(3,3,5,1)*dlll(5,3,3)-
     &                           cuuuu(3,4,1,1)*dlll(1,3,4)-
     &                           cuuuu(3,4,2,1)*dlll(2,3,4)-
     &                           cuuuu(3,4,3,1)*dlll(3,3,4)-
     &                           cuuuu(3,4,4,1)*dlll(4,3,4)-
     &                           cuuuu(3,4,5,1)*dlll(5,3,4)-
     &                           cuuuu(3,5,1,1)*dlll(1,3,5)-
     &                           cuuuu(3,5,2,1)*dlll(2,3,5)-
     &                           cuuuu(3,5,3,1)*dlll(3,3,5)-
     &                           cuuuu(3,5,4,1)*dlll(4,3,5)-
     &                           cuuuu(3,5,5,1)*dlll(5,3,5))
     &                              )


              !----------------------------------------------------------------
              ! computes diag. Jacobian of phi1_np1->L.phi1_np1 transformation
              ! by differentiating L.phi1=box.phi1-dV/dphi1 wrt. phi1_np1
              ! and remember: phi10=phi1*(1-x0**2)**3 
              ! 
              ! ddphi1_J_tx differs from ddphi_J due to forward/backward stencils
              ! at excision surfaces that affect the cross-derivative tx 
              ! (these are the only contributions, since the diag Jacobian is diff wrt. phi1_np1)
              !----------------------------------------------------------------
              dphi1_J=1/2/dt
              ddphi1_J=1/dt/dt

              if (i.eq.1.or.(chr(i-1).eq.ex)) then
                 if (i.le.(Nx-3)
     &               .and.((chr(i+1).ne.ex
     &               .and.chr(i+2).ne.ex
     &               .and.chr(i+3).ne.ex))) then
                    ddphi1_J_tx=-1/dt/dx
                 else if (i.le.(Nx-2)
     &                    .and.((chr(i+1).ne.ex
     &                    .and.chr(i+2).ne.ex))) then
                    ddphi1_J_tx=-3/4/dt/dx
                 else if (i.le.(Nx-1).and.chr(i+1).ne.ex) then
                    ddphi1_J_tx=-1/2/dt/dx
                 else
                    write(*,*) 'g_evo_opt: error in chr stencil (A)'
                    write(*,*) '    i,Nx,dx=',i,Nx,dx
                    write(*,*) '    (first error only)'
                    ddphi1_J_tx=0
                 end if
              else if (i.eq.Nx.or.(chr(i+1).eq.ex)) then
                 if (i.ge.4
     &               .and.((chr(i-1).ne.ex
     &               .and.chr(i-2).ne.ex
     &               .and.chr(i-3).ne.ex))) then
                    ddphi1_J_tx=1/dt/dx
                 else if (i.ge.3
     &                    .and.((chr(i-1).ne.ex
     &                    .and.chr(i-2).ne.ex))) then
                    ddphi1_J_tx=3/4/dt/dx
                 else if (i.ge.2.and.chr(i-1).ne.ex) then
                    ddphi1_J_tx=1/2/dt/dx
                 else
                    write(*,*) 'g_evo_opt: error in chr stencil (B)'
                    write(*,*) '    i,Nx,dx=',i,Nx,dx
                    write(*,*) '    (first error only)'
                    ddphi1_J_tx=0
                 end if
              else
                 if ((chr(i+1).ne.ex.and.chr(i-1).ne.ex)) then
                    ddphi1_J_tx=0
                 else
                    write(*,*) 'g_evo_opt: error in chr stencil (C)'
                    write(*,*) '    i,Nx,dx=',i,Nx,dx
                    write(*,*) '    (first error only)'
                    ddphi1_J_tx=0
                 end if
              end if

              phi1_J=    (
     &               g0_uu(1,1)*ddphi1_J
     &                         *(1-x0**2)**3
     &               -4*(3)*x0*g0_uu(1,2)*dphi1_J
     &                         *(1-x0**2)**(3)
     &               +2*g0_uu(1,2)*ddphi1_J_tx
     &                         *(1-x0**2)**(3)
     &                   )
     &              +
     &                dphi1_J*(1-x0**2)**3
     &                *(
     &                  g0_uu_x(1,1,1)+
     &                  g0_uu_x(2,1,2)+
     &                  g0_uu_x(3,1,3)+
     &                  g0_uu_x(4,1,4)+
     &                  g0_uu_x(5,1,5)
     &                 )
     &              +
     &                dphi1_J*(1-x0**2)**3
     &                *(
     &                  g0_uu(1,1)*gamma_ull(1,1,1)+
     &                  g0_uu(1,1)*gamma_ull(2,2,1)+
     &                  g0_uu(1,1)*gamma_ull(3,3,1)+
     &                  g0_uu(1,1)*gamma_ull(4,4,1)+
     &                  g0_uu(1,1)*gamma_ull(5,5,1)+
     &                  g0_uu(1,2)*gamma_ull(1,1,2)+
     &                  g0_uu(1,2)*gamma_ull(2,2,2)+
     &                  g0_uu(1,2)*gamma_ull(3,3,2)+
     &                  g0_uu(1,2)*gamma_ull(4,4,2)+
     &                  g0_uu(1,2)*gamma_ull(5,5,2)+
     &                  g0_uu(1,3)*gamma_ull(1,1,3)+
     &                  g0_uu(1,3)*gamma_ull(2,2,3)+
     &                  g0_uu(1,3)*gamma_ull(3,3,3)+
     &                  g0_uu(1,3)*gamma_ull(4,4,3)+
     &                  g0_uu(1,3)*gamma_ull(5,5,3)+
     &                  g0_uu(1,4)*gamma_ull(1,1,4)+
     &                  g0_uu(1,4)*gamma_ull(2,2,4)+
     &                  g0_uu(1,4)*gamma_ull(3,3,4)+
     &                  g0_uu(1,4)*gamma_ull(4,4,4)+
     &                  g0_uu(1,4)*gamma_ull(5,5,4)+
     &                  g0_uu(1,5)*gamma_ull(1,1,5)+
     &                  g0_uu(1,5)*gamma_ull(2,2,5)+
     &                  g0_uu(1,5)*gamma_ull(3,3,5)+
     &                  g0_uu(1,5)*gamma_ull(4,4,5)+
     &                  g0_uu(1,5)*gamma_ull(5,5,5)
     &                 )

              ! constraint damping terms added to efe,efe_J
              do a=1,5
                do b=1,5
                  cd_ll(a,b)=-kappa_cd*
     &                ( n_l(a)*c_l(a)+n_l(b)*c_l(a)-(1+rho_cd)*
     &                  g0_ll(a,b)*ndotc )
                end do
              end do

              cd_J_ll(1,1)=-kappa_cd*
     &            ( n_l(1)*g0_uu(1,1)*dgb_J*(1-x0**2)-(1+rho_cd)*
     &              g0_ll(1,1)*n_u(1)*0.5d0*g0_uu(1,1)*dgb_J*(1-x0**2) )
              cd_J_ll(1,2)=-kappa_cd*
     &            ( n_l(1)*0.5d0*g0_uu(1,1)*dgb_J*(1-x0**2)-(1+rho_cd)*
     &              g0_ll(1,2)*n_u(2)*g0_uu(1,1)*dgb_J*(1-x0**2)**2 )
              cd_J_ll(2,2)=-kappa_cd*
     &            ( 2*n_l(2)*g0_uu(1,2)*dgb_J*(1-x0**2)-(1+rho_cd)*
     &              g0_ll(2,2)*
     &              (-n_u(1)*0.5d0*g0_uu(2,2)*dgb_J*(1-x0**2)
     &               +n_u(2)*g0_uu(1,2)*dgb_J*(1-x0**2)) ) 
              cd_J_ll(3,3)=-kappa_cd*(1+rho_cd)*
     &              ( n_u(1)*0.5d0*g0_uu(3,3)*dgb_J*(1-x0**2)**2*x0**2 )
 
              if (kappa_cd.ne.0) then
                efe(1,1)=efe(1,1)+cd_ll(1,1)
                efe(1,2)=efe(1,2)+cd_ll(1,2)
                efe(2,2)=efe(2,2)+cd_ll(2,2)
                efe(3,3)=efe(3,3)+cd_ll(3,3)
                efe_J(1,1)=efe_J(1,1)+cd_J_ll(1,1)
                efe_J(1,2)=efe_J(1,2)+cd_J_ll(1,2)
                efe_J(2,2)=efe_J(2,2)+cd_J_ll(2,2)
                efe_J(3,3)=efe_J(3,3)+cd_J_ll(3,3)
              end if

              ! update gbars 
              if (is_nan(efe(1,1)).or.is_nan(efe_J(1,1)).or.
     &          efe_J(1,1).eq.0) then
                dump=.true.
              else
                gb_tt_np1(i)=gb_tt_np1(i)-efe(1,1)/efe_J(1,1)
              end if

              if (is_nan(efe(1,2)).or.is_nan(efe_J(1,2)).or.
     &          efe_J(1,2).eq.0) then
                dump=.true.
              else
                gb_tx_np1(i)=gb_tx_np1(i)-efe(1,2)/efe_J(1,2)
              end if

              if (is_nan(efe(2,2)).or.is_nan(efe_J(2,2)).or.
     &          efe_J(2,2).eq.0) then
                dump=.true.
              else
                gb_xx_np1(i)=gb_xx_np1(i)-efe(2,2)/efe_J(2,2)
              end if

              if (is_nan(efe(3,3)).or.is_nan(efe_J(3,3)).or.
     &          efe_J(3,3).eq.0) then
                dump=.true.
              else
                psi_np1(i)=psi_np1(i)-efe(3,3)/efe_J(3,3)
              end if

              ! update phi1, phi4 
              if (is_nan(phi1_res).or.is_nan(phi1_J)) then
                dump=.true.
              else
                phi1_np1(i)=phi1_np1(i)-phi1_res/phi1_J 
              end if

              ! save residuals
!              gb_res(i) =
!     &          max(abs(efe(1,1)/efe_J(1,1)),
!     &              abs(efe(1,2)/efe_J(1,2)),
!     &              abs(efe(2,2)/efe_J(2,2)),
!     &              abs(efe(3,3)/efe_J(3,3)))
              gb_res(i) =
     &          max(abs(efe(1,1)),
     &              abs(efe(1,2)),
     &              abs(efe(2,2)),
     &              abs(efe(3,3)))
              kg_res(i)=abs(phi1_res/phi1_J)

              ! check for NaNs
              if (dump.and.first_nan.or.ltrace) then
                first_nan=.false.
                write(*,*)
                write(*,*) 'g_evo_opt: Nan/zero at i,Nx=',
     &                                            i,Nx
                write(*,*) 'x=',x(i)
                write(*,*) 'dt,dx=',dt,dx
                write(*,*) 'x0',x0

                write(*,*) ' at tn:'
                write(*,*) ' gb_tt np1,n,nm1:',gb_tt_np1(i),
     &                 gb_tt_n(i),gb_tt_nm1(i)
                write(*,*) ' gb_tx np1,n,nm1:',gb_tx_np1(i),
     &                 gb_tx_n(i),gb_tx_nm1(i)
                write(*,*) ' gb_xx np1,n,nm1:',gb_xx_np1(i),
     &                 gb_xx_n(i),gb_xx_nm1(i)
                write(*,*) ' psi np1,n,nm1:',psi_np1(i),
     &                 psi_n(i),psi_nm1(i)
                write(*,*) ' g0_11 :',g0_ll(1,1)
                write(*,*) ' g0_12 :',g0_ll(1,2)
                write(*,*) ' g0_22 :',g0_ll(2,2)
                write(*,*) ' g0_33:',g0_ll(3,3)
                write(*,*) ' g0_44:',g0_ll(4,4)
                write(*,*) ' g0_55:',g0_ll(5,5)
                write(*,*) ' g0u_11 :',g0_uu(1,1)
                write(*,*) ' g0u_12 :',g0_uu(1,2)
                write(*,*) ' g0u_22 :',g0_uu(2,2)
                write(*,*) ' g0u_33:',g0_uu(3,3)
                write(*,*) ' g0u_44:',g0_uu(4,4)
                write(*,*) ' g0u_55:',g0_uu(5,5)
                write(*,*) ' phi np1,n,nm1:',phi1_np1(i),
     &                   phi1_n(i),phi1_nm1(i)
                write(*,*) ' phi1_t/x:',phi1_t,phi1_x
                write(*,*) ' phi1_ti..:',phi1_tt,phi1_tx,phi1_xx

                write(*,*) ' res J:'
                write(*,*) ' tt:',efe(1,1),efe_J(1,1)
                write(*,*) ' tx:',efe(1,2),efe_J(1,2)
                write(*,*) ' xx:',efe(2,2),efe_J(2,2)
                write(*,*) ' psi:',efe(3,3),efe_J(3,3)
                write(*,*) ' phi1:',phi1_res,phi1_J
              end if

            ! (REGION) non-interior points; set to zero prior to applying bcs 
            else 
              gb_tt_np1(i) = 0
              gb_tx_np1(i) = 0
              gb_xx_np1(i) = 0
              psi_np1(i)   = 0
              phi1_np1(i)  = 0 

            endif ! (near start of main loop)

          end do

        end do

        ! (REGION) x=0; impose regularity conditions 
        call axi_reg_Hb(Hb_t_np1,Hb_x_np1,chr,ex,L,x,Nx)
        call axi_reg_g(gb_tt_np1,gb_tx_np1,gb_xx_np1,
     &                            psi_np1,chr,ex,L,x,Nx)
        call axi_reg_phi(phi1_np1,chr,ex,L,x,Nx)

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
