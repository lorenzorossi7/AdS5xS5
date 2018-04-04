c----------------------------------------------------------------------
c NEED HARDCODED VERSION (that uses g_evo_opt.f instead of f_tt0_init.inc)
c
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
     &                      gb_xx_np1,gb_xx_n,gb_xx_nm1,gb_xx_t_n,
     &                      gb_yy_np1,gb_yy_n,gb_yy_nm1,gb_yy_t_n,
     &                      psi_np1,psi_n,psi_nm1,psi_t_n,
     &                      omega_np1,omega_n,omega_nm1,omega_t_n,
     &                      Hb_t_np1,Hb_t_n,Hb_t_nm1,Hb_t_t_n,
     &                      Hb_x_np1,Hb_x_n,Hb_x_nm1,Hb_x_t_n,
     &                      phi1_np1,phi1_n,phi1_nm1,phi1_t_n,
     &                      L,phys_bdy,x,dt,chr,ex,Nx)
        implicit none
        integer Nx
        integer phys_bdy(2)
        real*8 dt,ex,L
        real*8 chr(Nx)
        real*8 gb_tt_np1(Nx),gb_tt_n(Nx),gb_tt_nm1(Nx),gb_tt_t_n(Nx)
        real*8 gb_tx_np1(Nx),gb_tx_n(Nx),gb_tx_nm1(Nx),gb_tx_t_n(Nx)
        real*8 gb_xx_np1(Nx),gb_xx_n(Nx),gb_xx_nm1(Nx),gb_xx_t_n(Nx)
        real*8 gb_yy_np1(Nx),gb_yy_n(Nx),gb_yy_nm1(Nx),gb_yy_t_n(Nx)
        real*8 psi_np1(Nx),psi_n(Nx),psi_nm1(Nx),psi_t_n(Nx)
        real*8 omega_np1(Nx),omega_n(Nx),omega_nm1(Nx),omega_t_n(Nx)
        real*8 Hb_t_np1(Nx),Hb_t_n(Nx),Hb_t_nm1(Nx),Hb_t_t_n(Nx)
        real*8 Hb_x_np1(Nx),Hb_x_n(Nx),Hb_x_nm1(Nx),Hb_x_t_n(Nx)
        real*8 phi1_np1(Nx),phi1_n(Nx),phi1_nm1(Nx),phi1_t_n(Nx)

        real*8 x(Nx)

        real*8 lambda5
        real*8 tr_set

        logical is_nan

        !--------------------------------------------------------------
        ! the following are first and second time derivatives of *n*
        ! level variables, and as these are the only derivatives we
        ! use we drop any _n identifier
        !--------------------------------------------------------------
        real*8 gb_tt_t,gb_tt_tt
        real*8 gb_tx_t,gb_tx_tt
        real*8 gb_xx_t,gb_xx_tt
        real*8 gb_yy_t,gb_yy_tt
        real*8 Hb_t_t
        real*8 Hb_x_t
        real*8 phi1_t,phi1_tt

        real*8 h0_ll_tt(3,3),phi10_tt

        real*8 x0

        integer i,is,ie
        integer a,b,c,d

        real*8 dV1_dphi10 ! NEED TO HAVE THIS AS AN INPUT ARGUMENT

        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 dx

        logical ltrace
        parameter (ltrace=.false.)

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
        data i,is,ie/0,0,0/

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

        !---------------------------------------------------------------

        dx=x(2)-x(1)

        ! AdS5D cosmological constant
        !(lambda5=-(n-1)(n-2)/L^2) for n=5 dimensional AdS)
!        lambda5=-6/L/L
        ! NOTE: TEMPORARY CHECK
        lambda5=0.0d0

        do i=1,Nx
          phi1_nm1(i)=phi1_n(i)
        end do

        is=1
        ie=Nx
        if (phys_bdy(1).eq.1) is=2
        if (phys_bdy(2).eq.1) ie=Nx-1

        do i=is,ie
          x0=x(i)

          if (chr(i).ne.ex) then

          !-----------------------------------------------------------
          ! some other initializion, which needs to be done before
          ! temporal derivatives are calculated
          !-----------------------------------------------------------

          ! computes tensors at point i 
          call tensor_init(
     &            gb_tt_n,gb_tt_n,gb_tt_n,
     &            gb_tx_n,gb_tx_n,gb_tx_n,
     &            gb_xx_n,gb_xx_n,gb_xx_n,
     &            gb_yy_n,gb_yy_n,gb_yy_n,
     &            psi_n,psi_n,psi_n,
     &            omega_n,omega_n,omega_n,
     &            Hb_t_n,Hb_t_n,Hb_t_n,
     &            Hb_x_n,Hb_x_n,Hb_x_n,
     &            phi1_n,phi1_n,phi1_n,
     &            g0_ll,g0_uu,g0_ll_x,g0_uu_x,g0_ll_xx,
     &            gads_ll,gads_uu,gads_ll_x,gads_uu_x,gads_ll_xx,
     &            gA,gB,gA_x,gB_x,gA_xx,gB_xx,
     &            gAads,gBads,gAads_x,gBads_x,gAads_xx,gBads_xx,
     &            h0_ll,h0_uu,h0_ll_x,h0_uu_x,h0_ll_xx,
     &            A_l,A_l_x,Hads_l,Hads_l_x,
     &            gamma_ull,gamma_ull_x,
     &            riemann_ulll,ricci_ll,ricci_lu,ricci,
     &            einstein_ll,set_ll,f1_l,f2_ll,
     &            phi10_x,phi10_xx,
     &            x,dt,chr,L,ex,Nx,i)

          tr_set =set_ll(1,1)*g0_uu(1,1)+
     &            set_ll(2,2)*g0_uu(2,2)+
     &            set_ll(3,3)*g0_uu(3,3)+
     &         2*(set_ll(1,2)*g0_uu(1,2)+
     &            set_ll(1,3)*g0_uu(1,3)+
     &            set_ll(2,3)*g0_uu(2,3))

          ! initial first time derivatives; gb_ii_t_n,Hb_i_t_n,phi1_t_n were set in AdS5xS5_free_data()

          ! need this in h0_ll_tt,phi10_tt calculations
          phi10_x(1)    =phi1_t_n(i)*(1-x0**2)**3   
          h0_ll_x(1,1,1)=gb_tt_t_n(i)*(1-x0**2)  
          h0_ll_x(1,2,1)=gb_tx_t_n(i)*(1-x0**2)**2
          h0_ll_x(2,2,1)=gb_xx_t_n(i)*(1-x0**2)  
          h0_ll_x(3,3,1)=gb_yy_t_n(i)*(1-x0**2)  
          A_l_x(1,1)    =Hb_t_t_n(i)*(1-x0**2)**3
          A_l_x(2,1)    =Hb_x_t_n(i)*(1-x0**2)**2

          ! need this in gb_ii_nm1/np1,Hb_i_nm1/np1,phi1_nm1/np1 updates
          phi1_t =phi1_t_n(i)                
          gb_tt_t=gb_tt_t_n(i) 
          gb_tx_t=gb_tx_t_n(i) 
          gb_xx_t=gb_xx_t_n(i) 
          gb_yy_t=gb_yy_t_n(i) 
          Hb_t_t =Hb_t_t_n(i)
          Hb_x_t =Hb_x_t_n(i)

          ! 0 = efe_ab
          do a=1,3
            do b=a,3
              h0_ll_tt(a,b)=2/g0_uu(1,1)*
     &            ( 
     &                       -0.5d0*(                             
     &                          h0_uu(2,2)*h0_ll_xx(a,b,2,2)+
     &                          h0_uu(3,3)*h0_ll_xx(a,b,3,3)+
     &                       2*(h0_uu(1,2)*h0_ll_xx(a,b,1,2)+
     &                          h0_uu(1,3)*h0_ll_xx(a,b,1,3)+
     &                          h0_uu(2,3)*h0_ll_xx(a,b,2,3))
     &                       +
     &                          gads_uu(2,2)*h0_ll_xx(a,b,2,2)+
     &                          gads_uu(3,3)*h0_ll_xx(a,b,3,3)+
     &                       2*(gads_uu(1,2)*h0_ll_xx(a,b,1,2)+
     &                          gads_uu(1,3)*h0_ll_xx(a,b,1,3)+
     &                          gads_uu(2,3)*h0_ll_xx(a,b,2,3))
     &                       +
     &                          h0_uu(1,1)*gads_ll_xx(a,b,1,1)+
     &                          h0_uu(2,2)*gads_ll_xx(a,b,2,2)+
     &                          h0_uu(3,3)*gads_ll_xx(a,b,3,3)+
     &                       2*(h0_uu(1,2)*gads_ll_xx(a,b,1,2)+
     &                          h0_uu(1,3)*gads_ll_xx(a,b,1,3)+
     &                          h0_uu(2,3)*gads_ll_xx(a,b,2,3))
     &                       +
     &                          gads_uu(1,1)*gads_ll_xx(a,b,1,1)+
     &                          gads_uu(2,2)*gads_ll_xx(a,b,2,2)+
     &                          gads_uu(3,3)*gads_ll_xx(a,b,3,3)+
     &                       2*(gads_uu(1,2)*gads_ll_xx(a,b,1,2)+
     &                          gads_uu(1,3)*gads_ll_xx(a,b,1,3)+
     &                          gads_uu(2,3)*gads_ll_xx(a,b,2,3))
     &                              )
     &
     &                       -0.5d0*(
     &                          g0_uu_x(1,1,a)* g0_ll_x(b,1,1) +
     &                          g0_uu_x(1,2,a)*(g0_ll_x(b,1,2) +
     &                                          g0_ll_x(b,2,1))+
     &                          g0_uu_x(1,3,a)*(g0_ll_x(b,1,3) +
     &                                          g0_ll_x(b,3,1))+
     &                          g0_uu_x(2,2,a)* g0_ll_x(b,2,2) +
     &                          g0_uu_x(2,3,a)*(g0_ll_x(b,2,3) +
     &                                          g0_ll_x(b,3,2))+
     &                          g0_uu_x(3,3,a)* g0_ll_x(b,3,3) 
     &                            )
     &
     &                       -0.5d0*(
     &                          g0_uu_x(1,1,b)* g0_ll_x(a,1,1) +
     &                          g0_uu_x(1,2,b)*(g0_ll_x(a,1,2) +
     &                                          g0_ll_x(a,2,1))+
     &                          g0_uu_x(1,3,b)*(g0_ll_x(a,1,3) +
     &                                          g0_ll_x(a,3,1))+
     &                          g0_uu_x(2,2,b)* g0_ll_x(a,2,2) +
     &                          g0_uu_x(2,3,b)*(g0_ll_x(a,2,3) +
     &                                          g0_ll_x(a,3,2))+
     &                          g0_uu_x(3,3,b)* g0_ll_x(a,3,3) 
     &                            )
     &
     &                       -0.5d0*(Hads_l_x(a,b)+A_l_x(a,b))
     &
     &                       -0.5d0*(Hads_l_x(b,a)+A_l_x(b,a))
     &
     &                           +(
     &                          (Hads_l(1)+A_l(1))*gamma_ull(1,a,b)+
     &                          (Hads_l(2)+A_l(2))*gamma_ull(2,a,b)+
     &                          (Hads_l(3)+A_l(3))*gamma_ull(3,a,b)
     &                            )
     &
     &                           -(
     &                          gamma_ull(1,1,b)*gamma_ull(1,1,a)+
     &                          gamma_ull(1,2,b)*gamma_ull(2,1,a)+
     &                          gamma_ull(1,3,b)*gamma_ull(3,1,a)+
     &                          gamma_ull(2,1,b)*gamma_ull(1,2,a)+
     &                          gamma_ull(2,2,b)*gamma_ull(2,2,a)+
     &                          gamma_ull(2,3,b)*gamma_ull(3,2,a)+
     &                          gamma_ull(3,1,b)*gamma_ull(1,3,a)+
     &                          gamma_ull(3,2,b)*gamma_ull(2,3,a)+
     &                          gamma_ull(3,3,b)*gamma_ull(3,3,a)
     &                            )
     & 
     &                       -2*lambda5*g0_ll(a,b)/3
     &
     &                       -8*PI*(set_ll(a,b)-tr_set*g0_ll(a,b)/3) 
     &            )
            end do
          end do          

          ! 0 = g^ab phi1,ab - g^ab gamma^c_ab phi1,c 
          phi10_tt=-1/g0_uu(1,1)
     &              *(    
     &                    phi10_xx(2,2)*g0_uu(2,2)+
     &                    phi10_xx(3,3)*g0_uu(3,3)+
     &                 2*(phi10_xx(1,2)*g0_uu(1,2)+
     &                    phi10_xx(1,3)*g0_uu(1,3)+
     &                    phi10_xx(2,3)*g0_uu(2,3))
     &                -
     &                    phi10_x(1)*( gamma_ull(1,1,1)*g0_uu(1,1)+
     &                                 gamma_ull(1,2,2)*g0_uu(2,2)+
     &                                 gamma_ull(1,3,3)*g0_uu(3,3)+
     &                              2*(gamma_ull(1,1,2)*g0_uu(1,2)+
     &                                 gamma_ull(1,1,3)*g0_uu(1,3)+
     &                                 gamma_ull(1,2,3)*g0_uu(2,3)) )
     &                -
     &                    phi10_x(2)*( gamma_ull(2,1,1)*g0_uu(1,1)+
     &                                 gamma_ull(2,2,2)*g0_uu(2,2)+
     &                                 gamma_ull(2,3,3)*g0_uu(3,3)+
     &                              2*(gamma_ull(2,1,2)*g0_uu(1,2)+
     &                                 gamma_ull(2,1,3)*g0_uu(1,3)+
     &                                 gamma_ull(2,2,3)*g0_uu(2,3)) )
     &                -
     &                    phi10_x(3)*( gamma_ull(3,1,1)*g0_uu(1,1)+
     &                                 gamma_ull(3,2,2)*g0_uu(2,2)+
     &                                 gamma_ull(3,3,3)*g0_uu(3,3)+
     &                              2*(gamma_ull(3,1,2)*g0_uu(1,2)+
     &                                 gamma_ull(3,1,3)*g0_uu(1,3)+
     &                                 gamma_ull(3,2,3)*g0_uu(2,3)) )
     &                -
     &                    dV1_dphi10             
     &                                                                 )

          if (is_nan(h0_ll_tt(1,1)).or.is_nan(h0_ll_tt(1,2)) 
     &    .or.is_nan(h0_ll_tt(2,2)).or.is_nan(h0_ll_tt(3,3))) then
            write(*,*) 'h0_ll_tt(1,1)=',h0_ll_tt(1,1)
            write(*,*) 'h0_ll_tt(1,2)=',h0_ll_tt(1,2)
            write(*,*) 'h0_ll_tt(2,2)=',h0_ll_tt(2,2)
            write(*,*) 'h0_ll_tt(3,3)=',h0_ll_tt(3,3)
            stop
          end if

          ! initial second time derivatives
          gb_tt_tt=h0_ll_tt(1,1)/(1-x0**2) 
          gb_tx_tt=h0_ll_tt(1,2)/(1-x0**2)
          gb_xx_tt=h0_ll_tt(2,2)/(1-x0**2) 
          gb_yy_tt=h0_ll_tt(3,3)/(1-x0**2)**3 
          phi1_tt =phi10_tt/(1-x0**2)**3

          ! initialize past time level by O(h^3) expansion
           gb_tt_nm1(i)=gb_tt_n(i) - gb_tt_t*dt
     &                      + gb_tt_tt*dt**2/2
           gb_tx_nm1(i)=gb_tx_n(i) - gb_tx_t*dt
     &                      + gb_tx_tt*dt**2/2
           gb_xx_nm1(i)=gb_xx_n(i) - gb_xx_t*dt
     &                      + gb_xx_tt*dt**2/2
           gb_yy_nm1(i)=gb_yy_n(i) - gb_yy_t*dt
     &                      + gb_yy_tt*dt**2/2
!           psi_nm1(i)  =psi_n(i) - psi_t*dt  !NOTE: add these when you add psi subsector
!     &                      + psi_tt*dt**2/2
           Hb_t_nm1(i) =Hb_t_n(i) - Hb_t_t*dt
           Hb_x_nm1(i) =Hb_x_n(i) - Hb_x_t*dt
           phi1_nm1(i) =phi1_n(i) - phi1_t*dt
     &                      + phi1_tt*dt**2/2
        
          ! initialize future time level by O(h^3) expansion
           gb_tt_np1(i)=gb_tt_n(i) + gb_tt_t*dt
     &                      + gb_tt_tt*dt**2/2
           gb_tx_np1(i)=gb_tx_n(i) + gb_tx_t*dt
     &                      + gb_tx_tt*dt**2/2
           gb_xx_np1(i)=gb_xx_n(i) + gb_xx_t*dt
     &                      + gb_xx_tt*dt**2/2
           gb_yy_np1(i)=gb_yy_n(i) + gb_yy_t*dt
     &                      + gb_yy_tt*dt**2/2
!           psi_np1(i)  =psi_n(i) + psi_t*dt  !NOTE: !NOTE: add these when you add psi subsector
!     &                      + psi_tt*dt**2/2
           Hb_t_np1(i) =Hb_t_n(i) + Hb_t_t*dt
           Hb_x_np1(i) =Hb_x_n(i) + Hb_x_t*dt     
           phi1_np1(i) =phi1_n(i) + phi1_t*dt
     &                      + phi1_tt*dt**2/2
  
          end if

        end do

        return
        end
