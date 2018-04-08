c----------------------------------------------------------------------
c miscellaneous numerical routines for AdS5D1p1
c----------------------------------------------------------------------

c-----------------------------------------------------------------------
c specific x first derivative routine used by excision routines 
c below, for method 0,3,4 only
c
c-----------------------------------------------------------------------
        subroutine df1_int_x(f,f_x,dx,i,chr,ex,Nx)

        implicit none
        integer Nx,i
        real*8 f(Nx),chr(Nx),ex,f_x,dx
        logical first
        save first
        data first/.true./

        if (i.eq.1.or.chr(i-1).eq.ex) then
           if (i.le.(Nx-1).and.chr(i+1).ne.ex) then
              f_x=(-f(i)+f(i+1))/dx
!              write(*,*) 'df1_int_x: warning ... i=1 first order'
!              write(*,*) '    i,Nx,dx=',i,Nx,dx
           else
              if (first) then
                 first=.false.
                 write(*,*) 'df1_int_x: error in chr stencil (A)'
                 write(*,*) '    i,Nx,dx=',i,Nx,dx
                 write(*,*) '    (first error only)'
              end if
              f_x=0
              return
           end if
        else if (i.eq.Nx.or.chr(i+1).eq.ex) then
           if (i.ge.2.and.chr(i-1).ne.ex) then
              f_x=(f(i)-f(i-1))/dx
!              write(*,*) 'df1_int: warning ... i=Nx first order'
!              write(*,*) '    i,Nx,dx=',i,Nx,dx
           else
              if (first) then
                 first=.false.
                 write(*,*) 'df1_int: error in chr stencil (B)'
                 write(*,*) '    i,Nx,dx=',i,Nx,dx
                 write(*,*) '    (first error only)'
              end if
              f_x=0
              return
           end if
        else
           if (chr(i+1).ne.ex.and.chr(i-1).ne.ex) then
              f_x=(f(i+1)-f(i-1))/2/dx
           else
              if (first) then
                 first=.false.
                 write(*,*) 'df1_int: error in chr stencil (C)'
                 write(*,*) '    i,Nx,dx=',i,Nx,dx
                 write(*,*) '    (first error only)'
              end if
              f_x=0
              return
           end if
        end if

        return
        end

c----------------------------------------------------------------------
c the following computes all first derivatives of f,
c at a point i at time level n.
c
c stencil reduces to first order on excision surface
c----------------------------------------------------------------------
        subroutine df1_int(f_np1,f_n,f_nm1,f_t,f_x,
     &                     dx,dt,i,chr,ex,Nx,name)
        implicit none
        integer Nx,i
        real*8 f_np1(Nx),f_n(Nx),f_nm1(Nx)
        real*8 f_t,f_x,dx,dt,ex,chr(Nx)
        character*(*) name

        logical ltrace,first
        parameter (ltrace=.false.) 
        save first
        data first/.true./
        real*8 f_x_np1,f_x_nm1

        if (chr(i).eq.ex) then
         write(*,*) 'df1_int: error ... point excised'
         stop
        end if

        f_t=(f_np1(i)-f_nm1(i))/2/dt

        call df1_int_x(f_n,f_x,dx,i,chr,ex,Nx)

        if (ltrace) then
           write(*,*) 'df1_int for ',name
           write(*,*) ' f_t=',f_t
           write(*,*) ' f_x=',f_x
        end if

        return
        end

c----------------------------------------------------------------------
c same as df1_int above, but computes all second derivatives as 
c well 
c----------------------------------------------------------------------
        subroutine df2_int(f_np1,f_n,f_nm1,
     &                     f_t,f_x,
     &                     f_tt,f_tx,f_xx,
     &                     dx,dt,i,chr,ex,Nx,name)
        implicit none
        integer Nx,i
        real*8 f_np1(Nx),f_n(Nx),f_nm1(Nx)
        real*8 f_t,f_x,f_tt,f_tx,f_xx
        real*8 dx,dt,ex,chr(Nx)
        character*(*) name
        logical first
        save first
        data first/.true./

        logical ltrace
        parameter (ltrace=.false.)
        real*8 f_x_np1,f_x_nm1

        call df1_int(f_np1,f_n,f_nm1,f_t,f_x,
     &               dx,dt,i,chr,ex,Nx,name)

        f_tt=(f_np1(i)-2*f_n(i)+f_nm1(i))/dt/dt 

        f_xx=0

        if (chr(i).eq.ex) then
          write(*,*) 'df2_int: error ... point excised'
          stop
        end if

        call df1_int_x(f_np1,f_x_np1,dx,i,chr,ex,Nx)
        call df1_int_x(f_nm1,f_x_nm1,dx,i,chr,ex,Nx)
        f_tx=(f_x_np1-f_x_nm1)/2/dt

        !i

        if (i.eq.1.or.(chr(i-1).eq.ex)) then
          if (i.ge.(Nx-1).or.
     &        chr(i+1).eq.ex.or.chr(i+2).eq.ex) then
            if (first) then
              first=.false.
              write(*,*) 'df2_int: error in chr (A)'
              write(*,*) '    i,Nx,dx=',i,Nx,dx
              write(*,*) '    (first error only)'
            end if
            return
          end if
          f_xx=(f_n(i+2)-2*f_n(i+1)+f_n(i))/dx/dx 

        else if (i.eq.Nx.or.(chr(i+1).eq.ex)) then
          if (i.le.2.or.
     &        chr(i-1).eq.ex.or.chr(i-2).eq.ex) then
            if (first) then
              first=.false.
              write(*,*) 'df2_int: error in chr (B)'
              write(*,*) '    i,Nx,dx=',i,Nx,dx
              write(*,*) '    (first error only)'
            end if
            return
          end if
          f_xx=(f_n(i)-2*f_n(i-1)+f_n(i-2))/dx/dx 

        else if (chr(i+1).ne.ex.and.chr(i-1).ne.ex) then
          f_xx=(f_n(i+1)-2*f_n(i)+f_n(i-1))/dx/dx 

        else
          if (first) then
            first=.false.
            write(*,*) 'df2_int: error in chr (C)'
            write(*,*) '    i,Nx,dx=',i,Nx,dx
            write(*,*) '    (first error only)'
          end if
          return
        end if

        if (ltrace) then
           write(*,*) 'df2_int for ',name
           write(*,*) ' f_tt=',f_tt
           write(*,*) ' f_tx=',f_tx
           write(*,*) ' f_xx=',f_xx
        end if

        return
        end


c----------------------------------------------------------------------
c in polar coordinates t,x for x in [0,1]
c
c initializes f with a 2d gaussian-like profile ... multiplication
c by (1-x^2) is for correct asymptotics for AdS 
c
c f = amp*(1-x^2)*exp (- (r-r0)^2/delta^2) + A  , x > r0
c   = amp*(1-x^2) + A		                , x < r0
c
c where the profile is constructed to be phi1~(1-x^2)=(1-x)(1+x)
c for correct AdS asymptotics (1-x^2)^3*phi1~(1-x)^4, and regularity at x=0 
c----------------------------------------------------------------------
        subroutine gauss2d(f,amp,r0,delta,xu0,ax,
     &                     amp2,r02,delta2,xu02,ax2,
     &                     L,x,Nx)
        implicit none
        integer i,Nx
        real*8 f(Nx),x(Nx),L
        real*8 amp,r0,delta,ax,xu0,r
        real*8 amp2,r02,delta2,ax2,xu02,r2
        real*8 x0

        real*8 PI
        parameter (PI=3.141592653589793d0)

        !--------------------------------------------------------------

        do i=1,Nx
          f(i)=0
          x0=x(i)
          
          r=sqrt((x0-xu0)**2/ax**2)
          r2=sqrt((x0-xu02)**2/ax2**2)

          ! Gaussian phi=amp*exp(-r^2)/r^4 profile 
          ! remember that phi=phi1*(1-x^2)^3
          if (r.ge.1) then
             f(i)=0
          else if (r.gt.r0) then
             f(i)=amp*exp(-((r-r0)/delta)**2)*(1-x0**2)
          else
             f(i)=amp*(1-x0**2)
          end if

        end do

        return
        end

c-----------------------------------------------------------------------
c for variables with lin_zero_bnd ... zeros residual there
c-----------------------------------------------------------------------
        subroutine lin_zero_bnd_res(f,phys_bdy,all,Nx)
        implicit none
        integer Nx,all
        real*8 f(Nx)
        integer phys_bdy(2)

        integer i,is,ie

        ! initialize fixed-size variables
        data i,is,ie/0,0,0/

        !--------------------------------------------------------------

        if (phys_bdy(1).eq.1.or.all.eq.1) then
          f(2)=0
        end if

        if (phys_bdy(2).eq.1.or.all.eq.1) then
          f(Nx-1)=0
        end if

        return
        end

c----------------------------------------------------------------------
c initializes the metric to an exact black hole solution
c with radius parameter r0
c----------------------------------------------------------------------
        subroutine init_schw(gb_tt,gb_tx,gb_xx,gb_yy,psi,omega,
     &                       r0,L,phys_bdy,chr,ex,x,Nx)
        implicit none
        real*8 gb_tt(Nx),gb_tx(Nx),gb_xx(Nx)
        real*8 gb_yy(Nx),psi(Nx),omega(Nx)
        real*8 chr(Nx),x(Nx)
        real*8 r0,ex,L
        integer phys_bdy(2)
        integer i,Nx

        real*8 r_h,x_h,x0

        !--------------------------------------------------------------

        ! compute horizon global radius r_h and corresponding compactified x_h
        r_h=L*sqrt(2*sqrt(1+4*(r0/L)**2)-2)/2
        x_h=r_h/(1+r_h)

        ! initialize metric 
        do i=1,Nx
           if (chr(i).eq.ex) then
              gb_tt(i)=0
              gb_tx(i)=0
              gb_xx(i)=0
              psi(i)=0
           else
              x0=x(i)

              ! EF-like-near-horizon Schwarzschild-like-near-bdy coordinates
              gb_tt(i)=(r0/x0)**2/(1+x0)*(1-x0)
              gb_tx(i)=1/(1+x0)/(1-x_h)**4*(1-x0)
              gb_xx(i)=-((L**2*(-1 + x0)*((-1 + x0)**6*x0**4 + 
     &        L**2*((-1 + x0)**8*x0**2 - r0**2*(-1 + x_h)**10)))/
     &        ((1 + x0)*(L**2*(-1 + x0)**2 + x0**2)*
     &        (-x0**4 + L**2*(-1 + x0)**2*(r0**2*(-1 + x0)**2 - x0**2))*
     &        (-1 + x_h)**10)) 
              gb_yy(i)=0
              psi(i)  =0
              omega(i)=0

!              ! Schwarzschild coordinates (WARNING: not horizon-penetrating)
!              gb_tt(i)=(r0/x0)**2/(1+x0)*(1-x0)
!              gb_tx(i)=0
!              gb_xx(i)=(L**4*r0**2*(1 - x0))/
!     &                ((1 + (1 - x0)*(-1 + L**2*(1 - x0) - x0))*
!     &                (-1 - x0)*(L**2*r0**2*(1 - x0)**4 - 
!     &                (1 + (1 - x0)*(-1 + L**2*(1 - x0) - x0))*x0**2)
!     &                )
!              psi(i)=0

           end if
        end do

        ! (REGION) x=0; impose regularity conditions 
        call axi_reg_g(gb_tt,gb_tx,gb_xx,gb_yy,psi,chr,ex,L,x,Nx)

        return
        end

c----------------------------------------------------------------------
c calculates all the tensorial objects in x coordinates, at point i
c----------------------------------------------------------------------
        subroutine tensor_init(
     &                  gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                  gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                  gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                  gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &                  psi_np1,psi_n,psi_nm1,
     &                  omega_np1,omega_n,omega_nm1,
     &                  Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &                  Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &                  phi1_np1,phi1_n,phi1_nm1,
     &                  g0_ll,g0_uu,g0_ll_x,g0_uu_x,g0_ll_xx,
     &                  gads_ll,gads_uu,gads_ll_x,gads_uu_x,gads_ll_xx,
     &                  gA,gB,gA_x,gB_x,gA_xx,gB_xx,
     &                  gAads,gBads,gAads_x,gBads_x,gAads_xx,gBads_xx,
     &                  h0_ll,h0_uu,h0_ll_x,h0_uu_x,h0_ll_xx,
     &                  A_l,A_l_x,Hads_l,Hads_l_x,
     &                  gamma_ull,gamma_ull_x,
     &                  riemann_ulll,ricci_ll,ricci_lu,ricci,
     &                  s0_ll,t0_ll,f1_l,f2_ll,
     &                  phi10_x,phi10_xx,
     &                  x,dt,chr,L,ex,Nx,i)
        implicit none

        integer Nx
        integer i

        real*8 chr(Nx),ex
        real*8 x(Nx),dt,L

        real*8 gb_tt_np1(Nx),gb_tt_n(Nx),gb_tt_nm1(Nx)
        real*8 gb_tx_np1(Nx),gb_tx_n(Nx),gb_tx_nm1(Nx)
        real*8 gb_xx_np1(Nx),gb_xx_n(Nx),gb_xx_nm1(Nx)
        real*8 gb_yy_np1(Nx),gb_yy_n(Nx),gb_yy_nm1(Nx)
        real*8 psi_np1(Nx),psi_n(Nx),psi_nm1(Nx)
        real*8 omega_np1(Nx),omega_n(Nx),omega_nm1(Nx)
        real*8 Hb_t_np1(Nx),Hb_t_n(Nx),Hb_t_nm1(Nx)
        real*8 Hb_x_np1(Nx),Hb_x_n(Nx),Hb_x_nm1(Nx)
        real*8 phi1_np1(Nx),phi1_n(Nx),phi1_nm1(Nx)
        real*8 fb_t_np1(Nx),fb_t_n(Nx),fb_t_nm1(Nx)
        real*8 fb_x_np1(Nx),fb_x_n(Nx),fb_x_nm1(Nx)
        real*8 fb_y_np1(Nx),fb_y_n(Nx),fb_y_nm1(Nx)

        integer a,b,c,d,e,f,m,n,o,p

        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 dx,dy
        real*8 x0,y0   

        real*8 grad_phi1_sq

        integer dimA,dimB

        !--------------------------------------------------------------
        ! variables for tensor manipulations 
        !(indices are t,x,y,psi,omega)
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
        real*8 s0_ll(3,3),t0_ll(3,3)
        real*8 Hads_l(3),Hads_l_x(3,3),A_l(3),A_l_x(3,3)
        real*8 phi10_x(3),phi10_xx(3,3)
        real*8 f1_l(3),f2_ll(3,3)
        !NOTE: below I have implemented direct calculation
        real*8 levicivi3(3,3,3),vol(3,3,3),sqrtdetg
        real*8 riccibar_ll(3,3),riccibar_lu(3,3),riccibar

        real*8 efe(3,3),efest(3,3)
        real*8 term1(3,3),term2(3,3),term3(3,3),term4(3,3)
        real*8 term5(3,3),term6(3,3),term7(3,3),term8(3,3)
        real*8 term9(3,3)
        real*8 efe_norm,efest_norm

        !--------------------------------------------------------------
        ! the following are first and second time derivatives of *n*
        ! level variables, and as these are the only derivatives we
        ! use we drop any _n identifier
        !--------------------------------------------------------------
        real*8 gb_tt_t,gb_tt_x
        real*8 gb_tt_tt,gb_tt_tx,gb_tt_xx
        real*8 gb_tx_t,gb_tx_x
        real*8 gb_tx_tt,gb_tx_tx,gb_tx_xx
        real*8 gb_xx_t,gb_xx_x
        real*8 gb_xx_tt,gb_xx_tx,gb_xx_xx
        real*8 gb_yy_t,gb_yy_x
        real*8 gb_yy_tt,gb_yy_tx,gb_yy_xx
        real*8 psi_t,psi_x
        real*8 psi_tt,psi_tx,psi_xx
        real*8 omega_t,omega_x
        real*8 omega_tt,omega_tx,omega_xx
        real*8 phi1_t,phi1_x
        real*8 phi1_tt,phi1_tx,phi1_xx
  
        real*8 gb_tt0,gb_tx0,gb_xx0,gb_yy0,phi10

        real*8 g0_tt_ads0,g0_xx_ads0,g0_yy_ads0
        real*8 g0_tt_ads_x,g0_tt_ads_xx
        real*8 g0_xx_ads_x,g0_xx_ads_xx
        real*8 g0_yy_ads_x,g0_yy_ads_xx

        real*8 psi0,omega0

        real*8 gA_ads0,gB_ads0
        real*8 gA_ads_x,gA_ads_xx        
        real*8 gB_ads_yy

        real*8 Hb_t_t,Hb_t_x
        real*8 Hb_x_t,Hb_x_x

        real*8 Hb_t0,Hb_x0

        real*8 fb_t0,fb_x0,fb_y0
        real*8 fb_tx0,fb_ty0,fb_xy0

        real*8 f1_y_ads0
        real*8 f2_tx_ads0
!----------------------------------------------------------------------
        
        dx=(x(2)-x(1))

        x0=x(i)
        y0=L/2  !NOTE: change this to y0=y(j) when we add y-dependence

        ! set dimensions of S3 and S4 subspaces
        dimA=3
        dimB=4

        ! set gads values using sin(theta1)=sin(theta2)=1 w.l.o.g 
        !(considering theta1,theta2-independent case, so theta1=theta2=pi/2 slice will do)
        g0_tt_ads0 =-((1-x0)**2+x0**2)/(1-x0)**2
        g0_xx_ads0 =1/((1-x0)**2+x0**2)/(1-x0)**2
        g0_yy_ads0=PI**2

        ! set gAads values
        gA_ads0=x0**2/(1-x0)**2

        ! set gBads values
        gB_ads0=L**2*sin(PI*y0/L)**2

        ! set f1ads, f2ads values using sin(phi2)=sin(phi3)=sin(phi4)=1 w.l.o.g 
        !(considering phi2,phi3,phi4-independent case, so phi2=phi3=phi4=pi/2 slice will do)
        f1_y_ads0  = 4/L*PI*sin(PI*y0/L)**4
        f2_tx_ads0 = 4/L/(1-x0)**2

        ! set gbar values
        gb_tt0=gb_tt_n(i)
        gb_tx0=gb_tx_n(i)
        gb_xx0=gb_xx_n(i)
        gb_yy0=gb_yy_n(i)
        psi0  =psi_n(i) 
        omega0=omega_n(i)

        ! set fbar values
        fb_t0=fb_t_n(i)
        fb_x0=fb_x_n(i)
        fb_y0=fb_y_n(i)

        ! set hbar values
        Hb_t0=Hb_t_n(i)
        Hb_x0=Hb_x_n(i)

        ! set phi1 value
        phi10=phi1_n(i)

        ! set gads derivatives
        g0_tt_ads_x  =-2*x0/(1-x0)**3
        g0_tt_ads_xx =-(2+4*x0)/(1-x0)**4
        g0_xx_ads_x  =(4-2*(5-4*x0)*x0)
     &                /((1-x0)**2+x0**2)**2/(1-x0)**3
        g0_xx_ads_xx =(2*(3-10*(1-x0)*x0)*(3-6*x0+4*x0**2))
     &                /((1-x0)**2+x0**2)**3/(1-x0)**4
        g0_yy_ads_x =0
        g0_yy_ads_xx=0

        ! set gAads derivatives
        gA_ads_x =2*x0/(1-x0)**3
        gA_ads_xx=(2+4*x0)/(1-x0)**4
 
        ! set gBads derivatives
        gB_ads_yy=2*PI**2*cos(2*PI*y0/L)
 
        ! calculate gbar derivatives
        call df2_int(gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &       gb_tt_t,gb_tt_x,
     &       gb_tt_tt,gb_tt_tx,gb_tt_xx,
     &       dx,dt,i,chr,ex,Nx,'gb_tt')
        call df2_int(gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &       gb_tx_t,gb_tx_x,
     &       gb_tx_tt,gb_tx_tx,gb_tx_xx,
     &       dx,dt,i,chr,ex,Nx,'gb_tx')
        call df2_int(gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &       gb_xx_t,gb_xx_x,
     &       gb_xx_tt,gb_xx_tx,gb_xx_xx,
     &       dx,dt,i,chr,ex,Nx,'gb_xx')
        call df2_int(gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &       gb_yy_t,gb_yy_x,
     &       gb_yy_tt,gb_yy_tx,gb_yy_xx,
     &       dx,dt,i,chr,ex,Nx,'gb_yy')
        call df2_int(psi_np1,psi_n,psi_nm1,
     &       psi_t,psi_x,
     &       psi_tt,psi_tx,psi_xx,
     &       dx,dt,i,chr,ex,Nx,'psi')
        call df2_int(omega_np1,omega_n,omega_nm1,
     &       omega_t,omega_x,
     &       omega_tt,omega_tx,omega_xx,
     &       dx,dt,i,chr,ex,Nx,'omega')

        ! calculate hbar derivatives
        call df1_int(Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &       Hb_t_t,Hb_t_x,
     &       dx,dt,i,chr,ex,Nx,'Hb_t')
        call df1_int(Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &       Hb_x_t,Hb_x_x,
     &       dx,dt,i,chr,ex,Nx,'Hb_x')

        ! calculate phi1 derivatives
        call df2_int(phi1_np1,phi1_n,phi1_nm1,
     &       phi1_t,phi1_x,
     &       phi1_tt,phi1_tx,phi1_xx,
     &       dx,dt,i,chr,ex,Nx,'phi1')

        ! give values to the metric, using sin(theta1)=sin(theta2)=1 w.l.o.g 
        !(considering theta1,theta2-independent case, so theta1=theta2=pi/2 slice will do)
        g0_ll(1,1)=g0_tt_ads0+gb_tt0*(1-x0**2)
        g0_ll(1,2)=           gb_tx0*(1-x0**2)**2
        !g0_ll(1,3)=           gb_ty0*(1-x0**2)**2  !NOTE: add this when you add y-dependence
        g0_ll(2,2)=g0_xx_ads0+gb_xx0*(1-x0**2)
        g0_ll(3,3)=g0_yy_ads0+gb_yy0*(1-x0**2)*x0**2

        g0_uu(1,1)=
     &          (g0_ll(2,2)*g0_ll(3,3)-g0_ll(2,3)**2)
     &         /(-g0_ll(1,3)**2*g0_ll(2,2)
     &           +g0_ll(1,2)*g0_ll(1,3)*g0_ll(2,3)*2
     &           -g0_ll(1,1)*g0_ll(2,3)**2
     &           -g0_ll(1,2)**2*g0_ll(3,3)
     &           +g0_ll(1,1)*g0_ll(2,2)*g0_ll(3,3))
        g0_uu(1,2)=
     &         (g0_ll(1,3)*g0_ll(2,3)-g0_ll(1,2)*g0_ll(3,3))
     &         /(-g0_ll(1,3)**2*g0_ll(2,2)
     &           +g0_ll(1,2)*g0_ll(1,3)*g0_ll(2,3)*2
     &           -g0_ll(1,1)*g0_ll(2,3)**2
     &           -g0_ll(1,2)**2*g0_ll(3,3)
     &           +g0_ll(1,1)*g0_ll(2,2)*g0_ll(3,3))
        g0_uu(1,3)=
     &         (g0_ll(1,2)*g0_ll(2,3)-g0_ll(1,3)*g0_ll(2,2))
     &         /(-g0_ll(1,3)**2*g0_ll(2,2)
     &           +g0_ll(1,2)*g0_ll(1,3)*g0_ll(2,3)*2
     &           -g0_ll(1,1)*g0_ll(2,3)**2
     &           -g0_ll(1,2)**2*g0_ll(3,3)
     &           +g0_ll(1,1)*g0_ll(2,2)*g0_ll(3,3))
        g0_uu(2,2)=
     &         (g0_ll(1,1)*g0_ll(3,3)-g0_ll(1,3)**2)
     &         /(-g0_ll(1,3)**2*g0_ll(2,2)
     &           +g0_ll(1,2)*g0_ll(1,3)*g0_ll(2,3)*2
     &           -g0_ll(1,1)*g0_ll(2,3)**2
     &           -g0_ll(1,2)**2*g0_ll(3,3)
     &           +g0_ll(1,1)*g0_ll(2,2)*g0_ll(3,3))
        g0_uu(2,3)=
     &         (g0_ll(1,2)*g0_ll(1,3)-g0_ll(1,1)*g0_ll(2,3))
     &         /(-g0_ll(1,3)**2*g0_ll(2,2)
     &           +g0_ll(1,2)*g0_ll(1,3)*g0_ll(2,3)*2
     &           -g0_ll(1,1)*g0_ll(2,3)**2
     &           -g0_ll(1,2)**2*g0_ll(3,3)
     &           +g0_ll(1,1)*g0_ll(2,2)*g0_ll(3,3))
        g0_uu(3,3)=
     &         (g0_ll(1,1)*g0_ll(2,2)-g0_ll(1,2)**2)
     &         /(-g0_ll(1,3)**2*g0_ll(2,2)
     &           +g0_ll(1,2)*g0_ll(1,3)*g0_ll(2,3)*2
     &           -g0_ll(1,1)*g0_ll(2,3)**2
     &           -g0_ll(1,2)**2*g0_ll(3,3)
     &           +g0_ll(1,1)*g0_ll(2,2)*g0_ll(3,3))

        g0_ll_x(1,1,1)   =0
     &                   +gb_tt_t*(1-x0**2)
        g0_ll_x(1,1,2)   =g0_tt_ads_x
     &                   +gb_tt_x*(1-x0**2)
     &                   +gb_tt0*(-2*x0)
        g0_ll_xx(1,1,1,1)=0
     &                   +gb_tt_tt*(1-x0**2)
        g0_ll_xx(1,1,1,2)=0
     &                   +gb_tt_tx*(1-x0**2)
     &                   +gb_tt_t*(-2*x0)
        g0_ll_xx(1,1,2,2)=g0_tt_ads_xx
     &                   +gb_tt_xx*(1-x0**2)
     &                   +gb_tt_x*(-2*x0)
     &                   +gb_tt_x*(-2*x0)
     &                   +gb_tt0*(-2)

        g0_ll_x(1,2,1)   =0
     &                   +gb_tx_t*(1-x0**2)**2
        g0_ll_x(1,2,2)   =0
     &                   +gb_tx_x*(1-x0**2)**2
     &                   +gb_tx0*2*(1-x0**2)*(-2*x0)
        g0_ll_xx(1,2,1,1)=0
     &                   +gb_tx_tt*(1-x0**2)**2
        g0_ll_xx(1,2,1,2)=0
     &                   +gb_tx_tx*(1-x0**2)**2
     &                   +gb_tx_t*2*(1-x0**2)*(-2*x0)
        g0_ll_xx(1,2,2,2)=0
     &                   +gb_tx_xx*(1-x0**2)**2
     &                   +gb_tx_x*2*(1-x0**2)*(-2*x0)
     &                   +gb_tx_x*2*(1-x0**2)*(-2*x0)
     &                   +gb_tx0*2*(-2+6*x0**2)

!        g0_ll_x(1,3,1)   =0                            !NOTE: add this when you add y-dependence
!     &                   +gb_ty_t*(1-x0**2)**2
!        g0_ll_x(1,3,2)   =0
!     &                   +gb_ty_x*(1-x0**2)**2
!     &                   +gb_ty0*2*(1-x0**2)*(-2*x0)
!        g0_ll_xx(1,3,1,1)=0
!     &                   +gb_ty_tt*(1-x0**2)**2
!        g0_ll_xx(1,3,1,2)=0
!     &                   +gb_ty_tx*(1-x0**2)**2
!     &                   +gb_ty_t*2*(1-x0**2)*(-2*x0)
!        g0_ll_xx(1,3,2,2)=0
!     &                   +gb_ty_xx*(1-x0**2)**2
!     &                   +gb_ty_x*2*(1-x0**2)*(-2*x0)
!     &                   +gb_ty_x*2*(1-x0**2)*(-2*x0)
!     &                   +gb_ty0*2*(-2+6*x0**2)

        g0_ll_x(2,2,1)   =0
     &                   +gb_xx_t*(1-x0**2)
        g0_ll_x(2,2,2)   =g0_xx_ads_x
     &                   +gb_xx_x*(1-x0**2)
     &                   +gb_xx0*(-2*x0)
        g0_ll_xx(2,2,1,1)=0
     &                   +gb_xx_tt*(1-x0**2)
        g0_ll_xx(2,2,1,2)=0
     &                   +gb_xx_tx*(1-x0**2)
     &                   +gb_xx_t*(-2*x0)
        g0_ll_xx(2,2,2,2)=g0_xx_ads_xx
     &                   +gb_xx_xx*(1-x0**2)
     &                   +gb_xx_x*(-2*x0)
     &                   +gb_xx_x*(-2*x0)
     &                   +gb_xx0*(-2)

        g0_ll_x(3,3,1)   =0
     &                   +gb_yy_t*(1-x0**2)**3
        g0_ll_x(3,3,2)   =g0_yy_ads_x
     &                   +gb_yy_x*(1-x0**2)**3
     &                   +gb_yy0*3*(1-x0**2)**2*(-2*x0)
        g0_ll_xx(3,3,1,1)=0
     &                   +gb_yy_tt*(1-x0**2)**3
        g0_ll_xx(3,3,1,2)=0
     &                   +gb_yy_tx*(1-x0**2)**3
     &                   +gb_yy_t*3*(1-x0**2)**2*(-2*x0)
        g0_ll_xx(3,3,2,2)=g0_yy_ads_xx
     &                   +gb_yy_xx*(1-x0**2)**3
     &                   +gb_yy_x*3*(1-x0**2)**2*(-2*x0)
     &                   +gb_yy_x*3*(1-x0**2)**2*(-2*x0)
     &                   +gb_yy0*6*(-1+6*x0**2-5*x0**4)
        ! WARNING: from sin^2chi factor in pure ads term

        do a=1,2
          do b=a+1,3
            g0_ll(b,a)=g0_ll(a,b)
            g0_uu(b,a)=g0_uu(a,b) 
            do c=1,3
              g0_ll_x(b,a,c)=g0_ll_x(a,b,c)
            end do
          end do
        end do

        do a=1,3
          do b=1,3
            do c=1,3
              g0_uu_x(a,b,c)=0
              do d=1,3
                do e=1,3
                  g0_uu_x(a,b,c)=g0_uu_x(a,b,c)
     &                          -g0_ll_x(d,e,c)
     &                           *g0_uu(a,d)*g0_uu(b,e)
                end do
     &  
              end do
            end do
          end do
        end do

        do a=1,3
          do b=1,3
            do c=1,3
              do d=1,3
                g0_ll_xx(a,b,c,d)=
     &             g0_ll_xx(min(a,b),max(a,b),min(c,d),max(c,d))
              end do
            end do
          end do
        end do

        do a=1,3
          do b=1,3
            do c=1,3
              gamma_ull(a,b,c)=0
              do d=1,3
                gamma_ull(a,b,c)=gamma_ull(a,b,c)
     &                          +0.5d0*g0_uu(a,d)
     &                                *(g0_ll_x(c,d,b)
     &                                 -g0_ll_x(b,c,d)
     &                                 +g0_ll_x(d,b,c))
              end do
            end do
          end do
        end do

        gA=gA_ads0+psi0*(1-x0**2)
        gB=gB_ads0+omega0*(1-x0**2)**3

        gA_x(1)   =0
     &            +psi_t*(1-x0**2)
        gA_x(2)   =gA_ads_x
     &            +psi_x*(1-x0**2)
     &            +psi0*(-2*x0)
        gA_xx(1,1)=0
     &            +psi_tt*(1-x0**2)
        gA_xx(1,2)=0
     &            +psi_tx*(1-x0**2)
     &            +psi_t*(-2*x0)
        gA_xx(2,2)=gA_ads_xx
     &            +psi_xx*(1-x0**2)
     &            +psi_x*(-2*x0)
     &            +psi_x*(-2*x0)
     &            +psi0*(-2)

        gB_x(1)   =0
     &            +omega_t*(1-x0**2)**3
        gB_x(2)   =0
     &            +omega_x*(1-x0**2)**3
     &            +omega0*3*(1-x0**2)**2*(-2*x0)
        gB_xx(1,1)=0
     &            +omega_tt*(1-x0**2)**3
        gB_xx(1,2)=0
     &            +omega_tx*(1-x0**2)**3
     &            +omega_t*3*(1-x0**2)**2*(-2*x0)
        gB_xx(2,2)=0
     &            +omega_xx*(1-x0**2)**3
     &            +omega_x*3*(1-x0**2)**2*(-2*x0)
     &            +omega_x*3*(1-x0**2)**2*(-2*x0)
     &            +omega0*6*(-1+6*x0**2-5*x0**4)
        gB_xx(3,3)=gB_ads_yy  !NOTE: add the rest of the terms when you add y-dependence

        do a=1,2
          do b=a+1,3
            gA_xx(b,a)=gA_xx(a,b)
            gB_xx(b,a)=gB_xx(a,b)
          end do
        end do

        ! give values to the ads metric, using sin(theta1)=sin(theta2)=1 w.l.o.g 
        !(considering theta1,theta2-independent case, so theta1=theta2=pi/2 slice will do)
        gads_ll(1,1)=g0_tt_ads0
        gads_ll(2,2)=g0_xx_ads0
        gads_ll(3,3)=g0_yy_ads0

        gads_uu(1,1)=1/g0_tt_ads0
        gads_uu(2,2)=1/g0_xx_ads0
        gads_uu(3,3)=1/g0_yy_ads0

        gads_ll_x(1,1,2)   =g0_tt_ads_x  
        gads_ll_xx(1,1,2,2)=g0_tt_ads_xx 
        gads_ll_x(2,2,2)   =g0_xx_ads_x  
        gads_ll_xx(2,2,2,2)=g0_xx_ads_xx 
        gads_ll_x(3,3,2)   =g0_yy_ads_x 
        gads_ll_xx(3,3,2,2)=g0_yy_ads_xx
        ! WARNING: from sin^2theta factor in pure ads term

        do a=1,2
          do b=a+1,3
            gads_ll(b,a)=gads_ll(a,b)
            gads_uu(b,a)=gads_uu(a,b)
            do c=1,3
              gads_ll_x(b,a,c)=gads_ll_x(a,b,c)
            end do
          end do
        end do

        do a=1,3
          do b=1,3
            do c=1,3
              gads_uu_x(a,b,c)=
     &              -gads_ll_x(1,1,c)*gads_uu(a,1)*gads_uu(b,1)
     &              -gads_ll_x(1,2,c)*(gads_uu(a,1)*gads_uu(b,2)
     &                               +gads_uu(a,2)*gads_uu(b,1))
     &              -gads_ll_x(1,3,c)*(gads_uu(a,1)*gads_uu(b,3)
     &                               +gads_uu(a,3)*gads_uu(b,1))
     &              -gads_ll_x(2,2,c)*gads_uu(a,2)*gads_uu(b,2)
     &              -gads_ll_x(2,3,c)*(gads_uu(a,2)*gads_uu(b,3)
     &                               +gads_uu(a,3)*gads_uu(b,2))
     &              -gads_ll_x(3,3,c)*gads_uu(a,3)*gads_uu(b,3)
            end do
          end do
        end do

        gAads=gA_ads0

        gAads_x(2)   =gA_ads_x 
        gAads_xx(2,2)=gA_ads_xx

        gBads=gB_ads0

        gBads_xx(3,3)=gB_ads_yy

        ! give values to the metric deviation, using sin(theta1)=sin(theta2)=1 w.l.o.g 
        !(considering theta1,theta2-independent case, so theta1=theta2=pi/2 will do)
        h0_ll(1,1)=gb_tt0*(1-x0**2)
        h0_ll(1,2)=gb_tx0*(1-x0**2)**2
        !h0_ll(1,3)=gb_ty0*(1-x0**2)**2  !NOTE: add this when you add y-dependence
        h0_ll(2,2)=gb_xx0*(1-x0**2)
        h0_ll(3,3)=gb_yy0*(1-x0**2)**3
        
        h0_uu(1,1)=g0_uu(1,1)-gads_uu(1,1)
        h0_uu(1,2)=g0_uu(1,2)
        h0_uu(2,2)=g0_uu(2,2)-gads_uu(2,2)
        h0_uu(3,3)=g0_uu(3,3)-gads_uu(3,3)

        h0_ll_x(1,1,1)   =g0_ll_x(1,1,1)-gads_ll_x(1,1,1)
        h0_ll_x(1,1,2)   =g0_ll_x(1,1,2)-gads_ll_x(1,1,2)
        h0_ll_xx(1,1,1,1)=g0_ll_xx(1,1,1,1)-gads_ll_xx(1,1,1,1)
        h0_ll_xx(1,1,1,2)=g0_ll_xx(1,1,1,2)-gads_ll_xx(1,1,1,2)
        h0_ll_xx(1,1,2,2)=g0_ll_xx(1,1,2,2)-gads_ll_xx(1,1,2,2)

        h0_ll_x(1,2,1)   =g0_ll_x(1,2,1)-gads_ll_x(1,2,1)
        h0_ll_x(1,2,2)   =g0_ll_x(1,2,2)-gads_ll_x(1,2,2)
        h0_ll_xx(1,2,1,1)=g0_ll_xx(1,2,1,1)-gads_ll_xx(1,2,1,1)
        h0_ll_xx(1,2,1,2)=g0_ll_xx(1,2,1,2)-gads_ll_xx(1,2,1,2)
        h0_ll_xx(1,2,2,2)=g0_ll_xx(1,2,2,2)-gads_ll_xx(1,2,2,2)

        h0_ll_x(2,2,1)   =g0_ll_x(2,2,1)-gads_ll_x(2,2,1)
        h0_ll_x(2,2,2)   =g0_ll_x(2,2,2)-gads_ll_x(2,2,2)
        h0_ll_xx(2,2,1,1)=g0_ll_xx(2,2,1,1)-gads_ll_xx(2,2,1,1)
        h0_ll_xx(2,2,1,2)=g0_ll_xx(2,2,1,2)-gads_ll_xx(2,2,1,2)
        h0_ll_xx(2,2,2,2)=g0_ll_xx(2,2,2,2)-gads_ll_xx(2,2,2,2)

        h0_ll_x(3,3,1)   =g0_ll_x(3,3,1)-gads_ll_x(3,3,1)
        h0_ll_x(3,3,2)   =g0_ll_x(3,3,2)-gads_ll_x(3,3,2)
        h0_ll_xx(3,3,1,1)=g0_ll_xx(3,3,1,1)-gads_ll_xx(3,3,1,1)
        h0_ll_xx(3,3,1,2)=g0_ll_xx(3,3,1,2)-gads_ll_xx(3,3,1,2)
        h0_ll_xx(3,3,2,2)=g0_ll_xx(3,3,2,2)-gads_ll_xx(3,3,2,2)
        ! WARNING: from sin^2theta factor in pure ads term

        h0_uu_x(1,1,1)=g0_uu_x(1,1,1)-gads_uu_x(1,1,1)
        h0_uu_x(1,1,2)=g0_uu_x(1,1,2)-gads_uu_x(1,1,2)

        h0_uu_x(1,2,1)=g0_uu_x(1,2,1)-gads_uu_x(1,2,1)
        h0_uu_x(1,2,2)=g0_uu_x(1,2,2)-gads_uu_x(1,2,2)

        h0_uu_x(2,2,1)=g0_uu_x(2,2,1)-gads_uu_x(2,2,1)
        h0_uu_x(2,2,2)=g0_uu_x(2,2,2)-gads_uu_x(2,2,2)

        h0_uu_x(3,3,1)=g0_uu_x(3,3,1)-gads_uu_x(3,3,1)
        h0_uu_x(3,3,2)=g0_uu_x(3,3,2)-gads_uu_x(3,3,2)

        do a=1,2
          do b=a+1,3
            h0_ll(b,a)=h0_ll(a,b)
            h0_uu(b,a)=h0_uu(a,b)
            do c=1,3
              h0_ll_x(b,a,c)=h0_ll_x(a,b,c)
              h0_uu_x(b,a,c)=h0_uu_x(a,b,c)
            end do
          end do
        end do

        ! give values to the gh source functions
        A_l(1)=Hb_t0*(1-x0**2)**3
        A_l(2)=Hb_x0*(1-x0**2)**2

        A_l_x(1,1)=Hb_t_t*(1-x0**2)**3
        A_l_x(1,2)=Hb_t_x*(1-x0**2)**3
     &            -2*x0*3*(1-x0**2)**2*Hb_t0

        A_l_x(2,1)=Hb_x_t*(1-x0**2)**2
        A_l_x(2,2)=Hb_x_x*(1-x0**2)**2
     &            -2*x0*2*(1-x0**2)*Hb_x0

        ! give values to the ads gh source functions
        Hads_l(2)=(-2+4*x0)/(1-2*x0*(1-x0))

        ! give values to the ads gh source function derivatives
        Hads_l_x(2,2)=8*x0*(1-x0)/(1-2*x0*(1-x0))**2

        ! give values to the scalar field
        phi10_x(1)=phi1_t*(1-x0**2)**3
        phi10_x(2)=phi1_x*(1-x0**2)**3
     &            +phi10*(-6*x0)*(1-x0**2)**2
        phi10_x(3)=0

        phi10_xx(1,1)=phi1_tt*(1-x0**2)**3
        phi10_xx(1,2)=phi1_tx*(1-x0**2)**3
     &               +phi1_t*(-6*x0)*(1-x0**2)**2
        phi10_xx(1,3)=0
        phi10_xx(2,2)=phi1_xx*(1-x0**2)**3
     &               +phi1_x*(2)*(-6*x0)*(1-x0**2)**2
     &               +phi10*(-6*(1-x0**2)**2+24*x0**2*(1-x0**2))
        phi10_xx(2,3)=0
        phi10_xx(3,3)=0                

        do a=1,2
          do b=a+1,3
            phi10_xx(b,a)=phi10_xx(a,b)
          end do
        end do

        ! define Levi Civita symbol
        do a=1,3
          do b=1,3
            do c=1,3
                if ((a.eq.1.and.b.eq.2.and.c.eq.3)
     &          .or.(a.eq.3.and.b.eq.1.and.c.eq.2)
     &          .or.(a.eq.2.and.b.eq.3.and.c.eq.1))
     &          then
                  levicivi3(a,b,c)=+1
                end if
                if ((a.eq.2.and.b.eq.1.and.c.eq.3)
     &          .or.(a.eq.1.and.b.eq.3.and.c.eq.2)
     &          .or.(a.eq.3.and.b.eq.2.and.c.eq.1))
     &          then
                  levicivi3(a,b,c)=-1
                end if
            end do
          end do
        end do

       ! define square root of determinant of metric in (t,x,y) subsector
       sqrtdetg=sqrt(abs(g0_ll(1,1)*g0_ll(2,2)*g0_ll(3,3)
     &                  +g0_ll(1,2)*g0_ll(2,3)*g0_ll(3,1)
     &                  +g0_ll(1,3)*g0_ll(2,1)*g0_ll(3,2)
     &                  -g0_ll(1,3)*g0_ll(2,2)*g0_ll(3,1)
     &                  -g0_ll(1,2)*g0_ll(2,1)*g0_ll(3,3)
     &                  -g0_ll(1,1)*g0_ll(2,3)*g0_ll(3,2)))

       ! define volume form in in (t,x,y) subsector
        do a=1,3
          do b=1,3
            do c=1,3
              vol(a,b,c)=levicivi3(a,b,c)*sqrtdetg
            end do
          end do
        end do

        ! give values to the field strength, using sin(phi2)=sin(phi3)=sin(phi4)=1 w.l.o.g 
        !(considering phi2,phi3,phi4-independent case, so phi2=phi3=phi4=pi/2 slice will do)
        f1_l(1)   =0        !+fb_t0
        f1_l(2)   =0        !+fb_x0 !NOTE: add this when you add f evolution
        f1_l(3)   =f1_y_ads0!+fb_y0
!        do a=1,3
!          do b=1,3
!            fb_tx0=-vol(1,2,a)*f1_l(b)*g0_uu(a,b)-f2_tx_ads0   !check that this is zero for pure
!            fb_ty0=-vol(1,3,a)*f1_l(b)*g0_uu(a,b)
!            fb_xy0=-vol(2,3,a)*f1_l(b)*g0_uu(a,b)
!          end do
!        end do
        f2_ll(1,2)=f2_tx_ads0!+fb_tx0 
        f2_ll(1,3)=0         !+fb_ty0 !NOTE: add this when you add f evolution
        f2_ll(2,3)=0         !+fb_xy0

        do a=1,2
          do b=a+1,3
            f2_ll(b,a)=-f2_ll(a,b)
          end do
        end do

        ! calculate Christoffel symbol derivatives at point i
        !(gamma^a_bc,e = 1/2 g^ad_,e(g_bd,c  + g_cd,b  - g_bc,d)
        !              +   1/2 g^ad(g_bd,ce + g_cd,be - g_bc,de))
        !(WARNING: 
        ! this second derivative needs more info on chi,theta, so 
        ! cannot use the g0_ll_x, g0_uu_x, g0_ll_xx set by chi=pi/2,theta=pi/2, 
        ! only (from pure AdS piece):
        ! gamma_ull(4,4,3)=cotchi
        ! gamma_ull(5,5,3)=cotchi
        ! gamma_ull(2,4,4)=-1/L^2*x/(1-x)*(L^2*(1-x)^2+x^2)*sinchi^2
        ! gamma_ull(3,4,4)=-coschi*sinchi
        ! gamma_ull(4,3,4)=cotchi
        ! gamma_ull(5,5,4)=cottheta
        ! gamma_ull(2,5,5)=-1/L^2*x/(1-x)*(L^2*(1-x)^2+x^2)*sinchi^2*sintheta^2
        ! gamma_ull(3,5,5)=-coschi*sinchi*sintheta^2
        ! gamma_ull(4,5,5)=-costheta*sintheta
        ! gamma_ull(5,4,5)=cottheta
        ! gamma_ull(5,3,5)=cotchi
        ! are affected in this AdS5D case, so hardcoded these at the end)
        do a=1,3
          do b=1,3
            do c=1,3
              do e=1,3
                gamma_ull_x(a,b,c,e)=0
                do d=1,3
                  gamma_ull_x(a,b,c,e)=gamma_ull_x(a,b,c,e)
     &              +0.5d0*g0_uu_x(a,d,e)*(g0_ll_x(b,d,c)+
     &                     g0_ll_x(c,d,b)-g0_ll_x(b,c,d))
     &              +0.5d0*g0_uu(a,d)*(g0_ll_xx(b,d,c,e)+
     &                     g0_ll_xx(c,d,b,e)-g0_ll_xx(b,c,d,e))
                end do
              end do
            end do
          end do
        end do
        !gamma_ull_x(4,4,3,3)=-1 
        !gamma_ull_x(5,5,3,3)=-1 
        !gamma_ull_x(2,4,4,3)=0
        !gamma_ull_x(3,4,4,3)=1
        !gamma_ull_x(4,3,4,3)=-1
        !gamma_ull_x(5,5,4,4)=-1
        !gamma_ull_x(2,5,5,3)=0
        !gamma_ull_x(2,5,5,4)=0
        !gamma_ull_x(3,5,5,3)=1
        !gamma_ull_x(3,5,5,4)=0
        !gamma_ull_x(4,5,5,4)=1
        !gamma_ull_x(5,4,5,4)=-1
        !gamma_ull_x(5,3,5,3)=-1

        ! calculate Riemann tensor at point i
        !(R^a_bcd =gamma^a_bd,c - gamma^a_bc,d
        !          +gamma^a_ce gamma^e_bd - gamma^a_de gamma^e_bc)
        do a=1,3
          do b=1,3
            do c=1,3
              do d=1,3
                riemann_ulll(a,b,c,d)=
     &                gamma_ull_x(a,b,d,c)-gamma_ull_x(a,b,c,d)
                do e=1,3
                   riemann_ulll(a,b,c,d)=riemann_ulll(a,b,c,d)
     &               +gamma_ull(a,c,e)*gamma_ull(e,b,d)
     &               -gamma_ull(a,d,e)*gamma_ull(e,b,c)
                end do
              end do
            end do
          end do
        end do

        ! calculate Ricci tensor at point i
        !(R_bd = R^a_bad)
        do b=1,3
          do d=1,3
            riccibar_ll(b,d)=0
            do a=1,3
              riccibar_ll(b,d)=riccibar_ll(b,d)+riemann_ulll(a,b,a,d)
            end do
          end do
        end do

        ! calculate raised Ricci tensor at point i
        !(R_a^b = R_ad g^db)
        do a=1,3
          do b=1,3
            riccibar_lu(a,b)=0
            do d=1,3
              riccibar_lu(a,b)=riccibar_lu(a,b)
     &                        +riccibar_ll(a,d)*g0_uu(d,b)
            end do
          end do
        end do

        ! calculate Ricci scalar
        !(R = R_a^a)
        riccibar=0
        do a=1,3
          riccibar=riccibar+riccibar_lu(a,a)
        end do
  
        ! calculates stress-energy tensor at point i 
        !(T_ab = 2*phi1,a phi1,b - (phi1,c phi1,d) g^cd g_ab + ...)
        grad_phi1_sq=0
        do a=1,3
          do b=1,3
            grad_phi1_sq=grad_phi1_sq
     &                  +phi10_x(a)*phi10_x(b)*g0_uu(a,b)
          end do
        end do

        ! calculate s0_ll,t0_ll,ricci_ll
        do a=1,3
          do b=1,3
            s0_ll(a,b)=-dimA*( (gA_xx(a,b)
     &                         -gamma_ull(1,b,a)*gA_x(1)
     &                         -gamma_ull(2,b,a)*gA_x(2)
     &                         -gamma_ull(3,b,a)*gA_x(3)
     &                         )/(2*gA)
     &                       - (gA_x(a)*gA_x(b))/(4*gA**2) )
     &                 -dimB*( (gB_xx(a,b)
     &                         -gamma_ull(1,b,a)*gB_x(1)
     &                         -gamma_ull(2,b,a)*gB_x(2)
     &                         -gamma_ull(3,b,a)*gB_x(3)
     &                         )/(2*gB)
     &                       - (gB_x(a)*gB_x(b))/(4*gB**2) )

            t0_ll(a,b)=-(f1_l(a)*f1_l(b)
     &                  +g0_uu(1,1)*f2_ll(a,1)*f2_ll(b,1)
     &                  +g0_uu(1,2)*f2_ll(a,1)*f2_ll(b,2)
     &                  +g0_uu(1,3)*f2_ll(a,1)*f2_ll(b,3)
     &                  +g0_uu(2,1)*f2_ll(a,2)*f2_ll(b,1)
     &                  +g0_uu(2,2)*f2_ll(a,2)*f2_ll(b,2)
     &                  +g0_uu(2,3)*f2_ll(a,2)*f2_ll(b,3)
     &                  +g0_uu(3,1)*f2_ll(a,3)*f2_ll(b,1)
     &                  +g0_uu(3,2)*f2_ll(a,3)*f2_ll(b,2)
     &                  +g0_uu(3,3)*f2_ll(a,3)*f2_ll(b,3)
     &                  )/4

            ricci_ll(a,b)=riccibar_ll(a,b)+s0_ll(a,b)
          end do
        end do

!!TEST!
!              do a=1,3
!                do b=a,3
!
!                  term1(a,b)=-0.5d0*(
!     &                          g0_uu(1,1)*g0_ll_xx(a,b,1,1)+
!     &                          g0_uu(2,2)*g0_ll_xx(a,b,2,2)+
!     &                          g0_uu(3,3)*g0_ll_xx(a,b,3,3)+
!     &                       2*(g0_uu(1,2)*g0_ll_xx(a,b,1,2)+
!     &                          g0_uu(1,3)*g0_ll_xx(a,b,1,3)+
!     &                          g0_uu(2,3)*g0_ll_xx(a,b,2,3))
!     &                              )
!     &
!                  term2(a,b)=-0.5d0*(
!     &                          g0_uu_x(1,1,a)* g0_ll_x(b,1,1) +
!     &                          g0_uu_x(1,2,a)*(g0_ll_x(b,1,2) +
!     &                                          g0_ll_x(b,2,1))+
!     &                          g0_uu_x(1,3,a)*(g0_ll_x(b,1,3) +
!     &                                          g0_ll_x(b,3,1))+
!     &                          g0_uu_x(2,2,a)* g0_ll_x(b,2,2) +
!     &                          g0_uu_x(2,3,a)*(g0_ll_x(b,2,3) +
!     &                                          g0_ll_x(b,3,2))+
!     &                          g0_uu_x(3,3,a)* g0_ll_x(b,3,3)
!     &                            )
!     &
!                  term3(a,b)=-0.5d0*(
!     &                          g0_uu_x(1,1,b)* g0_ll_x(a,1,1) +
!     &                          g0_uu_x(1,2,b)*(g0_ll_x(a,1,2) +
!     &                                          g0_ll_x(a,2,1))+
!     &                          g0_uu_x(1,3,b)*(g0_ll_x(a,1,3) +
!     &                                          g0_ll_x(a,3,1))+
!     &                          g0_uu_x(2,2,b)* g0_ll_x(a,2,2) +
!     &                          g0_uu_x(2,3,b)*(g0_ll_x(a,2,3) +
!     &                                          g0_ll_x(a,3,2))+
!     &                          g0_uu_x(3,3,b)* g0_ll_x(a,3,3)
!     &                            )
!     &
!                  term4(a,b)=-0.5d0*(Hads_l_x(a,b)+A_l_x(a,b))
!     &
!                  term5(a,b)=-0.5d0*(Hads_l_x(b,a)+A_l_x(b,a))
!     &
!                  term6(a,b)=     (
!     &                          (Hads_l(1)+A_l(1))*gamma_ull(1,a,b)+
!     &                          (Hads_l(2)+A_l(2))*gamma_ull(2,a,b)+
!     &                          (Hads_l(3)+A_l(3))*gamma_ull(3,a,b)
!     &                            )
!     &
!                  term7(a,b)=    -(
!     &                          gamma_ull(1,1,b)*gamma_ull(1,1,a)+
!     &                          gamma_ull(1,2,b)*gamma_ull(2,1,a)+
!     &                          gamma_ull(1,3,b)*gamma_ull(3,1,a)+
!     &                          gamma_ull(2,1,b)*gamma_ull(1,2,a)+
!     &                          gamma_ull(2,2,b)*gamma_ull(2,2,a)+
!     &                          gamma_ull(2,3,b)*gamma_ull(3,2,a)+
!     &                          gamma_ull(3,1,b)*gamma_ull(1,3,a)+
!     &                          gamma_ull(3,2,b)*gamma_ull(2,3,a)+
!     &                          gamma_ull(3,3,b)*gamma_ull(3,3,a)
!     &                            )
!     &
!                  term8(a,b)=s0_ll(a,b)
!     &
!                  term9(a,b)=t0_ll(a,b)
!     &
!                  efe(a,b)=term1(a,b)+term2(a,b)+term3(a,b)+term4(a,b)
!     &                    +term5(a,b)+term6(a,b)+term7(a,b)+term8(a,b)
!     &                    +term9(a,b)
!     &
!                  efest(a,b)=riccibar_ll(a,b)+s0_ll(a,b)+t0_ll(a,b)
!
!                end do
!              end do
!
!              efe_norm=
!     &        max(abs(efe(1,1)),abs(efe(1,2)),
!     &            abs(efe(2,2)),abs(efe(3,3)))
!              efest_norm=
!     &        max(abs(efest(1,1)),abs(efest(1,2)),
!     &            abs(efest(2,2)),abs(efest(3,3)))
!
!!TEST!
!        if (x0.eq.0.5d0) then
!          write(*,*) 'efe_norm=',efe_norm
!          write(*,*) 'efest_norm=',efest_norm
!          write(*,*) 'g0_ll_x(2,2,1)=',g0_ll_x(2,2,1)
!          write(*,*) 'g0_ll_xx(2,2,1,1)=',g0_ll_xx(2,2,1,1)
!          write(*,*) 'riccibar_ll(1,1)=',riccibar_ll(1,1)
!          write(*,*) 'riccibar_ll(1,2)=',riccibar_ll(1,2)
!          write(*,*) 'riccibar_ll(2,2)=',riccibar_ll(2,2)
!          write(*,*) 'riccibar_ll(3,3)=',riccibar_ll(3,3)
!        end if
!!TEST!

        return
        end
