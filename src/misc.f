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
        subroutine init_schw(gb_tt,gb_tx,gb_xx,psi,r0,
     &                       L,phys_bdy,chr,ex,x,Nx)
        implicit none
        real*8 gb_tt(Nx),gb_tx(Nx),gb_xx(Nx),psi(Nx)
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
              gb_tx(i)=1/(1+x0)**2/(1-x_h)**5*(1-x0)
              gb_xx(i)=-((L**2*(1 - x0)*
     &                ((1 - x0)**6 - 4*(1 - x0)**7 + 
     &                (6 + L**2)*(1 - x0)**8 - 
     &                2*(2 + L**2)*(1 - x0)**9 + 
     &                (1 + L**2)*(1 - x0)**10 - 
     &                L**2*r0**2*(-1 + x_h)**10))/
     &                ((1 + (1 - x0)*(-1 + L**2*(1 - x0) - x0))*
     &                (-1 - x0)*(L**2*r0**2*(1 - x0)**4 - 
     &                (1 + (1 - x0)*(-1 + L**2*(1 - x0) - x0))*
     &                x0**2)*(-1 + x_h)**10))
              psi(i)=0

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
        call axi_reg_g(gb_tt,gb_tx,gb_xx,psi,chr,ex,L,x,Nx)

        return
        end

c----------------------------------------------------------------------
c calculates all the tensorial objects in x coordinates, at point i
c----------------------------------------------------------------------
        subroutine tensor_init(
     &                  gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                  gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                  gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                  psi_np1,psi_n,psi_nm1,
     &                  Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &                  Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &                  phi1_np1,phi1_n,phi1_nm1,
     &                  g0_ll,g0_uu,g0_ll_x,g0_uu_x,g0_ll_xx,
     &                  gads_ll,gads_uu,gads_ll_x,gads_uu_x,gads_ll_xx,
     &                  h0_ll,h0_uu,h0_ll_x,h0_uu_x,h0_ll_xx,
     &                  A_l,A_l_x,Hads_l,
     &                  gamma_ull,gamma_ull_x,
     &                  riemann_ulll,ricci_ll,ricci_lu,ricci,
     &                  einstein_ll,set_ll,
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
        real*8 psi_np1(Nx),psi_n(Nx),psi_nm1(Nx)
        real*8 Hb_t_np1(Nx),Hb_t_n(Nx),Hb_t_nm1(Nx)
        real*8 Hb_x_np1(Nx),Hb_x_n(Nx),Hb_x_nm1(Nx)
        real*8 phi1_np1(Nx),phi1_n(Nx),phi1_nm1(Nx)
        real*8 fb_t_np1(Nx),fb_t_n(Nx),fb_t_nm1(Nx)
        real*8 fb_x_np1(Nx),fb_x_n(Nx),fb_x_nm1(Nx)
        real*8 fb_y_np1(Nx),fb_y_n(Nx),fb_y_nm1(Nx)

        integer a,b,c,d,e,m,n,o,p

        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 dx,dy
        real*8 x0,y0   

        real*8 grad_phi1_sq

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
        real*8 f0_l(5),f0_ll(5,5),ff_ll(5,5)
        real*8 f_lllll(5,5,5,5,5),f_luuuu(5,5,5,5,5)
        real*8 fads_l(5),fads_ll(5,5)

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
        real*8 psi_t,psi_x
        real*8 psi_tt,psi_tx,psi_xx
        real*8 phi1_t,phi1_x
        real*8 phi1_tt,phi1_tx,phi1_xx
  
        real*8 gb_tt0,gb_tx0,gb_xx0,psi0,phi10
        real*8 g0_tt_ads0,g0_xx_ads0,g0_psi_ads0

        real*8 g0_tt_ads_x,g0_tt_ads_xx
        real*8 g0_xx_ads_x,g0_xx_ads_xx
        real*8 g0_psi_ads_x,g0_psi_ads_xx

        real*8 Hb_t_t,Hb_t_x
        real*8 Hb_x_t,Hb_x_x

        real*8 Hb_t0,Hb_x0

        real*8 fb_t0,fb_x0,fb_y0
        real*8 f0_tx_ads0,f0_y_ads0
!----------------------------------------------------------------------
        
        dx=(x(2)-x(1))

        x0=x(i)

        ! set gads values using sin(theta1)=sin(theta2)=1 w.l.o.g 
        !(considering theta1,theta2-independent case, so theta1=theta2=pi/2 slice will do)
        g0_tt_ads0 =-((1-x0)**2+x0**2)/(1-x0)**2
        g0_xx_ads0 =1/((1-x0)**2+x0**2)/(1-x0)**2
        g0_psi_ads0=x0**2/(1-x0)**2

        ! set fads values using sin(phi2)=sin(phi3)=sin(phi4)=1 w.l.o.g 
        !(considering phi2,phi3,phi4-independent case, so phi2=phi3=phi4=pi/2 slice will do)
        f0_tx_ads0 = 4/L*(-16*x0**3*(1+x0**2)/(1-x0**2)**5)
        f0_y_ads0  = 4/L*(PI*L**4*sin(PI*y0/L)**4)

        ! set gbar values
        gb_tt0=gb_tt_n(i)
        gb_tx0=gb_tx_n(i)
        gb_xx0=gb_xx_n(i)
        psi0  =psi_n(i)

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
        g0_psi_ads_x =2*x0/(1-x0)**3
        g0_psi_ads_xx=(2+4*x0)/(1-x0)**4
 
        ! calculate gbar derivatives
        call df2_int(gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &       gb_tt_t,gb_tt_x,
     &       gb_tt_tt,gb_tt_tx,gb_tt_xx,
     &       dx,dt,i,chr,ex,Nx,'phi1')
        call df2_int(gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &       gb_tx_t,gb_tx_x,
     &       gb_tx_tt,gb_tx_tx,gb_tx_xx,
     &       dx,dt,i,chr,ex,Nx,'phi1')
        call df2_int(gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &       gb_xx_t,gb_xx_x,
     &       gb_xx_tt,gb_xx_tx,gb_xx_xx,
     &       dx,dt,i,chr,ex,Nx,'phi1')
        call df2_int(psi_np1,psi_n,psi_nm1,
     &       psi_t,psi_x,
     &       psi_tt,psi_tx,psi_xx,
     &       dx,dt,i,chr,ex,Nx,'phi1')

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
        g0_ll(2,2)=g0_xx_ads0+gb_xx0*(1-x0**2)
        g0_ll(3,3)=g0_psi_ads0+psi0*(1-x0**2)*x0**2
        g0_ll(4,4)=g0_psi_ads0+psi0*(1-x0**2)*x0**2
        g0_ll(5,5)=g0_psi_ads0+psi0*(1-x0**2)*x0**2

        g0_uu(1,1)=
     &          (g0_ll(2,2)*g0_ll(3,3)*g0_ll(4,4)*g0_ll(5,5))
     &         /(g0_ll(1,1)*g0_ll(2,2)*g0_ll(3,3)*g0_ll(4,4)*g0_ll(5,5)
     &                  -g0_ll(1,2)**2*g0_ll(3,3)*g0_ll(4,4)*g0_ll(5,5))
        g0_uu(1,2)=
     &         (-g0_ll(1,2)*g0_ll(3,3)*g0_ll(4,4)*g0_ll(5,5))
     &         /(g0_ll(1,1)*g0_ll(2,2)*g0_ll(3,3)*g0_ll(4,4)*g0_ll(5,5)
     &                  -g0_ll(1,2)**2*g0_ll(3,3)*g0_ll(4,4)*g0_ll(5,5))
        g0_uu(2,2)=
     &          (g0_ll(1,1)*g0_ll(3,3)*g0_ll(4,4)*g0_ll(5,5))
     &         /(g0_ll(1,1)*g0_ll(2,2)*g0_ll(3,3)*g0_ll(4,4)*g0_ll(5,5)
     &                  -g0_ll(1,2)**2*g0_ll(3,3)*g0_ll(4,4)*g0_ll(5,5))
        g0_uu(3,3)=
     &          (g0_ll(1,1)*g0_ll(2,2)*g0_ll(4,4)*g0_ll(5,5)
     &                  -g0_ll(1,2)**2*g0_ll(4,4)*g0_ll(5,5))
     &         /(g0_ll(1,1)*g0_ll(2,2)*g0_ll(3,3)*g0_ll(4,4)*g0_ll(5,5)
     &                  -g0_ll(1,2)**2*g0_ll(3,3)*g0_ll(4,4)*g0_ll(5,5))
        g0_uu(4,4)=
     &          (g0_ll(1,1)*g0_ll(2,2)*g0_ll(3,3)*g0_ll(5,5)
     &                  -g0_ll(1,2)**2*g0_ll(3,3)*g0_ll(5,5))
     &         /(g0_ll(1,1)*g0_ll(2,2)*g0_ll(3,3)*g0_ll(4,4)*g0_ll(5,5)
     &                  -g0_ll(1,2)**2*g0_ll(3,3)*g0_ll(4,4)*g0_ll(5,5))
        g0_uu(5,5)=
     &          (g0_ll(1,1)*g0_ll(2,2)*g0_ll(3,3)*g0_ll(4,4)
     &                  -g0_ll(1,2)**2*g0_ll(3,3)*g0_ll(4,4))
     &         /(g0_ll(1,1)*g0_ll(2,2)*g0_ll(3,3)*g0_ll(4,4)*g0_ll(5,5)
     &                  -g0_ll(1,2)**2*g0_ll(3,3)*g0_ll(4,4)*g0_ll(5,5))

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
     &                   +psi_t*(1-x0**2)*x0**2
        g0_ll_x(3,3,2)   =g0_psi_ads_x
     &                   +psi_x*(1-x0**2)*x0**2
     &                   +psi0*(2*x0-4*x0**3)
        g0_ll_xx(3,3,1,1)=0
     &                   +psi_tt*(1-x0**2)*x0**2
        g0_ll_xx(3,3,1,2)=0
     &                   +psi_tx*(1-x0**2)*x0**2
     &                   +psi_t*(2*x0-4*x0**3)
        g0_ll_xx(3,3,2,2)=g0_psi_ads_xx
     &                   +psi_xx*(1-x0**2)*x0**2
     &                   +psi_x*(2*x0-4*x0**3)
     &                   +psi_x*(2*x0-4*x0**3)
     &                   +psi0*(2-12*x0**2)

        g0_ll_x(4,4,1)   =0
     &                   +psi_t*(1-x0**2)*x0**2
        g0_ll_x(4,4,2)   =g0_psi_ads_x
     &                   +psi_x*(1-x0**2)*x0**2
     &                   +psi0*(2*x0-4*x0**3)
        g0_ll_xx(4,4,1,1)=0
     &                   +psi_tt*(1-x0**2)*x0**2
        g0_ll_xx(4,4,1,2)=0
     &                   +psi_tx*(1-x0**2)*x0**2
     &                   +psi_t*(2*x0-4*x0**3)
        g0_ll_xx(4,4,2,2)=g0_psi_ads_xx
     &                   +psi_xx*(1-x0**2)*x0**2
     &                   +psi_x*(2*x0-4*x0**3)
     &                   +psi_x*(2*x0-4*x0**3)
     &                   +psi0*(2-12*x0**2)
        g0_ll_xx(4,4,3,3)=-2*g0_psi_ads0          ! WARNING: from sin^2chi factor in pure ads term

        g0_ll_x(5,5,1)   =0        
     &                   +psi_t*(1-x0**2)*x0**2
        g0_ll_x(5,5,2)   =g0_psi_ads_x
     &                   +psi_x*(1-x0**2)*x0**2
     &                   +psi0*(2*x0-4*x0**3)
        g0_ll_xx(5,5,1,1)=0
     &                   +psi_tt*(1-x0**2)*x0**2
        g0_ll_xx(5,5,1,2)=0
     &                   +psi_tx*(1-x0**2)*x0**2
     &                   +psi_t*(2*x0-4*x0**3)
        g0_ll_xx(5,5,2,2)=g0_psi_ads_xx
     &                   +psi_xx*(1-x0**2)*x0**2
     &                   +psi_x*(2*x0-4*x0**3)
     &                   +psi_x*(2*x0-4*x0**3)
     &                   +psi0*(2-12*x0**2)
        g0_ll_xx(5,5,3,3)=-2*g0_psi_ads0          ! WARNING: from sin^2chi factor in pure ads term
        g0_ll_xx(5,5,4,4)=-2*g0_psi_ads0          ! WARNING: from sin^2theta factor in pure ads term

        do a=1,4
          do b=a+1,5
            g0_ll(b,a)=g0_ll(a,b)
            g0_uu(b,a)=g0_uu(a,b) 
            do c=1,5
              g0_ll_x(b,a,c)=g0_ll_x(a,b,c)
            end do
          end do
        end do

        do a=1,5
          do b=1,5
            do c=1,5
              g0_uu_x(a,b,c)=0
              do d=1,5
                do e=1,5
                  g0_uu_x(a,b,c)=g0_uu_x(a,b,c)
     &                          -g0_ll_x(d,e,c)
     &                           *g0_uu(a,d)*g0_uu(b,e)
                end do
     &  
              end do
            end do
          end do
        end do

        do a=1,5
          do b=1,5
            do c=1,5
              do d=1,5
                g0_ll_xx(a,b,c,d)=
     &             g0_ll_xx(min(a,b),max(a,b),min(c,d),max(c,d))
              end do
            end do
          end do
        end do

        do a=1,5
          do b=1,5
            do c=1,5
              gamma_ull(a,b,c)=0
              do d=1,5
                gamma_ull(a,b,c)=gamma_ull(a,b,c)
     &                          +0.5d0*g0_uu(a,d)
     &                                *(g0_ll_x(c,d,b)
     &                                 -g0_ll_x(b,c,d)
     &                                 +g0_ll_x(d,b,c))
              end do
            end do
          end do
        end do

        ! give values to the ads metric, using sin(theta1)=sin(theta2)=1 w.l.o.g 
        !(considering theta1,theta2-independent case, so theta1=theta2=pi/2 slice will do)
        gads_ll(1,1)=g0_tt_ads0
        gads_ll(2,2)=g0_xx_ads0
        gads_ll(3,3)=g0_psi_ads0
        gads_ll(4,4)=g0_psi_ads0
        gads_ll(5,5)=g0_psi_ads0

        gads_uu(1,1)=1/g0_tt_ads0
        gads_uu(2,2)=1/g0_xx_ads0
        gads_uu(3,3)=1/g0_psi_ads0
        gads_uu(4,4)=1/g0_psi_ads0
        gads_uu(5,5)=1/g0_psi_ads0

        gads_ll_x(1,1,2)   =g0_tt_ads_x  
        gads_ll_xx(1,1,2,2)=g0_tt_ads_xx 
        gads_ll_x(2,2,2)   =g0_xx_ads_x  
        gads_ll_xx(2,2,2,2)=g0_xx_ads_xx 
        gads_ll_x(3,3,2)   =g0_psi_ads_x 
        gads_ll_xx(3,3,2,2)=g0_psi_ads_xx
        gads_ll_x(4,4,2)   =g0_psi_ads_x 
        gads_ll_xx(4,4,2,2)=g0_psi_ads_xx
        gads_ll_xx(4,4,3,3)=-2*g0_psi_ads0          ! WARNING: from sin^2chi factor in pure ads term
        gads_ll_x(5,5,2)   =g0_psi_ads_x 
        gads_ll_xx(5,5,2,2)=g0_psi_ads_xx
        gads_ll_xx(5,5,3,3)=-2*g0_psi_ads0          ! WARNING: from sin^2chi factor in pure ads term
        gads_ll_xx(5,5,4,4)=-2*g0_psi_ads0          ! WARNING: from sin^2theta factor in pure ads term
                
        do a=1,4
          do b=a+1,5
            gads_ll(b,a)=gads_ll(a,b)
            gads_uu(b,a)=gads_uu(a,b)
            do c=1,5
              gads_ll_x(b,a,c)=gads_ll_x(a,b,c)
            end do
          end do
        end do

        do a=1,5
          do b=1,5
            do c=1,5
              gads_uu_x(a,b,c)=
     &              -gads_ll_x(1,1,c)*gads_uu(a,1)*gads_uu(b,1)
     &              -gads_ll_x(1,2,c)*(gads_uu(a,1)*gads_uu(b,2)
     &                               +gads_uu(a,2)*gads_uu(b,1))
     &              -gads_ll_x(1,3,c)*(gads_uu(a,1)*gads_uu(b,3)
     &                               +gads_uu(a,3)*gads_uu(b,1))
     &              -gads_ll_x(1,4,c)*(gads_uu(a,1)*gads_uu(b,4)
     &                               +gads_uu(a,4)*gads_uu(b,1))
     &              -gads_ll_x(1,5,c)*(gads_uu(a,1)*gads_uu(b,5)
     &                               +gads_uu(a,5)*gads_uu(b,1))
     &              -gads_ll_x(2,2,c)*gads_uu(a,2)*gads_uu(b,2)
     &              -gads_ll_x(2,3,c)*(gads_uu(a,2)*gads_uu(b,3)
     &                               +gads_uu(a,3)*gads_uu(b,2))
     &              -gads_ll_x(2,4,c)*(gads_uu(a,2)*gads_uu(b,4)
     &                               +gads_uu(a,4)*gads_uu(b,2))
     &              -gads_ll_x(2,5,c)*(gads_uu(a,2)*gads_uu(b,5)
     &                               +gads_uu(a,5)*gads_uu(b,2))
     &              -gads_ll_x(3,3,c)*gads_uu(a,3)*gads_uu(b,3)
     &              -gads_ll_x(3,4,c)*(gads_uu(a,3)*gads_uu(b,4)
     &                               +gads_uu(a,4)*gads_uu(b,3))
     &              -gads_ll_x(3,5,c)*(gads_uu(a,3)*gads_uu(b,5)
     &                               +gads_uu(a,5)*gads_uu(b,3))
     &              -gads_ll_x(4,4,c)*gads_uu(a,4)*gads_uu(b,4)
     &              -gads_ll_x(4,5,c)*(gads_uu(a,4)*gads_uu(b,5)
     &                               +gads_uu(a,5)*gads_uu(b,4))
     &              -gads_ll_x(5,5,c)*gads_uu(a,5)*gads_uu(b,5)
            end do
          end do
        end do

        ! give values to the metric deviation, using sin(theta1)=sin(theta2)=1 w.l.o.g 
        !(considering theta1,theta2-independent case, so theta1=theta2=pi/2 will do)
        h0_ll(1,1)=gb_tt0*(1-x0**2)
        h0_ll(1,2)=gb_tx0*(1-x0**2)**2
        h0_ll(2,2)=gb_xx0*(1-x0**2)
        h0_ll(3,3)=psi0*x0**2*(1-x0**2)
        h0_ll(4,4)=psi0*x0**2*(1-x0**2)
        h0_ll(5,5)=psi0*x0**2*(1-x0**2) 
        
        h0_uu(1,1)=g0_uu(1,1)-gads_uu(1,1)
        h0_uu(1,2)=g0_uu(1,2)
        h0_uu(2,2)=g0_uu(2,2)-gads_uu(2,2)
        h0_uu(3,3)=g0_uu(3,3)-gads_uu(3,3)
        h0_uu(4,4)=g0_uu(4,4)-gads_uu(4,4)
        h0_uu(5,5)=g0_uu(4,4)-gads_uu(4,4)

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

        h0_ll_x(4,4,1)   =g0_ll_x(4,4,1)-gads_ll_x(4,4,1)
        h0_ll_x(4,4,2)   =g0_ll_x(4,4,2)-gads_ll_x(4,4,2)
        h0_ll_xx(4,4,1,1)=g0_ll_xx(4,4,1,1)-gads_ll_xx(4,4,1,1)
        h0_ll_xx(4,4,1,2)=g0_ll_xx(4,4,1,2)-gads_ll_xx(4,4,1,2)
        h0_ll_xx(4,4,2,2)=g0_ll_xx(4,4,2,2)-gads_ll_xx(4,4,2,2)
        h0_ll_xx(4,4,3,3)=g0_ll_xx(4,4,3,3)-gads_ll_xx(4,4,3,3)

        h0_ll_x(5,5,1)   =g0_ll_x(5,5,1)-gads_ll_x(5,5,1)
        h0_ll_x(5,5,2)   =g0_ll_x(5,5,2)-gads_ll_x(5,5,2)
        h0_ll_xx(5,5,1,1)=g0_ll_xx(5,5,1,1)-gads_ll_xx(5,5,1,1)
        h0_ll_xx(5,5,1,2)=g0_ll_xx(5,5,1,2)-gads_ll_xx(5,5,1,2)
        h0_ll_xx(5,5,2,2)=g0_ll_xx(5,5,2,2)-gads_ll_xx(5,5,2,2)
        h0_ll_xx(5,5,3,3)=g0_ll_xx(5,5,3,3)-gads_ll_xx(5,5,3,3)
        h0_ll_xx(5,5,4,4)=g0_ll_xx(5,5,4,4)-gads_ll_xx(5,5,4,4)

        h0_uu_x(1,1,1)=g0_uu_x(1,1,1)-gads_uu_x(1,1,1)
        h0_uu_x(1,1,2)=g0_uu_x(1,1,2)-gads_uu_x(1,1,2)

        h0_uu_x(1,2,1)=g0_uu_x(1,2,1)-gads_uu_x(1,2,1)
        h0_uu_x(1,2,2)=g0_uu_x(1,2,2)-gads_uu_x(1,2,2)

        h0_uu_x(2,2,1)=g0_uu_x(2,2,1)-gads_uu_x(2,2,1)
        h0_uu_x(2,2,2)=g0_uu_x(2,2,2)-gads_uu_x(2,2,2)

        h0_uu_x(3,3,1)=g0_uu_x(3,3,1)-gads_uu_x(3,3,1)
        h0_uu_x(3,3,2)=g0_uu_x(3,3,2)-gads_uu_x(3,3,2)

        h0_uu_x(4,4,1)=g0_uu_x(4,4,1)-gads_uu_x(4,4,1)
        h0_uu_x(4,4,2)=g0_uu_x(4,4,2)-gads_uu_x(4,4,2)

        h0_uu_x(5,5,1)=g0_uu_x(5,5,1)-gads_uu_x(5,5,1)
        h0_uu_x(5,5,2)=g0_uu_x(5,5,2)-gads_uu_x(5,5,2)

        do a=1,4
          do b=a+1,5
            h0_ll(b,a)=h0_ll(a,b)
            h0_uu(b,a)=h0_uu(a,b)
            do c=1,5
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
        Hads_l(1)=0
        Hads_l(2)=(3-4*x0*(2-x0)*(1-x0))/(x0*(1-x0)*(1-2*x0*(1-x0)))

        ! give values to the scalar field
        phi10_x(1)=phi1_t*(1-x0**2)**3
        phi10_x(2)=phi1_x*(1-x0**2)**3
     &            +phi10*(-6*x0)*(1-x0**2)**2
        phi10_x(3)=0
        phi10_x(4)=0
        phi10_x(5)=0

        phi10_xx(1,1)=phi1_tt*(1-x0**2)**3
        phi10_xx(1,2)=phi1_tx*(1-x0**2)**3
     &               +phi1_t*(-6*x0)*(1-x0**2)**2
        phi10_xx(1,3)=0
        phi10_xx(1,4)=0
        phi10_xx(1,5)=0
        phi10_xx(2,2)=phi1_xx*(1-x0**2)**3
     &               +phi1_x*(2)*(-6*x0)*(1-x0**2)**2
     &               +phi10*(-6*(1-x0**2)**2+24*x0**2*(1-x0**2))
        phi10_xx(2,3)=0
        phi10_xx(2,4)=0
        phi10_xx(2,5)=0
        phi10_xx(3,3)=0                
        phi10_xx(3,4)=0
        phi10_xx(3,5)=0
        phi10_xx(4,4)=0
        phi10_xx(4,5)=0
        phi10_xx(5,5)=0

        do a=1,4
          do b=a+1,5
            phi10_xx(b,a)=phi10_xx(a,b)
          end do
        end do

        ! give values to the field strength, using sin(phi2)=sin(phi3)=sin(phi4)=1 w.l.o.g 
        !(considering phi2,phi3,phi4-independent case, so phi2=phi3=phi4=pi/2 slice will do)
        f0_l(1)   =          fb_t0
        f0_l(2)   =          fb_x0
        f0_l(3)   =f0_y_ads0+fb_y0
        f0_ll(1,2)=0.0d0 ! NOTE: need to update as (7) in small_bh_AdS5xS5.tex
        f0_ll(1,3)=0.0d0 ! NOTE: need to update as (7) in small_bh_AdS5xS5.tex
        f0_ll(2,3)=0.0d0 ! NOTE: need to update as (7) in small_bh_AdS5xS5.tex

        ! give values to the ads field strength, using sin(phi2)=sin(phi3)=sin(phi4)=1 w.l.o.g 
        !(considering phi2,phi3,phi4-independent case, so phi2=phi3=phi4=pi/2 slice will do) 
        fads_ll(1,2)=f0_tx_ads0
        fads_l(3)   =f0_y_ads0

        ! calculate field strength F_abcde = 
        do a=1,5
          do b=1,5
            do c=1,5
              do d=1,5
                do e=1,5
                  f_lllll(a,b,c,d,e)=0.0d0 !need to update as (5) in small_bh_AdS5xS5.tex
                end do
              end do
            end do
          end do
        end do
        
        ! calculate raised field strength F_a^bcde = F_a
        do a=1,5
          do b=1,5
            do c=1,5
              do d=1,5
                do e=1,5
                  f_luuuu(a,b,c,d,e)=0.0d0
                  do m=1,5
                    do n=1,5
                      do o=1,5
                        do p=1,5
                          f_luuuu(a,b,c,d,e)=f_luuuu(a,b,c,d,e)
     &                                      +f_lllll(a,m,n,o,p)
     &                                      *g0_uu(b,m)*g0_uu(c,n)
     &                                      *g0_uu(d,o)*g0_uu(e,p)
                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do

        ! calculate contraction of field strength F_acdem F_b^cdem
        do a=1,5
          do b=1,5
            ff_ll(a,b)=0.0d0
            do c=1,5
              do d=1,5
                do e=1,5
                  do m=1,5
                    ff_ll(a,b)=ff_ll(a,b)
     &                        +f_lllll(a,c,d,e,m)
     ^                        *f_luuuu(b,c,d,e,m)
                  end do
                end do
              end do
            end do
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
        do a=1,5
          do b=1,5
            do c=1,5
              do e=1,5
                gamma_ull_x(a,b,c,e)=0
                do d=1,5
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
        gamma_ull_x(4,4,3,3)=-1 
        gamma_ull_x(5,5,3,3)=-1 
        gamma_ull_x(2,4,4,3)=0
        gamma_ull_x(3,4,4,3)=1
        gamma_ull_x(4,3,4,3)=-1
        gamma_ull_x(5,5,4,4)=-1
        gamma_ull_x(2,5,5,3)=0
        gamma_ull_x(2,5,5,4)=0
        gamma_ull_x(3,5,5,3)=1
        gamma_ull_x(3,5,5,4)=0
        gamma_ull_x(4,5,5,4)=1
        gamma_ull_x(5,4,5,4)=-1
        gamma_ull_x(5,3,5,3)=-1

        ! calculate Riemann tensor at point i
        !(R^a_bcd =gamma^a_bd,c - gamma^a_bc,d
        !          +gamma^a_ce gamma^e_bd - gamma^a_de gamma^e_bc)
        do a=1,5
          do b=1,5
            do c=1,5
              do d=1,5
                riemann_ulll(a,b,c,d)=
     &                gamma_ull_x(a,b,d,c)-gamma_ull_x(a,b,c,d)
                do e=1,5
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
        do b=1,5
          do d=1,5
            ricci_ll(b,d)=0
            do a=1,5
              ricci_ll(b,d)=ricci_ll(b,d)+riemann_ulll(a,b,a,d)
            end do
          end do
        end do

        ! calculate raised Ricci tensor at point i
        !(R_a^b = R_ad g^db)
        do a=1,5
          do b=1,5
            ricci_lu(a,b)=0
            do d=1,5
              ricci_lu(a,b)=ricci_lu(a,b)+ricci_ll(a,d)*g0_uu(d,b)
            end do
          end do
        end do

        ! calculate Ricci scalar
        !(R = R_a^a)
        ricci=0
        do a=1,5
          ricci=ricci+ricci_lu(a,a)
        end do
  
        ! calculates Einstein tensor at point i
        !(G_ab = R_ab - 1/2 R g_ab)
        do a=1,5
          do b=1,5
            einstein_ll(a,b)=ricci_ll(a,b)-0.5d0*ricci*g0_ll(a,b)
          end do
        end do

        ! calculates stress-energy tensor at point i 
        !(T_ab = 2*phi1,a phi1,b - (phi1,c phi1,d) g^cd g_ab + ...)
        grad_phi1_sq=0
        do a=1,5
          do b=1,5
            grad_phi1_sq=grad_phi1_sq
     &                  +phi10_x(a)*phi10_x(b)*g0_uu(a,b)
          end do
        end do

        do a=1,5
          do b=1,5
            set_ll(a,b)=
     &            phi10_x(a)*phi10_x(b)
     &           -g0_ll(a,b)*(grad_phi1_sq/2)
          end do
        end do

        return
        end
