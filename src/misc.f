c----------------------------------------------------------------------
c miscellaneous numerical routines for AdS5xS5
c----------------------------------------------------------------------

c-----------------------------------------------------------------------
c specific x first derivative routine used by excision routines 
c below, for method 0,3,4 only
c
c-----------------------------------------------------------------------
        subroutine df1_int_x(f,f_x,dx,i,j,chr,ex,Nx,Ny)

        implicit none
        integer Nx,Ny,i,j
        real*8 f(Nx,Ny),chr(Nx,Ny),ex,f_x,dx
        logical first
        save first
        data first/.true./

        if (i.eq.1.or.chr(i-1,j).eq.ex) then
           if (i.le.(Nx-1).and.chr(i+1,j).ne.ex) then
              f_x=(-f(i,j)+f(i+1,j))/dx
!              write(*,*) 'df1_int_x: warning ... i=1 first order'
!              write(*,*) '    i,Nx,dx=',i,Nx,dx
           else
              if (first) then
                 first=.false.
                 write(*,*) 'df1_int_x: error in chr stencil (A)'
                 write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
                 write(*,*) '    (first error only)'
              end if
              f_x=0
              return
           end if
        else if (i.eq.Nx.or.chr(i+1,j).eq.ex) then
           if (i.ge.2.and.chr(i-1,j).ne.ex) then
              f_x=(f(i,j)-f(i-1,j))/dx
!              write(*,*) 'df1_int_x: warning ... i=Nx first order'
!              write(*,*) '    i,Nx,dx=',i,Nx,dx
           else
              if (first) then
                 first=.false.
                 write(*,*) 'df1_int_x: error in chr stencil (B)'
                 write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
                 write(*,*) '    (first error only)'
              end if
              f_x=0
              return
           end if
        else
           if (chr(i+1,j).ne.ex.and.chr(i-1,j).ne.ex) then
              f_x=(f(i+1,j)-f(i-1,j))/2/dx
           else
              if (first) then
                 first=.false.
                 write(*,*) 'df1_int_x: error in chr stencil (C)'
                 write(*,*) '    i,j,Nx,Ny,dx=',i,j,Nx,Ny,dx
                 write(*,*) '    (first error only)'
              end if
              f_x=0
              return
           end if
        end if

        return
        end
 
!----------------------------------------------------------------------
        subroutine df1_int_y(f,f_y,dy,i,j,chr,ex,Nx,Ny)

        implicit none
        integer Nx,Ny,i,j
        real*8 f(Nx,Ny),chr(Nx,Ny),ex,f_y,dy
        logical first
        save first
        data first/.true./

        !--------------------------------------------------------------

        if ((j.eq.1).or.(chr(i,j-1).eq.ex)) then
           if (j.le.(Ny-1).and.chr(i,j+1).ne.ex) then
              f_y=(-f(i,j)+f(i,j+1))/dy
!              write(*,*) 'df1_int_y: warning ... j=1 first order'
!              write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
           else
              if (first) then
                 first=.false.
                 write(*,*) 'df1_int_y: error in chr stencil (D)'
                 write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
                 write(*,*) '    (first error only)'
              end if
              f_y=0
              return
           end if
        else if ((j.eq.Ny).or.(chr(i,j+1).eq.ex)) then
           if (j.ge.2.and.chr(i,j-1).ne.ex) then
              f_y=(f(i,j)-f(i,j-1))/dy
!              write(*,*) 'df1_int_y: warning ... j=Ny first order'
!              write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
           else
              if (first) then
                 first=.false.
                 write(*,*) 'df1_int_y: error in chr stencil (E)'
                 write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
                 write(*,*) '    (first error only)'
              end if
              f_y=0
              return
           end if
        else
           if ((chr(i,j+1).ne.ex.and.chr(i,j-1).ne.ex)) then
              f_y=(f(i,j+1)-f(i,j-1))/2/dy
           else
              if (first) then
                 first=.false.
                 write(*,*) 'df1_int_y: error in chr stencil (F)'
                 write(*,*) '    i,j,Nx,Ny,dy=',i,j,Nx,Ny,dy
                 write(*,*) '    (first error only)'
              end if
              f_y=0
              return
           end if
        end if

        f_y=0 !TEST!

        return
        end

c----------------------------------------------------------------------
c the following computes all first derivatives of f,
c at a point i at time level n.
c
c stencil reduces to first order on excision surface
c----------------------------------------------------------------------
        subroutine df1_int(f_np1,f_n,f_nm1,f_t,f_x,f_y,
     &                     dx,dy,dt,i,j,chr,ex,Nx,Ny,name)
        implicit none
        integer Nx,Ny,i,j
        real*8 f_np1(Nx,Ny),f_n(Nx,Ny),f_nm1(Nx,Ny)
        real*8 f_t,f_x,f_y,dx,dy,dt,ex,chr(Nx,Ny)
        character*(*) name

        logical ltrace,first
        parameter (ltrace=.false.) 
        save first
        data first/.true./
        real*8 f_x_np1,f_x_nm1

        if (chr(i,j).eq.ex) then
         write(*,*) 'df1_int: error ... point excised'
         stop
        end if

        f_t=(f_np1(i,j)-f_nm1(i,j))/2/dt

        call df1_int_x(f_n,f_x,dx,i,j,chr,ex,Nx,Ny)
        call df1_int_y(f_n,f_y,dy,i,j,chr,ex,Nx,Ny)

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
     &                     f_t,f_x,f_y,
     &                     f_tt,f_tx,f_ty,f_xx,f_xy,f_yy,
     &                     dx,dy,dt,i,j,chr,ex,Nx,Ny,name)
        implicit none
        integer Nx,Ny,i,j
        real*8 f_np1(Nx,Ny),f_n(Nx,Ny),f_nm1(Nx,Ny)
        real*8 f_t,f_x,f_y,f_tt,f_tx,f_ty,f_xx,f_xy,f_yy
        real*8 dx,dy,dt,ex,chr(Nx,Ny)
        character*(*) name
        logical first
        save first
        data first/.true./

        logical ltrace
        parameter (ltrace=.false.)
        real*8 f_x_np1,f_x_nm1,f_y_np1,f_y_nm1
        real*8 f_y_ip1,f_y_ip2
        real*8 f_y_im1,f_y_im2

        call df1_int(f_np1,f_n,f_nm1,f_t,f_x,f_y,
     &               dx,dy,dt,i,j,chr,ex,Nx,Ny,name)

        f_tt=(f_np1(i,j)-2*f_n(i,j)+f_nm1(i,j))/dt/dt 

        f_xx=0
        f_xy=0
        f_yy=0

        if (chr(i,j).eq.ex) then
          write(*,*) 'df2_int: error ... point excised'
          stop
        end if

        call df1_int_x(f_np1,f_x_np1,dx,i,j,chr,ex,Nx,Ny)
        call df1_int_x(f_nm1,f_x_nm1,dx,i,j,chr,ex,Nx,Ny)

        f_tx=(f_x_np1-f_x_nm1)/2/dt

        call df1_int_y(f_np1,f_y_np1,dy,i,j,chr,ex,Nx,Ny)
        call df1_int_y(f_nm1,f_y_nm1,dy,i,j,chr,ex,Nx,Ny)

        f_ty=(f_y_np1-f_y_nm1)/2/dt

        !i

        if (i.eq.1.or.(chr(i-1,j).eq.ex)) then
          if (i.ge.(Nx-1).or.chr(i+1,j).eq.ex.or.chr(i+2,j).eq.ex) then
            if (first) then
              first=.false.
              write(*,*) 'df2_int: error in chr (A)'
              write(*,*) '    i,j,Nx,Ny,dx,dy=',i,j,Nx,Ny,dx,dy
              write(*,*) '    (first error only)'
            end if
            return
          end if

          f_xx=(f_n(i+2,j)-2*f_n(i+1,j)+f_n(i,j))/dx/dx 

          call df1_int_y(f_n,f_y_ip1,dy,i+1,j,chr,ex,Nx,Ny)
          call df1_int_y(f_n,f_y_ip2,dy,i+2,j,chr,ex,Nx,Ny)
          f_xy=(-3*f_y+4*f_y_ip1-f_y_ip2)/2/dx

        else if (i.eq.Nx.or.(chr(i+1,j).eq.ex)) then
          if (i.le.2.or.
     &        chr(i-1,j).eq.ex.or.chr(i-2,j).eq.ex) then
            if (first) then
              first=.false.
              write(*,*) 'df2_int: error in chr (B)'
              write(*,*) '    i,j,Nx,Ny,dx,dy=',i,j,Nx,Ny,dx,dy
              write(*,*) '    (first error only)'
            end if
            return
          end if

          f_xx=(f_n(i,j)-2*f_n(i-1,j)+f_n(i-2,j))/dx/dx 

          call df1_int_y(f_n,f_y_im1,dy,i-1,j,chr,ex,Nx,Ny)
          call df1_int_y(f_n,f_y_im2,dy,i-2,j,chr,ex,Nx,Ny)
          f_xy=(3*f_y-4*f_y_ip1+f_y_ip2)/2/dx

        else if (chr(i+1,j).ne.ex.and.chr(i-1,j).ne.ex) then

          f_xx=(f_n(i+1,j)-2*f_n(i,j)+f_n(i-1,j))/dx/dx 

        else
          if (first) then
            first=.false.
            write(*,*) 'df2_int: error in chr (C)'
            write(*,*) '    i,j,Nx,Ny,dx,dy=',i,j,Nx,Ny,dx,dy
            write(*,*) '    (first error only)'
          end if
          return
        end if

        !j

        if (j.eq.1.or.(chr(i,j-1).eq.ex)) then
          if (j.ge.(Ny-1).or.chr(i,j+1).eq.ex.or.chr(i,j+2).eq.ex) then
            if (first) then
              first=.false.
              write(*,*) 'df2_int: error in chr (D)'
              write(*,*) '    i,j,Nx,Ny,dx,dy=',i,j,Nx,Ny,dx,dy
              write(*,*) '    (first error only)'
            end if
            return
          end if

          f_yy=(f_n(i,j+2)-2*f_n(i,j+1)+f_n(i,j))/dy/dy 

        else if (j.eq.Ny.or.(chr(i,j+1).eq.ex)) then
          if (j.le.2.or.
     &        chr(i,j-1).eq.ex.or.chr(i,j-2).eq.ex) then
            if (first) then
              first=.false.
              write(*,*) 'df2_int: error in chr (E)'
              write(*,*) '    i,j,Nx,Ny,dx,dy=',i,j,Nx,Ny,dx,dy
              write(*,*) '    (first error only)'
            end if
            return
          end if

          f_yy=(f_n(i,j)-2*f_n(i,j-1)+f_n(i,j-2))/dy/dy 

        else if (chr(i,j+1).ne.ex.and.chr(i,j-1).ne.ex) then

          f_yy=(f_n(i,j+1)-2*f_n(i,j)+f_n(i,j-1))/dy/dy 

        else
          if (first) then
            first=.false.
            write(*,*) 'df2_int: error in chr (F)'
            write(*,*) '    i,j,Nx,Ny,dx,dy=',i,j,Nx,Ny,dx,dy
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

        f_ty=0 !TEST!
        f_xy=0 !TEST!
        f_yy=0 !TEST!

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
        subroutine gauss2d(f,amp,r0,delta,xu0,yu0,ax,ay,
     &                     amp2,r02,delta2,xu02,yu02,ax2,ay2,
     &                     L,x,y,Nx,Ny)
        implicit none
        integer i,j,Nx,Ny
        real*8 f(Nx,Ny),x(Nx),y(Ny),L
        real*8 amp,r0,delta,ax,ay,xu0,yu0,r
        real*8 amp2,r02,delta2,ax2,ay2,xu02,yu02,r2
        real*8 x0,y0

        real*8 PI
        parameter (PI=3.141592653589793d0)

        !--------------------------------------------------------------

        do i=1,Nx
          do j=1,Ny
            f(i,j)=0
            x0=x(i)
            y0=y(j)
            
            r=sqrt((x0-xu0)**2/ax**2)
            r2=sqrt((x0-xu02)**2/ax2**2)

            ! Gaussian phi=amp*exp(-r^2)/r^4 profile 
            ! remember that phi=phi1*(1-x^2)^3
            if (r.ge.1) then
               f(i,j)=0
            else if (r.gt.r0) then
               f(i,j)=amp*exp(-((r-r0)/delta)**2)*(1-x0**2)
            else
               f(i,j)=amp*(1-x0**2)
            end if

          end do
        end do

        return
        end

c-----------------------------------------------------------------------
c for variables with lin_zero_bnd ... zeros residual there
c-----------------------------------------------------------------------
        subroutine lin_zero_bnd_res(f,phys_bdy,all,Nx,Ny)
        implicit none
        integer Nx,Ny,all
        real*8 f(Nx,Ny)
        integer phys_bdy(4)

        integer i,is,ie
        integer j,js,je

        ! initialize fixed-size variables
        data i,is,ie/0,0,0/
        data j,js,je/0,0,0/

        !--------------------------------------------------------------

        if (phys_bdy(1).eq.1.or.all.eq.1) then
          do j=2,Ny-1
            f(2,j)=0
          end do
        end if

        if (phys_bdy(2).eq.1.or.all.eq.1) then
          do j=2,Ny-1
            f(Nx-1,j)=0
          end do
        end if

        if (phys_bdy(3).eq.1.or.all.eq.1) then
          do i=2,Nx-1
            f(i,2)=0
          end do
        end if

        if (phys_bdy(4).eq.1.or.all.eq.1) then
          do i=2,Nx-1
            f(i,Ny-1)=0
          end do
        end if

        return
        end

c----------------------------------------------------------------------
c initializes the metric to an exact black hole solution
c with radius parameter r0
c----------------------------------------------------------------------
        subroutine init_schw(gb_tt,gb_tx,gb_ty,gb_xx,gb_xy,gb_yy,
     &                       psi,omega,r0,L,phys_bdy,chr,ex,x,y,Nx,Ny)
        implicit none
        real*8 gb_tt(Nx,Ny),gb_tx(Nx,Ny),gb_ty(Nx,Ny),gb_xx(Nx,Ny)
        real*8 gb_xy(Nx,Ny),gb_yy(Nx,Ny),psi(Nx,Ny),omega(Nx,Ny)
        real*8 chr(Nx,Ny),x(Nx),y(Ny)
        real*8 r0,ex,L
        integer phys_bdy(4)
        integer i,j,Nx,Ny

        real*8 r_h,x_h,x0

        !--------------------------------------------------------------

        ! compute horizon global radius r_h and corresponding compactified x_h
        r_h=L*sqrt(2*sqrt(1+4*(r0/L)**2)-2)/2
        x_h=r_h/(1+r_h)

        ! initialize metric 
        do i=1,Nx
          do j=1,Ny
           if (chr(i,j).eq.ex) then
              gb_tt(i,j)=0
              gb_tx(i,j)=0
              gb_xx(i,j)=0
              gb_yy(i,j)=0
              psi(i,j)=0
              omega(i,j)=0
           else
              x0=x(i)

              ! EF-like-near-horizon Schwarzschild-like-near-bdy coordinates
              gb_tt(i,j)=(r0/x0)**2/(1+x0)*(1-x0)
              gb_tx(i,j)=1/(1+x0)/(1-x_h)**4*(1-x0)
              gb_xx(i,j)=-((L**2*(-1 + x0)*((-1 + x0)**4*x0**4 + 
     &        L**2*((-1 + x0)**6*x0**2 - r0**2*(-1 + x_h)**8)))/
     &        ((1 + x0)*(L**2*(-1 + x0)**2 + x0**2)*
     &        (-x0**4 + L**2*(-1 + x0)**2*(r0**2*(-1 + x0)**2 - x0**2))*
     &        (-1 + x_h)**8)) 
              gb_yy(i,j)=0
              psi(i,j)  =0
              omega(i,j)=0

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
        end do

        ! (REGION) x=0; impose regularity conditions 
        call axi_reg_g(gb_tt,gb_tx,gb_ty,gb_xx,gb_xy,gb_yy,psi,omega,
     &                 chr,ex,L,x,y,Nx,Ny)

        return
        end

c----------------------------------------------------------------------
c calculates all the tensorial objects in x coordinates, at point i
c----------------------------------------------------------------------
        subroutine tensor_init(
     &                  gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                  gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                  gb_ty_np1,gb_ty_n,gb_ty_nm1,
     &                  gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                  gb_xy_np1,gb_xy_n,gb_xy_nm1,
     &                  gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &                  psi_np1,psi_n,psi_nm1,
     &                  omega_np1,omega_n,omega_nm1,
     &                  fb_t_np1,fb_t_n,fb_t_nm1,
     &                  fb_x_np1,fb_x_n,fb_x_nm1,
     &                  fb_y_np1,fb_y_n,fb_y_nm1,
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
     &                  A_l,A_l_x,Hads_l,Hads_l_x,
     &                  gamma_ull,gamma_ull_x,
     &                  riemann_ulll,ricci_ll,ricci_lu,ricci,
     &                  f1_l,f1_l_x,f2_ll,f2_ll_x,
     &                  f1ads_l,f1ads_l_x,f2ads_ll,f2ads_ll_x,
     &                  h1_l,h1_l_x,h2_ll,h2_ll_x,
     &                  sA,sB,tA,tB,
     &                  phi10_x,phi10_xx,
     &                  gammagg,gammahh,gammagh,gammahg,dlll,cuuuu,
     &                  x,y,dt,chr,L,ex,Nx,Ny,i,j)
        implicit none

        integer Nx,Ny
        integer i,j

        real*8 chr(Nx,Ny),ex
        real*8 x(Nx),y(Ny),dt,L

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
        real*8 fb_t_np1(Nx,Ny),fb_t_n(Nx,Ny),fb_t_nm1(Nx,Ny)
        real*8 fb_x_np1(Nx,Ny),fb_x_n(Nx,Ny),fb_x_nm1(Nx,Ny)
        real*8 fb_y_np1(Nx,Ny),fb_y_n(Nx,Ny),fb_y_nm1(Nx,Ny)

        integer a,b,c,d,e,f,m,n,o,p

        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 dx,dy
        real*8 x0,y0   

        real*8 grad_phi1_sq

        real*8 dimA,dimB

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
        real*8 levicivi3(3,3,3)
        real*8 vol(3,3,3),vol_x(3,3,3,3)
        real*8 volads(3,3,3),volads_x(3,3,3,3)
        real*8 volh(3,3,3),volh_x(3,3,3,3)
        real*8 sqrtdetg,sqrtdetg_x(3)
        real*8 sqrtdetgads,sqrtdetgads_x(3)
        real*8 sqrtdeth,sqrtdeth_x(3)
        real*8 riccibar_ll(3,3),riccibar_lu(3,3),riccibar

        real*8 gammagg(3,3,3),gammahh(3,3,3)
        real*8 gammagh(3,3,3),gammahg(3,3,3)
        real*8 cuuuu(3,3,3,3),dlll(3,3,3)

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
        real*8 gb_tt_t,gb_tt_x,gb_tt_y
        real*8 gb_tt_tt,gb_tt_tx,gb_tt_ty,gb_tt_xx,gb_tt_xy,gb_tt_yy
        real*8 gb_tx_t,gb_tx_x,gb_tx_y
        real*8 gb_tx_tt,gb_tx_tx,gb_tx_ty,gb_tx_xx,gb_tx_xy,gb_tx_yy
        real*8 gb_ty_t,gb_ty_x,gb_ty_y
        real*8 gb_ty_tt,gb_ty_tx,gb_ty_ty,gb_ty_xx,gb_ty_xy,gb_ty_yy
        real*8 gb_xx_t,gb_xx_x,gb_xx_y
        real*8 gb_xx_tt,gb_xx_tx,gb_xx_ty,gb_xx_xx,gb_xx_xy,gb_xx_yy
        real*8 gb_xy_t,gb_xy_x,gb_xy_y
        real*8 gb_xy_tt,gb_xy_tx,gb_xy_ty,gb_xy_xx,gb_xy_xy,gb_xy_yy
        real*8 gb_yy_t,gb_yy_x,gb_yy_y
        real*8 gb_yy_tt,gb_yy_tx,gb_yy_ty,gb_yy_xx,gb_yy_xy,gb_yy_yy
        real*8 psi_t,psi_x,psi_y
        real*8 psi_tt,psi_tx,psi_ty,psi_xx,psi_xy,psi_yy
        real*8 omega_t,omega_x,omega_y
        real*8 omega_tt,omega_tx,omega_ty,omega_xx,omega_xy,omega_yy
        real*8 phi1_t,phi1_x,phi1_y
        real*8 phi1_tt,phi1_tx,phi1_ty,phi1_xx,phi1_xy,phi1_yy
  
        real*8 gb_tt0,gb_tx0,gb_ty0,gb_xx0,gb_xy0,gb_yy0,phi10

        real*8 g0_tt_ads0,g0_xx_ads0,g0_yy_ads0
        real*8 g0_tt_ads_x,g0_tt_ads_xx
        real*8 g0_xx_ads_x,g0_xx_ads_xx
        real*8 g0_yy_ads_x,g0_yy_ads_xx

        real*8 psi0,omega0

        real*8 gA_ads0,gB_ads0
        real*8 gA_ads_x,gA_ads_xx        
        real*8 gB_ads_y,gB_ads_yy

        real*8 fb_t_t,fb_t_x,fb_t_y 
        real*8 fb_x_t,fb_x_x,fb_x_y
        real*8 fb_y_t,fb_y_x,fb_y_y

        real*8 Hb_t_t,Hb_t_x,Hb_t_y 
        real*8 Hb_x_t,Hb_x_x,Hb_x_y
        real*8 Hb_y_t,Hb_y_x,Hb_y_y

        real*8 Hb_t0,Hb_x0,Hb_y0

        real*8 fb_t0,fb_x0,fb_y0

        real*8 f1_y_ads0
!----------------------------------------------------------------------
        
        dx=(x(2)-x(1))
        dy=(y(2)-y(1))

        x0=x(i)
        y0=y(j)  

        ! set dimensions of S3 and S4 subspaces
        dimA=3d0 
        dimB=4d0

        ! set gads values using sin(theta1)=sin(theta2)=1 w.l.o.g 
        !(considering theta1,theta2-independent case, so theta1=theta2=pi/2 slice will do)
        g0_tt_ads0 =-((1-x0)**2+x0**2)/(1-x0)**2
        g0_xx_ads0 =1/((1-x0)**2+x0**2)/(1-x0)**2
        g0_yy_ads0=PI**2

        ! set gAads values
        gA_ads0=x0**2/(1-x0)**2

        ! set gBads values
        gB_ads0=L**2*sin(PI*y0/L)**2 !y-dependence

        ! set f1ads values using sin(phi2)=sin(phi3)=sin(phi4)=1 w.l.o.g 
        !(considering phi2,phi3,phi4-independent case, so phi2=phi3=phi4=pi/2 slice will do)
        f1_y_ads0  = 4/L*PI

        ! set gbar values
        gb_tt0=gb_tt_n(i,j)
        gb_tx0=gb_tx_n(i,j)
        gb_ty0=gb_ty_n(i,j)
        gb_xx0=gb_xx_n(i,j)
        gb_xy0=gb_xy_n(i,j)
        gb_yy0=gb_yy_n(i,j)
        psi0  =psi_n(i,j) 
        omega0=omega_n(i,j)

        ! set fbar values
        fb_t0=fb_t_n(i,j)
        fb_x0=fb_x_n(i,j)
        fb_y0=fb_y_n(i,j)

        ! set hbar values
        Hb_t0=Hb_t_n(i,j)
        Hb_x0=Hb_x_n(i,j)
        Hb_y0=Hb_y_n(i,j)

        ! set phi1 value
        phi10=phi1_n(i,j)

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
        gB_ads_y =L*PI*sin(2*PI*y0/L)
        gB_ads_yy=2*PI**2*cos(2*PI*y0/L) !y-dependence

        ! calculate gbar derivatives
        call df2_int(gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &       gb_tt_t,gb_tt_x,gb_tt_y,
     &       gb_tt_tt,gb_tt_tx,gb_tt_ty,
     &       gb_tt_xx,gb_tt_xy,gb_tt_yy,
     &       dx,dy,dt,i,j,chr,ex,Nx,Ny,'gb_tt')
        call df2_int(gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &       gb_tx_t,gb_tx_x,gb_tx_y,
     &       gb_tx_tt,gb_tx_tx,gb_tx_ty,
     &       gb_tx_xx,gb_tx_xy,gb_tx_yy,
     &       dx,dy,dt,i,j,chr,ex,Nx,Ny,'gb_tx')
        call df2_int(gb_ty_np1,gb_ty_n,gb_ty_nm1,
     &       gb_ty_t,gb_ty_x,gb_ty_y,
     &       gb_ty_tt,gb_ty_tx,gb_ty_ty,
     &       gb_ty_xx,gb_ty_xy,gb_ty_yy,
     &       dx,dy,dt,i,j,chr,ex,Nx,Ny,'gb_ty')
        call df2_int(gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &       gb_xx_t,gb_xx_x,gb_xx_y,
     &       gb_xx_tt,gb_xx_tx,gb_xx_ty,
     &       gb_xx_xx,gb_xx_xy,gb_xx_yy,
     &       dx,dy,dt,i,j,chr,ex,Nx,Ny,'gb_xx')
        call df2_int(gb_xy_np1,gb_xy_n,gb_xy_nm1,
     &       gb_xy_t,gb_xy_x,gb_xy_y,
     &       gb_xy_tt,gb_xy_tx,gb_xy_ty,
     &       gb_xy_xx,gb_xy_xy,gb_xy_yy,
     &       dx,dy,dt,i,j,chr,ex,Nx,Ny,'gb_xy')
        call df2_int(gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &       gb_yy_t,gb_yy_x,gb_yy_y,
     &       gb_yy_tt,gb_yy_tx,gb_yy_ty,
     &       gb_yy_xx,gb_yy_xy,gb_yy_yy,
     &       dx,dy,dt,i,j,chr,ex,Nx,Ny,'gb_yy')
        call df2_int(psi_np1,psi_n,psi_nm1,
     &       psi_t,psi_x,psi_y,
     &       psi_tt,psi_tx,psi_ty,
     &       psi_xx,psi_xy,psi_yy,
     &       dx,dy,dt,i,j,chr,ex,Nx,Ny,'psi')
        call df2_int(omega_np1,omega_n,omega_nm1,
     &       omega_t,omega_x,omega_y,
     &       omega_tt,omega_tx,omega_ty,
     &       omega_xx,omega_xy,omega_yy,
     &       dx,dy,dt,i,j,chr,ex,Nx,Ny,'omega')

        ! calculate fbar derivatives
        call df1_int(fb_t_np1,fb_t_n,fb_t_nm1,
     &       fb_t_t,fb_t_x,fb_t_y,
     &       dx,dy,dt,i,j,chr,ex,Nx,Ny,'fb_t')
        call df1_int(fb_x_np1,fb_x_n,fb_x_nm1,
     &       fb_x_t,fb_x_x,fb_x_y,
     &       dx,dy,dt,i,j,chr,ex,Nx,Ny,'fb_x')
        call df1_int(fb_y_np1,fb_y_n,fb_y_nm1,
     &       fb_y_t,fb_y_x,fb_y_y,
     &       dx,dy,dt,i,j,chr,ex,Nx,Ny,'fb_y')

        ! calculate hbar derivatives
        call df1_int(Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &       Hb_t_t,Hb_t_x,Hb_t_y,
     &       dx,dy,dt,i,j,chr,ex,Nx,Ny,'Hb_t')
        call df1_int(Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &       Hb_x_t,Hb_x_x,Hb_x_y,
     &       dx,dy,dt,i,j,chr,ex,Nx,Ny,'Hb_x')
        call df1_int(Hb_y_np1,Hb_y_n,Hb_y_nm1,
     &       Hb_y_t,Hb_y_x,Hb_y_y,
     &       dx,dy,dt,i,j,chr,ex,Nx,Ny,'Hb_y')

        ! calculate phi1 derivatives
        call df2_int(phi1_np1,phi1_n,phi1_nm1,
     &       phi1_t,phi1_x,phi1_y,
     &       phi1_tt,phi1_tx,phi1_ty,
     &       phi1_xx,phi1_xy,phi1_yy,
     &       dx,dy,dt,i,j,chr,ex,Nx,Ny,'phi1')

        ! give values to the metric, using sin(theta1)=sin(theta2)=1 w.l.o.g 
        !(considering theta1,theta2-independent case, so theta1=theta2=pi/2 slice will do)
        g0_ll(1,1)=g0_tt_ads0+gb_tt0*(1-x0**2)
        g0_ll(1,2)=           gb_tx0*(1-x0**2)
        g0_ll(1,3)=           gb_ty0*(1-x0**2)**2  
        g0_ll(2,2)=g0_xx_ads0+gb_xx0*(1-x0**2)
        g0_ll(2,3)=           gb_xy0*(1-x0**2)**2  
        g0_ll(3,3)=g0_yy_ads0+gb_yy0*(1-x0**2)**3

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
        g0_ll_x(1,1,3)   =0
     &                   +gb_tt_y*(1-x0**2)
        g0_ll_xx(1,1,1,1)=0
     &                   +gb_tt_tt*(1-x0**2)
        g0_ll_xx(1,1,1,2)=0
     &                   +gb_tt_tx*(1-x0**2)
     &                   +gb_tt_t*(-2*x0)
        g0_ll_xx(1,1,1,3)=0
     &                   +gb_tt_ty*(1-x0**2)
        g0_ll_xx(1,1,2,2)=g0_tt_ads_xx
     &                   +gb_tt_xx*(1-x0**2)
     &                   +gb_tt_x*(-2*x0)
     &                   +gb_tt_x*(-2*x0)
     &                   +gb_tt0*(-2)
        g0_ll_xx(1,1,2,3)=0
     &                   +gb_tt_xy*(1-x0**2)
     &                   +gb_tt_y*(-2*x0)
        g0_ll_xx(1,1,3,3)=0
     &                   +gb_tt_yy*(1-x0**2)

        g0_ll_x(1,2,1)   =0
     &                   +gb_tx_t*(1-x0**2)
        g0_ll_x(1,2,2)   =0
     &                   +gb_tx_x*(1-x0**2)
     &                   +gb_tx0*(-2*x0)
        g0_ll_x(1,2,3)   =0
     &                   +gb_ty_y*(1-x0**2)
        g0_ll_xx(1,2,1,1)=0
     &                   +gb_tx_tt*(1-x0**2)
        g0_ll_xx(1,2,1,2)=0
     &                   +gb_tx_tx*(1-x0**2)
     &                   +gb_tx_t*(-2*x0)
        g0_ll_xx(1,2,1,3)=0
     &                   +gb_tx_ty*(1-x0**2)
        g0_ll_xx(1,2,2,2)=0
     &                   +gb_tx_xx*(1-x0**2)
     &                   +gb_tx_x*(-2*x0)
     &                   +gb_tx_x*(-2*x0)
     &                   +gb_tx0*(-2)
        g0_ll_xx(1,2,2,3)=0
     &                   +gb_tx_xy*(1-x0**2)
     &                   +gb_tx_y*(-2*x0)
        g0_ll_xx(1,2,3,3)=0
     &                   +gb_tx_yy*(1-x0**2)

        g0_ll_x(1,3,1)   =0                            
     &                   +gb_ty_t*(1-x0**2)**2
        g0_ll_x(1,3,2)   =0
     &                   +gb_ty_x*(1-x0**2)**2
     &                   +gb_ty0*2*(1-x0**2)*(-2*x0)
        g0_ll_x(1,3,3)   =0                            
     &                   +gb_ty_y*(1-x0**2)**2
        g0_ll_xx(1,3,1,1)=0
     &                   +gb_ty_tt*(1-x0**2)**2
        g0_ll_xx(1,3,1,2)=0
     &                   +gb_ty_tx*(1-x0**2)**2
     &                   +gb_ty_t*2*(1-x0**2)*(-2*x0)
        g0_ll_xx(1,3,1,3)=0
     &                   +gb_ty_ty*(1-x0**2)**2
        g0_ll_xx(1,3,2,2)=0
     &                   +gb_ty_xx*(1-x0**2)**2
     &                   +gb_ty_x*2*(1-x0**2)*(-2*x0)
     &                   +gb_ty_x*2*(1-x0**2)*(-2*x0)
     &                   +gb_ty0*2*(-2+6*x0**2)
        g0_ll_xx(1,3,2,3)=0
     &                   +gb_ty_xy*(1-x0**2)**2
     &                   +gb_ty_y*2*(1-x0**2)*(-2*x0)
        g0_ll_xx(1,3,3,3)=0                            
     &                   +gb_ty_yy*(1-x0**2)**2

        g0_ll_x(2,2,1)   =0
     &                   +gb_xx_t*(1-x0**2)
        g0_ll_x(2,2,2)   =g0_xx_ads_x
     &                   +gb_xx_x*(1-x0**2)
     &                   +gb_xx0*(-2*x0)
        g0_ll_x(2,2,3)   =0
     &                   +gb_xx_y*(1-x0**2)
        g0_ll_xx(2,2,1,1)=0
     &                   +gb_xx_tt*(1-x0**2)
        g0_ll_xx(2,2,1,2)=0
     &                   +gb_xx_tx*(1-x0**2)
     &                   +gb_xx_t*(-2*x0)
        g0_ll_xx(2,2,1,3)=0
     &                   +gb_xx_ty*(1-x0**2)
        g0_ll_xx(2,2,2,2)=g0_xx_ads_xx
     &                   +gb_xx_xx*(1-x0**2)
     &                   +gb_xx_x*(-2*x0)
     &                   +gb_xx_x*(-2*x0)
     &                   +gb_xx0*(-2)
        g0_ll_xx(2,2,2,3)=0
     &                   +gb_xx_xy*(1-x0**2)
     &                   +gb_xx_y*(-2*x0)
        g0_ll_xx(2,2,3,3)=0
     &                   +gb_xx_yy*(1-x0**2)

        g0_ll_x(2,3,1)   =0                            
     &                   +gb_xy_t*(1-x0**2)**2
        g0_ll_x(2,3,2)   =0
     &                   +gb_xy_x*(1-x0**2)**2
     &                   +gb_xy0*2*(1-x0**2)*(-2*x0)
        g0_ll_x(2,3,3)   =0                            
     &                   +gb_xy_y*(1-x0**2)**2
        g0_ll_xx(2,3,1,1)=0
     &                   +gb_xy_tt*(1-x0**2)**2
        g0_ll_xx(2,3,1,2)=0
     &                   +gb_xy_tx*(1-x0**2)**2
     &                   +gb_xy_t*2*(1-x0**2)*(-2*x0)
        g0_ll_xx(2,3,1,3)=0
     &                   +gb_xy_ty*(1-x0**2)**2
        g0_ll_xx(2,3,2,2)=0
     &                   +gb_xy_xx*(1-x0**2)**2
     &                   +gb_xy_x*2*(1-x0**2)*(-2*x0)
     &                   +gb_xy_x*2*(1-x0**2)*(-2*x0)
     &                   +gb_xy0*2*(-2+6*x0**2)
        g0_ll_xx(2,3,2,3)=0
     &                   +gb_xy_xy*(1-x0**2)**2
     &                   +gb_xy_y*2*(1-x0**2)*(-2*x0)
        g0_ll_xx(2,3,3,3)=0
     &                   +gb_xy_yy*(1-x0**2)**2

        g0_ll_x(3,3,1)   =0
     &                   +gb_yy_t*(1-x0**2)**3
        g0_ll_x(3,3,2)   =g0_yy_ads_x
     &                   +gb_yy_x*(1-x0**2)**3
     &                   +gb_yy0*3*(1-x0**2)**2*(-2*x0)
        g0_ll_x(3,3,3)   =0
     &                   +gb_yy_y*(1-x0**2)**3
        g0_ll_xx(3,3,1,1)=0
     &                   +gb_yy_tt*(1-x0**2)**3
        g0_ll_xx(3,3,1,2)=0
     &                   +gb_yy_tx*(1-x0**2)**3
     &                   +gb_yy_t*3*(1-x0**2)**2*(-2*x0)
        g0_ll_xx(3,3,1,3)=0
     &                   +gb_yy_ty*(1-x0**2)**3
        g0_ll_xx(3,3,2,2)=g0_yy_ads_xx
     &                   +gb_yy_xx*(1-x0**2)**3
     &                   +gb_yy_x*3*(1-x0**2)**2*(-2*x0)
     &                   +gb_yy_x*3*(1-x0**2)**2*(-2*x0)
     &                   +gb_yy0*6*(-1+6*x0**2-5*x0**4)
        g0_ll_xx(3,3,2,3)=0
     &                   +gb_yy_xy*(1-x0**2)**3
     &                   +gb_yy_y*3*(1-x0**2)**2*(-2*x0)
        g0_ll_xx(3,3,3,3)=0
     &                   +gb_yy_yy*(1-x0**2)**3

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

        gA=gA_ads0+psi0*(1-x0**2)*x0**2
        gB=gB_ads0+omega0*(1-x0**2)**3*sin(PI*y0/L)**2 

        gA_x(1)   =0
     &            +psi_t*(1-x0**2)*x0**2
        gA_x(2)   =gA_ads_x
     &            +psi_x*(1-x0**2)*x0**2
     &            +psi0*(2*x0-4*x0**3)
        gA_x(3)   =0
     &            +psi_y*(1-x0**2)*x0**2
        gA_xx(1,1)=0
     &            +psi_tt*(1-x0**2)*x0**2
        gA_xx(1,2)=0
     &            +psi_tx*(1-x0**2)*x0**2
     &            +psi_t*(2*x0-4*x0**3)
        gA_xx(1,3)=0
     &            +psi_ty*(1-x0**2)*x0**2
        gA_xx(2,2)=gA_ads_xx
     &            +psi_xx*(1-x0**2)*x0**2
     &            +psi_x*(2*x0-4*x0**3)
     &            +psi_x*(2*x0-4*x0**3)
     &            +psi0*(2-12*x0**2)
        gA_xx(2,3)=0
     &            +psi_xy*(1-x0**2)*x0**2
     &            +psi_y*(2*x0-4*x0**3)
        gA_xx(3,3)=0
     &            +psi_yy*(1-x0**2)*x0**2

        gB_x(1)   =0
     &            +omega_t*(1-x0**2)**3*sin(PI*y0/L)**2 
        gB_x(2)   =0
     &            +omega_x*(1-x0**2)**3*sin(PI*y0/L)**2
     &            +omega0*3*(1-x0**2)**2*(-2*x0)*sin(PI*y0/L)**2
        gB_x(3)   =gB_ads_y
     &            +omega_y*(1-x0**2)**3*sin(PI*y0/L)**2
     &            +omega0*(1-x0**2)**3*sin(2*PI*y0/L)*PI/L
        gB_xx(1,1)=0
     &            +omega_tt*(1-x0**2)**3*sin(PI*y0/L)**2
        gB_xx(1,2)=0
     &            +omega_tx*(1-x0**2)**3*sin(PI*y0/L)**2
     &            +omega_t*3*(1-x0**2)**2*(-2*x0)*sin(PI*y0/L)**2
        gB_xx(1,3)=0
     &            +omega_ty*(1-x0**2)**3*sin(PI*y0/L)**2
     &            +omega_t*(1-x0**2)**3*sin(2*PI*y0/L)*PI/L
        gB_xx(2,2)=0
     &            +omega_xx*(1-x0**2)**3*sin(PI*y0/L)**2
     &            +omega_x*3*(1-x0**2)**2*(-2*x0)*sin(PI*y0/L)**2
     &            +omega_x*3*(1-x0**2)**2*(-2*x0)*sin(PI*y0/L)**2
     &            +omega0*6*(-1+6*x0**2-5*x0**4)*sin(PI*y0/L)**2
        gB_xx(2,3)=0
     &            +omega_xy*(1-x0**2)**3*sin(PI*y0/L)**2
     &            +omega_x*(1-x0**2)**3*sin(PI*y0/L)*cos(PI*y0/L)*2*PI/L
     &            +omega_y*3*(1-x0**2)**2*(-2*x0)*sin(PI*y0/L)**2
     &            +omega0*3*(1-x0**2)**2*(-2*x0)*sin(2*PI*y0/L)*PI/L
        gB_xx(3,3)=gB_ads_yy 
     &            +omega_yy*(1-x0**2)**3*sin(PI*y0/L)**2
     &            +omega_y*(1-x0**2)**3*sin(2*PI*y0/L)*PI/L
     &            +omega_y*(1-x0**2)**3*sin(2*PI*y0/L)*PI/L 
     &            +omega0*(1-x0**2)**3*cos(2*PI*y0/L)*2*PI**2/L**2

!TEST!
!       if (gB_x(1).ne.0) then
!         write(*,*) 'x,y=',x0,y0
!       end if
!TEST!
!       if (omega_np1(i,j)-omega_nm1(i,j).ne.0) then
!         write(*,*) 'np1,nm1=',omega_np1(i,j),omega_nm1(i,j)
!       end if

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

        gBads_x(3)   =gB_ads_y
        gBads_xx(3,3)=gB_ads_yy

        ! give values to the metric deviation, using sin(theta1)=sin(theta2)=1 w.l.o.g 
        !(considering theta1,theta2-independent case, so theta1=theta2=pi/2 will do)
        h0_ll(1,1)=gb_tt0*(1-x0**2)
        h0_ll(1,2)=gb_tx0*(1-x0**2)
        h0_ll(1,3)=gb_ty0*(1-x0**2)**2  
        h0_ll(2,2)=gb_xx0*(1-x0**2)
        h0_ll(2,3)=gb_xy0*(1-x0**2)**2  
        h0_ll(3,3)=gb_yy0*(1-x0**2)**3
        
        h0_uu(1,1)=g0_uu(1,1)-gads_uu(1,1)
        h0_uu(1,2)=g0_uu(1,2)
        h0_uu(1,3)=g0_uu(1,3)           
        h0_uu(2,2)=g0_uu(2,2)-gads_uu(2,2)
        h0_uu(2,3)=g0_uu(2,3)           
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

        h0_ll_x(1,3,1)   =g0_ll_x(1,3,1)-gads_ll_x(1,3,1) 
        h0_ll_x(1,3,2)   =g0_ll_x(1,3,2)-gads_ll_x(1,3,2)
        h0_ll_xx(1,3,1,1)=g0_ll_xx(1,3,1,1)-gads_ll_xx(1,3,1,1)
        h0_ll_xx(1,3,1,2)=g0_ll_xx(1,3,1,2)-gads_ll_xx(1,3,1,2)
        h0_ll_xx(1,3,2,2)=g0_ll_xx(1,3,2,2)-gads_ll_xx(1,3,2,2)

        h0_ll_x(2,2,1)   =g0_ll_x(2,2,1)-gads_ll_x(2,2,1)
        h0_ll_x(2,2,2)   =g0_ll_x(2,2,2)-gads_ll_x(2,2,2)
        h0_ll_xx(2,2,1,1)=g0_ll_xx(2,2,1,1)-gads_ll_xx(2,2,1,1)
        h0_ll_xx(2,2,1,2)=g0_ll_xx(2,2,1,2)-gads_ll_xx(2,2,1,2)
        h0_ll_xx(2,2,2,2)=g0_ll_xx(2,2,2,2)-gads_ll_xx(2,2,2,2)

        h0_ll_x(2,3,1)   =g0_ll_x(2,3,1)-gads_ll_x(2,3,1) 
        h0_ll_x(2,3,2)   =g0_ll_x(2,3,2)-gads_ll_x(2,3,2)
        h0_ll_xx(2,3,1,1)=g0_ll_xx(2,3,1,1)-gads_ll_xx(2,3,1,1)
        h0_ll_xx(2,3,1,2)=g0_ll_xx(2,3,1,2)-gads_ll_xx(2,3,1,2)
        h0_ll_xx(2,3,2,2)=g0_ll_xx(2,3,2,2)-gads_ll_xx(2,3,2,2)

        h0_ll_x(3,3,1)   =g0_ll_x(3,3,1)-gads_ll_x(3,3,1)
        h0_ll_x(3,3,2)   =g0_ll_x(3,3,2)-gads_ll_x(3,3,2)
        h0_ll_xx(3,3,1,1)=g0_ll_xx(3,3,1,1)-gads_ll_xx(3,3,1,1)
        h0_ll_xx(3,3,1,2)=g0_ll_xx(3,3,1,2)-gads_ll_xx(3,3,1,2)
        h0_ll_xx(3,3,2,2)=g0_ll_xx(3,3,2,2)-gads_ll_xx(3,3,2,2)

        h0_uu_x(1,1,1)=g0_uu_x(1,1,1)-gads_uu_x(1,1,1)
        h0_uu_x(1,1,2)=g0_uu_x(1,1,2)-gads_uu_x(1,1,2)

        h0_uu_x(1,2,1)=g0_uu_x(1,2,1)-gads_uu_x(1,2,1)
        h0_uu_x(1,2,2)=g0_uu_x(1,2,2)-gads_uu_x(1,2,2)

        h0_uu_x(1,3,1)=g0_uu_x(1,3,1)-gads_uu_x(1,3,1) 
        h0_uu_x(1,3,2)=g0_uu_x(1,3,2)-gads_uu_x(1,3,2)

        h0_uu_x(2,2,1)=g0_uu_x(2,2,1)-gads_uu_x(2,2,1)
        h0_uu_x(2,2,2)=g0_uu_x(2,2,2)-gads_uu_x(2,2,2)

        h0_uu_x(2,3,1)=g0_uu_x(2,3,1)-gads_uu_x(2,3,1) 
        h0_uu_x(2,3,2)=g0_uu_x(2,3,2)-gads_uu_x(2,3,2)

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

        hA=psi0*(1-x0**2)*x0**2
        hB=omega0*(1-x0**2)**3*sin(PI*y0/L)**2 

        hAu=1/gA-1/gAads
        hBu=1/gB-1/gBads

        hA_x(1)   =gA_x(1)-gAads_x(1) 
        hA_x(2)   =gA_x(2)-gAads_x(2)
        hA_x(3)   =gA_x(3)-gAads_x(3)
        hA_xx(1,1)=gA_xx(1,1)-gAads_xx(1,1) 
        hA_xx(1,2)=gA_xx(1,2)-gAads_xx(1,2)
        hA_xx(1,3)=gA_xx(1,3)-gAads_xx(1,3)
        hA_xx(2,2)=gA_xx(2,2)-gAads_xx(2,2)
        hA_xx(2,3)=gA_xx(2,3)-gAads_xx(2,3)
        hA_xx(3,3)=gA_xx(3,3)-gAads_xx(3,3)

        hB_x(1)   =gB_x(1)-gBads_x(1) 
        hB_x(2)   =gB_x(2)-gBads_x(2)
        hB_x(3)   =gB_x(3)-gBads_x(3)
        hB_xx(1,1)=gB_xx(1,1)-gBads_xx(1,1) 
        hB_xx(1,2)=gB_xx(1,2)-gBads_xx(1,2)
        hB_xx(1,3)=gB_xx(1,3)-gBads_xx(1,3)
        hB_xx(2,2)=gB_xx(2,2)-gBads_xx(2,2)
        hB_xx(2,3)=gB_xx(2,3)-gBads_xx(2,3)
        hB_xx(3,3)=gB_xx(3,3)-gBads_xx(3,3)

        ! give values to the gh source functions
        A_l(1)=Hb_t0*(1-x0**2)**2 
        A_l(2)=Hb_x0*(1-x0**2)**2
        A_l(3)=Hb_y0*(1-x0**2)**3

        A_l_x(1,1)=Hb_t_t*(1-x0**2)**2
        A_l_x(1,2)=Hb_t_x*(1-x0**2)**2
     &            -2*x0*2*(1-x0**2)*Hb_t0
        A_l_x(1,3)=Hb_t_y*(1-x0**2)**2 

        A_l_x(2,1)=Hb_x_t*(1-x0**2)**2
        A_l_x(2,2)=Hb_x_x*(1-x0**2)**2
     &            -2*x0*2*(1-x0**2)*Hb_x0
        A_l_x(2,3)=Hb_x_y*(1-x0**2)**2 

        A_l_x(3,1)=Hb_y_t*(1-x0**2)**3
        A_l_x(3,2)=Hb_y_x*(1-x0**2)**3
     &            -2*x0*3*(1-x0**2)**2*Hb_y0
        A_l_x(3,3)=Hb_y_y*(1-x0**2)**3

        ! give values to the ads gh source functions
        Hads_l(2)=(-2+4*x0)/(1-2*x0*(1-x0))

        ! give values to the ads gh source function derivatives
        Hads_l_x(2,2)=8*x0*(1-x0)/(1-2*x0*(1-x0))**2

        ! give values to the scalar field
        phi10_x(1)=phi1_t*(1-x0**2)**3
        phi10_x(2)=phi1_x*(1-x0**2)**3
     &            +phi10*(-6*x0)*(1-x0**2)**2
        phi10_x(3)=phi1_y*(1-x0**2)**3     

        phi10_xx(1,1)=phi1_tt*(1-x0**2)**3
        phi10_xx(1,2)=phi1_tx*(1-x0**2)**3
     &               +phi1_t*(-6*x0)*(1-x0**2)**2
        phi10_xx(1,3)=phi1_ty*(1-x0**2)**3
        phi10_xx(2,2)=phi1_xx*(1-x0**2)**3
     &               +phi1_x*(2)*(-6*x0)*(1-x0**2)**2
     &               +phi10*(-6*(1-x0**2)**2+24*x0**2*(1-x0**2))
        phi10_xx(2,3)=phi1_xy*(1-x0**2)**3 
     &               +phi1_y*(-6*x0)*(1-x0**2)**2
        phi10_xx(3,3)=phi1_yy*(1-x0**2)**3 

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
                else if ((a.eq.2.and.b.eq.1.and.c.eq.3)
     &          .or.(a.eq.1.and.b.eq.3.and.c.eq.2)
     &          .or.(a.eq.3.and.b.eq.2.and.c.eq.1))
     &          then
                  levicivi3(a,b,c)=-1
                else 
                  levicivi3(a,b,c)=0
                end if
            end do
          end do
        end do

       ! define sqrt(|detg|) of metric in (t,x,y) subsector
       sqrtdetg=sqrt(abs(g0_ll(1,1)*g0_ll(2,2)*g0_ll(3,3)
     &                  +g0_ll(1,2)*g0_ll(2,3)*g0_ll(3,1)
     &                  +g0_ll(1,3)*g0_ll(2,1)*g0_ll(3,2)
     &                  -g0_ll(1,3)*g0_ll(2,2)*g0_ll(3,1)
     &                  -g0_ll(1,2)*g0_ll(2,1)*g0_ll(3,3)
     &                  -g0_ll(1,1)*g0_ll(2,3)*g0_ll(3,2)))
       sqrtdetgads=sqrt(abs(gads_ll(1,1)*gads_ll(2,2)*gads_ll(3,3)
     &                     +gads_ll(1,2)*gads_ll(2,3)*gads_ll(3,1)
     &                     +gads_ll(1,3)*gads_ll(2,1)*gads_ll(3,2)
     &                     -gads_ll(1,3)*gads_ll(2,2)*gads_ll(3,1)
     &                     -gads_ll(1,2)*gads_ll(2,1)*gads_ll(3,3)
     &                     -gads_ll(1,1)*gads_ll(2,3)*gads_ll(3,2)))
       sqrtdeth=sqrtdetg-sqrtdetgads

       ! define derivatives sqrt(|detg|),x = 1/2/sqrt(|detg|) (detg/|detg|) detg,x
       do a=1,3
         sqrtdetg_x(a)= (g0_ll_x(1,1,a)*g0_ll(2,2)*g0_ll(3,3)
     &                  +g0_ll_x(1,2,a)*g0_ll(2,3)*g0_ll(3,1)
     &                  +g0_ll_x(1,3,a)*g0_ll(2,1)*g0_ll(3,2)
     &                  -g0_ll_x(1,3,a)*g0_ll(2,2)*g0_ll(3,1)
     &                  -g0_ll_x(1,2,a)*g0_ll(2,1)*g0_ll(3,3)
     &                  -g0_ll_x(1,1,a)*g0_ll(2,3)*g0_ll(3,2)
     &                  +g0_ll(1,1)*g0_ll_x(2,2,a)*g0_ll(3,3)
     &                  +g0_ll(1,2)*g0_ll_x(2,3,a)*g0_ll(3,1)
     &                  +g0_ll(1,3)*g0_ll_x(2,1,a)*g0_ll(3,2)
     &                  -g0_ll(1,3)*g0_ll_x(2,2,a)*g0_ll(3,1)
     &                  -g0_ll(1,2)*g0_ll_x(2,1,a)*g0_ll(3,3)
     &                  -g0_ll(1,1)*g0_ll_x(2,3,a)*g0_ll(3,2)
     &                  +g0_ll(1,1)*g0_ll(2,2)*g0_ll_x(3,3,a)
     &                  +g0_ll(1,2)*g0_ll(2,3)*g0_ll_x(3,1,a)
     &                  +g0_ll(1,3)*g0_ll(2,1)*g0_ll_x(3,2,a)
     &                  -g0_ll(1,3)*g0_ll(2,2)*g0_ll_x(3,1,a)
     &                  -g0_ll(1,2)*g0_ll(2,1)*g0_ll_x(3,3,a)
     &                  -g0_ll(1,1)*g0_ll(2,3)*g0_ll_x(3,2,a))
     &                 *(g0_ll(1,1)*g0_ll(2,2)*g0_ll(3,3)
     &                  +g0_ll(1,2)*g0_ll(2,3)*g0_ll(3,1)
     &                  +g0_ll(1,3)*g0_ll(2,1)*g0_ll(3,2)
     &                  -g0_ll(1,3)*g0_ll(2,2)*g0_ll(3,1)
     &                  -g0_ll(1,2)*g0_ll(2,1)*g0_ll(3,3)
     &                  -g0_ll(1,1)*g0_ll(2,3)*g0_ll(3,2))
     &              /abs(g0_ll(1,1)*g0_ll(2,2)*g0_ll(3,3)
     &                  +g0_ll(1,2)*g0_ll(2,3)*g0_ll(3,1)
     &                  +g0_ll(1,3)*g0_ll(2,1)*g0_ll(3,2)
     &                  -g0_ll(1,3)*g0_ll(2,2)*g0_ll(3,1)
     &                  -g0_ll(1,2)*g0_ll(2,1)*g0_ll(3,3)
     &                  -g0_ll(1,1)*g0_ll(2,3)*g0_ll(3,2))
     &              /2.0d0/sqrtdetg
         sqrtdetgads_x(a)= (gads_ll_x(1,1,a)*gads_ll(2,2)*gads_ll(3,3)
     &                     +gads_ll_x(1,2,a)*gads_ll(2,3)*gads_ll(3,1)
     &                     +gads_ll_x(1,3,a)*gads_ll(2,1)*gads_ll(3,2)
     &                     -gads_ll_x(1,3,a)*gads_ll(2,2)*gads_ll(3,1)
     &                     -gads_ll_x(1,2,a)*gads_ll(2,1)*gads_ll(3,3)
     &                     -gads_ll_x(1,1,a)*gads_ll(2,3)*gads_ll(3,2)
     &                     +gads_ll(1,1)*gads_ll_x(2,2,a)*gads_ll(3,3)
     &                     +gads_ll(1,2)*gads_ll_x(2,3,a)*gads_ll(3,1)
     &                     +gads_ll(1,3)*gads_ll_x(2,1,a)*gads_ll(3,2)
     &                     -gads_ll(1,3)*gads_ll_x(2,2,a)*gads_ll(3,1)
     &                     -gads_ll(1,2)*gads_ll_x(2,1,a)*gads_ll(3,3)
     &                     -gads_ll(1,1)*gads_ll_x(2,3,a)*gads_ll(3,2)
     &                     +gads_ll(1,1)*gads_ll(2,2)*gads_ll_x(3,3,a)
     &                     +gads_ll(1,2)*gads_ll(2,3)*gads_ll_x(3,1,a)
     &                     +gads_ll(1,3)*gads_ll(2,1)*gads_ll_x(3,2,a)
     &                     -gads_ll(1,3)*gads_ll(2,2)*gads_ll_x(3,1,a)
     &                     -gads_ll(1,2)*gads_ll(2,1)*gads_ll_x(3,3,a)
     &                     -gads_ll(1,1)*gads_ll(2,3)*gads_ll_x(3,2,a))
     &                    *(gads_ll(1,1)*gads_ll(2,2)*gads_ll(3,3)
     &                     +gads_ll(1,2)*gads_ll(2,3)*gads_ll(3,1)
     &                     +gads_ll(1,3)*gads_ll(2,1)*gads_ll(3,2)
     &                     -gads_ll(1,3)*gads_ll(2,2)*gads_ll(3,1)
     &                     -gads_ll(1,2)*gads_ll(2,1)*gads_ll(3,3)
     &                     -gads_ll(1,1)*gads_ll(2,3)*gads_ll(3,2))
     &                 /abs(gads_ll(1,1)*gads_ll(2,2)*gads_ll(3,3)
     &                     +gads_ll(1,2)*gads_ll(2,3)*gads_ll(3,1)
     &                     +gads_ll(1,3)*gads_ll(2,1)*gads_ll(3,2)
     &                     -gads_ll(1,3)*gads_ll(2,2)*gads_ll(3,1)
     &                     -gads_ll(1,2)*gads_ll(2,1)*gads_ll(3,3)
     &                     -gads_ll(1,1)*gads_ll(2,3)*gads_ll(3,2))
     &                 /2.0d0/sqrtdetgads
       end do

       ! define volume form in in (t,x,y) subsector
        do a=1,3
          do b=1,3
            do c=1,3
              vol(a,b,c)=levicivi3(a,b,c)*sqrtdetg
              volads(a,b,c)=levicivi3(a,b,c)*sqrtdetgads
              volh(a,b,c)=vol(a,b,c)-volads(a,b,c)
              do d=1,3
                vol_x(a,b,c,d)=levicivi3(a,b,c)*sqrtdetg_x(d)
                volads_x(a,b,c,d)=levicivi3(a,b,c)*sqrtdetgads_x(d)
                volh_x(a,b,c,d)=vol_x(a,b,c,d)-volads_x(a,b,c,d)
              end do
            end do
          end do
        end do

        ! give values to the field strength, using sin(phi2)=sin(phi3)=sin(phi4)=1 w.l.o.g 
        !(considering phi2,phi3,phi4-independent case, so phi2=phi3=phi4=pi/2 slice will do)
        f1_l(1)   =0        +fb_t0*(1-x0**2)**3
        f1_l(2)   =0        +fb_x0*(1-x0**2)**3
        f1_l(3)   =f1_y_ads0+fb_y0*(1-x0**2)**3
        f2_ll(1,2)=-vol(1,2,3)*f1_l(1)*g0_uu(3,1)
     &             -vol(1,2,3)*f1_l(2)*g0_uu(3,2)
     &             -vol(1,2,3)*f1_l(3)*g0_uu(3,3) 
        f2_ll(1,3)=-vol(1,3,2)*f1_l(1)*g0_uu(2,1)
     &             -vol(1,3,2)*f1_l(2)*g0_uu(2,2)
     &             -vol(1,3,2)*f1_l(3)*g0_uu(2,3)
        f2_ll(2,3)=-vol(2,3,1)*f1_l(1)*g0_uu(1,1)
     &             -vol(2,3,1)*f1_l(2)*g0_uu(1,2)
     &             -vol(2,3,1)*f1_l(3)*g0_uu(1,3)

        f1_l_x(1,1)=fb_t_t*(1-x0**2)**3
        f1_l_x(1,2)=fb_t_x*(1-x0**2)**3
     &             +fb_t0*3*(1-x0**2)**2*(-2*x0)
        f1_l_x(1,3)=fb_t_y*(1-x0**2)**3
        f1_l_x(2,1)=fb_x_t*(1-x0**2)**3
        f1_l_x(2,2)=fb_x_x*(1-x0**2)**3
     &             +fb_x0*3*(1-x0**2)**2*(-2*x0)
        f1_l_x(2,3)=fb_x_y*(1-x0**2)**3
        f1_l_x(3,1)=fb_y_t*(1-x0**2)**3
        f1_l_x(3,2)=fb_y_x*(1-x0**2)**3
     &             +fb_y0*3*(1-x0**2)**2*(-2*x0)
        f1_l_x(3,3)=fb_y_y*(1-x0**2)**3

        do c=1,3
          f2_ll_x(1,2,c)=-vol_x(1,2,3,c)*f1_l(1)*g0_uu(3,1)
     &                   -vol_x(1,2,3,c)*f1_l(2)*g0_uu(3,2)
     &                   -vol_x(1,2,3,c)*f1_l(3)*g0_uu(3,3)
     &                   -vol(1,2,3)*f1_l_x(1,c)*g0_uu(3,1)
     &                   -vol(1,2,3)*f1_l_x(2,c)*g0_uu(3,2)
     &                   -vol(1,2,3)*f1_l_x(3,c)*g0_uu(3,3)
     &                   -vol(1,2,3)*f1_l(1)*g0_uu_x(3,1,c)
     &                   -vol(1,2,3)*f1_l(2)*g0_uu_x(3,2,c)
     &                   -vol(1,2,3)*f1_l(3)*g0_uu_x(3,3,c)
          f2_ll_x(1,3,c)=-vol_x(1,3,2,c)*f1_l(1)*g0_uu(2,1)
     &                   -vol_x(1,3,2,c)*f1_l(2)*g0_uu(2,2)
     &                   -vol_x(1,3,2,c)*f1_l(3)*g0_uu(2,3)
     &                   -vol(1,3,2)*f1_l_x(1,c)*g0_uu(2,1)
     &                   -vol(1,3,2)*f1_l_x(2,c)*g0_uu(2,2)
     &                   -vol(1,3,2)*f1_l_x(3,c)*g0_uu(2,3)
     &                   -vol(1,3,2)*f1_l(1)*g0_uu_x(2,1,c)
     &                   -vol(1,3,2)*f1_l(2)*g0_uu_x(2,2,c)
     &                   -vol(1,3,2)*f1_l(3)*g0_uu_x(2,3,c)
          f2_ll_x(2,3,c)=-vol_x(2,3,1,c)*f1_l(1)*g0_uu(1,1)
     &                   -vol_x(2,3,1,c)*f1_l(2)*g0_uu(1,2)
     &                   -vol_x(2,3,1,c)*f1_l(3)*g0_uu(1,3)
     &                   -vol(2,3,1)*f1_l_x(1,c)*g0_uu(1,1)
     &                   -vol(2,3,1)*f1_l_x(2,c)*g0_uu(1,2)
     &                   -vol(2,3,1)*f1_l_x(3,c)*g0_uu(1,3)
     &                   -vol(2,3,1)*f1_l(1)*g0_uu_x(1,1,c)
     &                   -vol(2,3,1)*f1_l(2)*g0_uu_x(1,2,c)
     &                   -vol(2,3,1)*f1_l(3)*g0_uu_x(1,3,c)
        end do

        f1ads_l(3)   =f1_y_ads0

        f2ads_ll(1,2)=-volads(1,2,3)*f1ads_l(3)*gads_uu(3,3)
        f2ads_ll(1,3)=-volads(1,3,2)*f1ads_l(3)*gads_uu(2,3)
        f2ads_ll(2,3)=-volads(2,3,1)*f1ads_l(3)*gads_uu(1,3)

        do c=1,3
          f2ads_ll_x(1,2,c)=-volads_x(1,2,3,c)*f1ads_l(3)*gads_uu(3,3)
     &                      -volads(1,2,3)*f1ads_l(3)*gads_uu_x(3,3,c)
          f2ads_ll_x(1,3,c)=-volads_x(1,3,2,c)*f1ads_l(3)*gads_uu(2,3)
     &                      -volads(1,3,2)*f1ads_l(3)*gads_uu_x(2,3,c)
          f2ads_ll_x(2,3,c)=-volads_x(2,3,1,c)*f1ads_l(3)*gads_uu(1,3)
     &                      -volads(2,3,1)*f1ads_l(3)*gads_uu_x(1,3,c)
        end do

        h1_l(1)   =fb_t0*(1-x0**2)**3
        h1_l(2)   =fb_x0*(1-x0**2)**3
        h1_l(3)   =fb_y0*(1-x0**2)**3
        h2_ll(1,2)=f2_ll(1,2)-f2ads_ll(1,2) 
        h2_ll(1,3)=f2_ll(1,3)-f2ads_ll(1,3)
        h2_ll(2,3)=f2_ll(2,3)-f2ads_ll(2,3)

        h1_l_x(1,1)=fb_t_t*(1-x0**2)**3
        h1_l_x(1,2)=fb_t_x*(1-x0**2)**3
        h1_l_x(1,3)=fb_t_y*(1-x0**2)**3
        h1_l_x(2,1)=fb_x_t*(1-x0**2)**3
        h1_l_x(2,2)=fb_x_x*(1-x0**2)**3
        h1_l_x(2,3)=fb_x_y*(1-x0**2)**3
        h1_l_x(3,1)=fb_y_t*(1-x0**2)**3
        h1_l_x(3,2)=fb_y_x*(1-x0**2)**3
        h1_l_x(3,3)=fb_y_y*(1-x0**2)**3

        do c=1,3
          h2_ll_x(1,2,c)=f2_ll_x(1,2,c)-f2ads_ll_x(1,2,c) 
          h2_ll_x(1,3,c)=f2_ll_x(1,3,c)-f2ads_ll_x(1,3,c)
          h2_ll_x(2,3,c)=f2_ll_x(2,3,c)-f2ads_ll_x(2,3,c)
        end do

        do a=1,2
          do b=a+1,3
            f2_ll(b,a)=-f2_ll(a,b)
            f2ads_ll(b,a)=-f2ads_ll(a,b)
            h2_ll(b,a)=-h2_ll(a,b)
            do c=1,3
              f2_ll_x(b,a,c)=f2_ll_x(a,b,c)
              f2ads_ll_x(b,a,c)=f2ads_ll_x(a,b,c)
              h2_ll_x(b,a,c)=h2_ll_x(a,b,c)
            end do
          end do
        end do

        ! calculate Christoffel symbol derivatives at point i
        !(gamma^a_bc,e = 1/2 g^ad_,e(g_bd,c  + g_cd,b  - g_bc,d)
        !              +   1/2 g^ad(g_bd,ce + g_cd,be - g_bc,de))
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

!TEST!
!        gB_x(1)=0
!        gB_xx(1,3)=0

        ! calculate ricci_ll,sA,sB,tA,tB
        do a=1,3
          do b=1,3
            ricci_ll(a,b)=riccibar_ll(a,b)
     &                   -dimA*( (gA_xx(a,b)
     &                           -gamma_ull(1,b,a)*gA_x(1)
     &                           -gamma_ull(2,b,a)*gA_x(2)
     &                           -gamma_ull(3,b,a)*gA_x(3)
     &                           )/(2*gA)
     &                         - (gA_x(a)*gA_x(b))/(4*gA**2) )
     &                   -dimB*( (gB_xx(a,b)
     &                           -gamma_ull(1,b,a)*gB_x(1)
     &                           -gamma_ull(2,b,a)*gB_x(2)
     &                           -gamma_ull(3,b,a)*gB_x(3)
     &                           )/(2*gB)
     &                         - (gB_x(a)*gB_x(b))/(4*gB**2) )
          end do
        end do

        sA=-g0_uu(1,1)*gA_xx(1,1)/2
     &     -g0_uu(1,2)*gA_xx(1,2)/2
     &     -g0_uu(1,3)*gA_xx(1,3)/2
     &     -g0_uu(2,1)*gA_xx(2,1)/2
     &     -g0_uu(2,2)*gA_xx(2,2)/2
     &     -g0_uu(2,3)*gA_xx(2,3)/2
     &     -g0_uu(3,1)*gA_xx(3,1)/2
     &     -g0_uu(3,2)*gA_xx(3,2)/2
     &     -g0_uu(3,3)*gA_xx(3,3)/2
     &     +g0_uu(1,1)*gamma_ull(1,1,1)*gA_x(1)/2 
     &     +g0_uu(1,2)*gamma_ull(1,2,1)*gA_x(1)/2 
     &     +g0_uu(1,3)*gamma_ull(1,3,1)*gA_x(1)/2 
     &     +g0_uu(2,1)*gamma_ull(1,1,2)*gA_x(1)/2 
     &     +g0_uu(2,2)*gamma_ull(1,2,2)*gA_x(1)/2 
     &     +g0_uu(2,3)*gamma_ull(1,3,2)*gA_x(1)/2 
     &     +g0_uu(3,1)*gamma_ull(1,1,3)*gA_x(1)/2 
     &     +g0_uu(3,2)*gamma_ull(1,2,3)*gA_x(1)/2 
     &     +g0_uu(3,3)*gamma_ull(1,3,3)*gA_x(1)/2 
     &     +g0_uu(1,1)*gamma_ull(2,1,1)*gA_x(2)/2 
     &     +g0_uu(1,2)*gamma_ull(2,2,1)*gA_x(2)/2 
     &     +g0_uu(1,3)*gamma_ull(2,3,1)*gA_x(2)/2 
     &     +g0_uu(2,1)*gamma_ull(2,1,2)*gA_x(2)/2 
     &     +g0_uu(2,2)*gamma_ull(2,2,2)*gA_x(2)/2 
     &     +g0_uu(2,3)*gamma_ull(2,3,2)*gA_x(2)/2 
     &     +g0_uu(3,1)*gamma_ull(2,1,3)*gA_x(2)/2 
     &     +g0_uu(3,2)*gamma_ull(2,2,3)*gA_x(2)/2 
     &     +g0_uu(3,3)*gamma_ull(2,3,3)*gA_x(2)/2 
     &     +g0_uu(1,1)*gamma_ull(3,1,1)*gA_x(3)/2 
     &     +g0_uu(1,2)*gamma_ull(3,2,1)*gA_x(3)/2 
     &     +g0_uu(1,3)*gamma_ull(3,3,1)*gA_x(3)/2 
     &     +g0_uu(2,1)*gamma_ull(3,1,2)*gA_x(3)/2 
     &     +g0_uu(2,2)*gamma_ull(3,2,2)*gA_x(3)/2 
     &     +g0_uu(2,3)*gamma_ull(3,3,2)*gA_x(3)/2 
     &     +g0_uu(3,1)*gamma_ull(3,1,3)*gA_x(3)/2 
     &     +g0_uu(3,2)*gamma_ull(3,2,3)*gA_x(3)/2 
     &     +g0_uu(3,3)*gamma_ull(3,3,3)*gA_x(3)/2 
     &     -g0_uu(1,1)*(gA_x(1)*gA_x(1))/(4*gA)*(dimA-2d0)
     &     -g0_uu(1,2)*(gA_x(1)*gA_x(2))/(4*gA)*(dimA-2d0)
     &     -g0_uu(1,3)*(gA_x(1)*gA_x(3))/(4*gA)*(dimA-2d0)
     &     -g0_uu(2,1)*(gA_x(2)*gA_x(1))/(4*gA)*(dimA-2d0)
     &     -g0_uu(2,2)*(gA_x(2)*gA_x(2))/(4*gA)*(dimA-2d0)
     &     -g0_uu(2,3)*(gA_x(2)*gA_x(3))/(4*gA)*(dimA-2d0)
     &     -g0_uu(3,1)*(gA_x(3)*gA_x(1))/(4*gA)*(dimA-2d0)
     &     -g0_uu(3,2)*(gA_x(3)*gA_x(2))/(4*gA)*(dimA-2d0)
     &     -g0_uu(3,3)*(gA_x(3)*gA_x(3))/(4*gA)*(dimA-2d0)
     &     -g0_uu(1,1)*(gA_x(1)*gB_x(1))/(4*gB)*dimB
     &     -g0_uu(1,2)*(gA_x(1)*gB_x(2))/(4*gB)*dimB
     &     -g0_uu(1,3)*(gA_x(1)*gB_x(3))/(4*gB)*dimB
     &     -g0_uu(2,1)*(gA_x(2)*gB_x(1))/(4*gB)*dimB
     &     -g0_uu(2,2)*(gA_x(2)*gB_x(2))/(4*gB)*dimB
     &     -g0_uu(2,3)*(gA_x(2)*gB_x(3))/(4*gB)*dimB
     &     -g0_uu(3,1)*(gA_x(3)*gB_x(1))/(4*gB)*dimB
     &     -g0_uu(3,2)*(gA_x(3)*gB_x(2))/(4*gB)*dimB
     &     -g0_uu(3,3)*(gA_x(3)*gB_x(3))/(4*gB)*dimB

        sB=-g0_uu(1,1)*gB_xx(1,1)/2
     &     -g0_uu(1,2)*gB_xx(1,2)/2
     &     -g0_uu(1,3)*gB_xx(1,3)/2
     &     -g0_uu(2,1)*gB_xx(2,1)/2
     &     -g0_uu(2,2)*gB_xx(2,2)/2
     &     -g0_uu(2,3)*gB_xx(2,3)/2
     &     -g0_uu(3,1)*gB_xx(3,1)/2
     &     -g0_uu(3,2)*gB_xx(3,2)/2
     &     -g0_uu(3,3)*gB_xx(3,3)/2
     &     +g0_uu(1,1)*gamma_ull(1,1,1)*gB_x(1)/2 
     &     +g0_uu(1,2)*gamma_ull(1,2,1)*gB_x(1)/2 
     &     +g0_uu(1,3)*gamma_ull(1,3,1)*gB_x(1)/2 
     &     +g0_uu(2,1)*gamma_ull(1,1,2)*gB_x(1)/2 
     &     +g0_uu(2,2)*gamma_ull(1,2,2)*gB_x(1)/2 
     &     +g0_uu(2,3)*gamma_ull(1,3,2)*gB_x(1)/2 
     &     +g0_uu(3,1)*gamma_ull(1,1,3)*gB_x(1)/2 
     &     +g0_uu(3,2)*gamma_ull(1,2,3)*gB_x(1)/2 
     &     +g0_uu(3,3)*gamma_ull(1,3,3)*gB_x(1)/2 
     &     +g0_uu(1,1)*gamma_ull(2,1,1)*gB_x(2)/2 
     &     +g0_uu(1,2)*gamma_ull(2,2,1)*gB_x(2)/2 
     &     +g0_uu(1,3)*gamma_ull(2,3,1)*gB_x(2)/2 
     &     +g0_uu(2,1)*gamma_ull(2,1,2)*gB_x(2)/2 
     &     +g0_uu(2,2)*gamma_ull(2,2,2)*gB_x(2)/2 
     &     +g0_uu(2,3)*gamma_ull(2,3,2)*gB_x(2)/2 
     &     +g0_uu(3,1)*gamma_ull(2,1,3)*gB_x(2)/2 
     &     +g0_uu(3,2)*gamma_ull(2,2,3)*gB_x(2)/2 
     &     +g0_uu(3,3)*gamma_ull(2,3,3)*gB_x(2)/2 
     &     +g0_uu(1,1)*gamma_ull(3,1,1)*gB_x(3)/2 
     &     +g0_uu(1,2)*gamma_ull(3,2,1)*gB_x(3)/2 
     &     +g0_uu(1,3)*gamma_ull(3,3,1)*gB_x(3)/2 
     &     +g0_uu(2,1)*gamma_ull(3,1,2)*gB_x(3)/2 
     &     +g0_uu(2,2)*gamma_ull(3,2,2)*gB_x(3)/2 
     &     +g0_uu(2,3)*gamma_ull(3,3,2)*gB_x(3)/2 
     &     +g0_uu(3,1)*gamma_ull(3,1,3)*gB_x(3)/2 
     &     +g0_uu(3,2)*gamma_ull(3,2,3)*gB_x(3)/2 
     &     +g0_uu(3,3)*gamma_ull(3,3,3)*gB_x(3)/2 
     &     -g0_uu(1,1)*(gB_x(1)*gB_x(1))/(4*gB)*(dimB-2d0)
     &     -g0_uu(1,2)*(gB_x(1)*gB_x(2))/(4*gB)*(dimB-2d0)
     &     -g0_uu(1,3)*(gB_x(1)*gB_x(3))/(4*gB)*(dimB-2d0)
     &     -g0_uu(2,1)*(gB_x(2)*gB_x(1))/(4*gB)*(dimB-2d0)
     &     -g0_uu(2,2)*(gB_x(2)*gB_x(2))/(4*gB)*(dimB-2d0)
     &     -g0_uu(2,3)*(gB_x(2)*gB_x(3))/(4*gB)*(dimB-2d0)
     &     -g0_uu(3,1)*(gB_x(3)*gB_x(1))/(4*gB)*(dimB-2d0)
     &     -g0_uu(3,2)*(gB_x(3)*gB_x(2))/(4*gB)*(dimB-2d0)
     &     -g0_uu(3,3)*(gB_x(3)*gB_x(3))/(4*gB)*(dimB-2d0)
     &     -g0_uu(1,1)*(gA_x(1)*gB_x(1))/(4*gA)*dimA 
     &     -g0_uu(1,2)*(gA_x(1)*gB_x(2))/(4*gA)*dimA
     &     -g0_uu(1,3)*(gA_x(1)*gB_x(3))/(4*gA)*dimA
     &     -g0_uu(2,1)*(gA_x(2)*gB_x(1))/(4*gA)*dimA
     &     -g0_uu(2,2)*(gA_x(2)*gB_x(2))/(4*gA)*dimA
     &     -g0_uu(2,3)*(gA_x(2)*gB_x(3))/(4*gA)*dimA
     &     -g0_uu(3,1)*(gA_x(3)*gB_x(1))/(4*gA)*dimA
     &     -g0_uu(3,2)*(gA_x(3)*gB_x(2))/(4*gA)*dimA
     &     -g0_uu(3,3)*(gA_x(3)*gB_x(3))/(4*gA)*dimA

        tA=-gA*g0_uu(1,1)*g0_uu(1,1)*f2_ll(1,1)*f2_ll(1,1)/8
     &     -gA*g0_uu(1,2)*g0_uu(1,1)*f2_ll(1,1)*f2_ll(2,1)/8
     &     -gA*g0_uu(1,3)*g0_uu(1,1)*f2_ll(1,1)*f2_ll(3,1)/8
     &     -gA*g0_uu(2,1)*g0_uu(1,1)*f2_ll(2,1)*f2_ll(1,1)/8
     &     -gA*g0_uu(2,2)*g0_uu(1,1)*f2_ll(2,1)*f2_ll(2,1)/8
     &     -gA*g0_uu(2,3)*g0_uu(1,1)*f2_ll(2,1)*f2_ll(3,1)/8
     &     -gA*g0_uu(3,1)*g0_uu(1,1)*f2_ll(3,1)*f2_ll(1,1)/8
     &     -gA*g0_uu(3,2)*g0_uu(1,1)*f2_ll(3,1)*f2_ll(2,1)/8
     &     -gA*g0_uu(3,3)*g0_uu(1,1)*f2_ll(3,1)*f2_ll(3,1)/8
     &     -gA*g0_uu(1,1)*g0_uu(1,2)*f2_ll(1,1)*f2_ll(1,2)/8 
     &     -gA*g0_uu(1,2)*g0_uu(1,2)*f2_ll(1,1)*f2_ll(2,2)/8 
     &     -gA*g0_uu(1,3)*g0_uu(1,2)*f2_ll(1,1)*f2_ll(3,2)/8 
     &     -gA*g0_uu(2,1)*g0_uu(1,2)*f2_ll(2,1)*f2_ll(1,2)/8 
     &     -gA*g0_uu(2,2)*g0_uu(1,2)*f2_ll(2,1)*f2_ll(2,2)/8 
     &     -gA*g0_uu(2,3)*g0_uu(1,2)*f2_ll(2,1)*f2_ll(3,2)/8 
     &     -gA*g0_uu(3,1)*g0_uu(1,2)*f2_ll(3,1)*f2_ll(1,2)/8 
     &     -gA*g0_uu(3,2)*g0_uu(1,2)*f2_ll(3,1)*f2_ll(2,2)/8 
     &     -gA*g0_uu(3,3)*g0_uu(1,2)*f2_ll(3,1)*f2_ll(3,2)/8 
     &     -gA*g0_uu(1,1)*g0_uu(1,3)*f2_ll(1,1)*f2_ll(1,3)/8
     &     -gA*g0_uu(1,2)*g0_uu(1,3)*f2_ll(1,1)*f2_ll(2,3)/8
     &     -gA*g0_uu(1,3)*g0_uu(1,3)*f2_ll(1,1)*f2_ll(3,3)/8
     &     -gA*g0_uu(2,1)*g0_uu(1,3)*f2_ll(2,1)*f2_ll(1,3)/8
     &     -gA*g0_uu(2,2)*g0_uu(1,3)*f2_ll(2,1)*f2_ll(2,3)/8
     &     -gA*g0_uu(2,3)*g0_uu(1,3)*f2_ll(2,1)*f2_ll(3,3)/8
     &     -gA*g0_uu(3,1)*g0_uu(1,3)*f2_ll(3,1)*f2_ll(1,3)/8
     &     -gA*g0_uu(3,2)*g0_uu(1,3)*f2_ll(3,1)*f2_ll(2,3)/8
     &     -gA*g0_uu(3,3)*g0_uu(1,3)*f2_ll(3,1)*f2_ll(3,3)/8
     &     -gA*g0_uu(1,1)*g0_uu(2,1)*f2_ll(1,2)*f2_ll(1,1)/8
     &     -gA*g0_uu(1,2)*g0_uu(2,1)*f2_ll(1,2)*f2_ll(2,1)/8
     &     -gA*g0_uu(1,3)*g0_uu(2,1)*f2_ll(1,2)*f2_ll(3,1)/8
     &     -gA*g0_uu(2,1)*g0_uu(2,1)*f2_ll(2,2)*f2_ll(1,1)/8
     &     -gA*g0_uu(2,2)*g0_uu(2,1)*f2_ll(2,2)*f2_ll(2,1)/8
     &     -gA*g0_uu(2,3)*g0_uu(2,1)*f2_ll(2,2)*f2_ll(3,1)/8
     &     -gA*g0_uu(3,1)*g0_uu(2,1)*f2_ll(3,2)*f2_ll(1,1)/8
     &     -gA*g0_uu(3,2)*g0_uu(2,1)*f2_ll(3,2)*f2_ll(2,1)/8
     &     -gA*g0_uu(3,3)*g0_uu(2,1)*f2_ll(3,2)*f2_ll(3,1)/8
     &     -gA*g0_uu(1,1)*g0_uu(2,2)*f2_ll(1,2)*f2_ll(1,2)/8
     &     -gA*g0_uu(1,2)*g0_uu(2,2)*f2_ll(1,2)*f2_ll(2,2)/8
     &     -gA*g0_uu(1,3)*g0_uu(2,2)*f2_ll(1,2)*f2_ll(3,2)/8
     &     -gA*g0_uu(2,1)*g0_uu(2,2)*f2_ll(2,2)*f2_ll(1,2)/8
     &     -gA*g0_uu(2,2)*g0_uu(2,2)*f2_ll(2,2)*f2_ll(2,2)/8
     &     -gA*g0_uu(2,3)*g0_uu(2,2)*f2_ll(2,2)*f2_ll(3,2)/8
     &     -gA*g0_uu(3,1)*g0_uu(2,2)*f2_ll(3,2)*f2_ll(1,2)/8
     &     -gA*g0_uu(3,2)*g0_uu(2,2)*f2_ll(3,2)*f2_ll(2,2)/8
     &     -gA*g0_uu(3,3)*g0_uu(2,2)*f2_ll(3,2)*f2_ll(3,2)/8
     &     -gA*g0_uu(1,1)*g0_uu(2,3)*f2_ll(1,2)*f2_ll(1,3)/8
     &     -gA*g0_uu(1,2)*g0_uu(2,3)*f2_ll(1,2)*f2_ll(2,3)/8
     &     -gA*g0_uu(1,3)*g0_uu(2,3)*f2_ll(1,2)*f2_ll(3,3)/8
     &     -gA*g0_uu(2,1)*g0_uu(2,3)*f2_ll(2,2)*f2_ll(1,3)/8
     &     -gA*g0_uu(2,2)*g0_uu(2,3)*f2_ll(2,2)*f2_ll(2,3)/8
     &     -gA*g0_uu(2,3)*g0_uu(2,3)*f2_ll(2,2)*f2_ll(3,3)/8
     &     -gA*g0_uu(3,1)*g0_uu(2,3)*f2_ll(3,2)*f2_ll(1,3)/8
     &     -gA*g0_uu(3,2)*g0_uu(2,3)*f2_ll(3,2)*f2_ll(2,3)/8
     &     -gA*g0_uu(3,3)*g0_uu(2,3)*f2_ll(3,2)*f2_ll(3,3)/8
     &     -gA*g0_uu(1,1)*g0_uu(3,1)*f2_ll(1,3)*f2_ll(1,1)/8
     &     -gA*g0_uu(1,2)*g0_uu(3,1)*f2_ll(1,3)*f2_ll(2,1)/8
     &     -gA*g0_uu(1,3)*g0_uu(3,1)*f2_ll(1,3)*f2_ll(3,1)/8
     &     -gA*g0_uu(2,1)*g0_uu(3,1)*f2_ll(2,3)*f2_ll(1,1)/8
     &     -gA*g0_uu(2,2)*g0_uu(3,1)*f2_ll(2,3)*f2_ll(2,1)/8
     &     -gA*g0_uu(2,3)*g0_uu(3,1)*f2_ll(2,3)*f2_ll(3,1)/8
     &     -gA*g0_uu(3,1)*g0_uu(3,1)*f2_ll(3,3)*f2_ll(1,1)/8
     &     -gA*g0_uu(3,2)*g0_uu(3,1)*f2_ll(3,3)*f2_ll(2,1)/8
     &     -gA*g0_uu(3,3)*g0_uu(3,1)*f2_ll(3,3)*f2_ll(3,1)/8
     &     -gA*g0_uu(1,1)*g0_uu(3,2)*f2_ll(1,3)*f2_ll(1,2)/8
     &     -gA*g0_uu(1,2)*g0_uu(3,2)*f2_ll(1,3)*f2_ll(2,2)/8
     &     -gA*g0_uu(1,3)*g0_uu(3,2)*f2_ll(1,3)*f2_ll(3,2)/8
     &     -gA*g0_uu(2,1)*g0_uu(3,2)*f2_ll(2,3)*f2_ll(1,2)/8
     &     -gA*g0_uu(2,2)*g0_uu(3,2)*f2_ll(2,3)*f2_ll(2,2)/8
     &     -gA*g0_uu(2,3)*g0_uu(3,2)*f2_ll(2,3)*f2_ll(3,2)/8
     &     -gA*g0_uu(3,1)*g0_uu(3,2)*f2_ll(3,3)*f2_ll(1,2)/8
     &     -gA*g0_uu(3,2)*g0_uu(3,2)*f2_ll(3,3)*f2_ll(2,2)/8
     &     -gA*g0_uu(3,3)*g0_uu(3,2)*f2_ll(3,3)*f2_ll(3,2)/8
     &     -gA*g0_uu(1,1)*g0_uu(3,3)*f2_ll(1,3)*f2_ll(1,3)/8
     &     -gA*g0_uu(1,2)*g0_uu(3,3)*f2_ll(1,3)*f2_ll(2,3)/8
     &     -gA*g0_uu(1,3)*g0_uu(3,3)*f2_ll(1,3)*f2_ll(3,3)/8
     &     -gA*g0_uu(2,1)*g0_uu(3,3)*f2_ll(2,3)*f2_ll(1,3)/8
     &     -gA*g0_uu(2,2)*g0_uu(3,3)*f2_ll(2,3)*f2_ll(2,3)/8
     &     -gA*g0_uu(2,3)*g0_uu(3,3)*f2_ll(2,3)*f2_ll(3,3)/8
     &     -gA*g0_uu(3,1)*g0_uu(3,3)*f2_ll(3,3)*f2_ll(1,3)/8
     &     -gA*g0_uu(3,2)*g0_uu(3,3)*f2_ll(3,3)*f2_ll(2,3)/8
     &     -gA*g0_uu(3,3)*g0_uu(3,3)*f2_ll(3,3)*f2_ll(3,3)/8

!TEST! COMMENTING THIS OUT MAKES ONLY-GB_TY and ONLY-GB_XY UNSTABLE,
!      BUT tB, BELOW, DOES NOT APPEAR IN EITHER GB_TY OR GB_XY EQNS!
!        g0_uu(3,3)=
!     &         (g0_ll(1,1)*g0_ll(2,2))!-g0_ll(1,2)**2)
!     &         /(-g0_ll(1,3)**2*g0_ll(2,2)
!!     &           +g0_ll(1,2)*g0_ll(1,3)*g0_ll(2,3)*2
!     &           -g0_ll(1,1)*g0_ll(2,3)**2
!!     &           -g0_ll(1,2)**2*g0_ll(3,3)
!     &           +g0_ll(1,1)*g0_ll(2,2)*g0_ll(3,3))

        tB=-gB*g0_uu(1,1)*f1_l(1)*f1_l(1)/4
     &     -gB*g0_uu(1,2)*f1_l(1)*f1_l(2)/4
     &     -gB*g0_uu(1,3)*f1_l(1)*f1_l(3)/4
     &     -gB*g0_uu(2,1)*f1_l(2)*f1_l(1)/4
     &     -gB*g0_uu(2,2)*f1_l(2)*f1_l(2)/4
     &     -gB*g0_uu(2,3)*f1_l(2)*f1_l(3)/4
     &     -gB*g0_uu(3,1)*f1_l(3)*f1_l(1)/4
     &     -gB*g0_uu(3,2)*f1_l(3)*f1_l(2)/4
     &     -gB*g0_uu(3,3)*f1_l(3)*f1_l(3)/4

!TEST! RESETS BACK TO WHAT g0_UU(3,3) PHYSICALLY SHOULD BE
        g0_uu(3,3)=
     &         (g0_ll(1,1)*g0_ll(2,2)-g0_ll(1,2)**2)
     &         /(-g0_ll(1,3)**2*g0_ll(2,2)
     &           +g0_ll(1,2)*g0_ll(1,3)*g0_ll(2,3)*2
     &           -g0_ll(1,1)*g0_ll(2,3)**2
     &           -g0_ll(1,2)**2*g0_ll(3,3)
     &           +g0_ll(1,1)*g0_ll(2,2)*g0_ll(3,3))

        ! computes auxiliary objects at point i,j,k
        do a=1,3
          do b=1,3
            do c=1,3
              dlll(a,b,c)=
     &            g0_ll_x(b,c,a)-g0_ll_x(a,b,c)+g0_ll_x(c,a,b)
              gammagg(a,b,c)=0
              gammahh(a,b,c)=0
              gammagh(a,b,c)=0
              gammahg(a,b,c)=0
              do d=1,3
                gammagg(a,b,c)=gammagg(a,b,c)
     &                         +0.5d0*gads_uu(a,d)
     &                              *(gads_ll_x(c,d,b)
     &                               -gads_ll_x(b,c,d)
     &                               +gads_ll_x(d,b,c))
                gammahh(a,b,c)=gammahh(a,b,c)
     &                         +0.5d0*h0_uu(a,d)
     &                              *(h0_ll_x(c,d,b)
     &                               -h0_ll_x(b,c,d)
     &                               +h0_ll_x(d,b,c))
                gammagh(a,b,c)=gammagh(a,b,c)
     &                         +0.5d0*gads_uu(a,d)
     &                              *(h0_ll_x(c,d,b)
     &                               -h0_ll_x(b,c,d)
     &                               +h0_ll_x(d,b,c))
                gammahg(a,b,c)=gammahg(a,b,c)
     &                         +0.5d0*h0_uu(a,d)
     &                              *(gads_ll_x(c,d,b)
     &                               -gads_ll_x(b,c,d)
     &                               +gads_ll_x(d,b,c))
                cuuuu(a,b,c,d)=gads_uu(a,b)*gads_uu(c,d)+
     &                         h0_uu(a,b)*h0_uu(c,d)+
     &                         gads_uu(a,b)*h0_uu(c,d)+
     &                         h0_uu(a,b)*gads_uu(c,d)
              end do
            end do
          end do
        end do

!         if (x(i).eq.0.5d0) then
!           write(*,*) '(dimB-1)+sB+tB=',(dimB-1)+sB+tB
!           stop
!         end if
!TEST!
!        phi1_np1(i)=omega_t
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
!                    term8(a,b)=-dimA*( (gA_xx(a,b)
!     &                                 -gamma_ull(1,b,a)*gA_x(1)
!     &                                 -gamma_ull(2,b,a)*gA_x(2)
!     &                                 -gamma_ull(3,b,a)*gA_x(3)
!     &                                 )/(2*gA)
!     &                               - (gA_x(a)*gA_x(b))/(4*gA**2) )
!     &                         -dimB*( (gB_xx(a,b)
!     &                                 -gamma_ull(1,b,a)*gB_x(1)
!     &                                 -gamma_ull(2,b,a)*gB_x(2)
!     &                                 -gamma_ull(3,b,a)*gB_x(3)
!     &                                 )/(2*gB)
!     &                               - (gB_x(a)*gB_x(b))/(4*gB**2) )
!     &
!                    term9(a,b)=-(f1_l(a)*f1_l(b)
!     &                          +g0_uu(1,1)*f2_ll(a,1)*f2_ll(b,1)
!     &                          +g0_uu(1,2)*f2_ll(a,1)*f2_ll(b,2)
!     &                          +g0_uu(1,3)*f2_ll(a,1)*f2_ll(b,3)
!     &                          +g0_uu(2,1)*f2_ll(a,2)*f2_ll(b,1)
!     &                          +g0_uu(2,2)*f2_ll(a,2)*f2_ll(b,2)
!     &                          +g0_uu(2,3)*f2_ll(a,2)*f2_ll(b,3)
!     &                          +g0_uu(3,1)*f2_ll(a,3)*f2_ll(b,1)
!     &                          +g0_uu(3,2)*f2_ll(a,3)*f2_ll(b,2)
!     &                          +g0_uu(3,3)*f2_ll(a,3)*f2_ll(b,3)
!     &                          )/4
!     &
!                  efe(a,b)=term1(a,b)+term2(a,b)+term3(a,b)+term4(a,b)
!     &                    +term5(a,b)+term6(a,b)+term7(a,b)+term8(a,b)
!     &                    +term9(a,b)
!     &
!                  efest(a,b)=riccibar_ll(a,b)+term8(a,b)+term9(a,b)
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
