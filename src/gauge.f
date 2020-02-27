c----------------------------------------------------------------------
c in polar coordinates t==t, x==rho, y==chi/PI
c
c routines associated with the Harmonic source functions Hb_u
c----------------------------------------------------------------------

c----------------------------------------------------------------------
c smooth (C2) transition function from 0 at x=rho1, to 1 at x=rho2
c i.e. trans=0 everywhere for rho1=1.0,rho2=1.0, 
c      trans=1 everywhere for rho1=0.0,rho2=0.0
c----------------------------------------------------------------------
        real*8 function trans(x,rho1,rho2)
        implicit none
        real*8 x,rho1,rho2

        real*8 xb

        ! initialize fixed-size variables 
        data xb/0.0/

        !--------------------------------------------------------------
 
        xb=(rho2-x)/(rho2-rho1)

        if (x.ge.rho2) then
          trans=1
        else if (x.ge.rho1) then
          trans=1-xb**3*(6*xb**2-15*xb+10)
        else
          trans=0
        end if

        return
        end

c----------------------------------------------------------------------
c smooth (C-inf) transition function from 0 at x=rho1, to 1 at x=rho2
c i.e. trans2=0 everywhere for rho1=1.0,rho2=1.0, 
c      trans2=1 everywhere for rho1=0.0,rho2=0.0
c----------------------------------------------------------------------
        real*8 function trans2(x,rho1,rho2)
        implicit none
        real*8 x,rho1,rho2

        real*8 xb

        ! initialize fixed-size variables 
        data xb/0.0/

        !--------------------------------------------------------------
 
        xb=(x-rho1)/(rho2-rho1)

        if (x.ge.rho2) then
          trans2=1
        else if (x.ge.rho1) then
          trans2=exp(-1/xb)/(exp(-1/xb)+exp(-1/(1-xb)))
        else
          trans2=0
        end if

        return
        end

c----------------------------------------------------------------------
c Evolution of Hb is split into two routines, hb_t_evo and hb_i_evo.
c
c parameters:
c
c gauge : integer controlling which scheme
c
c rho1,rho2: real constants, gauge dependent
c
c current schemes:
c
c Hb_t:
c
c gauge = 0 : fixed gauge
c gauge = 1 : AdS-asymptotics gauge
c
c Hb_i:
c
c gauge = 0 : fixed gauge
c gauge = 1 : AdS-asymptotics gauge
c
c----------------------------------------------------------------------
        subroutine hb_t_evo(res,
     &                    gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                    gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                    gb_ty_np1,gb_ty_n,gb_ty_nm1,
     &                    gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                    gb_xy_np1,gb_xy_n,gb_xy_nm1,
     &                    gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &                    psi_np1,psi_n,psi_nm1,
     &                    omega_np1,omega_n,omega_nm1,
     &                    Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &                    Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &                    Hb_y_np1,Hb_y_n,Hb_y_nm1,
     &                    phi1_np1,phi1_n,phi1_nm1,
     &                    L,x,y,dt,chr,ex,
     &                    phys_bdy,ghost_width,Nx,Ny,
     &                    Hb_t_0,Hb_x_0,Hb_y_0,
     &                    gauge,t_n,rho1,rho2,xi1,xi2,
     &                    a1,a2)

        implicit none
        integer Nx,Ny,gauge,phys_bdy(4),ghost_width(4)
        real*8 res(Nx,Ny),t_n,t_np1
        real*8 chr(Nx,Ny),ex,L
        real*8 x(Nx),y(Ny),dt,rho1,rho2,xi1,xi2
        real*8 a1,a2
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

        real*8 Hb_t_0(Nx,Ny),Hb_x_0(Nx,Ny),Hb_y_0(Nx,Ny)

        integer i,j

        real*8 F_t_np1,F_x_np1,F_y_np1

        real*8 Hb_t0,Hb_x0,Hb_y0

        real*8 x0,y0

        real*8 trans

        real*8 f0,g0
        real*8 dx,dy

        real*8 PI
        parameter (PI=3.141592653589793d0)

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

        !--------------------------------------------------------------

        dx=x(2)-x(1)
        dy=y(2)-y(1)

        t_np1=t_n+dt

        do i=1,Nx
          do j=1,Ny
            res(i,j)=0
            if (chr(i,j).eq.ex) then
              Hb_t_np1(i,j)=0
            end if
          end do
        end do

        !--------------------------------------------------------------
        ! gauge 0
        !--------------------------------------------------------------

        if (gauge.eq.0) return

        !--------------------------------------------------------------
        ! gauge 1
        !--------------------------------------------------------------

        if (gauge.eq.1) then
          do i=2,Nx-1
            do j=2,Ny-1
              if (chr(i,j).ne.ex) then
                  x0=x(i)
                  y0=y(j)
                  Hb_t0=Hb_t_0(i,j)

                  ! calculate gbar derivatives for np1 time level
                  call df2_int(gb_tt_np1,gb_tt_np1,gb_tt_np1,
     &                 gb_tt_t,gb_tt_x,gb_tt_y,
     &                 gb_tt_tt,gb_tt_tx,gb_tt_ty,
     &                 gb_tt_xx,gb_tt_xy,gb_tt_yy,
     &                 dx,dy,dt,i,j,chr,ex,Nx,Ny,'gb_tt')
                  call df2_int(gb_tx_np1,gb_tx_np1,gb_tx_np1,
     &                 gb_tx_t,gb_tx_x,gb_tx_y,
     &                 gb_tx_tt,gb_tx_tx,gb_tx_ty,
     &                 gb_tx_xx,gb_tx_xy,gb_tx_yy,
     &                 dx,dy,dt,i,j,chr,ex,Nx,Ny,'gb_tx')
                  call df2_int(gb_ty_np1,gb_ty_np1,gb_ty_np1,
     &                 gb_ty_t,gb_ty_x,gb_ty_y,
     &                 gb_ty_tt,gb_ty_tx,gb_ty_ty,
     &                 gb_ty_xx,gb_ty_xy,gb_ty_yy,
     &                 dx,dy,dt,i,j,chr,ex,Nx,Ny,'gb_ty')
                  call df2_int(gb_xx_np1,gb_xx_np1,gb_xx_np1,
     &                 gb_xx_t,gb_xx_x,gb_xx_y,
     &                 gb_xx_tt,gb_xx_tx,gb_xx_ty,
     &                 gb_xx_xx,gb_xx_xy,gb_xx_yy,
     &                 dx,dy,dt,i,j,chr,ex,Nx,Ny,'gb_xx')
                  call df2_int(gb_xy_np1,gb_xy_np1,gb_xy_np1,
     &                 gb_xy_t,gb_xy_x,gb_xy_y,
     &                 gb_xy_tt,gb_xy_tx,gb_xy_ty,
     &                 gb_xy_xx,gb_xy_xy,gb_xy_yy,
     &                 dx,dy,dt,i,j,chr,ex,Nx,Ny,'gb_xy')
                  call df2_int(gb_yy_np1,gb_yy_np1,gb_yy_np1,
     &                 gb_yy_t,gb_yy_x,gb_yy_y,
     &                 gb_yy_tt,gb_yy_tx,gb_yy_ty,
     &                 gb_yy_xx,gb_yy_xy,gb_yy_yy,
     &                 dx,dy,dt,i,j,chr,ex,Nx,Ny,'gb_yy')
                  call df2_int(psi_np1,psi_np1,psi_np1,
     &                 psi_t,psi_x,psi_y,
     &                 psi_tt,psi_tx,psi_ty,
     &                 psi_xx,psi_xy,psi_yy,
     &                 dx,dy,dt,i,j,chr,ex,Nx,Ny,'psi')
                  call df2_int(omega_np1,omega_np1,omega_np1,
     &                 omega_t,omega_x,omega_y,
     &                 omega_tt,omega_tx,omega_ty,
     &                 omega_xx,omega_xy,omega_yy,
     &                 dx,dy,dt,i,j,chr,ex,Nx,Ny,'omega')
                  call df2_int(phi1_np1,phi1_np1,phi1_np1,
     &                 phi1_t,phi1_x,phi1_y,
     &                 phi1_tt,phi1_tx,phi1_ty,
     &                 phi1_xx,phi1_xy,phi1_yy,
     &                 dx,dy,dt,i,j,chr,ex,Nx,Ny,'phi1')

                  F_t_np1=gb_tx_np1(i,j)*2.0d0
     &                   +gb_tx_yy/10.0d0/PI**2*a1
     &                   +gb_tx_y  
     &                   *cos(PI*y0)/sin(PI*y0)
     &                   *2.0d0/5.0d0/PI*a2 !nonzero a2 induces instability at poles 

                  f0=trans(x0,rho1,rho2)
                  g0=(t_np1/(xi2*f0+xi1*(1-f0)))**4

                  if (xi2.le.1e-16) then
                    Hb_t_np1(i,j)=F_t_np1
                  else
                    Hb_t_np1(i,j)=F_t_np1+(Hb_t0-F_t_np1)*exp(-g0)
                  end if
              end if
            end do
          end do
          return
        end if

        !--------------------------------------------------------------
        ! otherwise
        !--------------------------------------------------------------

        write(*,*) 'hb_t_evo : error, gauge,',gauge,' unknown'

        return
        end

c-----------------------------------------------------------------------
        subroutine hb_i_evo(res,
     &                    gb_tt_np1,gb_tt_n,gb_tt_nm1,
     &                    gb_tx_np1,gb_tx_n,gb_tx_nm1,
     &                    gb_ty_np1,gb_ty_n,gb_ty_nm1,
     &                    gb_xx_np1,gb_xx_n,gb_xx_nm1,
     &                    gb_xy_np1,gb_xy_n,gb_xy_nm1,
     &                    gb_yy_np1,gb_yy_n,gb_yy_nm1,
     &                    psi_np1,psi_n,psi_nm1,
     &                    omega_np1,omega_n,omega_nm1,
     &                    Hb_t_np1,Hb_t_n,Hb_t_nm1,
     &                    Hb_x_np1,Hb_x_n,Hb_x_nm1,
     &                    Hb_y_np1,Hb_y_n,Hb_y_nm1,
     &                    phi1_np1,phi1_n,phi1_nm1,
     &                    L,x,y,dt,chr,ex,
     &                    phys_bdy,ghost_width,Nx,Ny,
     &                    Hb_t_0,Hb_x_0,Hb_y_0,
     &                    gauge,t_n,rho1,rho2,xi1,xi2,
     &                    b1,b2,b3,b4,b5,b6,b7,b8,b9,
     &                    b10,b11,b12,
     &                    c1,c2,c3,c4,c5,c6,c7,c8,c9,
     &                    c10,c11,c12,c13)

        implicit none
        integer Nx,Ny,gauge,phys_bdy(4),ghost_width(4)
        real*8 res(Nx,Ny),t_n,t_np1
        real*8 chr(Nx,Ny),ex,L
        real*8 x(Nx),y(Ny),dt,rho1,rho2,xi1,xi2
        real*8 b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12
        real*8 c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13
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

        real*8 Hb_t_0(Nx,Ny),Hb_x_0(Nx,Ny),Hb_y_0(Nx,Ny)

        integer i,j

        real*8 F_t_np1,F_x_np1,F_y_np1

        real*8 Hb_t0,Hb_x0,Hb_y0

        real*8 x0,y0

        real*8 trans

        real*8 f0,g0
        real*8 dx,dy

        real*8 PI
        parameter (PI=3.141592653589793d0) 

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

        !--------------------------------------------------------------

        dx=x(2)-x(1)
        dy=y(2)-y(1)

        t_np1=t_n+dt

        do i=1,Nx
          do j=1,Ny
            res(i,j)=0
            if (chr(i,j).eq.ex) then
              Hb_x_np1(i,j)=0
            end if
          end do
        end do

        !--------------------------------------------------------------
        ! gauge 0
        !--------------------------------------------------------------

        if (gauge.eq.0) return

        !--------------------------------------------------------------
        ! gauge 1
        !--------------------------------------------------------------

        if (gauge.eq.1) then 
          do i=2,Nx-1
            do j=2,Ny-1
              if (chr(i,j).ne.ex) then
                  x0=x(i)
                  y0=y(j)
                  Hb_x0=Hb_x_0(i,j)
                  Hb_y0=Hb_y_0(i,j)

                  ! calculate gbar derivatives for np1 time level
                  call df2_int(gb_tt_np1,gb_tt_np1,gb_tt_np1,
     &                 gb_tt_t,gb_tt_x,gb_tt_y,
     &                 gb_tt_tt,gb_tt_tx,gb_tt_ty,
     &                 gb_tt_xx,gb_tt_xy,gb_tt_yy,
     &                 dx,dy,dt,i,j,chr,ex,Nx,Ny,'gb_tt')
                  call df2_int(gb_tx_np1,gb_tx_np1,gb_tx_np1,
     &                 gb_tx_t,gb_tx_x,gb_tx_y,
     &                 gb_tx_tt,gb_tx_tx,gb_tx_ty,
     &                 gb_tx_xx,gb_tx_xy,gb_tx_yy,
     &                 dx,dy,dt,i,j,chr,ex,Nx,Ny,'gb_tx')
                  call df2_int(gb_ty_np1,gb_ty_np1,gb_ty_np1,
     &                 gb_ty_t,gb_ty_x,gb_ty_y,
     &                 gb_ty_tt,gb_ty_tx,gb_ty_ty,
     &                 gb_ty_xx,gb_ty_xy,gb_ty_yy,
     &                 dx,dy,dt,i,j,chr,ex,Nx,Ny,'gb_ty')
                  call df2_int(gb_xx_np1,gb_xx_np1,gb_xx_np1,
     &                 gb_xx_t,gb_xx_x,gb_xx_y,
     &                 gb_xx_tt,gb_xx_tx,gb_xx_ty,
     &                 gb_xx_xx,gb_xx_xy,gb_xx_yy,
     &                 dx,dy,dt,i,j,chr,ex,Nx,Ny,'gb_xx')
                  call df2_int(gb_xy_np1,gb_xy_np1,gb_xy_np1,
     &                 gb_xy_t,gb_xy_x,gb_xy_y,
     &                 gb_xy_tt,gb_xy_tx,gb_xy_ty,
     &                 gb_xy_xx,gb_xy_xy,gb_xy_yy,
     &                 dx,dy,dt,i,j,chr,ex,Nx,Ny,'gb_xy')
                  call df2_int(gb_yy_np1,gb_yy_np1,gb_yy_np1,
     &                 gb_yy_t,gb_yy_x,gb_yy_y,
     &                 gb_yy_tt,gb_yy_tx,gb_yy_ty,
     &                 gb_yy_xx,gb_yy_xy,gb_yy_yy,
     &                 dx,dy,dt,i,j,chr,ex,Nx,Ny,'gb_yy')
                  call df2_int(psi_np1,psi_np1,psi_np1,
     &                 psi_t,psi_x,psi_y,
     &                 psi_tt,psi_tx,psi_ty,
     &                 psi_xx,psi_xy,psi_yy,
     &                 dx,dy,dt,i,j,chr,ex,Nx,Ny,'psi')
                  call df2_int(omega_np1,omega_np1,omega_np1,
     &                 omega_t,omega_x,omega_y,
     &                 omega_tt,omega_tx,omega_ty,
     &                 omega_xx,omega_xy,omega_yy,
     &                 dx,dy,dt,i,j,chr,ex,Nx,Ny,'omega')
                  call df2_int(phi1_np1,phi1_np1,phi1_np1,
     &                 phi1_t,phi1_x,phi1_y,
     &                 phi1_tt,phi1_tx,phi1_ty,
     &                 phi1_xx,phi1_xy,phi1_yy,
     &                 dx,dy,dt,i,j,chr,ex,Nx,Ny,'phi1')

                  F_x_np1=gb_xx_np1(i,j)*2.0d0
     &                   +(gb_tt_yy/4.0d0/PI**2
     &                    +gb_tt_y*cos(PI*y0)/sin(PI*y0)/PI)*b1
     &                   +(gb_xx_yy/16.0d0/PI**2
     &                    +gb_xx_y*cos(PI*y0)/sin(PI*y0)/4.0d0/PI)*b2
     &                   +(-psi_yy/4.0d0/PI**2
     &                     -psi_y*cos(PI*y0)/sin(PI*y0)/PI)*b3
     &                   +(4.0d0*sin(PI*y0)*phi1_y/PI
     &                    +20.0d0*cos(PI*y0)*phi1_np1(i,j))*b4
     &                   +(-32.0d0*omega_np1(i,j))*b5
     &                   +(-8.0d0*gb_yy_np1(i,j)*PI**2)*b6
     &                   +(4.0d0*cos(PI*y0)/sin(PI*y0)
     &                     *gb_xy_np1(i,j)*PI)*b7
     &                   +(-sin(PI*y0)*phi1_y/PI
     &                     -5.0d0*cos(PI*y0)*phi1_np1(i,j))*b8
     &                   +(-1.5d0*psi_np1(i,j))*b9
     &                   +(8.0d0*omega_np1(i,j))*b10
     &                   +(2.0d0*gb_yy_np1(i,j)*PI**2)*b11
     &                   +(0.5d0*gb_tt_np1(i,j))*b12
                  F_y_np1=gb_xy_np1(i,j)*1.5d0
     &                   +gb_xy_yy/8.0d0/PI**2*c1
     &                   +(-omega_yy*tan(PI*y0)/2.0d0/PI**2
     &                     -omega_y*2.0d0/PI)/PI*c2
     &                   +(-sin(PI*y0)*phi1_yy/4.0d0/PI**2
     &                     -cos(PI*y0)*phi1_y*1.5d0/PI)/PI*c3
     &                   +(4.0d0*omega_y/PI**2)*c4
     &                   +(1.25d0*sin(PI*y0)*phi1_np1(i,j)/PI)*c5
     &                   +(0.375d0*psi_y/PI**2)*c6
     &                   +(2.0d0*omega_y/PI**2)*c7
     &                   +(-0.5d0*gb_yy_y)*c8
     &                   +(0.125d0*gb_xx_y/PI**2)*c9
     &                   +(-0.125d0*gb_tt_y/PI**2)*c10
     &                   +(6.0d0/sin(2.0d0*PI*y0)/PI*
     &                    (gb_yy_np1(i,j)/PI**2-omega_np1(i,j)))*c11
     &                   +(-2.0d0*sin(PI*y0)/PI
     &                     *(5.0d0*phi1_np1(i,j)
     &                      +tan(PI*y0)*phi1_y/PI))*c12
     &                   +(16.0d0/PI*tan(PI*y0)*omega_np1(i,j))*c13
                  f0=trans(x0,rho1,rho2)
                  g0=(t_np1/(xi2*f0+xi1*(1-f0)))**4

                  if (xi2.le.1e-16) then
                    Hb_x_np1(i,j)=F_x_np1
                    Hb_y_np1(i,j)=F_y_np1
                  else
                    Hb_x_np1(i,j)=F_x_np1+(Hb_x0-F_x_np1)*exp(-g0)
                    Hb_y_np1(i,j)=F_y_np1+(Hb_y0-F_y_np1)*exp(-g0)
                  end if
              end if
            end do
          end do
          return
        end if

        !--------------------------------------------------------------
        ! otherwise
        !--------------------------------------------------------------

        write(*,*) 'hb_i_evo : error, gauge,',gauge,' unknown'

        return
        end

c-----------------------------------------------------------------------
