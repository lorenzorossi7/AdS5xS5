c----------------------------------------------------------------------
c in polar coordinates t,x for x in [0,1]
c
c applies regularity conditions for the scalar field
c----------------------------------------------------------------------

        subroutine axi_reg_phi(phi,chr,ex,L,x,y,Nx,Ny)
        implicit none
        integer Nx,Ny
        real*8 phi(Nx,Ny)
        real*8 chr(Nx,Ny),ex,L
        real*8 x(Nx),y(Ny)

        integer i,j
        real*8 PI
        real*8 dx,dy
        parameter (PI=3.141592653589793d0)

        dx=x(2)-x(1)
        dy=y(2)-y(1)

        if (x(1).gt.dx/2.and.y(1).gt.dy/2.and.y(Ny).lt.1-dy/2) return

        if (x(1).le.dx/2) then
          do j=1,Ny
            if (chr(1,j).ne.ex) then
              if (chr(2,j).ne.ex.and.chr(3,j).ne.ex) then
                phi(1,j)=(4*phi(2,j)-phi(3,j))/3
              else
                write(*,*) 'WARNING axi_reg_phi'
              end if
            else
              phi(1,j)=0
            end if
          end do
        end if

        if (y(1).le.dy/2) then
          do i=1,Nx
            if (chr(i,1).ne.ex) then
              if (chr(i,2).ne.ex.and.chr(i,3).ne.ex) then
                phi(i,1)=(4*phi(i,2)-phi(i,3))/3
              else
                write(*,*) 'WARNING axi_reg_phi'
              end if
            else
              phi(i,1)=0
            end if
          end do
        end if

        if (y(Ny).ge.1-dy/2) then
          do i=1,Nx
            if (chr(i,Ny).ne.ex) then
              if (chr(i,Ny-1).ne.ex.and.chr(i,Ny-2).ne.ex) then
                phi(i,Ny)=(4*phi(i,Ny-1)-phi(i,Ny-2))/3
              else
                write(*,*) 'WARNING axi_reg_phi'
              end if
            else
              phi(i,Ny)=0
            end if
          end do
        end if

        return
        end

c----------------------------------------------------------------------
c in polar coordinates t,x for x in [0,1]
c
c applies regularity conditions for the metric
c----------------------------------------------------------------------

        subroutine axi_reg_g(gb_tt,gb_tx,gb_ty,gb_xx,gb_xy,gb_yy,
     &                       psi,omega,chr,ex,L,x,y,Nx,Ny)
        implicit none
        integer Nx,Ny
        real*8 gb_tt(Nx,Ny),gb_tx(Nx,Ny),gb_ty(Nx,Ny),gb_xx(Nx,Ny)
        real*8 gb_xy(Nx,Ny),gb_yy(Nx,Ny),psi(Nx,Ny),omega(Nx,Ny)
        real*8 chr(Nx,Ny),ex,L
        real*8 x(Nx),y(Ny)

        integer i,j
        real*8 PI,dx,dy
        parameter (PI=3.141592653589793d0)

        dx=x(2)-x(1)
        dy=y(2)-y(1)

        if (x(1).gt.dx/2.and.y(1).gt.dy/2.and.y(Ny).lt.1-dy/2) return

        if (x(1).le.dx/2) then
          do j=1,Ny
            if (chr(1,j).ne.ex) then
              if (chr(2,j).ne.ex.and.chr(3,j).ne.ex) then
                gb_tt(1,j)=(4*gb_tt(2,j)-gb_tt(3,j))/3
                gb_ty(1,j)=(4*gb_ty(2,j)-gb_ty(3,j))/3
                gb_xx(1,j)=(4*gb_xx(2,j)-gb_xx(3,j))/3
                gb_xy(1,j)=(4*gb_xy(2,j)-gb_xy(3,j))/3
                gb_yy(1,j)=(4*gb_yy(2,j)-gb_yy(3,j))/3
                psi(1,j)=(4*psi(2,j)-psi(3,j))/3
                omega(1,j)=(4*omega(2,j)-omega(3,j))/3
                gb_tx(1,j)=0
              else
                write(*,*) 'WARNING axi_reg_g'
              end if
            else
              gb_tt(1,j)=0
              gb_ty(1,j)=0
              gb_xx(1,j)=0
              gb_xy(1,j)=0
              gb_yy(1,j)=0
              gb_tx(1,j)=0
              psi(1,j)=0
              omega(1,j)=0
            end if
          end do
        end if

        if (y(1).le.dy/2) then
          do i=1,Nx
            if (chr(i,1).ne.ex) then
              if (chr(i,2).ne.ex.and.chr(i,3).ne.ex)
     &        then
                gb_tt(i,1)=(4*gb_tt(i,2)-gb_tt(i,3))/3
                gb_xx(i,1)=(4*gb_xx(i,2)-gb_xx(i,3))/3
                gb_yy(i,1)=(4*gb_yy(i,2)-gb_yy(i,3))/3
                gb_tx(i,1)=(4*gb_tx(i,2)-gb_tx(i,3))/3
                psi(i,1)=(4*psi(i,2)-psi(i,3))/3
                omega(i,1)=(4*omega(i,2)-omega(i,3))/3
                gb_ty(i,1)=0
                gb_xy(i,1)=0                
              else
                write(*,*) 'WARNING axi_reg_g'
              end if
            else
              gb_tt(i,1)=0
              gb_xx(i,1)=0
              gb_yy(i,1)=0
              gb_tx(i,1)=0
              gb_ty(i,1)=0
              gb_xy(i,1)=0
              psi(i,1)=0
              omega(i,1)=0
            end if
          end do
        end if

        if (y(Ny).ge.1-dy/2) then
          do i=1,Nx
            if (chr(i,Ny).ne.ex) then
              if (chr(i,Ny-1).ne.ex.and.chr(i,Ny-2).ne.ex)
     &        then
                gb_tt(i,Ny)=(4*gb_tt(i,Ny-1)-gb_tt(i,Ny-2))/3
                gb_xx(i,Ny)=(4*gb_xx(i,Ny-1)-gb_xx(i,Ny-2))/3
                gb_yy(i,Ny)=(4*gb_yy(i,Ny-1)-gb_yy(i,Ny-2))/3
                gb_tx(i,Ny)=(4*gb_tx(i,Ny-1)-gb_tx(i,Ny-2))/3
                gb_ty(i,Ny)=0
                gb_xy(i,Ny)=0
                psi(i,Ny)=0
                omega(i,Ny)=0
              else
                !write(*,*) 'WARNING axi_reg_g extras'
              end if
            else
              gb_tt(i,Ny)=0
              gb_xx(i,Ny)=0
              gb_yy(i,Ny)=0
              gb_tx(i,Ny)=0
              gb_ty(i,Ny)=0
              gb_xy(i,Ny)=0
              psi(i,Ny)=0
              omega(i,Ny)=0
            end if
          end do
        end if

        return
        end

c----------------------------------------------------------------------
c in polar coordinates t,x for x in [0,1]
c
c applies regularity conditions for the metric
c----------------------------------------------------------------------

        subroutine axi_reg_Hb(Hb_t,Hb_x,chr,ex,L,x,y,Nx,Ny)
        implicit none
        integer Nx,Ny
        real*8 Hb_t(Nx,Ny),Hb_x(Nx,Ny)
        real*8 chr(Nx,Ny),ex,L
        real*8 x(Nx),y(Ny)

        integer i,j
        real*8 PI,dx,dy
        parameter (PI=3.141592653589793d0)

        dx=x(2)-x(1)
        dy=y(2)-y(1)

        if (x(1).gt.dx/2.and.y(1).gt.dy/2.and.y(Ny).lt.1-dy/2) return

        if (x(1).le.dx/2) then
          do j=1,Ny
            if (chr(1,j).ne.ex) then
              if (chr(2,j).ne.ex.and.chr(3,j).ne.ex) then
                Hb_t(1,j)=(4*Hb_t(2,j)-Hb_t(3,j))/3
                Hb_x(1,j)=0
              else
                write(*,*) 'WARNING axi_reg_Hb'
              end if
            else
              Hb_t(1,j)=0
              Hb_x(1,j)=0
            end if
          end do
        end if

        if (y(1).le.dy/2) then
          do i=1,Nx
            if (chr(i,1).ne.ex) then
              if (chr(i,2).ne.ex.and.chr(i,3).ne.ex) then
                Hb_t(i,1)=(4*Hb_t(i,2)-Hb_t(i,3))/3
                Hb_x(i,1)=(4*Hb_x(i,2)-Hb_x(i,3))/3
              else
                write(*,*) 'WARNING axi_reg_Hb'
              end if
            else
              Hb_t(i,1)=0
              Hb_x(i,1)=0
            end if
          end do
        end if

        if (y(Ny).ge.1-dy/2) then
          do i=1,Nx
            if (chr(i,Ny).ne.ex) then
              if (chr(i,Ny-1).ne.ex.and.chr(i,Ny-2).ne.ex) then
                Hb_t(i,Ny)=(4*Hb_t(i,Ny-1)-Hb_t(i,Ny-2))/3
                Hb_x(i,Ny)=(4*Hb_x(i,Ny-1)-Hb_x(i,Ny-2))/3
              else
                write(*,*) 'WARNING axi_reg_Hb'
              end if
            else
              Hb_t(i,Ny)=0
              Hb_x(i,Ny)=0
            end if
          end do
        end if

        return
        end
