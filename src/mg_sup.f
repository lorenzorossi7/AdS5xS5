c----------------------------------------------------------------------
c in polar coordinates t==t, x==rho
c
c Relaxation, Lop, Residual routines for MG solver for zetab
c satisfying L.zetab=0, with induced MG right-hand side R giving L.zetab=R
c
c (fixed using correct g=gads+(1-rho^2)gb factorization)
c----------------------------------------------------------------------

c-----------------------------------------------------------------------
c if 
c 
c action = 1 : performs 1 RB relaxation sweep of zetab
c              (_lop, _res vars unused) ... sets norm
c action = 2 : calculates MG residual L.zetab-R --> _res vars ... sets norm
c              (_lop unused)
c action = 3 : calculates normal residual L.zetab --> _lop
c              (_res, _rhs unused)

c-----------------------------------------------------------------------      
        subroutine mg_sup(action,zetab,zetab_rhs,zetab_lop,zetab_res,
     &                    phi1,L,cmask,phys_bdy,chr,ex,x,norm,Nx)
        implicit none
        integer Nx,action
        integer phys_bdy(2)
        real*8 zetab(Nx),zetab_rhs(Nx),zetab_lop(Nx)
        real*8 zetab_res(Nx)
        real*8 phi1(Nx)
        real*8 cmask(Nx),chr(Nx)
        real*8 x(Nx),norm,ex,L

        real*8 trhoE_grad,trhoE_ptl
        real*8 ddzeta,ddzeta_Jac,grad_zetab_sq
        real*8 ddphi1,ddphi1_Jac,grad_phi1_sq

        real*8 phi10(Nx)

        real*8 g_ll(5,5)

        integer relax,lop,residual,cdiff_method
        parameter (relax=1,lop=3,residual=2)
        real*8 x0
        real*8 dx

        real*8 zetab0

        real*8 phi10_0

        real*8 Jac,res,rhs,new_rhs
        real*8 lambda5
        real*8 PI
        parameter (PI=3.141592653589793d0)

        integer i,pass,sum
        integer is

        logical first,do_red_black
        data first/.true./
!        data do_red_black/.false./
        data do_red_black/.true./
        save first

        ! initialize fixed-size variables
        data i,pass,sum/0,0,0/
        data is/0/

        data x0/0.0/
        data dx/0.0/

        data zetab0/0.0/

        data phi10_0/0.0/

        data Jac,res,rhs,new_rhs/0.0,0.0,0.0,0.0/
        data lambda5/0.0/

        !--------------------------------------------------------------

        first=.false.

        norm=0
        sum=0

        dx=(x(2)-x(1))

        rhs=0

        ! sets AdS5D cosmological constant
        lambda5=-6/L/L

        ! manually reconstruct phi10=phi1*(1-x^2)^3 
        do i=1,Nx
          x0=x(i)
          if (phi1(i).ne.0) phi10(i)=phi1(i)*(1-x0**2)**3
          if (phi1(i).eq.0) phi10(i)=0
        end do

        ! (REGION) interior, solve L.zetab=0 Hamiltonian constraint 
        do pass=0,1
          do i=2,Nx-1
            x0=x(i)

            if (
     &          ((do_red_black.and.mod(i+pass,2).eq.0)
     &           .or.(.not.do_red_black.and.pass.eq.0)) 
     &           .and.(cmask(i).eq.1)
     &           .and.(chr(i).ne.ex)      
     &         ) then

              ! fill in zetab 
              zetab0=zetab(i)

              ! fill in phi10_0
              phi10_0=phi10(i)

              ! computes initial energy density at i, with initial data
              ! time-symmetric so phi_t=0, and free scalar so V(phi)=0
              call df_int(phi10,ddphi1,ddphi1_Jac,
     &                    grad_phi1_sq,
     &                    x,i,chr,L,ex,Nx)
              trhoE_grad=grad_phi1_sq/2
              trhoE_ptl=0

              ! computes normal residual L.zetab
              !(NOTE: the physical energy density rhoE is such that
              ! rhoE_grad=trhoE_grad*zeta^(-2)
              ! where zeta=1+(1-x0**2)**3*zetab, and
              ! rhoE_ptl=trhoE_ptl)
              call df_int(zetab,ddzeta,ddzeta_Jac,
     &                    grad_zetab_sq,
     &                    x,i,chr,L,ex,Nx)
              res=ddzeta
     &            -lambda5*(1+(1-x0**2)**3*zetab0)/3
     &            +(lambda5+8*PI*trhoE_ptl)
     &             *((1+(1-x0**2)**3*zetab0)**3)/3
     &            +8*PI*trhoE_grad*(1+(1-x0**2)**3*zetab0)/3

              ! computes MG residual L.zetab-R
              rhs=res-zetab_rhs(i)

              Jac=ddzeta_Jac
     &            -( lambda5*( (1-x0**2)**3 )/3 )
     &            +(lambda5+8*PI*trhoE_ptl)
     &             *( (1+(1-x0**2)**3*zetab0)**2 )
     &             *(1-x0**2)**3
     &            +8*PI*trhoE_grad*( (1-x0**2)**3 )/3

              ! performs action
              if (action.eq.residual) then
                zetab_res(i)=rhs
              else if (action.eq.lop) then
                zetab_lop(i)=res
              else if (action.eq.relax) then
                zetab(i)=zetab(i)-rhs/Jac
              end if

              norm=norm+rhs**2
              sum=sum+1
              
            end if
          end do
        end do

        ! (REGION) x=0 axis, solve d/dx(zetab)=0  
        if (phys_bdy(1).ne.0) then

          ! computes normal residual d/dy(zetab)
          res=zetab(1)-(4*zetab(2)-zetab(3))/3

          ! computes MG residual d/dy(zetab)-R
          rhs=res-zetab_rhs(1)

          ! computes diag. Jacobian of zetab->L.zetab transformation
          ! by differentiating L.zetab wrt. z(1) diag. entries
          Jac=1

          if (action.eq.residual) then
            zetab_res(1)=rhs
          else if (action.eq.lop) then
            zetab_lop(1)=res
          else if (action.eq.relax) then
            zetab(1)=zetab(1)-rhs/Jac                  
          end if

          norm=norm+rhs**2
          sum=sum+1

        end if

        norm=sqrt(norm/sum)

        return
        end

c-----------------------------------------------------------------------
c in polar coordinates t==t, x==rho
c
c The following initializes the rest of the metric and Hb,
c given zetab
c
c NOTE: if we ever add gb_xy,gb_xz,gb_yz, must define them
c       as AMRD MG_cnst vars
c-----------------------------------------------------------------------
        subroutine init_ghb(zetab,phi1,gb_tt,gb_tx,gb_xx,psi,rhoa,rhob,
     &                      L,phys_bdy,chr,ex,x,Nx)
        implicit none
        integer Nx
        integer phys_bdy(2)
        real*8 zetab(Nx)
        real*8 phi1(Nx)
        real*8 gb_tt(Nx),gb_tx(Nx),gb_xx(Nx),psi(Nx)
        real*8 chr(Nx),ex,L
        real*8 x(Nx)
        real*8 rhoa,rhob
 
        real*8 dx

        real*8 zetab0
        
        integer i
        real*8 x0
        real*8 PI
        parameter (PI=3.141592653589793d0)

        real*8 g0_tt_ads0,g0_xx_ads0,g0_psi_ads0

        real*8 xb,trans

        ! initialize fixed-size variables 
        data zetab0/0.0/

        data i/0/
        data x0/0.0/

        data g0_tt_ads0,g0_xx_ads0,g0_psi_ads0/0.0,0.0,0.0/

        data xb,trans/0.0,0.0/

        !--------------------------------------------------------------

        dx=(x(2)-x(1))

        ! initialize gbars given zetab, using full metric expression
        ! g0_ij=g0_ij_ads*zetab^2, and 
        ! g-bar expression g0_ij=g0_ij_ads+gb_ij*(1-x0^2)
        ! psi expression g0_psi=g0_psi_ads+psi*(1-x0^2)*x0^2*Sin(PI*y0)^2
        do i=2,Nx-1
 
          x0=x(i)

          if (chr(i).ne.ex) then 

            zetab0=zetab(i)
            g0_tt_ads0=-((1-x0)**2+x0**2)/(1-x0)**2
            g0_xx_ads0=1/((1-x0)**2+x0**2)/(1-x0)**2
            g0_psi_ads0=x0**2/(1-x0)**2

            gb_tt(i)=0
            gb_tx(i)=0
            gb_xx(i)=g0_xx_ads0
     &                 *((1+(1-x0**2)**3*zetab0)**2-1)
     &                 /(1-x0**2)
            psi(i)=g0_psi_ads0
     &               *((1+(1-x0**2)**3*zetab0)**2-1)
     &               /(1-x0**2)
     &               /(x0**2)
            xb=(rhob-x0)/(rhob-rhoa)
            if (x0.ge.rhob) then
              trans=1
            else if (x0.ge.rhoa) then
              trans=1-xb**3*(6*xb**2-15*xb+10)
            else
              trans=0
            end if
            gb_tt(i)=(gb_xx(i)+3*psi(i))*trans

          endif

        end do

        ! x=0 regularization of gbars
        call axi_reg_g(gb_tt,gb_tx,gb_xx,psi,chr,ex,L,x,Nx)

        return
        end

c----------------------------------------------------------------------
c in polar coordinates t==t, x==rho
c
c the following computes ddf, the background Laplacian on f, and 
c grad_f_sq, the background squared gradient of f; 
c evaluated at point i and at the initial time.
c----------------------------------------------------------------------
        subroutine df_int(f,ddf,ddf_Jac,grad_f_sq,
     &                    x,i,chr,L,ex,Nx)
        implicit none
        integer Nx
        integer i
        real*8 chr(Nx),ex
        real*8 x(Nx),L
        real*8 f(Nx),fdot(Nx)
        real*8 f0,f0_x(5),f0_xx(5,5),ddf,ddf_Jac,grad_f_sq

        real*8 PI
        parameter (PI=3.141592653589793d0)

        integer a,b,c,d

        real*8 dx
        real*8 x0

        ! initialize fixed-size variables 
        data a,b,c,d/0,0,0,0/

        data dx/0.0/
        data x0/0.0/
 
        data f0/0.0/
        data f0_x/5*0.0/
        data f0_xx/25*0.0/

        !--------------------------------------------------------------

        dx=(x(2)-x(1))

        ! sets x0 to values at i
        x0=x(i)

        ! set first and second derivatives
        !(only the spatial part of the metric here, so don't need time derivatives)
        f0=f(i)
        f0_x(2)=(f(i+1)-f(i-1))/2/dx
        f0_xx(2,2)=(f(i+1)-2*f(i)+f(i-1))/dx/dx

        ! calculate ddzeta, background Laplacian acting on zeta
        ! and ddzeta_Jac, Jacobian of zeta->ddzeta transformation 
        ! (DDf = g^ab D_a D_b f)
        ! all in terms of zetab
        ! assuming zeta=1+(1-x0**2)**3*zetab
        ddf= 
     &       ( ((-1 + x0)*(2*(-2 + x0)*x0**2 + 
     &         L**2*(-1 + x0)**2*(-3 + 2*x0)))/(L**2*x0)  
     &       )*( -2*x0*3*(1-x0**2)**2*f0
     &           +(1-x0**2)**3*f0_x(2) )
     &      +( ((-1 + x0)**2*(L**2*(-1 + x0)**2 + x0**2))/L**2
     &       )*( (
     &            -2*3*(1-x0**2)**2
     &            +4*x0**2*3*2
     &             *(1-x0**2)
     &           )*f0
     &          -2*2*x0*3*(1-x0**2)**2*f0_x(2)
     &          +(1-x0**2)**3*f0_xx(2,2) )

        ddf_Jac=
     &           ( ((-1 + x0)*(2*(-2 + x0)*x0**2 + 
     &             L**2*(-1 + x0)**2*(-3 + 2*x0)))/(L**2*x0)  
     &           )*( -2*x0*3*(1-x0**2)**2 )
     &          +( ((-1 + x0)**2*(L**2*(-1 + x0)**2 + x0**2))/L**2
     &            )*( (
     &                 -2*3*(1-x0**2)**2
     &                 +4*x0**2*3*2
     &                  *(1-x0**2)
     &                )
     &               +(1-x0**2)**3*(-2/dx/dx) )

        ! calculate grad_f_sq, squared gradient of f
        !(Df^2 = g^ab D_a f D_b f)
        grad_f_sq=
     &            ( (-1 + x0)**2*(1 + x0*(-2 + x0 + x0/L**2))
     &            )*f0_x(2)*f0_x(2)

        return
        
        end
