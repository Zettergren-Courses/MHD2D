module MHD

use phys_consts
use numerical

implicit none


!INTERESTINGLY A LARGE VALUE OF THIS PARAMETER TENDS TO MESS UP SITUATIONS WITH SPHERICAL SYMMETRY, ESPECIALLY AT LOWER RESOLUTIONS
!real(8) :: xi=5d0    !converging double shock.
real(8) :: xi=3d0  


contains



  subroutine MHD_advection(dt,dx1,dx1i,dx2,dx2i, &
                           rhom,rhou1,rhou2,rhou3,rhoeps, &
                           B1,B2,B3)

    !------------------------------------------------------------
    !-------PARTIALLY UPDATE MHD VARIABLES BY ADVECTING THEM BY
    !-------TIME STEP DT.  ALL ARRAYS SHOULD HAVE GHOST CELLS.
    !------------------------------------------------------------

    real(8), intent(in) :: dt
    real(8), dimension(0:), intent(in) :: dx1
    real(8), dimension(:), intent(in) :: dx1i
    real(8), dimension(0:), intent(in) :: dx2
    real(8), dimension(:), intent(in) :: dx2i
    real(8), dimension(-1:,-1:), intent(inout) :: rhom,rhou1,rhou2,rhou3,rhoeps,B1,B2,B3

    real(8), dimension(1:size(rhom,1)-3,1:size(rhom,2)-4) :: v1i,v1izeros
    real(8), dimension(1:size(rhom,1)-4,1:size(rhom,2)-3) :: v2i,v2izeros

    integer :: ix1,ix2,lx1,lx2


    !SIZES WITHOUT GHOST CELLS
    lx1=size(rhom,1)-4; lx2=size(rhom,2)-4;


    !COMPUTE CELL INTERFACE VELOCITIES
    !$omp parallel do
    do ix2=1,lx2
      v1i(:,ix2)=0.5d0*(rhou1(0:lx1,ix2)/rhom(0:lx1,ix2)+rhou1(1:lx1+1,ix2)/rhom(1:lx1+1,ix2))    !CELL INTERFACE SPEEDS (LEFT WALL OF ITH CELL)
    end do
    !$omp parallel do
    do ix1=1,lx1
      v2i(ix1,:)=0.5d0*(rhou2(ix1,0:lx2)/rhom(ix1,0:lx2)+rhou2(ix1,1:lx2+1)/rhom(ix1,1:lx2+1))    !CELL INTERFACE SPEEDS (LEFT WALL OF ITH CELL)
    end do


    !NULL VELOCITIES FOR ADVECTING THE MAGNETIC FIELD
    v1izeros=0d0; v2izeros=0d0;


    !ADVECTION SUBSTEP  (parallelized inside function)
    rhom=advec2D_MC(rhom,v1i,v2i,dt,dx1,dx1i,dx2,dx2i)
    rhou1=advec2D_MC(rhou1,v1i,v2i,dt,dx1,dx1i,dx2,dx2i)
    rhou2=advec2D_MC(rhou2,v1i,v2i,dt,dx1,dx1i,dx2,dx2i)
    rhou3=advec2D_MC(rhou3,v1i,v2i,dt,dx1,dx1i,dx2,dx2i)
    rhoeps=advec2D_MC(rhoeps,v1i,v2i,dt,dx1,dx1i,dx2,dx2i)
    B1=advec2D_MC(B1,v1izeros,v2i,dt,dx1,dx1i,dx2,dx2i)
    B2=advec2D_MC(B2,v1i,v2izeros,dt,dx1,dx1i,dx2,dx2i)
    B3=advec2D_MC(B3,v1i,v2i,dt,dx1,dx1i,dx2,dx2i)

  end subroutine MHD_advection



  subroutine MHD_sources(dt,dx1,dx2,x1,x2,du1full,du2full, &
                         rhom,rhou1,rhou2,rhou3,rhoeps, &
                         B1,B2,B3, &
                         u1,u2,u3)

    !------------------------------------------------------------
    !-------PARTIALLY UPDATE MHD VARIABLES BY ADDING IN SOURCE
    !-------TERMS.  ALL STATE VARIABLES SHOULD INCLUDE GHOST
    !-------CELLS.
    !------------------------------------------------------------

    real(8), intent(in) :: dt
    real(8), dimension(0:), intent(in) :: dx1
    real(8), dimension(0:), intent(in) :: dx2
    real(8), dimension(-1:), intent(in) :: x1
    real(8), dimension(-1:), intent(in) :: x2
    real(8), dimension(:,:), intent(in) :: du1full,du2full
    real(8), dimension(-1:,-1:), intent(inout) :: rhom,rhou1,rhou2,rhou3,rhoeps,B1,B2,B3,u1,u2,u3

    real(8), dimension(1:size(rhom,1)-4,1:size(rhom,2)-4) :: Q11,Q22,p11,p22
    real(8), dimension(1:size(rhom,1)-4,1:size(rhom,2)-4) :: p,grad1B1u3,rhoepshalf,sourceterm
    real(8), dimension(1:size(rhom,1)-4,1:size(rhom,2)-4) :: deriv1,deriv2,deriv3

    integer :: ix1,ix2,lx1,lx2


    !SIZES WITHOUT GHOST CELLS
    lx1=size(rhom,1)-4; lx2=size(rhom,2)-4;


    !VON NEUMANN-RICHTMYER ARTIFICIAL VISCOSITY ('CLEANS UP' SOLUTIONS
    !WITH STRONG SHOCKS) [POTTER, 1970; STONE ET AL 1992]
    !THIS COMPUTES THE EQUIVALENT STRESS TENSOR TERMS.
    !$omp parallel do
    do ix2=1,lx2
      Q11(:,ix2)=0.25d0*xi**2*(min(du1full(:,ix2),0d0))**2*rhom(1:lx1,ix2)
      Q22(:,ix2)=0.25d0*xi**2*(min(du2full(:,ix2),0d0))**2*rhom(1:lx1,ix2)
    end do


    !SOURCE TERMS FOR 1-COMP OF MOMENTUM
    !$omp parallel do
    do ix2=1,lx2
      p11(:,ix2)=rhoeps(1:lx1,ix2)*(gammap-1d0) &
         +(B1(1:lx1,ix2)**2+B2(1:lx1,ix2)**2+B3(1:lx1,ix2)**2)/2d0/mu0 &
         +Q11(:,ix2)   !x1-direction pressure/stress
    end do
    deriv1=derivative2D1(p11,dx1(1:lx1))
!    deriv1=derivative2D1obj(p11,dx1(1:lx1),x1,x2)
    deriv2=derivative2D1(B1(1:lx1,1:lx2),dx1(1:lx1))
!    deriv2=derivative2D1obj(B1(1:lx1,1:lx2),dx1(1:lx1),x1,x2)
    deriv3=derivative2D2(B1(1:lx1,1:lx2),dx2(1:lx2))
!    deriv3=derivative2D2obj(B1(1:lx1,1:lx2),dx2(1:lx2),x1,x2)
    !$omp parallel do
    do ix2=1,lx2
      rhou1(1:lx1,ix2)=rhou1(1:lx1,ix2)+dt*((-1d0*deriv1(:,ix2))+ &    !pressure terms
                        1d0/mu0*(B1(1:lx1,ix2)*deriv2(:,ix2)+ &    !magnetic tension from dB1/dx1
                        B2(1:lx1,ix2)*deriv3(:,ix2))+ &    !magnetic tension from dB1/dx2
                        rhom(1:lx1,ix2)*g1(1:lx1,ix2))     !gravity
    end do


    !SOURCE TERMS FOR 2-COMP OF MOMENTUM
    !$omp parallel do
    do ix2=1,lx2
      p22(:,ix2)=rhoeps(1:lx1,ix2)*(gammap-1d0) &
         +(B1(1:lx1,ix2)**2+B2(1:lx1,ix2)**2+B3(1:lx1,ix2)**2)/2d0/mu0 &
         +Q22(:,ix2)   !x2-direction pressure/stress
    end do
    deriv1=derivative2D2(p22,dx2(1:lx2))
!    deriv1=derivative2D2obj(p22,dx2(1:lx1),x1,x2)
    deriv2=derivative2D1(B2(1:lx1,1:lx2),dx1(1:lx1))
!    deriv2=derivative2D1obj(B2(1:lx1,1:lx2),dx1(1:lx1),x1,x2)
    deriv3=derivative2D2(B2(1:lx1,1:lx2),dx2(1:lx2))
!    deriv3=derivative2D2obj(B2(1:lx1,1:lx2),dx2(1:lx2),x1,x2)
    !$omp parallel do
    do ix2=1,lx2
      rhou2(1:lx1,ix2)=rhou2(1:lx1,ix2)+dt*(-1d0*deriv1(:,ix2) &    !pressure terms
                       +1d0/mu0*(B1(1:lx1,ix2)*deriv2(:,ix2)+ &    !magnetic tension from dB2/dx1
                       B2(1:lx1,ix2)*deriv3(:,ix2))+ &    !magnetic tension from dB2/dx2
                       rhom(1:lx1,ix2)*g2(1:lx1,ix2))      !gravity
    end do
    
    
    !SOURCE TERMS FOR 3-COMP OF MOMENTUM
    deriv1=derivative2D1(B3(1:lx1,1:lx2),dx1(1:lx1))
    deriv2=derivative2D2(B3(1:lx1,1:lx2),dx2(1:lx2))
    !$omp parallel do
    do ix2=1,lx2
      rhou3(1:lx1,ix2)=rhou3(1:lx1,ix2)+ &
                         dt*(1d0/mu0)*(B1(1:lx1,ix2)*deriv1(:,ix2)+ &    !magnetic tension from dB3/dx1
                         B2(1:lx1,ix2)*deriv2(:,ix2))    !magnetic tension from dB3/dx2
    end do
  
  
    !PARTIALLY UPDATED FLOWS FROM MOMENTA
    !$omp parallel do
    do ix2=1,lx2
      u1(:,ix2)=rhou1(:,ix2)/rhom(:,ix2)
      u2(:,ix2)=rhou2(:,ix2)/rhom(:,ix2)
      u3(:,ix2)=rhou3(:,ix2)/rhom(:,ix2)
    end do
  
  
    !SOURCE TERMS FOR 2-COMP OF B-FIELD
    deriv1=derivative2D2(B2(1:lx1,1:lx2)*u1(1:lx1,1:lx2),dx2(1:lx2))
!    deriv1=derivative2D2obj(B2(1:lx1,1:lx2)*u1(1:lx1,1:lx2),dx2(1:lx2),x1,x2)
    !$omp parallel do
    do ix2=1,lx2
      B1(1:lx1,ix2)=B1(1:lx1,ix2)+ &
                       dt*deriv1(:,ix2)
    end do
  
  
    !SOURCE TERMS FOR 2-COMP OF B-FIELD
    deriv1=derivative2D1(B1(1:lx1,1:lx2)*u2(1:lx1,1:lx2),dx1(1:lx1))
!    deriv1=derivative2D1obj(B1(1:lx1,1:lx2)*u2(1:lx1,1:lx2),dx1(1:lx1),x1,x2)
    !$omp parallel do
    do ix2=1,lx2
      B2(1:lx1,ix2)=B2(1:lx1,ix2)+ &
                       dt*deriv1(:,ix2)        
    end do  
  
      
    !SOURCE TERMS FOR 3-COMP OF B-FIELD
  !  grad1B1u3=derivative(B1(1:lx1)*u3(1:lx1),dx1(1:lx1))
  !  grad1B1u3(1)=(B1(2)*u3(2)-B1(0)*u3(0))/(dx1(2)+dx1(1))
  !  grad1B1u3(lx1)=(B1(lx1+1)*u3(lx1+1)-B1(lx1-1)*u3(lx1-1))/(dx1(lx1+1)+dx1(lx1))
    deriv1=derivative2D1(B1(1:lx1,1:lx2)*u3(1:lx1,1:lx2),dx1(1:lx1))
    deriv2=derivative2D2(B2(1:lx1,1:lx2)*u3(1:lx1,1:lx2),dx2(1:lx2))
    !$omp parallel do
    do ix2=1,lx2
      B3(1:lx1,ix2)=B3(1:lx1,ix2)+ &
                dt*(deriv1(:,ix2)+ &
                deriv2(:,ix2))
    end do
  
  
    !SOURCE TERMS FOR INTERNAL ENERGY (RK2 STEPPING):  COMPRESSION/EXPANSION
    deriv1=derivative2D1(u1(1:lx1,1:lx2),dx1(1:lx1))
!    deriv1=derivative2D1obj(u1(1:lx1,1:lx2),dx1(1:lx1),x1,x2)
    deriv2=derivative2D2(u2(1:lx1,1:lx2),dx2(1:lx2))
!    deriv2=derivative2D2obj(u2(1:lx1,1:lx2),dx2(1:lx2),x1,x2)
    !$omp parallel do
    do ix2=1,lx2
      p11(:,ix2)=rhoeps(1:lx1,ix2)*(gammap-1d0)+Q11(:,ix2)
      p22(:,ix2)=rhoeps(1:lx1,ix2)*(gammap-1d0)+Q22(:,ix2)
      sourceterm(:,ix2)=-p11(:,ix2)*deriv1(:,ix2)-p22(:,ix2)*deriv2(:,ix2)
      rhoepshalf(:,ix2)=rhoeps(1:lx1,ix2)+dt/2d0*sourceterm(:,ix2)
  
      p11(:,ix2)=rhoepshalf(:,ix2)*(gammap-1d0)+Q11(:,ix2)
      p22(:,ix2)=rhoepshalf(:,ix2)*(gammap-1d0)+Q22(:,ix2)
      sourceterm(:,ix2)=-p11(:,ix2)*deriv1(:,ix2)-p22(:,ix2)*deriv2(:,ix2)
      rhoeps(1:lx1,ix2)=rhoeps(1:lx1,ix2)+dt*sourceterm(:,ix2)
    end do
  

  end subroutine MHD_sources



  subroutine visc_vel_diffs(u1,u2,du1full,du2full)

    !------------------------------------------------------------
    !-------LOAD UP VELOCITY DIFFERENCE ARRAYS FOR ART. VISCOSITY
    !-------INPUT VELOCS. HAVE GHOST CELLS.
    !------------------------------------------------------------

    real(8), dimension(-1:,-1:), intent(in) :: u1,u2
    real(8), dimension(:,:), intent(out) :: du1full,du2full

    integer :: ix1,ix2,lx1,lx2


    !SIZES WITHOUT GHOST CELLS
    lx1=size(u1,1)-4; lx2=size(u1,2)-4;


    !ART. VISCOSITY DIFFS OF U1 AND U2
    !$omp parallel do
    do ix2=1,lx2
      du1full(:,ix2)=u1(2:lx1+1,ix2)-u1(0:lx1-1,ix2)
    end do
    !$omp parallel do
    do ix1=1,lx1
      du2full(ix1,:)=u2(ix1,2:lx2+1)-u2(ix1,0:lx2-1)
    end do

  end subroutine visc_vel_diffs


end module MHD
