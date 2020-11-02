module phys_consts 

implicit none

real(8), parameter :: pi=3.14159d0
real(8), parameter :: amu=1.67d-27
real(8), parameter :: kb=1.38d-23 
real(8), parameter :: gammap=5d0/3d0
real(8), parameter :: mu0=4d0*pi*1d-7
real(8), parameter :: Mearth=5.972d24
real(8), parameter :: Msun=1.988d30
real(8), parameter :: magm_earth=-1d0*7.94d22
real(8), parameter :: Gconst=6.67e-11

real(8), parameter :: Rsun=6.955d8
real(8), parameter :: Rearth=6370d3
real(8), parameter :: Rmoon=1738d3

real(8), parameter :: robj=0d0          !a solid object to include in the flow
real(8), parameter :: x1ctr=0d0         !center of object
real(8), parameter :: x2ctr=0d0         !center of object

real(8), parameter :: B1obj=1d-9        !this block is to set state parameters inside object
real(8), parameter :: B2obj=0d0
real(8), parameter :: B3obj=0d0
real(8), parameter :: nobj=1d6
real(8), parameter :: Tobj=1d4

real(8), dimension(:,:), allocatable :: g1,g2


contains


  subroutine flow_obs(t,x1,x2,rhom,rhov1,rhov2,rhov3,rhoeps,B1,B2,B3,v1,v2,v3)

    !------------------------------------------------------------
    !-------CREATES AN OBSTRUCTION TO FLOW ON THE GRID.  HAVE
    !-------THIS SUBPROGRAM DO NOTHING FOR A `FREE' INTERIOR.
    !------------------------------------------------------------
  
    real(8), intent(in) :: t
    real(8), dimension(-1:), intent(in) :: x1
    real(8), dimension(-1:), intent(in) :: x2
    real(8), dimension(-1:,-1:), intent(inout) :: rhom,rhov1,rhov2,rhov3,rhoeps,B1,B2,B3,v1,v2,v3

    integer :: ix1,ix2,lx1,lx2

    real(8) :: Ti=1d4                             !temperature used to initialize internal energy
    real(8) :: r


    !SIZES WITHOUT GHOST CELLS
    lx1=size(rhom,1)-4; lx2=size(rhom,2)-4;


    !PICK A CIRCULAR REGION IN WHICH TO ZERO OUT VELOC.
    do ix2=1,lx2
      do ix1=1,lx1
        if (inside_obj(x1(ix1),x2(ix2))) then
          rhov1(ix1,ix2)=0d0
          rhov2(ix1,ix2)=0d0
          rhom(ix1,ix2)=nobj*amu
          rhoeps(ix1,ix2)=rhom(ix1,ix2)*kb*Tobj/amu/(gammap-1d0)

          B1(ix1,ix2)=B1obj
          B2(ix1,ix2)=B2obj
          B3(ix1,ix2)=B3obj

!          r=sqrt(x2(ix2)**2+x1(ix1)**2)
!          if (r<Rearth) then
!            r=Rearth
!            B2(ix1,ix2)=0d0
!            B3(ix1,ix2)=0d0
!            B1(ix1,ix2)=mu0*magm_earth/4d0/pi/(r**5)*(3*r**2-r**2)
!!            B1(ix1,ix2)=0d0
!          else
!            B2(ix1,ix2)=mu0*magm_earth/4d0/pi/(r**5)*3*x2(ix2)*x1(ix1)
!            B3(ix1,ix2)=0d0
!            B1(ix1,ix2)=mu0*magm_earth/4d0/pi/(r**5)*(3*x1(ix1)**2-r**2)
!          end if
        end if
      end do
    end do


    !RECOMPUTE VELOCITIES
    v1=rhov1/rhom
    v2=rhov2/rhom
    v3=rhov3/rhom

  end subroutine flow_obs


  function inside_obj(x1pos,x2pos)

    !------------------------------------------------------------
    !-------DETERMINES WHETHER OR NOT A POINT IS INSIDE THE OBJECT
    !------------------------------------------------------------


    real(8), intent(in) :: x1pos,x2pos

    real(8) :: r

    logical :: inside_obj


    r=sqrt((x1pos-x1ctr)**2+(x2pos-x2ctr)**2)
    if (r<robj) then
      inside_obj=.true.
    else
      inside_obj=.false.
    end if

  end function inside_obj


  subroutine gravity_allocstate(x1,x2)

    !------------------------------------------------------------
    !-------ALLOCATE SPACE TO STORE GRAVITATIONAL FIELD
    !------------------------------------------------------------

    real(8), dimension(-1:), intent(in) :: x1
    real(8), dimension(-1:), intent(in) :: x2
    integer :: ix1,ix2,lx1,lx2


    !SIZES WITHOUT GHOST CELLS
    lx1=size(x1,1)-4
    lx2=size(x2,1)-4

    if (allocated(g1)) then
      deallocate(g1,g2)
    else
      allocate(g1(-1:lx1+2,-1:lx2+2),g2(-1:lx1+2,-1:lx2+1))
      g1=0d0
      g2=0d0
    end if

  end subroutine gravity_allocstate


  subroutine gravity(x1,x2)

    !------------------------------------------------------------
    !-------CALCULATE GRAVITATIONAL FIELD
    !------------------------------------------------------------

    real(8), dimension(-1:), intent(in) :: x1
    real(8), dimension(-1:), intent(in) :: x2

    real(8) :: r
    integer :: ix1,ix2,lx1,lx2


    !SIZES WITHOUT GHOST CELLS
    lx1=size(x1,1)-4
    lx2=size(x2,1)-4


    !CONST ACCEL.
    g1=-9.8d0
    g2=0d0

!   !NEWTON'S LAW
!    do ix2=1,lx2
!      do ix1=1,lx1
!        r=x1(ix1)-x1(1)+Rsun
!        g1(ix1,ix2)=-Gconst*Msun/r**2
!        g2(ix1,ix2)=0d0
!      end do
!    end do

  end subroutine gravity


end module phys_consts
