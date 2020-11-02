module numerical 

use phys_consts                     !for differentiating around a object embedded in flow

implicit none

contains



  function derivative2D1(f,dx1)
    real(8), dimension(:,:), intent(in) :: f    !ghost cells not present...
    real(8), dimension(:), intent(in) :: dx1    !presumed backward diffs.

    integer :: lx1,lx2,ix2
    real(8), dimension(1:size(f,1),1:size(f,2)) :: derivative2D1

    lx1=size(f,1)
    lx2=size(f,2)

    !$omp parallel do
    do ix2=1,lx2
      derivative2D1(1,ix2)=(f(2,ix2)-f(1,ix2))/dx1(2)
      derivative2D1(2:lx1-1,ix2)=(f(3:lx1,ix2)-f(1:lx1-2,ix2))/(dx1(3:lx1)+dx1(2:lx1-1))
      derivative2D1(lx1,ix2)=(f(lx1,ix2)-f(lx1-1,ix2))/dx1(lx1)
    end do
  end function derivative2D1



  function derivative2D2(f,dx2)
    real(8), dimension(:,:), intent(in) :: f    !ghost cells not present...
    real(8), dimension(:), intent(in) :: dx2    !presumed backward diffs.

    integer :: lx1,lx2,ix1
    real(8), dimension(1:size(f,1),1:size(f,2)) :: derivative2D2

    lx1=size(f,1)
    lx2=size(f,2)

    !$omp parallel do
    do ix1=1,lx1
      derivative2D2(ix1,1)=(f(ix1,2)-f(ix1,1))/dx2(2)
      derivative2D2(ix1,2:lx2-1)=(f(ix1,3:lx2)-f(ix1,1:lx2-2))/(dx2(3:lx2)+dx2(2:lx2-1))
      derivative2D2(ix1,lx2)=(f(ix1,lx2)-f(ix1,lx2-1))/dx2(lx2)
    end do
  end function derivative2D2


  function derivative2D2p(f,dx2)

    !Derivative on a periodic domain

    real(8), dimension(:,:), intent(in) :: f    !ghost cells not present...
    real(8), dimension(:), intent(in) :: dx2    !presumed backward diffs.

    integer :: lx1,lx2,ix1
    real(8), dimension(1:size(f,1),1:size(f,2)) :: derivative2D2p

    lx1=size(f,1)
    lx2=size(f,2)

    !$omp parallel do
    do ix1=1,lx1
      derivative2D2p(ix1,1)=(f(ix1,2)-f(ix1,lx2))/(dx2(2)+dx2(1))
      derivative2D2p(ix1,2:lx2-1)=(f(ix1,3:lx2)-f(ix1,1:lx2-2))/(dx2(3:lx2)+dx2(2:lx2-1))
      derivative2D2p(ix1,lx2)=(f(ix1,1)-f(ix1,lx2-1))/(dx2(1)+dx2(lx2))
    end do

  end function derivative2D2p


  function derivative2D1obj(f,dx1,x1,x2)

    !DIFFERENTIATE AROUND AN OBJECT,x1

    real(8), dimension(:,:), intent(in) :: f    !ghost cells not present...
    real(8), dimension(:), intent(in) :: dx1    !presumed backward diffs.
    real(8), dimension(-1:), intent(in) :: x1
    real(8), dimension(-1:), intent(in) :: x2

    integer :: lx1,lx2,ix1,ix2
    real(8) :: ffwd,fback,dxval
    real(8), dimension(1:size(f,1),1:size(f,2)) :: derivative2D1obj


    lx1=size(f,1)
    lx2=size(f,2)

    do ix2=1,lx2
      derivative2D1obj(1,ix2)=(f(2,ix2)-f(1,ix2))/dx1(2)
      do ix1=2,lx1-1
        if (inside_obj(x1(ix1),x2(ix2))) then          !location of interest is inside object, do nothing
          derivative2D1obj(ix1,ix2)=0d0
        else
          if (inside_obj(x1(ix1+1),x2(ix2))) then        !point is outside obj, but need to check foward location.  If inside, back diff.
            ffwd=f(ix1,ix2)
!            fback=f(ix1-1,ix2)
            fback=ffwd    !force zero
            dxval=dx1(ix1)
          elseif (inside_obj(x1(ix1-1),x2(ix2))) then    !is the backward point inside object?  Use fwd diff.
            ffwd=f(ix1+1,ix2)
!            fback=f(ix1,ix2)
            fback=ffwd
            dxval=dx1(ix1+1)
          else                                           !neither is inside, so differentiate 2nd order centered
            ffwd=f(ix1+1,ix2)
            fback=f(ix1-1,ix2)
            dxval=dx1(ix1+1)+dx1(ix1)
          end if

          derivative2D1obj(ix1,ix2)=(ffwd-fback)/dxval
        end if
      end do
      derivative2D1obj(lx1,ix2)=(f(lx1,ix2)-f(lx1-1,ix2))/dx1(lx1)
    end do

  end function derivative2D1obj


  function derivative2D2obj(f,dx2,x1,x2)

    !DIFFERENTIATE AROUND AN OBJECT,x2

    real(8), dimension(:,:), intent(in) :: f    !ghost cells not present...
    real(8), dimension(:), intent(in) :: dx2    !presumed backward diffs.
    real(8), dimension(-1:), intent(in) :: x1
    real(8), dimension(-1:), intent(in) :: x2

    integer :: lx1,lx2,ix1,ix2
    real(8) :: ffwd,fback,dxval
    real(8), dimension(1:size(f,1),1:size(f,2)) :: derivative2D2obj


    lx1=size(f,1)
    lx2=size(f,2)

    do ix1=1,lx1
      derivative2D2obj(1,ix2)=(f(ix1,2)-f(ix1,1))/dx2(2)
      do ix2=2,lx2-1
        if (inside_obj(x1(ix1),x2(ix2))) then            !location of interest is inside object, do nothing
          derivative2D2obj(ix1,ix2)=0d0
        else
          if (inside_obj(x1(ix1),x2(ix2+1))) then        !point is outside obj, but need to check foward location.  If inside, back diff.
            ffwd=f(ix1,ix2)
!            fback=f(ix1,ix2-1)
            fback=ffwd
            dxval=dx2(ix2)
          elseif (inside_obj(x1(ix1),x2(ix2-1))) then    !is the backward point inside object?  Use fwd diff.
            ffwd=f(ix1,ix2+1)
!            fback=f(ix1,ix2)
            fback=ffwd
            dxval=dx2(ix2+1)
          else                                           !neither is inside, so differentiate 2nd order centered
            ffwd=f(ix1,ix2+1)
            fback=f(ix1,ix2-1)
            dxval=dx2(ix2+1)+dx2(ix2)
          end if

          derivative2D2obj(ix1,ix2)=(ffwd-fback)/dxval
        end if
      end do
      derivative2D2obj(ix1,lx2)=(f(ix1,lx2)-f(ix1,lx2-1))/dx2(lx2)
    end do

  end function derivative2D2obj


  function integral2D2(f,dx2)

    !------------------------------------------------------------
    !-------COMPUTE AN INTEGRAL ALONG THE 2-DIM.  IT IS EXPECTED THAT 
    !-------GHOST CELLS WILL HAVE BEEN TRIMMED FROM ARRAYS BEFORE
    !-------THEY ARE PASSED INTO THIS ROUTINE
    !------------------------------------------------------------

    real(8), dimension(:,:), intent(in) :: f
    real(8), dimension(1:size(f,2)), intent(in) :: dx2  !presumed backward diffs. (1:lx2)

    integer :: ix2,lx2
    real(8), dimension(1:size(f,1),1:size(f,2)) :: integral2D2

    lx2=size(f,2)

    integral2D2(:,1)=0d0
    do ix2=2,lx2
      integral2D2(:,ix2)=integral2D2(:,ix2-1)+0.5d0*(f(:,ix2)+f(:,ix2-1))*dx2(ix2)
    end do

  end function integral2D2


  function advec2D_MC(f,v1i,v2i,dt,dx1,dx1i,dx2,dx2i)
    real(8), dimension(-1:,-1:), intent(in) :: f    !f includes ghost cells
    real(8), dimension(:,:), intent(in) :: v1i
    real(8), dimension(0:), intent(in) :: dx1       !ith backward difference
    real(8), dimension(:), intent(in) :: dx1i
    real(8), dimension(:,:), intent(in) :: v2i
    real(8), dimension(0:), intent(in) :: dx2
    real(8), dimension(:), intent(in) :: dx2i
    real(8), intent(in) :: dt

    integer :: ix1,ix2,lx1,lx2
    real(8), dimension(-1:size(f,1)-2) :: fx1slice
    real(8), dimension(1:size(f,1)-3) :: v1slice    
    real(8), dimension(-1:size(f,2)-2) :: fx2slice
    real(8), dimension(1:size(f,2)-3) :: v2slice
    real(8), dimension(-1:size(f,1)-2,-1:size(f,2)-2) :: advec2D_MC

    lx1=size(f,1)-4
    lx2=size(f,2)-4

    !x2-sweep (it appears that fortran doesn't differentiate b/t row and column arrays, so no reshapes...)
    !omp parallel do private(fx2slice,v2slice)
    do ix1=1,lx1
      fx2slice=f(ix1,:)
      v2slice=v2i(ix1,:)
      fx2slice=advec1D_MC(fx2slice,v2slice,dt,dx2,dx2i)
      advec2D_MC(ix1,:)=fx2slice
    end do

    !copy x2 boundary conditions to partially updated variable for next sweep
    !advec2D_MC(:,-1:0)=f(:,-1:0);    
    !advec2D_MC(:,lx2+1:lx2+2)=f(:,lx2+1:lx2+2);  
    advec2D_MC(-1:0,:)=f(-1:0,:);    
    advec2D_MC(lx1+1:lx1+2,:)=f(lx1+1:lx1+2,:);  

    !x1-sweep
    !$omp parallel do private(fx1slice,v1slice)
    do ix2=1,lx2
      fx1slice=advec2D_MC(:,ix2); 
      v1slice=v1i(:,ix2);
      fx1slice=advec1D_MC(fx1slice,v1slice,dt,dx1,dx1i)
      advec2D_MC(:,ix2)=fx1slice;
    end do  

  end function advec2D_MC



  function advec1D_MC(f,v1i,dt,dx1,dx1i)
    real(8), dimension(-1:), intent(in) :: f    !f includes ghost cells
    real(8), dimension(:), intent(in) :: v1i
    real(8), dimension(0:), intent(in) :: dx1   !ith backward difference
    real(8), dimension(:), intent(in) :: dx1i
    real(8), intent(in) :: dt
    
    integer :: ix1,lx1
    real(8), dimension(size(v1i)) :: phi
    real(8), dimension(0:size(v1i)) :: slope         !slopes only need through first layer of ghosts
    real(8) :: lslope,rslope,cslope
    real(8), dimension(-1:size(f)-2) :: advec1D_MC


    lx1=size(f)-4

    !Slopes
    lslope=(f(0)-f(-1))/dx1(0)
    do ix1=0,lx1+1    
      rslope=(f(ix1+1)-f(ix1))/dx1(ix1+1)
      cslope=(f(ix1+1)-f(ix1-1))/(dx1(ix1)+dx1(ix1+1))
      slope(ix1)=minmod(cslope,minmod(2*lslope,2*rslope))

      lslope=rslope
    end do


    !Slope-limited flux at ***left*** wall of cell ix1. 
    !The treatment of slope here (ie the dx1(ix1)) assumes that the grid point isp centered within the cell
    do ix1=1,lx1+1
      if (v1i(ix1) < 0) then 
        phi(ix1)=f(ix1)*v1i(ix1) - 0.5*v1i(ix1)*(dx1(ix1)+v1i(ix1)*dt)*slope(ix1)
      else
        phi(ix1)=f(ix1-1)*v1i(ix1) + 0.5*v1i(ix1)*(dx1(ix1)-v1i(ix1)*dt)*slope(ix1-1)
      end if
    end do


    !flux differencing form
    advec1D_MC(1:lx1)=f(1:lx1)-dt*(phi(2:lx1+1)-phi(1:lx1))/dx1i


    !default to copying ghost cells
    advec1D_MC(-1:0)=f(-1:0)
    advec1D_MC(lx1+1:lx1+2)=f(lx1+1:lx1+2)
  end function advec1D_MC



  function minmod(a,b)
    real(8), intent(in) :: a,b
    real(8) :: minmod

    if (a*b <= 0.) then
      minmod=0.
    else if (abs(a) < abs(b)) then
      minmod=a
    else
      minmod=b
    end if
  end function minmod



end module numerical
