program MHD2D

!A 'CLASSICAL HYDRODYNAMIC' SOLVER (WAVES PROPAGATED THRU SOURCE TERMS) 
!FOR THE IDEAL MHD EQUATIONS IN 2 DIMENSIONS.  M. ZETTERGREN, T. SYMONS, K. CHAN
!
! LOOP LEVEL PARALLELIZATION HAS BEEN IMPLEMENTED USING OPENMP

use phys_consts
use numerical
use MHD

implicit none

integer, parameter :: outunit=42,method=1
!integer, parameter :: npts1=1001, npts2=1001
!integer, parameter :: npts1=501, npts2=501
integer, parameter :: npts1=251, npts2=251
real(8), parameter :: tcfl=0.75d0

!real(8), parameter :: stride1=2d3,stride2=2d3
!real(8), parameter :: stride1=4d3,stride2=4d3
real(8), parameter :: stride1=8d3,stride2=8d3

real(8), dimension(npts1) :: dx1i
real(8), dimension(-1:npts1+2) :: x1
real(8), dimension(npts1+1) :: x1i
real(8), dimension(0:npts1+2) :: dx1

real(8), dimension(npts2) :: dx2i
real(8), dimension(-1:npts2+2) :: x2
real(8), dimension(npts2+1) :: x2i
real(8), dimension(0:npts2+2) :: dx2

integer :: it,lx1,lx2,ix1,ix2

real(8), dimension(-1:npts1+2,-1:npts2+2) :: rhom,rhou1,rhou2,rhou3,B1,B2,B3,rhoeps,u1,u2,u3

real(8), dimension(npts1,npts2) :: du1full,du2full,vcon1,vcon2
real(8), dimension(npts1,npts2) :: p,vA,vsnd,vmatter

real(8) :: t=0d0,dt=1d-6,dtout=5d0,tdur=120d0
real(8) :: dt1,dt2
real(8) :: toutnext

real(8) :: tstart,tfin                           !for assessing performance of code blocks

real(8) :: Ti=8000d0                             !temperature used to initialize internal energy
real(8) :: vsw=15d3*5d0                          !initial condition solar wind speed



!GRID CONSTRUCTION
write(*,*) 'Creating grid of size:  ',npts1,npts2
x1=[ (ix1*stride1, ix1=-2,npts1+1) ]
lx1=size(x1)-4                             !exclude ghost cells in count
x1=x1-0.5d0*(x1(1)+x1(lx1))    !center the grid
dx1=x1(0:lx1+2)-x1(-1:lx1+1)                !backward diffs
x1i(1:lx1+1)=0.5d0*(x1(0:lx1)+x1(1:lx1+1))    !interface data
dx1i=x1i(2:lx1+1)-x1i(1:lx1)

x2=[ (ix2*stride2, ix2=-2,npts2+1) ]
lx2=size(x2)-4                             !exclude ghost cells in count
x2=x2-0.5d0*(x2(1)+x2(lx2))    !center
dx2=x2(0:lx2+2)-x2(-1:lx2+1)                !backward diffs
x2i(1:lx2+1)=0.5d0*(x2(0:lx2)+x2(1:lx2+1))    !interface data
dx2i=x2i(2:lx1+1)-x2i(1:lx2)


!PREP OUTPUT FILE (BINARY STREAM), WRITE GRID SIZES AND GRID
write(*,*) 'Opening output file...'
open(outunit,file='MHD2D.dat',status='replace',access='stream')
write(outunit) lx1
write(outunit) x1
write(outunit) lx2
write(outunit) x2


!INITIAL CONDITIONS
write(*,*) 'Setting initial conditions...'
rhom=2d9*amu
rhou1=0d0
rhou2=2d9*amu*vsw
!rhou2=0d0
rhou3=0d0
B1=0d0
B2=0d0
B3=0d0
rhoeps=rhom*kb*Ti/amu/(gammap-1)


!SET VELOCITIES
u1=rhou1/rhom
u2=rhou2/rhom
u3=rhou3/rhom


!MAIN LOOP
toutnext=0d0
do while (t<tdur)
  call cpu_time(tstart)


  !VELOCITY DIFFERENCES FOR CALCULATING ARTIFICAL VISCOSITY
  call visc_vel_diffs(u1,u2,du1full,du2full)


  !TIME STEP DETERMINATION (A BIT UGLY...)
  !$omp parallel do
  do ix2=1,lx2
    vA(:,ix2)=sqrt((B1(1:lx1,ix2)**2+B2(1:lx1,ix2)**2+B3(1:lx1,ix2)**2)/mu0/rhom(1:lx1,ix2))
    p(:,ix2)=rhoeps(1:lx1,ix2)*(gammap-1d0)
    vsnd(:,ix2)=sqrt(gammap*p(:,ix2)/rhom(1:lx1,ix2))
    vmatter(:,ix2)=sqrt(u1(1:lx1,ix2)**2+u2(1:lx1,ix2)**2)
    vcon1(:,ix2)=abs(0.5d0*xi**2*min(du1full(:,ix2),0d0))
    vcon2(:,ix2)=abs(0.5d0*xi**2*min(du2full(:,ix2),0d0))
  end do

  dt1=tcfl*minval(dx1i/(maxval(vA,2)+ & 
                      maxval(vsnd,2)+ &
                      maxval(vcon1,2)+ &
                      maxval(abs(u1(1:lx1,1:lx2)),2)))
  dt2=tcfl*minval(dx2i/(maxval(vA,1)+ &     
                      maxval(vsnd,1)+ &
                      maxval(vcon2,1)+ &
                      maxval(abs(u2(1:lx1,1:lx2)),1)))
  dt=min(dt1,dt2)
  t=t+dt                                        !time after this step
  write(*,*) 'Time:  ',t,'  Time step:  ',dt
  write(*,*) '  Alfven speed:  ',maxval(pack(vA,.true.))/1d3 
  write(*,*) '  Sound speed:  ',maxval(pack(vsnd,.true.))/1d3
  write(*,*) '  Matter speed:  ',maxval(pack(vmatter,.true.))/1d3
  write(*,*) '    viscous-1 speed:  ',maxval(pack(vcon1,.true.))/1d3, &
             ' viscous-2 speed:  ',maxval(pack(vcon2,.true.))/1d3


  !FILL GHOST CELLS
  call boundaries(t,x1,x2,rhom,rhou1,rhou2,rhou3,rhoeps,B1,B2,B3)


  !PICK A REGION IN WHICH TO INTRODUCE AN OBJECT INTO THE FLOW
  call flow_obs(t,x1,x2,rhom,rhou1,rhou2,rhou3,rhoeps,B1,B2,B3)


  !ADVECTION SUBSTEP
  call MHD_advection(dt,dx1,dx1i,dx2,dx2i,rhom,rhou1,rhou2,rhou3,rhoeps,B1,B2,B3)



  !SOURCES SUBSTEP
  call MHD_sources(dt,dx1,dx2,du1full,du2full,rhom,rhou1,rhou2,rhou3,rhoeps, &
                   B1,B2,B3,u1,u2,u3)



  !OUTPUT (RAW BINARY STREAM)
  if (t>toutnext) then
    write(*,*) 'Initiating output for time step:  ',t
    write(outunit) t
    write(outunit) rhom/amu
    write(outunit) u1
    write(outunit) u2
    write(outunit) u3
    write(outunit) B1
    write(outunit) B2
    write(outunit) B3
    write(outunit) rhoeps*(gammap-1)
    toutnext=toutnext+dtout
  end if


  call cpu_time(tfin)
  write(*,*) 'Time step completed in cpu time of:  ',tfin-tstart
end do

close(outunit)



contains

  subroutine boundaries(t,x1,x2,rhom,rhou1,rhou2,rhou3,rhoeps,B1,B2,B3)

    !------------------------------------------------------------
    !-------FILL THE GHOST CELLS TO ENFORCE BOUNDARY CONDITIONS
    !------------------------------------------------------------
  
    real(8), intent(in) :: t
    real(8), dimension(-1:), intent(in) :: x1
    real(8), dimension(-1:), intent(in) :: x2
    real(8), dimension(-1:,-1:), intent(inout) :: rhom,rhou1,rhou2,rhou3,rhoeps,B1,B2,B3  

    real(8) :: vsw


    !SIZES WITHOUT GHOST CELLS
    lx1=size(rhom,1)-4; lx2=size(rhom,2)-4;


    !SET THE SOLAR WIND SPEED AT BOUNDARY (COULD BE TIME-DEPENDENT)
    vsw=15d3*5d0


    !BOUNDARY CONDITIONS (THESE ARE FREE-FLOW)
    rhom(-1:0,1:lx2)=spread(rhom(1,1:lx2),1,2)              !bottom
    rhom(lx1+1:lx1+2,1:lx2)=spread(rhom(1,1:lx2),1,2)       !top
    rhom(1:lx1,-1:0)=2d9*amu                                !left
    rhom(1:lx1,lx2+1:lx2+2)=spread(rhom(1:lx1,lx2),2,2)     !right
  
    rhou1(-1:0,1:lx2)=spread(rhou1(1,1:lx2),1,2)
    rhou1(lx1+1:lx1+2,1:lx2)=spread(rhou1(lx1,1:lx2),1,2)
    rhou1(1:lx1,-1:0)=spread(rhou1(1:lx1,1),2,2)
    rhou1(1:lx1,lx2+1:lx2+2)=spread(rhou1(1:lx1,lx2),2,2)
  
    rhou2(-1:0,1:lx2)=spread(rhou2(1,1:lx2),1,2)
    rhou2(lx1+1:lx1+2,1:lx2)=spread(rhou2(1,1:lx2),1,2)
    rhou2(1:lx1,-1:0)=2d9*amu*vsw
    rhou2(1:lx1,lx2+1:lx2+2)=spread(rhou2(1:lx1,lx2),2,2)
  
    rhou3(-1:0,1:lx2)=spread(rhou3(1,1:lx2),1,2)
    rhou3(lx1+1:lx1+2,1:lx2)=spread(rhou3(lx1,1:lx2),1,2)
    rhou3(1:lx1,-1:0)=spread(rhou3(1:lx1,1),2,2)
    rhou3(1:lx1,lx2+1:lx2+2)=spread(rhou3(1:lx1,lx2),2,2)
  
    B1(-1:0,1:lx2)=spread(B1(1,1:lx2),1,2)
    B1(lx1+1:lx1+2,1:lx2)=spread(B1(lx1,1:lx2),1,2)
    B1(1:lx1,-1:0)=spread(B1(1:lx1,1),2,2)
!    B1(1:lx1,-1:0)=1d-9
    B1(1:lx1,lx2+1:lx2+2)=spread(B1(1:lx1,lx2),2,2)
  
    B2(-1:0,1:lx2)=spread(B2(1,1:lx2),1,2)
    B2(lx1+1:lx1+2,1:lx2)=spread(B2(lx1,1:lx2),1,2)    
    B2(1:lx1,-1:0)=spread(B2(1:lx1,1),2,2)
    B2(1:lx1,lx2+1:lx2+2)=spread(B2(1:lx1,lx2),2,2)
  
    B3(-1:0,1:lx2)=spread(B3(1,1:lx2),1,2)
    B3(lx1+1:lx1+2,1:lx2)=spread(B3(lx1,1:lx2),1,2)   
!    B3(1:lx1,-1:0)=spread(B3(1:lx1,1),2,2)
    B3(1:lx1,-1:0)=1d-9    !works well
    B3(1:lx1,lx2+1:lx2+2)=spread(B3(1:lx1,lx2),2,2)
  
    rhoeps(-1:0,1:lx2)=spread(rhoeps(1,1:lx2),1,2)
    rhoeps(lx1+1:lx1+2,1:lx2)=spread(rhoeps(lx1,1:lx2),1,2)
    rhoeps(1:lx1,-1:0)=spread(rhoeps(1:lx1,1),2,2)
    rhoeps(1:lx1,lx2+1:lx2+2)=spread(rhoeps(1:lx1,lx2),2,2)

  end subroutine boundaries


  subroutine flow_obs(t,x1,x2,rhom,rhou1,rhou2,rhou3,rhoeps,B1,B2,B3)

    !------------------------------------------------------------
    !-------CREATES AN OBSTRUCTION TO FLOW ON THE GRID.  HAVE
    !-------THIS SUBPROGRAM DO NOTHING FOR A `FREE' INTERIOR.
    !------------------------------------------------------------
  
    real(8), intent(in) :: t
    real(8), dimension(-1:), intent(in) :: x1
    real(8), dimension(-1:), intent(in) :: x2
    real(8), dimension(-1:,-1:), intent(inout) :: rhom,rhou1,rhou2,rhou3,rhoeps,B1,B2,B3 

    real(8) :: x1ctr,x2ctr,dist,objrad               !for solid object embedded in flow


    !SIZES WITHOUT GHOST CELLS
    lx1=size(rhom,1)-4; lx2=size(rhom,2)-4;


    !PICK A CIRCULAR REGION IN WHICH TO ZERO OUT VELOC.
    x1ctr=0d0
    x2ctr=0d0
    objrad=25d0*stride1/2d0
    do ix2=1,lx2
      do ix1=1,lx1
        dist=sqrt(x1(ix1)**2+x2(ix2)**2)
        if (dist<=objrad) then
          rhou1(ix1,ix2)=0d0
          rhou2(ix1,ix2)=0d0
          rhom(ix1,ix2)=2d9*amu
          rhoeps(ix1,ix2)=rhom(ix1,ix2)*kb*Ti/amu/(gammap-1)
          B1(ix1,ix2)=0d0
          B2(ix1,ix2)=0d0
          B3(ix1,ix2)=0d0
        end if
      end do
    end do

  end subroutine flow_obs


end program MHD2D
