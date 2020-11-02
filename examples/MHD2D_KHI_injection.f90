program MHD2D

!-------------------------------------------------------------------------------
!-------THE EP 410 MHD MODEL:
!-------A 'CLASSICAL HYDRODYNAMIC' SOLVER (WAVES PROPAGATED THRU SOURCE TERMS) 
!-------FOR THE IDEAL MHD EQUATIONS IN 2 DIMENSIONS.
!-------
!-------AUTHORS:  M. ZETTERGREN, T. SYMONS, K. CHAN
!-------
!-------LOOP LEVEL PARALLELIZATION HAS BEEN IMPLEMENTED USING OPENMP
!-------------------------------------------------------------------------------


use phys_consts
use numerical
use MHD

implicit none

integer, parameter :: outunit=42,method=1

integer, parameter :: npts1=501, npts2=501
real(8), parameter :: stride1=10d3,stride2=10d3

real(8), parameter :: tcfl=0.75d0

real(8), dimension(npts1) :: dx1i
real(8), dimension(-1:npts1+2) :: x1
real(8), dimension(npts1+1) :: x1i
real(8), dimension(0:npts1+2) :: dx1

real(8), dimension(npts2) :: dx2i
real(8), dimension(-1:npts2+2) :: x2
real(8), dimension(npts2+1) :: x2i
real(8), dimension(0:npts2+2) :: dx2

integer :: it,lx1,lx2,ix1,ix2

real(8), dimension(-1:npts1+2,-1:npts2+2) :: rhom,rhov1,rhov2,rhov3,B1,B2,B3,rhoeps,v1,v2,v3

real(8), dimension(npts1,npts2) :: dv1full,dv2full,vcon1,vcon2
real(8), dimension(npts1,npts2) :: p,vA,vsnd,vmatter

real(8) :: t=0d0,dt=1d-6
real(8) :: dtout=5d0,tdur=400d0                  !output rate and total runtime
real(8) :: dt1,dt2
real(8) :: toutnext

real(8) :: tstart,tfin                           !for assessing performance of code blocks

real(8) :: objrad=Rearth
real(8) :: r



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


!CALCULATE GRAVITY
call gravity_allocstate(x1,x2)
!call gravity(x1,x2)                    !comment out to run without gravity (probably most situations)


!PREP OUTPUT FILE (BINARY STREAM), WRITE GRID SIZES AND GRID
write(*,*) 'Opening output file...'
open(outunit,file='MHD2D.dat',status='replace',access='stream')
write(outunit) lx1
write(outunit) x1
write(outunit) lx2
write(outunit) x2


!INITIAL CONDITIONS
call initialization(t,x1,x2,rhom,rhov1,rhov2,rhov3,rhoeps,B1,B2,B3,v1,v2,v3)


!MAIN LOOP
it=1
toutnext=0d0
do while (t<tdur)
  call cpu_time(tstart)


  !FILL GHOST CELLS
  call boundaries(t,x1,x2,rhom,rhov1,rhov2,rhov3,rhoeps,B1,B2,B3)


  !PICK A REGION IN WHICH TO INTRODUCE AN OBJECT INTO THE FLOW
  call flow_obs(t,x1,x2,rhom,rhov1,rhov2,rhov3,rhoeps,B1,B2,B3,v1,v2,v3)


  !VELOCITY DIFFERENCES FOR CALCULATING ARTIFICAL VISCOSITY
  call visc_vel_diffs(v1,v2,dv1full,dv2full)


  !TIME STEP DETERMINATION (A BIT UGLY...)
  !$omp parallel do
  do ix2=1,lx2
    vA(:,ix2)=sqrt((B1(1:lx1,ix2)**2+B2(1:lx1,ix2)**2+B3(1:lx1,ix2)**2)/mu0/rhom(1:lx1,ix2))
    p(:,ix2)=rhoeps(1:lx1,ix2)*(gammap-1d0)
    vsnd(:,ix2)=sqrt(gammap*p(:,ix2)/rhom(1:lx1,ix2))
    vmatter(:,ix2)=sqrt(v1(1:lx1,ix2)**2+v2(1:lx1,ix2)**2)
    vcon1(:,ix2)=abs(0.5d0*xi**2*min(dv1full(:,ix2),0d0))
    vcon2(:,ix2)=abs(0.5d0*xi**2*min(dv2full(:,ix2),0d0))
  end do


!  !DON'T CALCULATE THE CHARACTERISTIC SPEEDS IN THE REGION
!  do ix2=1,lx2
!    do ix1=1,lx1
!      r=sqrt(x2(ix2)**2+x1(ix1)**2)
!      if (r<objrad) then
!        vA(ix1,ix2)=0d0
!      end if
!    end do
!  end do
!

  dt1=tcfl*minval(dx1i/(maxval(vA,2)+ & 
                      maxval(vsnd,2)+ &
                      maxval(vcon1,2)+ &
                      maxval(abs(v1(1:lx1,1:lx2)),2)))
  dt2=tcfl*minval(dx2i/(maxval(vA,1)+ &     
                      maxval(vsnd,1)+ &
                      maxval(vcon2,1)+ &
                      maxval(abs(v2(1:lx1,1:lx2)),1)))

  if (it==1) then                                                               !begin with a short time step
    dt=1d-6
  else
    dt=min(dt1,dt2)
  end if
  it=it+1

  t=t+dt                                                                        !time after this step
  write(*,*) 'Time:  ',t,'  Time step:  ',dt
  write(*,*) '  Alfven speed:  ',maxval(pack(vA,.true.))/1d3 
  write(*,*) '  Sound speed:  ',maxval(pack(vsnd,.true.))/1d3
  write(*,*) '  Matter speed:  ',maxval(pack(vmatter,.true.))/1d3
  write(*,*) '    viscous-1 speed:  ',maxval(pack(vcon1,.true.))/1d3, &
             '    viscous-2 speed:  ',maxval(pack(vcon2,.true.))/1d3


  !ADVECTION SUBSTEP
  call MHD_advection(dt,dx1,dx1i,dx2,dx2i,rhom,rhov1,rhov2,rhov3,rhoeps,B1,B2,B3)
  rhoeps=max(rhoeps,0d0)                                                       !for hypersonic flows the internal energy can be quite small...


  !SOURCES SUBSTEP
  call MHD_sources(dt,dx1,dx2,x1,x2,dv1full,dv2full,rhom,rhov1,rhov2,rhov3,rhoeps, &
                   B1,B2,B3,v1,v2,v3)
  rhoeps=max(rhoeps,0d0)


  !PICK A REGION IN WHICH TO INTRODUCE AN OBJECT INTO THE FLOW
  call flow_obs(t,x1,x2,rhom,rhov1,rhov2,rhov3,rhoeps,B1,B2,B3,v1,v2,v3)


  !OUTPUT (RAW BINARY STREAM)
  if (t>toutnext) then
    write(*,*) 'Initiating output for time step:  ',t
    write(outunit) t
    write(outunit) rhom/amu
    write(outunit) v1
    write(outunit) v2
    write(outunit) v3
    write(outunit) B1
    write(outunit) B2
    write(outunit) B3
    write(outunit) rhoeps*(gammap-1)
!    write(outunit) g1    

    toutnext=toutnext+dtout
  end if


  call cpu_time(tfin)
  write(*,*) 'Time step completed in cpu time of:  ',tfin-tstart
end do

close(outunit)
call gravity_allocstate(x1,x2)



contains


  subroutine initialization(t,x1,x2,rhom,rhov1,rhov2,rhov3,rhoeps,B1,B2,B3,v1,v2,v3)

    !------------------------------------------------------------
    !-------SET INITIAL CONDITIONS FOR THE SIMULATION
    !------------------------------------------------------------

    real(8), intent(in) :: t
    real(8), dimension(-1:), intent(in) :: x1
    real(8), dimension(-1:), intent(in) :: x2
    real(8), dimension(-1:,-1:), intent(inout) :: rhom,rhov1,rhov2,rhov3,rhoeps,B1,B2,B3,v1,v2,v3  

    integer :: ix1,ix2,lx1,lx2

    real(8) :: vsw                                !initial solar wind (background plasma) speed
    real(8) :: Ti=1d5                             !reference temperature used to initialize internal energy
    real(8) :: n0=1d9                             !reference number density value
    real(8) :: B0=0d0                             !reference value for magnetic field
    real(8) :: r,H,sigmax1,sigmax2


    !SIZES WITHOUT GHOST CELLS
    lx1=size(rhom,1)-4
    lx2=size(rhom,2)-4


    !COMPUTE THE INITIAL SOLAR WIND SPEED
    vsw=0d0


    !NOW STATE VARIABLE INITIAL CONDITIONS ARE COMPUTED
    rhom=n0*amu
!    sigmax1=250d3
!    sigmax2=250d3
!    do ix2=1,lx2
!      do ix1=1,lx1
!        rhom(ix1,ix2)=rhom(ix1,ix2)-0.5d0*n0*amu*exp(-(x1(ix1))**2/2d0/sigmax1**2)*exp(-(x2(ix2))**2/2d0/sigmax2**2)
!      end do
!    end do

    rhov1=0d0
    rhov2=rhom*vsw
    rhov3=0d0
    B1=B0
    B2=0d0
    B3=0d0
    rhoeps=rhom*kb*Ti/amu/(gammap-1d0)

  end subroutine initialization


  subroutine boundaries(t,x1,x2,rhom,rhov1,rhov2,rhov3,rhoeps,B1,B2,B3)

    !------------------------------------------------------------
    !-------FILL THE GHOST CELLS TO ENFORCE BOUNDARY CONDITIONS
    !------------------------------------------------------------
  
    real(8), intent(in) :: t
    real(8), dimension(-1:), intent(in) :: x1
    real(8), dimension(-1:), intent(in) :: x2
    real(8), dimension(-1:,-1:), intent(inout) :: rhom,rhov1,rhov2,rhov3,rhoeps,B1,B2,B3  

    integer :: ix1,ix2,lx1,lx2

    real(8) :: vsw                                !background solar wind speed
    real(8) :: Ti=1d5                             !temperature of plasma injected from the boundary
    real(8) :: n0=1d9                            !number density reference value
    real(8) :: B0=0d0                            !magnetic field reference value
    real(8) :: sigmax1,sigmax2


    !SIZES WITHOUT GHOST CELLS
    lx1=size(rhom,1)-4
    lx2=size(rhom,2)-4


    !SET THE BACKGROUND SPEED AT BOUNDARY (COULD BE TIME-DEPENDENT)
    vsw=0.5d0*sqrt(gammap*kb*Ti/amu)


    !BOUNDARY CONDITIONS FOR EACH OF THE MHD STATE VARIABLES (THIS EXAMPLE IS FREE-FLOW)
    rhom(-1:0,1:lx2)=spread(rhom(1,1:lx2),1,2)              !bottom
    rhom(lx1+1:lx1+2,1:lx2)=spread(rhom(lx1,1:lx2),1,2)     !top
    rhom(1:lx1,-1:0)=spread(rhom(1:lx1,1),2,2)              !left
    rhom(1:lx1,lx2+1:lx2+2)=spread(rhom(1:lx1,lx2),2,2)     !right
  
    rhov1(-1:0,1:lx2)=spread(rhov1(1,1:lx2),1,2)
    rhov1(lx1+1:lx1+2,1:lx2)=spread(rhov1(lx1,1:lx2),1,2)
!    sigmax2=200d3    !doesnt work for some reason
!    do ix2=1,lx2
!      rhov1(lx1+1:lx1+2,ix2)=-1d0*n0*amu*vsw*exp(-x2(ix2)**8/2d0/sigmax2**8)
!    end do
    rhov1(1:lx1,-1:0)=spread(rhov1(1:lx1,1),2,2)
    rhov1(1:lx1,lx2+1:lx2+2)=spread(rhov1(1:lx1,lx2),2,2)
  
    rhov2(-1:0,1:lx2)=spread(rhov2(1,1:lx2),1,2)
    rhov2(lx1+1:lx1+2,1:lx2)=spread(rhov2(1,1:lx2),1,2)
!    rhov2(1:lx1,-1:0)=spread(rhov2(1:lx1,1),2,2)
    sigmax1=200d3
    do ix1=1,lx1
      rhov2(ix1,-1:0)=1d0*n0*amu*vsw*exp(-x1(ix1)**8/2d0/sigmax1**8)
    end do
    rhov2(1:lx1,lx2+1:lx2+2)=spread(rhov2(1:lx1,lx2),2,2)
  
    rhov3(-1:0,1:lx2)=spread(rhov3(1,1:lx2),1,2)
    rhov3(lx1+1:lx1+2,1:lx2)=spread(rhov3(lx1,1:lx2),1,2)
    rhov3(1:lx1,-1:0)=spread(rhov3(1:lx1,1),2,2)
    rhov3(1:lx1,lx2+1:lx2+2)=spread(rhov3(1:lx1,lx2),2,2)
  
    B1(-1:0,1:lx2)=spread(B1(1,1:lx2),1,2)
    B1(lx1+1:lx1+2,1:lx2)=spread(B1(lx1,1:lx2),1,2)
    B1(1:lx1,-1:0)=spread(B1(1:lx1,1),2,2)
    B1(1:lx1,lx2+1:lx2+2)=spread(B1(1:lx1,lx2),2,2)
  
    B2(-1:0,1:lx2)=spread(B2(1,1:lx2),1,2)
    B2(lx1+1:lx1+2,1:lx2)=spread(B2(lx1,1:lx2),1,2)    
    B2(1:lx1,-1:0)=spread(B2(1:lx1,1),2,2)
    B2(1:lx1,lx2+1:lx2+2)=spread(B2(1:lx1,lx2),2,2)
  
    B3(-1:0,1:lx2)=spread(B3(1,1:lx2),1,2)
    B3(lx1+1:lx1+2,1:lx2)=spread(B3(lx1,1:lx2),1,2)   
    B3(1:lx1,-1:0)=spread(B3(1:lx1,1),2,2)
    B3(1:lx1,lx2+1:lx2+2)=spread(B3(1:lx1,lx2),2,2)
  
    rhoeps(-1:0,1:lx2)=spread(rhoeps(1,1:lx2),1,2)
    rhoeps(lx1+1:lx1+2,1:lx2)=spread(rhoeps(lx1,1:lx2),1,2)
    rhoeps(1:lx1,-1:0)=spread(rhoeps(1:lx1,1),2,2)
    rhoeps(1:lx1,lx2+1:lx2+2)=spread(rhoeps(1:lx1,lx2),2,2)

  end subroutine boundaries


end program MHD2D
