FC=gfortran
FL=gfortran
#FLAGS_OMP=-fopenmp -fmax-stack-var-size=1024 -fbounds-check
FLAGS_OMP=-fopenmp -fmax-stack-var-size=1024
FLAGS_OPTIM=-O3

all: MHD2D

MHD2D: MHD2D.o MHD.o numerical.o phys_consts.o
	$(FL) $(FLAGS_OPTIM) $(FLAGS_OMP) -o MHD2D phys_consts.o numerical.o MHD.o MHD2D.o

MHD2D.o: MHD2D.f90 MHD.o numerical.o phys_consts.o
	$(FC) $(FLAGS_OPTIM) $(FLAGS_OMP) -c MHD2D.f90

MHD.o: MHD.f90 numerical.o phys_consts.o
	$(FC) $(FLAGS_OPTIM) $(FLAGS_OMP) -c MHD.f90

numerical.o: numerical.f90 phys_consts.o
	$(FC) $(FLAGS_OPTIM) $(FLAGS_OMP) -c numerical.f90

phys_consts.o: phys_consts.f90
	$(FC) $(FLAGS_OPTIM) $(FLAGS_OMP) -c phys_consts.f90

clean:
	rm *.mod
	rm *.o
