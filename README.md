# MHD2D

A two-dimensional numerical solution to the ideal MHD equations.

## Installation

This directory includes source code and makefile for several numerical MHD demos including a blast wave, an external shock interacting with a blast wave, an example of magnetic tension, and a magnetosonic wave example.  Once built these examples may be run by invoking:

./MHD2D_blast
./MHD2D_shockblast
./MHD2D_tension
./MHD2D_waves 

respectively.  Once the run is done the output data are stored in the file MHD2D.dat.  A matlab script for reading in this output file is included as MHD2D.m.  When run from matlab, this will produce a sequence of plots that can be advanced in time by hitting any key.  To my knowledge, this plotting program also works with GNU/Octave 3.8.



The examples referenced above may be built, once pre-requisites are satisfied by:

make -f MHD2D.make

The compiler will complain about the manual stack size changes (necessary to avoid running out of memory and seg. faulting under Mac OS), but will produce correct output.  The result will be a multithreaded executable resulting from loop-level parallelization of the program.  It is approximately 30-40% faster than a single thread compile.  To turn off multithreading nullify the FLAG_OMP variable in the makefile.



For building on ubuntu-style linux, gfortran and make are required.  These may be installed by:

sudo apt-get install make
sudo apt-get install gfortran



For building under Mac OS, install Xcode, then Macports.  Then install the GNU compiler collection port (presently the latest version is gcc49):

sudo port install gcc49

Alternate versions of gcc could be used, as well.  You will need to link the shortcut gfortran (used by the makefile) to the actual executable.  E.g. for gcc49:

cd /usr/local/bin
sudo ln -s  /opt/local/bin/gfortran-mp-4.9 gfortran



These codes have never been installed to tested under a windows environment.
