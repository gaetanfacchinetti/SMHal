# SMHal [Simplified Models and (sub)Halos]

This software has been built to evaluate the properties of dark matter (sub)halos in a simplified dark matter model. The basic usage consist in setting the parameters of the Lagrangian in input and evaluating the dark matter abundance today, the minimal mass of halos, the number of subhalos in given specific host objects (taking into account tidal stripping effect) and their intrinsic luminosity. Dynamical effects within halos are evaluated through an update of the analytical model (Stref & Lavalle 2017). The properties of particles can be constrained from the observed abundance of dark matter in the Universe as well as direct detection and indirect detection. 

In order to be time-efficient all computations are made from scratch and do not rely on preexisting libraries. The only external inputs are from DDCalc ().

For more details please contact : gaetan.facchinetti@umontpellier.fr

## Installation

Requirement: In order to compile and run, this code needs the gnu scientific library gsl (https://www.gnu.org/software/gsl) and eign3 library (https://eigen.tuxfamily.org). Please make sure they are installed before going forward.

In order to download the software you can run: `git clone ... ` or use the download button. For the installation you need to open the Makefile in `/bin` and set the location of the gsl and eign3 libraries. Then run `make all` from the `/bin` folder.

The make command will first compile the DDCalc library before compiling the source of the software.

## Usage

