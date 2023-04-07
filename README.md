# SMHal [Simplified Models and (sub)Halos]

This software has been built (mainly during my PhD thesis) to evaluate the properties of dark matter (sub)halos in a simplified dark matter model. The basic usage consist in setting the parameters of the Lagrangian in input and evaluating the dark matter abundance today, the minimal mass of halos, the number of subhalos in given specific host objects (taking into account tidal stripping effect) and their intrinsic luminosity. Dynamical effects within halos are evaluated through an update of the analytical model (Stref & Lavalle 2017). The properties of particles can be constrained from the observed abundance of dark matter in the Universe as well as direct detection and indirect detection. 

In order to be time-efficient all computations are made from scratch and do not rely on preexisting libraries. The only external inputs are from the [DDCalc](https://github.com/GambitBSM/DDCalc) package.

gaetan.facchinetti@ulb.be


## Installation

_**WARNING**: the code is not user friendly in the current version and may not be before several updates. You can however use parts of the code provided that you refer to this repository and my work. In particular, please cite my PhD thesis:_
- Gaétan Facchinetti, Analytical study of particle dark matter structuring on small scales and implications for dark matter searches, Université de Montpellier, 2021[^1]

If you are interested in using it, following these steps should be the fist things to do: 
 
- [ ] To compile and run, this code needs the [GSL](https://www.gnu.org/software/gsl) and [eign3](https://eigen.tuxfamily.org) library. Please make sure they are installed before going forward.
- [ ] To download the software you can run: `git clone ... ` or use the download button. For the installation you need to open the Makefile in `/bin` and set the location of the gsl and eign3 libraries. Then run `make all` from the `/bin` folder. The make command will then first compile the DDCalc library before compiling the source of the software.


[^1]: [Link to my thesis](https://www.theses.fr/en/2021MONTS037)
