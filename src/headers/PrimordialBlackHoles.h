#ifndef PRIMORDIALBLACKHOLES
#define PRIMORDIALBLACKHOLES

#include "CppLib.h"
#include "DegreeFreedom.h"
#include "PowerSpectrum.h"

class PrimordialBlackHoles
{

public:
    PrimordialBlackHoles(PowerSpectrum n_ps, DegreeFreedom *n_degFree) : ps(n_ps), degFree(n_degFree){};
    ~PrimordialBlackHoles(){};

    double MassHorizon_vs_Tformation(double Tform); // Result in Msol
    double Tformation_vs_MassHorizon(double Mass); // Result in GeV for input in Msol
    double Radius_vs_MassHorizon(double Mass);  // Result in Mpc for input in Msol
    double MassHorizon_vs_Radius(double Radius);  // Result in Msol for input in Mpc 

    // Function for a spiked  curvature power spectrum P_R(k) = A_s k_s \delta(k-k_s)
    double SigmaR(double R, double As, double ks);
    double TransferFunction_RDera(double y);
    double Beta(double R, double alpha, double As, double ks);
    double fPBH(double Mass, double alpha, double As, double ks);
    double Max_fPBH(double alpha, double As, double ks);
    double Max_As(double max_f, double alpha, double ks);

    void plot_MassHorizon_vs_Radius();
    void plot_fPBH(double As, double ks);
    void plot_Max_As();


    static double CallBack_MassHorizon_vs_Tformation_toSolve(void *pt2Object, std::vector<double> xx) { return ((PrimordialBlackHoles *)pt2Object)->MassHorizon_vs_Tformation(xx[0]) - xx[1]; };
    static double CallBack_Radius_vs_MassHorizon_toSolve(void *pt2Object, std::vector<double> xx) { return ((PrimordialBlackHoles *)pt2Object)->Radius_vs_MassHorizon(xx[0]) - xx[1]; };
    static double CallBack_Max_fPBH_toSolve(void *pt2Object, std::vector<double> xx) { return ((PrimordialBlackHoles *)pt2Object)->Max_fPBH(xx[2], xx[0], xx[3]) - xx[1]; };

private :
    PowerSpectrum ps;
    DegreeFreedom *degFree;

};

#endif