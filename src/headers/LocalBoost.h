#ifndef DEF_LOCALBOOST
#define DEF_LOCALBOOST

#include "CppLib.h"
#include "FSLModel.h"
#include "DarkHalo.h"

class LocalBoost
{
    public:
    LocalBoost(FSLModel n_fslmod) : fslmod(n_fslmod), subModel(fslmod.get_subModel()){};
    ~LocalBoost(){};

    double Luminosity_OneSubhalo(double m200, double c200, double R);
    double DimensionlessLuminosity(double xt);

    // Luminosity of subhalos
    double fToIntOnc200_LuminositySubhalos(double c200, double m200, double R);
    double fToIntOnm200_LuminositySubhalos(double m200, double R);
    double LuminositySubhalos(double R);

    // Lumminosity of smooth and cross term
    double LuminositySmooth(double R);
    double LuminositySmoothSubhalos(double R);

    double TotalLuminosityWithSubhalos(double R);
    double TotalLuminosityWithoutSubhalos(double R);

    double Boost(double R);
    double TotalBoost();

    // CallBack functions
    static double CallBack_fToIntOnc200_LuminositySubhalos(void *pt2Object, std::vector<double> xx) { return ((LocalBoost *)pt2Object)->fToIntOnc200_LuminositySubhalos(xx[0], xx[1], xx[2]); };
    static double CallBack_fToIntOnm200_LuminositySubhalos(void *pt2Object, std::vector<double> xx) { return ((LocalBoost *)pt2Object)->fToIntOnm200_LuminositySubhalos(xx[0], xx[1]); };
    static double CallBack_ftoIntOnR_TotalLuminosityWithSubhalos(void *pt2Object, std::vector<double> xx) { return ((LocalBoost *)pt2Object)->TotalLuminosityWithSubhalos(xx[0]); };
    static double CallBack_ftoIntOnR_TotalLuminosityWithoutSubhalos(void *pt2Object, std::vector<double> xx) { return ((LocalBoost *)pt2Object)->TotalLuminosityWithoutSubhalos(xx[0]); };

    private:
    FSLModel fslmod;
    DarkHalo subModel;


};

#endif