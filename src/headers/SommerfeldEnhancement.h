#ifndef SOMMERFELDENHANCEMENT
#define SOMMERFELDENHANCEMENT

#include "CppLib.h"
#include "mymath.h"
#include "MyUnits.h"

class SommerfeldEnhancement
{

    public:
        SommerfeldEnhancement();
        ~SommerfeldEnhancement(){};

        std::vector<double> get_ephi_arr() const {return _ephi_arr;};
        std::vector<double> get_ephi_res_arr() const {return _ephi_res_arr;};
        std::vector<double> get_ephi_sat_arr() const {return _ephi_sat_arr;};

    ///

        double Sommerfeld_Hulthen_epsv(double ev, double ephi, double l, double alpha);
        double Sommerfeld_Hulthen_v(double v, double ephi, double l, double alpha) {return Sommerfeld_Hulthen_epsv(v/alpha, ephi, l, alpha);};

        void plot(double ev, double alpha);
        void plot_vs_v(double ephi, double alpha);

    private:
        std::vector<double> _ephi_res_arr, _ephi_sat_arr, _ephi_arr;

};


#endif