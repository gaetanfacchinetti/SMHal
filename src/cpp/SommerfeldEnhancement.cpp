#include "../headers/SommerfeldEnhancement.h"


SommerfeldEnhancement::SommerfeldEnhancement()
{
    std::vector<double> _ephi_peak_arr;

    _ephi_res_arr.resize(0);
    _ephi_sat_arr.resize(0);
    _ephi_peak_arr.resize(0);
    _ephi_arr.resize(0);


    int npts = 75, mpts = 15, n =0;
    double temp = 0, dephi = 0;
    
    // Construct ephi_res_arr and ephi_sat_arr
    for (int i = 0; i < npts; i++)
    {
        n = npts-i;
        _ephi_res_arr.push_back(6./PI/PI/n/n);
        _ephi_sat_arr.push_back(24./PI/PI/pow(2*n+1, 2));
        
       
        _ephi_peak_arr.push_back(24./PI/PI/pow(2*n+1, 2));
        _ephi_peak_arr.push_back(6./PI/PI/n/n);

        //std::cout << "We are here: " << _ephi_res_arr[i] << std::endl;
    }

    double ephi_min = 6./PI/PI/npts/npts;


    //std::cout  << "Size: " << _ephi_res_arr.size()-1 << std::endl;

    int ii = 0;
    for (int i = 0; i < _ephi_peak_arr.size()-1; i++)
    {
        double ephi_temp = 0;
        int k = 0;

        while (ephi_temp < _ephi_peak_arr[i+1])
        {
            ephi_temp = _ephi_peak_arr[i]*pow(10, k*5e-2);

            if(ephi_temp < _ephi_peak_arr[i+1])
            {
                _ephi_arr.push_back( ephi_temp );
            
                //std::cout << ii << "\t" << k << "\t" << _ephi_arr[ii] << "\t" << _ephi_peak_arr[i] << "\t" << _ephi_peak_arr[i+1] << std::endl;
                //ii += 1;
                k += 1;
            }
        }


    }


    dephi = log10( 100./_ephi_peak_arr[_ephi_peak_arr.size()-1])/25.;
    for (int j = 0; j < 25+1; j++)
        _ephi_arr.push_back( _ephi_res_arr[npts-1]*pow(10, j*dephi) );


    //for (int k = 0; k < _ephi_peak_arr.size(); k++)
    //    std::cout << _ephi_peak_arr[k] << " " << _ephi_peak_arr.size() << std::endl;

    //for (int k = 0; k < _ephi_arr.size(); k++)
    //    std::cout << _ephi_arr[k] << " " << _ephi_arr.size() << std::endl;

}

double SommerfeldEnhancement::Sommerfeld_Hulthen_epsv(double ev, double ephi, double l, double alpha)
{
    double evr = ev + pow(alpha, 3);
    double ephis = PI*PI/6.*ephi;
    double S0 = 0;


    if (evr < sqrt(ephis))
        S0 = PI/evr * sinh(2.*PI*evr/ephis) / (cosh(2*PI*evr/ephis) - cos(2*PI*sqrt(1/ephis - pow(evr/ephis, 2))));
    else
        S0 = PI/evr/(1-exp(-PI/evr));

    if (l == 0)
        return S0;
    
    if (l == 1)
        return S0*(pow(1.-ephis, 2) + 4*pow(evr, 2))/(pow(ephis, 2) + 4*pow(evr, 2));

    std::cout << "FATAL ERROR: Trying to access an unknwo value of l in " << __PRETTY_FUNCTION__ << std::endl; 

    exit(0);
}

void SommerfeldEnhancement::plot(double ev, double alpha)
{
    std::ofstream outfile;
    outfile.open("../output/CosmologicalBoost/SommerfeldBoost_epsv" + std::to_string(ev) + ".out");
    outfile.precision(8);

    double ephi =0;

    for (int i = 0; i < _ephi_arr.size(); i++)
        outfile << _ephi_arr[i] << "\t" << Sommerfeld_Hulthen_epsv(ev, _ephi_arr[i], 0, alpha) << "\t" << Sommerfeld_Hulthen_epsv(ev, _ephi_arr[i], 1, alpha) << std::endl;
    
    outfile.close();
}


void SommerfeldEnhancement::plot_vs_v(double ephi, double alpha)
{
    std::ofstream outfile;
    outfile.open("../output/CosmologicalBoost/SommerfeldBoost_ephi" + std::to_string(ephi) + ".out");
    outfile.precision(8);

    int Npts = 100;
    double ev_sat = ephi;
    double ev_min = 1e-6*ev_sat, ev_max = std::min(1., 1e+3*ev_sat);
    double vmin = alpha*ev_min, vmax = std::min(1., alpha*ev_max);

    double dv = log10(vmax/vmin)/(1.*Npts), v = 0;

    for (int i = 0; i < Npts+1; i++)
    {
        v = vmin*pow(10, i*dv);
        outfile << v*1e-3*C_LIGHT << "\t" << Sommerfeld_Hulthen_v(v, ephi, 0, alpha) << "\t" << v*v*Sommerfeld_Hulthen_v(v, ephi, 1, alpha) << "\t" << v*v << std::endl;
    }

    outfile.close();
}
