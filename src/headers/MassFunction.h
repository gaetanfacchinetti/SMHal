#ifndef DEF_PRESSSCHECHTER
#define DEF_PRESSSCHECHTER

#include "CppLib.h"
#include "Cosmology.h"
#include "PowerSpectrum.h"
#include "GaussLegendre.h"
#include "myspline.h"



template<typename Base, typename T>
inline bool instanceof(const T *ptr) {
   return dynamic_cast<const Base*>(ptr) != nullptr;
}


class MassFunction
{
  public:

    MassFunction(PowerSpectrum power_spectrum) : 
      _power_spectrum(power_spectrum), _cosmo(power_spectrum.get_Cosmology()), _window_type(power_spectrum.get_window()), _delta_c(1.686)
      {
        _extra_params.resize(0); 
        _rho_bar = 1e+9 * _cosmo.cosmological_density(0, Species::MATTER); // [Msol*Mpc^-3]
      };
    MassFunction(PowerSpectrum power_spectrum, std::vector<double> extra_params) : 
    _power_spectrum(power_spectrum), _cosmo(power_spectrum.get_Cosmology()), _window_type(power_spectrum.get_window()), _delta_c(1.686), _extra_params(extra_params){};
    
    virtual ~MassFunction(){};

    /** \brief Bryan and Norman virial overdensity (flat curvature)
     *  \details Taken from http://arxiv.org/abs/astro-ph/9710107
     *  \warning  Only valid here for E(z)*/
    double Overdensity_BryanNorman_1998(double z);

    /** \brief Peak height value in terms of the Mass
     *  \param (double) Mass \f$M~[\rm M_{\odot}]\f$
     *  \param (double) Redshift \f$z\f$ (dimensionless) 
     *  \return (double) Peak height \f$\nu_{M}\f$  (dimensionless)  */
    double nuM(double M, double z) {return _delta_c/_power_spectrum.SigmaM(M, z);};
    
    /** \brief Peak height value in terms of the Radius
    *  \param (double) Radius \f$R~[\rm kpc]\f$
    *  \param (double) Redshift \f$z\f$ (dimensionless) 
    *  \return (double) Peak height \f$\nu_{R}\f$  (dimensionless) */
    double nuR(double R, double z) {return _delta_c/_power_spectrum.SigmaR(R, z);};

   /** \brief Comoving mass function
    *  \details In between \f$M_{\rm min}\f$ and \f$M_{\rm max}\f$ at redshift \f$z\f$ 
    *            and for the Press-Schechter formalism : http://adsabs.harvard.edu/doi/10.1086/152650
    *  \param (double) Mass \f$M~[\rm M_\odot]\f$
    *  \param (double) Redshift\f$z\f$ (dimensionless) 
    *  \return (double) Mass function \f$[\rm M_\odot^{-1}~Mpc^{-3}]\f$ */
    double NumberDensity(double M, double z) {return Der_MassFraction(M, z) * _rho_bar / M;};

    /** \brief Comoving number density of halos 
    *  \details In between \f$M_{\rm min}\f$ and \f$M_{\rm max}\f$ at redshift \f$z\f$ 
    *            and for the Press-Schechter formalism : http://adsabs.harvard.edu/doi/10.1086/152650
    *  \param (double) Minimal mass \f$M_{\rm min}~[\rm M_\odot]\f$
    *  \param (double) Maximal mass \f$M_{\rm max}~[\rm M_\odot]\f$ 
    *  \param (double) Redshift\f$z\f$ (dimensionless) 
    *  \return (double) Number density \f$[\rm Mpc^{-3}]\f$ */
    double NumberHalos_unitsOfVolume(double Mmin, double Mmax, double z);

    // The mass of the Universe contained in Halos
    // Be careful that by default mmin = 1e-14
    double MassInHalos(double z, double Mmin = 1e-15);
    double MassFractionInHalos(double z, double Mmin = 1e-15){return MassInHalos(z, Mmin)/(1e+9*_cosmo.cosmological_density(z, Species::MATTER));};



    // Functions to be evaluated in a new instance
    virtual double Der_MassFraction(double M, double z)  
    {
      std::cout << "FATAL ERROR: cannot call " << __PRETTY_FUNCTION__ << " from an object of the parent class" << std::endl;
      exit(0); 
    };

    virtual double NumberDensity_spike(double As_spike, double ks_spike, double z)
    {
      std::cout << "FATAL ERROR: cannot call " << __PRETTY_FUNCTION__ << " from an object of the parent class" << std::endl;
      exit(0); 
    };
    virtual double Der_MassFraction_spike(double As_spike, double ks_spike, double z)
    {
      std::cout << "FATAL ERROR: cannot call " << __PRETTY_FUNCTION__ << " from an object of the parent class" << std::endl;
      exit(0); 
    };

    


    // General Functions related to the EST formalism
    double NumberProgenitors_EST(double M1, double M2, double z1, double z2);
    double CMF_Density_EST(double M2, double M1, double z); // The cumulative mass function in terms of density
    double fToIntOns_CMF_Density_EST(double S2, double S1, double z);
    double CMF_Approx_EST(double M2, double M1, double z); // Approximative number of structures in the halo of mass M0
    double MF_EST(double M2, double M1, double z1);
    void plot_MF_EST(double M1, double z1);
    void Interpolate_M2_MF_EST_vs_M(double M1, double z1);
    double Interp_M2_MF_EST_vs_M(double M2, double M1, double z1);

    // Setters and getter
    PowerSpectrum get_PowerSpectrum() const { return _power_spectrum; };
    double get_delta_c() const { return _delta_c; };
    Cosmology const get_cosmo() const { return _cosmo; };
    std::vector<double> get_extra_params() const {return _extra_params;};

    // Plotting functions
    void plot_NormalisedHaloMassFunction(double z);
    void plot(double z, bool in_units_of_h = false, std::string add = "");
    void plot_HaloMassFunction_vs_redshift(double M, std::string add);
    void plot_NumberProgenitors_EST(double M1, double z1, double z2);
    void plot_CMF_EST(double M1, double z1);
    void plot_DeltaS_EST(double M1, int version);

    // Static functions (CallBacks)
    static double CallBack_fToIntOns_CMF_Density_EST(void *pt2Object, std::vector<double> xx) { return ((MassFunction *)pt2Object)->fToIntOns_CMF_Density_EST(xx[0], xx[1], xx[2]); };
    static double CallBack_NumberDensity(void *pt2Object, std::vector<double> xx) { return ((MassFunction *)pt2Object)->NumberDensity(xx[0], xx[1]); };
    static double CallBack_NumberDensity_Times_M(void *pt2Object, std::vector<double> xx) { return xx[0]*((MassFunction *)pt2Object)->NumberDensity(xx[0], xx[1]); };

  protected:  
    Cosmology _cosmo;
    PowerSpectrum _power_spectrum;

    PSWindowType _window_type;
    my_spline::spline _spline_M2_MF_EST_log10_M;

    double _delta_c;
    double _rho_bar; // [Msol*Mpc^-3]

    std::vector<double> _extra_params; // extra_parameters of the mass function

};


class MassFunction_PS : public MassFunction
{
public:
  MassFunction_PS(PowerSpectrum power_spectrum) : MassFunction(power_spectrum){};
  virtual ~MassFunction_PS(){};

  virtual double Der_MassFraction(double M, double z);

  /** Comoving number density of halos in between \f$M_{\rm min}\f$ and \f$M_{\rm max}\f$ at redshift \f$z\f$ 
   *  when there is a SPIKE at \f$k_{\rm s}\f$ with amplitude \f$A_{\rm s}\f$ in the matter power spectrum 
   *  -> only relevant for the fourier space top hat window
   *  This function is in factor of a dirac distribution \f$\delta(M- M_{\rm s})\f$ at \f$M_{\rm s}\f$ corresponding to \f$k_{\rm s}\f$
   *  Input  : \f$M\f$ (M\f$_\odot\f$)
   *           \f$z\f$ (dimensionless) 
   *  Output : Number density per unit mass (\f${\rm Mpc}^{-3}/\delta(M-M_{\rm s})\f$) (M\f$_\odot^{-1}\f$) */
  virtual double NumberDensity_spike(double As_spike, double ks_spike, double z);
  virtual double Der_MassFraction_spike(double As_spike, double ks_spike, double z);
};




class MassFunction_ST : public MassFunction
{

  public:

    MassFunction_ST(PowerSpectrum power_spectrum, double small_a = 0.707, double A = 0.3222, double q = 0.3) : 
     MassFunction(power_spectrum), _small_a(small_a), _A(A), _q(q)
      {_extra_params.resize(3);
     _extra_params[0] = _small_a;
     _extra_params[1] = _A;
      _extra_params[2] = _q;};
    MassFunction_ST(PowerSpectrum power_spectrum, std::vector<double> extra_params) : 
      MassFunction(power_spectrum, extra_params), _small_a(extra_params[0]), _A(extra_params[1]), _q(extra_params[2]) {};
    virtual ~MassFunction_ST(){}; 

    virtual double Der_MassFraction(double M, double z);

    virtual double NumberDensity_spike(double As_spike, double ks_spike, double z)
    {
      std::cout << "FATAL ERROR: cannot call " << __PRETTY_FUNCTION__ << " with Sheth and Tormen mass function" << std::endl;
      exit(0); 
    };
    virtual double Der_MassFraction_spike(double As_spike, double ks_spike, double z)
    {
      std::cout << "FATAL ERROR: cannot call " << __PRETTY_FUNCTION__ << " with Sheth and Tormen mass function" << std::endl;
      exit(0); 
    };


  private:
    double _small_a, _A, _q;
  
};

#endif
