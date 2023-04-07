#ifndef DEF_COSMOLOGICALBOOST
#define DEF_COSMOLOGICALBOOST

#include "MassConcentrationModel.h"
#include "ExceptionHandler.h"
#include "DarkHalo.h"
#include "SommerfeldEnhancement.h"


enum class ProfileConfig
{
    ALL_NFW,
    ALL_MOORE,
    MIXED
};

class CosmologicalBoost 
{
    public:
    CosmologicalBoost(std::shared_ptr<MassFunction> mass_function);
    ~CosmologicalBoost(){};

    /** \brief Boost factor for a single halo
    *  \param (double) Mass \f$m~[\rm M_{\odot}]\f$
    *  \param (double) Redshift of formation \f$z_{\rm f}\f$ (dimensionless) 
    *  \param (double) Redshift z of observation \f$z\f$ (dimensionless) 
    *  \return (DensityProfile) Profile of the halo  */
    double OneHaloBoost(double m, double zf, double z, DensityProfile profile = DensityProfile::NFW);

    double fToIntOnm_for_TotalBoost(double m, double z, double n_prof, double distrib_model);
    double fToIntOnzf_for_TotalBoost(double m, double zf, double z, double n_prof, double distrib_model);
    double HaloDistribution(double m, double zf, double z, double distrib_model);
    double TotalBoost(double z, double mmin, ProfileConfig prof_config, double distrib_model);
    double MinimalBound_TotalBoost(double z, double mmin);
    
    double AverageRedshiftOfCollapse_OnSpike(double z);
    double fToIntOnmu_AverageRedshiftOfCollapse_OnSpike(double nu, double z);

    /** \brief Formation reshift of a halo of mass m at z (simple model)
    *  \param (double) Mass \f$m~[\rm M_{\odot}]\f$
    *  \param (double) Redshift z of observation \f$z\f$ (dimensionless) */
    double SimpleModel_FormationRedshift_Delta(double m, double z);

    double BBKSFunction(double x);

    double MassFractionInSpike(double mmin, double z);

    /** \brief function to integrate over in mass to get the approximate sommerfeld enhanced boost factor */
    double fToIntOnm_for_TotalBoost_with_Sommerfeld_approx(double m, double z, double ephi, double l, double alpha);
    double TotalBoost_with_Sommerfeld_approx(double z, double ephi, double l, double alpha, double mmin);
    void plot_TotalBoost_with_Sommerfeld_approx(double l, double ephi, double alpha, double mmin, std::string add_file_name="");

    void plot_TotalBoost(double mmin, double distrib_config, std::string add_file_name = "");
    void plot_FormationRedshiftDistribution(std::string name, double distrib_model);

    static double CallBack_fToIntOnm_for_TotalBoost(void *pt2Object, std::vector<double> xx) { return ((CosmologicalBoost *)pt2Object)->fToIntOnm_for_TotalBoost(xx[0], xx[1], xx[2], xx[3]); };
    static double CallBack_fToIntOnmu_AverageRedshiftOfCollapse_OnSpike(void *pt2Object, std::vector<double> xx) { return ((CosmologicalBoost *)pt2Object)->fToIntOnmu_AverageRedshiftOfCollapse_OnSpike(xx[0], xx[1]); };
    static double CallBack_fToIntOnzf_for_TotalBoost(void *pt2Object, std::vector<double> xx) { return ((CosmologicalBoost *)pt2Object)->fToIntOnzf_for_TotalBoost(xx[0], xx[1], xx[2], xx[3], xx[4]); };
    static double CallBack_fToIntOnm_for_TotalBoost_with_Sommerfeld_approx(void *pt2Object, std::vector<double> xx) { return ((CosmologicalBoost *)pt2Object)->fToIntOnm_for_TotalBoost_with_Sommerfeld_approx(xx[0], xx[1], xx[2], xx[3], xx[4]); };

    private:
    PowerSpectrum _power_spectrum, _power_spectrum_without_spikes;
    MassConcentrationModel _mass_concentration_model, _mass_concentration_model_without_spikes;
    std::shared_ptr<MassFunction> _mass_function;
    SommerfeldEnhancement _sommerfeld_enhancement;

    double _rhom0; // [Msol * Mpc^{-3}]

};



#endif