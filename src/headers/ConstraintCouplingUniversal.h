#ifndef DEF_CONSTRAINTCOUPLINGUNIVERSAL
#define DEF_CONSTRAINTCOUPLINGUNIVERSAL

#include "CrossSections/CrossSection_ccgaga.h"
#include "CrossSections/CrossSection_ccff.h"
#include "CrossSections/CrossSection_ccss.h"
#include "CrossSections/CrossSection_ccpp.h"
#include "CrossSections/CrossSection_ccsp.h"
#include "CrossSections/CrossSectionMajorana_cccc.h"
#include "CrossSections/CrossSectionDirac_cccc.h"
#include "CrossSections/CrossSectionDirac_ccbccb.h"
#include "CrossSections/CrossSection_cfcf.h"
#include "SpecialIntegrations.h"
#include "StandardModel.h"
#include "DegreeFreedom.h"
#include "ChemicalDecoupling.h"
#include "KineticDecoupling.h"
#include "MinimalMass.h"
#include "CPOddHiggsLowMass.h"
#include "InputReader.h"
#include "DecayRates.h"
#include "gsl/gsl_complex.h"
#include "gsl/gsl_complex_math.h"
#include "ExceptionHandler.h"
#include "FSLModel.h"
#include <ctime>
#include "DirectDetection.h"
//#include "DirectDetection.h"

enum class CouplingType
{
    UNIVERSAL_CONSTANT,
    UNIVERSAL_YUKAWA,
    SINGLE_CONSTANT,
    SINGLE_YUKAWA
};


class ConstraintCouplingUniversal
{
public:
    ConstraintCouplingUniversal(CouplingType n_cptype = CouplingType::UNIVERSAL_CONSTANT, Cosmology n_Cosmo = Cosmology(CosmoModel::PLANCKONLY)){n_SM = StandardModel::getInstance()->get_n_couples_SM_ferm(); cptype = n_cptype; Cosmo = n_Cosmo;};
    ~ConstraintCouplingUniversal(){};

    double Coupling_correctAbundance(std::string const &name_input, std::string const &name_output);
    double TkdAndMinMass(std::string const &name_model_file_input, std::string const &name_coupling_file_input, std::string const &name_output);
    double UnevolvedNsubhalos(std::string const &name_minmass_file_input, std::string const &name_output);
    double SelfInteractingCS(std::string const &name_model_file_input, std::string const &name_coupling_file_input, std::string const &name_output);
    double DecayTime(std::string const &name_model_file_input, std::string const &name_coupling_file_input, std::string const &name_output);

    std::vector<double> ConstraintOneModel(DarkSectorModel &model, DegreeFreedom *degFree, double _omega_C_h2);
    std::vector<double> TemperatureKdOneModel(DarkSectorModel &model, double lambda, DegreeFreedom *degFree);
    bool SetNewUniversalCoupling(DarkSectorModel &model, double lambda_SM, double lambda_DM);
    bool SetNewUniversalYukawaCoupling(DarkSectorModel &model, double lambda_SM, double lambda_DM);
    bool SetNewSingleCoupling(DarkSectorModel &model, double lambda_SM, double lambda_DM, int coupl);
    bool SetCoupling(DarkSectorModel &model, double lambda_SM, double lambda_DM, int coupl=3); // By default coupling to electron/positron only if nothing specified
    void Write_input_MinimalMassHalos_vs_massDM_vs_massProp(std::string const &name_input, Proptype type, int Npts, double massMin, double massMax);

    double MaximalValueCouplingDM(DarkSectorModel &model, int coupl = 3);

    // getter funcitons
    CouplingType get_coupling_type() const {return cptype;}; 

    // Reader functions
    std::vector<std::vector<double>> ReadCouplings(std::string const &name);
    std::vector<std::vector<double> > ReadMinimalMass(std::string const &name);
    
    double WriteMassesFromModel(std::string const &name_input);
    double plot_approximate_Yf();

private:
    int n_SM;  
    Cosmology Cosmo;
    CouplingType cptype;
};

#endif