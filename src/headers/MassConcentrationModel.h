#ifndef DEF_MASSCONCENTRATIONMODEL
#define DEF_MASSCONCENTRATIONMODEL

#include "CppLib.h"
#include "Cosmology.h"
#include "MassFunction.h"
#include "PowerSpectrum.h"

enum class MassFunctionType
{
  PRESS_SCHECHTER,
  SETH_TORMEN
};

class MassConcentrationModel
{
private:
  MassConcentrationModel(PowerSpectrum power_spectrum, std::shared_ptr<MassFunction> mass_function) : 
    _power_spectrum(power_spectrum), _mass_function(mass_function) {};

public: 
  MassConcentrationModel(){};
  static MassConcentrationModel MassConcentrationModel_from_PowerSpectrum(PowerSpectrum power_spectrum)
  // This function creates an instance of the class MassConcentrationModel with the default Mass Function (Press Schechter)
  {
    std::shared_ptr<MassFunction> mass_function = std::make_shared<MassFunction_PS>(power_spectrum);
    return MassConcentrationModel(power_spectrum, mass_function);
  }

  static MassConcentrationModel MassConcentrationModel_from_MassFunction(std::shared_ptr<MassFunction> mass_function)
  // This function creates an instance of the class MassConcentrationModel from a mass function only
  {
    PowerSpectrum power_spectrum = mass_function->get_PowerSpectrum();
    return MassConcentrationModel(power_spectrum, mass_function);
  }

  ~MassConcentrationModel(){};

  double RedshiftOfCollapse_vs_Mass(double Mass);
  double Mass_vs_RedshiftOfCollapse(double zc);
  double fToDichotomieFor_RedshiftOfCollapse_vs_Mass(std::vector<double> zM);

  double fToIntOnM2_RedshiftOfFormation_LaceyCole(double St, double omFt, double M1, double S1, double Sh);
  void plot_fToIntOnM2_RedshiftOfFormation_LaceyCole(double omFt, double M1);
  double ProbaRedshiftOfFormation_LaceyCole_vs_zF(double zF, double M1, double z1);
  double ProbaRedshiftOfFormation_LaceyCole_vs_omFt(double omFt, double M1);
  void plot_ProbaRedshiftOfFormation_LaceyCole_vs_zF(double M1, double z1);
  void plot_ProbaRedshiftOfFormation_LaceyCole_vs_omFt(double M1);

  double fToIntOnM2_IntegratedProbaRedshiftOfFormation_LaceyCole(double St, double omFt, double M1, double S1, double Sh);
  double IntegratedProbaRedshiftOfFormation_LaceyCole_vs_omFt(double omFt, double M1);
  void plot_IntegratedProbaRedshiftOfFormation_LaceyCole_vs_omFt(double M1);

  // norm = true if normalisation to 1 on positive redshifts only
  // norm = false if normalisation to 1 on redshifts on [-1, \infty]
  double ProbaRedshiftOfCollapse(double zc, double M, bool norm = true);
  double CumulativeRedshiftOfCollapse(double zc, double M);
  double MeanRedshiftOfCollapse(double M);

  double fToIntOnM_for_dndzc(double zc, double M, double z);
  double dndzc(double zc, double z);
  void plot_dndzc(double z);


  // Concentration models
  double ConcentrationMaccio(double zc, double z);
  double ConcentrationMaccioMass(double z, double M);
  double MeanConcentrationMaccio(double z, double M);
  double SigmaConcentrationMaccio(double z, double M);

  double ConcentrationOkoli(double z, double m200);
  double A2_Okoli(double z, double R);
  double B2_Okoli(double z, double R);
  double AB_Okoli(double z, double R);
  double Window_B_Okoli(double kR);
  double fToIntFor_AB_Okoli(std::vector<double> lnkRz);
  double fToIntFor_B2_Okoli(std::vector<double> lnkRz);

  double ConcentrationFacchinetti(double zc, double z);
  double ConcentrationFacchinetti_vs_m(double m, double z);
  double RedshiftOfCollapse_Facchinetti(double c, double m, double z);
  double fForBissection_ConcentrationFacchinetti(double c, double z, double zc);
  double fForBissection_RedshiftOfCollapse_Facchinetti(double zc, double c, double m, double z);
  double fForBissection_ConcentrationFacchinetti_vs_m(double c, double m, double z);

  double ConcentrationSanchezCondePrada(double m200);



  // Plot functions
  void plot_ProbaRedshiftOfCollapse_vs_zc(double M);
  void plot_ProbaRedshiftOfCollapse_vs_M(double zc);
  void plot_CumulativeRedshiftOfCollapse_vs_zc(double M);
  void plot_concentration(double z);
  void plot_B_over_A_Okoli(double z);
  void plot_fToIntFor_AB_Okoli(double z, double R);
  //void plot_ProbaRedshiftOfCollapse_vs_k(double zc);

  static double CallBack_fToDichotomieFor_RedshiftOfCollapse_vs_Mass(void *pt2Object, std::vector<double> xx) { return ((MassConcentrationModel *)pt2Object)->fToDichotomieFor_RedshiftOfCollapse_vs_Mass(xx); };

  static double CallBack_fForBissection_ConcentrationFacchinetti(void *pt2Object, std::vector<double> xx) { return ((MassConcentrationModel *)pt2Object)->fForBissection_ConcentrationFacchinetti(xx[0], xx[1], xx[2]); };
  static double CallBack_fForBissection_RedshiftOfCollapse_Facchinetti(void *pt2Object, std::vector<double> xx) { return ((MassConcentrationModel *)pt2Object)->fForBissection_RedshiftOfCollapse_Facchinetti(xx[0], xx[1], xx[2], xx[3]); };
  static double CallBack_fForBissection_ConcentrationFacchinetti_vs_m(void *pt2Object, std::vector<double> xx) { return ((MassConcentrationModel *)pt2Object)->fForBissection_ConcentrationFacchinetti_vs_m(xx[0], xx[1], xx[2]); };

  static double CallBack_fToIntFor_AB_Okoli(void *pt2Object, std::vector<double> xx) { return ((MassConcentrationModel *)pt2Object)->fToIntFor_AB_Okoli(xx); };
  static double CallBack_fToIntFor_B2_Okoli(void *pt2Object, std::vector<double> xx) { return ((MassConcentrationModel *)pt2Object)->fToIntFor_B2_Okoli(xx); };
  static double CallBack_fToIntFor_NormalisationProbazc(void *pt2Object, std::vector<double> zcM) { return ((MassConcentrationModel *)pt2Object)->ProbaRedshiftOfCollapse(zcM[0], zcM[1]); };
  static double CallBack_fToIntFor_MomentProbazc(void *pt2Object, std::vector<double> zcM) { return pow(zcM[0], zcM[2]) * ((MassConcentrationModel *)pt2Object)->ProbaRedshiftOfCollapse(zcM[0], zcM[1]); };
  static double CallBack_fToIntFor_MeanConcentrationMaccio(void *pt2Object, std::vector<double> zcMz) 
    { return ((MassConcentrationModel *)pt2Object)->ConcentrationMaccio(zcMz[0], zcMz[2]) * ((MassConcentrationModel *)pt2Object)->ProbaRedshiftOfFormation_LaceyCole_vs_zF(zcMz[0], zcMz[1], 0); };
  static double CallBack_fToIntFor_Log10ConcentrationNMaccio(void *pt2Object, std::vector<double> zcMzn) 
    { return pow(log10(((MassConcentrationModel *)pt2Object)->ConcentrationMaccio(zcMzn[0], zcMzn[2])), zcMzn[3]) * ((MassConcentrationModel *)pt2Object)->ProbaRedshiftOfFormation_LaceyCole_vs_zF(zcMzn[0], zcMzn[1], 0); };
 
  static double CallBack_fToIntFor_Concentration2Maccio(void *pt2Object, std::vector<double> zcMz) 
    { return pow(((MassConcentrationModel *)pt2Object)->ConcentrationMaccio(zcMz[0], zcMz[2]), 2) * ((MassConcentrationModel *)pt2Object)->ProbaRedshiftOfCollapse(zcMz[0], zcMz[1]); };
  
  static double CallBack_fToIntOnM2_RedshiftOfFormation_LaceyCole(void *pt2Object, std::vector<double> xx) { return ((MassConcentrationModel *)pt2Object)->fToIntOnM2_RedshiftOfFormation_LaceyCole(xx[0], xx[1], xx[2], xx[3], xx[4]); };
  static double CallBack_fToIntOnM2_IntegratedProbaRedshiftOfFormation_LaceyCole(void *pt2Object, std::vector<double> xx) 
    { return ((MassConcentrationModel *)pt2Object)->fToIntOnM2_IntegratedProbaRedshiftOfFormation_LaceyCole(xx[0], xx[1], xx[2], xx[3], xx[4]); };

 static double CallBack_fToIntOnM_for_dndzc(void *pt2Object, std::vector<double> xx) { return ((MassConcentrationModel *)pt2Object)->fToIntOnM_for_dndzc(xx[0], xx[1], xx[2]); };

private:
  Cosmology _cosmo;
  std::shared_ptr<MassFunction> _mass_function;
  TransferFunction_EH98 _transfer_function;
  PowerSpectrum _power_spectrum;
};

#endif