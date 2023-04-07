#ifndef DARKSECTORMODEL
#define DARKSECTORMODEL

#include "Particle.h"
#include "StandardModel.h"


// Attention here a_VEC are the vector coupliens b_VEC the axial vector ones
class DarkSectorModel
{

public:
  DarkSectorModel();
  ~DarkSectorModel() {};

  // Getter
  Particle get_chi(int i) const { return chi[i]; };
  Particle get_phis(int i) const { return phis[i]; };
  Particle get_phip(int i) const { return phip[i]; };
  Particle get_X(int i) const { return X[i]; };

  Particle *get_chi_ptr(int i) {return &chi[i];};
  Particle *get_phis_ptr(int i) {return &phis[i];};
  Particle *get_phip_ptr(int i) {return &phip[i];};
  Particle *get_X_ptr(int i) {return &X[i];};

  std::vector<Particle> get_chi() const { return chi; };
  std::vector<Particle> get_phis() const { return phis; };
  std::vector<Particle> get_phip() const { return phip; };
  std::vector<Particle> get_X() const { return X; };
  std::vector<Particle> get_DS_Part() const { return DS_Part; };

  int get_n_DS_DM() const { return n_DS_DM; };
  int get_n_DS_phis() const { return n_DS_phis; };
  int get_n_DS_phip() const { return n_DS_phip; };
  int get_n_DS_X() const { return n_DS_X; };

  double get_lambda_S_SM(int prop_s, int couple) const { return lambda_S_SM[prop_s][couple]; };
  std::vector<double> get_lambda_S_SM(int prop_s) const { return lambda_S_SM[prop_s]; };
  std::vector<std::vector<double>> get_lambda_S_SM() const { return lambda_S_SM; };
  double get_lambda_PS_SM(int prop_ps, int couple) const { return lambda_PS_SM[prop_ps][couple]; };
  std::vector<double> get_lambda_PS_SM(int prop_p) const { return lambda_PS_SM[prop_p]; };
  double get_a_VEC_SM(int prop_vec, int couple) const { return a_VEC_SM[prop_vec][couple]; };
  std::vector<double> get_a_VEC_SM(int prop_vec) const { return a_VEC_SM[prop_vec]; };
  double get_b_VEC_SM(int prop_vec, int couple) const { return b_VEC_SM[prop_vec][couple]; };
  std::vector<double> get_b_VEC_SM(int prop_vec) const { return b_VEC_SM[prop_vec]; };
 
  double get_lambda_S_DM(int prop_s, int dm_1, int dm_2) const { return lambda_S_DM[prop_s][dm_1][dm_2]; };
  std::vector<std::vector<double>> get_lambda_S_DM(int prop_s) const { return lambda_S_DM[prop_s]; };
  std::vector<std::vector<std::vector<double>>> get_lambda_S_DM() const { return lambda_S_DM; };
  double get_lambda_PS_DM(int prop_ps, int dm_1, int dm_2) const { return lambda_PS_DM[prop_ps][dm_1][dm_2]; };
  std::vector<std::vector<double>> get_lambda_PS_DM(int prop_p) const { return lambda_PS_DM[prop_p]; };
  std::vector<std::vector<std::vector<double>>> get_lambda_PS_DM() const { return lambda_PS_DM; };
  
  /// Vector coupling between the mediator and the DM particles return (double) input (id_med, id_dm1, id_dm2)
  double get_a_VEC_DM(int prop_vec, int dm_1, int dm_2) const { return a_VEC_DM[prop_vec][dm_1][dm_2]; };
  /// Vector coupling between the mediator and the DM particles return (std::vector<std::vector<double>>) input (id_med)
  std::vector<std::vector<double>> get_a_VEC_DM(int prop_vec) const { return a_VEC_DM[prop_vec]; };
  /// Vector coupling between the mediator and the DM particles return (std::vector<std::vector<std::vector<double>>>) input ()
  std::vector<std::vector<std::vector<double>>> get_a_VEC_DM() const { return a_VEC_DM; };

  /// Axial-Vector coupling between the mediator and the DM particles return (double) input (id_med, id_dm1, id_dm2)
  double get_b_VEC_DM(int prop_vec, int dm_1, int dm_2) const { return b_VEC_DM[prop_vec][dm_1][dm_2]; };
  /// Axial-Vector coupling between the mediator and the DM particles return (std::vector<std::vector<double>>) input (id_med)
  std::vector<std::vector<double>> get_b_VEC_DM(int prop_vec) const { return b_VEC_DM[prop_vec]; };
   /// Axial-Vector coupling between the mediator and the DM particles return (std::vector<std::vector<std::vector<double>>>) input ()
  std::vector<std::vector<std::vector<double>>> get_b_VEC_DM() const { return b_VEC_DM; };
  
  double get_c_sss(int prop_scalar_1, int prop_scalar_2, int prop_scalar_3) const { return c_sss[prop_scalar_1][prop_scalar_2][prop_scalar_3]; };
  std::vector<std::vector<double>> get_c_sss(int prop_scalar_1) const { return c_sss[prop_scalar_1]; };
  std::vector<std::vector<std::vector<double>>> get_c_sss() const { return c_sss; };
  double get_d_spp(int prop_scalar_1, int prop_pseudoscalar_1, int prop_pseudoscalar_2) const { return d_spp[prop_scalar_1][prop_pseudoscalar_1][prop_pseudoscalar_2]; };
  std::vector<std::vector<double>> get_d_spp(int prop_scalar_1) const { return d_spp[prop_scalar_1]; };
  std::vector<std::vector<std::vector<double>>> get_d_spp() const { return d_spp; };
  double get_d_psp(int prop_pseudoscalar_1, int prop_scalar_1, int prop_pseudoscalar_2) const { return d_spp[prop_pseudoscalar_1][prop_scalar_1][prop_pseudoscalar_2]; };
  std::vector<std::vector<double>> get_d_psp(int prop_pseudoscalar_1) const { return d_psp[prop_pseudoscalar_1]; };
  std::vector<std::vector<std::vector<double>>> get_d_psp() const { return d_psp; };
  double get_g_sXX(int prop_scalar_1, int prop_vector_1, int prop_vector_2) const { return g_sXX[prop_scalar_1][prop_vector_1][prop_vector_2]; };
  std::vector<std::vector<double>> get_g_sXX(int prop_scalar_1) const { return g_sXX[prop_scalar_1]; };


  //std::vector<int> get_indexToDSChiIndex() const { return indexToDSChiIndex; };

  void add_propagator(Proptype type, int dof, std::string name, double mass, double width);
  void add_darkmatter(double mass, double width, Fermiontype ftype = Fermiontype::majorana); // DM particles are set to Majorana by default

  void set_coupling_DS_S_SM(int prop_s, int couple, double const& lambda){lambda_S_SM[prop_s][couple] = lambda;};
  void set_coupling_DS_PS_SM(int prop_p, int couple, double const& lambda){lambda_PS_SM[prop_p][couple] = lambda;};
  void set_coupling_DS_a_VEC_SM(int prop_vec, int couple, double const& lambda){a_VEC_SM[prop_vec][couple] = lambda; };
  void set_coupling_DS_b_VEC_SM(int prop_vec, int couple, double const& lambda){b_VEC_SM[prop_vec][couple] = lambda; };

  void set_coupling_DS_S_DM(int prop_s, int dm_1, int dm_2, double const& lambda){lambda_S_DM[prop_s][dm_1][dm_2] = lambda;};
  void set_coupling_DS_PS_DM(int prop_p, int dm_1, int dm_2, double const& lambda){lambda_PS_DM[prop_p][dm_1][dm_2] = lambda;};
  void set_coupling_DS_a_VEC_DM(int prop_vec, int dm_1, int dm_2, double const& lambda){a_VEC_DM[prop_vec][dm_1][dm_2] = lambda; };
  void set_coupling_DS_b_VEC_DM(int prop_vec, int dm_1, int dm_2, double const& lambda){b_VEC_DM[prop_vec][dm_1][dm_2] = lambda; };

  void set_couplings_DS_S_SM(std::vector<std::vector<double>> const& n_lambda_S_SM) { lambda_S_SM = n_lambda_S_SM; };
  void set_couplings_DS_PS_SM(std::vector<std::vector<double>> const& n_lambda_PS_SM) { lambda_PS_SM = n_lambda_PS_SM; };
  void set_couplings_DS_VEC_SM(std::vector<std::vector<double>> const&  n_a_VEC_SM, std::vector<std::vector<double>> n_b_VEC_SM)
  {
    a_VEC_SM = n_a_VEC_SM;
    b_VEC_SM = n_b_VEC_SM;
  };
  void set_couplings_DS_S_DM(std::vector<std::vector<std::vector<double>>> const& n_lambda_S_DM) { lambda_S_DM = n_lambda_S_DM; };
  void set_couplings_DS_PS_DM(std::vector<std::vector<std::vector<double>>> const& n_lambda_PS_DM) { lambda_PS_DM = n_lambda_PS_DM; };
  void set_couplings_DS_VEC_DM(std::vector<std::vector<std::vector<double>>> const& n_a_VEC_DM, std::vector<std::vector<std::vector<double>>> n_b_VEC_DM)
  {
    a_VEC_DM = n_a_VEC_DM;
    b_VEC_DM = n_b_VEC_DM;
  };
  void set_couplings_c_sss(std::vector<std::vector<std::vector<double>>> const& n_c_sss) { c_sss = n_c_sss; };
  void set_couplings_g_sXX(std::vector<std::vector<std::vector<double>>> const& n_g_sXX) { g_sXX = n_g_sXX; };
  void set_couplings_d_spp(std::vector<std::vector<std::vector<double>>> const& n_d_spp) 
  { 
    d_spp = n_d_spp; 
    
    for (int i = 0; i < d_spp.size(); i++)
      for (int j = 0; j < d_spp[i].size(); j++)
        for (int k = 0; k < d_spp[i][j].size(); k++)
            d_psp[j][i][k] = d_spp[i][j][k];
  };

  //void set_width_phis(int prop_s, n_width) { phis[0]->set_width(n_width);}

  // Private functions
  void Initialise();
  void ForceInitCouplings();

private:
  // List of couplings
  //double cSXX, cSSS, cSPP;
  std::vector<double> cSff, cPff, AXff, BXff;
  std::vector<std::vector<double>> cSCC, cPCC, AXCC, BXCC;

  // Additional vectors and constants
  int n_DS_DM, n_DS_phis, n_DS_phip, n_DS_X;

  std::vector<Particle> chi;
  std::vector<Particle> phis;
  std::vector<Particle> phip;
  std::vector<Particle> X;

  std::vector<bool> isPresentBeforeQCDPT, isPresentAfterQCDPT;
  std::vector<Particle> DS_Elem_Part, DS_Prop, DS_Part;

  //std::vector<int> indexToSMFermIndex;
  //std::vector<int> indexToDSChiIndex;

  // Coupling tables
  std::vector<std::vector<double>> lambda_S_SM, lambda_PS_SM, a_VEC_SM, b_VEC_SM;
  std::vector<std::vector<std::vector<double>>> lambda_S_DM, lambda_PS_DM, a_VEC_DM, b_VEC_DM;
  std::vector<std::vector<std::vector<double>>> c_sss, d_spp, d_psp, g_sXX;

  int n_SM_couples;
  int n_SM_part;
  int part_DS_index;
};

#endif
