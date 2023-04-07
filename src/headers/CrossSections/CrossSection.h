/**
* This class is one of the most important of the code.  
* It defines all the function that cross-section subclasses inherit.  
* It allows to compute the standard cross-section or
* the transfer cross-section associated to any process.
*/

#ifndef CROSSSECTION
#define CROSSSECTION

#include "../Particle.h"
#include "../SpecialIntegrations.h"
#include "../StandardModel.h"
#include "../DarkSectorModel.h"


class CrossSection
{
public:
  /// Constructor for interactions \f$p_1 + p_2 \rightarrow p_3 + p_4\f$ 
  CrossSection(Particle const &p_1, Particle const &p_2, Particle const &p_3, Particle const &p_4);
  virtual ~CrossSection(){};

  // =============================
  // Main functions
  // =============================

  /// Standard cross-section multiplied by \f$p_1^2\f$ (in the center of mass)
  double EvalSigmaP1cm2(double s);
  /// Standard cross-section multiplied by \f$p_1^2\f$ (in the center of mass) evaluated from numerical integration
  double EvalSigmaP1cm2_fromIntegral(double s);
  /// Standard cross-section
  double EvalSigma(double s) { return EvalSigmaP1cm2(s) * pow(f_p1cm(s), -2); };
  /// Standard cross-section evaluated from numerical integration
  double EvalSigma_fromIntegral(double s) { return EvalSigmaP1cm2_fromIntegral(s) * pow(f_p1cm(s), -2); };
  /// Transfer cross-section multiplied by \f$p_1^4\f$ (in the center of mass)
  double EvalSigmaTransferP1cm4(double s);
  /// Transfer cross-section
  double EvalSigmaTransfer(double s) { return EvalSigmaTransferP1cm4(s) * pow(f_p1cm(s), -4); };
  /// Transfer cross-section as defined in 10.1103/PhysRevD.96.115010 multiplied by \f$p_1^4\f$ (in the center of mass)
  double EvalSigmaTransferBringmannP1cm4(double omega);
  /// Transfer cross-section as defined in 10.1103/PhysRevD.96.115010 multiplied by \f$p_1^4\f$ (in the center of mass) evaluated from numerical integration
  double EvalSigmaTransferBringmannP1cm4_fromIntegral(double omega);

  /// Standard cross section evaluated in term of the relative velocity
  double EvalSigma_fVrel(double vrel) { return EvalSigma(s_vs_vrel(vrel)); };
  /// Standard cross section evaluated in term of the relative velocity and from the numerical integration
  double EvalSigma_fromIntegral_fVrel(double vrel) { return EvalSigma_fromIntegral(s_vs_vrel(vrel)); };
  /// Transfer Cross-section evaluated in term of the relative velocity
  double EvalSigmaTransfer_fVrel(double vrel) { return EvalSigmaTransfer(s_vs_vrel(vrel)); };

protected:
  // ================================================
  // Protected functions (for daughter class)
  // ================================================

  /// Initialise the particle content of the interaction and give it a name
  void Initialise(std::vector<Particle> &, std::vector<Particle> &, std::vector<Particle> &, std::string name);

  /** Initialisation of the virtual functions for the coefficient polynomials
   * The variable t is here in order to account for depencies more complicated than just rational fractions */ 
  virtual dcomp Q(ParticlePair const &partPair, int n, double s, double t = 0) = 0;
  /// Initialisation of order of the coefficient polynomial
  virtual int Q_order(ParticlePair const &partPair) = 0;

  // Redefinition of function pow in Power for the interface with mathematica
  template <typename T1, typename T2>
  T1 Power(T1 const &n1, T2 const &n2) { return pow(n1, n2); };

  // Error message function for the access to the coefficient of the polynomial decomposition of the amplitude
  void printErrorMessage(std::string func_name)
  {
    std::cout << "FATAL ERROR : " << func_name << " -> trying to access non defined coeff." << std::endl;
    exit(0);
  };

private:
  // ====================================================
  // Private functions used to evaluate the integrals
  // ====================================================

  /// Integrals (functionals) for the cross section
  dcomp IntLPtPt(ParticlePair const &partPair, double s);
  dcomp IntLPtPu(ParticlePair const &partPair, double s);
  dcomp IntLPuPu(ParticlePair const &partPair, double s);
  dcomp IntLP(ParticlePair const &partPair, double s);
  dcomp IntL1(ParticlePair const &partPair, double s);

  // Integrals (functionals) for the s-wave term computation
  dcomp IntLHPtPt(ParticlePair const &partPair, double s);
  dcomp IntLHPtPu(ParticlePair const &partPair, double s);
  dcomp IntLHPuPu(ParticlePair const &partPair, double s);
  dcomp IntLHP(ParticlePair const &partPair, double s);
  dcomp IntLH1(ParticlePair const &partPair, double s);

  // Integrals (functionals) for the transfer cross-section
  dcomp IntRPtPt(ParticlePair const &partPair, double s);
  dcomp IntRPtPu(ParticlePair const &partPair, double s);
  dcomp IntRPuPu(ParticlePair const &partPair, double s);
  dcomp IntRP(ParticlePair const &partPair, double s);
  dcomp IntR1(ParticlePair const &partPair, double s);

  // Function to compute integral over (only use for the full numerical integration)
  double fToIntOnt_for_EvalSigmaP1cm2_fromIntegral(double t, double s);
  double fToIntOnt_for_EvalSigmaTransferBringmannP1cm4_fromIntegral(double t, double s);

  // The value of the Center Of Mass momenta
  // and bounds of Mandelstam t variable
  double f_p1cm(double s);
  double f_p3cm(double s);
  double tmin(double s);
  double tmax(double s);

public:

  /// Gives the s-wave term
  double EvalSWaveTerm();

  double EvalSWaveTermInterfChannelS();

  /// Mandelstam variable \f$s\f$ v.s. relative velocity
  double s_vs_vrel(double vrel);

  /// Getter for the particle \f$p_1\f$
  Particle const get_particle_input_1() const { return _p_1; };
  /// Getter for the particle \f$p_2\f$
  Particle const get_particle_input_2() const { return _p_2; };
  /// Getter for the particle \f$p_3\f$
  Particle const get_particle_input_3() const { return _p_3; };
  /// Getter for the particle \f$p_4\f$
  Particle const get_particle_input_4() const { return _p_4; };

  /**  Plot functions (overriden to take into account several configurations)  
  * If var_type = "vrel" plot is made w.r.t. relative velocity (linear scale)  
  * Othersiwe plot is made w.r.t. the Mandelsam variable s (logatmic scale) */
  void plot(double var_max, std::string name, std::string var_type = "s");
  void plot(std::string name, std::string var_type = "s")
  {
    if (var_type == "vrel")
      plot(1, name, var_type);
    else
      plot(100 * pow(_m1 + _m2, 2), name, var_type);
  };
  
  /**  Plot functions (overriden to take into account several configurations)  
  * If var_type = "vrel" plot is made w.r.t. relative velocity (linear scale)  
  * Othersiwe plot is made w.r.t. the Mandelsam variable s (logatmic scale) 
  * -> Here the numerical integration method is used and not the analytic decomposition */
  void plot_fromIntegral (double var_max, std::string name, std::string var_type = "s");
  void plot_fromIntegral(std::string name, std::string var_type = "s")
  {
    if (var_type == "vrel")
      plot_fromIntegral(1, name, var_type);
    else
      plot_fromIntegral(10 * pow(_m1 + _m2, 2), name, var_type);
  };

  void plot_sigmaTransfer(double s_max, std::string name, std::string var_type = "s");
  void plot_sigmaTransfer_fromIntegral(double s_max, std::string name, std::string var_type = "s");

  void plot_tminAndtmax(double s_max);
  void plot_L1(ParticlePair const &partPair, double s_max);
  void plot_RealQ(ParticlePair const &partPair, int order, double s_max);

  // ====================================================
  // Functions to test all the terms separately
  // ====================================================

  // Normal Cross-section
  double EvalSigmaInterfChannelSP1cm2(double s);
  double EvalSigmaInterfChannelS(double s) { return EvalSigmaInterfChannelSP1cm2(s) / (f_p1cm(s) * f_p1cm(s)); };
  double EvalSigmaInterfChannelTP1cm2(double s);
  double EvalSigmaInterfChannelT(double s) { return EvalSigmaInterfChannelTP1cm2(s) / (f_p1cm(s) * f_p1cm(s)); };
  double EvalSigmaInterfChannelUP1cm2(double s);
  double EvalSigmaInterfChannelU(double s) { return EvalSigmaInterfChannelUP1cm2(s) / (f_p1cm(s) * f_p1cm(s)); };
  double EvalSigmaInterfChannelSTP1cm2(double s, bool turnOnForSameParticle);
  double EvalSigmaInterfChannelST(double s, bool turnOnForSameParticle) { return EvalSigmaInterfChannelSTP1cm2(s, turnOnForSameParticle) / (f_p1cm(s) * f_p1cm(s)); };
  double EvalSigmaInterfChannelSUP1cm2(double s, bool turnOnForSameParticle);
  double EvalSigmaInterfChannelSU(double s, bool turnOnForSameParticle) { return EvalSigmaInterfChannelSUP1cm2(s, turnOnForSameParticle) / (f_p1cm(s) * f_p1cm(s)); };
  double EvalSigmaInterfChannelTUP1cm2(double s, bool turnOnForSameParticle);
  double EvalSigmaInterfChannelTU(double s, bool turnOnForSameParticle) { return EvalSigmaInterfChannelTUP1cm2(s, turnOnForSameParticle) / (f_p1cm(s) * f_p1cm(s)); };
  double EvalSigmaSquareChannelSP1cm2(double s);
  double EvalSigmaSquareChannelS(double s) { return EvalSigmaSquareChannelSP1cm2(s) / (f_p1cm(s) * f_p1cm(s)); };
  double EvalSigmaSquareChannelTP1cm2(double s);
  double EvalSigmaSquareChannelT(double s) { return EvalSigmaSquareChannelTP1cm2(s) / (f_p1cm(s) * f_p1cm(s)); };
  double EvalSigmaSquareChannelUP1cm2(double s);
  double EvalSigmaSquareChannelU(double s) { return EvalSigmaSquareChannelUP1cm2(s) / (f_p1cm(s) * f_p1cm(s)); };

  double EvalSigmaTransferSquareChannelSP1cm4(double s);
  double EvalSigmaTransferSquareChannelS(double s) { return EvalSigmaTransferSquareChannelSP1cm4(s) * pow(f_p1cm(s), -4); };
  double EvalSigmaTransferSquareChannelS_vrel(double vrel)
  {
    double s = pow(_m1 - _m2, 2) + 2 * _m1 * _m2 * (1. + (1. / sqrt(1 - vrel * vrel)));
    return EvalSigmaTransferSquareChannelSP1cm4(s) * pow(f_p1cm(s), -4);
  };
  double EvalSigmaTransferSquareChannelTP1cm4(double s);
  double EvalSigmaTransferSquareChannelT(double s) { return EvalSigmaTransferSquareChannelTP1cm4(s) * pow(f_p1cm(s), -4); };
  double EvalSigmaTransferSquareChannelT_vrel(double vrel)
  {
    double s = pow(_m1 - _m2, 2) + 2 * _m1 * _m2 * (1. + (1. / sqrt(1 - vrel * vrel)));
    return EvalSigmaTransferSquareChannelTP1cm4(s) * pow(f_p1cm(s), -4);
  };
  double EvalSigmaTransferSquareChannelUP1cm4(double s);
  double EvalSigmaTransferSquareChannelU(double s) { return EvalSigmaTransferSquareChannelUP1cm4(s) * pow(f_p1cm(s), -4); };
  double EvalSigmaTransferSquareChannelU_vrel(double vrel)
  {
    double s = pow(_m1 - _m2, 2) + 2 * _m1 * _m2 * (1. + (1. / sqrt(1 - vrel * vrel)));
    return EvalSigmaTransferSquareChannelUP1cm4(s) * pow(f_p1cm(s), -4);
  };

  double EvalSigmaTransferInterfChannelSTP1cm4(double s);
  double EvalSigmaTransferInterfChannelST(double s) { return EvalSigmaTransferInterfChannelSTP1cm4(s) * pow(f_p1cm(s), -4); };
  double EvalSigmaTransferInterfChannelST_vrel(double vrel)
  {
    double s = pow(_m1 - _m2, 2) + 2 * _m1 * _m2 * (1. + (1. / sqrt(1 - vrel * vrel)));
    return EvalSigmaTransferInterfChannelSTP1cm4(s) * pow(f_p1cm(s), -4);
  };
  double EvalSigmaTransferInterfChannelSUP1cm4(double s);
  double EvalSigmaTransferInterfChannelSU(double s) { return EvalSigmaTransferInterfChannelSUP1cm4(s) * pow(f_p1cm(s), -4); };
  double EvalSigmaTransferInterfChannelSU_vrel(double vrel)
  {
    double s = pow(_m1 - _m2, 2) + 2 * _m1 * _m2 * (1. + (1. / sqrt(1 - vrel * vrel)));
    return EvalSigmaTransferInterfChannelSUP1cm4(s) * pow(f_p1cm(s), -4);
  };
  double EvalSigmaTransferInterfChannelTUP1cm4(double s);
  double EvalSigmaTransferInterfChannelTU(double s) { return EvalSigmaTransferInterfChannelTUP1cm4(s) * pow(f_p1cm(s), -4); };
  double EvalSigmaTransferInterfChannelTU_vrel(double vrel)
  {
    double s = pow(_m1 - _m2, 2) + 2 * _m1 * _m2 * (1. + (1. / sqrt(1 - vrel * vrel)));
    return EvalSigmaTransferInterfChannelTUP1cm4(s) * pow(f_p1cm(s), -4);
  };

  static double CallBack_fToIntOnt_for_EvalSigmaP1cm2_fromIntegral(void *pt2Object, std::vector<double> xx) { return ((CrossSection *)pt2Object)->fToIntOnt_for_EvalSigmaP1cm2_fromIntegral(xx[0], xx[1]); };
  static double CallBack_fToIntOnt_for_EvalSigmaTransferBringmannP1cm4_fromIntegral(void *pt2Object, std::vector<double> xx) { return ((CrossSection *)pt2Object)->fToIntOnt_for_EvalSigmaTransferBringmannP1cm4_fromIntegral(xx[0], xx[1]); };

  /*



  */

private:
  SpecialIntegrations speInt;

  Particle _p_1, _p_2, _p_3, _p_4;

  std::vector<ParticlePair> ParticlePairsSSSq;
  std::vector<ParticlePair> ParticlePairsTTSq;
  std::vector<ParticlePair> ParticlePairsUUSq;
  std::vector<ParticlePair> ParticlePairsSSInt;
  std::vector<ParticlePair> ParticlePairsTTInt;
  std::vector<ParticlePair> ParticlePairsUUInt;
  std::vector<ParticlePair> ParticlePairsSTInt;
  std::vector<ParticlePair> ParticlePairsSUInt;
  std::vector<ParticlePair> ParticlePairsTUInt;

protected:
  // std::vector<ff> polyQ;
  // std::vector<unsigned int> polyQ_order;

  double _m1, _m2, _m3, _m4;
  int _g1, _g2;
  double _sym;

  int ns, nt, nu;

  double mandSum;
};

#endif
