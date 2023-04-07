#include "../headers/ChemicalDecoupling.h"

// ################### Constructor #####################

ChemicalDecoupling::ChemicalDecoupling(DarkSectorModel const &_DS, DegreeFreedom *_degFree, CS_TYPE n_cs_type, double n_simpleSigmaV)
{

  chi = _DS.get_chi();
  int n_chi = chi.size();

  if (n_chi > 0)
    is_CoAnnihilation = true;
  else
    is_CoAnnihilation = false;

  // Initialisation of the simple values
  simpleSigmaV = n_simpleSigmaV;
  // Variable to know if it is necessary to save the Y resolution results
  saveYToPrint = false;

  cs_type = n_cs_type;

  // Table containing the value of the s-wave [0] and p-wave [1] terms
  coeffSigmaV.resize(2);
  coeffSigmaV[0] = 0;
  coeffSigmaV[1] = 0;

  if (cs_type == CS_TYPE::S_WAVE_ONLY)
    coeffSigmaV[0] = simpleSigmaV;
  if (cs_type == CS_TYPE::P_WAVE_ONLY)
    coeffSigmaV[1] = simpleSigmaV;

  m_chi.resize(n_chi);
  g_chi.resize(n_chi);

  crossSection.resize(n_chi);
  for (int i = 0; i < n_chi; i++)
  {
    m_chi[i] = chi[i].get_mass();
    g_chi[i] = chi[i].get_degFree();
    crossSection[i].resize(n_chi);
  }

  std::vector<Particle> phi_s = _DS.get_phis();
  std::vector<Particle> phi_p = _DS.get_phip();
  int n_phis = phi_s.size();
  int n_phip = phi_p.size();

  // Masses and  widths of the propagators to adapt the integration scheme
  m_prop.resize(n_phis + n_phip);
  w_prop.resize(n_phis + n_phip);



  for (int i = 0; i < n_phis; i++)
  {
    m_prop[i] = phi_s[i].get_mass();
    w_prop[i] = phi_s[i].get_width();
  }

  for (int i = n_phis; i < n_phis + n_phip; i++)
  {
    m_prop[i] = phi_p[i].get_mass();
    w_prop[i] = phi_p[i].get_width();
  }

  // Make a loop on all fermions that count in the model
  n_crossSections_ccff = StandardModel::getInstance()->get_n_couples_SM_ferm();

  /** We sum over all possible cross-sections
    * Note that crossSection[i][j] = crossSection[j][i] therefore we use this symmetry
    * For Dirac DM we also have crossSection[i][i] = 0 (no terms on the diagonal)
    * Therefore we also use this symmetry to reduce computation time */

  int jmax = 0;
  for (int i = 0; i < n_chi; i++)
  {

    jmax = i;
    if (chi[0].get_fermiontype() == Fermiontype::dirac)
      jmax = i - 1;

    for (int j = 0; j <= jmax; j++)
    {
      for (int k = 0; k < n_crossSections_ccff; k++)
        crossSection[i][j].push_back(new CrossSection_ccff(_DS, chi[i], chi[j], StandardModel::getInstance()->get_couples_SM_ferm(k)));

      for (int k = 0; k < n_phis; k++)
        for (int l = 0; l < n_phis; l++)
          crossSection[i][j].push_back(new CrossSection_ccss(_DS, chi[i], chi[j], phi_s[k], phi_s[l]));

      for (int k = 0; k < n_phip; k++)
        for (int l = 0; l < n_phip; l++)
          crossSection[i][j].push_back(new CrossSection_ccpp(_DS, chi[i], chi[j], phi_p[k], phi_p[l]));

      for (int k = 0; k < n_phis; k++)
        for (int l = 0; l < n_phip; l++)
          crossSection[i][j].push_back(new CrossSection_ccsp(_DS, chi[i], chi[j], phi_s[k], phi_p[l]));
    }
  }

  /*
  for(int i = 0; i < crossSection.size(); i++)
    for(int j = 0; j < crossSection[i].size(); j++)
      for(int k = 0; k < crossSection[i][j].size(); k++)
      {
        std::cout << i << " " << j << " " << k << " " << crossSection[i][j][k]->get_particle_input_3().get_name() << std::endl;
      }
  */
  //exit(0);

  // Get the spline from degFree
  spline_gSmooth = _degFree->get_gSmoothInterp();
  spline_hSmooth = _degFree->get_hSmoothInterp();
  spline_derLogHSmooth = _degFree->get_derLogHSmoothInterp();
  spline_gStarHalf = _degFree->get_gStarHalfInterp();

  // Get the temperature of QCD phase transition
  T_transQCD = _degFree->get_T_transQCD();

  /** Initialise the values of the vectors that can
   * be printed to follow the evolution of the resolution */
  xPrint.resize(0);
  YeqPrint.resize(0);
  YIEPrint.resize(0);
  YRK6Print.resize(0);

  // workspace for the integration with the gsl library
  n_points_int = 1000;
  workspace = gsl_integration_workspace_alloc(n_points_int);

  if (cs_type != CS_TYPE::P_WAVE_ONLY && cs_type != CS_TYPE::S_WAVE_ONLY)
  {
 
    // Computation of the s-wave term
    Initialise_s_wave_term();

    // In order to make the computation faster we integrate sigv
    InterpolateSigmaV();
  }
}

void ChemicalDecoupling::Initialise_s_wave_term()
{

  double s_wave, s_wave_ij = 0, s_wave_ijk = 0;

  // Search for all chi particles with the minimal mass
  // We use the fact that the DM particles are ordered
  // from the lowest mass to the largest.
  int imax = 0, jmax = 0;
  for (int i = 0; i < chi.size(); i++)
  {
    if (m_chi[i] == m_chi[0])
      imax = i;
  }

  // Regulator : Theta(log10(T)-log10(T_QCD_GEV)) smoothed
  //double reg = regulator(log10(m_chi[0]), 1.0 / 0.3, log10(T_QCD_GeV), 1, 0);

  //std::cout << "Number of DM particles of minimal mass : " << imax << std::endl;

  for (int i = 0; i <= imax; i++)
  {
    jmax = i;
    if (chi[0].get_fermiontype() == Fermiontype::dirac)
      jmax = i - 1;

    for (int j = 0; j <= jmax; j++)
    {
      for (int k = 0; k < crossSection[i][j].size(); k++)
      {
        s_wave_ijk = crossSection[i][j][k]->EvalSWaveTerm();

        //std::cout << i << " " << j << " " << k << " " << s_wave_i

        // Checking if the s-wave term is not nan
        if (s_wave_ij != s_wave_ij)
          std::cout << " WARNING : s-wave computation failed for cross-section : " << k << " in " << __PRETTY_FUNCTION__ << std::endl;

        // We do not consider quarks as at low temperature their are no longer in the plasma
        if (crossSection[i][j][k]->get_particle_input_3().get_QCDType() != QCDTYPE::QUARK)
        {

          s_wave_ij += s_wave_ijk;
        }
        //else
        // We pultiply by a factor of 3 to take into account the 3 colors per quarks
        // Indeed, this is not taken into account in the computation of the cross-section
        // s_wave_ij += 3 * reg * s_wave_ijk;
      }

      if (i != j)
        s_wave_ij = 2 * s_wave_ij;
      // Multiply by two in order to take into account s_wave_ij = s_wave_ji

      s_wave += s_wave_ij;
    }
  }

  coeffSigmaV[0] = s_wave;
}

/*
void ChemicalDecoupling::InterpolateCrossSectionP1cm2()
{

  spline_crossSection.resize(n_crossSections);

  int numPoints = 500000;

  std::vector<double> logs;
  std::vector<double> crossSec;
  logs.resize(numPoints);
  crossSec.resize(numPoints, 0);

  double s = 0;

  // Be carefull with these values
  double s_min = 4 * pow(chi[0].get_mass(), 2);
  double s_max = pow(chi[0].get_mass() * (2 + 30 * 1e+4), 2);

  double dellogs = (log10(s_max) - log10(s_min)) / numPoints;

  for (int j = 0; j < n_crossSections; j++)
  {
    for (int i = 0; i < numPoints; i++)
    {
      s = s_min * pow(10, i * dellogs);
      logs[i] = log10(s_min) + i * dellogs;

      crossSec[i] = crossSection[j]->EvalSigmaP1cm2(s);
      //std::cout << i << std::endl;
    }

    spline_crossSection[j].set_points(logs, crossSec, true);
  }
}*/

void ChemicalDecoupling::InterpolateSigmaV()
{

  /** Note that the reference "noquarks" means that quarks are removed from the
  * particles that can create DM particles and be created by DM particles after
  * the QCD phase transition */

  int numPoints = 250;

  std::vector<double> logT;                //, logT_noquarks;
  std::vector<double> sigmaV_vec_noquarks; // sigmaV_vec_all
  logT.resize(numPoints+1);
  //logT_noquarks.resize(0);
  //sigmaV_vec_all.resize(numPoints);
  sigmaV_vec_noquarks.resize(numPoints+1);

  double T = 0;

  // Be carefull with these values
  double xf = 0, T_max = 0, T_min = 0;

  // If we are in a freeze-out configuration we interpolate sigmaV
  try
  {
    xf = DeterminationXfMethod2();
    T_max = chi[0].get_mass();
    T_min = chi[0].get_mass() / 299.;
  }
  catch (const ExceptionHandler &e)
  {
    // Probably to be checked or re-written
    if (e.get_exception_type() == ExceptionType::FreezeIn)
    {
      T_max = chi[0].get_mass() / 0.001;
      T_min = chi[0].get_mass() / 299.;
    }

    if (e.get_exception_type() == ExceptionType::GSL_error)
    {
      ExceptionHandler ee((char *)"WARNING :: GSL exception", __LINE__, (char *)__PRETTY_FUNCTION__, ExceptionType::GSL_error);
      throw ee;
    }
  }

  //std::cout << SigmaV(T_min) << " aaa" << std::endl;

  double dellogT = (log10(T_max / T_min)) / numPoints;

  double xmin = m_chi[0] / T_max;
  double xmax = m_chi[0] / T_min;
  double x = 0, dx = log10(xmax / xmin) / numPoints;
  double reg = 0;

  try
  {
    for (int i = 0; i < numPoints + 1; i++)
    {
      //std::cout << "here " << T_min << std::endl;
      x = xmax * pow(10, -i * dx);
      logT[i] = log10(m_chi[0] / x);
      T = m_chi[0] / x;
      //sigmaV_vec_all[i] = SigmaV(T, true);
      //logT_noquarks.push_back(logT[i]);

      // We remove the annihilation/creation into/from quarks when the temperature is too low
      // We do that with a regulator following the phase transition as in DegreeOfFreedom class
      sigmaV_vec_noquarks[i] = SigmaV(T, false);
    }
  }
  catch (ExceptionHandler const &e)
  {
    if (e.get_exception_type() == ExceptionType::GSL_error)
    {
      ExceptionHandler ee((char *)"WARNING :: GSL exception", __LINE__, (char *)__PRETTY_FUNCTION__, ExceptionType::GSL_error);
      throw ee;

      /*for (int i = 0; i < sigmaV_vec_all.size(); i++)
      {
        x = xmax * pow(10, -i * dx);
        logT[i] = log10(m_chi[0] / x);
        sigmaV_vec_all[i] = 0;
      }*/
      for (int i = 0; i < sigmaV_vec_noquarks.size(); i++)
      {
        x = xmax * pow(10, -i * dx);
        logT[i] = log10(m_chi[0] / x);
        sigmaV_vec_noquarks[i] = 0;
      }
    }
    else
    {
      std::cout << e.what() << std::endl;
      exit(0);
    }
  }

  //std::cout << "here" << std::endl;

  //std::cout << "Interpolation of <sigma v_Mol> completed" << std::endl;

  //spline_sigmaV_all.set_points(logT, sigmaV_vec_all, false);
  spline_sigmaV_noquarks.set_points(logT, sigmaV_vec_noquarks, false);
}


/*






*/
// #############################################################################
// ###################### Evalutation of xf for model 1 #######################
// ############################################################################

double ChemicalDecoupling::DeterminationXfMethod1()
{

  double xMin = 1;
  double xMax = 200;

  int compt = 0, comptMax = 10;

  double result = Dichotomie_Static(this, CallBackfToSolveXfModel1Method1, xMin, xMax, 1e-10, 0);

  while (result == 0 && compt < comptMax)
  {
    xMin += 2.;
    result = Dichotomie_Static(this, CallBackfToSolveXfModel1Method1, xMin, xMax, 1e-10, 0);
    compt++;
  }

  std::cout << "Results method 1 : xf = " << result << std::endl;

  return result;
}

double ChemicalDecoupling::DeterminationXfMethod2()
{
  double xMin = 1.;
  double xMax = 200.;

  int compt = 0, comptMax = 10;
  std::vector<double> xx{1.};

  double result = 0;

  try
  {
    result = Dichotomie_Static(0, 1, xx, this, CallBackfToSolveXfModel1Method2, xMin, xMax, 1e-5, -1);
  }
  catch (ExceptionHandler const &e)
  {
    //std::cout << "Here !!" << std::endl;
    if (e.get_exception_type() == ExceptionType::GSL_error)
    {
      ExceptionHandler ee((char *)"WARNING :: FreezeIn exception", __LINE__, (char *)__PRETTY_FUNCTION__, ExceptionType::GSL_error);
      throw ee;
    }
  }

  xx[0] = sqrt(xMax * xMin);

  if (result == -1 && fToSolveXfModel1Method2(xx[0]) < 0)
  {
    //std::cout << "Bah ici" << std::endl;
    //std::cout << fToSolveXfModel1Method2(xMax) << " " << fToSolveXfModel1Method2(xMin) << std::endl;
    // Means that we are likely to be in a freeze-in scenario
    ExceptionHandler e((char *)"WARNING :: FreezeIn exception", __LINE__, (char *)__PRETTY_FUNCTION__, ExceptionType::FreezeIn);
    throw e;
  }

  return result;
}

double ChemicalDecoupling::fToSolveXfModel1Method1(double const &x)
{
  double cteF = g_chi[0] * sqrt(45 / 2) * ((M_PLANCK * m_chi[0]) / (4 * PI * PI * PI));
  double sqrtGeff = sqrt(spline_gSmooth(log10(m_chi[0] / x)));
  double sigmv = SigmaV(m_chi[0] / x);

  // Return the efficient cross section development
  return (cteF * sigmv * pow(x, 0.5) * exp(-x) / sqrtGeff - 1.0);
}

double ChemicalDecoupling::fToSolveXfModel1Method2(double const &x)
{
  double delta = 1.5; // According to Gondolo and Gelmini

  //std::cout << m_chi[0] << std::endl;

  double gStarHalf = spline_gStarHalf(log10(m_chi[0] / x));
  double heff = spline_hSmooth(log10(m_chi[0] / x));
  double sigmv = 0;
  double derH = spline_derLogHSmooth(log10(m_chi[0] / x));

  sigmv = SigmaV(m_chi[0] / x);

  double cteF = sqrt(PI / 45) * 45 * g_chi[0] * m_chi[0] * M_PLANCK * delta * (delta + 2) / (4 * pow(PI, 4));

  return cteF * BesselK2(x) * sigmv * gStarHalf / heff - BesselK1(x) / BesselK2(x) + derH / x;
}

void ChemicalDecoupling::OutputGammaAndH()
{
  std::ofstream out_GammaAndH;
  out_GammaAndH.open("../output/GammaAndH.out");
  out_GammaAndH.precision(8);

  std::ofstream out_GammaComplete;
  out_GammaComplete.open("../output/GammaComplete.out");
  out_GammaComplete.precision(8);

  double H;
  double Gamma1, Gamma2, Gamma3, GammaTot, GammaReal;
  double Geff;
  double sigma1 = 3e-26 / GeVm2_to_cm3persec;
  double sigma2 = 3e-26 / GeVm2_to_cm3persec;
  double sigma3 = 3e-26 / GeVm2_to_cm3persec;
  double sigmaReal;
  double x;

  std::vector<double> var;
  var.resize(2);

  int Npts = 1000;
  double xMin = 1;
  double xMax = 200;
  double dx = (xMax - xMin) / Npts;

  for (int i = 0; i < Npts + 1; i++)
  {
    x = xMin + dx * i;

    // test << x << "\t" << fToSolveXfModel1Method1(x) << std::endl;
    // test2 << x << " " <<  DegLibComp->ComputeGSmoothOnly(m_chi[0]/x) << std::endl;
    Geff = spline_gSmooth(log10(m_chi[0] / x));
    H = 2 * pow(PI, 1.5) * m_chi[0] * m_chi[0] * sqrt(Geff) / (sqrt(45) * M_PLANCK * x * x);

    Gamma1 = sigma1 * g_chi[0] * pow(m_chi[0], 3) * BesselK2(x) / (2 * PI * PI * x);
    Gamma2 = sigma2 / x * g_chi[0] * pow(m_chi[0], 3) * BesselK2(x) / (2 * PI * PI * x);
    Gamma3 = sigma3 / (x * x) * g_chi[0] * pow(m_chi[0], 3) * BesselK2(x) / (2 * PI * PI * x);

    sigmaReal = coeffSigmaV[0] + coeffSigmaV[1] / x;
    GammaReal = sigmaReal * g_chi[0] * pow(m_chi[0], 3) * BesselKn(2, x) / (2 * PI * PI * x);
    GammaTot = Gamma1 + Gamma2 + Gamma3;

    out_GammaAndH << x << " " << H << " " << GammaReal << " " << Gamma1 << " " << Gamma2 << " " << Gamma3 << " " << GammaTot << std::endl;

    out_GammaComplete << x << "\t" << H << "\t" << SigmaV(m_chi[0] / x) * g_chi[0] * pow(m_chi[0], 3) * BesselKn(2, x) / (2 * PI * PI * x) << std::endl;
  }

  out_GammaAndH.close();
  out_GammaComplete.close();
}

double ChemicalDecoupling::FunctionTest(std::vector<double> x)
{
  return exp(-(x[0] - 10) * (x[0] - 10));
}

/*






*/

double ChemicalDecoupling::fToIntegrateForSigmaV_noquarks(std::vector<double> variables)
{
  // WARNING :: WE COMPUTE SIGMA*(s-4m^2) IN ORDER TO PREVENT INFINIES
  // IN THE COMPUTATION LATER ON AT s=4m^2.
  double sigmaSm4M2 = 0, sigmaSm4M2_ij = 0;

  double z = variables[0]; // z = sqrt(s)*x_chi/m_chi
  double x = variables[1];
  double s = z * z * m_chi[0] * m_chi[0] / (x * x);
  double T = m_chi[0] / x;

  // Regulator : Theta(log10(T)-log10(T_QCD_GEV)) smoothed
  double reg = regulatorTQCD(log10(T), 1.0 / 0.3, log10(T_QCD_GeV), 1, 0);

  int jmax, n_chi = chi.size();
  double pij2, p002 = 1. / 4. * (s - 4 * pow(m_chi[0], 2));


  for (int i = 0; i < n_chi; i++)
  {
    jmax = i;
    if (chi[0].get_fermiontype() == Fermiontype::dirac)
      jmax = i - 1;

    for (int j = 0; j <= jmax; j++)
    {
      // Initialisation of the sum for a given initial particle pair
      sigmaSm4M2_ij = 0;
      for (int k = 0; k < crossSection[i][j].size(); k++)
      {
        if (crossSection[i][j][k]->get_particle_input_3().get_QCDType() != QCDTYPE::QUARK)
        {

          //std::cout << i << " " << j << " " << k << " " << 4*crossSection[i][j][k]->EvalSigmaP1cm2(s) << std::endl;
          //if(crossSection[i][j][k]->get_particle_input_3().get_name() == "electron")
          sigmaSm4M2_ij += 4 * crossSection[i][j][k]->EvalSigmaP1cm2(s);
          // Moreover we multply by 3 in order to take into account the 3 colors
        }
        else if (crossSection[i][j][k]->get_particle_input_3().get_QCDType() == QCDTYPE::QUARK)
        {
          // When dealing with quarks in the out states we only consider them if the temperature
          // is above the temperature of the QCD phase transition
          // Regulator : Theta(log10(T)-log10(T_QCD_GEV)) smoothed
          //std::cout << i << " " << j << " " << k << " " << reg * 3 * 4 * crossSection[i][j][k]->EvalSigmaP1cm2(s) << std::endl;
          sigmaSm4M2_ij += reg * 3 * 4 * crossSection[i][j][k]->EvalSigmaP1cm2(s);
          //*< The factor of 4 fromes from the fact that \f$p_00^2 = (s-4m_\chi[0]^2)/4\f$
        }
      }

      if (i != j)
        sigmaSm4M2_ij *= 2;
      // We multiply by two to take into account the fact that we also need to add sigma_ji = sigma_ij

      pij2 = 1. / (4. * s) * (s - pow(m_chi[i] + m_chi[j], 2)) * (s - pow(m_chi[i] - m_chi[j], 2));
      sigmaSm4M2_ij = g_chi[i] * g_chi[j] / pow(g_chi[0], 2) * pij2 / p002 * sigmaSm4M2_ij;

      sigmaSm4M2 = sigmaSm4M2 + sigmaSm4M2_ij;
    }
  }

  //std::cout  << " ------ " << z << " " << sigmaSm4M2 << std::endl;
  //std::cout << sigmaSm4M2 << std::endl;

  return  sigmaSm4M2 * (2 * s * m_chi[0] / x) * BesselK1(z);
}

double ChemicalDecoupling::fToIntegrateForSigmaV_all(std::vector<double> variables)
{
  // WARNING :: WE COMPUTE SIGMA*(s-4m^2) IN ORDER TO PREVENT INFINIES
  // IN THE COMPUTATION LATER ON AT s=4m^2.
  double sigmaSm4M2 = 0, sigmaSm4M2_ij = 0;

  double z = variables[0]; // z = sqrt(s)*x_chi/m_chi
  double x = variables[1];
  double s = z * z * m_chi[0] * m_chi[0] / (x * x);
  double T = m_chi[0] / x;

  int jmax, n_chi = chi.size();
  double pij2, p002 = 1. / 4. * (s - 4 * pow(m_chi[0], 2));
  for (int i = 0; i < n_chi; i++)
  {
    jmax = i;
    if (chi[0].get_fermiontype() == Fermiontype::dirac)
      jmax = i - 1;

    for (int j = 0; j <= jmax; j++)
    {
      // Initialisation of the sum for a given initial particle pair
      sigmaSm4M2_ij = 0;
      for (int k = 0; k < crossSection[i][j].size(); k++)
      {
        if (crossSection[i][j][k]->get_particle_input_3().get_QCDType() == QCDTYPE::QUARK)
        {
          // We multply by 3 in order to take into account the 3 colors
          sigmaSm4M2_ij += 3 * 4 * crossSection[i][j][k]->EvalSigmaP1cm2(s);
          //*< The factor of 4 fromes from the fact that \f$p_00^2 = (s-4m_\chi[0]^2)/4\f$
        }
        else
          sigmaSm4M2_ij += 4 * crossSection[i][j][k]->EvalSigmaP1cm2(s);
      }

      if (i != j)
        sigmaSm4M2_ij *= 2;
      // We multiply by two to take into account the fact that we also need to add sigma_ji = sigma_ij

      pij2 = 1. / (4. * s) * (s - pow(m_chi[i] + m_chi[j], 2)) * (s - pow(m_chi[i] - m_chi[j], 2));
      sigmaSm4M2_ij = g_chi[i] * g_chi[j] / pow(g_chi[0], 2) * pij2 / p002 * sigmaSm4M2_ij;

      sigmaSm4M2 = sigmaSm4M2 + sigmaSm4M2_ij;
    }
  }

  return sigmaSm4M2 * (2 * s * m_chi[0] / x) * BesselK1(z);
}

/*
double ChemicalDecoupling::FunctionToIntegrateForSigmaV(double z, void* par)
{ 
  gsl_p = (gsl_f_params *)par;
  return FunctionToIntegrateForSigmaV(par->parameters);
};*/

void ChemicalDecoupling:: plotFuncToIntegrateSigmaV(double T)
{

  double z = 0; // z= sqrt(s)*x/m_chi[0]
  double x = (1. * m_chi[0]) / (1. * T);

  std::vector<double> variables;
  variables.resize(2);
  variables[0] = z;
  variables[1] = x;

  double zMin = 2 * x;
  double zMax = zMin + 30; // To be precise in the computation at 1e-4

  int nPoints = 10000;
  double delta = log10(zMax / zMin) / (1.0 * nPoints);

  std::ofstream testEncore;
  testEncore.open("../output/FToIntSigmaVAverage.out");
  testEncore.precision(8);

  for (int i = 0; i < nPoints + 1; i++)
  {
    variables[0] = zMin * pow(10, delta * i);
    testEncore << variables[0] << " " << fToIntegrateForSigmaV_noquarks(variables) << " " << BesselK1(variables[0]) << std::endl;
  }

  testEncore.close();
}

double ChemicalDecoupling::SigmaV(double T, bool with_quarks_after_TQCD)
{
  double z = 0; // z= sqrt(s)*x/m_chi[0]
  double x = (1. * m_chi[0]) / (1. * T);

  // Be Carefull SigmaV(m_chi[0]/x) i only possible for x < 300 approximately
  // Otherwise Bessel functions in the expression make divergences appear

  // Condition if one has to use a simplified version
  if (cs_type == CS_TYPE::SIMPLIFIED_ONLY)
  {
    double sigmaV = 0;
    std::vector<double> coeff;
    coeff.resize(2);
    for (int i = 0; i < n_crossSections_ccff; i++)
    {
      sigmaV += coeffSigmaV[0] + coeffSigmaV[1] / x;
    }
    return sigmaV;
  }
  if (cs_type == CS_TYPE::S_WAVE_ONLY)
    return simpleSigmaV;
  if (cs_type == CS_TYPE::P_WAVE_ONLY)
    return simpleSigmaV / x;
  if (cs_type == CS_TYPE::FULL_SIMPLIFIED && x > 300)
  {
    double sigmaV = 0;
    std::vector<double> coeff;
    coeff.resize(2);
    for (int i = 0; i < n_crossSections_ccff; i++)
    {
      sigmaV += coeffSigmaV[0] + coeffSigmaV[1] / x;
    }
    return sigmaV;
  }
  if (cs_type == CS_TYPE::FULL_ONLY && x > 300)
  {
    std::cout << "Not possible to evaluate numerically SigmaV() over x > 300 we have to use the approximation given by the expension if possible" << std::endl;
    std::cout << "This is due to Bessel function decreasing too fast. However, this expension for the complete model is available only if cs_type is FULL_SIMPLIFIED" << std::endl;
    std::cout << "Please modify this parameter if you want to access values of x > 300" << std::endl;

    exit(0);
  }

  // double sMin = 4 * m_chi[0] * m_chi[0];
  // double sMax = 3e+6;

  double zMin = 2 * x;
  double zMax = zMin + 30; // To be precise in the computation at 1e-4

  /*
  std::vector<double> z_mass_prop, z_interm_mass_prop;
  z_mass_prop.resize(m_prop.size());

  if (m_prop.size() > 0)
    z_interm_mass_prop.resize(m_prop.size() - 1);
  else
    z_interm_mass_prop.resize(0);

  for (int i = 0; i < m_prop.size(); i++)
    z_mass_prop[i] = m_prop[i] * x / m_chi[0];

  for (int i = 1; i < m_prop.size(); i++)
    z_interm_mass_prop[i - 1] = pow(10, log10(z_mass_prop[i] * z_mass_prop[i - 1]) / 2);
*/
  double integral = 0;

  std::vector<double> variables;
  variables.resize(2);
  variables[0] = z;
  variables[1] = x;

  //std::cout << zMin << " " << z_mass_prop[0] << " " << z_mass_prop[z_mass_prop.size()-1] << " " << zMax << std::endl;

  // Integrate the function over z
  //std::cout << z_mass_prop.size() << " "  << z_interm_mass_prop.size() << std::endl;
  /*if(zMin < z_mass_prop[0])
    integral += GaussLegendre_IntegralInverseLnLn_Static(0, 2, gl500, zMin, z_mass_prop[0], variables, this, CallBackFunctionToIntegrateForSigmaV);
  else
    intgral += GaussLegendre_IntegralLn_Static(0, 2, gl500, zMin, z_mass_prop[0], variables, this, CallBackFunctionToIntegrateForSigmaV);
  
  if(zMax > z_mass_prop[z_mass_prop.size()-1])
    integral += GaussLegendre_IntegralLnLn_Static(0, 2, gl500, z_mass_prop[z_mass_prop.size()-1], zMax, variables, this, CallBackFunctionToIntegrateForSigmaV);

  for (int i = 1; i < m_prop.size() - 1; i++)
  {
      integral += GaussLegendre_IntegralLnLn_Static(0, 2, gl500, z_mass_prop[i - 1], z_interm_mass_prop[i-1], variables, this, CallBackFunctionToIntegrateForSigmaV);
      integral += GaussLegendre_IntegralInverseLnLn_Static(0, 2, gl500, z_interm_mass_prop[i-1], z_mass_prop[i], variables, this, CallBackFunctionToIntegrateForSigmaV);
  }*/

  //std::cout << "Houra : " << T << std::endl;

  double err;

  gsl_f_params *gsl_p = new gsl_f_params();

  gsl_p->pt_MyClass = this;
  gsl_p->_x = x;

  gsl_function F;
  F.params = gsl_p;

  if (with_quarks_after_TQCD == true)
    F.function = &gslClassWrapperLn_fToIntegrateForSigmaV_all;
  else
    F.function = &gslClassWrapperLn_fToIntegrateForSigmaV_noquarks;

  int n_singularities = m_prop.size(); // Position of the pseudo-singularities

  std::vector<double> sing_vec;
  sing_vec.resize(0);
  sing_vec.push_back(log(zMin));

  for (int i = 1; i < n_singularities + 1; i++)
  {
    if (log(m_prop[i - 1] * x / m_chi[0]) > log(zMin) && log(m_prop[i - 1] * x / m_chi[0]) < log(zMax))
    {
      sing_vec.push_back(log(m_prop[i - 1] * x / m_chi[0]));
    }
  }

  sing_vec.push_back(log(zMax));

  if (sing_vec.size() == 2) // no singularity in the integration range
  {
    if (with_quarks_after_TQCD == true)
      integral = GaussLegendre_IntegralLn_Static(0, 2, gl500, zMin, zMax, variables, this, CallBackfToIntegrateForSigmaV_all);
    else
      integral = GaussLegendre_IntegralLn_Static(0, 2, gl500, zMin, zMax, variables, this, CallBackfToIntegrateForSigmaV_noquarks);
  }
  else
  {
    double eps = 0;
    for (int i = 0; i < sing_vec.size() - 1; i++)
    {
      double integ = 0;
      gsl_integration_qag(&F, sing_vec[i], sing_vec[i + 1], 0, 1e-7, n_points_int, 1, workspace, &integ, &err);
      integral += integ;
    }
  }

  //gsl_integration_qag(&F, log(zMin), log(zMax), 0, 1e-7, n_points_int, 2, workspace, &integral, &err);
  //gsl_integration_qags(&F, log(zMin), log(zMax), 0, 1e-7, 1000, workspace, &integral, &err);
  //gsl_integration_qagp(&F, sing, sing_vec.size(), 0, 1e-2, n_points_int, workspace, &integral, &err);

  // WARNING : Do not forget to delete the pointer or memory leaks can occur
  delete gsl_p;

  //integral = GaussLegendre_IntegralLn_Static(0, 2, gl500, zMin, zMax, variables, this, CallBackFunctionToIntegrateForSigmaV);
  double first_factor = 0;
  //std::cout << "here : " << chi.size() << std::endl;
  for (int i = 0; i < chi.size(); i++)
  {
    first_factor += g_chi[i] * pow(m_chi[i], 2) * BesselK2(m_chi[i] / T);
  }
  first_factor = 1 / pow(first_factor, 2);
  first_factor *= pow(g_chi[0], 2) / (8 * T);

  //std::cout << first_factor << std::endl;

  return first_factor * integral;
  //return (x * integral) / (8 * pow(m_chi[0], 5) * pow(BesselK2(x), 2));
}

void ChemicalDecoupling::plotSigmaV()
{
  double x = 0;

  double xMin = 1;
  double xMax = 299; // To be precise in the computation at 1e-4

  int nPoints = 100;
  double delta = log10(xMax / xMin) / (1.0 * nPoints);

  std::ofstream outfile;
  outfile.open("../output/SigmaV.out");

  outfile << "# x | sigv(x) | s-wave term | interpolated sigv(x)" << std::endl;

  //crossSection[0]->plot(1e+10, "sigmaV");

  for (int i = 0; i < nPoints + 1; i++)
  {
    x = xMin * pow(10, delta * i);
    outfile << x << " " << SigmaV(m_chi[0] / x, true)  << " " << m_chi[0]/x << std::endl; // << " " << SigmaV(m_chi[0] / x, false) << " " << coeffSigmaV[0] << " " << spline_sigmaV_all(log10(m_chi[0] / x)) << " " << spline_sigmaV_noquarks(log10(m_chi[0] / x)) << std::endl;
  }

  outfile.close();

  /* std::vector<double> coeff;
    coeff.resize(2);
  
    coeff = ((CrossSection_ccff *)crossSection[0])->AppSigmaVChanSPropSPS();
    std::cout << coeff[0]  << " " << coeff[1] << std::endl;*/
}

/*






*/
// ####################################################################################
// ###################### Abundance equation resolution ###############################
// ####################################################################################

void ChemicalDecoupling::DerivativeOfComovingDensityNonStiff(std::vector<double> variables, double W, double &deriv)
{
  double x = variables[0];
  double Yeq = variables[1];

  deriv = sqrt(PI / 45) * M_PLANCK * m_chi[0] * spline_gStarHalf(log10(m_chi[0] / x)) *
          SigmaV(m_chi[0] / x) * (Yeq * Yeq * exp(-W) - exp(W)) / (x * x);

  return;
}

std::vector<double> ChemicalDecoupling::SolveFreezeOut()
{
  std::vector<double> variables;
  variables.resize(2);

  double x;
  double Y, Yeq, W; // Where W will represents ln(Y) : W = ln(Y)


  double xf = DeterminationXfMethod2();

     
  if (xf == -2)
  {
    variables[0] = -2;
    variables[1] = 0;

    return variables;
  }
  else if (xf == -1)
  {
    std::cout << "WARNING : Impossible to solve an estimation of xf in " << __PRETTY_FUNCTION__ << std::endl;
    return variables;
  }

  // double xMin = xf - 15;
  double xMin = 1; // Freeze-in study
  double xMax = xf + 150;
  double dx;

  double xMinBF = 0.9; // xMin for the part Before the Freeze-out
  double nPointsBF = 100;
  double delta = abs(xMin - xMinBF) / (1.0 * nPointsBF);

  int compt = 0;
  int comptMax = 500000;

  double normY = 1, normN;

  // Output of this function used to compute the abundance today
  std::vector<double> xYInitForAbundanceToday_xcd;
  xYInitForAbundanceToday_xcd.resize(3);

  xPrint.resize(0);
  YRK6Print.resize(0);
  YIEPrint.resize(0);
  YeqPrint.resize(0);

  /* WARNING : The Runge-Kutta resolution of the differential equation only works well
  if we start close to x of the decoupling otherwise it numerically diverges
  Before the freeze-out Y follow Yeq, therefore one uses that Y=Yeq before to complete
  the curves of Y before the freeze out.*/

  if (saveYToPrint == true)
  {

    normY = 45 * g_chi[0] * BesselK2(1) / (4 * pow(PI, 4) * spline_hSmooth(log10(m_chi[0])));

    for (int i = 0; i < nPointsBF + 1; i++)
    {
      x = xMinBF + delta * i;
      Yeq = 45 * g_chi[0] * x * x * BesselK2(x) / (4 * pow(PI, 4) * spline_hSmooth(log10(m_chi[0] / x)));

      xPrint.push_back(x);
      YRK6Print.push_back(Yeq / normY);
      YIEPrint.push_back(Yeq / normY);
      YeqPrint.push_back(Yeq / normY);
    }
  }

  // Initlizing x, W and dx
  x = xMin;
  dx = 0.02;
  Yeq = 45 * g_chi[0] * x * x * BesselK2(x) / (4 * pow(PI, 4) * spline_hSmooth(log10(m_chi[0] / x)));
  W = log(Yeq);
  Y = Yeq;

  bool varstep = false;
  double dxmin = dx, prec = 1e-3;
  bool okprec;
  bool xcd_defined = false;

  // ---- Loop for the abundance computation ----
  while (x < xMax && compt < comptMax)
  {

    Yeq = 45 * g_chi[0] * x * x * BesselK2(x) / (4 * pow(PI, 4) * spline_hSmooth(log10(m_chi[0] / x)));

    variables[0] = x;
    variables[1] = Yeq;

    //std::cout << "Verbose is " << Verbose << std::endl;

    if (compt % 100 == 0 && Verbose == true)
    {
      std::cout << "Numerical integration comoving number density : " << std::endl;
      std::cout << x << " " << Y << " " << Yeq << std::endl;
    }

    // Runge-Kutta resolution to go one step forward (resolution on W)
    // W = myRungeKutta4_Static(1, W, variables, dx, this, CallBackDerivativeOfAbundanceNonStiff);
    // W = myRungeKutta6_Static(0, variables, W, this, CallBackDerivativeOfAbundanceNonStiff, varstep, dx, dxmin, prec, okprec);


    // Impicit Euler resolutio to go one step forward (resolution on Y)
    Y = ImplicitEulerSolverSquares_Static(0, variables, Y, this, CallBackFunctionForImplicitEulerA,
                                             CallBackFunctionForImplicitEulerB,
                                             CallBackFunctionForImplicitEulerZ, dx);

    //std::cout << Y << std::endl;

    // Saving the results
    if (saveYToPrint)
    {
      xPrint.push_back(x);
      YRK6Print.push_back(exp(W) / normY);
      YIEPrint.push_back(Y / normY);
      YeqPrint.push_back(Yeq / normY);
    }

    x += dx; // Incrementing for the next step
    compt++; // Incrementing the compter too

    // std::cout << dx << std::endl;

    // Define the value of x at the decoupling xcd
    if (abs(Yeq / Y) < 1e-1 && xcd_defined == false)
    {
      xYInitForAbundanceToday_xcd[2] = x;
      xcd_defined = true;
    }

    // write the output
    // BE CAREFULL THIS HAS BEEN CHANGED
    if ((x < xMax && x > xMax - 1)) // || (abs(Yeq / Y) < 1e-3 && coeffSigmaV[0] != 0))
    {
      xYInitForAbundanceToday_xcd[0] = x;
      xYInitForAbundanceToday_xcd[1] = Y;

      //std::cout << "ici " << coeffSigmaV[0] << " " << x << " " << xMax << std::endl;

      //if (x < xMax && x > xMax - 1 && coeffSigmaV[0] != 0)
      //  std::cout << "WARNING : Comoving number density solver did not converge !!!" << std::endl;

      return xYInitForAbundanceToday_xcd;
    }
  }

  if (compt == comptMax)
    std::cout << "Error  : ChemicalDecoupling::SolveAbundance -> Not enough step to reach the final point" << std::endl;

  std::vector<double> outputError{-1};

  return outputError;
}

double ChemicalDecoupling::ComputeAbundanceToday(std::vector<double> const &xYInitForAbundanceToday)
{

  double xinit = xYInitForAbundanceToday[0];
  double Yinit = xYInitForAbundanceToday[1];
  double xmax = xinit * 1e+4;

  // std::cout << "The resolution of abundance today start at x=" << xinit << " for Y=" << Yinit << std::endl;

  std::vector<double> variables;
  variables.resize(1);

  // Plot of the function to integrate
  double nPoints = 100000;
  double delta = abs(xmax - xinit) / (1.0 * nPoints);

  double integral, Y0;

  if (coeffSigmaV[0] != 0)
  {
    integral = GaussLegendre_Integral_Static(0, 1, gl200, xinit, xmax, variables, this, CallBackFunctionToIntForAbundanceToday);
    //integral = Simpson_Integral1_Static(0, 1, 0, nPoints, xinit, xmax, variables, this, CallBackFunctionToIntForAbundanceToday, err);
    Y0 = Yinit / (1 + Yinit * pow(PI / 45, 0.5) * m_chi[0] * M_PLANCK * integral);
  }
  else
    Y0 = Yinit;

  return Y0;
}

double ChemicalDecoupling::FunctionToIntForAbundanceToday(std::vector<double> variables)
{
  double x = variables[0];
  return (spline_gStarHalf(log10(m_chi[0] / x)) / (x * x)) * (coeffSigmaV[0] + coeffSigmaV[1] / x);
}

void ChemicalDecoupling::SolveFreezeIn()
{
  double Yfi, x;
  std::vector<double> variables;
  variables.resize(2);
  double xMin = 0.0001;
  double xMax = 100.;
  double xMinBF = 0.00005;
  int compt = 0, comptMax = 50000;
  double Yeq;
  double normY;
  double du = 0.001;
  double u; // u = ln(x)

  double nPointsBF = 1000;
  double delta = abs(xMin - xMinBF) / (1.0 * nPointsBF);

  xPrint.resize(0);
  YFIPrint.resize(0);
  YeqPrint.resize(0);

  if (saveYToPrint == true)
  {
    for (int i = 0; i < nPointsBF + 1; i++)
    {
      x = xMinBF + delta * i;
      Yeq = 45 * g_chi[0] * x * x * BesselK2(x) / (4 * pow(PI, 4) * spline_hSmooth(log10(m_chi[0] / x)));

      if (i == 0)
        normY = Yeq;

      xPrint.push_back(x);
      YFIPrint.push_back(0);
      YeqPrint.push_back(Yeq / normY);
    }
  }

  Yfi = 0;
  x = xMin;
  u = log(x);

  while (x < xMax && compt < comptMax)
  {
    Yeq = 45 * g_chi[0] * x * x * BesselK2(x) / (4 * pow(PI, 4) * spline_hSmooth(log10(m_chi[0] / x)));

    variables[0] = x;
    variables[1] = Yeq;

    if (compt % 200 == 0)
      std::cout << x << std::endl;

    // Impicit Euler resolutio to go one step forward (resolution on Yfi)
    Yfi = ImplicitEulerSolverSquares_Static(0, variables, Yfi, this, CallBackFunctionForImplicitEulerLnA, CallBackFunctionForImplicitEulerZ, du);

    // Saving the results
    if (saveYToPrint)
    {
      xPrint.push_back(x);
      YFIPrint.push_back(Yfi / normY);
      YeqPrint.push_back(Yeq / normY);
    }

    u += du; // Incrementing for the next step
    x = exp(u);

    compt++; // Incrementing the compter too
  }
}

//
//

// Zone of tested functions
/*
double ChemicalDecoupling::ComputeYFreezeInApprox(double x)
{

  Y = sqrt(PI / 45) * M_PLANCK * m_chi[0];
  integral = Simpson_Integral1_Static(0, 1, 0, 1000, 1/x, xmax, variables, this, CallBackFunctionToIntForAbundanceToday, err);
}

double ChemicalDecoupling::FunctionToIntFreezeInApprox(std::vector<double> var)
{
  double x = var[0];
  double Yeq = var[1];

  return spline_gStarHalf(log10(m_chi[0] / x)) * SigmaV(m_chi[0] / x) * pow(Yeq, 2) / (x * x);
}*/

void ChemicalDecoupling::PlotYFO()
{
  // Output file
  std::ofstream out_Y, out_Y_norm;

  std::string name_file = "../output/Y";
  std::ostringstream strc;
  strc << m_chi[0];
  name_file += strc.str();

  out_Y.open(name_file + ".out");
  out_Y_norm.open(name_file + "_norm.out");
  out_Y.precision(8);
  out_Y_norm.precision(8);

  double yeq1 = 45 * g_chi[0] * BesselK2(1) / (4 * pow(PI, 4) * spline_hSmooth(log10(m_chi[0])));

  for (int i = 0; i < xPrint.size(); i++)
  {
    out_Y << xPrint[i] << "\t" << YIEPrint[i] * yeq1 << "\t" << YeqPrint[i] * yeq1 << "\t" << YRK6Print[i] * yeq1 << std::endl;
    out_Y_norm << xPrint[i] << "\t" << YIEPrint[i] << "\t" << YeqPrint[i] << "\t" << YRK6Print[i] << std::endl;
  }

  out_Y.close();
  out_Y_norm.close();
}

void ChemicalDecoupling::PlotYFI()
{
  // Output file
  std::ofstream out_Y;

  std::string name_file = "..output/Y_massDM.out";
  std::ostringstream strc;
  strc << m_chi[0];
  name_file += strc.str();

  out_Y.open("../output/" + name_file);
  out_Y.precision(8);

  for (int i = 0; i < xPrint.size(); i++)
  {
    out_Y << xPrint[i] << "\t" << YFIPrint[i] << "\t" << YeqPrint[i] << std::endl;
  }

  out_Y.close();
}

double ChemicalDecoupling::FunctionForImplicitEulerA(std::vector<double> variables)
{
  double x = variables[0];

  
  if(cs_type == CS_TYPE::FULL_ONLY)
    return sqrt(PI / 45) * M_PLANCK * m_chi[0] * spline_gStarHalf(log10(m_chi[0] / x)) * spline_sigmaV_noquarks(log10(m_chi[0] / x)) / (x * x);
  else
    return sqrt(PI / 45) * M_PLANCK * m_chi[0] * spline_gStarHalf(log10(m_chi[0] / x)) * SigmaV(m_chi[0] / x, false) / (x * x);

  // -> Use this return command instead in order to consider that all particles are present in the bath
  // return sqrt(PI / 45) * M_PLANCK * m_chi[0] * spline_gStarHalf(log10(m_chi[0] / x)) * spline_sigmaV_all(log10(m_chi[0] / x)) / (x * x);
}

double ChemicalDecoupling::FunctionForImplicitEulerB(std::vector<double> variables)
{
  return FunctionForImplicitEulerA(variables);
  //double x = variables[0];
  //return sqrt(PI / 45) * M_PLANCK * m_chi[0] * spline_gStarHalf(log10(m_chi[0] / x)) * spline_sigmaV_all(log10(m_chi[0] / x)) / (x * x);
}

double ChemicalDecoupling::FunctionForImplicitEulerLnA(std::vector<double> variables)
{
  double x = variables[0];
  return sqrt(PI / 45) * M_PLANCK * m_chi[0] * spline_gStarHalf(log10(m_chi[0] / x)) * spline_sigmaV_all(log10(m_chi[0] / x)) / (x);
}

double ChemicalDecoupling::FunctionForImplicitEulerZ(std::vector<double> variables)
{
  // Simply returns Yeq contained in variables[2]
  return variables[1];
}
