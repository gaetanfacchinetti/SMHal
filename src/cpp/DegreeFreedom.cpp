#include "../headers/DegreeFreedom.h"

DegreeFreedom::DegreeFreedom(double nT_transQCD, DarkSectorModel *_DS)
{

  // Define a shorter name for Parameters singleton
  StandardModel *_SM = StandardModel::getInstance();
  T_transQCD = nT_transQCD;

  bool includeDS;

  if (_DS != NULL)
    includeDS = true;
  else
    includeDS = false;

  nBessel = 100;    // Number of terms in sums
  numPoints = 1000; // Number of discretisation points
  Tmax = 5000;      // Maximal temperature (en GeV)
  Tmin = 0.00001;   // Minimal temperature (en GeV)

  logTmax = log10(Tmax);
  logTmin = log10(Tmin);
  deltaLogT = (logTmax - logTmin) / (1.0 * numPoints);

  all_Part.resize(0);
  std::vector<Particle> _SM_Part = _SM->get_SM_Part();
  all_Part.insert(std::end(all_Part), std::begin(_SM_Part), std::end(_SM_Part));

  int n_SM_Part = _SM_Part.size();

  if (includeDS)
  {
    std::vector<Particle> _DS_Part = _DS->get_DS_Part();
    all_Part.insert(std::end(all_Part), std::begin(_DS_Part), std::end(_DS_Part));
  }

  numSpecies = all_Part.size();

  degFree.resize(numSpecies);
  mass.resize(numSpecies);
  stat.resize(numSpecies);
  isPresentBeforeQCDPT.resize(numSpecies);
  isPresentAfterQCDPT.resize(numSpecies);

  int index = 0;

  for (int i = 0; i < _SM->get_SM_Part().size(); i++)
  {
    mass[_SM->get_SM_Part()[i].get_index()] = _SM->get_SM_Part()[i].get_mass();
    degFree[_SM->get_SM_Part()[i].get_index()] = _SM->get_SM_Part()[i].get_degFree();
    stat[_SM->get_SM_Part()[i].get_index()] = _SM->get_SM_Part()[i].get_stat();
    isPresentBeforeQCDPT[_SM->get_SM_Part()[i].get_index()] = _SM->get_SM_Part()[i].is_PresentBeforeQCDPT();
    isPresentAfterQCDPT[_SM->get_SM_Part()[i].get_index()] = _SM->get_SM_Part()[i].is_PresentAfterQCDPT();
    index++;
  }

  if (includeDS)
  {
    for (int i = 0; i < _DS->get_DS_Part().size(); i++)
    {
      mass[_DS->get_DS_Part()[i].get_index()] = _DS->get_DS_Part()[i].get_mass();
      degFree[_DS->get_DS_Part()[i].get_index()] = _DS->get_DS_Part()[i].get_degFree();
      stat[_DS->get_DS_Part()[i].get_index()] = _DS->get_DS_Part()[i].get_stat();
      isPresentBeforeQCDPT[_DS->get_DS_Part()[i].get_index()] = _DS->get_DS_Part()[i].is_PresentBeforeQCDPT();
      isPresentAfterQCDPT[_DS->get_DS_Part()[i].get_index()] = _DS->get_DS_Part()[i].is_PresentAfterQCDPT();
      index++;
    }
  }

  g_vec.resize(numPoints + 1);
  h_vec.resize(numPoints + 1);
  gSmooth_vec.resize(numPoints + 1);
  hSmooth_vec.resize(numPoints + 1);
  logH_vec.resize(numPoints + 1);
  logHSmooth_vec.resize(numPoints + 1);
  logT_vec.resize(numPoints + 1);

  n_neutrinos = _SM->get_SM_neutrinos().size();
  SM_neutrinos.resize(n_neutrinos);
  SM_neutrinos = _SM->get_SM_neutrinos();

  neutrinos_index.resize(n_neutrinos);
  for (int i = 0; i < n_neutrinos; i++)
    neutrinos_index[i] = SM_neutrinos[i].get_index();

  if (NUM_THREADS == 1)
  {
    std::clog << std::endl;
    std::clog << "===============================================" << std::endl;
    std::clog << "EFFECTIVE D.O.F. INITIALIZATION " << std::endl;
    std::clog << std::endl;
    std::clog << " -- Physical parameters : Particle content ----------------------------------------" << std::endl;
    std::clog.precision(4);

    std::clog << all_Part;

    std::clog << " ---------------------------------------------------------------------------------------- " << std::endl;

    std::clog.precision(8);
  }

  // Computation to check results
  float sumInit = 0;
  for (int i = 0; i < n_SM_Part; i++)
    if (isPresentBeforeQCDPT[i])
    {
      if (stat[i] == +1)
        sumInit += degFree[i];
      else
        sumInit += (7.0 / 8.0) * degFree[i];
    }

  if (NUM_THREADS == 1)
  {
    std::clog << " ---------------------------------------------------------------------------------------- " << std::endl;
    std::clog << " -- QCD phase transition & test computation -----------------------------------" << std::endl;
    std::clog << "| QCD phase transition is set at : " << T_transQCD << " GeV                                    |" << std::endl;
    std::clog << "| The number of degree of freedom at high temperature from SM must be : " << sumInit << " |" << std::endl;
    std::clog << " ------------------------------------------------------------------------------ " << std::endl;
  }

  // ComputeGAndH();
  ComputeGAndHSmooth();
  ComputeGStarHalf();
  isInterpolateCalled = false;
  Interpolate();
  InitialiseLattice();

  if (NUM_THREADS == 1)
    std::cout << "--> Degrees of freedom initialised : see log file" << std::endl;
}

//
//

void DegreeFreedom::ComputeGAndH()
{
  double logTransQCD = log10(T_transQCD);
  double logT = 0, g = 0, h = 0, x;

  for (int ilogT = 0; ilogT < numPoints + 1; ilogT++)
  {
    // New temperature
    logT = logTmax - ilogT * deltaLogT;
    logT_vec[ilogT] = logT;

    // Initialisation of the variable g
    g = 0;
    h = 0;

    // Loop on all the species
    for (int esp = 0; esp < numSpecies; esp++)
    {
      // Relevant variable x
      x = mass[esp] / pow(10, logT);

      // Only the relevant species at a given temperature are taken into account
      if (((isPresentBeforeQCDPT[esp] == true) && logT > logTransQCD) || ((isPresentAfterQCDPT[esp] == true) && logT < logTransQCD))
      {
        // Formula for massives species different from neutrino
        if (mass[esp] > 0 && ((find(neutrinos_index.begin(), neutrinos_index.end(), esp) == neutrinos_index.end()) || logT > log10(M_e)))
        {
          for (int n = 1; n < nBessel; n++)
          {
            g += (15.0 / pow(PI, 4)) * degFree[esp] * pow(stat[esp], n + 1) * ((6 * x * pow(n, -3) + pow(x, 3) / n) * BesselK1(x * n) + (3 * pow(x, 2) * pow(n, -2)) * BesselK0(x * n));
            h += (45.0 / (4 * pow(PI, 4))) * degFree[esp] * pow(stat[esp], n + 1) * ((8 * x * pow(n, -3) + pow(x, 3) / n) * BesselK1(x * n) + (4 * pow(x, 2) * pow(n, -2)) * BesselK0(x * n));
          }
        }
        // formula for massless species
        else if ((find(neutrinos_index.begin(), neutrinos_index.end(), esp) == neutrinos_index.end()) || logT > log10(M_e))
        {
          if (stat[esp] == +1) // massless bosons
          {
            g += degFree[esp];
            h += degFree[esp];
          }
          else // massless fermions
          {
            g += (7.0 / 8.0) * degFree[esp];
            h += (7.0 / 8.0) * degFree[esp];
          }
        }
        // peculiar case of neutrinos after e+/e- non relativistic
        else if ((find(neutrinos_index.begin(), neutrinos_index.end(), esp) != neutrinos_index.end()) && (logT < log10(M_e)))
        {
          //std::cout << ilogT << " Passage" << std::endl;
          g += (7.0 / 8.0) * 2 * (N_EFF / 6) * pow(4.0 / 11.0, 4.0 / 3.0);
          h += (7.0 / 8.0) * 2 * (N_EFF / 6) * (4.0 / 11.0);
        }
      }
    }

    // Writing the results in vectors
    g_vec[ilogT] = g;
    h_vec[ilogT] = h;
    logH_vec[ilogT] = log10(h);
  }
}

void DegreeFreedom::ComputeGAndHSmooth()
{

  double logTransQCD = log10(T_transQCD);

  double smoothnessQCD = 1.0 / 0.3; // smoothness of the QCD phase transition
  double smoothnessNeut = 1.0 / 0.3;
  double logT = 0, g = 0, h = 0, x; // loop variable

  // Neutrino "transition" parameters (neutrinos are species 13)
  double maxGNeut = (7.0 / 8.0) * degFree[neutrinos_index[0]];
  double maxHNeut = (7.0 / 8.0) * degFree[neutrinos_index[0]];
  double minGNeut = (7.0 / 8.0) * 2 * N_EFF / 6 * pow(4.0 / 11.0, 4.0 / 3.0);
  double minHNeut = (7.0 / 8.0) * 2 * N_EFF / 6 * (4.0 / 11.0);

  double regQCD;
  double regGNeut;
  double regHNeut;

  regQCD = regulatorTQCD(logT, smoothnessQCD, logTransQCD, 1.0, 0.0);
  regGNeut = regulatorTQCD(logT, smoothnessNeut, log10(M_e), maxGNeut, minGNeut);
  regHNeut = regulatorTQCD(logT, smoothnessNeut, log10(M_e), maxHNeut, minHNeut);

  // Initialisation of the variable g
  g = 0;
  h = 0;

  for (int ilogT = 0; ilogT < numPoints + 1; ilogT++)
  {
    // New temperature
    logT = logTmax - ilogT * deltaLogT;
    logT_vec[ilogT] = logT;

    regQCD = regulatorTQCD(logT, smoothnessQCD, logTransQCD, 1.0, 0.0);
    regGNeut = regulatorTQCD(logT, smoothnessNeut, log10(M_e), maxGNeut, minGNeut);
    regHNeut = regulatorTQCD(logT, smoothnessNeut, log10(M_e), maxHNeut, minHNeut);

    // Initialisation of the variable g
    g = 0;
    h = 0;

    // Loop on all the species
    for (int esp = 0; esp < numSpecies; esp++)
    {
      // Relevant variable x
      x = mass[esp] / pow(10, logT);

      // Only the relevant species at a given temperature are taken into account

      if (mass[esp] > 0 && (find(neutrinos_index.begin(), neutrinos_index.end(), esp) == neutrinos_index.end())) // Formula for massives species
      {
        for (int n = 1; n < nBessel; n++)
        {
          if (isPresentAfterQCDPT[esp] == true && isPresentBeforeQCDPT[esp] == true)
          {
            g += (15.0 / pow(PI, 4)) * degFree[esp] * pow(stat[esp], n + 1) * ((6 * x * pow(n, -3) + pow(x, 3) / n) * BesselK1(x * n) + (3 * pow(x, 2) * pow(n, -2)) * BesselK0(x * n));
            h += (45.0 / (4 * pow(PI, 4))) * degFree[esp] * pow(stat[esp], n + 1) * ((8 * x * pow(n, -3) + pow(x, 3) / n) * BesselK1(x * n) + (4 * pow(x, 2) * pow(n, -2)) * BesselK0(x * n));
          }
          else if (isPresentAfterQCDPT[esp] == false && isPresentBeforeQCDPT[esp] == true)
          {
            g += (15.0 / pow(PI, 4)) * degFree[esp] * regQCD * pow(stat[esp], n + 1) * ((6 * x * pow(n, -3) + pow(x, 3) / n) * BesselK1(x * n) + (3 * pow(x, 2) * pow(n, -2)) * BesselK0(x * n));
            h += (45.0 / (4 * pow(PI, 4))) * degFree[esp] * regQCD * pow(stat[esp], n + 1) * ((8 * x * pow(n, -3) + pow(x, 3) / n) * BesselK1(x * n) + (4 * pow(x, 2) * pow(n, -2)) * BesselK0(x * n));
          }
          else
          {
            g += (15.0 / pow(PI, 4)) * degFree[esp] * (1 - regQCD) * pow(stat[esp], n + 1) * ((6 * x * pow(n, -3) + pow(x, 3) / n) * BesselK1(x * n) + (3 * pow(x, 2) * pow(n, -2)) * BesselK0(x * n));
            h += (45.0 / (4 * pow(PI, 4))) * degFree[esp] * (1 - regQCD) * pow(stat[esp], n + 1) * ((8 * x * pow(n, -3) + pow(x, 3) / n) * BesselK1(x * n) + (4 * pow(x, 2) * pow(n, -2)) * BesselK0(x * n));
          }
        }
      }
      // formula for massless species
      else if ((mass[esp] == 0) && (find(neutrinos_index.begin(), neutrinos_index.end(), esp) == neutrinos_index.end()))
      {
        if (isPresentAfterQCDPT[esp] == true && isPresentBeforeQCDPT[esp] == true)
        {
          if (stat[esp] == +1)
          { // massless bosons
            g += degFree[esp];
            h += degFree[esp];
          }
          else
          { // massless fermion
            g += (7.0 / 8.0) * degFree[esp];
            h += (7.0 / 8.0) * degFree[esp];
          }
        }
        else if (isPresentAfterQCDPT[esp] == false && isPresentBeforeQCDPT[esp] == true)
        {
          if (stat[esp] == +1)
          { // massless bosons
            g += degFree[esp] * regQCD;
            h += degFree[esp] * regQCD;
          }
          else
          { // massless fermion
            g += (7.0 / 8.0) * degFree[esp] * regQCD;
            h += (7.0 / 8.0) * degFree[esp] * regQCD;
          }
        }
        else
        {
          if (stat[esp] == +1)
          { // massless bosons
            g += degFree[esp] * (1 - regQCD);
            h += degFree[esp] * (1 - regQCD);
          }
          else
          { // massless fermion
            g += (7.0 / 8.0) * degFree[esp] * (1 - regQCD);
            h += (7.0 / 8.0) * degFree[esp] * (1 - regQCD);
          }
        }
      }
      // peculiar case of neutrinos after e+/e- non relativistic
      else if ((find(neutrinos_index.begin(), neutrinos_index.end(), esp) != neutrinos_index.end()))
      {
        //std::cout << ilogT << " Passage" << std::endl;
        g += regGNeut;
        h += regHNeut;
      }
    }

    // Writing the results in vectors
    gSmooth_vec[ilogT] = g;
    hSmooth_vec[ilogT] = h;
    logHSmooth_vec[ilogT] = log10(h);
  }
}



// Function that compute g_star^{1/2}, must me called after gSmooth_vec and
// hSmooth_vec have been computed in DegreeFreedom::ComputeGAndHSmooth()
void DegreeFreedom::ComputeGStarHalf()
{
  Der = new Derivateur();

  // Keeping in memory that we computed gStarHalf_vec and derLogHSmooth_vec
  isComputeGStarHalfCalled = true;

  gStarHalf_vec.resize(numPoints + 1);
  derLogHSmooth_vec.resize(numPoints + 1);

  // Putting by hand the derivatives on the edges since they may diverge
  derLogHSmooth_vec[numPoints + 1] = 0;
  derLogHSmooth_vec[numPoints] = 0;
  derLogHSmooth_vec[0] = 0;
  derLogHSmooth_vec[1] = 0;

  for (int i = 2; i < numPoints - 1; i++)
  {
    derLogHSmooth_vec[i] = -Der->derivePrm(logHSmooth_vec, i, numPoints + 1) / deltaLogT;
    //derLogHSmooth_Lattice_ver[i] = -Der->derivePrm(logHSmooth_Lattice_vec, i, numPoints + 1) / deltaLogT;
  }

  for (int i = 0; i < numPoints + 1; i++)
  {
    gStarHalf_vec[i] = hSmooth_vec[i] * (1 + derLogHSmooth_vec[i] / 3.0) / sqrt(gSmooth_vec[i]);
  }
}

// This function uses spline to intepolate the diffrent vectors
// Must be called after the computation of the vectors to interpolate
int DegreeFreedom::Interpolate()
{
  if (isComputeGStarHalfCalled == false)
  {
    std::cout << "WARNING : Try to interpolate but the vectors to interpolate have not been computed" << std::endl;
    std::cout << "          Please run before the function DegreeFreedom::ComputeGStarHalf()" << std::endl;
    return -1;
  }

  std::clog << std::endl;
  std::clog << "|| INFO.  : (degFreeComp.) Effective degree of freedom functions are interpolated !" << std::endl;

  // This function can be called only once
  if (isInterpolateCalled == true)
    return -1;

  isInterpolateCalled = true;

  // We will interpolate derLogHSmooth_vec by a spline_derLogHSmooth
  // Before that we need to reverse the vectors since spline works only
  // for the "x" vectors ordered from min(x) to max(x)
  std::vector<double> logTReverse_vec;
  std::vector<double> gStarHalfReverse_vec, gSmoothReverse_vec, hSmoothReverse_vec, derLogHSmoothReverse_vec;

  gSmoothReverse_vec.resize(gSmooth_vec.size());
  hSmoothReverse_vec.resize(hSmooth_vec.size());
  gStarHalfReverse_vec.resize(gStarHalf_vec.size());
  derLogHSmoothReverse_vec.resize(derLogHSmooth_vec.size());
  logTReverse_vec.resize(logT_vec.size());

  derLogHSmoothReverse_vec = reverseVector(derLogHSmooth_vec);
  gSmoothReverse_vec = reverseVector(gSmooth_vec);
  hSmoothReverse_vec = reverseVector(hSmooth_vec);
  gStarHalfReverse_vec = reverseVector(gStarHalf_vec);
  logTReverse_vec = reverseVector(logT_vec);

  spline_derLogHSmooth.set_points(logTReverse_vec, derLogHSmoothReverse_vec, true);
  spline_gSmooth.set_points(logTReverse_vec, gSmoothReverse_vec, true);
  spline_hSmooth.set_points(logTReverse_vec, hSmoothReverse_vec, true);
  spline_gStarHalf.set_points(logTReverse_vec, gStarHalfReverse_vec, true);

  return 0;
}

std::vector<double> DegreeFreedom::reverseVector(std::vector<double> const &vec)
{
  int size = vec.size();
  std::vector<double> vecBis;
  vecBis.resize(size);

  for (int i = 0; i < size; i++)
  {
    vecBis[i] = vec[size - 1 - i];
  }

  return vecBis;
}

// Output function, write the ouptu files for g_eff and h_eff
void DegreeFreedom::plot_GAndH()
{
  // Output file

  out_GAndH.open("../output/DegFreedom/GAndH.out");
  out_GAndH.precision(8);

  out_GAndHSmooth.open("../output/DegFreedom/GAndHSmooth.out");
  out_GAndHSmooth.precision(8);

  out_GStarHalf.open("../output/DegFreedom/GStarHalf.out");
  out_GStarHalf.precision(8);

  for (int i = 0; i < numPoints + 1; i++)
  {
    out_GAndH << pow(10, logT_vec[i]) << "\t" << g_vec[i] << "\t" << h_vec[i] << std::endl;
  }

  for (int i = 0; i < numPoints + 1; i++)
  {
    //std::cout << i << std::endl;
    out_GAndHSmooth << pow(10, logT_vec[i]) << "\t" << gSmooth_vec[i] << "\t" << hSmooth_vec[i] << "\t" << pow(gStarHalf_vec[i], 2)
                    << "\t" << spline_gLattice(logT_vec[i]) << "\t" << spline_hLattice(logT_vec[i]) << "\t" << pow(spline_gStarHalfLattice(logT_vec[i]), 2) << std::endl;
  }

  /*
  for (int i = 0; i < numPoints + 1; i++)
  {
    out_GStarHalf << pow(10, logT_vec[i]) << "\t" << gStarHalf_vec[i] << " "
                  << pow(gStarHalf_vec[i], 2) << " " << derLogHSmooth_vec[i] << std::endl;
  }*/
}

void DegreeFreedom::InitialiseLattice()
{
  std::vector<double> log10T_MeV = {0., 0.5, 1., 1.25, 1.60, 2.00, 2.15, 2.20, 2.40, 2.50, 3.00, 4.00, 4.30, 4.60, 5.00, 5.45};
  std::vector<double> g = {10.71, 10.74, 10.76, 11.09, 13.68, 17.61, 24.07, 29.84, 47.83, 53.04, 73.48, 83.10, 84.56, 91.97, 102.17, 104.98};
  std::vector<double> g_h = {1.00228, 1.00029, 1.00048, 1.00505, 1.02159, 1.02324, 1.05423, 1.07578, 1.06118, 1.04690, 1.01778, 1.00123, 1.00389, 1.00887, 1.00750, 1.00023};

  std::vector<double> h, log10T_GeV;
  h.resize(g_h.size());
  log10T_GeV.resize(log10T_MeV.size());

  for (int i = 0; i < g_h.size(); i++)
    h[i] = g[i] / g_h[i];

  for (int i = 0; i < log10T_MeV.size(); i++)
    log10T_GeV[i] = log10T_MeV[i] - 3;

  spline_gLattice.set_points(log10T_GeV, g, true);
  spline_hLattice.set_points(log10T_GeV, h, true);

  std::vector<double> gs, logHSmooth_Lattice_vec, derLogHSmooth_Lattice_vec;
  gs.resize(numPoints + 1);
  logHSmooth_Lattice_vec.resize(numPoints + 1);
  derLogHSmooth_Lattice_vec.resize(numPoints + 1);
  // Setting the initial and final values
  derLogHSmooth_Lattice_vec[numPoints + 1] = 0;
  derLogHSmooth_Lattice_vec[numPoints] = 0;
  derLogHSmooth_Lattice_vec[0] = 0;
  derLogHSmooth_Lattice_vec[1] = 0;

  for (int i = 0; i < numPoints + 1; i++)
    logHSmooth_Lattice_vec[i] = log10(spline_hLattice(logT_vec[i]));

  for (int i = 2; i < numPoints - 1; i++)
  {
    derLogHSmooth_Lattice_vec[i] = -Der->derivePrm(logHSmooth_Lattice_vec, i, numPoints + 1) / deltaLogT;
    //derLogHSmooth_Lattice_vec[i] = (logHSmooth_Lattice_vec[i+1] - logHSmooth_Lattice_vec[i])/ (logT_vec[i+1] - logT_vec[i]);
  }

  for (int i = 0; i < numPoints + 1; i++)
    gs[i] = spline_hLattice(logT_vec[i]) * (1 + derLogHSmooth_Lattice_vec[i] / 3.0) / sqrt(spline_gLattice(logT_vec[i]));

  std::vector<double> gs_reverse = reverseVector(gs), logT_reverse = reverseVector(logT_vec);

  spline_gStarHalfLattice.set_points(logT_reverse, gs_reverse, true);
}