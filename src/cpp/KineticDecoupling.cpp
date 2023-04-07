#include "../headers/KineticDecoupling.h"

KineticDecoupling::KineticDecoupling(DarkSectorModel const &_DS, DegreeFreedom *_degFree, int index_chi)
{
    chi = _DS.get_chi(index_chi);

    double n_scatter_SM = 24; // BE CAREFUL JUST CHANGED
    double n_scatter_phi = 0;

    n_scatter = n_scatter_SM;
    crossSectionScatt.resize(0);
    mass_scattered_part.resize(0);
    spin_scattered_part.resize(0);

    // Here we put the n_scatter to 24 as we can diffuse over the 12 particles of the SM
    // as well as the 12 antiparticles associated -> For a simple t-channel the diffusion
    // against the particle and the anti-particle have the same cross-section
    for (int i = 0; i < n_scatter_SM; i++)
    {
        //std::cout << "Check for the segmentation fault beg : " << i << std::endl;
        crossSectionScatt.push_back(new CrossSection_cfcf(_DS, chi, StandardModel::getInstance()->get_couples_SM_ferm(i % 12)));
        crossSectionScatt.push_back(new CrossSection_cfcf(_DS, chi, StandardModel::getInstance()->get_couples_SM_ferm(i % 12)));
        //std::cout << "Check for the segmentation fault : " << i << std::endl;
        mass_scattered_part.push_back(StandardModel::getInstance()->get_couples_SM_ferm(i % 12).get_first().get_mass());
        mass_scattered_part.push_back(StandardModel::getInstance()->get_couples_SM_ferm(i % 12).get_second().get_mass());
        spin_scattered_part.push_back(StandardModel::getInstance()->get_couples_SM_ferm(i % 12).get_first().get_spin());
        spin_scattered_part.push_back(StandardModel::getInstance()->get_couples_SM_ferm(i % 12).get_second().get_spin());
        //std::cout << "Check for the segmentation fault end : " << i << std::endl;
    }

    // Cross-section scattering with photons -> One loop order
    crossSectionScatt_photons = new CrossSection_cgacga(_DS, chi);

    // If mediators are less massive than the DM mass they can be produced
    // we assume that they acquire an equilibrium distribution in the plasma
    // therefore we add their contribution to the kinetic equilibrium
    /*std::vector<Particle> phis = _DS.get_phis();
    for (int i = 0; i < phis.size(); i++)
        if (phis[i].get_mass() < chi.get_mass())
        {
            n_scatter_phi++;
            crossSectionScatt.push_back(new CrossSection_cscs(_DS, chi, phis[i], chi, phis[i]));
            mass_scattered_part.push_back(phis[i].get_mass());
            spin_scattered_part.push_back(phis[i].get_spin());
        }*/
    // Attention : in practice we could also scatter on two different scalars
    // In which case we may have to revise the expression for the relaxation momentum rate


    spline_gSmooth = _degFree->get_gSmoothInterp();
    spline_hSmooth = _degFree->get_hSmoothInterp();
    spline_derLogHSmooth = _degFree->get_derLogHSmoothInterp();
    spline_gStarHalf = _degFree->get_gStarHalfInterp();

    //std::cout << "Introduction : " << spline_gSmooth(-4) << std::endl;

    T_transQCD = _degFree->get_T_transQCD();

    InterpolateGammaTot();
}


double KineticDecoupling::gammaTot(double T)
{
    double res = 0;

    // Regulator : Theta(log10(T)-log10(T_QCD_GEV)) smoothed
    double reg = regulatorTQCD(log10(T), 1.0 / 0.3, log10(T_transQCD), 1, 0);

    for (int i = 0; i < n_scatter; i++)

        if (crossSectionScatt[i]->get_particle_input_2().get_QCDType() == QCDTYPE::QUARK)
        {
            // When dealing with quarks in the out states we only consider them if the temperature
            // is above the temperature of the QCD phase transition -> regulated by the regulator function
            // Moreover we multply by 3 in order to take into account the 3 colors
            res += reg * 3 * gamma(i, T);
        }
        else if (crossSectionScatt[i]->get_particle_input_2().get_QCDType() != QCDTYPE::QUARK)
        {
            res += gamma(i, T);
        }

    return res;
}

double KineticDecoupling::gamma(int index_cs, double T)
{
    std::vector<double> omiT = {0, (double)index_cs, T};

    double m2 = mass_scattered_part[index_cs];

    if (m2 > 100 * T)
        return 0;

    return GaussLegendre_IntegralLn_Static(0, 3, gl500, m2, 100 * T, omiT, this, CallBack_fToIntForGamma) / (3 * PI * PI * T * chi.get_degFree() * chi.get_mass());
}



double KineticDecoupling::fToIntForGamma(double om, double index_cs, double T)
{
    int eps = 0;

    if (spin_scattered_part[(int)index_cs] == floor(spin_scattered_part[(int)index_cs])) // check if it is a boson
        eps = -1;
    else // if it is a fermion
        eps = +1;

    //std::cout << eps << std::endl;

    return exp(om / T) * pow(exp(om / T) + eps, -2) * crossSectionScatt[(int)index_cs]->EvalSigmaTransferBringmannP1cm4_fromIntegral(om);
    //return  exp(om / T) * pow(exp(om / T) + eps, -2);
}

double KineticDecoupling::gamma_1Loop_photons(double T)
// Leading one loop corrections to the value of gamma
{
    std::vector<double> omiT = {0, T};

    return GaussLegendre_IntegralLn_Static(0, 2, gl1000, 0.01*T, 100 * T, omiT, this, CallBack_fToIntForGamma_1Loop_photons) / (3 * PI * PI * T * chi.get_degFree() * chi.get_mass());
}

double KineticDecoupling::fToIntForGamma_1Loop_photons(double om,  double T)
{
    int eps = -1; // Here we deal with hotons that are bosons
    // Here we use the computation from the numerical integral in order to account for the dependance in $t$ in the coupling constant
    return exp(om / T) * pow(exp(om / T) + eps, -2) * crossSectionScatt_photons->EvalSigmaTransferBringmannP1cm4_fromIntegral(om);
}

void KineticDecoupling::InterpolateGammaTot()
{

    int numPoints = 250;

    std::vector<double> log10x;
    std::vector<double> gamma_vec;
    log10x.resize(numPoints + 1);
    gamma_vec.resize(numPoints + 1);

    double m_chi = chi.get_mass();
    double xmin = 1.;
    double xmax = 1e+9;

    double x = 0, dx = log10(xmax / xmin) / numPoints;

    for (int i = 0; i < numPoints + 1; i++)
    {
        //std::cout << "here " << T_min << std::endl;
        x = xmin * pow(10, i * dx);
        log10x[i] = log10(x);
        gamma_vec[i] = gammaTot(m_chi / x);
        //std::cout <<  log10x[i] << " " << m_chi/x << " " << gamma_vec[i]  << std::endl;
    }

    //std::cout << "Interpolation of <sigma v_Mol> completed" << std::endl;

    spline_gammaTot.set_points(log10x, gamma_vec, false);
}

double KineticDecoupling::SolveTemperatureEI()
{

    double z, wbis, yeq; // w = ln(TDM/T), z = ln(x) = ln(massDM/T), wbis = ln(massDM TDM/s^2/3)
    double yEIp1, yEI;   // variables for Implicit Euler method
    std::vector<double> varEI;

    varEI.resize(2);

    xSaved.clear();
    ySavedEIL.clear();
    yeqSaved.clear();

    double dz = 0.005;
    double massDM = chi.get_mass();
    double Tkd = EstimateDecouplingTemperature();

    if (Tkd == -2)
        return -2; // means that it is not possible to solve
    else
        std::cout << "Ici Tkd = " << Tkd << std::endl;

    double xmax = std::min(1e+9, massDM / 1e-5);

    double zMax = std::min(log(massDM / Tkd) + 4, log(xmax));
    double zMin = log(massDM / Tkd) - 4;

    int compt = 0, comptMax = 100000;

    z = zMin;

    yeq = exp(z) * pow(spline_hSmooth(log10(massDM * exp(-z))), -2. / 3.) * pow(45 / (2 * PI * PI), 2. / 3.);

    // Initialization of temperature
    //wbis = log(yeq);
    yEI = yeq;

    xSaved.push_back(exp(z));
    ySavedEIL.push_back(yEI);

    yeqSaved.push_back(yeq);

    bool problemDetected;
    int numTry = 10;

    double wbisp1;
    double approxDeriv;

    // std::cout << pow((2 * PI * PI / 45.), 2. / 3.) << std::endl;
    for (int j = 0; j < numTry; j++)
    {
        problemDetected = false;

        while (z < zMax && compt < comptMax && problemDetected == false)
        {
            z += dz;
            compt++;

            varEI[0] = z;

            yeq = exp(z) * pow(spline_hSmooth(log10(massDM * exp(-z))), -2. / 3.) * pow(45 / (2 * PI * PI), 2. / 3.);

            varEI[1] = yeq;

            // Second using Euler implicit
            yEIp1 = ImplicitEulerSolverLinear_Static(0, varEI, yEI, this, CallBackFunctionForImplicitEulerLinearA, CallBackFunctionForImplicitEulerLinearZ, dz);

            approxDeriv = 0; //abs(yEIp1 - yEI) / dz;
            //std::cout << abs(yEIp1 - yEI) / dz << std::endl;
            yEI = yEIp1;

            // Saving the results for plot if needed
            xSaved.push_back(exp(z));
            ySavedEIL.push_back(yEI);
            yeqSaved.push_back(yeq);

            // Check if there is a problem in the resolution
            //if ((exp(wbis) != exp(wbis) || isinf(exp(wbis)) == true || approxDeriv > 100) && j < numTry - 1)
            if ((yEI != yEI || isinf(yEI) == true || approxDeriv > 1000) && j < numTry - 1)
            {

                //exit(0);
                problemDetected = true;
                compt = 0;
                zMin += 0.5;

                z = zMin;

                xSaved.clear();
                ySavedEIL.clear();
                yeqSaved.clear();

                yeq = exp(z) * pow(spline_hSmooth(log10(massDM * exp(-z))), -2. / 3.) * pow(45 / (2 * PI * PI), 2. / 3.);

                // Initialization of temperature
                wbis = log(yeq);
                yEI = yeq;

                xSaved.push_back(exp(z));
                ySavedEIL.push_back(yEI);
                yeqSaved.push_back(yeq);

                std::cout << "WARNING : Problem detected in resolution for thermal decoupling" << std::endl;
                std::cout << "          We try again with  log10(xMin) = " << zMin << " -> xMin = " << exp(zMin) << std::endl;
            }
            // else if (exp(wbis) != exp(wbis) || isinf(exp(wbis)) || approxDeriv > 100) == true)
            else if ((yEI != yEI || isinf(yEI) || approxDeriv > 1000) == true)
            {
                problemDetected = true;
                std::cout << "ERROR : Impossible to solve thermal decoupling" << std::endl;
            }
        }
    }

    yinf = yEI;

    double TkdApproxNew = Dichotomie_Static(this, CallBackFToCheckDecTem, massDM * exp(-zMin), massDM * exp(-zMax), 1e-3, 0);

    std::cout << "Approximate x dec (final eval.) : " << massDM / TkdApproxNew << std::endl;

    return TkdApproxNew;
}

double KineticDecoupling::FunctionForImplicitEulerLinearA(std::vector<double> var)
{

    double z = var[0];
    double massDM = chi.get_mass();
    double T = massDM * exp(-z);

    double fact = (1 + spline_derLogHSmooth(log10(T)) / 3.0);
    double H = 2 * sqrt((pow(PI, 3) / 45) * spline_gSmooth(log10(T))) * pow(T, 2) / M_PLANCK;

    return fact * spline_gammaTot(log10(massDM / T)) / H;
}

double KineticDecoupling::FunctionForImplicitEulerLinearZ(std::vector<double> var)
{
    return var[1];
}

double KineticDecoupling::SolveTemperature()
{

    double z, wbis, yeq;                      // w = ln(TDM/T), z = ln(x) = ln(massDM/T), wbis = ln(massDM TDM/s^2/3)
    double yAMp1, yAM, yAMm1, yAMm2, yAMm3;   // variables for Adam moulton method of order 5
    double yBDp1, yBDm1, yBDm2, yBDm3, yBDm4; // variables for Backward differential method of order 5
    double yEIp1, yEI, yBD;                   // variables for Implicit Euler method
    std::vector<double> varEI, varAM, varBD;
    std::vector<double> varAM0, varAMm1, varAMm2, varAMm3;

    varEI.resize(2);

    varAM.resize(2);
    varAM0.resize(2);
    varAMm1.resize(2);
    varAMm2.resize(2);
    varAMm3.resize(2);

    varBD.resize(2);

    xSaved.clear();
    //ySavedRK6.clear();
    ySavedBDL.clear();
    ySavedAML.clear();
    ySavedEIL.clear();
    yeqSaved.clear();

    double dz = 0.005;
    double massDM = chi.get_mass();
    double Tkd = EstimateDecouplingTemperature();
    double zMax = log(massDM / Tkd) + 2;
    double zMin = log(massDM / Tkd) - 10;

    int compt = 0, comptMax = 100000;

    z = zMin;

    yeq = exp(z) * pow(spline_hSmooth(log10(massDM * exp(-z))), -2. / 3.) *
          pow(45 / (2 * PI * PI), 2. / 3.);

    // Initialization of temperature
    //wbis = log(yeq);
    yEI = yeq;
    yAM = yeq;
    yBD = yeq;

    yAMm1 = exp(z - dz) * pow(spline_hSmooth(log10(massDM * exp(-z + dz))), -2. / 3.) *
            pow(45 / (2 * PI * PI), 2. / 3.);
    yAMm2 = exp(z - 2 * dz) * pow(spline_hSmooth(log10(massDM * exp(-z + 2 * dz))), -2. / 3.) *
            pow(45 / (2 * PI * PI), 2. / 3.);
    yAMm3 = exp(z - 3 * dz) * pow(spline_hSmooth(log10(massDM * exp(-z + 3 * dz))), -2. / 3.) *
            pow(45 / (2 * PI * PI), 2. / 3.);

    yBDm1 = exp(z - dz) * pow(spline_hSmooth(log10(massDM * exp(-z + dz))), -2. / 3.) *
            pow(45 / (2 * PI * PI), 2. / 3.);
    yBDm2 = exp(z - 2 * dz) * pow(spline_hSmooth(log10(massDM * exp(-z + 2 * dz))), -2. / 3.) *
            pow(45 / (2 * PI * PI), 2. / 3.);
    yBDm3 = exp(z - 3 * dz) * pow(spline_hSmooth(log10(massDM * exp(-z + 3 * dz))), -2. / 3.) *
            pow(45 / (2 * PI * PI), 2. / 3.);
    yBDm4 = exp(z - 4 * dz) * pow(spline_hSmooth(log10(massDM * exp(-z + 4 * dz))), -2. / 3.) *
            pow(45 / (2 * PI * PI), 2. / 3.);

    xSaved.push_back(exp(z));
    //ySavedRK6.push_back(exp(wbis));
    ySavedAML.push_back(yAM);
    ySavedEIL.push_back(yEI);
    ySavedBDL.push_back(yBD);
    yeqSaved.push_back(yeq);

    bool problemDetected;
    int numTry = 10;

    double approxDeriv;

    // std::cout << pow((2 * PI * PI / 45.), 2. / 3.) << std::endl;
    for (int j = 0; j < numTry; j++)
    {
        problemDetected = false;

        while (z < zMax && compt < comptMax && problemDetected == false)
        {
            z += dz;
            compt++;

            varEI[0] = z;
            varAM[0] = z;
            varBD[0] = z;

            // Strange but 0 is one step before the step +1 currntly calculated
            // Therefore 0 correspond to the "z" : z - dz and so on and so forth
            // This is because z += dz is above
            varAM0[0] = z - dz;
            varAMm1[0] = z - 2. * dz;
            varAMm2[0] = z - 3. * dz;
            varAMm3[0] = z - 4. * dz;

            yeq = exp(z) * pow(spline_hSmooth(log10(massDM * exp(-z))), -2. / 3.) * pow(45 / (2 * PI * PI), 2. / 3.);

            varEI[1] = yeq;
            varAM[1] = yeq;
            varBD[1] = yeq;

            varAM0[1] = exp(z - dz) * pow(spline_hSmooth(log10(massDM * exp(-z + dz))), -2. / 3.) *
                        pow(45 / (2 * PI * PI), 2. / 3.);
            varAMm1[1] = exp(z - 2 * dz) * pow(spline_hSmooth(log10(massDM * exp(-z + 2 * dz))), -2. / 3.) *
                         pow(45 / (2 * PI * PI), 2. / 3.);
            varAMm2[1] = exp(z - 3 * dz) * pow(spline_hSmooth(log10(massDM * exp(-z + 3 * dz))), -2. / 3.) *
                         pow(45 / (2 * PI * PI), 2. / 3.);
            varAMm3[1] = exp(z - 4 * dz) * pow(spline_hSmooth(log10(massDM * exp(-z + 4 * dz))), -2. / 3.) *
                         pow(45 / (2 * PI * PI), 2. / 3.);

            // First method using RK6
            //wbisp1 = myRungeKutta6_Static(0, var, wbis, this, CallBackDerivativeOfTemperature,
            //                                 varstep, dz, dzmin, prec, okprec);

            //approxDeriv = abs(wbisp1 - wbis) / dz;
            //wbis = wbisp1;

            // Second Method using Euler implicit
            yBDp1 = ImplicitBackwardDiffSolver5Linear_Static(0, varBD, yBD, yBDm1, yBDm2, yBDm3, yBDm4, this, CallBackFunctionForImplicitEulerLinearA, CallBackFunctionForImplicitEulerLinearZ, dz);
            yAMp1 = ImplicitAdamsMoultonSolver5Linear_Static(0, varAM, varAM0, varAMm1, varAMm2, varAMm3, yAM, yAMm1, yAMm2, yAMm3, this, CallBackFunctionForImplicitEulerLinearA, CallBackFunctionForImplicitEulerLinearZ, dz);
            yEIp1 = ImplicitEulerSolverLinear_Static(0, varEI, yEI, this, CallBackFunctionForImplicitEulerLinearA, CallBackFunctionForImplicitEulerLinearZ, dz);

            yAMm3 = yAMm2;
            yAMm2 = yAMm1;
            yAMm1 = yAM;
            yAM = yAMp1;

            yBDm4 = yBDm3;
            yBDm3 = yBDm2;
            yBDm2 = yBDm1;
            yBDm1 = yBD;
            yBD = yBDp1;

            approxDeriv = abs(yEIp1 - yEI) / dz;
            yEI = yEIp1;

            // Saving the results for plot if needed
            xSaved.push_back(exp(z));
            //ySavedRK6.push_back(exp(wbis));
            ySavedEIL.push_back(yEI);
            ySavedBDL.push_back(yBD);
            ySavedAML.push_back(yAM);
            yeqSaved.push_back(yeq);

            // Check if there is a problem in the resolution
            //if ((exp(wbis) != exp(wbis) || isinf(exp(wbis)) == true || approxDeriv > 100) && j < numTry - 1)
            if ((yEI != yEI || isinf(yEI) == true || approxDeriv > 100) && j < numTry - 1)
            {
                problemDetected = true;
                compt = 0;
                zMin += 0.5;

                z = zMin;

                xSaved.clear();
                // ySavedRK6.clear();
                ySavedEIL.clear();
                ySavedBDL.clear();
                ySavedAML.clear();
                yeqSaved.clear();

                yeq = exp(z) * pow(spline_hSmooth(log10(massDM * exp(-z))), -2. / 3.) *
                      pow(45 / (2 * PI * PI), 2. / 3.);

                // Initialization of temperature
                wbis = log(yeq);
                yEI = yeq;
                yAM = yeq;
                yBD = yeq;

                xSaved.push_back(exp(z));
                // ySavedRK6.push_back(exp(wbis));
                ySavedEIL.push_back(yEI);
                ySavedBDL.push_back(yBD);
                ySavedAML.push_back(yAM);
                yeqSaved.push_back(yeq);

                std::cout << "WARNING : Problem detected in resolution for thermal decoupling" << std::endl;
                std::cout << "          We try again with  log10(xMin) = " << zMin << " -> xMin = " << exp(zMin) << std::endl;
            }
            // else if (exp(wbis) != exp(wbis) || isinf(exp(wbis)) || approxDeriv > 100) == true)
            else if ((yEI != yEI || isinf(yEI) || approxDeriv > 100) == true)
            {
                problemDetected = true;
                std::cout << "ERROR : Impossible to solve thermal decoupling" << std::endl;
            }
        }
    }

    yinf = yEI;

    double TkdApproxNew = Dichotomie_Static(this, CallBackFToCheckDecTem, massDM * exp(-zMin), massDM * exp(-zMax), 1e-3, 0);

    std::cout << "Approximate x dec (final eval.) : " << massDM / TkdApproxNew << std::endl;

    return TkdApproxNew;
}

void KineticDecoupling::DerivativeOfTemperature(std::vector<double> varBis, double w, double &deriv)
{
    double z = varBis[0];
    double yeq = varBis[1];
    double massDM = chi.get_mass();
    double T = massDM * exp(-z);
    double H = 2 * sqrt((pow(PI, 3) / 45) * spline_gSmooth(log10(T))) * pow(T, 2) / M_PLANCK;

    //     std::cout << spline_gSmooth(log10(T)) << std::endl;

    double fact = (1 + spline_derLogHSmooth(log10(T)) / 3.0);

    //deriv = 1 + ((gammaTot(T) / H) * (exp(-w) - 1) - 2) * fact;
    deriv = -fact * (gammaTot(T) / H) * (1 - yeq * exp(-w));
}

void KineticDecoupling::plotTemperatureEvolution()
{
    std::ofstream out_TDM;
    out_TDM.open("../output/TDM");
    out_TDM.precision(8);

    std::ofstream out_y;

    std::string file_name = "yChiVsxChiTestEIL";
    std::ostringstream strc;
    //strc << "_LSDM" << 100 * couplingScalarDM << "_LSF" << 100 * couplingScalarFermions << "_LPSDM"
    //     << 100 * couplingPScalarDM << "_LPSF" << 100 * couplingPScalarFermions << "_massDM" << massDM;
    file_name += strc.str();

    out_y.open("../output/" + file_name);

    out_y.precision(8);

    int nsize = ySavedEIL.size();
    for (int i = 0; i < nsize; i++)
    {
        out_y << xSaved[i] << "\t " << ySavedBDL[i] << "\t" << ySavedEIL[i] << "\t" << ySavedBDL[i] << "\t" << yeqSaved[i] << std::endl;
    }
}

void KineticDecoupling::plotTemperatureEvolutionEI(std::string name)
{

    std::ofstream out_y;
    out_y.open("../output/TemperatureEvolution/yChiVsxChi_" + name + ".out");
    out_y.precision(8);

    out_y << "# x | y | yeq | gammaTot [s^{-1}] | H [s^{-1}]" << std::endl;

    int nsize = ySavedEIL.size();

    double HH = 0, mchi = chi.get_mass(), T = 0;
    for (int i = 0; i < nsize; i++)
    {
        T = mchi / xSaved[i];
        HH = 2 * sqrt((pow(PI, 3) / 45) * spline_gSmooth(log10(T))) * pow(T, 2) / M_PLANCK;
        out_y << xSaved[i] << "\t "
              << "\t" << ySavedEIL[i]
              << "\t" << yeqSaved[i]
              << "\t" << spline_gammaTot(log10(xSaved[i])) * GeV_to_sm1
              << "\t" << HH * GeV_to_sm1 << std::endl;
    }
}

/*
// Be carefull, need to have created a vector xSaved before using this function
void KineticDecoupling::PlotGammaKineticVsTemperature()
{
    std::ofstream out_Gamma;

    std::string file_name = "GammaKineticVsTemperature";
    out_Gamma.open("../output/Gamma/" + file_name);
    out_Gamma.precision(8);

    double H, Gamma, T;

    int nsize = xSaved.size();
    
    for (int i = 0; i < nsize; i++)
    {
        T = massDM / xSaved[i];
        Gamma = gammaTot(T);
        H = 2 * sqrt((pow(PI, 3) / 45) * spline_gSmooth(log10(T))) * pow(T, 2) / M_PLANCK;
        out_Gamma << xSaved[i] << "\t " << Gamma << "\t" << H << "\t" << Gamma / H << "\t" << (1 + spline_derLogHSmooth(log10(T)) / 3.0) * Gamma / H << std::endl;
    }
}*/

double KineticDecoupling::EstimateDecouplingTemperature()
{
    double TkdApprox;
    double massDM = chi.get_mass();
    double Tmin = std::max(massDM / 1e+9, 1e-5);
    double Tmax = massDM / 1.;
    std::vector<double> TT{1.};

    TkdApprox = Dichotomie_Static(0, 1, TT, this, CallBackFToEstimateDecTem, Tmin, Tmax, 1e-3, -1);

    TT[0] = sqrt(Tmin * Tmax);

    //std::cout << "There : " << FToEstimateDecTem(TT[0]) << std::endl;

    //double m_chi = chi.get_mass();
    //double HH = 2 * sqrt((pow(PI, 3) / 45) * spline_gSmooth(log10(TT[0]))) * pow(TT[0], 2) / M_PLANCK; // in GeV

    //std::cout << TT[0] << " " << m_chi/TT[0] << " " << HH*GeV_to_sm1 << " " << spline_gammaTot(log10(m_chi / TT[0]))*GeV_to_sm1 << std::endl;

    if (TkdApprox == -1 && FToEstimateDecTem(TT[0]) < 0)
    {
        // Throw an exception if it is impossible to evalate the decouping temperature because the kinetic decoupling
        ExceptionHandler e((char *)"WARNING :: KineticDecoupling gamma > H exception", __LINE__, (char *)__PRETTY_FUNCTION__, ExceptionType::KineticDec_gammaGreatThanH);
        throw e;
    }
    else if (TkdApprox == -1 && FToEstimateDecTem(TT[0]) > 0)
    {
        ExceptionHandler e((char *)"WARNING :: KineticDecoupling gamma < H exception", __LINE__, (char *)__PRETTY_FUNCTION__, ExceptionType::KineticDec_gammaLessThanH);
        throw e;
    }
    else
        std::cout << "Approximate x dec (first eval.) : " << massDM / TkdApprox << std::endl;

    return TkdApprox;
}

double KineticDecoupling::FToEstimateDecTem(double T)
{
    //std::cout << T << std::endl;
    //std::cout << spline_gSmooth(-3) << std::endl;
    double m_chi = chi.get_mass();
    double H = 2 * sqrt((pow(PI, 3) / 45) * spline_gSmooth(log10(T))) * pow(T, 2) / M_PLANCK; // in GeV

    return H / spline_gammaTot(log10(m_chi / T)) - 1; // in GeV
}

double KineticDecoupling::FToCheckDecTem(double T)
{

    double massDM = chi.get_mass();
    double a = pow((2 * PI * PI / 45.) * spline_hSmooth(log10(T)), 2. / 3.);

    return massDM / T - yinf * a;
}

//
//
// Plotting functions
//
//

void KineticDecoupling::plotGammaTot()
{
    std::ofstream outfile;
    std::string nameTot = "../output/Gamma/GammaTot.out";
    outfile.open(nameTot);
    outfile.precision(8);

    double T_min = std::max(chi.get_mass() * 1e-8, 1e-6);
    double T_max = chi.get_mass();

    //std::cout << chi.get_mass() << std::endl;

    int Npts = 100;
    double logdel_T = (log10(T_max) - log10(T_min)) / (1.0 * Npts);

    double T, H;

    outfile << "# mchi = " << chi.get_mass() << "GeV" << std::endl;
    outfile << "# T [GeV] \t Gamma(T) [GeV] | H [GeV] | x = m_chi/T" << std::endl;
    outfile << std::endl;

    for (int i = 0; i < Npts + 1; i++)
    {
        T = T_min * pow(10, i * logdel_T);
        H = 2 * sqrt((pow(PI, 3) / 45) * spline_gSmooth(log10(T))) * pow(T, 2) / M_PLANCK;
        outfile << T << "\t" << gammaTot(T) << "\t" << "\t" << H << "\t" << chi.get_mass() / T << std::endl; //spline_gammaTot(log10(chi.get_mass() / T))<< std::endl;
    }

    outfile.close();
}


void KineticDecoupling::plotGamma_1Loop_photons()
{
    std::ofstream outfile;
    std::string nameTot = "../output/Gamma/Gamma_1Loop_photons.out";
    outfile.open(nameTot);
    outfile.precision(8);

    double T_min = std::max(chi.get_mass() * 1e-8, 1e-6);
    double T_max = chi.get_mass();

    //std::cout << chi.get_mass() << std::endl;

    int Npts = 100;
    double logdel_T = (log10(T_max) - log10(T_min)) / (1.0 * Npts);

    double T, H;

    outfile << "# mchi = " << chi.get_mass() << "GeV" << std::endl;
    outfile << "# T [GeV] \t Gamma(T) [GeV] | H [GeV] | x = m_chi/T" << std::endl;
    outfile << std::endl;

    for (int i = 0; i < Npts + 1; i++)
    {
        T = T_min * pow(10, i * logdel_T);
        H = 2 * sqrt((pow(PI, 3) / 45) * spline_gSmooth(log10(T))) * pow(T, 2) / M_PLANCK;
        outfile << T << "\t" << gamma_1Loop_photons(T) << "\t" << "\t" << H << "\t" << chi.get_mass() / T << std::endl;
    }

    outfile.close();
}


// Plotting function to plot the cross section with s
void KineticDecoupling::plotfToIntForGamma(int index_cs, double T)
{
    std::ofstream outfile;
    std::string nameTot = "../output/fToIntForGamma.out";
    outfile.open(nameTot);
    outfile.precision(15);

    double om_min = mass_scattered_part[index_cs];

    if (om_min == 0)
        om_min = 1e-10;

    double om_max = 100 * T;

    if (om_min > om_max)
        std::cout << "impossible to plot fToIntForGamma" << std::endl;

    int Npts = 1000;
    double logdel_om = (log10(om_max) - log10(om_min)) / (1.0 * Npts);

    double om, s;

    outfile << "# om [GeV] | fToIntForGamma(om) [GeV^{-1} s^{-1}] | s [GeV^2]" << std::endl;
    outfile << std::endl;

    for (int i = 0; i < Npts; i++)
    {
        om = om_min * pow(10, i * logdel_om);
        s = pow(chi.get_mass(), 2) + 2 * om * chi.get_mass() + pow(mass_scattered_part[index_cs], 2);
        outfile << om << "\t" << fToIntForGamma(om, (double)index_cs, T) << "\t" << s << std::endl;
    }

    outfile.close();
}

void KineticDecoupling::plotfToEstimateDecTem()
{
    std::ofstream outfile;
    std::string nameTot = "../output/fToEstimateDecTem.out";
    outfile.open(nameTot);
    outfile.precision(15);

    double T_min = std::max(chi.get_mass() / 10000, 1e-5);
    double T_max = chi.get_mass() / 1.;

    int Npts = 1000;
    double logdel_T = (log10(T_max) - log10(T_min)) / (1.0 * Npts);

    double T;

    outfile << "# T [GeV] | fToEstimateDecTem" << std::endl;
    outfile << std::endl;

    for (int i = 0; i < Npts; i++)
    {
        T = T_min * pow(10, i * logdel_T);
        outfile << T << "\t" << FToEstimateDecTem(T) << std::endl;
    }

    outfile.close();
}
