#include "../headers/StandardModel.h"

StandardModel::StandardModel()
{

    n_SM_Elem = 30;

    // ---- SM particles ----

    // Mediator boson
    photon = Particle(1., 2, "photon", 0, 0, 0, 0);
    gluon = Particle(1., 16, "gluon", 0, 0, 0, 1, QCDTYPE::GLUON);
    higgs = Particle(0., 1, "higgs", M_H, W_H, 0, 2);
    boson_Z = Particle(1., 3, "boson_Z", M_Z, W_Z, 0, 3);
    boson_Wp = Particle(1., 3, "boson_Wp", M_W, W_W, +1., 4);
    boson_Wm = Particle(1., 3, "boson_Wm", M_W, W_W, -1., 5);

    //Neutrinos
    neutrino_e = Particle(0.5, 1, "neut_e", M_nue, W_nue, 0, 6);
    anti_neutrino_e = Particle(0.5, 1, "anti_neut_e", M_nue, W_nue, 0, 7);
    neutrino_mu = Particle(0.5, 1, "neut_mu", M_numu, W_numu, 0, 8);
    anti_neutrino_mu = Particle(0.5, 1, "anti_neut_mu", M_numu, W_numu, 0, 9);
    neutrino_tau = Particle(0.5, 1, "neut_tau", M_nutau, W_nutau, 0, 10);
    anti_neutrino_tau = Particle(0.5, 1, "anti_neut_tau", M_nutau, W_nutau, 0, 11);

    // Leptons
    electron = Particle(0.5, 2, "electron", M_e, W_e, -1., 12);
    positron = Particle(0.5, 2, "positron", M_e, W_e, +1., 13);
    muon = Particle(0.5, 2, "muon", M_MU, W_MU, -1., 14);
    anti_muon = Particle(0.5, 2, "anti_muon", M_MU, W_MU, +1., 15);
    tau = Particle(0.5, 2, "tau", M_TAU, W_TAU, -1., 16);
    anti_tau = Particle(0.5, 2, "anti_tau", M_TAU, W_TAU, +1., 17);

    // Quarks
 
    quark_down = Particle(0.5, 6, "quark_down", M_QUARK_DOWN, 0, -1./3., 18, QCDTYPE::QUARK);
    anti_quark_down = Particle(0.5, 6, "anti_quark_down", M_QUARK_DOWN, 0, +1./3., 19, QCDTYPE::QUARK);
    quark_up = Particle(0.5, 6, "quark_up", M_QUARK_UP, 0, +2./3., 20, QCDTYPE::QUARK);
    anti_quark_up = Particle(0.5, 6, "anti_quark_up", M_QUARK_UP, 0, -2./3., 21, QCDTYPE::QUARK);
    quark_strange = Particle(0.5, 6, "quark_str", M_QUARK_STRANGE, 0, -1./3., 22, QCDTYPE::QUARK);
    anti_quark_strange = Particle(0.5, 6, "anti_quark_str", M_QUARK_STRANGE, 0, +1./3., 23, QCDTYPE::QUARK);
    quark_charm = Particle(0.5, 6, "quark_charm", M_QUARK_CHARM, 0, +2./3., 24, QCDTYPE::QUARK);
    anti_quark_charm = Particle(0.5, 6, "anti_quark_charm", M_QUARK_CHARM, 0, -2./3., 25, QCDTYPE::QUARK);
    quark_bottom = Particle(0.5, 6, "quark_bottom", M_QUARK_BOTTOM, 0, -1./3., 26, QCDTYPE::QUARK);
    anti_quark_bottom = Particle(0.5, 6, "anti_quark_bottom", M_QUARK_BOTTOM, 0, +1./3., 27, QCDTYPE::QUARK);
    quark_top = Particle(0.5, 6, "quark_top", M_QUARK_TOP, 0, 2./3.,28, QCDTYPE::QUARK);
    anti_quark_top = Particle(0.5, 6, "anti_quark_top", M_QUARK_TOP, 0, -2./3.,29, QCDTYPE::QUARK);

    // Composite particles

    n_SM_Poly = 11;

    pion0 = Particle(0., 1, "pion0", M_PI_0, W_PI_0, 0, 30, QCDTYPE::MESON);
    pionp = Particle(0., 1, "pionp", M_PI_C, W_PI_C, +1., 31, QCDTYPE::MESON);
    pionm = Particle(0., 1, "pionm", M_PI_C, W_PI_C, -1., 32, QCDTYPE::MESON);
    eta = Particle(0., 1, "eta", M_ETA, W_ETA, 0, 33, QCDTYPE::MESON);
    eta_prime = Particle(0., 1, "eta_prime", M_ETA_P, W_ETA_P, 0, 34, QCDTYPE::MESON);
    kaonp = Particle(0., 1, "kaonp", M_K_C, W_K_C, +1., 35, QCDTYPE::MESON);
    kaonm = Particle(0., 1, "kaonm", M_K_C, W_K_C, -1., 36, QCDTYPE::MESON);
    kaon0_short = Particle(0., 1, "kaon0_short", M_K_0, W_K_0_S, 0, 36, QCDTYPE::MESON);
    kaon0_long = Particle(0., 1, "kaon0_long", M_K_C, W_K_C, 0, 37, QCDTYPE::MESON);

    proton = Particle(0.5, 2, "proton", M_PROTON, 0, +1, 38, QCDTYPE::BARYON);
    anti_proton = Particle(0.5, 2, "anti_proton", M_PROTON, 0, -1, 39, QCDTYPE::BARYON);
    neutron = Particle(0.5, 2, "neutron", M_NEUTRON, 0, 0, 40, QCDTYPE::BARYON);

    n_SM_Part = n_SM_Elem + n_SM_Poly;

    // Vector containing pointers to elementary particles
    SM_Elem_Part.resize(n_SM_Elem);

    SM_Elem_Part[0] = photon;
    SM_Elem_Part[1] = gluon;
    SM_Elem_Part[2] = higgs;
    SM_Elem_Part[3] = boson_Z;
    SM_Elem_Part[4] = boson_Wp;
    SM_Elem_Part[5] = boson_Wm;
    SM_Elem_Part[6] = neutrino_e;
    SM_Elem_Part[7] = anti_neutrino_e;
    SM_Elem_Part[8] = neutrino_mu;
    SM_Elem_Part[9] = anti_neutrino_mu;
    SM_Elem_Part[10] = neutrino_tau;
    SM_Elem_Part[11] = anti_neutrino_tau;
    SM_Elem_Part[12] = electron;
    SM_Elem_Part[13] = positron;
    SM_Elem_Part[14] = muon;
    SM_Elem_Part[15] = anti_muon;
    SM_Elem_Part[16] = tau;
    SM_Elem_Part[17] = anti_tau;
    SM_Elem_Part[18] = quark_down;
    SM_Elem_Part[19] = anti_quark_down;
    SM_Elem_Part[20] = quark_up;
    SM_Elem_Part[21] = anti_quark_up;
    SM_Elem_Part[22] = quark_strange;
    SM_Elem_Part[23] = anti_quark_strange;
    SM_Elem_Part[24] = quark_charm;
    SM_Elem_Part[25] = anti_quark_charm;
    SM_Elem_Part[26] = quark_bottom;
    SM_Elem_Part[27] = anti_quark_bottom;
    SM_Elem_Part[28] = quark_top;
    SM_Elem_Part[29] = anti_quark_top;

    SM_neutrinos.resize(6);
    SM_neutrinos[0] = neutrino_e;
    SM_neutrinos[1] = anti_neutrino_e;
    SM_neutrinos[2] = neutrino_mu;
    SM_neutrinos[3] = anti_neutrino_mu;
    SM_neutrinos[4] = neutrino_tau;
    SM_neutrinos[5] = anti_neutrino_tau;

    // Vector containing vectors to all particles
    SM_Part.resize(0);
    SM_Part.insert(std::end(SM_Part), std::begin(SM_Elem_Part), std::end(SM_Elem_Part));

    SM_Part.resize(n_SM_Elem + n_SM_Poly);
    SM_Part[n_SM_Elem + 0] = pion0;
    SM_Part[n_SM_Elem + 1] = pionp;
    SM_Part[n_SM_Elem + 2] = pionm;
    SM_Part[n_SM_Elem + 3] = eta;
    SM_Part[n_SM_Elem + 4] = eta_prime;
    SM_Part[n_SM_Elem + 5] = kaonp;
    SM_Part[n_SM_Elem + 6] = kaonm;
    SM_Part[n_SM_Elem + 7] = kaon0_short;
    SM_Part[n_SM_Elem + 8] = kaon0_long;

    SM_Part[n_SM_Elem + 9] = proton;
    SM_Part[n_SM_Elem + 10] = anti_proton;
    SM_Part[n_SM_Elem + 11] = neutron;

    SM_Elem_Fermions.resize(0);

    for (int i = 0; i < SM_Elem_Part.size(); i++)
        if (!(floor(SM_Elem_Part[i].get_spin()) == SM_Elem_Part[i].get_spin()))
            SM_Elem_Fermions.push_back(SM_Elem_Part[i]);

    SM_Fermions.resize(0);
    n_SM_Fermions = 0;
    indexToSMFermIndex.resize(n_SM_Part);

    for (int i = 0; i < SM_Part.size(); i++)
        if (!(floor(SM_Part[i].get_spin()) == SM_Part[i].get_spin()))
        {
            SM_Fermions.push_back(SM_Part[i]);
            indexToSMFermIndex[i] = n_SM_Fermions;
            n_SM_Fermions++;
        }

    isPresentBeforeQCDPT.resize(n_SM_Part, false);
    isPresentAfterQCDPT.resize(n_SM_Part, false);

    for (int i = 0; i < SM_Part.size(); i++)
    {
        if (SM_Part[i].get_QCDType() != QCDTYPE::MESON && SM_Part[i].get_QCDType() != QCDTYPE::BARYON)
            isPresentBeforeQCDPT[i] = true;

        if (SM_Part[i].get_QCDType() != QCDTYPE::GLUON && SM_Part[i].get_QCDType() != QCDTYPE::QUARK)
            isPresentAfterQCDPT[i] = true;
    }

    std::cout << "--> Standard Model Particles initalised : see log file" << std::endl;
    std::clog << "## Standard Model Particles :" << std::endl;
    std::clog << SM_Part;

    // Pairs of fermion / anti_fermions
    nu_e_pair = ParticlePair(neutrino_e, anti_neutrino_e, 0);
    nu_mu_pair = ParticlePair(neutrino_mu, anti_neutrino_mu, 1);
    nu_tau_pair = ParticlePair(neutrino_tau, anti_neutrino_tau, 2);
    e_pair = ParticlePair(electron, positron, 3);
    muon_pair = ParticlePair(muon, anti_muon, 4);
    tau_pair = ParticlePair(tau, anti_tau, 5);
    q_u_pair = ParticlePair(quark_up, anti_quark_up, 6);
    q_d_pair = ParticlePair(quark_down, anti_quark_down, 7);
    q_s_pair = ParticlePair(quark_strange, anti_quark_strange, 8);
    q_c_pair = ParticlePair(quark_charm, anti_quark_charm, 9);
    q_b_pair = ParticlePair(quark_bottom, anti_quark_bottom, 10);
    q_t_pair = ParticlePair(quark_top, anti_quark_top, 11);
    //p_pair = ParticlePair(&proton, &anti_proton, 12);

    n_couples_SM_ferm = 12;
    couples_SM_ferm.resize(n_couples_SM_ferm);


    couples_SM_ferm[0] = nu_e_pair;
    couples_SM_ferm[1] = nu_mu_pair;
    couples_SM_ferm[2] = nu_tau_pair;
    couples_SM_ferm[3] = e_pair;
    couples_SM_ferm[4] = muon_pair;
    couples_SM_ferm[5] = tau_pair;
    couples_SM_ferm[6] = q_u_pair;
    couples_SM_ferm[7] = q_d_pair;
    couples_SM_ferm[8] = q_s_pair;
    couples_SM_ferm[9] = q_c_pair;
    couples_SM_ferm[10] = q_b_pair;
    couples_SM_ferm[11] = q_t_pair;
    //couples_SM_ferm[12] = p_pair;
}
