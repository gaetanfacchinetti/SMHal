#ifndef STANDARDMODEL
#define STANDARDMODEL

#include "CppLib.h"
#include "Particle.h"

class StandardModel
{

  public:
    // Fonctions qui en font un singleton
    static StandardModel *getInstance()
    {
        if (NULL == _SM)
            _SM = new StandardModel;

        return _SM;
    }

    static void kill()
    {
        if (NULL != _SM)
        {
            delete _SM;
            _SM = NULL;
        }
    }

    std::vector<int> get_indexToSMFermIndex() const {return indexToSMFermIndex;};
    Particle get_SM_Part(int i) const {return SM_Part[i];};
    std::vector<Particle> get_SM_Part() const  {return SM_Part;};
    ParticlePair get_couples_SM_ferm(int i) const  {return couples_SM_ferm[i];};
    std::vector<ParticlePair> get_couples_SM_ferm() const  {return couples_SM_ferm;};
    std::vector<Particle> get_SM_neutrinos() const {return SM_neutrinos;};

    int get_n_SM_Part() const {return n_SM_Part;};
    int get_n_couples_SM_ferm() const {return n_couples_SM_ferm;};

    Particle get_electron() const {return electron;};
    Particle * get_photon_ptr() {return &photon;};
    // Particle get_neutrino() const {return neutrino;};


  private:
    StandardModel();
    ~StandardModel(){};

    static StandardModel *_SM;

    Particle photon, gluon, higgs, boson_Z, boson_Wp, boson_Wm;
    Particle neutrino_e, neutrino_mu, neutrino_tau;
    Particle anti_neutrino_e, anti_neutrino_mu, anti_neutrino_tau;
    Particle electron, positron, muon, anti_muon, tau, anti_tau;
    Particle quark_up, quark_down, quark_strange, quark_charm, quark_bottom, quark_top;
    Particle anti_quark_up, anti_quark_down, anti_quark_strange, anti_quark_charm, anti_quark_bottom, anti_quark_top;
    Particle pion0, pionp, pionm, eta, eta_prime, kaon0_short, kaon0_long, kaonp, kaonm, proton, anti_proton, neutron;

    ParticlePair e_pair, muon_pair, tau_pair, nu_e_pair, nu_mu_pair, nu_tau_pair, q_u_pair, q_d_pair, q_s_pair, q_c_pair, q_b_pair, q_t_pair, p_pair;
    std::vector<ParticlePair> couples_SM_ferm;
    int n_couples_SM_ferm;

    int n_SM_Part, n_SM_Elem, n_SM_Poly;

    std::vector<Particle> SM_Elem_Part, SM_Elem_Fermions, SM_Part, SM_Fermions;
    std::vector<Particle> SM_neutrinos;
    std::vector<bool> isPresentBeforeQCDPT, isPresentAfterQCDPT;

    std::vector<int> indexToSMFermIndex;
    int n_SM_Fermions;
};

#endif