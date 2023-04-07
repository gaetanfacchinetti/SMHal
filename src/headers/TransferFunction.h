#ifndef DEF_TRANSFERFUNCTION
#define DEF_TRANSFERFUNCTION

#include "CppLib.h"
#include "MyUnits.h"
#include "Cosmology.h"

// Content of the class provided by M. Stref

class TransferFunction
{
    public : 
        TransferFunction(){_cosmology = Cosmology();};
        TransferFunction(Cosmology cosmology) : _cosmology(cosmology){};
        virtual ~TransferFunction(){};

        Cosmology get_cosmo() const {return _cosmology;};
        void plot();

        void test(){std::cout << "test of the function" << std::endl;};

        /// Pure virtual function making the class abstract
        /// k should be in Mpc^{-1}
        virtual double evaluate(double k) {return 1.;};

    protected:
        Cosmology _cosmology;

};


// Transfer function class from Green Hofmann and Schwarz (2005)
class TransferFunction_WIMP_GHS05 : public TransferFunction
{

    public:
        /// x=Tkd and y=mchi in GeV if Tm = True
        /// x=kfs and y = kd in Mpc if Tm = False 
        TransferFunction_WIMP_GHS05(Cosmology n_cosmo, double x, double y, bool Tm = false);   
        TransferFunction_WIMP_GHS05() : TransferFunction(){};
        virtual ~TransferFunction_WIMP_GHS05(){};

        static std::shared_ptr<TransferFunction_WIMP_GHS05> TransferFunction_WIMP_GHS05_from_Tempmchi(Cosmology cosmology, double Tkd, double mchi)
            {return std::make_shared<TransferFunction_WIMP_GHS05>(cosmology, Tkd, mchi, true);};

        static std::shared_ptr<TransferFunction_WIMP_GHS05> TransferFunction_WIMP_GHS05_from_kfskd(Cosmology cosmology, double kfs, double kd)
            {return std::make_shared<TransferFunction_WIMP_GHS05>(cosmology, kfs, kd, false);};

        virtual double evaluate(double k); 

    private: 
        double kfs, kd; // Mpc^{-1}
};



// Transfer function class from Eisenstein and Hu 1998
class TransferFunction_EH98 : public TransferFunction
{
    
public:

    TransferFunction_EH98() : TransferFunction(){};
    TransferFunction_EH98(Cosmology n_cosmo, bool with_baryons = true);
    virtual ~TransferFunction_EH98(){};
    

    double R_ratio(double z);
    double T0_tilde(double q, double alphac, double betac);
    double shape_parameter(double k);
    double T_cdm(double k);
    double G_func(double y);
    double s_tilde(double k);
    double T_baryons(double k);

    /**
     * \brief Transfer function
     * \param (double) Comoving scale \f$k~[\rm Mpc^{-1}]\f$ 
     * \return (double) Transfer function [dimensionless] */
    virtual double evaluate(double k);
    double effective_CDM_transfer_function(double k);
    

    
private:
    bool _with_baryons;
    double Omega_c_h2, Omega_b_h2, Omega_m_h2;
    double z_eq, k_eq, z_drag, sound_horizon, alpha_c, beta_c, k_silk, alpha_b, beta_b;
    double Theta27 = T0_CMB_K/2.7; // Normalised CMB temperature
    
};



#endif
