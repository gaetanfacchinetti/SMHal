#include "../headers/DarkSectorModel.h"

DarkSectorModel::DarkSectorModel()
{
    n_DS_DM = 0;
    n_DS_phis = 0;
    n_DS_phip = 0;
    n_DS_X = 0;
    n_SM_part = StandardModel::getInstance()->get_n_SM_Part();
    part_DS_index = n_SM_part;
    n_SM_couples = StandardModel::getInstance()->get_n_couples_SM_ferm();
}

// This function has to be called once all the dark sector
// particles have been declared. It makes lists of all particles
void DarkSectorModel::Initialise()
{
    d_psp.resize(n_DS_phip);

    for (int i = 0; i < n_DS_phip; i++)
    {
        d_psp[i].resize(n_DS_phis);

        for (int j = 0; j < n_DS_phis; j++)
            d_psp[i][j].resize(n_DS_phip);
    }

    for (int i = 0; i < n_DS_phip; i++)
        for (int j = 0; j < n_DS_phis; j++)
            for (int k = 0; k < n_DS_phip; k++)
                d_psp[i][j][k] = d_spp[j][i][k];

    DS_Prop.resize(0);
    DS_Prop.insert(std::end(DS_Prop), std::begin(phis), std::end(phis));
    DS_Prop.insert(std::end(DS_Prop), std::begin(phip), std::end(phip));
    DS_Prop.insert(std::end(DS_Prop), std::begin(X), std::end(X));

    DS_Elem_Part.resize(0);
    DS_Elem_Part.insert(std::end(DS_Elem_Part), std::begin(chi), std::end(chi));
    DS_Elem_Part.insert(std::end(DS_Elem_Part), std::begin(DS_Prop), std::end(DS_Prop));

    DS_Part.resize(0);
    DS_Part.insert(std::end(DS_Part), std::begin(chi), std::end(chi));
    DS_Part.insert(std::end(DS_Part), std::begin(DS_Prop), std::end(DS_Prop));

    isPresentBeforeQCDPT.resize(DS_Part.size());
    isPresentAfterQCDPT.resize(DS_Part.size());

    for (int i = 0; i < DS_Part.size(); i++)
    {
        if (DS_Part[i].get_QCDType() != QCDTYPE::MESON && DS_Part[i].get_QCDType() != QCDTYPE::BARYON)
            isPresentBeforeQCDPT[i] = true;

        if (DS_Part[i].get_QCDType() != QCDTYPE::GLUON && DS_Part[i].get_QCDType() != QCDTYPE::QUARK)
            isPresentAfterQCDPT[i] = true;
    }

    // Writting the particles of the dark sector
    std::clog << DS_Part;
    std::clog << " ---------------------------------------------------------------------------------------------------------" << std::endl;
}

void DarkSectorModel::add_darkmatter(double mass, double width, Fermiontype ftype)
{
    // By default all new DM particles are Majorana particles
    chi.push_back(Particle(0.5, 2, "chi", mass, width, 0, part_DS_index, Proptype::dm, n_DS_DM, ftype));
    part_DS_index++;
    n_DS_DM++;
}

void DarkSectorModel::add_propagator(Proptype type, int dof, std::string name, double mass, double width)
{

    double spin = 0;
    if (type == Proptype::scalar)
    {
        phis.push_back(Particle(0, dof, name, mass, width, 0, part_DS_index, type, n_DS_phis));
        n_DS_phis++;
    }
    if (type == Proptype::pseudoscalar)
    {
        phip.push_back(Particle(0, dof, name, mass, width, 0, part_DS_index, Proptype::pseudoscalar, n_DS_phip));
        n_DS_phip++;
    }
    if (type == Proptype::vector)
    {
        X.push_back(Particle(1, dof, name, mass, width, 0, part_DS_index, Proptype::vector, n_DS_X));
        n_DS_X++;
    }
    
    part_DS_index++;
}

// This function is used to resize all the coupling tables if
// we want to set them by hand afterward
void DarkSectorModel::ForceInitCouplings()
{

    //. Resizing the couplings to SM
    lambda_S_SM.resize(n_DS_phis);
    lambda_PS_SM.resize(n_DS_phip);
    a_VEC_SM.resize(n_DS_X);
    b_VEC_SM.resize(n_DS_X);

    // Resizing the couplings to DM
    lambda_S_DM.resize(n_DS_phis);
    lambda_PS_DM.resize(n_DS_phip);
    a_VEC_DM.resize(n_DS_X);
    b_VEC_DM.resize(n_DS_X);

    // Resizing the coupling between mediators
    c_sss.resize(n_DS_phis);
    d_spp.resize(n_DS_phis);
    d_psp.resize(n_DS_phip);
    g_sXX.resize(n_DS_phis);

    for (int i = 0; i < n_DS_phis; i++)
    {
        lambda_S_DM[i].resize(n_DS_DM);
        lambda_S_SM[i].resize(n_SM_couples);
        for (int j = 0; j < n_DS_DM; j++)
            lambda_S_DM[i][j].resize(n_DS_DM, 0);
    }

    for (int i = 0; i < n_DS_phip; i++)
    {
        lambda_PS_DM[i].resize(n_DS_DM);
        lambda_PS_SM[i].resize(n_SM_couples);
        for (int j = 0; j < n_DS_DM; j++)
            lambda_PS_DM[i][j].resize(n_DS_DM, 0);

        //std::cout << "coucou" << " " << lambda_PS_DM.size() << " " << lambda_PS_DM[0].size() << " " << lambda_PS_DM[0][0].size() << std::endl;
    }

    for (int i = 0; i < n_DS_X; i++)
    {
        a_VEC_DM[i].resize(n_DS_DM);
        b_VEC_DM[i].resize(n_DS_DM);
        a_VEC_SM[i].resize(n_SM_couples);
        b_VEC_SM[i].resize(n_SM_couples);

        for (int j = 0; j < n_DS_DM; j++)
        {
            a_VEC_DM[i][j].resize(n_DS_DM, 0);
            b_VEC_DM[i][j].resize(n_DS_DM, 0);
        }
    }


    for(int i = 0; i < n_DS_phis; i++)
    {
        c_sss[i].resize(n_DS_phis);
        for(int j = 0; j < n_DS_phis; j++)
            c_sss[i][j].resize(n_DS_phis, 0);
    }

    for(int i = 0; i < n_DS_phis; i++)
    {
        d_spp[i].resize(n_DS_phip);
        for(int j = 0; j < n_DS_phip; j++)
            d_spp[i][j].resize(n_DS_phip, 0);
    }

    for(int i = 0; i < n_DS_phip; i++)
    {
        d_psp[i].resize(n_DS_phis);
        for(int j = 0; j < n_DS_phis; j++)
            d_psp[i][j].resize(n_DS_phip, 0);
    }

    for(int i = 0; i < n_DS_phis; i++)
    {
        g_sXX[i].resize(n_DS_X);
        for(int j = 0; j < n_DS_X; j++)
            g_sXX[i][j].resize(n_DS_X, 0);
    }




}
