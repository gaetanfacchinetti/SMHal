#include "../headers/InputReader.h"

InputReader::InputReader()
{
    mass_DM.resize(0);
    //mass_DM[0].resize(0);

    mass_prop_scalar.resize(0);
    //mass_prop_scalar[0].resize(0);

    mass_prop_pseudoscalar.resize(0);
    //mass_prop_pseudoscalar[0].resize(0);

    mass_prop_vector.resize(0);

    lambda_S_SM.resize(0);
    lambda_PS_SM.resize(0);
    a_VEC_SM.resize(0);
    b_VEC_SM.resize(0);

    lambda_S_DM.resize(0);
    lambda_PS_DM.resize(0);
    a_VEC_DM.resize(0);
    b_VEC_DM.resize(0);

    c_sss.resize(0);
    d_spp.resize(0);
    g_sXX.resize(0);

    prop_type.resize(0);
    glob_ind_to_scalar_ind.resize(0);
    glob_ind_to_pseudoscalar_ind.resize(0);
    glob_ind_to_vector_ind.resize(0);
}

//
//
//
//

std::vector<DarkSectorModel> InputReader::Read(std::string name_file)
{
    std::string line;
    int n_line = 0;
    int n_line_written = 0;

    dof_prop_scalar.resize(0, 0);
    dof_prop_pseudoscalar.resize(0, 0);
    dof_prop_vector.resize(0, 0);

    glob_ind_to_scalar_ind.resize(0);
    glob_ind_to_pseudoscalar_ind.resize(0);
    glob_ind_to_vector_ind.resize(0);

    prop_type.resize(0, Proptype::scalar);
    ftype_DM.resize(0, Fermiontype::majorana); // By default Majorana particles
    n_dm_part = 0;

    std::string delimiter_0 = "#";
    std::string delimiter_1 = "::";

    std::string token;
    size_t pos0, pos1;

    n_prop_scalar = 0;
    n_prop_pseudoscalar = 0;
    n_prop_vector = 0;

    std::vector<DarkSectorModel> _DS_vec;

    _DS_vec.resize(0);
    mass_DM.resize(0);
    mass_prop_scalar.resize(0);
    mass_prop_pseudoscalar.resize(0);
    mass_prop_vector.resize(0);

    lambda_S_SM.resize(0);
    lambda_PS_SM.resize(0);
    a_VEC_SM.resize(0);
    b_VEC_SM.resize(0);

    lambda_S_DM.resize(0);
    lambda_PS_DM.resize(0);
    a_VEC_DM.resize(0);
    b_VEC_DM.resize(0);

    c_sss.resize(0);
    d_spp.resize(0);
    g_sXX.resize(0);

    std::fstream infile(name_file);
    if (!infile.is_open())
        std::cout << "WARNING :: File " << name_file << " not found in " << __PRETTY_FUNCTION__ << std::endl;
    else
    {
        // Loop on all the lines
        while (getline(infile, line))
        {
            pos1 = line.find(delimiter_1);

            if (pos1 != std::string::npos)
            {
                token = line.substr(0, pos1);

                if (trim(token) == "Npoints")
                {
                    line.erase(0, pos1 + delimiter_1.length());
                    Npoints = stoi(trim(line));
                    _DS_vec.resize(Npoints);
                    mass_DM.resize(Npoints);
                    mass_prop_scalar.resize(Npoints);
                    mass_prop_pseudoscalar.resize(Npoints);
                    mass_prop_vector.resize(Npoints);

                    // Initialise all models to majorana fermions
                    ftype_DM.resize(Npoints, Fermiontype::majorana);

                    lambda_S_SM.resize(Npoints);
                    lambda_PS_SM.resize(Npoints);
                    a_VEC_SM.resize(Npoints);
                    b_VEC_SM.resize(Npoints);

                    lambda_S_DM.resize(Npoints);
                    lambda_PS_DM.resize(Npoints);
                    a_VEC_DM.resize(Npoints);
                    b_VEC_DM.resize(Npoints);

                    c_sss.resize(Npoints);
                    d_spp.resize(Npoints);
                    g_sXX.resize(Npoints);
                }
                if (trim(token) == "DSCONT")
                {
                    ReadDSContent(line, n_line);
                }
            }

            n_line++;
        }

        // We reread the file a second time (for the DS content and masses)
        n_line = 0;

        infile.clear();
        infile.seekg(0, std::ios::beg);

        if (infile.is_open())
        {
            while (getline(infile, line))
            {
                pos1 = line.find(delimiter_1);

                if (pos1 != std::string::npos)
                {
                    token = line.substr(0, pos1);

                    //std::cout << token << std::endl;

                    if (trim(token) == "TYPEDM")
                        ReadTypeDM(line, n_line);
                    if (trim(token) == "MDMDS")
                        ReadMassDMDS(line, n_line);
                }
                //else
                //    std::cout << line << std::endl;

                n_line++;
            }
        }

        //std::cout << "Number particles " << n_dm_part << std::endl;

        // Out of the loop reading the file
        for (int i = 0; i < Npoints; i++)
        {
            lambda_S_SM[i].resize(n_prop_scalar);
            lambda_PS_SM[i].resize(n_prop_pseudoscalar);
            a_VEC_SM[i].resize(n_prop_vector);
            b_VEC_SM[i].resize(n_prop_vector);

            for (int j = 0; j < n_prop_scalar; j++)
                lambda_S_SM[i][j].resize(20, 0);

            for (int j = 0; j < n_prop_pseudoscalar; j++)
                lambda_PS_SM[i][j].resize(20, 0);

            for (int j = 0; j < n_prop_vector; j++)
            {
                a_VEC_SM[i][j].resize(20, 0);
                b_VEC_SM[i][j].resize(20, 0);
            }
        }

        for (int i = 0; i < Npoints; i++)
        {
            lambda_S_DM[i].resize(n_prop_scalar);
            lambda_PS_DM[i].resize(n_prop_pseudoscalar);
            a_VEC_DM[i].resize(n_prop_vector);
            b_VEC_DM[i].resize(n_prop_vector);

            c_sss[i].resize(n_prop_scalar);
            d_spp[i].resize(n_prop_scalar);
            g_sXX[i].resize(n_prop_scalar);

            for (int j = 0; j < n_prop_scalar; j++)
            {
                lambda_S_DM[i][j].resize(n_dm_part);
                c_sss[i][j].resize(n_prop_scalar);
                d_spp[i][j].resize(n_prop_pseudoscalar);
                g_sXX[i][j].resize(n_prop_vector);

                for (int k = 0; k < n_dm_part; k++)
                    lambda_S_DM[i][j][k].resize(n_dm_part, 0);

                for (int k = 0; k < n_prop_scalar; k++)
                    c_sss[i][j][k].resize(n_prop_scalar, 0);

                for (int k = 0; k < n_prop_pseudoscalar; k++)
                    d_spp[i][j][k].resize(n_prop_pseudoscalar, 0);

                for (int k = 0; k < n_prop_vector; k++)
                    g_sXX[i][j][k].resize(n_prop_vector, 0);
            }

            for (int j = 0; j < n_prop_pseudoscalar; j++)
            {
                lambda_PS_DM[i][j].resize(n_dm_part);
                for (int k = 0; k < n_dm_part; k++)
                    lambda_PS_DM[i][j][k].resize(n_dm_part, 0);
            }

            for (int j = 0; j < n_prop_vector; j++)
            {
                a_VEC_DM[i][j].resize(n_dm_part);
                b_VEC_DM[i][j].resize(n_dm_part);

                for (int k = 0; k < n_dm_part; k++)
                {
                    a_VEC_DM[i][j][k].resize(n_dm_part, 0);
                    b_VEC_DM[i][j][k].resize(n_dm_part, 0);
                }
            }
        }
    }

    //std::cout << lambda_PS_DM[0][0][1][0] << std::endl;

    // We read the file a third time (for the interactions)
    n_line = 0;

    infile.clear();
    infile.seekg(0, std::ios::beg);

    if (infile.is_open())
    {
        while (getline(infile, line))
        {
            pos1 = line.find(delimiter_1);

            if (pos1 != std::string::npos)
            {
                token = line.substr(0, pos1);

                //std::cout << token << std::endl;
                if (trim(token) == "IntDSSM")
                    ReadIntDSSM(line, n_line);
                if (trim(token) == "IntDSDM")
                    ReadIntDSDM(line, n_line);
                if (trim(token) == "IntDSDS")
                    ReadIntDSDS(line, n_line);
            }
            //else
            //    std::cout << line << std::endl;

            n_line++;
        }
    }

    infile.close();

    //std::cout << "jdozij 4 : " << dof_prop_scalar[0] << " " << dof_prop_scalar[1] << " " << dof_prop_scalar[2] << std::endl;
    //std::cout << "Expected " << Npoints << " points in this file" << std::endl;

    // Creation of the particles the table of DarkSectorModels
    for (int i = 0; i < Npoints; i++)
    {
        _DS_vec[i] = DarkSectorModel();

        for (int j = 0; j < n_dm_part; j++)
            _DS_vec[i].add_darkmatter(mass_DM[i][j], 0, ftype_DM[i]);

        for (int j = 0; j < n_prop_scalar; j++)
            _DS_vec[i].add_propagator(Proptype::scalar, dof_prop_scalar[j], "phis", mass_prop_scalar[i][j], 0);

        for (int j = 0; j < n_prop_pseudoscalar; j++)
            _DS_vec[i].add_propagator(Proptype::pseudoscalar, dof_prop_pseudoscalar[j], "phip", mass_prop_pseudoscalar[i][j], 0);

        for (int j = 0; j < n_prop_vector; j++)
            _DS_vec[i].add_propagator(Proptype::vector, dof_prop_vector[j], "X", mass_prop_vector[i][j], 0);

        for (int j = 0; j < n_prop_scalar; j++)
        {
            Symetrise2(lambda_S_DM[i][j]);
            Symetrise2(d_spp[i][j]);
            Symetrise2(g_sXX[i][j]);
        }

        Symetrise3(c_sss[i]);

        for (int j = 0; j < n_prop_pseudoscalar; j++)
            Symetrise2(lambda_PS_DM[i][j]);

        for (int j = 0; j < n_prop_vector; j++)
        {
            Symetrise2(a_VEC_DM[i][j]);
            Symetrise2(b_VEC_DM[i][j]);
        }

        /* try{
            std::cout << "A la lecture : " << b_VEC_SM[0][0][3] << std::endl;}
            catch(const std::exception & e){std::cerr << e.what(); }; */

        _DS_vec[i].set_couplings_DS_S_SM(lambda_S_SM[i]);
        _DS_vec[i].set_couplings_DS_PS_SM(lambda_PS_SM[i]);
        _DS_vec[i].set_couplings_DS_VEC_SM(a_VEC_SM[i], b_VEC_SM[i]);

        _DS_vec[i].set_couplings_DS_S_DM(lambda_S_DM[i]);
        _DS_vec[i].set_couplings_DS_PS_DM(lambda_PS_DM[i]);
        _DS_vec[i].set_couplings_DS_VEC_DM(a_VEC_DM[i], b_VEC_DM[i]);

        _DS_vec[i].set_couplings_c_sss(c_sss[i]);
        _DS_vec[i].set_couplings_d_spp(d_spp[i]);
        _DS_vec[i].set_couplings_g_sXX(g_sXX[i]);

        std::clog << "## Dark sector particles : model " << i << std::endl;
        _DS_vec[i].Initialise();

        //std::cout << "here " << lambda_PS_DM[0][0][0][0] << std::endl;

        /*std::cout << "lambda_S_SM" << std::endl;
        for(int i = 0; i<n_prop_scalar; i++)
            for(int j = 0; j<12; j++)
             {
                 std::cout << i << " " << j << " " << lambda_S_SM[0][i][j] << std::endl; 
            } */
    }

    if (Verbose == true)
    {
        std::cout << "--> DarkSector Model initalised :  see log file" << std::endl;
    }
    return _DS_vec;
}

void InputReader::ReadDSContent(std::string &content, int n_line)
{

    content = trim(content);

    size_t pos;
    std::string token1, token2, token3, token4, token5;
    std::string delimiter_1 = "::";
    std::string delimiter_2 = ",";

    std::string particle_prop;

    try
    {
        pos = content.find(delimiter_1);
        token1 = content.substr(0, pos);
        content.erase(0, pos + delimiter_1.length());

        pos = content.find(delimiter_1);
        token2 = content.substr(0, pos);
        content.erase(0, pos + delimiter_1.length());
        token3 = content;
    }
    catch (std::exception const &e)
    {
        std::cerr << "SYNTAX ERROR :: line " << n_line << " -> " << e.what() << std::endl;
    }

    if (trim(token2) == "NDM")
    {
        //std::cout << "ici " << stoi(trim(token3)) << std::endl;
        n_dm_part = stoi(trim(token3));
    }

    if (trim(token2) == "DSPart")
        particle_prop = trim(token3);

    //std::cout << particle_prop << std::endl;

    pos = particle_prop.find(delimiter_2);
    token4 = particle_prop.substr(0, pos);
    particle_prop.erase(0, pos + delimiter_2.length());
    token5 = particle_prop;

    //std::cout << ltrim(token4, "(") << std::endl;

    if (ltrim(token4, "(") == "scalar")
    {

        dof_prop_scalar.push_back(stoi(rtrim(token5, ")")));
        prop_type.push_back(Proptype::scalar);
        glob_ind_to_scalar_ind.push_back(n_prop_scalar);
        glob_ind_to_pseudoscalar_ind.push_back(-1);
        glob_ind_to_vector_ind.push_back(-1);
        n_prop_scalar++;
        //std::cout << "dof : " << dof_prop_scalar[dof_prop_scalar.size() - 1] << std::endl;
        //if(prop_type[prop_type.size()-1] == Proptype::scalar)
        //    std::cout << static_cast<std::underlying_type<Proptype>::type>(prop_type[prop_type.size()-1]) << std::endl;
    }
    if (ltrim(token4, "(") == "pseudoscalar")
    {

        dof_prop_pseudoscalar.push_back(stoi(rtrim(token5, ")")));
        prop_type.push_back(Proptype::pseudoscalar);
        glob_ind_to_scalar_ind.push_back(-1);
        glob_ind_to_pseudoscalar_ind.push_back(n_prop_pseudoscalar);
        glob_ind_to_vector_ind.push_back(-1);
        n_prop_pseudoscalar++;
    }
    if (ltrim(token4, "(") == "vector")
    {
        dof_prop_vector.push_back(stoi(rtrim(token5, ")")));
        prop_type.push_back(Proptype::vector);
        glob_ind_to_scalar_ind.push_back(-1);
        glob_ind_to_pseudoscalar_ind.push_back(-1);
        glob_ind_to_vector_ind.push_back(n_prop_vector);
        n_prop_vector++;
    }
}

//
//

void InputReader::ReadTypeDM(std::string &content, int n_line)
{
    content = trim(content);

    size_t pos;
    std::string token1, token2, token3, token4, token5, token6, token7;
    std::string delimiter_1 = "::";
    std::string delimiter_2 = ",";
    std::string delimiter_3 = " ";
    int index;

    try
    {
        pos = content.find(delimiter_1);
        token1 = content.substr(0, pos);
        content.erase(0, pos + delimiter_1.length());

        pos = content.find(delimiter_1);
        token2 = content.substr(0, pos);
        content.erase(0, pos + delimiter_1.length());
        token3 = content;
    }
    catch (std::exception const &e)
    {
        std::cerr << "SYNTAX ERROR :: line " << n_line << " -> " << e.what() << std::endl;
    }

    try
    {
        //std::cout << n_line << " " << trim(token2) << std::endl;
        index = stoi(trim(token2));
    }
    catch (std::exception const &e)
    {
        std::cerr << "SYNTAX ERROR :: line " << n_line << " -> " << e.what() << std::endl;
    }

    token3 = trim(token3);

    //std::cout << token3 << std::endl;

    if (index > Npoints - 1)
    {
        std::cerr << "Number of declared models (Npoints) is too low in comparison of the content of the file at line :: " << n_line << std::endl;
        exit(0);
    }

    if (token3 == "Majorana")
    {
        ftype_DM[index] = Fermiontype::majorana;
        return;
    }
    if (token3 == "Dirac")
    {
        ftype_DM[index] = Fermiontype::dirac;
        n_dm_part = 2 * n_dm_part; // We distinguish between the particle and its anti-particle
        return;
    }

    std::cerr << "SYNTAX ERROR :: line " << n_line << " -> Please enter a valid Fermion type (Majorana/Dirac)" << std::endl;
}

void InputReader::ReadMassDMDS(std::string &content, int n_line)
{
    content = trim(content);

    size_t pos;
    std::string token1, token2, token3, token4, token5, token6, token7;
    std::string delimiter_1 = "::";
    std::string delimiter_2 = ",";
    std::string delimiter_3 = " ";
    int index;

    try
    {
        pos = content.find(delimiter_1);
        token1 = content.substr(0, pos);
        content.erase(0, pos + delimiter_1.length());

        pos = content.find(delimiter_1);
        token2 = content.substr(0, pos);
        content.erase(0, pos + delimiter_1.length());
        token3 = content;
    }
    catch (std::exception const &e)
    {
        std::cerr << "SYNTAX ERROR :: line " << n_line << " -> " << e.what() << std::endl;
    }

    try
    {
        //std::cout << n_line << " " << trim(token2) << std::endl;
        index = stoi(trim(token2));
    }
    catch (std::exception const &e)
    {
        std::cerr << "SYNTAX ERROR :: line " << n_line << " -> " << e.what() << std::endl;
    }

    token3 = trim(token3);

    pos = token3.find(delimiter_2);
    token4 = token3.substr(0, pos);
    token3.erase(0, pos + delimiter_2.length());
    token5 = token3;

    //std::cout << "We have:" << token3 << "|" << token4 << "|" << token5 << std::endl;

    if (index > Npoints - 1)
    {
        std::cerr << "Number of declared models (Npoints) is too low in comparison of the content of the file at line :: " << n_line << std::endl;
        exit(0);
    }

    mass_DM[index].resize(0);
    //std::cout << "Number of DM particles : " << n_dm_part << std::endl;

    int n_dm_species = n_dm_part;
    if (ftype_DM[index] == Fermiontype::dirac)
        n_dm_species = (int)(n_dm_part / 2.);

    for (int i = 0; i < n_dm_species; i++)
    {
        if (i < n_dm_species - 1)
        {
            pos = token4.find(delimiter_3);
            token6 = token4.substr(0, pos);
            mass_DM[index].push_back(stod(trim(token6)));
            token4.erase(0, pos + delimiter_3.length());
        }
        else
        {
            //std::cout << "here" << std::endl;
            token6 = token4;
            mass_DM[index].push_back(stod(trim(token6)));
            try
            {
                token4.erase(0, pos + delimiter_3.length());
            }
            catch (std::exception const &e)
            {
                std::cerr << "SYNTAX ERROR :: line " << n_line << " -> " << e.what() << std::endl;
            }
        }

        if (ftype_DM[index] == Fermiontype::dirac)
            mass_DM[index].push_back(mass_DM[index][mass_DM[index].size() - 1]);
    }

    // Here we check that we correclty initialised all DM particles
    if (mass_DM[index].size() != n_dm_part || trim(token4) != "")
        std::cout << "WARNING : Number of DM mass given do not correspond to the number of DM mass introduced" << std::endl;

    // We order the mass index to place the smallest first
    order1(mass_DM[index]);

    // Initialise the vectors containing the mass of the propagators
    for (int i = 0; i < Npoints; i++)
    {
        mass_prop_scalar[i].resize(n_prop_scalar);
        mass_prop_pseudoscalar[i].resize(n_prop_pseudoscalar);
        mass_prop_vector[i].resize(n_prop_vector);
    }

    int i_scalar = 0, i_pseudoscalar = 0, i_vector = 0;

    token5 = trim(token5);
    //std::cout << token5 << std::endl;

    for (int i = 0; i < n_prop_scalar + n_prop_pseudoscalar + n_prop_vector; i++)
    {
        if (i < n_prop_scalar + n_prop_pseudoscalar + n_prop_vector - 1)
        {
            pos = token5.find(delimiter_3);
            token7 = token5.substr(0, pos);
        }
        else
            token7 = token5;

        //std::cout << "ici " << i << " " << token7 << std::endl;
        //std::cout << "pff " << static_cast<std::underlying_type<Proptype>::type>(prop_type[i]) << std::endl;

        if (prop_type[i] == Proptype::scalar)
        {
            mass_prop_scalar[index][i_scalar] = stod(trim(token7));
            i_scalar++;
        }
        if (prop_type[i] == Proptype::pseudoscalar)
        {
            mass_prop_pseudoscalar[index][i_pseudoscalar] = stod(trim(token7));
            i_pseudoscalar++;
        }
        if (prop_type[i] == Proptype::vector)
        {
            mass_prop_vector[index][i_vector] = stod(trim(token7));
            i_vector++;
        }

        if (i < n_prop_scalar + n_prop_pseudoscalar + n_prop_vector - 1)
            token5.erase(0, pos + delimiter_3.length());
    }
}

void InputReader::ReadIntDSSM(std::string &interactions, int n_line)
{
    interactions = trim(interactions);

    size_t pos, pos2;
    std::string token1, token2, token3, token4, token5, token6, token7;
    std::string delimiter_1 = "::";
    std::string delimiter_2 = ",";
    std::string delimiter_3 = " ";
    std::string delimiter_4 = "=";
    int index;

    try
    {
        pos = interactions.find(delimiter_1);
        token1 = interactions.substr(0, pos);
        interactions.erase(0, pos + delimiter_1.length());

        pos = interactions.find(delimiter_1);
        token2 = interactions.substr(0, pos);
        interactions.erase(0, pos + delimiter_1.length());

        pos = interactions.find(delimiter_1);
        token3 = interactions.substr(0, pos);
        interactions.erase(0, pos + delimiter_1.length());

        token4 = interactions;
    }
    catch (std::exception const &e)
    {
        std::cerr << "SYNTAX ERROR :: line " << n_line << " -> " << e.what() << std::endl;
    }

    int id_prop = 0;

    try
    {

        index = stoi(trim(token2));
        token3 = trim(token3);
        id_prop = stoi(token3);

        //std::cout << n_line << " " << trim(token2)  << " " << id_prop << std::endl;
    }
    catch (std::exception const &e)
    {
        std::cerr << "SYNTAX ERROR :: line " << n_line << " -> " << e.what() << std::endl;
    }

    double val = 0;

    token4 = trim(token4);

    bool end_of_line_reached = false;

    do
    {
        pos = token4.find(delimiter_2);

        if (pos == std::string::npos)
        {
            token5 = token4;
            end_of_line_reached = true;
        }
        else
        {
            //std::cout << "token4 = " << token4 << " " << pos << std::endl;
            //std::cout << delimiter_2.length() << std::endl;
            //std::cout << trim(token5) << std::endl;

            token5 = token4.substr(0, pos);
            token4.erase(0, pos + delimiter_2.length());
            token4 = trim(token4);
        }

        //std::cout << "token4 = " << token4 << std::endl;
        //std::cout << "token5 = " << token5 << std::endl;

        pos2 = token5.find(delimiter_4);
        token6 = token5.substr(0, pos2);
        token5.erase(0, pos2 + delimiter_4.length());
        token7 = token5;

        token6 = trim(token6);
        token7 = trim(token7);

        try
        {
            //std::cout << "oui : " << token7 << " " << token6 << " | " << token4 << " | " << std::endl;
            val = stod(token7);
        }
        catch (std::exception const &e)
        {
            std::cerr << "SYNTAX ERROR :: line " << token7 << " " << n_line << " -> " << e.what() << std::endl;
        }

        //std::cout << stod(token7) << std::endl;

        //token4.erase(0, pos + delimiter_2.length());

        if (index > Npoints - 1)
        {
            std::cerr << "Number of declared models (Npoints) is too low in comparison of the content of the file at line :: " << n_line << std::endl;
            exit(0);
        }

        if (prop_type[id_prop] == Proptype::scalar)
        {
            //std::cout << token6 << " " << val << " " << index << " " << num_s << std::endl;
            CompleteCouplingTables_lambda_S_SM(token6, val, index, glob_ind_to_scalar_ind[id_prop]);
        }

        if (prop_type[id_prop] == Proptype::pseudoscalar)
            CompleteCouplingTables_lambda_PS_SM(token6, val, index, glob_ind_to_pseudoscalar_ind[id_prop]);

        if (prop_type[id_prop] == Proptype::vector)
            CompleteCouplingTables_ab_VEC_SM(token6, val, index, glob_ind_to_vector_ind[id_prop]);

        //std::cout << token6 << " " << token7 << std::endl;
    } while (!end_of_line_reached);
}

void InputReader::CompleteCouplingTables_lambda_S_SM(std::string const &name, double val, int id, int id_part)
{

    //std::cout << name << " " << id_part << " " << val << " " << std::endl;

    //std::cout << lambda_S_SM.size() << std::endl;
    //exit(0);

    if (name == "l_nu_e")
        lambda_S_SM[id][id_part][0] = val;
    if (name == "l_nu_mu")
        lambda_S_SM[id][id_part][1] = val;
    if (name == "l_nu_tau")
        lambda_S_SM[id][id_part][2] = val;
    if (name == "l_e")
        lambda_S_SM[id][id_part][3] = val;
    if (name == "l_mu")
        lambda_S_SM[id][id_part][4] = val;
    if (name == "l_tau")
        lambda_S_SM[id][id_part][5] = val;
    if (name == "l_u")
        lambda_S_SM[id][id_part][6] = val;
    if (name == "l_d")
        lambda_S_SM[id][id_part][7] = val;
    if (name == "l_s")
        lambda_S_SM[id][id_part][8] = val;
    if (name == "l_c")
        lambda_S_SM[id][id_part][9] = val;
    if (name == "l_b")
        lambda_S_SM[id][id_part][10] = val;
    if (name == "l_t")
        lambda_S_SM[id][id_part][11] = val;

    if (name == "l_ga")
        lambda_S_SM[id][id_part][12] = val;
}

void InputReader::CompleteCouplingTables_lambda_PS_SM(std::string const &name, double val, int id, int id_part)
{
    if (name == "l_nu_e")
        lambda_PS_SM[id][id_part][0] = val;
    if (name == "l_nu_mu")
        lambda_PS_SM[id][id_part][1] = val;
    if (name == "l_nu_tau")
        lambda_PS_SM[id][id_part][2] = val;
    if (name == "l_e")
        lambda_PS_SM[id][id_part][3] = val;
    if (name == "l_mu")
        lambda_PS_SM[id][id_part][4] = val;
    if (name == "l_tau")
        lambda_PS_SM[id][id_part][5] = val;
    if (name == "l_u")
        lambda_PS_SM[id][id_part][6] = val;
    if (name == "l_d")
        lambda_PS_SM[id][id_part][7] = val;
    if (name == "l_s")
        lambda_PS_SM[id][id_part][8] = val;
    if (name == "l_c")
        lambda_PS_SM[id][id_part][9] = val;
    if (name == "l_b")
        lambda_PS_SM[id][id_part][10] = val;
    if (name == "l_t")
        lambda_PS_SM[id][id_part][11] = val;
    if (name == "l_ga")
        lambda_PS_SM[id][id_part][12] = val;
}

void InputReader::CompleteCouplingTables_ab_VEC_SM(std::string const &name, double val, int id, int id_part)
{
    if (name == "a_nu_e")
        a_VEC_SM[id][id_part][0] = val;
    if (name == "a_nu_mu")
        a_VEC_SM[id][id_part][1] = val;
    if (name == "a_nu_tau")
        a_VEC_SM[id][id_part][2] = val;
    if (name == "a_e")
        a_VEC_SM[id][id_part][3] = val;
    if (name == "a_mu")
        a_VEC_SM[id][id_part][4] = val;
    if (name == "a_tau")
        a_VEC_SM[id][id_part][5] = val;
    if (name == "a_u")
        a_VEC_SM[id][id_part][6] = val;
    if (name == "a_d")
        a_VEC_SM[id][id_part][7] = val;
    if (name == "a_s")
        a_VEC_SM[id][id_part][8] = val;
    if (name == "a_c")
        a_VEC_SM[id][id_part][9] = val;
    if (name == "a_b")
        a_VEC_SM[id][id_part][10] = val;
    if (name == "a_t")
        a_VEC_SM[id][id_part][11] = val;

    if (name == "b_nu_e")
        b_VEC_SM[id][id_part][0] = val;
    if (name == "b_nu_mu")
        b_VEC_SM[id][id_part][1] = val;
    if (name == "b_nu_tau")
        b_VEC_SM[id][id_part][2] = val;
    if (name == "b_e")
        b_VEC_SM[id][id_part][3] = val;
    if (name == "b_mu")
        b_VEC_SM[id][id_part][4] = val;
    if (name == "b_tau")
        b_VEC_SM[id][id_part][5] = val;
    if (name == "b_u")
        b_VEC_SM[id][id_part][6] = val;
    if (name == "b_d")
        b_VEC_SM[id][id_part][7] = val;
    if (name == "b_s")
        b_VEC_SM[id][id_part][8] = val;
    if (name == "b_c")
        b_VEC_SM[id][id_part][9] = val;
    if (name == "b_b")
        b_VEC_SM[id][id_part][10] = val;
    if (name == "b_t")
        b_VEC_SM[id][id_part][11] = val;
}

//
//

void InputReader::ReadIntDSDM(std::string &interactions, int n_line)
{

    interactions = trim(interactions);

    size_t pos, pos2;
    std::string token1, token2, token3, token4, token5, token6, token7;
    std::string delimiter_1 = "::";
    std::string delimiter_2 = ",";
    std::string delimiter_3 = " ";
    std::string delimiter_4 = "=";
    std::string delimiter_5 = "_";
    int index;

    try
    {
        pos = interactions.find(delimiter_1);
        token1 = interactions.substr(0, pos);
        interactions.erase(0, pos + delimiter_1.length());

        pos = interactions.find(delimiter_1);
        token2 = interactions.substr(0, pos);
        interactions.erase(0, pos + delimiter_1.length());

        pos = interactions.find(delimiter_1);
        token3 = interactions.substr(0, pos);
        interactions.erase(0, pos + delimiter_1.length());

        token4 = interactions;
    }
    catch (std::exception const &e)
    {
        std::cerr << "SYNTAX ERROR :: line " << n_line << " -> " << e.what() << std::endl;
    }

    try
    {
        //std::cout << n_line << " " << trim(token2) << std::endl;
        index = stoi(trim(token2));
    }
    catch (std::exception const &e)
    {
        std::cerr << "SYNTAX ERROR :: line " << n_line << " -> " << e.what() << std::endl;
    }

    token3 = trim(token3);
    int id_prop = stoi(token3);

    int num_s = 0;
    int num_ps = 0;
    int num_vec = 0;

    double val = 0;

    token4 = trim(token4);

    //std::cout << token4 << std::endl;

    bool end_of_line_reached = false;

    do
    {
        pos = token4.find(delimiter_2);

        if (pos == std::string::npos)
        {
            token5 = token4;
            end_of_line_reached = true;
        }
        else
        {
            //std::cout << "token4 = " << token4 << " " << pos << std::endl;
            //std::cout << delimiter_2.length() << std::endl;
            //std::cout << trim(token5) << std::endl;

            token5 = token4.substr(0, pos);
            token4.erase(0, pos + delimiter_2.length());
            token4 = trim(token4);
        }
        //std::cout << "token4 = " << token4 << std::endl;
        //std::cout << "token5 = " << token5 << std::endl;

        pos2 = token5.find(delimiter_4);
        token6 = token5.substr(0, pos2);
        token5.erase(0, pos2 + delimiter_4.length());
        token7 = token5;

        token6 = trim(token6);

        pos2 = token6.find(delimiter_5);
        std::string coupling_type = token6.substr(0, pos2);
        token6.erase(0, pos2 + delimiter_5.length());
        coupling_type = trim(coupling_type);

        pos2 = token6.find(delimiter_5);
        std::string ind_dm_1_str = token6.substr(0, pos2);
        int ind_dm_1 = stoi(trim(ind_dm_1_str));
        token6.erase(0, pos2 + delimiter_5.length());
        std::string ind_dm_2_str = token6;
        int ind_dm_2 = stoi(trim(ind_dm_2_str));

        token7 = trim(token7);

        try
        {
            //std::cout << "oui : " << ind_dm_1 << " " << token6 << " | " << token4 << " | " << std::endl;
            val = stod(token7);
        }
        catch (std::exception const &e)
        {
            std::cerr << "SYNTAX ERROR :: line " << token7 << " " << n_line << " -> " << e.what() << std::endl;
        }

        if (ftype_DM[index] == Fermiontype::dirac && ind_dm_1 == ind_dm_2 && val != 0)
            std::cout << "WARNING : Trying to couple a Dirac fermion to itself | from " << __PRETTY_FUNCTION__ << std::endl;
        else
        {
            if (prop_type[id_prop] == Proptype::scalar)
                lambda_S_DM[index][glob_ind_to_scalar_ind[id_prop]][ind_dm_1][ind_dm_2] = val;

            if (prop_type[id_prop] == Proptype::pseudoscalar)
                lambda_PS_DM[index][glob_ind_to_pseudoscalar_ind[id_prop]][ind_dm_1][ind_dm_2] = val;

            if (prop_type[id_prop] == Proptype::vector)
            {
                if (coupling_type == "a")
                    a_VEC_DM[index][glob_ind_to_vector_ind[id_prop]][ind_dm_1][ind_dm_2] = val;
                if (coupling_type == "b")
                    b_VEC_DM[index][glob_ind_to_vector_ind[id_prop]][ind_dm_1][ind_dm_2] = val;
            }
        }

    } while (!end_of_line_reached);
}

void InputReader::ReadIntDSDS(std::string &interactions, int n_line)
{
    //std::cout << "here!" << std::endl;
    interactions = trim(interactions);

    size_t pos, pos2;
    std::string token1, token2, token3, token4, token5, token6, token7;
    std::string delimiter_1 = "::";
    std::string delimiter_2 = ",";
    std::string delimiter_3 = " ";
    std::string delimiter_4 = "=";
    std::string delimiter_5 = "_";
    int index;

    try
    {
        pos = interactions.find(delimiter_1);
        token1 = interactions.substr(0, pos);
        interactions.erase(0, pos + delimiter_1.length());

        pos = interactions.find(delimiter_1);
        token2 = interactions.substr(0, pos);
        interactions.erase(0, pos + delimiter_1.length());
        token3 = interactions;

        //std::cout << "|" << interactions << "|" << std::endl;
        /*
        std::cout << interactions << std::endl;
        
        pos = interactions.find(delimiter_1);
        token1 = interactions.substr(0, pos);
        interactions.erase(0, pos + delimiter_1.length());

        std::cout << interactions << std::endl;

        pos = interactions.find(delimiter_1);
        token2 = interactions.substr(0, pos);
        interactions.erase(0, pos + delimiter_1.length());

        std::cout << "|"  << interactions << "|" << std::endl;

        //pos = interactions.find(delimiter_1);
        token3 = interactions.erase(0, pos+2);

        std::cout << token3 << std::endl;
        interactions.erase(0, pos + delimiter_1.length());

        token4 = interactions;
        std::cout << interactions << std::endl;*/
    }
    catch (std::exception const &e)
    {
        std::cerr << "SYNTAX ERROR :: line " << n_line << " -> " << e.what() << std::endl;
    }

    try
    {
        //std::cout << n_line << " " << trim(token2) << std::endl;
        index = stoi(trim(token2));
    }
    catch (std::exception const &e)
    {
        std::cerr << "SYNTAX ERROR :: line " << n_line << " -> " << e.what() << std::endl;
    }

    //std::cout << index << std::endl;

    //token3 = trim(token3);
    /*int ind_prop_1 = stoi(token3);

    if (prop_type[ind_prop_1] != Proptype::scalar)
    {
        std::cerr << "SYNTAX ERROR :: line " << n_line << " -> "
                  << "propagator must be a scalar" << std::endl;
        exit(0);
    }*/

    int num_s = 0;
    int num_ps = 0;
    int num_vec = 0;

    double val = 0;

    //int ind_prop_1 = 0;
    //token3 = trim(token3);

    //std::cout << token3 << std::endl;

    bool end_of_line_reached = false;

    do
    {
        pos = token3.find(delimiter_2);

        if (pos == std::string::npos)
        {
            token4 = token3;
            end_of_line_reached = true;
        }
        else
        {
            //std::cout << "token4 = " << token4 << " " << pos << std::endl;
            //std::cout << delimiter_2.length() << std::endl;
            //std::cout << trim(token5) << std::endl;

            token4 = token3.substr(0, pos);
            token3.erase(0, pos + delimiter_2.length());
            token3 = trim(token3);
            //std::cout << token3 << " | " << token4 << std::endl;
        }

        //std::cout << "token4 = " << token4 << std::endl;
        //std::cout << "token5 = " << token5 << std::endl;

        pos2 = token4.find(delimiter_4);
        token5 = token4.substr(0, pos2);
        token4.erase(0, pos2 + delimiter_4.length());
        token6 = token4;

        //std::cout << token5 << " " << token6 << " " << std::endl;

        token5 = trim(token5);

        pos2 = token5.find(delimiter_5);
        std::string coupling_type = token5.substr(0, pos2);
        token5.erase(0, pos2 + delimiter_4.length());
        coupling_type = trim(coupling_type);

        //std::cout << token5 << std::endl;

        pos2 = token5.find(delimiter_5);
        std::string ind_prop_1_str = token5.substr(0, pos2);
        int ind_prop_1 = stoi(trim(ind_prop_1_str));
        token5.erase(0, pos2 + delimiter_5.length());

        //std::cout << token5 << std::endl;

        pos2 = token5.find(delimiter_5);
        std::string ind_prop_2_str = token5.substr(0, pos2);
        int ind_prop_2 = stoi(trim(ind_prop_2_str));
        token5.erase(0, pos2 + delimiter_5.length());

        pos2 = token5.find(delimiter_5);
        std::string ind_prop_3_str = token5;
        int ind_prop_3 = stoi(trim(ind_prop_3_str));

        token6 = trim(token6); // value of the coupling

        //std::cout << ind_prop_1_str << " " << ind_prop_2_str << " " << ind_prop_3_str << " " << token6 << std::endl;

        try
        {
            //std::cout << "oui : " << token7 << " " << token6 << " | " << token4 << " | " << std::endl;
            val = stod(token6);
        }
        catch (std::exception const &e)
        {
            std::cerr << "SYNTAX ERROR :: line " << n_line << " -> " << e.what() << std::endl;
        }

        if (coupling_type == "csss" && prop_type[ind_prop_1] == Proptype::scalar && prop_type[ind_prop_2] == Proptype::scalar && prop_type[ind_prop_3] == Proptype::scalar)
        {
            //std::cout << "Here : " << val << std::endl;
            c_sss[index][glob_ind_to_scalar_ind[ind_prop_1]][glob_ind_to_scalar_ind[ind_prop_2]][glob_ind_to_scalar_ind[ind_prop_3]] = val;
        }
        else if (coupling_type == "csss" && (prop_type[ind_prop_1] != Proptype::scalar || prop_type[ind_prop_2] != Proptype::scalar || prop_type[ind_prop_3] != Proptype::scalar))
        {
            std::cerr << "SYNTAX ERROR :: line " << n_line << " -> "
                      << "with csss propagators must all be scalars" << std::endl;
            exit(0);
        }

        if (coupling_type == "dspp" && prop_type[ind_prop_1] == Proptype::scalar && prop_type[ind_prop_2] == Proptype::pseudoscalar && prop_type[ind_prop_3] == Proptype::pseudoscalar)
        {
            d_spp[index][glob_ind_to_scalar_ind[ind_prop_1]][glob_ind_to_pseudoscalar_ind[ind_prop_2]][glob_ind_to_pseudoscalar_ind[ind_prop_3]] = val;
        }
        else if (coupling_type == "dspp" && (prop_type[ind_prop_1] != Proptype::scalar || prop_type[ind_prop_2] != Proptype::pseudoscalar || prop_type[ind_prop_3] != Proptype::pseudoscalar))
        {
            //std::cout << coupling_type << " " << ind_prop_2 << " " << ind_prop_3 << std::endl;
            std::cerr << "SYNTAX ERROR :: line " << n_line << " -> "
                      << "with dspp the last two propagators must be pseudoscalars" << std::endl;
            exit(0);
        }

        if (coupling_type == "gsXX" && prop_type[ind_prop_2] == Proptype::vector && prop_type[ind_prop_3] == Proptype::vector)
        {
            g_sXX[index][glob_ind_to_scalar_ind[ind_prop_1]][glob_ind_to_vector_ind[ind_prop_2]][glob_ind_to_vector_ind[ind_prop_3]] = val;
        }
        else if (coupling_type == "gsXX" && (prop_type[ind_prop_2] != Proptype::vector || prop_type[ind_prop_3] != Proptype::vector))
        {
            std::cerr << "SYNTAX ERROR :: line " << n_line << " -> "
                      << "with gsXX the last two propagators must be vectors" << std::endl;
            exit(0);
        }
    } while (!end_of_line_reached);

    //std::cout << "reached the end" << std::endl;
}

//
//
//
//

void InputReader::Symetrise2(std::vector<std::vector<double> > &vec)
{

    for (int i = 0; i < vec.size(); i++)
    {
        if (vec[i].size() != vec.size())
        {
            std::cerr << "FATAL ERROR : " << __PRETTY_FUNCTION__ << " -> Check length" << std::endl;
            exit(0);
        }

        for (int j = i + 1; j < vec[i].size(); j++)
        {
            if (vec[i][j] != 0)
                vec[j][i] = vec[i][j];

            if (vec[j][i] != 0)
                vec[i][j] = vec[j][i];
        }
    }
}

void InputReader::Symetrise3(std::vector<std::vector<std::vector<double> > > &vec)
{

    for (int i = 0; i < vec.size(); i++)
    {
        if (vec[i].size() != vec.size())
        {
            std::cerr << "FATAL ERROR : " << __PRETTY_FUNCTION__ << " -> Check length" << std::endl;
            exit(0);
        }

        for (int j = 0; j < vec[i].size(); j++)
        {

            if (vec[i][j].size() != vec[i].size())
            {
                std::cerr << "FATAL ERROR : " << __PRETTY_FUNCTION__ << " -> Check length" << std::endl;
                exit(0);
            }

            for (int k = 0; k < vec[i][j].size(); k++)
            {

                if (vec[i][j][k] != 0)
                {
                    vec[k][i][j] = vec[i][j][k];
                    vec[j][k][i] = vec[i][j][k];
                    vec[j][i][k] = vec[i][j][k];
                    vec[k][j][i] = vec[i][j][k];
                    vec[i][k][j] = vec[i][j][k];
                }

                if (vec[k][i][j] != 0)
                {
                    vec[i][j][k] = vec[k][i][j];
                    vec[j][k][i] = vec[k][i][j];
                    vec[j][i][k] = vec[k][i][j];
                    vec[k][j][i] = vec[k][i][j];
                    vec[i][k][j] = vec[k][i][j];
                }

                if (vec[j][k][i] != 0)
                {
                    vec[i][j][k] = vec[j][k][i];
                    vec[k][i][j] = vec[j][k][i];
                    vec[j][i][k] = vec[j][k][i];
                    vec[k][j][i] = vec[j][k][i];
                    vec[i][k][j] = vec[j][k][i];
                }

                if (vec[j][i][k] != 0)
                {
                    vec[i][j][k] = vec[j][i][k];
                    vec[k][i][j] = vec[j][i][k];
                    vec[j][k][i] = vec[j][i][k];
                    vec[k][j][i] = vec[j][i][k];
                    vec[i][k][j] = vec[j][i][k];
                }

                if (vec[k][j][i] != 0)
                {
                    vec[i][j][k] = vec[k][j][i];
                    vec[k][i][j] = vec[k][j][i];
                    vec[j][k][i] = vec[k][j][i];
                    vec[j][i][k] = vec[k][j][i];
                    vec[i][k][j] = vec[k][j][i];
                }

                if (vec[i][k][j] != 0)
                {
                    vec[i][j][k] = vec[i][k][j];
                    vec[k][i][j] = vec[i][k][j];
                    vec[j][k][i] = vec[i][k][j];
                    vec[j][i][k] = vec[i][k][j];
                    vec[k][j][i] = vec[i][k][j];
                }
            }
        }
    }
}

void InputReader::order1(std::vector<double> &vec)
{
    int n = vec.size();
    int min = 0;
    double trans;

    //std::cout << n << "" " << std::endl;

    for (int i = 0; i < n - 1; i++)
    {
        min = i;
        for (int j = i + 1; j < n; j++)
        {
            if (vec[j] < vec[min])
                min = j;
        }

        if (min != i)
        {
            trans = vec[i];
            vec[i] = vec[min];
            vec[min] = trans;
        }

        //for(int j = 0; j<n; j++)
        //    std::cout << i << " " << vec[j] << std::endl;
    }
}

std::string &InputReader::ltrim(std::string &str, const std::string &chars)
{
    str.erase(0, str.find_first_not_of(chars));
    return str;
}

std::string &InputReader::rtrim(std::string &str, const std::string &chars)
{
    str.erase(str.find_last_not_of(chars) + 1);
    return str;
}

std::string &InputReader::trim(std::string &str, const std::string &chars)
{
    return ltrim(rtrim(str, chars), chars);
}
