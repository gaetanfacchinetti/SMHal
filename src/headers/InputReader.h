#ifndef DEF_INPUTREADER
#define DEF_INPUTREADER

#include "CppLib.h"
#include "DarkSectorModel.h"

class InputReader
{
  public:
    static InputReader *getInstance()
    {
        if (NULL == _IR)
            _IR = new InputReader;

        return _IR;
    }

    static void kill()
    {
        if (NULL != _IR)
        {
            delete _IR;
            _IR = NULL;
        }
    }

    std::vector<DarkSectorModel> Read(std::string name_file);

  private:

    InputReader();
    ~InputReader(){};

    void ReadDSContent(std::string &content, int n_line);
    void ReadTypeDM(std::string &content, int n_line);
    void ReadMassDMDS(std::string &content, int n_line);
    void ReadIntDSSM(std::string &interactions, int n_line);
    void ReadIntDSDM(std::string &interactions, int n_line);
    void ReadIntDSDS(std::string &interactions, int n_line);

    void CompleteCouplingTables_lambda_S_SM(std::string const &name, double val, int id, int id_part);
    void CompleteCouplingTables_lambda_PS_SM(std::string const &name, double val, int id, int id_part);
    void CompleteCouplingTables_ab_VEC_SM(std::string const &name, double val, int id, int id_part);

    void order1(std::vector<double> & vec);
    void Symetrise2(std::vector<std::vector<double>> & vec);
    void Symetrise3(std::vector<std::vector<std::vector<double>>> & vec);
    std::string &ltrim(std::string &str, const std::string &chars = "\t\n\v\f\r ");
    std::string &rtrim(std::string &str, const std::string &chars = "\t\n\v\f\r ");
    std::string &trim(std::string &str, const std::string &chars = "\t\n\v\f\r ");

    static InputReader *_IR;

    int Npoints;
    int n_dm_part;
    int n_prop_scalar;
    int n_prop_pseudoscalar;
    int n_prop_vector;

    std::vector<int> dof_prop_scalar;
    std::vector<int> dof_prop_pseudoscalar;
    std::vector<int> dof_prop_vector;

    std::vector<Fermiontype> ftype_DM;

    std::vector<std::vector<double>> mass_DM;
    std::vector<std::vector<double>> mass_prop_scalar;
    std::vector<std::vector<double>> mass_prop_pseudoscalar;
    std::vector<std::vector<double>> mass_prop_vector;

    std::vector<std::vector<std::vector<double>>> lambda_S_SM, lambda_PS_SM, a_VEC_SM, b_VEC_SM;
    std::vector<std::vector<std::vector<std::vector<double>>>> lambda_S_DM, lambda_PS_DM, a_VEC_DM, b_VEC_DM;
    std::vector<std::vector<std::vector<std::vector<double>>>> c_sss, d_spp, g_sXX;

    std::vector<Proptype> prop_type;

    std::vector<int> glob_ind_to_scalar_ind, glob_ind_to_pseudoscalar_ind, glob_ind_to_vector_ind;

    //std::vector<DarkSectorModel> _DS_vec;
};

#endif