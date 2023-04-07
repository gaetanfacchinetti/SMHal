#include "../headers/Particle.h"

Particle::Particle(double n_spin, int n_degFree, std::string n_name, double n_mass, double n_width, double n_electric_charge, int n_ind, QCDTYPE n_qtype)
{
    spin = n_spin;
    qtype = n_qtype;
    name = n_name;

    mass = n_mass;
    width = n_width;

    degFree = n_degFree;

    electric_charge = n_electric_charge;

    specific_ind = 0; 

    if (floor(spin) == spin)
        stat = +1;
    else
        stat = -1;

    if (qtype != QCDTYPE::MESON && qtype != QCDTYPE::BARYON)
        isPresentBeforeQCDPT = true;
    else
        isPresentBeforeQCDPT = false;

    if (qtype != QCDTYPE::GLUON && qtype != QCDTYPE::QUARK)
        isPresentAfterQCDPT = true;
    else
        isPresentAfterQCDPT = false;

    ind = n_ind;
}

Particle::Particle(double n_spin, int n_degFree, std::string n_name, double n_mass,  double n_width, double n_electric_charge, int n_ind, Proptype n_ptype, int n_specific_ind, Fermiontype n_ftype, QCDTYPE n_qtype)
{
    spin = n_spin;
    qtype = n_qtype;
    name = n_name;

    mass = n_mass;
    width = n_width;

    degFree = n_degFree;

    electric_charge = n_electric_charge;

    specific_ind = n_specific_ind;

    if (floor(spin) == spin)
        stat = +1;
    else
        stat = -1;

    if (qtype != QCDTYPE::MESON && qtype != QCDTYPE::BARYON)
        isPresentBeforeQCDPT = true;
    else
        isPresentBeforeQCDPT = false;

    if (qtype != QCDTYPE::GLUON && qtype != QCDTYPE::QUARK)
        isPresentAfterQCDPT = true;
    else
        isPresentAfterQCDPT = false;

    ptype = n_ptype;

    ftype = n_ftype;

    ind = n_ind;
}

/*
bool Particle::isEqualTo(Particle const& b) const
{
 
    if (stat == b.stat && mass == b.mass && name == b.name &&width == b.width && spin == b.spin &&  type == b.type && degFree == b.degFree && ind == b.ind 
        && specific_ind == b.specific_ind && isPresentAfterQCDPT == b.isPresentAfterQCDPT && isPresentBeforeQCDPT == b.isPresentBeforeQCDPT && qtype == b.qtype && ptype == b.ptype)
        return true;
    else
        return false;
}*/

ParticlePair::ParticlePair(Particle & n_part1, Particle & n_part2, int n_index, INT_TYPE n_int_type)
{
    part1 = n_part1;
    part2 = n_part2;
    index = n_index;

    int_type = n_int_type;
}


std::ostream &operator<<(std::ostream &os, const ParticlePair pp)
{

    os << "ParticlePair index = " << pp.ind() << " | "; 

    if (pp.get_type() == INT_TYPE::SS)
        os << "type : SS";
    if (pp.get_type() == INT_TYPE::TT)
        os << "type : TT";
    if (pp.get_type() == INT_TYPE::UU)
        os << "type : UU";
    if (pp.get_type() == INT_TYPE::ST)
        os << "type : ST";
    if (pp.get_type() == INT_TYPE::SU)
        os << "type : SU";
    if (pp.get_type() == INT_TYPE::TU)
        os << "type : TU";
    
    return os;
}



std::ostream &operator<<(std::ostream &os, const std::vector<ParticlePair *> &v)
{
    os << "[ ";

    if (v.size() > 0)
        for (int i = 0; i < v.size(); i++)
        {
            os << "(" << v[i]->part1.get_name() << "," << v[i]->part2.get_name() << " | ind : " << v[i]->index << ")";
            if (i != v.size() - 1)
                os << " ; ";
        }
    else
    {
        os << "void";
    }

    os << " ]";

    return os;
}

std::ostream &operator<<(std::ostream &os, const std::vector<ParticlePair> &v)
{
    os << "[ ";

    if (v.size() > 0)
        for (int i = 0; i < v.size(); i++)
        {
            os << "(" << v[i].part1.get_name() << "," << v[i].part2.get_name() << " | ind : " << v[i].index << ")";
            if (i != v.size() - 1)
                os << " ; ";
        }
    else
    {
        os << "void";
    }

    os << " ]";

    return os;
}

std::ostream &operator<<(std::ostream &os, const std::vector<Particle *> &v)
{
    for (int i = 0; i < v.size(); i++)
    {
        os.precision(2);
        os << " | Spin : " << v[i]->get_spin();
        os.precision(4);
        os << std::scientific;
        os << " \t |  mass : "
           << v[i]->get_mass() << " \t | width : " << v[i]->get_width()
           << "\t | degFree : " << v[i]->get_degFree();
        os << std::defaultfloat;
        os << "  \t | index : " << v[i]->get_index() << " \t | -> " << v[i]->get_name() << std::endl;
    }

    return os;
}

std::ostream &operator<<(std::ostream &os, const std::vector<Particle> &v)
{
    for (int i = 0; i < v.size(); i++)
    {
        os.precision(2);
        os << " | Spin : " << v[i].get_spin();
        os.precision(4);
        os << std::scientific;
        os << " \t |  mass : "
           << v[i].get_mass() << " \t | width : " << v[i].get_width()
           << "\t | degFree : " << v[i].get_degFree();
        os << std::defaultfloat;
        os << "  \t | index : " << v[i].get_index() << " \t | . " << v[i].get_name() << std::endl;
    }

    return os;
}



/*
bool operator==(Particle const& a, Particle const& b)
{
    return a.isEqualTo(b);
}

bool operator!=(Particle const & a, Particle const &b)
{
  return !(a.isEqualTo(b));
}*/

