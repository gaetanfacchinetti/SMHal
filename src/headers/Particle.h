#ifndef PARTICLE
#define PARTICLE

#include "CppLib.h"
#include "MyUnits.h"

enum class QCDTYPE
{
  QUARK,
  GLUON,
  MESON,
  BARYON,
  NOTHING
};

enum class INT_TYPE
{
  SS,
  ST,
  SU,
  TT,
  TU,
  UU,
  NON
};

enum class Proptype
{
  scalar,
  pseudoscalar,
  vector,
  dm
};

enum class Fermiontype
{
  dirac,
  majorana,
  nothing
};

class Particle
{

public:
  Particle(double spin, int n_degFree, std::string name, double mass, 
           double width, double n_electric_charge, int ind, QCDTYPE qtype = QCDTYPE::NOTHING);
  Particle(double spin, int n_degFree, std::string name, double mass, 
           double width, double n_electric_charge, int ind, Proptype n_ptype, int n_specific_ind = 0, Fermiontype n_ftype = Fermiontype::nothing, QCDTYPE n_qtype = QCDTYPE::NOTHING);

  Particle(){};
  ~Particle(){};

  double get_mass() const { return mass; };
  double get_width() const { return width; };
  double get_spin() const { return spin; };
  QCDTYPE get_QCDType() const { return qtype; };
  Proptype get_proptype() const { return ptype; };
  Fermiontype get_fermiontype() const { return ftype; };
  std::string get_fermiontype_str() const
  {
    if (ftype == Fermiontype::majorana)
      return "Majorana";
    if (ftype == Fermiontype::dirac)
      return "Dirac";
    if (ftype == Fermiontype::nothing)
      return "Nothing";
    return "";
  }
  int get_index() const { return ind; };
  int get_specific_index() const { return specific_ind; };
  int get_degFree() const { return degFree; };
  double get_stat() const { return stat; };
  double get_electric_charge() const { return electric_charge; };
  bool is_PresentBeforeQCDPT() const { return isPresentBeforeQCDPT; };
  bool is_PresentAfterQCDPT() const { return isPresentAfterQCDPT; };
  std::string get_name() const { return name; };

  void set_name(std::string const &n_name) { name = n_name; };
  void set_mass(double const &n_mass) { mass = n_mass; };
  void set_degFree(int const &n_degFree) { degFree = n_degFree; };
  void set_width(double const &n_width) { width = n_width; };
  void set_ftype(Fermiontype const &n_ftype) { ftype = n_ftype; };
  friend std::ostream &operator<<(std::ostream &os, const std::vector<Particle *> &v);
  friend std::ostream &operator<<(std::ostream &os, const std::vector<Particle> &v);

  //bool isEqualTo(Particle const& b) const;

private:
  double mass, width, spin, electric_charge;
  std::string type, name;
  int degFree, ind, specific_ind;
  double stat;
  bool isPresentBeforeQCDPT, isPresentAfterQCDPT;
  QCDTYPE qtype;
  Proptype ptype;
  Fermiontype ftype;
};

class ParticlePair
{
public:
  ParticlePair(){};
  ParticlePair(Particle &part1, Particle &part2, int n_index, INT_TYPE n_int_type = INT_TYPE::NON);
  friend std::ostream &operator<<(std::ostream &os, const std::vector<ParticlePair *> &v);
  friend std::ostream &operator<<(std::ostream &os, const std::vector<ParticlePair> &v);
  friend std::ostream &operator<<(std::ostream &os, const ParticlePair pp);

  ~ParticlePair(){};

  Particle get_first() const { return part1; };
  Particle get_second() const { return part2; };

  Particle get_particle(int selec) const
  {
    if (selec == 0)
      return part1;
    if (selec == 1)
      return part2;
    else
    {
      std::cout << "FATAL ERROR : Trying to access a particle in a particle pair with an indicator greater than 1" << std::endl;
      exit(0);
    }
  }

  INT_TYPE get_type() const { return int_type; };

  int ind() const { return index; };

private:
  Particle part1, part2;
  INT_TYPE int_type;
  int index;
};

/*
class DSParticle : public Particle
{
public:
  DSParticle(double spin, int n_degFree, std::string name, double mass,
             double width, int ind, Proptype n_ptype = Proptype::DM, QCDTYPE qtype = QCDTYPE::NOTHING);
  DSParticle();

  Proptype get_proptype() const { return ptype; };

private:
  Proptype ptype;
};*/

#endif
