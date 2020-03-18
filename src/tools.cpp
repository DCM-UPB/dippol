#include "tools.hpp"

constexpr unsigned int str2int(const char* str, int h = 0)
{
  return !str[h] ? 5381 : (str2int(str, h+1) * 33) ^ str[h];
}

double deduceMass(const std::string name)
{
  std::string name_uc = name;  
  std::for_each(name_uc.begin(), name_uc.end(), [](char & c){c = ::toupper(c);});
  double mass = 0.0;
  switch(str2int(name_uc.c_str()))
    {
    case str2int("H"): mass=1.0079;
      break;
    case str2int("O"): mass=15.9994;
      break;
    case str2int("C"): mass=12.011;
      break;
    case str2int("N"): mass=14.007;
      break;
    case str2int("P"): mass=30.973762;
      break;
    case str2int("CL"): mass=35.45;
      break;
    case str2int("NA"): mass=22.98976928;
      break;
    case str2int("S"): mass=32.06;
      break;
    case str2int("K"): mass=39.0983;
      break;
    case str2int("CA"): mass=40.078;
      break;
    case str2int("MN"): mass=54.938045;
      break;
    case str2int("FE"): mass=55.845;
      break;
    case str2int("MG"): mass=24.3050;
      break;
    case str2int("ZN"): mass=65.38;
      break;
    case str2int("HE"): mass=4.002602;
      break;
    case str2int("LI"): mass=6.94;
      break;
    case str2int("BE"): mass=9.012182;
      break;
    case str2int("B"): mass=10.81;
      break;
    case str2int("F"): mass=18.9984032;
      break;
    case str2int("NE"): mass=20.1797;
      break;
    case str2int("AL"): mass=26.9815386;
      break;
    case str2int("SI"): mass=28.085;
      break;
    case str2int("AR"): mass=39.948;
      break;
    case str2int("SC"): mass=44.955912;
      break;
    case str2int("TI"): mass=47.867;
      break;
    case str2int("V"): mass=50.9415;
      break;
    case str2int("CR"): mass=51.9961;
      break;
    case str2int("NI"): mass=58.6934;
      break;
    case str2int("CO"): mass=58.933195;
      break;
    case str2int("CU"): mass=63.546;
      break;
    case str2int("GA"): mass=69.723;
      break;
    case str2int("GE"): mass=72.63;
      break;
    case str2int("AS"): mass=74.92160;
      break;
    case str2int("SE"): mass=78.96;
      break;
    case str2int("BR"): mass=79.904;
      break;
    case str2int("KR"): mass=83.798;
      break;
    case str2int("RB"): mass=85.4678;
      break;
    case str2int("SR"): mass=87.62;
      break;
    case str2int("Y"): mass=88.90585;
      break;
    case str2int("ZR"): mass=91.224;
      break;
    case str2int("NB"): mass=92.90638;
      break;
    case str2int("MO"): mass=95.96;
      break;
    case str2int("TC"): mass=98;
      break;
    case str2int("RU"): mass=101.07;
      break;
    case str2int("RH"): mass=102.90550;
      break;
    case str2int("PD"): mass=106.42;
      break;
    case str2int("AG"): mass=107.8682;
      break;
    case str2int("CD"): mass=112.411;
      break;
    case str2int("IN"): mass=114.818;
      break;
    case str2int("SN"): mass=118.710;
      break;
    case str2int("SB"): mass=121.760;
      break;
    case str2int("I"): mass=126.90447;
      break;
    case str2int("TE"): mass=127.60;
      break;
    case str2int("XE"): mass=131.293;
      break;
    case str2int("CS"): mass=132.9054519;
      break;
    case str2int("BA"): mass=137.327;
      break;
    case str2int("LA"): mass=138.90547;
      break;
    case str2int("CE"): mass=140.116;
      break;
    case str2int("PR"): mass=140.90765;
      break;
    case str2int("ND"): mass=144.242;
      break;
    case str2int("PM"): mass=145;
      break;
    case str2int("SM"): mass=150.36;
      break;
    case str2int("EU"): mass=151.964;
      break;
    case str2int("GD"): mass=157.25;
      break;
    case str2int("TB"): mass=158.92535;
      break;
    case str2int("DY"): mass=162.500;
      break;
    case str2int("HO"): mass=164.93032;
      break;
    case str2int("ER"): mass=167.259;
      break;
    case str2int("TM"): mass=168.93421;
      break;
    case str2int("YB"): mass=173.054;
      break;
    case str2int("LU"): mass=174.9668;
      break;
    case str2int("HF"): mass=178.49;
      break;
    case str2int("TA"): mass=180.94788;
      break;
    case str2int("W"): mass=183.84;
      break;
    case str2int("RE"): mass=186.207;
      break;
    case str2int("OS"): mass=190.23;
      break;
    case str2int("IR"): mass=192.217;
      break;
    case str2int("PT"): mass=195.084;
      break;
    case str2int("AU"): mass=196.966569;
      break;
    case str2int("HG"): mass=200.59;
      break;
    case str2int("TL"): mass=204.38;
      break;
    case str2int("PB"): mass=207.2;
      break;
    case str2int("BI"): mass=208.98040;
      break;
    case str2int("PO"): mass=209;
      break;
    case str2int("AT"): mass=210;
      break;
    case str2int("RN"): mass=222;
      break;
    case str2int("FR"): mass=223;
      break;
    case str2int("RA"): mass=226;
      break;
    case str2int("AC"): mass=227;
      break;
    case str2int("PA"): mass=231.03588;
      break;
    case str2int("TH"): mass=232.03806;
      break;
    case str2int("NP"): mass=237;
      break;
    case str2int("U"): mass=238.02891;
      break;
    case str2int("AM"): mass=243;
      break;
    case str2int("PU"): mass=244;
      break;  
    case str2int("CM"): mass=247;
      break;
    } 
  return mass;
}

Eigen::Vector3d gcoord_to_eigenv(const loos::GCoord in)
{
   Eigen::Vector3d out(in.x(),in.y(),in.z());
   return out;
}

Eigen::Matrix3d gmatrix_to_eigenm(const loos::GMatrix in)
{
  Eigen::Matrix3d eigenm;
  eigenm << in(0,0), in(0,1), in(0,2),
    in(1,0), in(1,1), in(1,2),
    in(2,0), in(2,1), in(2,2);
  return eigenm;
}

Eigen::Matrix3d inertia_tensor(const loos::AtomicGroup group)
{
  Eigen::Matrix3d toi = Eigen::Matrix3d::Zero();
  for (int i=0; i<group.size(); i++)
  {
    const loos::GCoord r = group[i]->coords() - group.centerOfMass();
    const double x = r.x();
    const double y = r.y();
    const double z = r.z();
    const double m = group[i]->mass();
    toi(0,0) += (m * (y*y + z*z));
    toi(1,1) += (m * (x*x + z*z));
    toi(2,2) += (m * (x*x + y*y));
    toi(0,1) += (m * x * y);
    toi(0,2) += (m * x * z);
    toi(1,2) += (m * y * z);
  }
  toi(0,1) = -1*toi(0,1);
  toi(0,2) = -1*toi(0,2);
  toi(1,2) = -1*toi(1,2);
  toi(1,0) = toi(0,1);
  toi(2,0) = toi(0,2);
  toi(2,1) = toi(1,2);
  return toi;
}

void writexyz(const loos::AtomicGroup group, std::ostream &output)
{
  output << group.size() << std::endl << std::endl;
  loos::AtomicGroup::Iterator it(group);
  loos::pAtom current;
  while(current = it())
  {
    output.width(3);
    output << current->name();
    output.width(21);
    output << current->coords().x();
    output.width(21);
    output << current->coords().y();
    output.width(21);
    output << current->coords().z()
           << std::endl;  
  }
}


// Eigen::Matrix3d align(loos::AtomicGroup ref, loos::AtomicGroup mol, Eigen::Matrix<double,Eigen::Dynamic,3> w)
// {
//     if(ref.size() != mol.size())
//     {
//         cout << "Not possible to align two molecules of different sizes" << endl;
//         exit(0);
//     }
//     
//     if (w.rows() != ref.size())
//     {
//         cout << "Weight matrix should have the same size as molecules to be aligned!" << endl;
//         exit(0);
//     }
//     
//     Eigen::Map<Eigen::MatrixX<double,ref.size(),3>>
//     
//     Eigen::MatrixXd<ref.size(),3> M;
//     // M.resize(ref.size(),3);
//     
//     
//     Eigen::Matrix3d R;
//     
//     return R;
// }
