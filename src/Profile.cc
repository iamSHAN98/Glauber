#define Pi M_PI
#include <cmath>
#include "Profile.h"

// Spherical harmonics
double Y20(double x){
  return sqrt(5./(16*Pi))*(3*pow(cos(x), 2) - 1.);
}

double Y40(double x){
  return (3./(16*sqrt(Pi)))*(35*pow(cos(x), 4) - 30*pow(cos(x), 2) + 3.);
}

namespace EventGen{

  double ChargeProfile :: GetProfile(double *R){
    switch (N.Profile) {
      case KeyChargeProfile::Fermi    : return Fermi(R);
      case KeyChargeProfile::Deformed : return Deformed(R);
      case KeyChargeProfile::Gaussian : return Gaussian(R);
      case KeyChargeProfile::Hulthen  : return Hulthen(R);
      default :
        std::cout << "Can't find given nucleus, exiting \n";
        exit(1);
    }
  }

  void ChargeProfile :: Normalize(){
    Rho = 1.;
    double Min[3] = {}, Max[] = {2*N.R, Pi, 2*Pi};
    double T = Int.Integrate3D<ChargeProfile>
               (this, &ChargeProfile::GetProbability, Min, Max);
    Rho = N.A/T;
  }

  double ChargeProfile :: GetProfile(Position R){
    double Ri[] = {R.Norm(), R.Theta(), R.Phi()};
    return GetProfile(Ri);
  }

  double ChargeProfile :: GetProbability(double *R){
    return R[0]*R[0]*sin(R[1])*GetProfile(R);
  }

  double ChargeProfile :: GetProbability(Position R){
    double Ri[] = {R.Norm(), R.Theta(), R.Phi()};
    return GetProbability(Ri);
  }

  // Definitions from arXiv:1408.2549 [nucl-ex]
  double ChargeProfile :: Fermi(double *R){
    double A = 1. + N.w*pow(R[0]/N.R, 2);
    double B = 1. + exp((R[0] - N.R)/N.a);
    return Rho*A/B;
  }

  double ChargeProfile :: Deformed(double *R){
    double Rd = N.R*(1. + N.b2*Y20(R[1]) + N.b4*Y40(R[1]));
    double B = (1. + exp((R[0] - Rd)/N.a));
    return Rho/B;
  }

  double ChargeProfile :: Gaussian(double *R){
    double A = 1. + N.w*pow(R[0]/N.R, 2);
    double B = (R[0]*R[0] - N.R*N.R)/(N.a*N.a);
    return Rho*A/(1. + exp(B));
  }

  double ChargeProfile :: Hulthen(double *R){
    double r = R[0];
    double B = exp(-N.a*r) - exp(-N.b2*r);
    return Rho*(B*B)/(r*r + 1e-10);
  }

}
