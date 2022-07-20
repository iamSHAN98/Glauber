#define Pi M_PI
#include "Glauber.h"
using namespace Utility;

namespace EventGen{

  void Glauber :: Configure(YAML::Node Param){
    f = Param["f"].as<double>();
    Npp = Param["Npp"].as<double>();
  }

  bool Glauber :: GenerateEvent(){
    ObtainImpactVector();
    ObtainRotationAngle();

    CalculateNPartNColl();
    if(Npart[0] + Npart[1] < 2.) return false;

    CalculateMultiplicity();
    for(int n = 0; n < Ord; n++)
      CalculateEccentricity(n + 2);

    return true;
  }

  void Glauber :: CalculateNPartNColl(){
    double F = 0.5*(1. - f);

    Npart[0] = Npart[1] = Ncoll = 0.;
    Rcm = {};

    Position R;
    double Tm, Tp, Wa, Wb, Wc;
    double Ya, Yb, Yc, Yx, Yy;

    for(int i = 0; i < Fa.Nq; i++){
      R.x = R0*Arr[i].Node;

      Ya = 0., Yb = 0., Yc = 0.;    // Y integrand for Npart & Ncoll
      Yx = 0., Yy = 0.;             // Y integrand for CoM

      for(int j = 0; j < Fa.Nq; j++){
        R.y = R0*Arr[j].Node;

        Thickness(R, Tm, Tp);
        Wa = Tm*(1. - exp(-SgNN*Tp));
        Wb = Tp*(1. - exp(-SgNN*Tm));
        Wc = SgNN*Tm*Tp;

        // Y integral
        Ya += Wa*Arr[j].Weight;
        Yb += Wb*Arr[j].Weight;
        Yc += Wc*Arr[j].Weight;
        Yx += (F*(Wa + Wb) + f*Wc)*Arr[j].Weight;
        Yy += R.y*(F*(Wa + Wb) + f*Wc)*Arr[j].Weight;
      }

      // X integral
      Npart[0] += Ya*Arr[i].Weight;
      Npart[1] += Yb*Arr[i].Weight;
      Ncoll += Yc*Arr[i].Weight;
      Rcm.x += R.x*Yx*Arr[i].Weight;
      Rcm.y += Yy*Arr[i].Weight;
    }

    // Dimension correction
    Npart[0] *= R0*R0, Npart[1] *= R0*R0, Ncoll *= R0*R0;
    Rcm *= R0*R0/(F*(Npart[0] + Npart[1]) + f*Ncoll + 1e-10);

    // Spectator neutron count
    Spec[0] = (Na.A - Npart[0])*(1. - Na.Z/Na.A);
    Spec[1] = (Nb.A - Npart[1])*(1. - Nb.Z/Nb.A);
  }

  void Glauber :: Thickness(Position R, double& Tm, double& Tp){
    Position Rm, Rp;

    Tm = 0., Tp = 0.;
    for(int i = 0; i < Fa.Nq; i++){
      Rm.z = Rp.z = R0*Arr[i].Node;

      // Nucleus A
      Rm.x = R.x - 0.5*B.x, Rm.y = R.y - 0.5*B.y;
      EulerRotation(Rm, Theta[0], Phi[0]);

      // Nucleus B
      Rp.x = R.x + 0.5*B.x, Rp.y = R.y + 0.5*B.y;
      EulerRotation(Rp, Theta[1], Phi[1]);

      Tm += Fa.GetProfile(Rm)*Arr[i].Weight;
      Tp += Fb.GetProfile(Rp)*Arr[i].Weight;
    }

    // Dimension correction
    Tm *= R0, Tp *= R0;
  }

  void Glauber :: EulerRotation(Position& R, double Th, double Ph){
    double x = R.x, y = R.y, z = R.z;
    R.x = x*cos(Th)*cos(Ph) + z*cos(Ph)*sin(Th) - y*sin(Ph);
    R.y = y*cos(Ph) + x*cos(Th)*sin(Ph) + z*sin(Th)*sin(Ph);
    R.z = z*cos(Th) - x*sin(Th);
  }

  void Glauber :: CalculateMultiplicity(){
    Nch = (0.5*(1. - f)*(Npart[0] + Npart[1]) + f*Ncoll)*Npp;
  }

  void Glauber :: CalculateEccentricity(int n){
    double e, r, phi;
    double NumC = 0., NumS = 0., Denom = 0.;

    Position R;                           // Shift by -Rcm
    for(int i = 0; i < Fa.Nq; i++){
      R.x = R0*Arr[i].Node - Rcm.x;

      double SumC = 0., SumS = 0., SumD = 0.;

      for(int j = 0; j < Fa.Nq; j++){
        R.y = R0*Arr[j].Node - Rcm.y;

        e = GetTransverseProfile(R), r = R.TransverseNorm(), phi = R.Phi();

        SumS += e*pow(r, n)*sin(n*phi)*Arr[j].Weight;
        SumC += e*pow(r, n)*cos(n*phi)*Arr[j].Weight;
        SumD += e*pow(r, n)*Arr[j].Weight;
      }

      NumS += SumS*Arr[i].Weight;
      NumC += SumC*Arr[i].Weight;
      Denom += SumD*Arr[i].Weight;
    }

    Ecc[n - 2] = sqrt(NumC*NumC + NumS*NumS)/(Denom + 1e-10);
    Psi[n - 2] = (1./n)*atan2(NumS, NumC);
  }

  double Glauber :: GetTransverseProfile(Position R){
    double Tm, Tp;
    Thickness(R, Tm, Tp);

    double ePart = Tm*(1. - exp(-SgNN*Tp)) + Tp*(1. - exp(-SgNN*Tm));
    double eColl = SgNN*Tm*Tp;
    return 0.5*(1. - f)*ePart + f*eColl;
  }

  void Glauber :: SetOutputFormat(int Nev){
    // Initialize DataStream
    Output.Open(FileName.str(), DataStream::Write);

    // Datasets
    Output.Add("Glauber/b", &b, DataStream::Double);
    Output.Add("Glauber/Theta", &Theta, DataStream::Double, {2});
    Output.Add("Glauber/Phi", &Phi, DataStream::Double, {2});
    Output.Add("Glauber/Npart", &Npart, DataStream::Double, {2});
    Output.Add("Glauber/Spec", &Spec, DataStream::Double, {2});
    Output.Add("Glauber/Ncoll", &Ncoll, DataStream::Double);
    Output.Add("Glauber/Nch", &Nch, DataStream::Double);
    Output.Add("Glauber/Ecc", &Ecc, DataStream::Double, {Ord});
    Output.Add("Glauber/Psi", &Psi, DataStream::Double, {Ord});

    Output.Configure(Nev);

    // Model parameters
    Output.SetAttribute("Glauber", "sNN", &sNN, DataStream::Double);
    Output.SetAttribute("Glauber", "f", &f, DataStream::Double);
    Output.SetAttribute("Glauber", "Npp", &Npp, DataStream::Double);

    // Nucleus parameters
    DataStream::MetaData InfoN(DataStream::Compound, {1}, "Nucleus");
    InfoN.AddMember(&NucleusData::Name, "Name", DataStream::String);
    InfoN.AddMember(&NucleusData::A, "A", DataStream::Integer);
    InfoN.AddMember(&NucleusData::Z, "Z", DataStream::Integer);  

    Output.SetAttribute("Glauber", "NucleusA", &Na, InfoN);
    Output.SetAttribute("Glauber", "NucleusB", &Nb, InfoN);
  }

  void Glauber :: ObtainRotationAngle(){
    for(int i = 0; i < 2; i++)
      Theta[i] = Rng.Sine(), Phi[i] = Rng.Uniform(0., 2*Pi);  
  }

  // Expression from arXiv:1710.07098 [nucl-ex]  
  double Glauber :: CrossSection(double x){
    return 0.1*(29.8 + 0.038*pow(log(x*x), 2.43));
  }

  void Glauber :: ObtainImpactVector(){
    b = Rng.Linear(Bmin, Bmax);
    double ph = Rng.Uniform(0., 2*Pi);
    double x = b*cos(ph), y = b*sin(ph);
    B = Position{x, y, 0.};
  }

  double Glauber :: GetTheta(int i){
    if((i < 1) || (i > 2)) return 0.;
    return Theta[i - 1];
  }

  double Glauber :: GetPhi(int i){
    if((i < 1) || (i > 2)) return 0.;
    return Phi[i - 1];
  }

  double Glauber :: GetNPart(int i){
    if((i < 1) || (i > 2)) return Npart[0] + Npart[1];
    return Npart[i - 1];
  }

  double Glauber :: GetSpecNeutron(int i){
    if((i < 1) || (i > 2)) return Spec[0] + Spec[1];
    return Spec[i - 1];
  }

  double Glauber :: GetEccentricity(int i){
    if((i < 2) || (i > Ord + 2)) return 0.;
    return Ecc[i - 2];
  }

  double Glauber :: GetPlaneAngle(int i){
    if((i < 2) || (i > Ord + 2)) return 0.;
    return Psi[i - 2];
  }

}
