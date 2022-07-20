#define Eps 1e-10
#define Pi M_PI
#include "GlauberMC.h"

namespace EventGen{

  void GlauberMC :: Configure(YAML::Node Param){
    Glauber::Configure(Param);       // Base parameters

    // Redefine node
    auto Node = Param["Monte Carlo"];

    NegBinFlag = Node["Negative Binomial"]["Flag"].as<bool>();
    if(NegBinFlag)
      k = Node["Negative Binomial"]["k"].as<double>();

    GammaFlag = Node["Gamma"]["Flag"].as<bool>();
    if(GammaFlag)
      Tpp = Node["Gamma"]["Tpp"].as<double>();
  }

  void GlauberMC :: ShiftMassCenter(){
    double F = 0.5*(1. - f);
    double D = F*(Npart[0] + Npart[1]) + f*Ncoll;
    if(D < Eps) return;

    // Obtain center of mass
    Position Rp = {}, Rc = {};
    for(int i = 0; i < (int)Npart[0]; i++) Rp += PartA[i];
    for(int i = 0; i < (int)Npart[1]; i++) Rp += PartB[i];
    for(int i = 0; i < (int)Ncoll; i++)    Rc += Coll[i];

    // Shift to origin
    auto R = (F*Rp + f*Rc)/D;               
    for(int i = 0; i < (int)Npart[0]; i++) PartA[i] += R;
    for(int i = 0; i < (int)Npart[1]; i++) PartB[i] += R;
    for(int i = 0; i < (int)Ncoll; i++)    Coll[i] += R;
  }

  bool GlauberMC :: GenerateEvent(){
    ObtainImpactVector();
    ObtainRotationAngle();

    NucleusA.GenerateNucleons();
    NucleusA.ShiftNucleus();
    NucleusA.RotateNucleus(Theta[0], Phi[0]);
    NucleusA.ShiftNucleus(-0.5*B);

    NucleusB.GenerateNucleons();
    NucleusB.ShiftNucleus();
    NucleusB.RotateNucleus(Theta[1], Phi[1]);
    NucleusB.ShiftNucleus(0.5*B);

    CalculateNPartNColl();
    if(Npart[0] + Npart[1] < 2.) return false;

    ShiftMassCenter();

    CalculateMultiplicity();
    for(int n = 0; n < Ord; n++)
      CalculateEccentricity(n + 2);

    return true;
  }

  void GlauberMC :: CalculateNPartNColl(){
    int NpA = 0, NpB = 0, Nc = 0;
    Position Ra, Rb, R;

    for(int i = 0; i < Na.A; i++)
      for(int j = 0; j < Nb.A; j++){

        Ra = NucleusA.GetPosition(i);
        Rb = NucleusB.GetPosition(j);
        R = Ra - Rb;
        double dT = pow(R.TransverseNorm(), 2);

        if(SgNN >= Pi*dT){

          Coll[Nc++] = (Na.A*Ra + Nb.A*Rb)/(Na.A + Nb.A);

          if(!NucleusA.IsParticipant(i)){
            NucleusA.TagAsParticipant(i);
            PartA[NpA++] = Ra;
          }

          if(!NucleusB.IsParticipant(j)){
            NucleusB.TagAsParticipant(j);
            PartB[NpB++] = Rb;
          }
        }
      }

    Npart[0] = NpA, Npart[1] = NpB, Ncoll = Nc;

    // Spectator neutron count
    Spec[0] = 0, Spec[1] = 0;

    for(int i = 0; i < Na.A; i++)
      if(!NucleusA.IsParticipant(i))
        if(NucleusA.IsNeutron(i))
          Spec[0] += 1;

    for(int i = 0; i < Nb.A; i++)
      if(!NucleusB.IsParticipant(i))
        if(NucleusB.IsNeutron(i))
          Spec[1] += 1;
  }

  void GlauberMC :: CalculateMultiplicity(){
    if(!NegBinFlag) Glauber::CalculateMultiplicity();
    else{
      double F = 0.5*(1. - f), Weight;
      Nch = 0.;

      for(int n = 0; n < (int)(Npart[0] + Npart[1]); n++){
        Weight = Rng.NegativeBinomial(Npp, k);
        Nch += F*Weight;
      }

      for(int n = 0; n < (int)Ncoll; n++){
        Weight = Rng.NegativeBinomial(Npp, k);
        Nch += f*Weight;
      }
    }
  }

  void GlauberMC :: CalculateEccentricity(int n){
    double F = 0.5*(1. - f);
    double NumC = 0., NumS = 0., Denom = 0.;
    double rT, phi, c;

    for(int i = 0; i < (int)Npart[0]; i++){
      rT = PartA[i].TransverseNorm(), phi = PartA[i].Phi();
      c = F*pow(rT, n);
      Denom += c, NumC += c*cos(n*phi), NumS += c*sin(n*phi);
    }

    for(int i = 0; i < (int)Npart[1]; i++){
      rT = PartB[i].TransverseNorm(), phi = PartB[i].Phi();
      c = F*pow(rT, n);
      Denom += c, NumC += c*cos(n*phi), NumS += c*sin(n*phi);
    }

    for(int i = 0; i < (int)Ncoll; i++){
      rT = Coll[i].TransverseNorm(), phi = Coll[i].Phi();
      c = f*pow(rT, n);
      Denom += c, NumC += c*cos(n*phi), NumS += c*sin(n*phi);
    }

    Ecc[n - 2] = sqrt(NumC*NumC + NumS*NumS)/(Denom + 1e-10);
    Psi[n - 2] = (1./n)*atan2(NumS, NumC);
  }

  double GlauberMC :: GetTransverseProfile(Position R){
    double F = 0.5*(1. - f), Weight = 1.;
    double ThP = F/Tpp, ThC = f/Tpp;

     // Participant contribution
    double ePart = 0.;
    for(int i = 0; i < (int)Npart[0]; i++){
      if(GammaFlag) Weight = Rng.Gamma(1., ThP);
      ePart += Weight*Thickness(R - PartA[i]);
    }

    for(int i = 0; i < (int)Npart[1]; i++){
      if(GammaFlag) Weight = Rng.Gamma(1., ThP);
      ePart += Weight*Thickness(R - PartB[i]);
    }

    // Binary collisions contribution
    double eColl = 0.;
    for(int i = 0; i < (int)Ncoll; i++){
      if(GammaFlag) Weight = Rng.Gamma(1., ThC);
      eColl += Weight*Thickness(R - Coll[i]);
    }

    return F*ePart + f*eColl;
  }

  // Definition from arXiv:1409.8164 [nucl-th]
  double GlauberMC :: Thickness(Position R){
    double Sg = SgNN/(8*Pi);
    double rT = pow(R.TransverseNorm(), 2);    // Transverse distance only
    return (1./(2*Pi*Sg))*exp(-0.5*rT/Sg);
  }

  void GlauberMC :: SetOutputFormat(int Nev){
    Glauber::SetOutputFormat(Nev);
    NucleusA.SetOutputFormat(Output, Nev, "Nucleus/NucleusA");
    NucleusB.SetOutputFormat(Output, Nev, "Nucleus/NucleusB");
  }

}
