#define Pi M_PI
#include "Nucleus.h"
using namespace Utility;

namespace EventGen{

  Position Nucleus :: GetMassCenter(){
    Position R = {};

    for(int i = 0; i < Data.A; i++)
      R += Arr[i].R;

    return R/Data.A;
  }

  void Nucleus :: ShiftNucleus(){
    Position R = GetMassCenter();
    ShiftNucleus(-R);
  }

  void Nucleus :: ShiftNucleus(Position dR){
    for(int i = 0; i < Data.A; i++)
      Arr[i].R += dR;
  }

  void Nucleus :: RotateNucleus(double th, double ph){
    Position R;
    double x, y, z;
    for(int i = 0; i < Data.A; i++){
      x = Arr[i].R.x, y = Arr[i].R.y, z = Arr[i].R.z;

      // Euler rotation
      R.x = x*cos(th)*cos(ph) + z*cos(ph)*sin(th) - y*sin(ph);
      R.y = y*cos(ph) + x*cos(th)*sin(ph) + z*sin(th)*sin(ph);
      R.z = z*cos(th) - x*sin(th);

      Arr[i].R = R;
    }
  }

  void Nucleus :: GenerateNucleons(){
    Reset();
    Position R;
    int Ntry = 0;

    int i = 0;
    while(i < Data.A){
      if(SamplePosition(R)){
        if(!CheckOverlap(i, R)){
          Arr[i++].R = R;
          Ntry = 0;
        }
        else Ntry += 1;
      }
      else if(Ntry == NMax){
        Arr[i++].R = R;
        Ntry = 0;
      }
      else Ntry += 1;
    }

    // Randomly tag a nucleon as a neutron
    int Index[Data.A];
    std::iota(Index, Index + Data.A, 0);
    Rng->Shuffle(Index, Data.A);

    for(int n = 0; n < Data.A - Data.Z; n++){
      i = Index[n];
      Arr[i].Neutron = true;
    }
  }

  bool Nucleus :: SamplePosition(Position& R){
    double u[4];
    for(int j = 0; j < 4; j++)
      u[j] = Rng->Uniform(0., 1.);

    double r = 3*M*u[0];
    double th = Pi*u[1];
    double phi = 2*Pi*u[2];

    double Ri[] = {r, th, phi};
    double Prob = Prof->GetProbability(Ri);

    // Rejection sampling
    if(M*u[3] <= Prob){
      R.x = r*sin(th)*cos(phi);
      R.y = r*sin(th)*sin(phi);
      R.z = r*cos(th);
      return true;
    }

    return false;
  }

  bool Nucleus :: CheckOverlap(int i, Position R){
    double dR;
    for(int j = 0; j < i; j++){
      dR = (Arr[j].R - R).Norm();
      if(dR < R0) return true;
    }
    return false;
  }

  void Nucleus :: Reset(){
    for(int i = 0; i < Data.A; i++){
      Arr[i].R = {};
      Arr[i].Neutron = false;
      Arr[i].Participant = false;
    }
  }

  void Nucleus :: SetOutputFormat(DataStream::File& Obj, int N, std::string Path){

    // MetaData
    DataStream::MetaData InfoR(DataStream::Compound, {1}, "Position");
    InfoR.AddMember(&Position::x, "x", DataStream::Double);
    InfoR.AddMember(&Position::y, "y", DataStream::Double);
    InfoR.AddMember(&Position::y, "z", DataStream::Double);

    DataStream::MetaData InfoN(DataStream::Compound, {Data.A}, "Nucleon");
    InfoN.AddMember(&Nucleon::Participant, "Participant", DataStream::Bool);
    InfoN.AddMember(&Nucleon::Neutron, "Neutron", DataStream::Bool);
    InfoN.AddMember(&Nucleon::R, InfoR);

    // Dataset
    Obj.Add(Path, Arr.data(), InfoN);
    
    Obj.Configure(N);
  }

}
