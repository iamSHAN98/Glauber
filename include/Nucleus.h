#pragma once
#include <vector>
#include "Random.h"
#include "Profile.h"
#include "DataStream.h"

namespace EventGen{

  struct Nucleon{
    bool Participant;
    bool Neutron;                     // Neutron : True, Proton : False
    Position R;
  };

  class Nucleus{

    // Inputs
    NucleusData Data;

    // Output
    std::vector<Nucleon> Arr;

    // Utilities
    ChargeProfile *Prof;
    Utility::Random *Rng;

    // Constants
    static constexpr double R0 = 0.5; // Hard sphere radius for nucleons (fm)
    static constexpr double M = 10.;  // Constant : a) Normalization for uniform probability
                                      //            b) Spatial range for sampling
    static const int NMax = 100;      // Max number of interation for sampling each nucleon

    // Operations
    void Reset();
    bool SamplePosition(Position&);
    bool CheckOverlap(int, Position);

    public :

      Nucleus() = default;
      Nucleus(NucleusData N){ Initialize(N); }
      Nucleus(ChargeProfile *F){ Initialize(F); }
      ~Nucleus() = default;

      void Initialize(NucleusData N){
        Data = N;
        Prof->Initialize(Data);
        Arr.resize(Data.A);
      }

      void Initialize(ChargeProfile *F){
        Prof = F;
        Data = Prof->N;
        Arr.resize(Data.A);
      }

      // Setter
      void SetRandomGenerator(Utility::Random *Gen){ Rng = Gen; }
      void TagAsParticipant(int i){ Arr[i].Participant = true; }

      // Getter
      Position GetMassCenter();
      Position GetPosition(int i){ return Arr[i].R; }
      bool IsParticipant(int i){ return Arr[i].Participant; }
      bool IsNeutron(int i){ return Arr[i].Neutron; }

      // Operations
      void GenerateNucleons();
      void ShiftNucleus();                  // Shifts to origin
      void ShiftNucleus(Position);
      void RotateNucleus(double, double);

      // Output
      void SetOutputFormat(DataStream::File&, int, std::string);
  };

}
