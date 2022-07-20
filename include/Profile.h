#pragma once
#include "Integration.h"
#include "NucleusData.h"
#include "Position.h"

namespace EventGen{

  class ChargeProfile{

    NucleusData N;
    double Rho;

    static const int Nq = 31;                      // Quadrature length
    Utility::Integration Int;                      // Integration
    void Normalize();

    friend class Glauber;
    friend class Nucleus;

    public :

      ChargeProfile() = default;
      
      ChargeProfile(NucleusData Data){
        Initialize(Data);
      }

      ~ChargeProfile() = default;

      void Initialize(NucleusData Data){
        N = Data;
        Int.Initialize(Nq);
        Normalize();
      }

      // Getter
      double GetProfile(double*);
      double GetProfile(Position);
      double GetProbability(double*);
      double GetProbability(Position);
      double GetNormalization(){ return Rho; }

      // Charge density profiles
      double Fermi(double*);
      double Deformed(double*);
      double Gaussian(double*);
      double Hulthen(double*);

  };

}
