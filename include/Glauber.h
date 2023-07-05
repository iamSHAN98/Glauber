#pragma once
#include <cmath>
#include <sstream>
#include "yaml-cpp/yaml.h"
#include "DataStream.h"
#include "Random.h"
#include "Profile.h"

namespace EventGen{

  class Glauber{

    private :

      std::vector<Utility::Quadrature> Arr;          // Gauss-Legendre quadratures
      double R0;                                     // Integration limit (fm)
      Position Rcm;                                  // Center of mass of sources

      // Helper functions
      void Thickness(Position, double& Ta, double& Tb);
      void EulerRotation(Position&, double, double);

    protected :

      // Model inputs (user)
      NucleusData Na, Nb;
      double sNN, f, Npp;
      double Bmin, Bmax;

      // Model inputs (sampled)
      double b, Theta[2], Phi[2];
      Position B;

      // Model inputs (computed)
      double SgNN;
      ChargeProfile Fa, Fb;

      // Model outputs
      double Npart[2], Spec[2], Ncoll, Nch;
      static const int Ord = 5;                      // Max harmonic order
      double Ecc[Ord], Psi[Ord];

      // Utilities
      Utility::Random Rng;                           // Random number generation

      // Output
      DataStream::File Output;
      std::stringstream Name;
      std::stringstream FileName;

      // Operations
      virtual void CalculateNPartNColl();
      virtual void CalculateMultiplicity();
      virtual void CalculateEccentricity(int);

      // Helper functions
      void ObtainImpactVector();
      void ObtainRotationAngle();
      double CrossSection(double);           // Inelastic pp cross section (in GeV, out fm^2)

    public :

      Glauber() = default;

      Glauber(std::string A, std::string B){
        Na = ListNucleus[MapNucleus[A]];
        Nb = ListNucleus[MapNucleus[B]];

        Name << "Optical Glauber, " << Na.Name << "+" << Na.Name;
        FileName << "opt_glauber_" << Na.Name << "_" << Na.Name << ".h5";

        Fa.Initialize(Na);
        Fb.Initialize(Nb);

        R0 = 2*std::max(Na.R, Nb.R);
        Arr.resize(Fa.Nq);
        Fa.Int.GetQuadrature(Arr.data());
      }

      virtual ~Glauber() = default;

      // Setter
      void SetImpactRange(double Min, double Max){
        Bmin = Min, Bmax = Max;
      }
      void SetCollisionEnergy(double SNN){
        sNN = SNN;
        SgNN = CrossSection(sNN);
        Name << " at ECoM = " << sNN << " GeV";
      }

      // Operations
      virtual void Configure(YAML::Node);
      virtual bool GenerateEvent();

      // Getter
      double GetImpactParameter(){ return b; }
      double GetTheta(int);
      double GetPhi(int);
      double GetNPart(int i = 0);
      double GetSpecNeutron(int i = 0);
      double GetNColl(){ return Ncoll; }
      double GetMultiplicity(){ return Nch; }
      double GetEccentricity(int);
      double GetPlaneAngle(int);
      virtual double GetTransverseProfile(Position);

      // Output
      std::string GetInfo(){ return Name.str(); }
      virtual void SetOutputFormat(int); 
      void WriteOutput(){ Output.Write(); }
      void SaveOutput(){ Output.Close(); }

  };

}
