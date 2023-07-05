#pragma once
#include "Glauber.h"
#include "Nucleus.h"

namespace EventGen{

  class GlauberMC : public Glauber{

    Nucleus NucleusA, NucleusB;

    // Source positions
    std::vector<Position> PartA, PartB;
    std::vector<Position> Coll;

    // Flags
    bool NegBinFlag, GammaFlag;

    // Model inputs (user)
    double k, Tpp;

    // Operations
    void CalculateNPartNColl();
    void CalculateMultiplicity();
    void CalculateEccentricity(int);
    void ShiftMassCenter();

    // Helper functions
    double Thickness(Position);

    public :

      GlauberMC() : Glauber(){}

      GlauberMC(std::string A, std::string B) : Glauber(A, B){
        Name.str("");
        Name << "Monte Carlo Glauber, " << Na.Name << "+" << Na.Name;

        FileName.str("");
        FileName << "mc_glauber_" << Na.Name << "_" << Na.Name << ".h5";

        NucleusA.Initialize(&Fa);
        NucleusA.SetRandomGenerator(&Rng);

        NucleusB.Initialize(&Fb);
        NucleusB.SetRandomGenerator(&Rng);

        // Initial estimate for Npart & Ncoll
        PartA.resize(Na.A), PartB.resize(Nb.A), Coll.resize(5000);
      }

      ~GlauberMC() = default;

      // Operations
      void Configure(YAML::Node);
      bool GenerateEvent();

      // Getter
      double GetTransverseProfile(Position);

      // Output
      void SetOutputFormat(int Nev);

  };

}
