#pragma once
#include <string>
#include <map>

namespace EventGen{

  // Nuclear charge distributions
  enum KeyChargeProfile{Fermi = 0, Deformed, Gaussian, Hulthen};

  // Nucleus
  enum KeyNucleus{Au = 0, Au2, U, Pb, S, Cu, d};
  static std::map<std::string, KeyNucleus> MapNucleus{
    {"Au", KeyNucleus::Au},
    {"Au2", KeyNucleus::Au2},
    {"U", KeyNucleus::U},
    {"Pb", KeyNucleus::Pb},
    {"S", KeyNucleus::S},
    {"Cu", KeyNucleus::Cu},
    {"d", KeyNucleus::d}
  };

  // Energy units
  enum KeyUnit{MeV = -3, GeV = 0, TeV = 3};
  static std::map<std::string, KeyUnit> MapUnit{
    {"MeV", KeyUnit::MeV},
    {"GeV", KeyUnit::GeV},
    {"TeV", KeyUnit::TeV}
  };

}
