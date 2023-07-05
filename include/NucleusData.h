#pragma once
#include "KeyMap.h"

namespace EventGen{

  struct NucleusData{
    char Name[10];
    uint A, Z;
    double R, a, w, b2, b4;
    KeyChargeProfile Profile;
  };

  // Values from arXiv:1408.2549 [nucl-ex]
  const NucleusData ListNucleus[] = {
    {"Au", 197, 79, 6.38, 0.535, 0., 0., 0., KeyChargeProfile::Fermi},
    {"Au2", 197, 79, 6.38, 0.535, 0., -0.131, -0.031, KeyChargeProfile::Deformed},
    {"U", 238, 92, 6.81, 0.44, 0., 0.28, 0.093, KeyChargeProfile::Deformed},
    {"Pb", 207, 82, 6.62, 0.546, 0., 0., 0., KeyChargeProfile::Fermi},
    {"S", 32, 16, 2.54, 2.191, 0.16, 0., 0., KeyChargeProfile::Gaussian},
    {"Cu", 63, 29, 4.2, 0.596, 0.0, 0., 0., KeyChargeProfile::Fermi},
    {"d", 2, 1, 0.01, 0.228, 0., 1.177, 0., KeyChargeProfile::Hulthen},
    {"He", 4, 2, 0.964, 0.322, 0.517, 0., 0., KeyChargeProfile::Fermi},
    {"O", 16, 8, 2.608, 0.513, -0.051, 0., 0., KeyChargeProfile::Fermi},
  };

}
