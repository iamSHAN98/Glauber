#pragma once
#include <array>
#include "yaml-cpp/yaml.h"
#include "NucleusData.h"

namespace EventGen{

  class Config{

    YAML::Node ConfigFile;

    public :

      Config() = default;
      Config(std::string File){ LoadFile(File); }
      ~Config() = default;

      void LoadFile(std::string File){
        ConfigFile = YAML::LoadFile(File);
      }

      YAML::Node GetNode(std::string);      // Primary nodes

      // Model independent inputs
      bool CheckMonteCarlo();
      std::array<std::string, 2> GetCollidingNuclei();
      std::array<double, 2> GetImpactRange();
      double GetCollisionEnergy();
      int GetEventNumber();
      int GetStatusStep();

  };
}
