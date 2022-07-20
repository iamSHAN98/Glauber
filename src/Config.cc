#include <cmath>
#include "Config.h"

namespace EventGen{

  YAML::Node Config :: GetNode(std::string Key){
    return ConfigFile[Key];
  }

  bool Config :: CheckMonteCarlo(){
    return ConfigFile["Monte Carlo"].as<bool>();
  }

  std::array<std::string, 2> Config :: GetCollidingNuclei(){
    std::array<std::string, 2> List;
    List[0] = ConfigFile["Nucleus"]["A"].as<std::string>();
    List[1] = ConfigFile["Nucleus"]["B"].as<std::string>();
    return List;
  }

  std::array<double, 2> Config :: GetImpactRange(){
    std::array<double, 2> Range;
    Range[0] = ConfigFile["b"]["Min"].as<double>();
    Range[1] = ConfigFile["b"]["Max"].as<double>();
    return Range;
  }

  double Config :: GetCollisionEnergy(){
    double Val = ConfigFile["sNN"]["Value"].as<double>();
    auto Unit = ConfigFile["sNN"]["Unit"].as<std::string>();
    return Val*pow(10., MapUnit[Unit]);  // GeV
  }

  int Config :: GetEventNumber(){
    return ConfigFile["Events"].as<int>();
  }

  int Config:: GetStatusStep(){
    return ConfigFile["Status"].as<int>();
  }
  
}
