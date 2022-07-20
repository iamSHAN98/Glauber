#include <chrono>
#include "Config.h"
#include "GlauberMC.h"

int main(int argc, char** argv){

  // Read configuration
  auto ConfigFile = EventGen::Config(argv[1]);
  
  auto Nuclei = ConfigFile.GetCollidingNuclei();
  auto Range = ConfigFile.GetImpactRange();
  auto sNN  = ConfigFile.GetCollisionEnergy();
  auto Param = ConfigFile.GetNode("Parameter");
  int Nev = ConfigFile.GetEventNumber();
  int Nf = ConfigFile.GetStatusStep();

  // Configure event generator
  EventGen::Glauber *Obj;

  if(ConfigFile.CheckMonteCarlo())
    Obj = new EventGen::GlauberMC(Nuclei[0], Nuclei[1]);
  else
    Obj = new EventGen::Glauber(Nuclei[0], Nuclei[1]);
  Obj->SetImpactRange(Range[0], Range[1]);
  Obj->SetCollisionEnergy(sNN);
  Obj->Configure(Param);
  Obj->SetOutputFormat(Nev);

  // Put additional inputs from command line (using argv) here 

  // Print information
  std::cout << "Configuration done \n";
  std::cout << "Using " << Obj->GetInfo() << "\n";
  std::cout << "Starting event generation for " << Nev << " event(s) \n";
  auto Start = std::chrono::system_clock::now();

  int n = 0;
  while(n < Nev){
    bool EventFlag = Obj->GenerateEvent();

    if(EventFlag){

      /*
        Put additional flags here, as follows in the pseudo-code
        if (!condition) continue;
      */

      Obj->WriteOutput();

      // Do additional work here; Certainly, this requires the event output

      if((n + 1)%Nf == 0)
        std::cout << (n + 1) << " events generated \n";

      n += 1;
    }
  }

  auto End = std::chrono::system_clock::now();
  double Time = std::chrono::duration_cast<std::chrono::duration<double>>(End - Start).count();
  std::cout << "Time : " << Time << "\n";

  Obj->SaveOutput();
  delete Obj;
  return 0;
}
