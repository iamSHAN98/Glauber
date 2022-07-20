#pragma once
#include <random>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

namespace Utility{

  class Random{

    private :

      std::random_device seed;
      gsl_rng *gen;

    public :

      Random();
      ~Random();

      template <typename T>
      void Shuffle(T *Arr, int Size){
        gsl_ran_shuffle(gen, Arr, Size, sizeof(T));
      }

      double Uniform(double, double);
      int Uniform(int, int);

      double Linear(double, double);
      double Sine();

      int NegativeBinomial(double, double);
      double Gamma(double, double);

      int Poisson(double);
      double Exponential(double);
      double Gaussian(double, double);

  };
}
