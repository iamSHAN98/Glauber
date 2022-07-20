#include <cmath>
#include "Random.h"

namespace Utility{

  Random :: Random(){
    gen = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(gen, seed());
  }

  Random :: ~Random(){
    gsl_rng_free(gen);
  }

  double Random :: Uniform(double Min, double Max){
    double r = (Max - Min)*gsl_rng_uniform(gen);
    return r + Min;
  }

  int Random :: Uniform(int Min, int Max){
    int r = gsl_rng_uniform_int(gen, Max - Min + 1);
    return r + Min;
  }

  int Random :: NegativeBinomial(double Mu, double K){
    double p = K/(Mu + K);      // Success probability
    return gsl_ran_negative_binomial(gen, p, K);
  }

  double Random :: Gamma(double Mu, double Th){
    double k = Mu/Th;
    return gsl_ran_gamma(gen, k, Th);
  }

  int Random :: Poisson(double Mu){
    return gsl_ran_poisson(gen, Mu);
  }

  double Random :: Exponential(double Mu){
    return gsl_ran_exponential(gen, Mu);
  }

  double Random :: Gaussian(double Mu, double Sg){
    double x = gsl_ran_gaussian(gen, Sg);
    return x + Mu;
  }

  // Inverse sampling
  double Random :: Linear(double Min, double Max){
    double r = gsl_rng_uniform(gen);
    return sqrt(r*(Max*Max - Min*Min) + Min*Min);
  }

  double Random :: Sine(){
    double r = gsl_rng_uniform(gen);
    return acos(2*r - 1.);
  }

}
