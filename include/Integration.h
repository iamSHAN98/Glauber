#pragma once
#include "gsl/gsl_integration.h"

/*
    Note : Quadratures are generated on the symmetric interval [-1, 1].
    Using parametrization x = 1/2 [ (b - a) t + (a + b) ], t in [-1, 1]
    summation is applied for arbitrary interval [a, b].
*/

namespace Utility{

  struct Quadrature{
    double Node, Weight;
  };

  class Integration{

    private :

      gsl_integration_fixed_workspace *w;

    public :

      Integration() = default;

      Integration(int N){
        Initialize(N);
      }

      ~Integration(){
        gsl_integration_fixed_free(w);
      }

      void Initialize(int N){
        auto type = gsl_integration_fixed_legendre;
        w = gsl_integration_fixed_alloc(type, N, -1., 1., 0., 0.); 
      }

      template <class T>
      double Integrate1D(T *Obj, double (T::*F)(double*), double *Min, double *Max){

        double M = Max[0] - Min[0], P = Max[0] + Min[0];
        double R[] = {0., Min[1], Min[2]};

        double Sum = 0.;
        for(int i = 0; i < w->n; i++){
          R[0] = 0.5*(M*w->x[i] + P);
          Sum += (*Obj.*F)(R)*w->weights[i];
        }

        return 0.5*M*Sum;
      }

      template <class T>
      double Integrate2D(T *Obj, double (T::*F)(double*), double *Min, double *Max){

        double M = Max[1] - Min[1], P = Max[1] + Min[1];
        double Ri[] = {Min[0], 0., Min[2]};
        double Rf[] = {Max[0], 0., Max[2]};

        double Sum = 0.;
        for(int i = 0; i < w->n; i++){
          Ri[1] = Rf[1] = 0.5*(M*w->x[i] + P);
          Sum += Integrate1D(Obj, F, Ri, Rf)*w->weights[i];
        }

        return 0.5*M*Sum;
      }

      template <class T>
      double Integrate3D(T *Obj, double (T::*F)(double*), double *Min, double *Max){

        double M = Max[2] - Min[2], P = Max[2] + Min[2];
        double Ri[] = {Min[0], Min[1], 0.};
        double Rf[] = {Max[0], Max[1], 0.};

        double Sum = 0.;
        for(int i = 0; i < w->n; i++){
          Ri[2] = Rf[2] = 0.5*(M*w->x[i] + P);
          Sum += Integrate2D(Obj, F, Ri, Rf)*w->weights[i];
        }

        return 0.5*M*Sum;
      }

      int GetQuadratureLength(){ return w->n; }

      void GetQuadrature(Quadrature *Arr){
        for(int i = 0; i < w->n; i++){
          Arr[i].Node = w->x[i];
          Arr[i].Weight = w->weights[i];
        }
      }

  };
}
