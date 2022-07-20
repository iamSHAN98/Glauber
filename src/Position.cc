#include <iomanip>
#include <cmath>
#include "Position.h"

namespace EventGen{

  double Position :: Norm(){
    return sqrt(x*x + y*y + z*z);
  }

  double Position :: TransverseNorm(){
    return sqrt(x*x + y*y);
  }

  double Position :: Theta(){
    return acos(z/(Norm() + 1e-10));
  }

  double Position :: Phi(){
    return atan2(y, x);
  }

  double Position :: Dot(const Position& v){
    return x*v.x + y*v.y + z*v.z;
  }

  Position Position :: Cross(const Position& v){
    double a = y*v.z - z*v.y;
    double b = z*v.x - x*v.z;
    double c = x*v.y - y*v.x;
    return Position{a, b, c};
  }

  Position Position :: operator - () const{
    return Position{-x, -y, -z};
  }

  Position Position :: operator * (double a) const{
    return Position{a*x, a*y, a*z};
  }

  Position Position :: operator / (double a) const{
    if(a == 0)
      throw std::runtime_error("Zero division error");
    return Position{x/a, y/a, z/a};
  }

  Position Position :: operator + (const Position& u) const{
    return Position{x + u.x, y + u.y, z + u.z};
  }

  Position Position :: operator - (const Position& u) const{
    return Position{x - u.x, y - u.y, z - u.z};
  }

  Position& Position :: operator += (const Position& u){
    x += u.x, y += u.y, z += u.z;
    return *this;
  }

  Position& Position :: operator -= (const Position& u){
    x -= u.x, y -= u.y, z -= u.z;
    return *this;
  }

  Position& Position :: operator *= (double a){
    x *= a, y *= a, z *= a;
    return *this;
  }

  Position operator * (double a, const Position& u){
    return Position{a*u.x, a*u.y, a*u.z};
  }

  std::istream& operator >> (std::istream& In, Position& V) {
    In >> V.x >> V.y >> V.z;
    return In;
  }

  std::ostream& operator << (std::ostream& Out, const Position& V) {
    Out << std::left << std::setw(15) << V.x;
    Out << std::left << std::setw(15) << V.y;
    Out << std::left << std::setw(15) << V.z;
    return Out;
  }

}
