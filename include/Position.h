#pragma once
#include <iostream>

namespace EventGen{

  struct Position{
    double x, y, z;

    double Norm();
    double TransverseNorm();
    double Theta();
    double Phi();
    double Dot(const Position&);
    Position Cross(const Position&);

    Position operator - () const;
    Position operator * (double) const;
    Position operator / (double) const;
    Position operator + (const Position&) const;
    Position operator - (const Position&) const;
    Position& operator += (const Position&);
    Position& operator -= (const Position&);
    Position& operator *= (double);
  };

  Position operator * (double, const Position&);
  std::istream& operator >> (std::istream&, Position&);
  std::ostream& operator << (std::ostream&, const Position&);

}
