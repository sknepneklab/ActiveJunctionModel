#ifndef __CONSTRAINT_RADIAL_HPP__  
#define __CONSTRAINT_RADIAL_HPP__ 

#include "constraint.hpp"

namespace AJM
{

  class ConstraintRadial : public Constraint
  {
    public: 

      ConstraintRadial() { }
      Vec apply(const VertexHandle<Property>& vh, const Vec& f) override
      {
        const Vec r = vh->r.unit();
        double f_dot_r = dot(f, r);
        return Vec(f_dot_r * r.x, f_dot_r * r.y);
      }
      Vec apply(const Vec& f) override
      {
        return Vec(0.0, 0.0);
      }
  };

}
#endif