#ifndef __CONSTRAINT_X_HPP__  
#define __CONSTRAINT_X_HPP__ 

#include "constraint.hpp"

namespace AJM
{

  class ConstraintX : public Constraint
  {
    public: 

      ConstraintX() { }
      Vec apply(const VertexHandle<Property>& vh, const Vec& f) override
      {
        return Vec(f.x, 0.0);
      }
      Vec apply(const Vec& f) override
      {
        return Vec(f.x, 0.0);
      }
  };

}
#endif