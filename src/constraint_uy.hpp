#ifndef __CONSTRAINT_UY_HPP__  
#define __CONSTRAINT_UY_HPP__ 

#include "constraint.hpp"

namespace AJM
{

  class ConstraintUY : public Constraint
  {
    public: 

      ConstraintUY() { }
      Vec apply(const VertexHandle<Property>& vh, const Vec& f) override
      {
        const Vec fv = vh->data().force;
        return Vec(0.0, (f - fv).y);
      }
      Vec apply(const Vec& f) override
      {
        return Vec(0.0, 0.0);
      }
  };

}
#endif