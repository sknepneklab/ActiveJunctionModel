#ifndef __CONSTRAINT_UX_HPP__  
#define __CONSTRAINT_UX_HPP__ 

#include "constraint.hpp"

namespace AJM
{

  class ConstraintUX : public Constraint
  {
    public: 

      ConstraintUX() { }
      Vec apply(const VertexHandle<Property>& vh, const Vec& f) override
      {
        const Vec fv =  vh->data().force;
        return Vec((f - fv).x, 0.0);
      }
      Vec apply(const Vec& f) override
      {
        return Vec(0.0, 0.0);
      }
  };

}
#endif