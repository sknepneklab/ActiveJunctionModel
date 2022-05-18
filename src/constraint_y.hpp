#ifndef __CONSTRAINT_Y_HPP__  
#define __CONSTRAINT_Y_HPP__ 

#include "constraint.hpp"

namespace AJM
{

  class ConstraintY : public Constraint
  {
    public: 

      ConstraintY() { }
      Vec apply(const VertexHandle<Property>& vh, const Vec& f) override
      {
        return Vec(0.0, f.y);
      }
      Vec apply(const Vec& f) override
      {
        return Vec(0.0, f.y);
      }
  };

}
#endif