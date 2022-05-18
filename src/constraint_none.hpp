#ifndef __CONSTRAINT_NONE_HPP__  
#define __CONSTRAINT_NONE_HPP__ 

#include "constraint.hpp"

namespace AJM
{

  class ConstraintNone : public Constraint
  {
    public: 

      ConstraintNone() { }
      Vec apply(const VertexHandle<Property>& vh, const Vec& f) override
      {
        return Vec(f);
      }
      Vec apply(const Vec& f) override
      {
        return Vec(f);
      }

  };

}
#endif