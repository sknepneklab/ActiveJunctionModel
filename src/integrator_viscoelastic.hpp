/* ***************************************************************************
 *
 *  Copyright (C) 2017 University of Dundee
 *  All rights reserved. 
 *
 *  This file is part of AJM (Active Junction Model) program.
 *
 *  AJM is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  AJM is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * ****************************************************************************/

/*!
 * \file integrator_brownian.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 20-May-2019
 * \brief IntegratorBrownian class 
*/

#ifndef __INTEGRATOR_VISCOELASTIC_HPP__
#define __INTEGRATOR_VISCOELASTIC_HPP__

#include "integrator.hpp"

#include <chrono>
#include <utility>
#include <map>
#include <memory>


namespace AJM
{

  using ConstrainerType = unique_ptr<Constrainer>;

  class IntegratorViscoelastic : public Integrator 
  {

    public:

      IntegratorViscoelastic(System& sys, ForceCompute& fc) : Integrator{sys, fc},
                                                              _tau_l{1.0},
                                                              _tau_s{1.0}
                                                                    
      { 
        
      }

      void step() override;
      void set_params(const params_type& params) override 
      { 
        for (auto& p : params)
        {
          if (p.first == "tau_l")
            _tau_l = p.second;
          if (p.first == "tau_s")
            _tau_s = p.second;
        }
      };
      void set_type_params(const string& type, const params_type& params) override { }
      void set_external_force(const string& vtype, const Vec& f) override { }
      void set_radial_force(const string& vtype, double f) override  { }
      void set_flag(const string& flag) override {  }
      bool converged() override { return true; }

      
      
    private:

      double _tau_l;
      double _tau_s;  // spokes relaxation time scale
  };

}

#endif

