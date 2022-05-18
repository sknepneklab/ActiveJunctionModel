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
 * \file integrator_myosin_diffusion.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Aug-2019
 * \brief IntegratorMyosinDiffusion class 
*/

#ifndef __INTEGRATOR_MYOSIN_DIFFUSION_HPP__
#define __INTEGRATOR_MYOSIN_DIFFUSION_HPP__


#include "integrator.hpp"

#include <chrono>
#include <utility>
#include <map>
#include <memory>

using namespace std::chrono;
using std::map;


namespace AJM
{

  
  class IntegratorMyosinDiffusion : public Integrator 
  {

    public:

      IntegratorMyosinDiffusion(System& sys, ForceCompute& fc, int seed) : Integrator{sys, fc},
                                                                           _D{1.0},
                                                                           _use_spokes{false}
      { 
        
      }

      void step() override;
      void set_params(const params_type& params) override 
      { 
        for (auto p : params)
        {
          if (p.first == "D")
            _D = p.second;
          else 
            throw runtime_error("Myosin diffusion integrator: Unknown parameter " + p.first + ".");
        } 
        
      };

      void set_type_params(const string& type, const params_type& params) override 
      { 
        if (_sys.cell_types().find(type) == _sys.cell_types().end())
          throw runtime_error("Cell type " + type + " is not defined.");
        if (params.find("D") == params.end())
          runtime_error("Myosin diffusion integrator: parameter D is required for cell type "+type+".");
        try
        {
          for (FaceHandle<Property> fh = _sys.mesh().faces().begin(); fh != _sys.mesh().faces().end(); fh++)
            if (fh->data().face_type == _sys.cell_types()[type])
              fh->data().D = params.at("D");
        }
        catch(const exception& e)
        {
          cerr << "Problem setting type paramters for myosin diffusion integrator. Exception : " << e.what() << '\n';
        }
      }
      
      void set_external_force(const string& vtype, const Vec& f) override { }
      void set_radial_force(const string& vtype, double f) override  { }
      void set_flag(const string& flag) override 
      {  
        if (flag == "use_spokes")
          _use_spokes = true;
        else
          throw runtime_error("Myosin diffussion integrator: Unknown flag : " + flag + ".");
      }
      bool converged() override { return true; }

      
      
    private:

      double _D; // Diffusion constant
      bool _use_spokes;  // If true, diffuse between junctions and spokes
  };

}

#endif

