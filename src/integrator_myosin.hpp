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
 * \file integrator_myosin.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 20-May-2019
 * \brief IntegratorMyosin class 
*/

#ifndef __INTEGRATOR_MYOSIN_HPP__
#define __INTEGRATOR_MYOSIN_HPP__


#include "integrator.hpp"

#include <chrono>
#include <utility>
#include <map>
#include <memory>

using namespace std::chrono;
using std::map;


namespace AJM
{

  
  class IntegratorMyosin : public Integrator 
  {

    public:

      IntegratorMyosin(System& sys, ForceCompute& fc, int seed) : Integrator(sys, fc),
                                                                  _rng(RNG((seed >= 0) ? seed : system_clock::now().time_since_epoch().count())),
                                                                  _conserve_myosin{false},
                                                                  _use_differential_tension{false},
                                                                  _use_myo_pool{false},
                                                                  _myo_fluc{0.0},
                                                                  _myomax{1.0},
                                                                  _myo_cell{3.0},
                                                                  _k_on{0.1},
                                                                  _t_star{0.3},
                                                                  _k_min{0.1},
                                                                  _k_0{5.0},
                                                                  _J{0.0},
                                                                  _Q{1.0},
                                                                  _myo_diff{0.0},
                                                                  _myo_diff_vertex{0.0}
      { 
        
      }

      void step() override;
      void set_params(const params_type& params) override 
      { 
        for (auto p : params)
        {
          if (p.first == "kon")
            _k_on = p.second;
          else if (p.first == "kmin")
            _k_min = p.second;
          else if (p.first == "k0")
            _k_0 = p.second;
          else if (p.first == "J")
            _J = p.second;
          else if (p.first == "t*")
            _t_star = p.second;
          else if (p.first == "myosin_diffusion")
            _myo_diff = p.second;
          else if (p.first == "myosin_vertex_diffusion")
            _myo_diff_vertex = p.second;
          else if (p.first == "Q")
            _Q = p.second;
          else if (p.first == "cell_myosin")
            _myo_cell = p.second;
          else if (p.first == "myosin_fluctuation")
            _myo_fluc = p.second;
          else if (p.first == "max_myosin")
            _myomax = p.second;
          else if (p.first == "r")
            _k_on = p.second;
          else if (p.first == "p")
            _k_min = p.second;
          else
            throw runtime_error("Myosin integrator: Unknown parameter " + p.first + ".");
        } 
        
      };
      void set_type_params(const string& type, const params_type& params) override 
      { 
        if (_sys.cell_types().find(type) == _sys.cell_types().end())
          throw runtime_error("Cell type " + type + " is not defined.");
        if (params.find("cell_myosin") == params.end())
          runtime_error("Myosin integrator: parameter cell_myosin is required for cell type "+type+".");
        try
        {
          for (FaceHandle<Property> fh = _sys.mesh().faces().begin(); fh != _sys.mesh().faces().end(); fh++)
            if (fh->data().face_type == _sys.cell_types()[type])
              fh->data().cell_myo = params.at("cell_myosin");
        }
        catch(const exception& e)
        {
          cerr << "Problem setting type paramters for myosin integrator. Exception : " << e.what() << '\n';
        }
      }
      void set_external_force(const string& vtype, const Vec& f) override { }
      void set_radial_force(const string& vtype, double f) override  { }
      void set_flag(const string& flag) override 
      {  
        if (flag == "conserve_myosin")
          _conserve_myosin = true;
        else if (flag == "differential_tension")
          _use_differential_tension = true;
        else if (flag == "fixed_myosin")
          _use_myo_pool = true;
        else
          throw runtime_error("Myosin integrator: Unknown flag : " + flag + ".");

        if (_conserve_myosin && _use_differential_tension)
          throw runtime_error("One cannot simultaneously use conserved myosin and differential tension.");
      }
      bool converged() override { return true; }

      
      
    private:

      RNG  _rng;
      bool _conserve_myosin;
      bool _use_differential_tension;
      bool _use_myo_pool;  // If true, assume that there is a fixed pool of myosin in each cell
      double _k_on;
      double _k_min;
      double _k_0;
      double _J;
      double _t_star;
      double _myo_diff;
      double _myo_diff_vertex;
      double _Q;
      double _myo_cell;
      double _myomax;
      double _myo_fluc;

      // Auxiliary private functions
      double _myo_dyn(HEHandle<Property>&);
      double _myo_dyn(VertexHandle<Property>&, HEHandle<Property>&);


      
  };

}

#endif

