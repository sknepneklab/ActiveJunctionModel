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
 * \file force_line_tension.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 12-Jun-2019
 * \brief ForceLineTension class 
 */ 

#ifndef __FORCE_LINE_TENSION_HPP__
#define __FORCE_LINE_TENSION_HPP__

#include "force.hpp"


namespace AJM
{

  // Force on a vertex
  class ForceLineTension : public Force
  {
    public:
      // Orphaned model options
      // const param_type& gammaLen, const param_type& l0, const param_type& counter
      // _gammaLen(gammaLen), _l0(l0), _counter(counter)
                                                                                                                                
      ForceLineTension(System& sys) : Force(sys), 
                                      _double_boundary_constants{false} 
                                    { 
                                      _lambda.resize(_sys.cell_types().size(), 0.0);
                                    }
      virtual ~ForceLineTension() { }
        
      // computes force on vertex by a given edge
      Vec compute(const VertexHandle<Property>&, const HEHandle<Property>&) override;  
      
      // compute edge tension for handling 4-vertices (and moves that computation out of integrator)
      // takes care of correct type of force law that way
      double tension(const HEHandle<Property>&, double, double) override;

      // Stress calculation not implemented yet
      void stress(FaceHandle<Property>& fh) override 
      {
        //throw runtime_error("Stress calculation has not been implemented for line tension force.");
      }

      // We do not compute face force for this force type
      Vec compute_face_force(FaceHandle<Property> &fh) override
      {
        return Vec(0, 0);
      }

      // Energy calculation 
      double energy(const FaceHandle<Property>&) override;
      
      // set all paramters for a given type
      void set_params(const string& cell_type, const params_type& params) override
      {
        for (auto p : params)
          if (p.first != "lambda")
            throw runtime_error("Unknown parameter "+p.first+".");

        if (params.find("lambda") == params.end())
          throw runtime_error("LineTension force requires parameter k.");
       
        try 
        {
          if (_sys.cell_types().find(cell_type) == _sys.cell_types().end())
            throw runtime_error("Fore line tension: Cell type " + cell_type + " is not defined.");
          int ct = _sys.cell_types()[cell_type];
          if (ct + 1 > _lambda.size())
            _lambda.push_back(params.at("lambda"));
          else
            _lambda[ct] = params.at("lambda");
        } 
        catch(const exception& e)
        {
          cerr << "Problem with setting line tension force paramters. Exception: " << e.what() << endl;
          throw;
        }
      }

      // set all vector-valued paramters for a given type
      void set_vec_params(const string& cell_type, const vec_params_type& params) override { };

      void set_relative_paramters(const string& cell_type, const multi_params_type& params) override { };

      void set_flag(const string& flag) override 
      { 
        if (flag == "use_cell_type")  
        {
          _use_cell_type = true;
          cout << "Warning! Setting use_cell_type flag in line tension force will override all parameters read from the input JSON file (if the read_params flag in the read_input function was set to True)." << endl;
        }
        else
          throw runtime_error("Unknown flag : " + flag + ".");
      }

      void copy_type_param_to_cell() override
      {
        for (FaceHandle<Property> fh = _sys.mesh().faces().begin(); fh != _sys.mesh().faces().end(); fh++)
          fh->data().lambda  = _lambda[fh->data().face_type];
      }
       
    
    private:

      vector<double> _lambda;
      bool _double_boundary_constants;   // if true, lambda for cells bordering a hole are doubled
      
  };

  
}

#endif
