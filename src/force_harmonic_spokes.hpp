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
 * \file force_harmonic_spokes.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Aug-2019
 * \brief ForceHarmonicSpokes class 
 */ 

#ifndef __FORCE_HARMONIC_SPOKES_HPP__
#define __FORCE_HARMONIC_SPOKES_HPP__

#include "force.hpp"


namespace AJM
{

  // Force on a vertex
  class ForceHarmonicSpokes : public Force
  {
    public:
                                                                                                                              
      ForceHarmonicSpokes(System& sys) : Force(sys) 
                                      { 
                                        _k.resize(_sys.cell_types().size(), 0.0);
                                      }
      virtual ~ForceHarmonicSpokes() { }
        
      // computes force on vertex by a given edge
      Vec compute(const VertexHandle<Property>&, const HEHandle<Property>&) override;  
      
      // compute edge tension for handling 4-vertices (and moves that computation out of integrator)
      // takes care of correct type of force law that way
      double tension(const HEHandle<Property>&, double, double) override;

      // Stress calculation not implemented yet
      void stress(FaceHandle<Property>& fh) override 
      {
        //throw runtime_error("Stress calculation has not been implemented for harmonic spokes force.");
      }

      // We do not compute face force for this force type
      Vec compute_face_force(FaceHandle<Property> &fh) override
      {
        return Vec(0, 0);
      }

      // Energy of harmonic spokes is currently not defined, so we return 0
      double energy(const FaceHandle<Property>& fh) override
      {
        return 0.0;   
      }
      
      // set all paramters for a given type
      void set_params(const string& cell_type, const params_type& params) override
      {
        for (auto p : params)
          if (p.first != "k")
            throw runtime_error("Unknown parameter "+p.first+".");

        if (params.find("k") == params.end())
          throw runtime_error("Harmonic force on spokes requires parameter k (spoke stiffness).");
      
        try 
        {
          if (_sys.cell_types().find(cell_type) == _sys.cell_types().end())
            throw runtime_error("Fore harmonic spokes: Cell type " + cell_type + " is not defined.");
          int ct = _sys.cell_types()[cell_type];
          if (ct + 1 > _k.size())
            _k.push_back(params.at("k"));
          else
            _k[ct] = params.at("k");
        } 
        catch(const exception& e)
        {
          cerr << "Problem with setting harmonic spokes force paramters. Exception: " << e.what() << '\n';
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
        cout << "Warning! Setting use_cell_type flag in harmonic spokes force will override all parameters read from the input JSON file (if the read_params flag in the read_input function was set to True)." << endl;
        }
        else
        throw runtime_error("Uknown flag : "+ flag + ".");
      }

      void copy_type_param_to_cell() override
      {
        for (FaceHandle<Property> fh = _sys.mesh().faces().begin(); fh != _sys.mesh().faces().end(); fh++)
          fh->data().k_spoke = _k[fh->data().face_type];
      }

    
    private:

      vector<double> _k;  // stiffness of spokes (defined by the cell type)
            
  };

  
}

#endif
