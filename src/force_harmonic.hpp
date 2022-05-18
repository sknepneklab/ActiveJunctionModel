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
 * \file force_harmonic.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 19-May-2019
 * \brief ForceHarmonic class 
 */ 

#ifndef __FORCE_HARMONIC_HPP__
#define __FORCE_HARMONIC_HPP__

#include "force.hpp"


namespace AJM
{

  // Force on a vertex
  class ForceHarmonic : public Force
  {
    public:
      // Orphaned model options
      // const param_type& gammaLen, const param_type& l0, const param_type& counter
      // _gammaLen(gammaLen), _l0(l0), _counter(counter)
                                                                                                                                
      ForceHarmonic(System& sys) : Force(sys) 
                                    { 
                                      _k.resize(_sys.cell_types().size(), 0.0);
                                    }
      virtual ~ForceHarmonic() { }
        
      // computes force on vertex by a given edge
      Vec compute(const VertexHandle<Property>&, const HEHandle<Property>&) override;  
      
      // compute edge tension for handling 4-vertices (and moves that computation out of integrator)
      // takes care of correct type of force law that way
      double tension(const HEHandle<Property>&, double, double) override;

      // Stress calculation not implemented yet
      void stress(FaceHandle<Property>& fh) override 
      {
        //throw runtime_error("Stress calculation has not been implemented for harmonic force.");
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
          if (p.first != "k")
            throw runtime_error("Unknown parameter "+p.first+".");

        if (params.find("k") == params.end())
          throw runtime_error("Harmonic force requires parameter k.");
       
        try 
        {
          if (_sys.cell_types().find(cell_type) == _sys.cell_types().end())
            throw runtime_error("Fore harmonic: Cell type " + cell_type + " is not defined.");
          int ct = _sys.cell_types()[cell_type];
          if (ct + 1 > _k.size())
            _k.push_back(params.at("k"));
          else
            _k[ct] = params.at("k");
        } 
        catch(const exception& e)
        {
          cerr << "Problem with setting harmonic force paramters. Exception: " << e.what() << endl;
          throw;
        }
      }

      // set all vector-valued paramters for a given type
      void set_vec_params(const string& cell_type, const vec_params_type& params) override { };

      void set_relative_paramters(const string& cell_type, const multi_params_type& params) override
      {
        if (_use_cell_type)
          throw runtime_error("Using the \"use_cell_type\" flag is not compatible with setting relative values of per-cell paramters.");
        for (auto p : params)
          if (p.first != "k" && p.first != "seed")
            throw runtime_error("Unknown parameter " + p.first + ".");
        
        bool has_k = params.find("k") != params.end();
        
        try 
        {
          unsigned int seed = system_clock::now().time_since_epoch().count();
          if (params.find("seed") != params.end())
            seed = static_cast<unsigned int>(params.at("seed")[0]);
          RNG rng(seed);
          bool setall = false; 
          if (cell_type == "all")
            setall = true;
          else
          {
            if (_sys.cell_types().find(cell_type) == _sys.cell_types().end())
              throw runtime_error("Fore area: Cell type " + cell_type + " is not defined.");
          }
          
          for (FaceHandle<Property> fh = _sys.mesh().faces().begin(); fh != _sys.mesh().faces().end(); fh++)
          {
            if (has_k)
            {
              double mean_k = params.at("k")[0];
              double std_k = 0.0;
              if (params.at("k").size() > 1) std_k = params.at("k")[1];
              double scale_k = (fabs(std_k) < 1e-6) ? mean_k : rng.normal(mean_k, std_k);
              if (setall)
                fh->data().k *= scale_k;
              else
                if (fh->data().type_name == cell_type)
                  fh->data().k *= scale_k;
            }
          }
        } 
        catch(const exception& e)
        {
          cerr << "Problem with setting relative paramters for harmonic force. Exception: " << e.what() << '\n';
          throw;
        }  

      }

      void set_flag(const string& flag) override 
      { 
        if (flag == "use_cell_type") 
        { 
          _use_cell_type = true;
          cout << "Warning! Setting use_cell_type flag in harmonic force will override all parameters read from the input JSON file (if the read_params flag in the read_input function was set to True)." << endl;
        }
        else
          throw runtime_error("Unknown flag : " + flag + ".");
      }
       
      void copy_type_param_to_cell() override
      {
        for (FaceHandle<Property> fh = _sys.mesh().faces().begin(); fh != _sys.mesh().faces().end(); fh++)
          fh->data().k  = _k[fh->data().face_type];
      }
      
    
    private:

      vector<double> _k;
      
  };

  
}

#endif
