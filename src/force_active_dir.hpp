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
 * \file force_active_dir.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 14-Jul-2019
 * \brief ForceActiveDir class 
 */ 

#ifndef __FORCE_ACTIVE_DIR_HPP__
#define __FORCE_ACTIVE_DIR_HPP__

#include "force.hpp"


namespace AJM
{

  // Force on a vertex
  class ForceActiveDir : public Force
  {
    public:
      // Orphaned model options
      // const param_type& gammaLen, const param_type& l0, const param_type& counter
      // _gammaLen(gammaLen), _l0(l0), _counter(counter)
                                                                                                                                
      ForceActiveDir(System& sys) : Force(sys) 
                                    { 
                                      _beta.resize(_sys.cell_types().size(), 0.0);
                                      _alpha.resize(_sys.cell_types().size(), 1.0);
                                      _psi.resize(_sys.cell_types().size(), 2);
                                    }
      virtual ~ForceActiveDir() { }
        
      // computes force on vertex by a given edge
      Vec compute(const VertexHandle<Property>&, const HEHandle<Property>&) override;  
      
      // compute edge tension for handling 4-vertices (and moves that computation out of integrator)
      // takes care of correct type of force law that way
      double tension(const HEHandle<Property>&, double, double) override;

      // Stress calculation not implemented yet
      void stress(FaceHandle<Property>& fh) override 
      {
        //throw runtime_error("Stress calculation has not been implemented for active dir force.");
      }

      // We do not compute face force for this force type
      Vec compute_face_force(FaceHandle<Property> &fh) override
      {
        return Vec(0, 0);
      }

      // Energy of Active Dir force is currently not defined, so we return 0
      double energy(const FaceHandle<Property>& fh) override
      {
        return 0.0;   
      }
      
      // set all paramters for a given type
      void set_params(const string& cell_type, const params_type& params) override
      {
        for (auto p : params)
          if (p.first != "alpha" && p.first != "beta" && p.first != "psi")
            throw runtime_error("Unknown parameter "+p.first+".");

        if (params.find("beta") == params.end())
          throw runtime_error("Active force requires parameter beta.");
        if (params.find("alpha") == params.end())
          throw runtime_error("Active force requires parameter alpha.");
        if (params.find("psi") == params.end())
          throw runtime_error("Active force requires parameter psi.");

        try 
        {
          if (_sys.cell_types().find(cell_type) == _sys.cell_types().end())
            throw runtime_error("Fore active direction: Cell type " + cell_type + " is not defined.");
          int ct = _sys.cell_types()[cell_type];
          if (ct + 1 > _beta.size())
            _beta.push_back(params.at("beta"));
          else
            _beta[ct] = params.at("beta");
          if (ct + 1 > _alpha.size())
            _alpha.push_back(params.at("alpha"));
          else
            _alpha[ct] = params.at("alpha");
          if (ct + 1 > _psi.size())
            _psi.push_back(static_cast<int>(params.at("psi")));
          else
            _psi[ct] = static_cast<int>(params.at("psi"));
        } 
        catch(const exception& e)
        {
          cerr << "Problem with setting active force paramters. Exception: " << e.what() << endl;
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
        cout << "Warning! Setting use_cell_type flag in active dir will override all parameters read from the input JSON file (if the read_params flag in the read_input function was set to True)." << endl;
        }
        else 
        throw runtime_error("Unknown flag " + flag + ".");
      }

      void copy_type_param_to_cell() override
      {
        for (FaceHandle<Property> fh = _sys.mesh().faces().begin(); fh != _sys.mesh().faces().end(); fh++)
        {
          fh->data().alpha = _alpha[fh->data().face_type];
          fh->data().beta  = _beta[fh->data().face_type];
          fh->data().psi = _psi[fh->data().face_type];
        }
      }

    
    private:

      vector<double> _beta;
      vector<double> _alpha;
      vector<int> _psi;
      
  };

  
}

#endif
