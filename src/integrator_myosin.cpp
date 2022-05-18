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
 * \file integrator_myosin.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 20-May-2019
 * \brief IntegratorMyosin class 
*/

#include "integrator_myosin.hpp"

namespace AJM
{
  void IntegratorMyosin::step()
  {
    double sqrt_dt = sqrt(_dt);
        
    for (FaceHandle<Property> fh = _sys.mesh().faces().begin(); fh != _sys.mesh().faces().end(); fh++)
    {
      double alpha = fh->data().alpha;
      // update average tension per cell
      if (_use_differential_tension || _use_myo_pool)
      {
        HEHandle<Property> first = fh->he();
        HEHandle<Property> he = fh->he();
        double T = 0.0;
        double active_myo = 0.0;
        do
        {
          if (_use_differential_tension)
            T += he->data().tension;
          if (_use_myo_pool)
            active_myo += he->data().myo;
          he = he->next();
        } while (he != first);
        if (_use_differential_tension)
          fh->data().avg_tension = T/fh->nsides;
        if (_use_myo_pool)
          fh->data().active_myo = active_myo; 
      }

      // no myosin conservation if we don't set it or if we are on the outer face
      if (!_conserve_myosin || fh->outer) 
      {
        HEHandle<Property> first = fh->he();
        HEHandle<Property> he_run = first;
        do
        {
          if (!he_run->edge()->erased)
            he_run->data().myo += _dt*this->_myo_dyn(he_run) + sqrt_dt*_myo_fluc*_rng.gauss_rng();		  
          he_run = he_run->next();
        } while (he_run != first);
      }
      else
      {
        // Myosin dynamics including conservation law
        // This can be done using a Lagrange multiplier (see notes)
        // lambda = -1/z (\sum _myo_dyn - Q*kmin (M-M0))
        // where Q is akin the thermal coupling parameter that you get in numerical thermostats
        // a myostat, so to say ...
        HEHandle<Property> first = fh->he(); //select half-edge	that belongs to the face fh
        HEHandle<Property> he_run = first;
        double fconst = 0.0;                 //sum_j dot(m_j) + sum_v dot(m_v) 
        int ze = 0;                          //number of edges + 4-vertices
        double myotot = 0.0;
        
        do 
        {
          double le = (he_run->from()->r - he_run->to()->r).len();
          fconst += this->_myo_dyn(he_run)*(1 + alpha*le);
          myotot += he_run->data().myo*(1 + alpha*le);
          if (he_run->from()->coordination > 3) // if collapsed junction vertex
          {
            VertexHandle<Property> vh = he_run->from();
            myotot += _sys.get_myostore(vh).get_myo(fh->id);
            fconst += _myo_dyn(vh,he_run); // myosin update rule for collapsed junctions 
            ze++; //  count vertex
          }
          he_run = he_run->next();
          ze++; // count edges
        } while (he_run != first);
        fconst += _Q*_k_on*(myotot - fh->data().cell_myo*(1 + alpha));
        fconst /= ze + alpha*_sys.mesh().perim(fh);

        do
        {
          he_run->data().myo += _dt*(this->_myo_dyn(he_run) - fconst) + sqrt_dt*_myo_fluc*_rng.gauss_rng(); 

          if (he_run->data().myo < 0) he_run->data().myo = 0.0;

          if (he_run->from()->coordination > 3) // if collapsed junction vertex
          {
            VertexHandle <Property> vh = he_run->from();
            double mv = _sys.get_myostore(vh).get_myo(he_run->face()->id);
            // unidirectional diffusion to the two half-edges
            _sys.get_myostore(vh).get_myo(he_run->face()->id) += _dt*_myo_dyn(vh,he_run); 
            he_run->data().myo += _dt*_myo_diff_vertex*mv;
            he_run->prev()->data().myo += _dt*_myo_diff_vertex*mv;
            _sys.get_myostore(vh).get_l0() -= _dt*0.5*_k_on*_Q*_sys.get_myostore(vh).get_l0(); // slow decay to 0            
          }
          he_run = he_run->next();
        } while (he_run != first );
      }
    }
  }

  double IntegratorMyosin::_myo_dyn(HEHandle<Property>& he)
  {
    double T = he->data().tension;
    double inhibition = 0.0;

    double k_on = _k_on;
    if (_J > 0)
    {
      if (_sys.symmetric_myosin())
        inhibition = he->prev()->data().myo + he->next()->data().myo - 1.0;
      else
        inhibition = he->prev()->data().myo + he->next()->data().myo;  
      k_on *= exp(-_J*inhibition);
    }

    double add_diff = 0.0;
    if (_myo_diff > 0)
      add_diff = _myo_diff*(he->prev()->data().myo - 2*he->data().myo + he->next()->data().myo);
    
    if (_use_differential_tension)
    {
      double dT = T - he->face()->data().avg_tension;
      return (k_on - _k_min*he->data().myo*(1.0/_myomax + exp(-_k_0*(dT - _t_star))) + add_diff);
    }
    else if (_use_myo_pool)
    {
      FaceHandle<Property> fh = he->face();
      double M = fh->data().cell_myo != 0.0 ? fh->data().cell_myo : _myo_cell;
      double m = fh->data().active_myo;
      double z = static_cast<double>(fh->nsides);
      double r = k_on; 
      double p = _k_min;
      if (fh->outer)
        return 0.0;
      else
        return r*((M-m) - p*z*he->data().myo*(exp(-_k_0*(T - _t_star)) + 1.0/_myomax));
    }
    else
      return (k_on - _k_min*he->data().myo*(1.0/_myomax + exp(-_k_0*(T - _t_star))) + add_diff);
  }

  
  double IntegratorMyosin::_myo_dyn(VertexHandle<Property>& vh, HEHandle<Property>& he)
  { 
    return -2.0*_myo_diff_vertex*_sys.get_myostore(vh).get_myo(he->face()->id);
  }  
}
