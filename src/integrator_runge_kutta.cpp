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
 * \file integrator_runge_kutta.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com 
 * \date 20-Aug-2020
 * \brief IntegratorRungeKutta class 
*/

#include "integrator_runge_kutta.hpp"

namespace AJM
{
  void IntegratorRungeKutta::step()
  {
    double mu = 1.0 / _gamma;    // mobility 
    // Compute force on each vertex
    for (VertexHandle<Property> vh = _sys.mesh().vertices().begin(); vh != _sys.mesh().vertices().end(); vh++)
      if (!vh->erased)
        _force_compute.compute(vh);

    // If we are computing stress, add the term due to friction
    if (_force_compute.stress_compute())
      for (FaceHandle<Property> fh = _sys.mesh().faces().begin(); fh != _sys.mesh().faces().end(); fh++)
      {
        if (!(fh->erased || fh->outer)) 
        {
          _force_compute.compute_stress(fh);
          HEHandle<Property> he = fh->he();
          HEHandle<Property> first = fh->he();
          Vec rc = _sys.mesh().get_face_centre(fh);
          vector<double> stress(4,0.0);
          do
          {
            double z = he->from()->coordination;
            Vec ri = he->from()->r - rc;
            Vec vi = mu*he->from()->data().force;
            stress[0] -= ri.x*vi.x/z;
            stress[1] -= 0.5*(ri.x*vi.y + ri.y*vi.x)/z;
            stress[2] -= 0.5*(ri.x*vi.y + ri.y*vi.x)/z;
            stress[3] -= ri.y*vi.y/z;
            he = he->next();
          } while (he != first);
          double fact = _gamma/_sys.mesh().area(fh);
          transform(stress.begin(),stress.end(),stress.begin(),[fact](double v) -> double {return fact*v; });
          transform(fh->data().stress.begin(),fh->data().stress.end(), stress.begin(), fh->data().stress.begin(), std::plus<double>());
          transform(fh->data().stress_v.begin(),fh->data().stress_v.end(), stress.begin(), fh->data().stress_v.begin(), std::plus<double>());
        }
      }
    
    // This is actual integrator 
    vector<Vec> k1(_sys.mesh().vertices().size(), Vec(0.0, 0.0)); 
    vector<Vec> k2(_sys.mesh().vertices().size(), Vec(0.0, 0.0)); 
    vector<Vec> k3(_sys.mesh().vertices().size(), Vec(0.0, 0.0)); 
    vector<Vec> k4(_sys.mesh().vertices().size(), Vec(0.0, 0.0)); 

    // Step 1:

    for (VertexHandle<Property> vh = _sys.mesh().vertices().begin(); vh != _sys.mesh().vertices().end(); vh++)
    {
      if (!vh->erased)
      {
        // add external force 
        Vec f = vh->data().force + _constant_force[vh->data().vert_type];

        if (_has_radial_force)
          f = f + _radial_force_magnitude[vh->data().vert_type] * vh->r.unit();

        // apply constraint
        if (_constraint_enabled) 
          f = _constrainer->apply_vertex(vh, f);

        k1[vh->id] = _dt*mu*f;
      } 
    }

    // Virtual move of vertices 
    for (VertexHandle<Property> vh = _sys.mesh().vertices().begin(); vh != _sys.mesh().vertices().end(); vh++)
      if (!vh->erased)
        vh->r += 0.5*k1[vh->id];

    // Compute force
    for (VertexHandle<Property> vh = _sys.mesh().vertices().begin(); vh != _sys.mesh().vertices().end(); vh++)
      if (!vh->erased)
        _force_compute.compute(vh);

    // Move back
    // Virtual move of vertices 
    for (VertexHandle<Property> vh = _sys.mesh().vertices().begin(); vh != _sys.mesh().vertices().end(); vh++)
      if (!vh->erased)
        vh->r -= 0.5*k1[vh->id];

    // Step 2:
    for (VertexHandle<Property> vh = _sys.mesh().vertices().begin(); vh != _sys.mesh().vertices().end(); vh++)
    {
      if (!vh->erased)
      {
        // add external force 
        Vec f = vh->data().force + _constant_force[vh->data().vert_type];

        if (_has_radial_force)
          f = f + _radial_force_magnitude[vh->data().vert_type] * vh->r.unit();

        // apply constraint
        if (_constraint_enabled) 
          f = _constrainer->apply_vertex(vh, f);

        k2[vh->id] = _dt*mu*f;
      } 
    }

    // Virtual move of vertices 
    for (VertexHandle<Property> vh = _sys.mesh().vertices().begin(); vh != _sys.mesh().vertices().end(); vh++)
      if (!vh->erased)
        vh->r += 0.5*k2[vh->id];

    // Compute force
    for (VertexHandle<Property> vh = _sys.mesh().vertices().begin(); vh != _sys.mesh().vertices().end(); vh++)
      if (!vh->erased)
        _force_compute.compute(vh);

    // Move back
    // Virtual move of vertices 
    for (VertexHandle<Property> vh = _sys.mesh().vertices().begin(); vh != _sys.mesh().vertices().end(); vh++)
      if (!vh->erased)
        vh->r -= 0.5*k2[vh->id];

    // Step 3:
    for (VertexHandle<Property> vh = _sys.mesh().vertices().begin(); vh != _sys.mesh().vertices().end(); vh++)
    {
      if (!vh->erased)
      {
        // add external force 
        Vec f = vh->data().force + _constant_force[vh->data().vert_type];

        if (_has_radial_force)
          f = f + _radial_force_magnitude[vh->data().vert_type] * vh->r.unit();

        // apply constraint
        if (_constraint_enabled) 
          f = _constrainer->apply_vertex(vh, f);

        k3[vh->id] = _dt*mu*f;
      } 
    }

    // Virtual move of vertices 
    for (VertexHandle<Property> vh = _sys.mesh().vertices().begin(); vh != _sys.mesh().vertices().end(); vh++)
      if (!vh->erased)
        vh->r += k3[vh->id];

    // Compute force
    for (VertexHandle<Property> vh = _sys.mesh().vertices().begin(); vh != _sys.mesh().vertices().end(); vh++)
      if (!vh->erased)
        _force_compute.compute(vh);

    // Move back
    // Virtual move of vertices 
    for (VertexHandle<Property> vh = _sys.mesh().vertices().begin(); vh != _sys.mesh().vertices().end(); vh++)
      if (!vh->erased)
        vh->r -= k3[vh->id];

    // Step 4:
    for (VertexHandle<Property> vh = _sys.mesh().vertices().begin(); vh != _sys.mesh().vertices().end(); vh++)
    {
      if (!vh->erased)
      {
        // add external force 
        Vec f = vh->data().force + _constant_force[vh->data().vert_type];

        if (_has_radial_force)
          f = f + _radial_force_magnitude[vh->data().vert_type] * vh->r.unit();

        // apply constraint
        if (_constraint_enabled) 
          f = _constrainer->apply_vertex(vh, f);

        k4[vh->id] = _dt*mu*f;
      } 
    }

    // Final step
    for (VertexHandle<Property> vh = _sys.mesh().vertices().begin(); vh != _sys.mesh().vertices().end(); vh++)
      if (!vh->erased)
        vh->r += (1.0/6.0)*(k1[vh->id] + 2*k2[vh->id] + 2*k3[vh->id] + k4[vh->id]);


  }

}