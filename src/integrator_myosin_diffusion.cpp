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
 * \file integrator_myosin_diffusion.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Aug-2019
 * \brief IntegratorMyosinDiffusion class 
*/

#include "integrator_myosin_diffusion.hpp"

namespace AJM
{
  void IntegratorMyosinDiffusion::step()
  {
    for (FaceHandle<Property> fh = _sys.mesh().faces().begin(); fh != _sys.mesh().faces().end(); fh++)
    {
      HEHandle<Property> first = fh->he();
      HEHandle<Property> he = first;
      vector<double> he_myo;
      vector<double> spoke_myo;
      // Store myosin before one can update it in diffusion
      do
      {
        he_myo.push_back(he->data().myo);
        if (!fh->outer && _use_spokes)
          spoke_myo.push_back(fh->data().spoke[he->from()->id].myo);
      } while ((he = he->next()) != first);


      int i = 0;
      do
      {
        if (!fh->outer && _use_spokes)  // for inner edges we have diffusion between spokes and half-edge
        {
          int j = (i == (fh->nsides - 1)) ? 0 : i + 1;
          VertexHandle<Property> vh_from = he->from();
          VertexHandle<Property> vh_to   = he->to();
          he->data().myo += _dt*_D*(spoke_myo[i] - 2*he_myo[i] + spoke_myo[j]);
          fh->data().spoke[he->to()->id].myo += _dt * _D * (he_myo[i] - 2 * spoke_myo[j] + he_myo[j]);
          if (fh->data().spoke[he->to()->id].myo < 0.0) fh->data().spoke[he->to()->id].myo = 0.0;  
        }
        else  // for outer face or is spokes are note used, myosin diffuses to neighbours
        {
          int j = (i == (fh->nsides - 1)) ? 0 : i + 1;
          int k = (i == 0) ? fh->nsides - 1 : i - 1;
          he->data().myo += _dt * _D * (he_myo[k] - 2 * he_myo[i] + he_myo[j]);
        }
        if (he->data().myo < 0.0) he->data().myo = 0.0;  
        i++;
      } while ((he = he->next()) != first);
     }
  }

}
