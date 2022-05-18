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
 * \file force_active_spokes.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 20-May-2019
 * \brief ForceActiveSpokes class 
*/ 

#include "force_active_spokes.hpp"

namespace AJM 
{
  Vec ForceActiveSpokes::compute(const VertexHandle<Property>& vh, const HEHandle<Property>& he)
  {
    FaceHandle<Property> fh   = he->face();         // face to the right of the half edge; we only need one face to avoid double counting
                                                    // when looping around the vertex
  
    if (fh->outer)
      return Vec(0,0);

    double beta;
    
    if (_use_cell_type)
      beta = (fh->outer)  ? 0.0 : _beta[fh->data().face_type];
    else
      beta = (fh->outer)  ? 0.0 : fh->data().beta_a;

    double myo_1 = (_sys.symmetric_myosin()) ? (he->data().myo - 0.5) : he->data().myo;
    double myo_2 = (_sys.symmetric_myosin()) ? (he->prev()->data().myo - 0.5) : he->prev()->data().myo;
    
    Vec dr = _sys.mesh().get_face_centre(fh) - vh->r;
    Vec fs = beta*(myo_1 + myo_2)*dr.unit();
    return fs; 
  }

  double ForceActiveSpokes::tension(const HEHandle<Property>& he, double myo, double l0)
  {
    return 0.0;
  }
}