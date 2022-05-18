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
 * \file force_harmonic_spokes.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 20-May-2019
 * \brief ForceHarmonicSpokes class 
*/ 

#include "force_harmonic_spokes.hpp"

namespace AJM 
{
  Vec ForceHarmonicSpokes::compute(const VertexHandle<Property>& vh, const HEHandle<Property>& he)
  {
    FaceHandle<Property> fh   = he->face();         // cell to the right of the half edge
    FaceHandle<Property> fh_p = he->pair()->face(); // pair cell (opposite side of the same junction)
  
    double k_1, k_2;
    
    if (_use_cell_type)
    {
      k_1 = (fh->outer)    ? 0.0 : _k[fh->data().face_type];
      k_2 = (fh_p->outer)  ? 0.0 : _k[fh_p->data().face_type];
    }
    else
    {
      k_1  = (fh->outer)   ? 0.0 : fh->data().k_spoke;
      k_2  = (fh_p->outer) ? 0.0 : fh_p->data().k_spoke;
    }
    
    Vec fs(0,0);
    if (!fh->outer)
    {
      Vec l = _sys.mesh().get_face_centre(fh) - vh->r;
      fs = k_1 * (l.len() - fh->data().spoke[vh->id].l0) * l.unit();
    }
    if (!fh_p->outer)
    {
      Vec l = _sys.mesh().get_face_centre(fh_p) - vh->r;
      fs = k_2 * (l.len() - fh_p->data().spoke[vh->id].l0) * l.unit();
    }  
    return fs; 
  }

  double ForceHarmonicSpokes::tension(const HEHandle<Property>& he, double myo, double l0)
  {
    return 0.0;
  }

}