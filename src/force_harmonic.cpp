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
 * \file force_harmonic.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 20-May-2019
 * \brief ForceHarmonic class 
*/ 

#include "force_harmonic.hpp"

namespace AJM 
{
  Vec ForceHarmonic::compute(const VertexHandle<Property>& vh, const HEHandle<Property>& he)
  {
    Vec l = he->to()->r - vh->r;                    // vector along the junction pointing away from the vertex
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
      k_1  = (fh->outer)   ? 0.0 : fh->data().k;
      k_2  = (fh_p->outer) ? 0.0 : fh_p->data().k;
    }
    k_1  = (fh_p->outer) ? 2*k_1 : k_1;
    k_2  = (fh->outer)   ? 2*k_2 : k_2;
    
    double fharmonic = (k_1 + k_2)*(l.len() - he->edge()->data().l0);

    return fharmonic*l.unit();
  }

  double ForceHarmonic::tension(const HEHandle<Property>& he, double myo, double l0)
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
      k_1  = (fh->outer)   ? 0.0 : fh->data().k;
      k_2  = (fh_p->outer) ? 0.0 : fh_p->data().k;
    }
    k_1  = (fh_p->outer) ? 2*k_1 : k_1;
    k_2  = (fh->outer)   ? 2*k_2 : k_2;
  
    return -(k_1 + k_2)*l0;  // Note: Assumption here is that edge has zero length.
  }

  double ForceHarmonic::energy(const FaceHandle<Property>& fh)
  {
    HEHandle<Property>& he = fh->he();
    HEHandle<Property>& first = fh->he();
    if (fh->outer || fh->erased)
      return 0.0;
    
    double k;
    if (_use_cell_type)
      k = _k[fh->data().face_type];
    else
      k = fh->data().k;

    double eng = 0.0;
    do 
    {
      Vec r_from = he->from()->r;
      Vec r_to = he->to()->r;
      double dl = (r_to - r_from).len() - he->edge()->data().l0;
      eng += 0.5*k*dl*dl;
      he = he->next();
    } while (he != first);
    return eng;
  }

}