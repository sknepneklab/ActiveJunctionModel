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
 * \file force_active_area.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 20-May-2019
 * \brief ForceActiveArea class 
*/ 

#include "force_active_area.hpp"

namespace AJM 
{
  Vec ForceActiveArea::compute(const VertexHandle<Property>& vh, const HEHandle<Property>& he)
  {
    Vec l = he->to()->r - vh->r;                    // vector along the junction pointing away from the vertex
    FaceHandle<Property> fh   = he->face();         // cell to the right of the half edge
    FaceHandle<Property> fh_p = he->pair()->face(); // pair cell (opposite side of the same junction)
  
    double beta_1, beta_2;
    
    if (_use_cell_type)
    {
      beta_1 = (fh->outer)    ? 0.0 : _beta[fh->data().face_type];
      beta_2 = (fh_p->outer)  ? 0.0 : _beta[fh_p->data().face_type];
    }
    else
    {
      beta_1  = (fh->outer)   ? 0.0 : fh->data().beta_a;
      beta_2  = (fh_p->outer) ? 0.0 : fh_p->data().beta_a;
    }

    beta_2 = fabs(beta_1) < SMALL_NUMBER ? 0.0 : beta_2;
    beta_1 = fabs(beta_2) < SMALL_NUMBER ? 0.0 : beta_1;
    
    double factive;
    
    if (_sys.symmetric_myosin())
      factive = beta_1*(he->data().myo - 0.5) - beta_2*(he->pair()->data().myo - 0.5);
    else
      factive = beta_1*(he->data().myo) - beta_2*(he->pair()->data().myo);

    return factive*l.ez_cross_v().unit();
  }

  double ForceActiveArea::tension(const HEHandle<Property>& he, double myo, double l0)
  {
    return 0.0;
  }
}