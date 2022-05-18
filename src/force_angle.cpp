/* ***************************************************************************
 *
 *  Copyright (C) 2017 University of Dundee
 *  All rights reserved. 
 *
 *  This file is part of AJM (Angle Junction Model) program.
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
 * \file force_angle.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 14-Jun-2019
 * \brief ForceAngle class 
*/ 

#include "force_angle.hpp"

namespace AJM 
{
  Vec ForceAngle::compute(const VertexHandle<Property>& vh, const HEHandle<Property>& he)
  {
    double k;
    if (_use_cell_type)
      k = _k[he->face()->data().face_type];
    else
      k  = he->face()->data().k_angle;
    
    VertexHandle<Property> vh_i  = he->from();
    VertexHandle<Property> vh_ip = he->to();
    VertexHandle<Property> vh_im = he->prev()->from();
    Vec r1 = vh_ip->r - vh_i->r;
    Vec r2 = vh_im->r - vh_i->r;
    double cross = r1.x*r2.y - r2.x*r1.y;
    
    if (cross < 0.0)
    {
      double c = dot(r1,r2)/(r1.len()*r2.len());

      if (c > 1.0) c = 1.0;
      else if (c < -1.0) c = -1.0;
      
      double s = sqrt(1.0 - c*c);
      if (s < 1e-7) s = 1e7;
      else s = 1.0/s;
      
      double dtheta = acos(c) - M_PI;
      double tk = k*dtheta;

      double aa = -2.0 * tk * s;
      double a11 = aa*c / r1.len2();
      double a12 = -aa / (r1.len()*r2.len());
      double a22 = aa*c / r2.len2();

      Vec fi_p = a11*r1 + a12*r2;
      Vec fi_m = a12*r1 + a22*r2;

      return -(fi_p+fi_m);
    }
    else
      return Vec(0,0);
  }

  double ForceAngle::tension(const HEHandle<Property>& he, double myo, double l0)
  {
    return 0.0;
  }

  // NOTE: This is clearly well defined and we need to implement it
  double ForceAngle::energy(const FaceHandle<Property>& fh)
  {
    return 0.0;
  }
}