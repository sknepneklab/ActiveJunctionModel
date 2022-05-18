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
 * \file force_line_tension.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 12-Jun-2019
 * \brief ForceLineTension class 
*/ 

#include "force_line_tension.hpp"

namespace AJM 
{
  Vec ForceLineTension::compute(const VertexHandle<Property>& vh, const HEHandle<Property>& he)
  {
    Vec l = he->to()->r - vh->r;                    // vector along the junction pointing away from the vertex
    FaceHandle<Property> fh   = he->face();         // cell to the right of the half edge
    FaceHandle<Property> fh_p = he->pair()->face(); // pair cell (opposite side of the same junction)
    double lambda_1, lambda_2;
    
    
    if (_use_cell_type)
    {
      lambda_1 = (fh->outer)    ? 0.0 : _lambda[fh->data().face_type];
      lambda_2 = (fh_p->outer)  ? 0.0 : _lambda[fh_p->data().face_type];
    }
    else
    {
      lambda_1  = (fh->outer)   ? 0.0 : fh->data().lambda;
      lambda_2  = (fh_p->outer) ? 0.0 : fh_p->data().lambda;
    }
    if (_double_boundary_constants)
    {
      lambda_1  = (fh_p->outer) ? 2*lambda_1 : lambda_1;
      lambda_2  = (fh->outer)   ? 2*lambda_2 : lambda_2;
    }

    double lambda = lambda_1 + lambda_2;

    return -lambda*l.unit();
  }

  double ForceLineTension::tension(const HEHandle<Property>& he, double myo, double l0)
  {
    FaceHandle<Property> fh   = he->face();         // cell to the right of the half edge
    FaceHandle<Property> fh_p = he->pair()->face(); // pair cell (opposite side of the same junction)

    double lambda_1, lambda_2;
    
    if (_use_cell_type)
    {
      lambda_1 = (fh->outer)    ? 0.0 : _lambda[fh->data().face_type];
      lambda_2 = (fh_p->outer)  ? 0.0 : _lambda[fh_p->data().face_type];
    }
    else
    {
      lambda_1  = (fh->outer)   ? 0.0 : fh->data().lambda;
      lambda_2  = (fh_p->outer) ? 0.0 : fh_p->data().lambda;
    }
    if (_double_boundary_constants)
    {
      lambda_1  = (fh_p->outer) ? 2*lambda_1 : lambda_1;
      lambda_2  = (fh->outer)   ? 2*lambda_2 : lambda_2;
    }
    double lambda = lambda_1 + lambda_2;

    return -lambda;  // Note: Assumption here is that edge has zero length.
  }

  double ForceLineTension::energy(const FaceHandle<Property>& fh)
  {
    HEHandle<Property>& he = fh->he();
    HEHandle<Property>& first = fh->he();
    if (fh->outer || fh->erased)
      return 0.0;
    
    double lambda;
    if (_use_cell_type)
      lambda = _lambda[fh->data().face_type];
    else
      lambda = fh->data().lambda;

    double eng = 0.0;
    do 
    {
      Vec r_from = he->from()->r;
      Vec r_to = he->to()->r;
      double l = (r_to - r_from).len();
      eng += -lambda*l;
      he = he->next();
    } while (he != first);
    return eng;
  }

}