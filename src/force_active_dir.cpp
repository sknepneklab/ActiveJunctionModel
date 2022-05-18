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
 * \file force_active_dir.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 14-Jul-2019
 * \brief ForceActiveDir class 
*/ 

#include "force_active_dir.hpp"

namespace AJM 
{
  Vec ForceActiveDir::compute(const VertexHandle<Property>& vh, const HEHandle<Property>& he)
  {
    Vec l = he->to()->r - vh->r;                    // vector along the junction pointing away from the vertex
    FaceHandle<Property> fh   = he->face();         // cell to the right of the half edge
    FaceHandle<Property> fh_p = he->pair()->face(); // pair cell (opposite side of the same junction)
  
    double beta_1, beta_2;
    double alpha_1, alpha_2;
    int psi_1, psi_2;
    
    if (_use_cell_type)
    {
      beta_1 = (fh->outer)    ? 0.0 : _beta[fh->data().face_type];
      beta_2 = (fh_p->outer)  ? 0.0 : _beta[fh_p->data().face_type];
      alpha_1 = (fh->outer)   ? 0.0 : _alpha[fh->data().face_type];
      alpha_2 = (fh_p->outer) ? 0.0 : _alpha[fh_p->data().face_type];
      psi_1 = (fh->outer)     ? 0   : _psi[fh->data().face_type];
      psi_2 = (fh_p->outer)   ? 0   : _psi[fh_p->data().face_type];
    }
    else
    {
      beta_1  = (fh->outer)   ? 0.0 : fh->data().beta;
      beta_2  = (fh_p->outer) ? 0.0 : fh_p->data().beta;
      alpha_1 = (fh->outer)   ? 0.0 : fh->data().alpha;
      alpha_2 = (fh_p->outer) ? 0.0 : fh_p->data().alpha;
      psi_1   = (fh->outer)   ? 0   : fh->data().psi;
      psi_2   = (fh_p->outer) ? 0   : fh_p->data().psi;
    }

    if ((psi_1 % 2 != 0) || (psi_2 % 2 != 0))
      throw runtime_error("Power in the active direction coupling has to be even.");
    
    double alpha = alpha_1 + alpha_2;
    double factive;

    Vec n_1  = _sys.mesh().get_face_direction(fh);
    Vec n_2  = _sys.mesh().get_face_direction(fh_p);

    Vec ul = l.unit();
    double n_dot_l_1 = n_1.dot(ul);
    double n_dot_l_2 = n_2.dot(ul);

    if (psi_1 == 2)
      beta_1 = (1.0 + alpha_1*(n_dot_l_1*n_dot_l_1 - 1.0));
    else
      beta_1 = (1.0 + alpha_1*(pow(n_dot_l_1,psi_1) - 1.0));

    if (psi_2 == 2)
      beta_2 = (1.0 + alpha_2*(n_dot_l_2*n_dot_l_2 - 1.0));
    else
      beta_2 = (1.0 + alpha_2*(pow(n_dot_l_2,psi_2) - 1.0));
    
    if (_sys.symmetric_myosin())
      factive = 0.5*(beta_1*(he->data().myo - 0.5) + beta_2*(he->pair()->data().myo - 0.5));
    else
      factive = 0.5*(beta_1*(he->data().myo) + beta_2*(he->pair()->data().myo));

    return factive*l.unit();
  }

  double ForceActiveDir::tension(const HEHandle<Property>& he, double myo, double l0)
  {
    Vec l = he->to()->r - he->from()->r;
    FaceHandle<Property> fh   = he->face();         // cell to the right of the half edge
    FaceHandle<Property> fh_p = he->pair()->face(); // pair cell (opposite side of the same junction)
  
    double beta_1, beta_2;
    double alpha_1, alpha_2;
    int psi_1, psi_2;
    
    if (_use_cell_type)
    {
      beta_1 = (fh->outer)    ? 0.0 : _beta[fh->data().face_type];
      beta_2 = (fh_p->outer)  ? 0.0 : _beta[fh_p->data().face_type];
      alpha_1 = (fh->outer)   ? 0.0 : _alpha[fh->data().face_type];
      alpha_2 = (fh_p->outer) ? 0.0 : _alpha[fh_p->data().face_type];
      psi_1 = (fh->outer)     ? 0   : _psi[fh->data().face_type];
      psi_2 = (fh_p->outer)   ? 0   : _psi[fh_p->data().face_type];
    }
    else
    {
      beta_1  = (fh->outer)   ? 0.0 : fh->data().beta;
      beta_2  = (fh_p->outer) ? 0.0 : fh_p->data().beta;
      alpha_1 = (fh->outer)   ? 0.0 : fh->data().alpha;
      alpha_2 = (fh_p->outer) ? 0.0 : fh_p->data().alpha;
      psi_1   = (fh->outer)   ? 0   : fh->data().psi;
      psi_2   = (fh_p->outer) ? 0   : fh_p->data().psi;
    }

    Vec n_1  = _sys.mesh().get_face_direction(fh);
    Vec n_2  = _sys.mesh().get_face_direction(fh_p);

    Vec ul = l.unit();
    double n_dot_l_1 = n_1.dot(ul);
    double n_dot_l_2 = n_2.dot(ul);

    if (psi_1 == 2)
      beta_1 = (1.0 + alpha_1*(n_dot_l_1*n_dot_l_1 - 1.0));
    else
      beta_1 = (1.0 + alpha_1*(pow(n_dot_l_1,psi_1) - 1.0));

    if (psi_2 == 2)
      beta_2 = (1.0 + alpha_2*(n_dot_l_2*n_dot_l_2 - 1.0));
    else
      beta_2 = (1.0 + alpha_2*(pow(n_dot_l_2,psi_2) - 1.0));

    double beta = 0.5*(beta_1 + beta_2);
    if (_sys.symmetric_myosin())
      return beta*(myo - 0.5);
    else
      return beta*myo;
  }
}