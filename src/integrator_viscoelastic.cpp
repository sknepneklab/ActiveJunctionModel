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
 * \file integrator_viscoelastic.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 20-May-2019
 * \brief IntegratorViscoelastic class 
*/

#include "integrator_viscoelastic.hpp"

namespace AJM
{
  void IntegratorViscoelastic::step()
  {
    for (EdgeHandle<Property> eh = _sys.mesh().edges().begin(); eh != _sys.mesh().edges().end(); eh++)
    {
      double l = (eh->he()->to()->r - eh->he()->from()->r).len();
      double l0 = eh->data().l0;
      eh->data().l0 += _dt*(l - l0)/_tau_l;
    }
    if (_sys.has_spokes())
      for (FaceHandle<Property> fh = _sys.mesh().faces().begin(); fh != _sys.mesh().faces().end(); fh++)
      {
        Vec rc = _sys.mesh().get_face_centre(fh);
        HEHandle<Property> he = fh->he();
        HEHandle<Property> first = fh->he();
        do
        {
          VertexHandle<Property> vh = he->from();
          Vec l = rc - vh->r;
          Spoke& s = fh->data().spoke[vh->id];
          s.l0 += _dt * (l.len() - s.l0) / _tau_s;
          he = he->next();
        } while (he != first);
      }
  }
}
