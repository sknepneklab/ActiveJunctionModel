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
 * \file myostore.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 13-Mar-2019
 * \brief MyoStore class 
*/

#ifndef __MYOSTORE_HPP__
#define __MYOSTORE_HPP__

#include "system.hpp"

#include <map>
#include <vector>
#include <algorithm>
#include <iostream>

using std::map;
using std::vector;
using std::min;
using std::cout;
using std::endl;

namespace AJM
{
  template<typename Property>
  class MyoStore 
  {
    public:

      MyoStore() { }
      MyoStore(HEHandle<Property>& he, double myo) : _myo(myo), _myo_dep(4,0.0)
      {
        HEHandle<Property> hep = he->pair();
        _myosin[he->face()->id] = he->data().myo;
        _myosin[hep->face()->id] = hep->data().myo;
        _myosin[he->prev()->pair()->face()->id] = 0.0;
        _myosin[hep->prev()->pair()->face()->id] = 0.0;
        _faces.push_back(he->face()->id);
        _faces.push_back(hep->face()->id);
        _faces.push_back(he->prev()->pair()->face()->id);
        _faces.push_back(hep->prev()->pair()->face()->id);
        _l0 = he->data().l0;
      }

      void spread_backward(HEHandle<Property>& he)
      {
        he->data().myo = this->_myosin[he->face()->id];
        he->next()->pair()->data().myo += 0.5*this->_myosin[he->next()->pair()->face()->id];
        he->prev()->pair()->data().myo += 0.5*this->_myosin[he->prev()->pair()->face()->id];
        HEHandle<Property> hep = he->pair();
        hep->data().myo = this->_myosin[hep->face()->id];
        hep->next()->pair()->data().myo += 0.5*this->_myosin[hep->next()->pair()->face()->id];
        hep->prev()->pair()->data().myo += 0.5*this->_myosin[hep->prev()->pair()->face()->id];
      }

      void spread_forward(HEHandle<Property>& he)
      {
        he->data().myo = _myo;
        he->next()->data().myo -= min(he->next()->data().myo, 0.5*_myo);
        he->prev()->data().myo -= min(he->prev()->data().myo, 0.5*_myo);
        he->prev()->pair()->data().myo += 0.5*_myosin[he->prev()->pair()->face()->id];
        he->next()->pair()->data().myo += 0.5*_myosin[he->next()->pair()->face()->id];
        HEHandle<Property> hep = he->pair();
        hep->data().myo = _myo;
        hep->next()->data().myo -= min(hep->next()->data().myo, 0.5*_myo);
        hep->prev()->data().myo -= min(hep->prev()->data().myo, 0.5*_myo);
        hep->prev()->pair()->data().myo += 0.5*_myosin[hep->prev()->pair()->face()->id];
        hep->next()->pair()->data().myo += 0.5*_myosin[hep->next()->pair()->face()->id];
      }
  
      void virtual_spread(HEHandle<Property>& he, bool sign) 
      {
        HEHandle<Property>& hep = he->pair();
        double sgn = sign ? 1.0 : -1.0;
        if (sign)
        {
          _myo_dep[0] = min(he->prev()->data().myo, 0.5*_myo);
          _myo_dep[1] = min(he->next()->data().myo, 0.5*_myo);
          _myo_dep[2] = min(hep->prev()->data().myo, 0.5*_myo);
          _myo_dep[3] = min(hep->next()->data().myo, 0.5*_myo);
        }
        
        he->next()->pair()->data().myo += 0.5*sgn*_myosin[he->next()->pair()->face()->id];
        he->prev()->pair()->data().myo += 0.5*sgn*_myosin[he->prev()->pair()->face()->id];

        // We also take part of the myosin from previous edges and assign them to the new edge
        
        he->prev()->data().myo -=  sgn*_myo_dep[0];
        he->next()->data().myo -=  sgn*_myo_dep[1];

        // Same thing for the pair
        hep->next()->pair()->data().myo += 0.5*sgn*_myosin[hep->next()->pair()->face()->id];
        hep->prev()->pair()->data().myo += 0.5*sgn*_myosin[hep->prev()->pair()->face()->id];
        
        hep->prev()->data().myo -=  sgn*_myo_dep[2];
        hep->next()->data().myo -=  sgn*_myo_dep[3];
      }
   
      double& get_myo(int f) 
      { 
        if (_myosin.find(f) == _myosin.end())
          throw runtime_error("Face does not exist in myo.");
        return _myosin[f]; 
      }
      double& get_l0() { return _l0; }
      int get_face(int f) { return _faces[f]; }

    private:
      
      vector<int> _faces;       // constans all faces around the 
      map<int,double> _myosin;  // myosin for all half edges
      vector<double> _myo_dep;  
      double _l0;                // l0 before collapse
      double _myo;               // myosin base level

  };
}


#endif
