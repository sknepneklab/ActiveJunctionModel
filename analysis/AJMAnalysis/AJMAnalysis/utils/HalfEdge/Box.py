###########################################################################
#
#  Copyright (C) 2017, 2018, 2019, 2020 University of Dundee
#  All rights reserved.
#
#  This file is part of AJM (Active Junction Model) program.
#
#  AJM is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  AJM is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ############################################################################

#
# \file Box.py
# \author Rastko Sknepnek, sknepnek@gmail.com
# \date 11-May-2020
# Handles simulation box data
#

import numpy as np


class Box:

    def __init__(self, a, b=None):
        """
        Construct simulation box.
          Parameters
          ----------
          a : list 
             Components of the a vector of the simulation box a = [ax, ay]
          b : list
             Components of the a vector of the simulation box a = [ax, ay]
          Note
          ----
            Simulation box is centred at (0,0).
        """
        has_b = None
        if type(a) is not list:
            raise ValueError('a vector has to be a list.')
        if len(a) != 2:
            raise ValueError('a vector has to be give with two components.')
        if b is not None:
            if (type(b) is list) or (len(b) == 2):
                has_b = True

        self.a = np.array(a)
        if has_b:
            self.b = np.array(b)
        else:
            self.b = np.array([a[1], a[0]])

        self.h = np.vstack((self.a, self.b)).T
        self.invh = np.linalg.inv(self.h)
