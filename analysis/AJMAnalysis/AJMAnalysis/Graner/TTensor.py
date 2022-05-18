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
# \file TTensor.py
# \author Rastko Sknepnek, sknepnek@gmail.com
# \date 12-May-2020
# Computes the T tensor (Graner statistics) and produces a VTP file for it
#

import json
import os

import numpy as np

from ..utils.HalfEdge import *
from .Tensor import Tensor


class TTensor:

    def __init__(self, frame_1, frame_2, dt=1.0):
        ext_1 = os.path.splitext(frame_1)[1]
        ext_2 = os.path.splitext(frame_2)[1]
        if (ext_1 != '.json' and ext_1 != '.bz2') or (ext_2 != '.json' and ext_2 != '.bz2'):
            raise Exception(
                'JSON files with two consecutive frames of the mesh state have to be provided.')
        self.mesh_1 = Mesh()
        self.mesh_2 = Mesh()
        self.mesh_1.read(frame_1)
        self.mesh_2.read(frame_2)
        self.T = Tensor(self.mesh_2.num_inner_faces)
        if self.mesh_1.system['time'] and self.mesh_2.system['time']:
            dt = self.mesh_2.system['time'] - self.mesh_1.system['time']
        fidx = 0
        # Implicit assumption that faces are not relabeled here
        for (f_1, f_2) in zip(self.mesh_1.faces, self.mesh_2.faces):
            if not (f_1.outer or f_2.outer):
                neigh_1, neigh_2 = f_1.neighbours(), f_2.neighbours()
                # compute l(t) and l(t+dt)
                nv_1 = f_1.neigh_vectors(self.mesh_1.box)
                nv_2 = f_2.neigh_vectors(self.mesh_2.box)
                # links that appear after a T1 transition
                appear = list(set(neigh_2) - set(neigh_1))
                # links that dissapear in T1 transition
                disappear = list(set(neigh_1) - set(neigh_2))
                ma = np.zeros((2, 2))
                md = np.zeros((2, 2))
                if len(appear) > 0:
                    for a in appear:
                        idx = neigh_2.index(a)
                        l = nv_2[idx]
                        ma += np.array(l.outer(l))
                if len(disappear) > 0:
                    for d in disappear:
                        idx = neigh_1.index(d)
                        l = nv_1[idx]
                        md += np.array(l.outer(l))
                # Do not keep the time here for now
                # Use central average of contact numbers again
                zav = 0.5*(len(neigh_1) + len(neigh_2))
                self.T.T[fidx, :, :] = (1/dt)*(ma - md)/zav
                fidx += 1

    def plot_vtk_tensor(self, filename):
        self.T.plot_vtk_tensor(filename, self.mesh_2, 'T_tensor')
