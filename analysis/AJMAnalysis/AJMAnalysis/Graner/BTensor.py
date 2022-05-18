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
# \file BTensor.py
# \author Rastko Sknepnek, sknepnek@gmail.com
# \date 11-May-2020
# Computes the B tensor (Graner statistics) and produces a VTP file for it
#

import json
import os

import numpy as np

from ..utils.HalfEdge import *
from .Tensor import Tensor


class BTensor:

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
        self.B = Tensor(self.mesh_2.num_inner_faces)
        self.C = Tensor(self.mesh_2.num_inner_faces)
        if self.mesh_1.system['time'] and self.mesh_2.system['time']:
            dt = self.mesh_2.system['time'] - self.mesh_1.system['time']

        idx = 0
        for f in range(self.mesh_2.num_inner_faces):
            f_2 = self.mesh_2.faces[f]
            if (self.mesh_1.num_inner_faces == self.mesh_2.num_inner_faces) or ("original_id" not in self.mesh_2.faces[f].params):
                f_1 = self.mesh_1.faces[f]
            else:
                f_1 = self.mesh_1.faces[self.mesh_2.faces[f].params["original_id"]]
            if not (f_1.outer or f_2.outer):
                # Number of links of either cells
                neigh_1, neigh_2 = f_1.neighbours(), f_2.neighbours()
                # Conserved links
                neigh = list(set(neigh_1) & set(neigh_2))
                nv_1 = f_1.neigh_vectors(self.mesh_1.box)
                nv_2 = f_2.neigh_vectors(self.mesh_2.box)
                for n in neigh:
                    idx_1, idx_2 = neigh_1.index(n), neigh_2.index(n)
                    # l vectors at two consecutive time steps
                    l_1, l_2 = nv_1[idx_1], nv_2[idx_2]
                    l = 0.5*(l_1 + l_2)
                    dl = l_2 - l_1
                    # Adding up the conserved links
                    self.C.T[idx, :, :] += np.array(l.outer(dl))
                # Need to normalise by the total number of links. Per cell, that comes out to (symmetrised) 0.5(len(neigh_1) + len(neigh_2))
                zav = 0.5*(len(neigh_1) + len(neigh_2))
                self.C.T[idx, :, :] = (1/dt)*self.C.T[idx, :, :]/zav
                self.B.T[idx, :, :] = self.C.T[idx, :, :] + \
                    self.C.T[idx, :, :].T
                idx += 1

    def plot_vtk_tensor(self, filename):
        self.B.plot_vtk_tensor(filename, self.mesh_2, 'B_tensor')
