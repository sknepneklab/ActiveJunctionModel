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
# \file Face.py
# \author Rastko Sknepnek, sknepnek@gmail.com
# \date 11-May-2020
# Handles faces in the half-edge implementation of the mesh

from .Vec import Vec


class Face:

    def __init__(self, idx):
        self.idx = idx
        self.he = None
        self.outer = False
        self.type = None
        self.params = {}

    def area(self):
        # Use shoe-lace formula to compute area
        # By defintion, outer face will have negative area
        first = self.he
        he = self.he
        area = 0.0
        r0 = he.vfrom.r
        while True:
            vi = he.vfrom
            vj = he.vto
            ri = vi.r - r0
            rj = vj.r - r0
            area += ri.r[0]*rj.r[1] - ri.r[1]*rj.r[0]
            he = he.next
            if he.idx == first.idx:
                break
        return 0.5*abs(area)

    def perimeter(self):
        first = self.he
        he = self.he
        perim = 0.0
        while True:
            ri = he.vfrom.r
            rj = he.vto.r
            dr = rj - ri
            perim += dr.length()
            he = he.next
            if he.idx == first.idx:
                break
        return perim

    def rc(self, box=None):
        first = self.he
        he = self.he
        Rc = Vec([0.0, 0.0])
        N = 0
        while True:
            dr = he.vfrom.r - first.vfrom.r
            Rc += Vec(dr.r)
            N += 1
            he = he.next
            if he.idx == first.idx:
                break
        Rc = (1.0/N)*Rc + first.vfrom.r
        return Vec(Rc.r, box)

    def neighbours(self):
        neigh = []
        first = self.he
        he = self.he
        while True:
            if not he.pair.face.outer:
                neigh.append(he.pair.face.idx)
            he = he.next
            if he.idx == first.idx:
                break
        return neigh

    def neigh_vectors(self, box=None):
        nv = []
        first = self.he
        he = self.he
        Rc = self.rc(box)
        while True:
            if not he.pair.face.outer:
                Rci = he.pair.face.rc(box)
                dR = Rci - Rc
                nv.append(dR)
            he = he.next
            if he.idx == first.idx:
                break
        return nv

    def edge_vectors(self):
        ev = []
        first = self.he
        he = self.he
        while True:
            dR = (he.vto.r - he.vfrom.r).unit()
            ev.append(dR)
            he = he.next
            if he.idx == first.idx:
                break
        return ev

    def num_sides(self):
        first = self.he
        he = self.he
        N = 0
        while True:
            N += 1
            he = he.next
            if he.idx == first.idx:
                break
        return N

    def centroid(self):
        first = self.he
        he = self.he
        r0 = he.vfrom.r
        Rc = Vec([0.0, 0.0])
        while True:
            ri = he.vfrom.r - r0
            rj = he.vto.r - r0
            fact = ri.r[0]*rj.r[1] - ri.r[1]*rj.r[0]
            Rc.r[0] += (ri.r[0] + rj.r[0])*fact
            Rc.r[1] += (ri.r[1] + rj.r[1])*fact
            he = he.next
            if he.idx == first.idx:
                break
        return (1.0/(6*self.area()))*Rc + r0
