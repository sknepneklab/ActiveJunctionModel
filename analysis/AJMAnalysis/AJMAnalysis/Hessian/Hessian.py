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
# \file Hessian.py
# \author Rastko Sknepnek, sknepnek@gmail.com
# \date 29-Jun-2020

from ..utils.HalfEdge import *
import numpy as np
import os.path
import vtk


class Hessian:

    def __init__(self, mesh, params=None):
        ext = os.path.splitext(mesh)[1][1:]
        if ext != 'json':
            raise Exception(
                'Hessian calculation needs a mesh file in JSON format.')
        self.mesh = Mesh()
        self.mesh.read(mesh)
        self.N = len(self.mesh.vertices)
        self.H = np.zeros((2*self.N, 2*self.N))
        self.params = params
        self.evaluated = False
        self.has_eigvals = False

    def __H13(self):
        for f in self.mesh.faces:
            if not f.outer:
                if self.params is None:
                    Kappa = f.params['kappa']
                    Gamma = f.params['gamma']
                else:
                    Kappa = self.params[f.type]['kappa']
                    Gamma = self.params[f.type]['gamma']
                he = f.he
                first = f.he
                while True:
                    k = he.vfrom.idx
                    rk = he.vfrom.r
                    rkp = he.vto.r
                    rkm = he.prev.vfrom.r
                    ak = (rkp - rkm).ez_cross()
                    pk = (rk - rkm).unit() - (rkp - rk).unit()
                    hem = he
                    while True:
                        m = hem.vfrom.idx
                        rm = hem.vfrom.r
                        rmp = hem.vto.r
                        rmm = hem.prev.vfrom.r
                        am = (rmp - rmm).ez_cross()
                        pm = (rm - rmm).unit() - (rmp - rm).unit()
                        self.H[2*k, 2*m] += 0.25*Kappa*ak.r[0] * \
                            am.r[0] + Gamma*pk.r[0]*pm.r[0]
                        self.H[2*k, 2*m + 1] += 0.25*Kappa * \
                            ak.r[0]*am.r[1] + Gamma*pk.r[0]*pm.r[1]
                        self.H[2*k + 1, 2*m] += 0.25*Kappa * \
                            ak.r[1]*am.r[0] + Gamma*pk.r[1]*pm.r[0]
                        self.H[2*k + 1, 2*m + 1] += 0.25*Kappa * \
                            ak.r[1]*am.r[1] + Gamma*pk.r[1]*pm.r[1]
                        hem = hem.next
                        if hem.idx == he.idx:
                            break
                    he = he.next
                    if he.idx == first.idx:
                        break

    def __H2(self):
        for f in self.mesh.faces:
            if not f.outer:
                if self.params is None:
                    Kappa = f.params['kappa']
                else:
                    Kappa = self.params[f.type]['kappa']
                A0 = f.params['A0']
                fact = 0.5*Kappa*(f.area() - A0)
                he = f.he
                first = f.he
                while True:
                    k = he.vfrom.idx
                    m1 = he.vto.idx
                    m2 = he.prev.vfrom.idx
                    self.H[2*k, 2*m1 + 1] += fact
                    self.H[2*k, 2*m2 + 1] -= fact
                    self.H[2*k + 1, 2*m1] -= fact
                    self.H[2*k + 1, 2*m2] += fact
                    he = he.next
                    if he.idx == first.idx:
                        break

    def __H4(self):
        for f in self.mesh.faces:
            if not f.outer:
                if self.params is None:
                    Gamma = f.params['gamma']
                else:
                    Gamma = self.params[f.type]['gamma']
                fact = Gamma*f.perimeter()
                he = f.he
                first = f.he
                while True:
                    k = he.vfrom.idx
                    kp = he.vto.idx
                    km = he.prev.vfrom.idx
                    rkkm = he.vfrom.r - he.prev.vfrom.r
                    rkpk = he.vto.r - he.vfrom.r
                    lrkkm = rkkm.length()
                    lrkpk = rkpk.length()
                    lrkkm3 = lrkkm*lrkkm*lrkkm
                    lrkpk3 = lrkpk*lrkpk*lrkpk
                    for alpha in range(2):
                        for beta in range(2):
                            self.H[2*k + alpha, 2*k + beta] -= fact*(
                                rkkm.r[alpha]*rkkm.r[beta]/lrkkm3 + rkpk.r[alpha]*rkpk.r[beta]/lrkpk3)
                            self.H[2*k + alpha, 2*km + beta] += fact * \
                                rkkm.r[alpha]*rkkm.r[beta]/lrkkm3
                            self.H[2*k + alpha, 2*kp + beta] += fact * \
                                rkpk.r[alpha]*rkpk.r[beta]/lrkpk3
                            if alpha == beta:
                                self.H[2*k + alpha, 2*k + beta] += fact * \
                                    (1/lrkkm + 1/lrkpk)
                                self.H[2*k + alpha, 2*km + beta] -= fact/lrkkm
                                self.H[2*k + alpha, 2*kp + beta] -= fact/lrkpk
                    he = he.next
                    if he.idx == first.idx:
                        break

    def __H5(self):
        for f in self.mesh.faces:
            if not f.outer:
                if self.params is None:
                    Lambda = f.params['lambda']
                else:
                    Lambda = self.params[f.type]['lambda']
                he = f.he
                first = f.he
                while True:
                    i = he.vfrom.idx
                    ip = he.vto.idx
                    ripi = he.vto.r - he.vfrom.r
                    lripi = ripi.length()
                    lripi3 = lripi*lripi*lripi
                    f1 = Lambda/lripi3
                    f2 = Lambda/lripi
                    for alpha in range(2):
                        for beta in range(2):
                            self.H[2*i + alpha, 2*i + beta] += f1 * \
                                ripi.r[alpha]*ripi.r[beta]
                            self.H[2*i + alpha, 2*ip + beta] -= f1 * \
                                ripi.r[alpha]*ripi.r[beta]
                            self.H[2*ip + alpha, 2*i + beta] -= f1 * \
                                ripi.r[alpha]*ripi.r[beta]
                            self.H[2*ip + alpha, 2*ip + beta] += f1 * \
                                ripi.r[alpha]*ripi.r[beta]
                            if alpha == beta:
                                self.H[2*i + alpha, 2*i + beta] -= f2
                                self.H[2*i + alpha, 2*ip + beta] += f2
                                self.H[2*ip + alpha, 2*i + beta] += f2
                                self.H[2*ip + alpha, 2*ip + beta] -= f2
                    he = he.next
                    if he.idx == first.idx:
                        break

    def compute(self):
        self.__H13()
        self.__H2()
        self.__H4()
        self.__H5()
        self.evaluated = True

    def eig(self):
        if not self.evaluated:
            raise Exception(
                'Need to compute Hessian before its eigenvalues and eigenvectors can be found.')
        self.eval, self.evec = np.linalg.eigh(self.H)
        self.has_eigvals = True

    def dump_modes(self, filename, modes, draw_periodic=False):
        if not self.has_eigvals:
            raise Exception(
                'Need to compute Hessian eigenvalues and eigenvectors before producing output.')
        points = vtk.vtkPoints()
        faces = vtk.vtkCellArray()
        face = vtk.vtkPolygon()
        vecs = {}
        for m in modes:
            vecs[m] = vtk.vtkDoubleArray()
            vecs[m].SetName("mode_{:06d}".format(m))
            vecs[m].SetNumberOfComponents(3)
        i = 0
        idx = 0
        for v in self.mesh.vertices:
            points.InsertNextPoint([v.r.r[0], v.r.r[1], 0])
            for m in modes:
                vecs[m].InsertNextTuple(
                    [self.evec[i, m], self.evec[i + 1, m], 0])
            i += 2
            idx += 1

        for f in self.mesh.faces:
            omit = False
            crosses_pbc = False
            face_verts = []
            if self.mesh.box != None:
                he = f.he
                first = f.he
                s0 = np.dot(self.mesh.box.invh, he.vfrom.r.r)
                while True:
                    s = np.dot(self.mesh.box.invh, he.vfrom.r.r)
                    ds = s - s0
                    if np.abs(np.rint(ds[0])) != 0 or np.abs(np.rint(ds[1])) != 0:
                        if not draw_periodic:
                            omit = True
                            break
                        else:
                            crosses_pbc = True
                            r = np.dot(self.mesh.box.h, np.array(
                                [s[0] - np.rint(ds[0]), s[1] - np.rint(ds[1])]))
                            points.InsertNextPoint([r[0], r[1], 0])
                            for m in modes:
                                vecs[m].InsertNextTuple(
                                    [self.evec[2*he.vfrom.idx, m], self.evec[2*he.vfrom.idx + 1, m], 0])
                            face_verts.append(idx)
                            idx += 1
                    else:
                        face_verts.append(he.vfrom.idx)
                    he = he.next
                    if he.idx == first.idx:
                        break
            if not (f.outer or omit):
                face.GetPointIds().SetNumberOfIds(f.num_sides())
                he = f.he
                first = f.he
                i = 0
                while True:
                    if not crosses_pbc:
                        face.GetPointIds().SetId(i, he.vfrom.idx)
                    else:
                        face.GetPointIds().SetId(i, face_verts[i])
                    i += 1
                    he = he.next
                    if he.idx == first.idx:
                        break
                faces.InsertNextCell(face)
        polyData = vtk.vtkPolyData()
        polyData.SetPoints(points)
        for m in modes:
            polyData.GetPointData().AddArray(vecs[m])
        polyData.SetPolys(faces)
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(filename)
        if vtk.VTK_MAJOR_VERSION <= 5:
            writer.SetInput(polyData)
        else:
            writer.SetInputData(polyData)
        # writer.SetDataModeToAscii()
        writer.Write()
