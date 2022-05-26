###########################################################################
#
#  Copyright (C) 2017, 2018 University of Dundee
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
# \file Tensor.py
# \author Rastko Sknepnek, sknepnek@gmail.com
# \date 10-May-2020
# \brief Stores and processes a sequence of 2x2 matrices
#

import numpy as np
import vtk


class Tensor:

    def __init__(self, N):
        self.N = N
        self.T = np.zeros((N, 2, 2))
        self.eigvals = np.zeros((N, 2))
        self.eigvecs = np.zeros((N, 2, 2))
        self.type = np.zeros(N, dtype=np.int)
        self.has_eigenvals = False

    def compute_eigvals(self):
        for i in range(self.T.shape[0]):
            e, v = np.linalg.eigh(self.T[i, :, :])
            self.eigvals[i, :] = e
            self.eigvecs[i, :, :] = v
            if np.all(e >= 0):
                self.type[i] = 1
            elif np.all(e < 0):
                self.type[i] = 3
            else:
                self.type[i] = 2
        self.compute_eigvals = True

    def mean(self):
        Tavg = np.mean(self.T, axis=0)
        eigval, eigvec = np.linalg.eigh(Tavg)
        vx, vy = eigvec[:, 1]
        alpha = np.arctan2(vy, vx)
        return (eigval, alpha)

    def histogram(self, nbin=25):
        if not self.has_eigenvals:
            self.compute_eigvals()
        anis = self.eigvals[:, 1]/self.eigvals[:, 0]
        disbin = np.linspace(np.min(anis), np.max(anis), nbin + 1)
        anishist, edges = np.histogram(anis, disbin, density=True)
        return anishist, disbin

    def rose(self, nbin=90):
        if not self.has_eigenvals:
            self.compute_eigvals()
        rosebin = np.linspace(-0.5*np.pi, 0.5*np.pi, nbin + 1)
        vx = self.eigvecs[:, 0, 1]
        vy = self.eigvecs[:, 1, 1]
        alpha = np.arctan2(vy, vx)
        hmm = np.where(alpha < -0.5*np.pi)
        alpha[hmm] += np.pi
        hmm2 = np.where(alpha > 0.5*np.pi)
        alpha[hmm2] -= np.pi
        alphahist, edges = np.histogram(alpha, rosebin, normed=True)
        return alphahist, rosebin

    def plot_vtk_tensor(self, filename, mesh, tensorname):
        if not self.has_eigenvals:
            self.compute_eigvals()
        celltypedict = {}
        ct = 1
        for f in mesh.faces:
            if not f.outer:
                if not f.type in celltypedict:
                    celltypedict[f.type] = ct
                    ct += 1
        points = vtk.vtkPoints()
        types = vtk.vtkIntArray()
        celltypes = vtk.vtkIntArray()
        types.SetNumberOfComponents(1)
        types.SetName("EigType")
        celltypes.SetNumberOfComponents(1)
        celltypes.SetName("CellType")
        tensor = vtk.vtkDoubleArray()
        tensor.SetNumberOfComponents(9)
        tensor.SetName(tensorname)
        for f in mesh.faces:
            if not f.outer:
                rc = f.rc(mesh.box)
                points.InsertNextPoint([rc.r[0], rc.r[1], 0.0])
                types.InsertNextValue(self.type[f.idx])
                celltypes.InsertNextValue(celltypedict[f.type])
                tensor.InsertNextTuple([self.T[f.idx, 0, 0], self.T[f.idx, 0, 1],
                                       0.0, self.T[f.idx, 1, 0], self.T[f.idx, 1, 1], 0.0, 0.0, 0.0, 0.0])

        polyData = vtk.vtkPolyData()
        polyData.SetPoints(points)
        polyData.GetPointData().AddArray(types)
        polyData.GetPointData().AddArray(celltypes)
        polyData.GetPointData().AddArray(tensor)
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(filename)
        if vtk.VTK_MAJOR_VERSION <= 5:
            writer.SetInput(polyData)
        else:
            writer.SetInputData(polyData)
        # writer.SetDataModeToAscii()
        writer.Write()
