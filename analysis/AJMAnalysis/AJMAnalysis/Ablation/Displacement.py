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
# \file Displacement.py
# \author Rastko Sknepnek, sknepnek@gmail.com
# \date 27-Oct-2020
# Calculates radial and tangential displacement after ablation 

from ..utils.HalfEdge import *
import numpy as np
import os.path
import vtk

class Displacement:

    def __init__(self, mesh_1, mesh_2, verts):
        self.init_mash  = Mesh()
        self.final_mash = Mesh()
        self.init_mash.read(mesh_1)
        self.final_mash.read(mesh_2)
        if len(verts) != 2:
            raise ValueError("Two vertices of the ablated junction have to be specified.")
        self.v1 = self.init_mash.vertices[verts[0]]
        self.v2 = self.init_mash.vertices[verts[1]]
        rc = 0.5*(self.v1.r + self.v2.r)
        dr = self.v2.r - self.v1.r
        self.n = dr.unit()
        self.cos = np.zeros(len(self.init_mash.vertices))
        self.sin = np.zeros(len(self.init_mash.vertices))
        self.theta = np.zeros(len(self.init_mash.vertices))
        self.dr = np.zeros(len(self.init_mash.vertices))
        for i in range(len(self.init_mash.vertices)):
            v = self.init_mash.vertices[i]
            ri = v.r - rc 
            theta = np.arctan2(self.n.r[0]*ri.r[1]-self.n.r[1]*ri.r[0], self.n.dot(ri))
            self.theta[i] = theta
            self.cos[i] = np.cos(theta)
            self.sin[i] = np.sin(theta)
            self.dr[i] = ri.length()
    
    def compute_displacement(self):
        self.Dr = np.zeros(len(self.init_mash.vertices))
        self.Dt = np.zeros(len(self.init_mash.vertices))
        self.D = []
        for i in range(len(self.init_mash.vertices)):
            ri = self.init_mash.vertices[i].r
            rf = self.final_mash.vertices[i].r
            D = rf - ri 
            Dx = D.dot(self.n)
            Dy = (D - Dx*self.n).length()
            self.Dr[i] =  self.cos[i]*Dx + self.sin[i]*Dy
            self.Dt[i] = -self.sin[i]*Dx + self.cos[i]*Dy
            self.D.append(D)


    def dump_displacement(self, filename):
        points = vtk.vtkPoints()
        disp = vtk.vtkDoubleArray()
        disp.SetName("displacement")
        disp.SetNumberOfComponents(3)
        for i in range(len(self.init_mash.vertices)):
            v = self.init_mash.vertices[i]
            D = self.D[i]
            points.InsertNextPoint([v.r.r[0], v.r.r[1], 0])
            disp.InsertNextTuple([D.r[0], D.r[1], 0]) 
        polyData = vtk.vtkPolyData()
        polyData.SetPoints(points)
        polyData.GetPointData().AddArray(disp)
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(filename)
        if vtk.VTK_MAJOR_VERSION <= 5:
            writer.SetInput(polyData)
        else:
            writer.SetInputData(polyData)
        writer.Write()
        