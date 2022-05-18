from triangular_lattice import *
from make_mesh import *

t = TriangularLattice(21)
t.build_periodic_lattice()
m = MakeMesh(t)
m.make_initial_configuration()
A0 = m.polygon_area(m.mesh_faces[0], m.mesh_vertices)
m.json_out("hexagonal_grow.json", params={'kappa': 1.0, 'gamma': 0.25, 'maxA0': 1.5*A0})
