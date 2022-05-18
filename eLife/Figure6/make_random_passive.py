import sys as system 
system.path.insert(1,'/home/sh18581/Documents/AJM_Quarantine/AJM/makeconf/box')

from voro_box import VoroBox
from make_mesh import MakeMesh
from shapely.geometry.polygon import Polygon
import numpy as np

Lx = 40
Ly = 40
avarea = 2.598076211353316
N = int(Lx*Ly/avarea)
print(N)
minimal_interparticle_distance = 1.0
# Lower to make more disordered packings
max_iterations = 2000
# Stopping criterion
tolerance = 5e-5
optimised_coords = 'seeds.dat'
from_file = False  # Set this to True to build from a relaxed configuration
output = 'random_40x40_patch.json'

# Optimise random cell packing
if not from_file:
  v = VoroBox(N, Lx, Ly, minimal_interparticle_distance)
  v.optimise(tol = tolerance, max_iter= max_iterations)
  np.savetxt(optimised_coords, v.random_box.get_coords())

# Build the mesh based on the optimal packing
m = MakeMesh(Lx, Ly)
if from_file:
  seeds = np.loadtxt(optimised_coords)
else:
  seeds = v.random_box.get_coords()

m.construct_random(seeds)
m.mark_sides({'left' : ('left', False), 'bottom': ('bottom', False), 'right': ('right', False), 'top': ('top', False)})

m.set_elastic_parameters({'kappa': 1.0, 'gamma': 0.5, 'k': 0.5})
m.set_cell_myo(cell_myo = 6.0)

# Set active region
active_region = [[-18,-18], [18,-18], [18,18], [-18,18]]
m.set_cell_types(active_region, 'active')
m.json_out(output)


