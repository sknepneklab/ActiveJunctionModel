from honeycomb_lattice import *

# Genrate hexagonal lattice in a rectangular box [-10,10] x [-10,10]
h = HoneycombLattice(100.0,100.0,1.0)

# build lattice 
h.build(R=20)

# Set vertex types
h.set_vertex_type([17, 22], 'radial')

# Set constraints
h.set_constraint('radial', 'radial')

# record mesh
h.json_out('disk.json')
