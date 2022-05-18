from honeycomb_lattice import *

# Genrate hexagonal lattice in a penatog
h = HoneycombLattice(200.0, 200.0, 1.0)

# build lattice 
h.build(poly=[[25,0],[8,-24],[-20,-15],[-20,15],[8,24]]) 

# record mesh
h.json_out('pentagon.json')
