from honeycomb_lattice import *

# Genrate hexagonal lattice in a rectangular box [-10,10] x [-10,10]
h = HoneycombLattice(60.0,60.0,1.0)

# build lattice 
h.build()
h.minmax()

# mark top row of vertices
top =  [[1.1*h.minx,0.95*h.maxy], 
         [1.1*h.minx,1.2*h.maxy],
         [1.1*h.maxx,1.2*h.maxy],
         [1.1*h.maxx,0.95*h.maxy]]

# mark bottom row of vertices
bottom  =  [[1.1*h.minx,0.95*h.miny], 
         [1.1*h.minx,1.2*h.miny],
         [1.1*h.maxx,1.2*h.miny],
         [1.1*h.maxx,0.95*h.miny]]


# Set type of top vertices
h.set_vertex_type(top, "top")
# Set type of bottom vertices
h.set_vertex_type(bottom, "bottom")

# Fix top and bottom vertices, so forces don't act on them
h.set_constraint('top','fixed')
h.set_constraint('bottom','fixed')

# Set types oce top and bottom cells
h.set_cell_type(h.get_cells_with_vert_type('top'), 'top')
h.set_cell_type(h.get_cells_with_vert_type('bottom'), 'bottom')

# Print everything out as JSON
h.json_out("large_hexagonal_shear.json")
