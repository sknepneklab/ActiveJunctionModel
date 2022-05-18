from honeycomb_lattice import *

# Generate hexagonal lattice in a rectangular box [-10,10] x [-10,10]
lx = 20.0
ly = 20.0

h = HoneycombLattice(lx, ly, 1.0)

# build lattice 
h.build()
h.minmax()

active = [77, 78, 93, 94]
buffer = [76, 61, 62, 63, 79, 95, 109, 108, 107, 92]

h.set_cell_type(active, 'active')
h.set_cell_type(buffer, 'buffer')

h.set_l0()


# Set a rectangle around the left side of the system
left = [[1.2*h.minx,0.95*h.miny], 
        [0.9*h.minx,0.95*h.miny],
        [0.9*h.minx,0.95*h.maxy],
        [1.2*h.minx,0.95*h.maxy]]

# Set a rectangle around the right side of the system
right = [[0.9*h.maxx,0.95*h.miny], 
         [1.2*h.maxx,0.95*h.miny],
         [1.2*h.maxx,0.95*h.maxy],
         [0.9*h.maxx,0.95*h.maxy]]

# Set boundary vertices on the left as type 1
h.set_vertex_type(left, "left", True)
# Set boundary vertices on the right as type 2
h.set_vertex_type(right, "right", True)

# Set constraint on the left and right types to move only along x-axis
#h.set_constraint('left','x')
#h.set_constraint('right','x')

# Print everything out as JSON
h.json_out("honeycomb_active_buffer.json")