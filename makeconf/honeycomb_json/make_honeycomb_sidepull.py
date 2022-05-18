from honeycomb_lattice import *

# Genrate hexagonal lattice in a rectangular box [-10,10] x [-10,10]
h = HoneycombLattice(20.0,40.0,1.0)

# build lattice 
h.build()

# Generate a strip in the middle that is active
width = 0.0
height = 3.0
centstrip = [[-width, -height], 
        [-width,height],
        [1.1*h.maxx,height],
        [1.1*h.maxx,-height]]

h.set_cell_type(centstrip, "active") 

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


# Set a rectangle around the top side of the system
top = [[1.2*h.minx,0.95*h.maxy], 
       [1.2*h.minx,1.2*h.maxy],
       [1.2*h.maxx,1.2*h.maxy],
       [1.2*h.maxx,0.95*h.maxy]]

# Set a rectangle around the bottom side of the system
bottom = [[1.2*h.minx,0.95*h.miny], 
          [1.2*h.minx,1.2*h.miny],
          [1.2*h.maxx,1.2*h.miny],
          [1.2*h.maxx,0.95*h.miny]]



# Set boundary vertices on the left as type 1
h.set_vertex_type(left, "left", True)
# Set boundary vertices on the right as type 2
h.set_vertex_type(right, "right", True)

h.set_vertex_type(top, "top", True)
h.set_vertex_type(bottom, "bottom", True)

# Set constraint on the left and right types to move only along x-axis
h.set_constraint('left','fixed')
h.set_constraint('right','ux')
h.set_constraint('top','x')
h.set_constraint('bottom','x')

# Print everything out as JSON
h.json_out("honeycomb_sidepull.json")