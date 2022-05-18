from honeycomb_lattice import *

# Generate hexagonal lattice in a rectangular box [-20,20] x [-20,20]
lx = 40.0
ly = 40.0

h = HoneycombLattice(lx, ly, 1.0)

# build lattice 
h.build()

# Generate a lace pattern of active 4-cell pieces
widthx = 1.6
widthy = 1.3
distx = 6.0
disty = 5.2

nx = 2*int(lx/(2*distx))
print('nx = ', nx)
ny = 2*int(ly/(2*disty))
print('ny = ', ny)

for k in range(nx + 1):
    for l in range(ny + 1):
        x = distx*(-nx/2 + k)
        y = -0.5 + disty*(-ny/2 + l)
        box = [[x-widthx, y-widthy], [x+widthx, y-widthy],[x+widthx, y+widthy],[x+widthx, y-widthy]]
        h.set_cell_type(box, "active") 

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
h.json_out("honeycomb_lace.json")