# ###############################################################
#
#  Performs simulations of a random patch of particles
#
#  Author: Rastko Sknepnek, (c) 2021
#  Date: 07-Nov-2021
#
# ###############################################################


# ###############################################################
#
# Load standard Python modules
#
# ###############################################################
import sys as s
from os.path import exists
from os import remove
import argparse 
from voro_box import VoroBox
from make_mesh import MakeMesh
from shapely.geometry.polygon import Polygon
import numpy as np

# ###############################################################
#
# Load AJM modules
#
# ###############################################################
from ajm import *
from AJMAnalysis import *

# ###############################################################
#
# Read command line arguments 
#
# ###############################################################
parser = argparse.ArgumentParser()
parser.add_argument('--fpull', dest = 'fpull', type = float, default = 0.05, help = 'pulling force on left and right sides')
parser.add_argument('--beta', dest = 'beta', type = float, default = 0.1, help = 'activity')
parser.add_argument('--tpassive', dest = 'tpassive', type = int, default = 1000 , help = 'number of time units of the passive pull')
parser.add_argument('--tactive', dest = 'tactive', type = int, default = 800 , help = 'number of time units of the active pull')
parser.add_argument('--dt', dest = 'dt', type = float, default = 0.1, help = 'timestep')
parser.add_argument('--seed', dest = 'seed', type = int, default = None, help = 'random number generator seed')
parser.add_argument('--dumpfreq', dest = 'dumpfreq', type = int, default = 100, help = 'how often to produce output')
args = parser.parse_args()

# ###############################################################
#
# Set parameters
#
# ###############################################################
freq = int(round(1.0/args.dt))   # This makes sure that we output data once per unit of time

# Cell types
cell_types = ['active', 'passive']
      
# passive cells and active cells have different elastic and myosin related properties
beta = {} 
beta['passive'] = 0.0 
beta['active']  = args.beta

# sets P0 since lambda = P0*Gamma 
lam = {} 
lam['passive'] = 3.0
lam['active']  = 3.0

#total myosin
myocell = {}
myocell['passive'] = 6.0
myocell['active']  = 6.0

# Cell mechanics parameters
kappa = 1.0     # area stiffness
gamma = 0.5     # perimeter stiffness
k = 0.5         # junction spring constant
myosat = 10.0   # saturation value of myosin (this the parameter 1/alpha in the equation for F(T))
tau_v = 20.0    # viscose relaxation timescale 

# Myosin kinetics parameters
M = 6.0          # total cell myosin (M in eq. 5 in the draft)
myo = 0.5        # myosin reference level m_0  (m_0 in eq. 4 in the draft)
myofluc = 0.01   # mysosin fluctuations (alpha in the paper in eq. 1)
r = 0.01         # myosin associations rate (this is 1/\tau_m)
p = 1.0          # myosin dissociation rate (this is the prefactor that multiplies the exponential term in F(T))
tstar = 0.3      # T^* constant 
k0 = 2.0/tstar   # scaling factor in the exponential 

# ################################################################
#
# Set up simulation objects
#
# ################################################################

if args.seed != None:
    seed = args.seed
else:
    seed = -1        # Use current time as the seed

tissue  = Tissue()                                               # initialise mesh
sim_sys = System(tissue)                                         # base object for the system
forces = Force(sim_sys)                                          # handles all types of forces
integrators = Integrate(sim_sys, forces, seed)                   # handles all integrators
topology = Topology(sim_sys, forces)                             # handles all topology changes (T1, division, ingression)
dumps = Dump(sim_sys, forces)                                    # handles all data output 
simulation = Simulation(sim_sys, integrators, forces, topology)  # simulation object

# #################################################################
#
# Create the initial configuration and read it
#
# #################################################################

input_file = 'beta_{:.4f}_fpull_{:.4f}_seed_{:}.json'.format(args.beta, args.fpull, seed)
if not exists(input_file):
    print('Building random initial configuration...')
    N = 600
    Lx = 40
    Ly = 40
    minimal_interparticle_distance = 1.0
    max_iterations = 20000
    tolerance = 1e-6
    v = VoroBox(N, Lx, Ly, minimal_interparticle_distance)
    v.optimise(tol=tolerance, max_iter=max_iterations)
    # Build the mesh based on the optimal packin
    m = MakeMesh(Lx, Ly)
    seeds = v.random_box.get_coords()
    m.construct_random(seeds)
    m.mark_sides({'left': ('left', True), 'right': ('right', True)})
    m.set_elastic_parameters({'kappa': 2.0, 'gamma': 0.25,
                             'lambda': 1.0, 'k': 0.5, 'l0': 1.0}, stress_free=True)
    m.set_cell_myo(cell_myo=6.0)
    # Set active region
    active_region = [[-18, -18], [18, -18], [18, 18], [-18, 18]]
    m.set_cell_types(active_region, 'active')
    m.json_out(input_file)

sim_sys.read_input(input_file, read_params = True)           # read input configuration


# #################################################################
#
# Add forces to the system
#
# #################################################################

forces.set_update(True)    # enable proper forces summation for detailed output
forces.add('area')         # add area force form term E = 0.5*kappa*(A-A0)^2
forces.add('perimeter')    # add perimeter force term from E = 0.5*gamma*P^2 + lambda*P (maybe -?)
forces.add('active')       # add active term due to myosin on junction
forces.add('harmonic')     # add active term due to myosin on junction

# All force terms will be using paramters defines based on cell type (default is that parateters are given in the JSON file)
forces.set_flag('area', 'use_cell_type')
forces.set_flag('perimeter', 'use_cell_type')
forces.set_flag('active', 'use_cell_type')
forces.set_flag('harmonic', 'use_cell_type')

# Set parameters for each cell type
for tp in cell_types :
   forces.set_params('area', tp , {'kappa' : kappa})
   forces.set_params('perimeter', tp ,  {'gamma': gamma, 'lambda': lam[tp]})
   forces.set_params('harmonic', tp, {'k': k})


# #################################################################
#
# Set conditions for the T1 transition
#
# #################################################################

topology.set_params({'min_edge_len': 0.05, 'new_edge_len': 0.055, 'myosin': 0.5}) 


# #################################################################
#
# Add Brownian integrator that will handle mechanical part
#
# #################################################################

integrators.add('brownian')    

integrators.set_external_force('brownian', 'right', Vec(args.fpull,0.0))  # pulling on the right-most column of vertices
integrators.set_external_force('brownian', 'left', Vec(-args.fpull,0.0))  # pulling on the left-most column of vertices

# #################################################################
#
# Simulation starts here
#
# #################################################################

integrators.set_dt(10*args.dt) # set time step

step = 0       # Step counter in terms of time units

print('Passive pull for {:d} steps'.format(args.tpassive))
# Pulling on the passive system
for i in range(args.tpassive):
    if i % args.dumpfreq == 0:
        dumps.dump_junctions('junctions_fpull_{:.4f}_beta_{:.3f}_step_{:08d}.vtp'.format(args.fpull, args.beta, i))
        dumps.dump_cells('cells_fpull_{:.4f}_beta_{:.3f}_step_{:08d}.vtp'.format(args.fpull, args.beta, i))
    simulation.run(int(round(freq)))
    step += 1 

# Turn on activity
integrators.add('myosin')      # add myosin integrator that will handle myosin dynamics (note: viscoelastic term is missing)

# specify parameters and flags for each integrator
integrators.set_params('myosin', {'r': r, 'p': p, 't*': tstar ,'k0': k0})
integrators.set_flag('myosin', 'fixed_myosin')                   #conserve cell myosin
integrators.set_params('myosin',{"myosin_fluctuation": myofluc}) #myosin noise
integrators.set_params('myosin',{"max_myosin": myosat})          #setting saturated myosin value

sim_sys.set_myomax(0.5)            #sets initial  myosin value
sim_sys.set_symmetric_myosin(True) # active force is of the form beta*(m-m*)

# Add viscoelasticity 
integrators.add('viscoelastic')                           # add Maxwell-like viscoelastic integrator on junctions
integrators.set_params('viscoelastic', {'tau_l': tau_v})  #enabling viscoelasticity

# Set activity ON
for tp in cell_types :  
   forces.set_params('active', tp , {'beta': beta[tp], 'alpha': 0.0})

integrators.set_dt(args.dt) # set time step

#Vfile = 'V_fpull_{:.4f}_beta_{:.3f}_seed_{:d}.dat'.format(args.fpull, args.beta, seed)
Mfile = 'M_fpull_{:.4f}_beta_{:.3f}_seed_{:d}.dat'.format(args.fpull, args.beta, seed)
Mmyofile = 'Mmyo_fpull_{:.4f}_beta_{:.3f}_seed_{:d}.dat'.format(args.fpull, args.beta, seed)
Tensionfile = 'Tens_fpull_{:.4f}_beta_{:.3f}_seed_{:d}.dat'.format(args.fpull, args.beta, seed)

print('Simulating active pull for {:d} steps'.format(args.tactive))
with open(Mfile, 'w') as Mout, open(Mmyofile,'w') as Mmyoout, open(Tensionfile, 'w') as Tout:
    first_iter = True
    for i in range(step, step + args.tactive):
        if i % args.dumpfreq == 0:
            dumps.dump_junctions('junctions_fpull_{:.4f}_beta_{:.3f}_step_{:08d}.vtp'.format(args.fpull, args.beta, i))
            dumps.dump_cells('cells_fpull_{:.4f}_beta_{:.3f}_step_{:08d}.vtp'.format(args.fpull, args.beta, i))
        dumps.dump_mesh('mesh_fpull_{:.4f}_beta_{:.3f}_step_{:08d}.json'.format(args.fpull, args.beta, i))
        if first_iter:
            mesh = utils.HalfEdge.Mesh()
            mesh.read('mesh_fpull_{:.4f}_beta_{:.3f}_step_{:08d}.json'.format(args.fpull, args.beta, i))
            active = []
            for face in mesh.faces:
                if face.type == 'active':
                    active.append(face.idx)
            first_iter = False 
        else:
            #V = Graner.VTensor('mesh_fpull_{:.4f}_beta_{:.3f}_step_{:08d}.json'.format(args.fpull, args.beta, i-1), 'mesh_fpull_{:.4f}_beta_{:.3f}_step_{:08d}.json'.format(args.fpull, args.beta, i))
            M = Graner.MTensor('mesh_fpull_{:.4f}_beta_{:.3f}_step_{:08d}.json'.format(args.fpull, args.beta, i))
            Mmyo = Graner.MyoTensor('mesh_fpull_{:.4f}_beta_{:.3f}_step_{:08d}.json'.format(args.fpull, args.beta, i))
            Tension = Graner.TensionTensor('mesh_fpull_{:.4f}_beta_{:.3f}_step_{:08d}.json'.format(args.fpull, args.beta, i))
            #meanV = np.mean(V.V.T[active,:,:], axis=0).flatten()
            meanM = np.mean(M.M.T[active,:,:], axis=0).flatten()
            meanMmyo = np.mean(Mmyo.Myo.T[active,:,:], axis=0).flatten()
            meanTens = np.mean(Tension.Tension.T[active,:,:], axis=0).flatten() 
            #print('{:d}  {:.16e}  {:.16e}  {:.16e}  {:.16e}'.format(step-args.tpassive, meanV[0], meanV[1], meanV[2], meanV[3]), file=Vout)
            print('{:d}  {:.16e}  {:.16e}  {:.16e}  {:.16e}'.format(step-args.tpassive, meanM[0], meanM[1], meanM[2], meanM[3]), file=Mout)
            print('{:d}  {:.16e}  {:.16e}  {:.16e}  {:.16e}'.format(step-args.tpassive, meanMmyo[0], meanMmyo[1], meanMmyo[2], meanMmyo[3]), file=Mmyoout)
            print('{:d}  {:.16e}  {:.16e}  {:.16e}  {:.16e}'.format(step-args.tpassive, meanTens[0], meanTens[1], meanTens[2], meanTens[3]), file=Tout)
            remove('mesh_fpull_{:.4f}_beta_{:.3f}_step_{:08d}.json'.format(args.fpull, args.beta, i-1))
        simulation.run(int(round(freq)))
        step += 1 

















