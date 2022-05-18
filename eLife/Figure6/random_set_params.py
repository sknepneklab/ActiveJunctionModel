import sys as system 
import argparse
import os 
# AJM import and also condfiguration builder which is in the same structure
system.path.insert(1,'/home/sh18581/Documents/AJM_Quarantine/AJM/build')
from ajm import *
import math
import numpy as np
import bz2

parser = argparse.ArgumentParser(description="""Run passive systems with varying hybrid contributions, p0 and applied force""")
parser.add_argument("--forcex",type=float, default=0.1, help="applied force in x-direction")
parser.add_argument("--forcey",type=float, default=0.0, help="applied force in y-direction")
parser.add_argument("-f", "--freq", type=int, default=200, help="cycle duration (in number of timesteps)")
parser.add_argument("--tauv",type=float, default=20.0, help="viscous time scale")
parser.add_argument("--beta",type=float, default=0.6, help="beta myosin force coupling")
parser.add_argument("--r",type=float, default=0.01, help="myosin on rate")
parser.add_argument("--inifile",type=str, default = "random_40x40_patch.json",  help="name of input file")
parser.add_argument("--dir", type = str, default="./")
parser.add_argument("--graner", action="store_true")
parser.add_argument("--mel", type=float,default="0.05",help="collapse length")
parser.add_argument("--seed", type=int,default=0,help="fluctuation seed")
args = parser.parse_args()


# Ilyas' stuff
def compress_bz2(fname):
   with open(fname,"rb") as file :
      with open('{}.bz2'.format(fname),"wb") as file_compressed :
         data_compressed = bz2.compress(file.read())
         file_compressed.write(data_compressed)

####### Setting simulation parameters
ctypes = ["passive","active"]    
betas={}
betas["passive"] = 0.0
betas["active"] = args.beta

beta_as={}
beta_as["passive"] = 0.0
beta_as["active"] = 0.0
#beta_as["active"] = 0.5*args.beta


myocells={}
myocells["passive"] = 6.0
myocells["active"] = 6.0



# Cell mechanics parameters - fix most
# set in input file now
kappa=1.0
gamma=0.5
alpha=0.0
k=0.5 # harmonic spring constant ( = gamma usually )
tau_l = args.tauv # viscous time scale (start at not viscous)

lams={}
lams["passive"] = 6*gamma
lams["active"] = 6*gamma


# Myosin paramters
#total myosin
cellmyo = 6.0 
# Myosin kinetic parameters
myo = 0.5 
myosat = 10.0 #saturating myosin value
r = args.r #myosin on rate
myofluc = args.r
p = 1.0 # ratio of activaeble to total myosin
tstar = 0.3 # tstar constant np.log(2-1/myosat) / k0 #fixed point value
k0 = 2.0/tstar #scaling factor in the exponential


###### Start of simulation proper
tissue = Tissue()           # initialise mesh
sys = System(tissue)        # base object for the system
f = Force(sys)              # handles all types of forces
integ = Integrate(sys,f,0)  # handles all integrators
mes = Measurement(tissue,sys)
t = Topology(sys, f)        # handles all topolgy changes (T1, division, ingression)
d = Dump(sys, f)            # handles all data output 
sim = Simulation(sys, integ, f, t)  # simulation object

#sys.read_input(args.inifile)           # read input configuration
# If flag True is set, read all the local parameters from json
sys.read_input(args.inifile) 


#----------- Forces 
f.set_update(True) # enable proper forces summation (see commit ef8cd1477c9be0df3c32cb72127dcca936bc6467 on force_compute.hpp)
f.add('area')         # add area force form term E = 0.5*kappa*(A-A0)^2
f.add('perimeter')    # add perimeter force term from E = 0.5*gamma*P^2 + lambda*P (maybe -?)
f.add('harmonic')       # add harmonic part of potential
f.add('active')       # add active term due to myosin on junction
f.add('active_area')       # add active term due to area activity

   

######## For non j-son defined loacl paramters

# All force terms will be using paramters defines based on cell type (default is that parateters are given in the JSON file)
f.set_flag('area', 'use_cell_type')
f.set_flag('perimeter', 'use_cell_type')
f.set_flag('harmonic', 'use_cell_type')

##### This has to always be done locally
f.set_flag('active', 'use_cell_type')
f.set_flag('active_area', 'use_cell_type')


# set paramters for each force term
# general syntax is:
# 1st argument: force term (type of force)
# 2nd argument: cell type 
# 3rd argument: dictionary with the paramter,value pair for each paramter

# Set in input file now
for tp in ctypes :
   f.set_params('area', tp , {'kappa' : kappa})
   f.set_params('perimeter', tp ,  {'gamma': gamma, 'lambda': lams[tp]})
   f.set_params('harmonic', tp, {'k': k})

# DO NOT SET CELL PROPERTIES set (some) cells  properties 
for c in tissue.cells():
   c.property().A0 = 2.598076211353316    # hexagon native area of side length 1 is ~ 2.598 (3/2*\sqrt(3))
   c.property().fa = 0.0 # self-propulsion activity
###########  

#----------- Integrators 
integ.add('brownian')    # add Brownian integrator that handles all vercecx movements
# No viscoelasticity
#integ.add('viscoelastic')# add Maxwell-like viscoelastic integrator on junctions
integ.add('myosin')      # add myosin integrator that will handle myosin dynamics (note: viscoelastic term is missing)

# specify parameters and flags for each integrator
# for myossin
# integ.set_params('myosin', {'kon': kon, 'kmin': kmin, 't*': tstar ,'k0':k0,'Q': Q,'J': J})
integ.set_flag('myosin', 'fixed_myosin')
integ.set_params('myosin', {'t*': tstar
                            ,"r":r,"p":p})
# integ.set_flag('myosin', 'conserve_myosin') # conserve cell myosin
integ.set_params('myosin',{"myosin_fluctuation": myofluc}) #myosin noise
integ.set_params('myosin',{"max_myosin": myosat}) #setting staturating myosin
# integ.set_params('myosin',{"myosin_diffusion": D }) # myosin diffusion within cell
for tp in ctypes : # set total cell myosin per cell types
   integ.set_type_params('myosin', tp , {'cell_myosin': myocells[tp]})
sys.set_myomax(0.5) # sets initial  myosin value
sys.set_symmetric_myosin(True) # active force is of the form beta*(m-m*)

# specify paramters for topology change
t1_resolution_params={'min_edge_len': args.mel, 'new_edge_len': 1.1*args.mel, 'myosin': 0.5}
t.set_params(t1_resolution_params) # /!\ need to set new edges myosin value



#----------------------------------------------------------- We can now stretch -------------------------------------------------------------------------------------
# Start the simulation protocol
step=0
dt_passive=1e-1
dt_active =1e-2
freq_passive=200
freq_active=args.freq
N_active=800
N_passive=50

# Passive part (non-viscous)
print("Stretching passively (build up myosin) - not viscous")
for tp in ctypes : #turn active forces OFF (just to make sure)
   f.set_params('active', tp, {'beta': 0.0, 'alpha': 0.0})
   f.set_params('active_area', tp, {'beta': 0.0})
integ.set_external_force('brownian', 'right', Vec(args.forcex,0))  # pulling on the right-most column of vertices
integ.set_external_force('brownian', 'left', Vec(-args.forcex,0.0))
integ.set_external_force('brownian', 'top', Vec(0,args.forcey))  # pulling on the right-most column of vertices
integ.set_external_force('brownian', 'bottom', Vec(0,-args.forcey))

#no viscoelasticity at the begining
#saving frequency and timestep

integ.set_dt(dt_passive) #set time step
for i in range(step,step+N_passive):
   d.dump_junctions('junctions_{:08d}.vtp'.format(i),True)
   d.dump_cells('cells_{:08d}.vtp'.format(i),True)
   if args.graner == True :
      d.dump_mesh('mesh_{:08d}.json'.format(i)) #dump mesh in json format
      compress_bz2('mesh_{:08d}.json'.format(i))
      os.remove('mesh_{:08d}.json'.format(i))
   sim.run(int(round(freq_passive)))
   step+=1

print("Turning on activity and viscoelasticity")
for tp in ctypes : #turn activity ON
   f.set_params('active', tp , {'beta': betas[tp], 'alpha': alpha})
   f.set_params('active_area', tp , {'beta': beta_as[tp]})

integ.add('viscoelastic')# add Maxwell-like viscoelastic integrator on junctions
integ.set_params('viscoelastic',{'tau_l': tau_l}) #enabling viscoelasticity
integ.set_external_force('brownian', 'right', Vec(args.forcex,0))  # pulling on the right-most column of vertices
integ.set_external_force('brownian', 'left', Vec(-args.forcex,0.0))
integ.set_external_force('brownian', 'top', Vec(0,args.forcey))  # pulling on the right-most column of vertices
integ.set_external_force('brownian', 'bottom', Vec(0,-args.forcey))

#saving frequency and timestep

integ.set_dt(dt_active) #smaller timestep for the active part is generally a good thing
for i in range(step,step+N_active):
   d.dump_junctions('junctions_{:08d}.vtp'.format(i),True)
   d.dump_cells('cells_{:08d}.vtp'.format(i),True)
   if args.graner == True :
      d.dump_mesh('mesh_{:08d}.json'.format(i)) #dump mesh in json format
      compress_bz2('mesh_{:08d}.json'.format(i))
      os.remove('mesh_{:08d}.json'.format(i))
   sim.run(int(round(freq_active)))
   step+=1













