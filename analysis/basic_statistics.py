# AJM python running script to compute basic statistics, including p0 and myosin distributions
# (c) Silke Henkes and Rastko Sknepnek, 2018
# Aberdeen University, Dundee University

import sys
import glob
import argparse
import numpy as np
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser()

parser.add_argument("-path","--path", type=str, default="/home/luke/AJM18/AJM/analysis/Texture/", help="path for ajm files (absolute)")
parser.add_argument("-confdir","--confdir", type=str, default='/home/luke/AJM18/GranerTensors/', help="configuration path (absolute)")
parser.add_argument("-skip","--skip", type=int, default=0, help="number of files to skip in analysis")
parser.add_argument("-nbin","--nbin", type=int, default=25, help="number of bins for the histograms")

args = parser.parse_args()

ajmpath=args.path
sys.path.insert(1,ajmpath)

areafiles = sorted(glob.glob(args.confdir + 'area_*.dat'))[args.skip:]
myosinfiles = sorted(glob.glob(args.confdir + 'myosin_*.dat'))[args.skip:]

nfiles=len(areafiles)

p0bin=np.linspace(3.5,4.5,args.nbin+1)
dp0=p0bin[2]-p0bin[1]
myobin=np.linspace(0.0,1.25,args.nbin+1)
dmy=myobin[1]=myobin[0]

p0hist=np.zeros((nfiles,args.nbin))
myohist=np.zeros((nfiles,args.nbin))
u=0
for ar in areafiles:
    data = np.loadtxt(ar)
    area = data[:,0]
    perimeter = data[:,1]
    p0hist[u,:],edges=np.histogram(perimeter/np.sqrt(area),bins=p0bin,normed=True)
    u+=1
    
u=0
for my in myosinfiles:
    myosin = np.loadtxt(my)
    myohist[u,:],edges=np.histogram(myosin,bins=myobin,normed=True)
    u+=1
    
    
plt.figure()
for u in range(0,nfiles,10):
    plt.plot(p0bin[1:]-dp0/2.0,p0hist[u,:],'.-')
    
plt.figure()
for u in range(0,nfiles,10):
    plt.plot(myobin[1:]-dmy/2.0,myohist[u,:],'.-')



plt.show()
    