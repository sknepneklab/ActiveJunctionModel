# AJM python running script to compute the sequential Graner tensors on already existing data
# (c) Silke Henkes and Rastko Sknepnek, 2018
# Aberdeen University, Dundee University

import sys
import glob
import argparse
import numpy as np
import matplotlib.pyplot as plt
import pickle


parser = argparse.ArgumentParser()

parser.add_argument("-path","--path", type=str, default="/home/luke/AJM18/AJM/analysis/Texture/", help="path for ajm files (absolute)")
parser.add_argument("-confdir","--confdir", type=str, default='/home/luke/AJM18/GranerTensors/', help="configuration path (absolute)")
parser.add_argument("-skip","--skip", type=int, default=0, help="number of files to skip in analysis")
parser.add_argument("-nbin","--nbin", type=int, default=35, help="number of bins for the histograms")

args = parser.parse_args()

ajmpath=args.path
#confpath=args.confdir
sys.path.insert(1,ajmpath)
#sys.path.insert(1,confpath)

from texture import Texture

plotting=False
vtkoutput=False
savepickle=True

areafiles = sorted(glob.glob(args.confdir + 'area_*.dat'))[args.skip:]
myosinfiles = sorted(glob.glob(args.confdir + 'myosin_*.dat'))[args.skip:]

nfiles=len(areafiles)

p0bin=np.linspace(3.5,5.5,args.nbin+1)
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
    

# locate the files to carry out the Graner computations
texturefiles = sorted(glob.glob(args.confdir + 'texture_*.dat'))[args.skip:]
nfiles=len(texturefiles)

u=0
meigval=np.zeros((nfiles,2))
alpha=np.zeros((nfiles,))
Lx=np.zeros((nfiles,))
Ly=np.zeros((nfiles,))

anisdis=np.zeros((nfiles,args.nbin))
rosedis=np.zeros((nfiles,args.nbin))
for tx in texturefiles:
    #print tx
    txtr = Texture(tx)
    txtr.compute_eigvals()
    m,a = txtr.get_MeanTexture()
    Lx[u],Ly[u]=txtr.get_Dimensions()
    anisdis[u,:],anibin=txtr.get_Distribution(args.nbin)
    rosedis[u,:],rosebin=txtr.get_Rose(args.nbin)
    # note nematic symmetry
    if a<-np.pi/2:
        a+=np.pi
    if a>np.pi/2:
        a-=np.pi
    meigval[u,:],alpha[u] = m,a
    if vtkoutput:
        txtr.compute_ellipses(10,0.2)
        txtvtk = args.confdir + 'texture' + str(u) + '.vtp' 
        txtr.plot_vtk(txtvtk)
    u+=1
    
# Aaand ... the same thing for the myosin 
myofiles = sorted(glob.glob(args.confdir + 'myotexture_*.dat'))[args.skip:]
nfiles=len(myofiles)

u=0
myoeigval=np.zeros((nfiles,2))
malpha=np.zeros((nfiles,))
manisdis=np.zeros((nfiles,args.nbin))
mrosedis=np.zeros((nfiles,args.nbin))
for tx in myofiles:
    #print tx
    txtr = Texture(tx)
    txtr.compute_eigvals()
    
    m,a = txtr.get_MeanTexture()
    manisdis[u,:],anibin=txtr.get_Distribution(args.nbin)
    mrosedis[u,:],rosebin=txtr.get_Rose(args.nbin)
    # note nematic symmetry
    if a<-np.pi/2:
        a+=np.pi
    if a>np.pi/2:
        a-=np.pi
    myoeigval[u,:],malpha[u] = m,a
    if vtkoutput:
        txtr.compute_ellipses(10,1.0)
        txtvtk = args.confdir + 'myotexture' + str(u) + '.vtp' 
        txtr.plot_vtk(txtvtk)
    u+=1

  
timeval=np.linspace(0,nfiles,nfiles)
    
if savepickle:
    filename=args.confdir +'statistics.p'
    data={'confdir':args.confdir,'skip':args.skip,'p0bin':p0bin,'p0hist':p0hist,'myobin':myobin,'myohist':myohist,'timeval':timeval,'meigval':meigval,'alpha':alpha,'Lx':Lx,'Ly':Ly,'myoeigval':myoeigval,'malpha':malpha,'anibin':anibin,'anisdis':anisdis,'manisdis':manisdis,'rosebin':rosebin,'rosedis':rosedis,'mrosedis':mrosedis}
    pickle.dump(data,open(filename,'wb'))
                
if plotting:
    
    plt.figure()
    plt.plot(timeval,meigval[:,0],'-k',label='min eigval',lw=2)
    plt.plot(timeval,meigval[:,1],'-r',label='max eigval',lw=2)
    plt.plot(timeval,alpha,'-b',label='angle',lw=2)
    plt.xlabel('time')
    plt.ylabel('mean texture')
    plt.legend()

    plt.figure()
    plt.plot(timeval,Lx,'-k',label='Lx',lw=2)
    plt.plot(timeval,Ly,'-r',label='Ly',lw=2)
    plt.xlabel('time')
    plt.ylabel('size')
    plt.legend()

    plt.figure()
    plt.plot(timeval,myoeigval[:,0],'-k',label='min eigval',lw=2)
    plt.plot(timeval,myoeigval[:,1],'-r',label='max eigval',lw=2)
    plt.plot(timeval,malpha,'-b',label='angle',lw=2)
    plt.xlabel('time')
    plt.ylabel('mean myosin texture')
    plt.legend()


    da=anibin[1]-anibin[0]

    plt.figure()
    for u in range(0,nfiles,10):
        plt.plot(anibin[1:]-da/2.0,anisdis[u,:],'.-')
    plt.xlabel('anisotropy')
    plt.title('Texture distribution')
        
    plt.figure()
    for u in range(0,nfiles,10):
        plt.plot(anibin[1:]-da/2.0,manisdis[u,:],'.-')
    plt.xlabel('anisotropy')
    plt.title('Myosin texture distribution')

    dr=rosebin[1]-rosebin[0]
    
    plt.figure()
    for u in range(0,nfiles,10):
        plt.plot(rosebin[1:]-dr/2.0,rosedis[u,:],'.-')
    plt.xlabel('angle (radians)')
    plt.title('Texture angle distribution')
    
    plt.figure()
    for u in range(0,nfiles,10):
        plt.plot(rosebin[1:]-dr/2.0,mrosedis[u,:],'.-')
    plt.xlabel('angle (radians)')
    plt.title('Myosin angle distribution')
        
    plt.figure()
    for u in range(0,nfiles,10):
        plt.plot(p0bin[1:]-dp0/2.0,p0hist[u,:],'.-')
    plt.xlabel('p0')
    plt.title('p0 distribution')
        
    plt.figure()
    for u in range(0,nfiles,10):
        plt.plot(myobin[1:]-dmy/2.0,myohist[u,:],'.-')
    plt.xlabel('myosin')
    plt.title('Myosin distribution')



    plt.show()
    
    