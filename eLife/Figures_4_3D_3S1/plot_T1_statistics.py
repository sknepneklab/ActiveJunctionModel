import numpy as np  
import pickle
import matplotlib.pyplot as plt


# Seeds for our simulations (maximum 32)
seeds = range(1,33)

# Values of active myosin coupling connstant
betaval=['0.0','0.2','0.4','0.6','0.8','1.0','1.2','1.4']
# Values of external pulling force
fval=['0.0','0.025','0.05','0.075','0.1','0.125','0.15','0.175','0.2','0.225','0.25','0.275','0.3']

# Clunky but correct: actual number of seeds run. Note that for the nseeds = 16, no actual T1s happen at all
# and simulation is very nearly deterministic. N = 16 is plenty for averaging here.
nseeds = np.zeros((len(betaval),len(fval))).astype(int)
# and then just list them
nseeds[0,:] = np.array([16,16,16,16,16,16,16,16,16,16,16,16,16])
nseeds[1,:] = np.array([16,16,16,16,16,16,16,16,16,16,16,16,16])
nseeds[2,:] = np.array([16,16,16,16,16,16,16,16,16,16,16,16,16])
nseeds[3,:] = np.array([16,16,16,16,32,32,32,32,16,16,16,16,16])
nseeds[4,:] = np.array([16,16,16,32,32,32,32,32,32,32,16,16,16])
nseeds[5,:] = np.array([16,16,32,32,32,32,32,32,32,32,32,32,32])
nseeds[6,:] = np.array([16,32,32,32,32,32,32,32,32,32,32,32,32])
nseeds[7,:] = np.array([32,32,32,32,32,32,32,32,32,32,32,32,32])

# We now have a handful of crashed simulations, all at beta 1.4 
# label: bet idx, f idx, seed idx
badpts = [[7,2,27],[7,2,28],[7,3,26],[7,3,29],[7,9,28]]
nbad=np.zeros((len(betaval),len(fval)))
nbad[7,:]=[0,0,2,2,0,0,0,0,0,1,0,0,0]

# All simulations are 899 time units long, with a spacing of 1 time unit between snapshots
# after the first 100 snapshots (equilibration). We discard those.
Nsnap = 899
dt = 1

# Graner coarse-grained U, P and V tensors for the 14-cell central active region
Uval = np.zeros((len(betaval),len(fval),len(seeds),Nsnap,2,2))
Pval = np.zeros((len(betaval),len(fval),len(seeds),Nsnap-1,2,2))
Vval = np.zeros((len(betaval),len(fval),len(seeds),Nsnap-1,2,2))

# Mean tension and myosin for the central acive region
Tensavg = np.zeros((len(betaval),len(fval),len(seeds),Nsnap,2,2))
Myoavg = np.zeros((len(betaval),len(fval),len(seeds),Nsnap,2,2))

# Statistics on area, perimeter and shape factor
avArea = np.zeros((len(betaval),len(fval),len(seeds),Nsnap))
avPerim = np.zeros((len(betaval),len(fval),len(seeds),Nsnap))
avp0 = np.zeros((len(betaval),len(fval),len(seeds),Nsnap))
nbin = 40

# Distributions of shape parameter, myosin and tension
p0hist = np.zeros((len(betaval),len(fval),len(seeds),Nsnap,nbin))
myohist = np.zeros((len(betaval),len(fval),len(seeds),Nsnap,nbin))
tenshist = np.zeros((len(betaval),len(fval),len(seeds),Nsnap,nbin))

# total number of T1s in the simulation (from Graner statistics)
NT1tot = np.zeros((len(betaval),len(fval),len(seeds)))
# for the angle histograms. This appears to be -pi to pi indiscriminately
# wrap in half
nbin =11
anglebin = np.linspace(0,np.pi,nbin+1)
# make a second one for the first T1, with finer spacing
nbinfirst =59
anglebinfirst = np.linspace(0,np.pi,nbinfirst+1)
anglebinfirst_VAR = np.linspace(-np.pi/2,np.pi/2,nbinfirst+1)

# Histogram of distributions of T1s, separated into 'all' and 'first'
T1plus_all_hist = np.zeros((len(betaval),len(fval),nbin))
T1plus_first_hist = np.zeros((len(betaval),len(fval),nbinfirst))
T1min_all_hist = np.zeros((len(betaval),len(fval),nbin))
T1min_first_hist = np.zeros((len(betaval),len(fval),nbinfirst))
T1min_first_hist_VAR = np.zeros((len(betaval),len(fval),nbinfirst))

# first T1, and second T1, from Graner statistics (i.e. not necessarily central)
isT1 = np.zeros((len(betaval),len(fval),len(seeds)))
isT12 = np.zeros((len(betaval),len(fval),len(seeds)))
# start of active phase in frames
activestart = 100
firstT1time = -1*np.ones((len(betaval),len(fval),len(seeds)))
secondT1time = -1*np.ones((len(betaval),len(fval),len(seeds)))

# Dynamical T1 data: From the central junction analysis (tensor_timeseries pickle files)
# These existence and times of first and second T1 are for the central T1, ultimately what we chose to show
isT1_Dyn = np.zeros((len(betaval),len(fval),len(seeds)))
isT12_Dyn = np.zeros((len(betaval),len(fval),len(seeds)))

firstT1time_Dyn = -1*np.ones((len(betaval),len(fval),len(seeds)))
secondT1time_Dyn = -1*np.ones((len(betaval),len(fval),len(seeds)))

# We now have a handful of crashed simulations, all at beta 1.4 
# label: bet idx, f idx, seed idx
badpts = [[7,2,27],[7,2,28],[7,3,26],[7,3,29],[7,9,28]]

# Read in all of our pickle files and fill in statistics above
for b in range(len(betaval)):
    for f in range(len(fval)):
        # collecting Graner T1 statistics from time series. 'f' is for first, 'a' is for all, 'p' is for appearing  junction, 'm' is for diseappearing junction
        T1pflist = []
        T1mflist = []
        T1palist = []
        T1malist = []
        for s in range(nseeds[b,f]):
            if not [b,f,seeds[s]] in badpts:

                ##### Central junction data #######
                # Dynamics information from the bespoke analysis of the central junctions
                picklefile2 = './picklefiles2/beta' + betaval[b] +'/fpull' + fval[f]+'/junction_dynamics_timeseries_' + str(seeds[s]) + '.p'
                print(picklefile2)
                # data.update({'isdata':isdata,'firstT1':firstT1,'secondT1':secondT1,'myosins':myosins,'tensions':tensions,'l0s':l0s,'lengths':lengths})
                data2  = pickle.load(open(picklefile2,'rb'))
                firstT1 = data2['firstT1']
                secondT1 = data2['secondT1']

                # Existence and time of first and second central T1. Take away equilibration period from time count.
                if firstT1>0:
                    isT1_Dyn[b,f,s]=1
                    firstT1time_Dyn[b,f,s]=tval[firstT1]-tval[activestart]
                if secondT1>0:
                    isT12_Dyn[b,f,s]=1
                    secondT1time_Dyn[b,f,s]=tval[secondT1]-tval[activestart]

                ######## Graner statistics averaged data #####
                # Graner statistics output pickle files
                picklefile = './picklefiles2/beta' + betaval[b] + '/fpull' + fval[f]+'/tensor_timeseries_' + str(seeds[s]) + '.p'
                print(picklefile)
                data  = pickle.load(open(picklefile,'rb'))
                tval = data['tval']
                
                # Averaged tensors
                Uval[b,f,s,:,:,:] = data["U"]
                Pval[b,f,s,:,:,:] = data["P"]
                Vval[b,f,s,:,:,:] = data["V"]

                # T1 data
                nt1 = data['NT1']
                t1neg = data['T1angle_pos']
                t1pos = data['T1angle_neg']

                # Every junction counts 4 times (2 disappear from each cell, 2 appear). This is not an integer if one of the cells
                # involved is part of the boundary
                NT1tot[b,f,s] = np.sum(nt1)/4.0

                # T1 statistics from averaged data
                # locate first two nonzero elements
                hasT1 = np.nonzero(nt1)[0]
                if len(hasT1)>0:
                    firstT1time[b,f,s]=tval[hasT1[0]]-tval[activestart]
                    # might fail since it could be a junction half in half out
                    try:
                        T1pflist.append(t1pos[0])
                        T1mflist.append(t1neg[0])
                        isT1[b,f,s]=1
                    except:
                        pass
                if len(hasT1)>1:
                    secondT1time[b,f,s]=tval[hasT1[1]]-tval[activestart]
                    isT12[b,f,s]=1
                # put all of the angles into the running list of T1s
                T1palist.extend(t1pos)
                T1malist.extend(t1neg)

                # We are collecting all the tissue statistics here
                Tensavg[b,f,s,:,:,:] = data["Tensavg"]
                Myoavg[b,f,s,:,:,:] = data["Myoavg"]

                # Global statistics
                avArea[b,f,s,:] = data['avArea']
                avPerim[b,f,s,:] = data['avPerim']
                avp0[b,f,s,:] = data['avp0']

                # p0, myosin and tension distributions
                # bins are always the same, simply replace
                p0bin = data['p0bin']
                p0hist[b,f,s,:,:]=data['p0hist']

                myobin = data['myobin']
                myohist[b,f,s,:,:]=data['myohist']

                tensbin = data['tensbin']
                tenshist[b,f,s,:,:]=data['tenshist']

        # wrap T1 angles into upper two quadrants
        # [f(x) if condition else g(x) for x in sequence]
        T1pflist2 = [x if x>0 else (x+np.pi)  for x in T1pflist]
        T1mflist2 = [x if x>0 else (x+np.pi)  for x in T1mflist]
        #T1mflist2_VAR = [x if x>0 else (x+np.pi)  for x in T1mflist]
        T1palist2 = [x if x>0 else (x+np.pi)  for x in T1palist]
        T1malist2 = [x if x>0 else (x+np.pi)  for x in T1malist]
        print(T1mflist)
        print(T1mflist2)
    
        # Actual histograms. Note the 'VAR' is simply playing around with the binning
        T1plus_first_hist[b,f,:],edges = np.histogram(T1pflist2,bins=anglebinfirst)
        T1min_first_hist[b,f,:],edges = np.histogram(T1mflist2,bins=anglebinfirst)
        T1min_first_hist_VAR[b,f,:],edges = np.histogram(T1mflist,bins=anglebinfirst_VAR)
        T1plus_all_hist[b,f,:],edges = np.histogram(T1palist2,bins=anglebin)
        T1min_all_hist[b,f,:],edges = np.histogram(T1malist2,bins=anglebin)

        # Normalisation: Divide by number of configurations, except for the 'all' counts where every T1 is in the list twice
        T1plus_first_hist[b,f,:]=T1plus_first_hist[b,f,:]/(nseeds[b,f]-nbad[b,f])
        T1min_first_hist[b,f,:]=T1min_first_hist[b,f,:]/(nseeds[b,f]-nbad[b,f])
        T1min_first_hist_VAR[b,f,:]=T1min_first_hist_VAR[b,f,:]/(nseeds[b,f]-nbad[b,f])
        T1plus_all_hist[b,f,:]=T1plus_all_hist[b,f,:]/(2*(nseeds[b,f]-nbad[b,f]))
        T1min_all_hist[b,f,:]=T1min_all_hist[b,f,:]/(2*(nseeds[b,f]-nbad[b,f]))


# Constructing appropriate edges for the two dimensional colour phase diagram plots 
fs = []
for f in range(len(fval)):
    fs.append(float(fval[f]))

fedges=[]
df = fs[1]-fs[0]
for f in range(len(fval)):
    fedges.append(fs[f]-df/2.0)
fedges.append(fs[f]+df/2.0)

betas = np.zeros((len(betaval),))
for b in range(len(betaval)):
    betas[b]=float(betaval[b])

betaedges=[]
db = betas[1]-betas[0]
for b in range(len(betaval)):
    betaedges.append(betas[b]-db/2.0)
betaedges.append(betas[b]+db/2.0)

fs=np.array(fs)
betas=np.array(betas)


# raw Number of T1s as a checkerboard plot
plt.figure()
T1data = np.average(NT1tot[:,:,:],axis=2)
plt.pcolormesh(fedges,betaedges,T1data,vmin=0,vmax=50)
plt.colorbar()
plt.xlabel('pulling force')
plt.ylabel(r'$\beta$')
plt.title('Number of T1s') 

# probability of having a first T1 (a la Ilyas), i.e. any T1 at all. 
# Resolution set by number of seeds for that point in parameter space
plt.figure()
T1plot = np.sum(isT1[:,:,:],axis=2)/len(seeds)
plt.pcolormesh(fedges,betaedges,T1data,vmin=0,vmax=1)
plt.colorbar()
plt.xlabel('pulling force')
plt.ylabel(r'$\beta$')
plt.title('Probability of first T1') 

### In publication
# Dynamics data: Probability of having a first T1 being on the central junction
plt.figure()
T1plot_Dyn = np.sum(isT1_Dyn[:,:,:],axis=2)/len(seeds)
plt.pcolormesh(fedges,betaedges,T1plot_Dyn,vmin=0,vmax=1)
plt.colorbar()
# Including the 25% contour i the probability of any T1 at that point.
plt.contour(fs,betas,T1plot,[0.25],color='red',linewidths=2)

plt.xlabel('pulling force')
plt.ylabel(r'$\beta$')
plt.title('Probability of first T1 (junction method)')

### In Publication
# Cleaned up first T1 transition plot, only where at least 25% of first T1s  happen
# Final version: Error bar is error of the mean, i.e. divided by sqrt(number of samples)
# Data point is the mean (not median)
plt.figure()
contraction = np.zeros((len(betaval),len(fval)))
d_contraction = np.zeros((len(betaval),len(fval)))
expansion = np.zeros((len(betaval),len(fval)))
d_expansion = np.zeros((len(betaval),len(fval)))
color=plt.cm.nipy_spectral(np.linspace(0,1,len(betaval)+1))
for b,c in zip(range(len(betaval)),color):
    isdata = []
    for f in range(len(fval)):
        if sum(isT1_Dyn[b,f,:])>0.25*nseeds[b,f]:
            isdata.append(f)
            hasT1=np.nonzero(isT1_Dyn[b,f,:])[0]
            contraction[b,f]=np.average(firstT1time_Dyn[b,f,hasT1])
            #contraction[b,f]=np.median(firstT1time_Dyn[b,f,hasT1])
            d_contraction[b,f]=np.std(firstT1time_Dyn[b,f,hasT1])/np.sqrt(len(hasT1))
    plt.errorbar(fs[isdata],contraction[b,isdata],d_contraction[b,isdata],marker='o',markersize=8,lw=2,linestyle='--',color=c,label=betaval[b])
plt.xlabel('pulling force')
plt.ylabel('time') 
plt.title('first T1 time (contraction time) - junction method') 
plt.ylim(0,600)
plt.xlim(0.0,0.325)
plt.legend()

# First and second T1 transitions from Graner data only (i.e. not controlled for being central)
# time of first and second T1s, and hence contraction and expansion phase
# be a bit more careful, there are places where none happen. Don't average or plot over those
contraction = np.zeros((len(betaval),len(fval)))
d_contraction = np.zeros((len(betaval),len(fval)))
expansion = np.zeros((len(betaval),len(fval)))
d_expansion = np.zeros((len(betaval),len(fval)))
color=plt.cm.nipy_spectral(np.linspace(0,1,len(betaval)+1))
plt.figure()
for b,c in zip(range(len(betaval)),color):
    isdata = []
    for f in range(len(fval)):
        if sum(isT1[b,f,:])>0:
            isdata.append(f)
            hasT1=np.nonzero(isT1[b,f,:])[0]
            contraction[b,f]=np.average(firstT1time[b,f,hasT1])
            d_contraction[b,f]=np.std(firstT1time[b,f,hasT1])/np.sqrt(len(hasT1))
    plt.errorbar(fs[isdata],contraction[b,isdata],d_contraction[b,isdata],marker='o',markersize=8,linestyle='--',color=c,label=betaval[b])
plt.xlabel('pulling force')
plt.ylabel('time') 
plt.title('first T1 time (contraction time)') 
plt.ylim(0,600)
plt.legend()

plt.figure()
# if there is a second T1, there was a first T1
for b,c in zip(range(len(betaval)),color):
    isdata = []
    for f in range(len(fval)):
        if sum(isT12[b,f,:])>0:
            isdata.append(f)
            hasT1=np.nonzero(isT12[b,f,:])[0]
            expansion[b,f]=np.average(secondT1time[b,f,hasT1]-firstT1time[b,f,hasT1])
            d_expansion[b,f]=np.std(secondT1time[b,f,hasT1]-firstT1time[b,f,hasT1])/np.sqrt(len(hasT1))
    plt.errorbar(fs[isdata],expansion[b,isdata],d_expansion[b,isdata],marker='o',markersize=8,linestyle='--',color=c,label=betaval[b])
plt.xlabel('pulling force')
plt.ylabel('time') 
plt.title('2nd T1 - 1st T1 time (expansion time)') 
plt.ylim(0,600)
plt.legend()


# Graner statistics T1 orientation plots. Does include non-central T1s for the first T1 transition
# Switch to turn off as it's a lot of plots
plotT1orient=False
if plotT1orient:
    ## Supplementary?
    ## For all T1s (these are reasonably nice)

    ## Line plots
    # for f in range(len(fval)):
    #     plt.figure()
    #     color2=plt.cm.nipy_spectral(np.linspace(0,1,len(betaval)+1))
    #     for b,c in zip(range(len(betaval)),color2):
    #         if sum(isT1[b,f,:])>0:
    #             dang=anglebin[1]-anglebin[0]
    #             plt.plot(anglebin[1:]-dang/2.0,T1plus_all_hist[b,f,:],marker='o',linestyle='-',color=c,label=betaval[b])
    #             plt.plot(anglebin[1:]-dang/2.0,-T1min_all_hist[b,f,:],marker='o',linestyle='-',color=c)
    #     plt.xlabel('angle')
    #     plt.ylabel('All T1 orientation')
    #     plt.title('force ' + fval[f])
    #     plt.legend()

    ## Polar histograms (fiddly, still not convincing)
    # for f in range(len(fval)):
    #     plt.figure()
    #     ax = plt.subplot(111, projection='polar')
    #     color2=plt.cm.nipy_spectral(np.linspace(0,1,len(betaval)+1))
    #     for b,c in zip(range(len(betaval)),color2):
    #         if sum(isT1[b,f,:])>0:
    #             dang=anglebin[1]-anglebin[0]
    #             anglebin2 = anglebin[1:]-dang/2.0+np.pi
    #             #plt.plot(anglebin[1:]-dang/2.0,T1plus_all_hist[b,f,:],marker='.',linestyle='-',color=c,label=betaval[b])
    #             #plt.plot(anglebin2,T1min_all_hist[b,f,:],marker='.',linestyle='-',color=c)
    #             plt.bar(anglebin[1:]-dang/2.0,T1plus_all_hist[b,f,:],width=0.27,lw=2,edgecolor=c,color=c,alpha=0.75,label=betaval[b],zorder=len(betaval)-b)
    #             plt.bar(anglebin2,T1min_all_hist[b,f,:],width=0.27,lw=2,edgecolor=c,color=c,alpha=0.75,zorder=len(betaval)-b)
    #     plt.xlabel('angle')
    #     plt.ylabel('All T1 orientation')
    #     plt.title('force ' + fval[f])
    #     plt.legend()

    ##  First T1 (Graner)
    # Line plots
    for f in range(len(fval)):
        plt.figure()
        color2=plt.cm.nipy_spectral(np.linspace(0,1,len(betaval)+1))
        for b,c in zip(range(len(betaval)),color2):
            if sum(isT1[b,f,:])>0:
                dang=anglebinfirst[1]-anglebinfirst[0]
                plt.plot(anglebinfirst[1:]-dang/2.0,T1plus_first_hist[b,f,:],marker='.',linestyle='-',color=c,label=betaval[b])
                #plt.plot(anglebinfirst[1:]-dang/2.0,-T1min_first_hist[b,f,:],marker='.',linestyle='-',color=c)
                plt.plot(anglebinfirst_VAR[1:]-dang/2.0,-T1min_first_hist_VAR[b,f,:],marker='.',linestyle='-',color=c)
        plt.xlabel('angle')
        plt.ylabel('First T1 orientation')
        plt.title('force ' + fval[f])
        plt.legend()

    # Polar plots
    # for f in range(len(fval)):
    #     plt.figure()
    #     ax = plt.subplot(111, projection='polar')
    #     color2=plt.cm.nipy_spectral(np.linspace(0,1,len(betaval)+1))
    #     sum1 = np.zeros((len(anglebinfirst[1:]),2))
    #     for b,c in zip(range(len(betaval)),color2):
    #         if sum(isT1[b,f,:])>0:
    #             dang=anglebinfirst[1]-anglebinfirst[0]
    #             anglebin2 = anglebinfirst[1:]-dang/2.0+np.pi
    #             plt.bar(anglebinfirst[1:]-dang/2.0,sum1[:,0] + T1plus_first_hist[b,f,:],width=0.03,lw=0,edgecolor=c,color=c,alpha=1,label=betaval[b],zorder=len(betaval)-b)
    #             plt.bar(anglebin2,sum1[:,0] + T1plus_first_hist[b,f,:],width=0.03,lw=0,edgecolor=c,color=c,alpha=1,zorder=len(betaval)-b)
    #             plt.bar(anglebin2,sum1[:,1] + T1min_first_hist[b,f,:],width=0.03,lw=0,edgecolor=c,color=c,alpha=1,zorder=len(betaval)-b)
    #             plt.bar(anglebinfirst[1:]-dang/2.0,sum1[:,1] + T1min_first_hist[b,f,:],width=0.03,lw=0,edgecolor=c,color=c,alpha=1,zorder=len(betaval)-b)
    #             sum1[:,0] += T1plus_first_hist[b,f,:]
    #             sum1[:,1] += T1min_first_hist[b,f,:]
    #     plt.xlabel('angle')
    #     plt.ylabel('First T1 orientation')
    #     plt.title('force ' + fval[f])
    #     plt.legend()

### Paper plot
# At Kees' request: Orientation of first T1s. See how much info can go into a single plot
# Final version: f = 0.125, beta = 1
# NOPE: final version f = 0.15, beta = 0.8
fidx = 6
bidx = 4
plt.figure()
ax = plt.subplot(111, projection='polar')
if sum(isT1[bidx,fidx,:])>0:
    dang=anglebinfirst[1]-anglebinfirst[0]
    anglebin2 = anglebinfirst[1:]-dang/2.0+np.pi
    plt.bar(anglebinfirst[1:]-dang/2.0,T1plus_first_hist[bidx,fidx,:],width=0.03,lw=0,edgecolor='g',color='g',alpha=0.75,label='appear')
    plt.bar(anglebin2,T1plus_first_hist[bidx,fidx,:],width=0.03,lw=0,edgecolor='g',color='g',alpha=0.75)
    plt.bar(anglebin2,T1min_first_hist[bidx,fidx,:],width=0.03,lw=0,edgecolor='r',color='r',alpha=0.75,label='disappear')
    plt.bar(anglebinfirst[1:]-dang/2.0,T1min_first_hist[bidx,fidx,:],width=0.03,lw=0,edgecolor='r',color='r',alpha=0.75)
plt.xlabel('angle')
plt.ylabel('First T1 orientation')
plt.title('force ' + fval[fidx] + ' beta ' + betaval[bidx])
plt.legend()



### In the paper, plotting results for time 400 in the end, i.e. 300 in real time units since starting activity   
# Integral, if we can identify V with dot epsilontot
# do only over plastic time frame
start = 99
dt = 1
convext = np.zeros((len(betaval),len(fval),5,3))
dconvext = np.zeros((len(betaval),len(fval),5,3))
times = [150,200,400,600,800]
dav=5
plotthis=False
for b in range(len(betaval)):
    if plotthis:
        plt.figure()
    color=plt.cm.nipy_spectral(np.linspace(0,1,len(fval)+1))
    for f,c in zip(range(len(fval)),color):
        epsxx0 = np.cumsum(Vval[b,f,:,start:,0,0],axis=1)*dt
        epsyy0 = np.cumsum(Vval[b,f,:,start:,1,1],axis=1)*dt
        epsxy0 = np.cumsum(Vval[b,f,:,start:,0,1],axis=1)*dt
        epsxx = np.sum(epsxx0[:nseeds[b,f],:],axis=0)/(nseeds[b,f]-nbad[b,f])
        depsxx = np.std(epsxx0[:nseeds[b,f],:],axis=0)/np.sqrt(nseeds[b,f]-nbad[b,f])
        epsxy = np.sum(epsxy0[:nseeds[b,f],:],axis=0)/(nseeds[b,f]-nbad[b,f])
        depsxy = np.std(epsxy0[:nseeds[b,f],:],axis=0)/np.sqrt(nseeds[b,f]-nbad[b,f])
        epsyy = np.sum(epsyy0[:nseeds[b,f],:],axis=0)/(nseeds[b,f]-nbad[b,f])
        depsyy = np.std(epsyy0[:nseeds[b,f],:],axis=0)/np.sqrt(nseeds[b,f]-nbad[b,f])
        if plotthis:
            tvals = dt*np.arange(start,Nsnap-1)
            plt.plot(tvals,epsxx,'-',color=c,label=fval[f])
            #plt.plot(tvals,epsxy,':',color=c)
            plt.plot(tvals,epsyy,'--',color=c)

        # Locate point of maximum convergence-extension and average around it
        #maxpt = np.argmax(epsyy-epsxx)
        for t in range(len(times)):
            convext[b,f,t,0] = np.average(epsxx[(times[t]-start-dav):(times[t]-start+dav)])
            dconvext[b,f,t,0] = np.average(depsxx[(times[t]-start-dav):(times[t]-start+dav)])
            convext[b,f,t,1] = np.average(epsxy[(times[t]-start-dav):(times[t]-start+dav)])
            dconvext[b,f,t,1] = np.average(depsxy[(times[t]-start-dav):(times[t]-start+dav)])
            convext[b,f,t,2] = np.average(epsyy[(times[t]-start-dav):(times[t]-start+dav)])
            dconvext[b,f,t,2] = np.average(depsyy[(times[t]-start-dav):(times[t]-start+dav)])


    if plotthis:   
        plt.plot(tvals,0*tvals,'--k',lw=2)
        plt.xlabel('time')
        plt.ylabel('Total strain from V')
        plt.title('Beta =' + betaval[b])
        plt.legend()


# Convergence-extension at time 400
### Paper plot
plt.figure()
color3=plt.cm.nipy_spectral(np.linspace(0,1,len(betaval)+1))
for b,c in zip(range(len(betaval)),color3):
    plt.errorbar(fs,convext[b,:,2,0]-convext[b,:,2,2],dconvext[b,:,2,0]+dconvext[b,:,2,2],marker='o',linestyle='--',color=c,label=betaval[b])
plt.xlabel('pulling force')
plt.ylabel('shear strain') 
plt.title('Convergence-extension xx-yy, time 400') 
plt.legend()


# plt.figure()
# for b,c in zip(range(len(betaval)),color3):
#     plt.plot(fs,convext[b,:,4,0]-convext[b,:,4,2],marker='o',linestyle='--',color=c,label=betaval[b])
# plt.xlabel('pulling force')
# plt.ylabel('shear strain') 
# plt.title('Convergence-extension xx-yy, time 800') 
# plt.legend()

# Compression part 
plt.figure()
for b,c in zip(range(len(betaval)),color3):
    plt.plot(fs,convext[b,:,2,0]+convext[b,:,2,2],dconvext[b,:,2,0]+convext[b,:,2,2],marker='o',linestyle='--',color=c,label=betaval[b])
plt.xlabel('pulling force')
plt.ylabel('compression') 
plt.title('Convergence-extension xx+yy, time 400') 
plt.legend()


# plt.figure()
# for b,c in zip(range(len(betaval)),color3):
#     plt.plot(fs,convext[b,:,4,0]+convext[b,:,4,2],marker='o',linestyle='--',color=c,label=betaval[b])
# plt.xlabel('pulling force')
# plt.ylabel('compression') 
# plt.title('Convergence-extension xx+yy, time 800') 
# plt.legend()


# Finish with some pretty plots of the U, V and P
## Final version: f = 0.125, beta = 1.0
# NOPE: final version f = 0.15, beta = 0.8
fidx = 6
bidx = 4
for bidx in range(bidx,bidx+1):


    # Look for the averages and integrals
    # Mean strain U since activity was turned on 
    # epsxx0 = np.cumsum(Vval[b,f,:,start:,0,0],axis=1)*dt
    start=99
    Uav = np.average(Uval[bidx,fidx,:,start:,:,:],axis=0)
    dUav = np.std(Uval[bidx,fidx,:,start:,:,:],axis=0)
    plt.figure()
    tvals = dt*np.arange(start,Nsnap)
    plt.plot(tvals,Uav[:,0,0],'-r',label='exx')
    plt.fill_between(tvals, Uav[:,0,0]-dUav[:,0,0], Uav[:,0,0]+dUav[:,0,0], color='r', alpha=.2)
    plt.plot(tvals,Uav[:,0,1],'-g',label='exy')
    plt.fill_between(tvals, Uav[:,0,1]-dUav[:,0,1], Uav[:,0,1]+dUav[:,0,1], color='g', alpha=.2)
    plt.plot(tvals,Uav[:,1,1],'-k',label='eyy')
    plt.fill_between(tvals, Uav[:,1,1]-dUav[:,1,1], Uav[:,1,1]+dUav[:,1,1], color='k', alpha=.2)
    plt.xlabel('time')
    plt.ylabel('Total strain from U')
    plt.title('beta ' +betaval[bidx] + ', force ' + fval[fidx])
    plt.legend()


    # Integral, if we can identify V with dot epsilontot
    # do only over plastic time frame
    start = 99
    epsxx0 = np.cumsum(Vval[bidx,fidx,:,start:,0,0],axis=1)*dt
    epsyy0 = np.cumsum(Vval[bidx,fidx,:,start:,1,1],axis=1)*dt
    epsxy0 = np.cumsum(Vval[bidx,fidx,:,start:,0,1],axis=1)*dt
    epsxx = np.average(epsxx0,axis=0)
    depsxx = np.std(epsxx0,axis=0)
    epsxy = np.average(epsxy0,axis=0)
    depsxy = np.std(epsxy0,axis=0)
    epsyy = np.average(epsyy0,axis=0)
    depsyy = np.std(epsyy0,axis=0)
    plt.figure()
    tvals = dt*np.arange(start,Nsnap-1)
    plt.plot(tvals,epsxx,'-r',label='exx')
    plt.fill_between(tvals, epsxx-depsxx, epsxx+depsxx, color='r', alpha=.2)
    plt.plot(tvals,epsxy,'-g',label='exy')
    plt.fill_between(tvals, epsxy-depsxy, epsxy+depsxy, color='g', alpha=.2)
    plt.plot(tvals,epsyy,'-k',label='eyy')
    plt.fill_between(tvals, epsyy-depsyy, epsyy+depsyy, color='k', alpha=.2)
    plt.xlabel('time')
    plt.ylabel('Total strain from V')
    plt.title('beta ' +betaval[bidx] + ', force ' + fval[fidx])
    plt.ylim(-0.1,0.1)
    plt.legend()


    start = 99
    epsxx0 = np.cumsum(Pval[bidx,fidx,:,start:,0,0],axis=1)*dt
    epsyy0 = np.cumsum(Pval[bidx,fidx,:,start:,1,1],axis=1)*dt
    epsxy0 = np.cumsum(Pval[bidx,fidx,:,start:,0,1],axis=1)*dt
    epsxx = np.average(epsxx0,axis=0)
    depsxx = np.std(epsxx0,axis=0)
    epsxy = np.average(epsxy0,axis=0)
    depsxy = np.std(epsxy0,axis=0)
    epsyy = np.average(epsyy0,axis=0)
    depsyy = np.std(epsyy0,axis=0)

    plt.figure()
    tvals = dt*np.arange(start,Nsnap-1)
    plt.plot(tvals,epsxx,'-r',label='exx')
    plt.fill_between(tvals, epsxx-depsxx, epsxx+depsxx, color='r', alpha=.2)
    plt.plot(tvals,epsxy,'-g',label='exy')
    plt.fill_between(tvals, epsxy-depsxy, epsxy+depsxy, color='g', alpha=.2)
    plt.plot(tvals,epsyy,'-k',label='eyy')
    plt.fill_between(tvals, epsyy-depsyy, epsyy+depsyy, color='k', alpha=.2)
    plt.xlabel('time')
    plt.ylabel('Total strain from P')
    plt.title('beta ' +betaval[bidx] + ', force ' + fval[fidx])
    plt.legend()


# And, to finish, some actual myosin and tension distributions on junctions

# Now look at the actual histograms to make sure we understand what's going on.

# measure at 3 points: end of passive phase, couple of steps into active phase, far into simulation
# take average over 5 points either side

# Note: Turned out to be far from clear, not including in the publication. Need to approach issue differently
plotMyoTension = False
if plotMyoTension:
    nbin=40
    tmeasure = [90,130,600]
    myohist_markers = np.zeros((len(betaval),len(fval),3,nbin))
    myohist_markers2 = np.zeros((len(betaval),len(fval),3,3))
    for b in range(len(betaval)):
        for f in range(len(fval)):
            for t in range(len(tmeasure)):
                myohist_markers[b,f,t,:] = np.average(np.average(myohist[b,f,:,(tmeasure[t]-5):(tmeasure[t]+5),:],axis=0),axis=0)


    dm = myobin[1]-myobin[0]
    color=plt.cm.nipy_spectral(np.linspace(0,1,len(fval)+1))
    for b in range(bidx,bidx+1):
        plt.figure()
        for f,c in zip(range(len(fval)),color):
            #plt.plot(myobin[1:]-dm/2.0,myohist_markers[b,f,0,:],'.--',color=c)
            plt.plot(myobin[1:]-dm/2.0,myohist_markers[b,f,1,:],'.-',color=c,label=fval[f])
            #plt.plot(myobin[1:]-dm/2.0,myohist_markers[b,f,1,:],'.-',color=c,label=fval[f])
        plt.axvline(0.5,color='k')
        plt.axvline(1,color='k')
        plt.xlabel('myosin')
        plt.ylabel('probability')
        plt.legend()
        plt.title('Myosin dist, beta=' + betaval[b])
        plt.xlim(0,1.5)
        plt.ylim(0,10)


    dm = myobin[1]-myobin[0]
    color=plt.cm.nipy_spectral(np.linspace(0,1,len(betaval)+1))
    for f in range(fidx,fidx+1):
        plt.figure()
        for b,c in zip(range(len(betaval)),color):
            #plt.plot(myobin[1:]-dm/2.0,myohist_markers[b,f,0,:],'.-',color=c,label=betaval[b])
            plt.plot(myobin[1:]-dm/2.0,myohist_markers[b,f,1,:],'.-',color=c,label=betaval[b])
            #plt.plot(myobin[1:]-dm/2.0,myohist_markers[b,f,1,:],'.-',color=c,label=betaval[b])
        plt.axvline(0.5,color='k')
        #plt.axvline(1,color='k')
        plt.xlabel('myosin')
        plt.ylabel('probability')
        plt.legend()
        plt.title('Myosin dist, f=' + fval[f])
        plt.xlim(0,1.5)
        plt.ylim(0,10)



    nbin=40
    tmeasure = [90,150,600]
    tenshist_markers = np.zeros((len(betaval),len(fval),3,nbin))
    for b in range(len(betaval)):
        for f in range(len(fval)):
            for t in range(len(tmeasure)):
                tenshist_markers[b,f,t,:] = np.average(np.average(tenshist[b,f,:,(tmeasure[t]-5):(tmeasure[t]+5),:],axis=0),axis=0)


    dt = tensbin[1]-tensbin[0]
    color=plt.cm.nipy_spectral(np.linspace(0,1,len(fval)+1))
    for b in range(bidx,bidx+1):
        plt.figure()
        for f,c in zip(range(len(fval)),color):
            #plt.plot(myobin[1:]-dm/2.0,myohist_markers[b,f,0,:],'.--',color=c)
            plt.plot(tensbin[1:]-dt/2.0,tenshist_markers[b,f,1,:],'.-',color=c,label=fval[f])
            #plt.plot(myobin[1:]-dm/2.0,myohist_markers[b,f,1,:],'.-',color=c,label=fval[f])
        plt.axvline(0.3,color='k')
        #plt.axvline(1,color='k')
        plt.xlabel('tension')
        plt.ylabel('probability')
        plt.legend()
        plt.title('Tension dist, beta=' + betaval[b])


    color=plt.cm.nipy_spectral(np.linspace(0,1,len(betaval)+1))
    for f in range(fidx,fidx+1):
        plt.figure()
        for b,c in zip(range(len(betaval)),color):
            #plt.plot(myobin[1:]-dm/2.0,myohist_markers[b,f,0,:],'.-',color=c,label=betaval[b])
            plt.plot(tensbin[1:]-dm/2.0,tenshist_markers[b,f,1,:],'.-',color=c,label=betaval[b])
            #plt.plot(myobin[1:]-dm/2.0,myohist_markers[b,f,1,:],'.-',color=c,label=betaval[b])
        plt.axvline(0.3,color='k')
        plt.xlabel('tension')
        plt.ylabel('probability')
        plt.legend()
        plt.title('Tension dist, f=' + fval[f])



plt.show()

