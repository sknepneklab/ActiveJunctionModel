import numpy as np  
import pickle
import matplotlib.pyplot as plt


# Seeds for our simulations (maximum 32)
seeds = [1,2,3,4,5]
#seeds = [5]
Nactive = 520 # 520 cells in our lone disordered configuration, currently

# Values of active myosin coupling connstant
betaval=['0.5']
# Values of external pulling force
fval=['0.0','0.05','0.1','0.15','0.2','0.3','0.4']

# Clunky but correct: actual number of seeds run. Note that for the nseeds = 16, no actual T1s happen at all
# and simulation is very nearly deterministic. N = 16 is plenty for averaging here.
nseeds = 5*np.ones((len(betaval),len(fval))).astype(int)
#nseeds = 1*np.ones((len(betaval),len(fval))).astype(int)

# We now have a handful of crashed simulations, all at beta 1.4 
# label: bet idx, f idx, seed idx
badpts=[[]]
nbad=np.zeros((len(betaval),len(fval)))

# All simulations are 800 time units long, with a spacing of 2 time unit between snapshots
# after the first 50 snapshots (equilibration). We discard those.
Nsnap = 849
#Tmax = 300 # up to where to analyse T1s (experimental)
Tmax = Nsnap
# crashes, correct normalisation
NsnapReal = np.zeros((len(betaval),len(fval),len(seeds))).astype(int)
dt = 2

# Graner coarse-grained U, P and V tensors for the whole region
Uval = np.zeros((len(betaval),len(fval),len(seeds),Nsnap,2,2))
Pval = np.zeros((len(betaval),len(fval),len(seeds),Nsnap-1,2,2))
Vval = np.zeros((len(betaval),len(fval),len(seeds),Nsnap-1,2,2))

# Mean tension and myosin for the whole region
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
nbin =21
anglebin = np.linspace(0,np.pi,nbin+1)

# Histogram of distributions of T1s, separated into 'all' and 'first'
T1plus_all_hist = np.zeros((len(betaval),len(fval),nbin))
T1min_all_hist = np.zeros((len(betaval),len(fval),nbin))


# first T1, existence
isT1 = np.zeros((len(betaval),len(fval),len(seeds)))
# start of active phase in frames
activestart = 50
firstT1time = -1*np.ones((len(betaval),len(fval),len(seeds)))


# Read in all of our pickle files and fill in statistics above
for b in range(len(betaval)):
    for f in range(len(fval)):
        # collecting Graner T1 statistics from time series. 'f' is for first, 'a' is for all, 'p' is for appearing  junction, 'm' is for diseappearing junction
        T1palist = []
        T1malist = []
        T1pflist = []
        T1mflist = []
        for s in range(nseeds[b,f]):
            if not [b,f,seeds[s]] in badpts:

                ######## Graner statistics averaged data #####
                # Graner statistics output pickle files
                picklefile = './picklefiles/beta' + betaval[b] + '/fpull' + fval[f]+'/tensor_timeseries_' + str(seeds[s]) + '_Utest.p'
                print(picklefile)
                data  = pickle.load(open(picklefile,'rb'))
                tval = data['tval']
                
                # Averaged tensors
                # data missing for crashed simulations (and that region is unstable)
                ngood = len(data["U"])
                Uval[b,f,s,:ngood,:,:] = data["U"]
                Pval[b,f,s,:ngood-1,:,:] = data["P"]
                Vval[b,f,s,:ngood-1,:,:] = data["V"]
                NsnapReal[b,f,s]=ngood
                # # Polarisation

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
                        print("half junction!")
                # put all of the angles into the running list of T1s
                # alternative: stop at some time
                currT1s = np.cumsum(nt1)
                if Tmax < ngood:
                    idmax = int(currT1s[Tmax])
                else:
                    idmax = ngood
                T1palist.extend(t1pos[:idmax])
                T1malist.extend(t1neg[:idmax])
                print(len(T1palist))
                print(len(t1pos))
                print(len(t1neg))

                # We are collecting all the tissue statistics here
                Tensavg[b,f,s,:ngood,:,:] = data["Tensavg"]
                Myoavg[b,f,s,:ngood,:,:] = data["Myoavg"]

                # Global statistics
                avArea[b,f,s,:ngood] = data['avArea']
                avPerim[b,f,s,:ngood] = data['avPerim']
                avp0[b,f,s,:ngood] = data['avp0']

                # always the same, simply replace
                p0bin = data['p0bin']
                p0hist[b,f,s,:ngood,:]=data['p0hist']

                myobin = data['myobin']
                myohist[b,f,s,:ngood,:]=data['myohist']

                tensbin = data['tensbin']
                tenshist[b,f,s,:ngood,:]=data['tenshist']

        # wrap T1 angles into upper two quadrants
        # [f(x) if condition else g(x) for x in sequence]
        T1pflist2 = [x if x>0 else (x+np.pi)  for x in T1pflist]
        T1mflist2 = [x if x>0 else (x+np.pi)  for x in T1mflist]
        #T1mflist2_VAR = [x if x>0 else (x+np.pi)  for x in T1mflist]
        T1palist2 = [x if x>0 else (x+np.pi)  for x in T1palist]
        T1malist2 = [x if x>0 else (x+np.pi)  for x in T1malist]
        #print(T1mflist)
        #print(T1mflist2)
        if b==5 and f==4:
            plt.figure()
            plt.plot(T1palist2,'.b')
            plt.plot(T1malist2,'.r')
    
        # Actual histograms. 

        T1plus_all_hist[b,f,:],edges = np.histogram(T1palist2,bins=anglebin)
        T1min_all_hist[b,f,:],edges = np.histogram(T1malist2,bins=anglebin)

        # Normalisation: Divide by number of configurations, except for the 'all' counts where every T1 is in the list twice
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


# Graner statistics T1 orientation plots. Does include non-central T1s for the first T1 transition
# Switch to turn off as it's a lot of plots
plotT1orient=False
if plotT1orient:
    ## Supplementary?
    ## For all T1s (these are reasonably nice)

    # Line plots
    for f in range(len(fval)):
        plt.figure()
        color2=plt.cm.nipy_spectral(np.linspace(0,1,len(betaval)+1))
        for b,c in zip(range(len(betaval)),color2):
            if sum(isT1[b,f,:])>0:
                dang=anglebin[1]-anglebin[0]
                plt.plot(anglebin[1:]-dang/2.0,T1plus_all_hist[b,f,:],marker='o',linestyle='-',color=c,label=betaval[b])
                plt.plot(anglebin[1:]-dang/2.0,-T1min_all_hist[b,f,:],marker='o',linestyle='-',color=c)
        plt.xlabel('angle')
        plt.ylabel('All T1 orientation')
        plt.title('force ' + fval[f])
        plt.legend()


### Paper plot
# At Kees' request: Orientation of all T1s
# Final version: f = 0.2, beta = 0.5
fidx = 4
bidx = 0
#fidx = 4
#bidx = 6
plt.figure()
ax = plt.subplot(111, projection='polar')
if sum(isT1[bidx,fidx,:])>0:
    dang=anglebin[1]-anglebin[0]
    anglebin2 = anglebin[1:]-dang/2.0+np.pi
    plt.bar(anglebin[1:]-dang/2.0,T1plus_all_hist[bidx,fidx,:],width=0.15,lw=0,edgecolor='b',color='b',alpha=0.5,label='appear')
    plt.bar(anglebin2,T1plus_all_hist[bidx,fidx,:],width=0.15,lw=0,edgecolor='b',color='b',alpha=0.5)
    plt.bar(anglebin2,T1min_all_hist[bidx,fidx,:],width=0.15,lw=0,edgecolor='r',color='r',alpha=0.5,label='disappear')
    plt.bar(anglebin[1:]-dang/2.0,T1min_all_hist[bidx,fidx,:],width=0.15,lw=0,edgecolor='r',color='r',alpha=0.5)
    plt.plot(anglebin[1:]-dang/2.0,T1plus_all_hist[bidx,fidx,:],color='b',marker='o',linestyle='none')
    plt.plot(anglebin2,T1plus_all_hist[bidx,fidx,:],color='b',marker='o',linestyle='none')
    plt.plot(anglebin2,T1min_all_hist[bidx,fidx,:],color='r',marker='*',linestyle='none')
    plt.plot(anglebin[1:]-dang/2.0,T1min_all_hist[bidx,fidx,:],marker='*',linestyle='none',color='r')
    
plt.xlabel('angle')
plt.ylabel('All T1 orientation')
plt.title('force ' + fval[fidx] + ' beta ' + betaval[bidx])
plt.legend()



### In the paper, plotting results for time 400 in the end, i.e. 300 in real time units since starting activity   
# Integral, if we can identify V with dot epsilontot
# do only over plastic time frame
start = 49
dt = 2
convext = np.zeros((len(betaval),len(fval),5,3))
dconvext = np.zeros((len(betaval),len(fval),5,3))
times = [150,200,400,600,800]
dav=5
plotthis=True
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
plt.title('Convergence-extension xx-yy, time 700') 
plt.legend()



# Compression part 
plt.figure()
for b,c in zip(range(len(betaval)),color3):
    plt.errorbar(fs,convext[b,:,2,0]+convext[b,:,2,2],dconvext[b,:,2,0]+dconvext[b,:,2,2],marker='o',linestyle='--',color=c,label=betaval[b])
plt.xlabel('pulling force')
plt.ylabel('compression') 
plt.title('Convergence-extension xx+yy, time 700') 
plt.legend()



# Finish with some pretty plots of the U, V and P
## Final version: f = 0.125, beta = 1.0
# NOPE: final version f = 0.15, beta = 0.8
fidx = 4
bidx = 0
#fidx = 4
#bidx = 6
for bidx in range(bidx,bidx+1):


    # Look for the averages and integrals
    # Mean strain U since activity was turned on 
    # epsxx0 = np.cumsum(Vval[b,f,:,start:,0,0],axis=1)*dt
    start=49
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
    start = 49
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
    #plt.ylim(-0.1,0.1)
    plt.legend()


    start = 49
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


# Polarisation tensors
plotPolarisation = True
dt = 2
tenspol = np.zeros((len(betaval),len(fval),5,3))
myopol = np.zeros((len(betaval),len(fval),5,3))
shapepol = np.zeros((len(betaval),len(fval),5,3))
dtenspol = np.zeros((len(betaval),len(fval),5,3))
dmyopol = np.zeros((len(betaval),len(fval),5,3))
dshapepol = np.zeros((len(betaval),len(fval),5,3))
times = [150,200,400,600,800]
dav=5
start = 50
for b in range(len(betaval)):
	if plotPolarisation:
		plt.figure()
	color=plt.cm.nipy_spectral(np.linspace(0,1,len(fval)+1))
	for f,c in zip(range(len(fval)),color):
		Tav = np.average(Tensavg[b,f,:,:,:,:],axis=0)
		dTav = np.std(Tensavg[b,f,:,:,:,:],axis=0)
		for t in range(len(times)):
			tenspol[b,f,t,0] = np.average(Tav[(times[t]-start-dav):(times[t]-start+dav),0,0])
			tenspol[b,f,t,1] = np.average(Tav[(times[t]-start-dav):(times[t]-start+dav),0,1])
			tenspol[b,f,t,2] = np.average(Tav[(times[t]-start-dav):(times[t]-start+dav),1,1])
			dtenspol[b,f,t,0] = np.average(dTav[(times[t]-start-dav):(times[t]-start+dav),0,0])
			dtenspol[b,f,t,1] = np.average(dTav[(times[t]-start-dav):(times[t]-start+dav),0,1])
			dtenspol[b,f,t,2] = np.average(dTav[(times[t]-start-dav):(times[t]-start+dav),1,1])
		if plotPolarisation:
			plt.plot(tvals,Tav[start:,0,0],'-',color=c,label=fval[f])
			plt.plot(tvals,Tav[start:,0,1],':',color=c)
			plt.plot(tvals,Tav[start:,1,1],'--',color=c)
	if plotPolarisation:
		plt.plot(tvals,0*tvals,'--k',lw=2)
		plt.xlabel('time')
		plt.ylabel('Tension tensor')
		plt.title('Beta =' + betaval[b])
		plt.legend()
	
for b in range(len(betaval)):
	if plotPolarisation:
		plt.figure()
	color=plt.cm.nipy_spectral(np.linspace(0,1,len(fval)+1))
	for f,c in zip(range(len(fval)),color):
		Mav = np.average(Myoavg[b,f,:,:,:,:],axis=0)
		dMav = np.std(Myoavg[b,f,:,:,:,:],axis=0)
		for t in range(len(times)):
			myopol[b,f,t,0] = np.average(Mav[(times[t]-start-dav):(times[t]-start+dav),0,0])
			myopol[b,f,t,1] = np.average(Mav[(times[t]-start-dav):(times[t]-start+dav),0,1])
			myopol[b,f,t,2] = np.average(Mav[(times[t]-start-dav):(times[t]-start+dav),1,1])
			dmyopol[b,f,t,0] = np.average(dMav[(times[t]-start-dav):(times[t]-start+dav),0,0])
			dmyopol[b,f,t,1] = np.average(dMav[(times[t]-start-dav):(times[t]-start+dav),0,1])
			dmyopol[b,f,t,2] = np.average(dMav[(times[t]-start-dav):(times[t]-start+dav),1,1])
		if plotPolarisation:
			plt.plot(tvals,Mav[start:,0,0],'-',color=c,label=fval[f])
			plt.plot(tvals,Mav[start:,0,1],':',color=c)
			plt.plot(tvals,Mav[start:,1,1],'--',color=c)
	if plotPolarisation:
		plt.plot(tvals,0*tvals,'--k',lw=2)
		plt.xlabel('time')
		plt.ylabel('Myosin tensor')
		plt.title('Beta =' + betaval[b])
		plt.legend()
		
for b in range(len(betaval)):
	if plotPolarisation:
		plt.figure()
	color=plt.cm.nipy_spectral(np.linspace(0,1,len(fval)+1))
	for f,c in zip(range(len(fval)),color):
		#Upol = np.average(Uav[b,f,:,:,:,:],axis=0)
		#Uav = np.average(Uval[bidx,fidx,:,start:,:,:],axis=0)
		Upol = np.average(Uval[b,f,:,:,:,:],axis=0)
		dUpol = np.std(Uval[b,f,:,:,:,:],axis=0)
		for t in range(len(times)):
			shapepol[b,f,t,0] = np.average(Upol[(times[t]-start-dav):(times[t]-start+dav),0,0])
			shapepol[b,f,t,1] = np.average(Upol[(times[t]-start-dav):(times[t]-start+dav),0,1])
			shapepol[b,f,t,2] = np.average(Upol[(times[t]-start-dav):(times[t]-start+dav),1,1])
			dshapepol[b,f,t,0] = np.average(dUpol[(times[t]-start-dav):(times[t]-start+dav),0,0])
			dshapepol[b,f,t,1] = np.average(dUpol[(times[t]-start-dav):(times[t]-start+dav),0,1])
			dshapepol[b,f,t,2] = np.average(dUpol[(times[t]-start-dav):(times[t]-start+dav),1,1])
		if plotPolarisation:
			plt.plot(tvals,Upol[start:,0,0],'-',color=c,label=fval[f])
			plt.plot(tvals,Upol[start:,0,1],':',color=c)
			plt.plot(tvals,Upol[start:,1,1],'--',color=c)
	if plotPolarisation:
		plt.plot(tvals,0*tvals,'--k',lw=2)
		plt.xlabel('time')
		plt.ylabel('Shape tensor')
		plt.title('Beta =' + betaval[b])
		plt.legend()
		

# A single summary figure at beta = 0.5 ...
bidx = 0
plt.figure()
plt.errorbar(fs,myopol[bidx,:,2,0],dmyopol[bidx,:,2,0],marker='o',linestyle='-',color='g',label='myosin')
plt.errorbar(fs,myopol[bidx,:,2,2],dmyopol[bidx,:,2,2],marker='x',linestyle='--',color='g')
plt.errorbar(fs,tenspol[bidx,:,2,0],dtenspol[bidx,:,2,0],marker='o',linestyle='-',color='r',label='stress')
plt.errorbar(fs,tenspol[bidx,:,2,2],dtenspol[bidx,:,2,2],marker='x',linestyle='--',color='r')
plt.errorbar(fs,shapepol[bidx,:,2,0],dshapepol[bidx,:,2,0],marker='o',linestyle='-',color='k',label='shape')
plt.errorbar(fs,shapepol[bidx,:,2,2],dshapepol[bidx,:,2,2],marker='x',linestyle='--',color='k')
plt.xlabel('pulling force')
plt.ylabel('polarisation (xx,yy)')
plt.title('Polarisations at time 700, beta = ' + betaval[bidx])
plt.legend()






plt.show()

