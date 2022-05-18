import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat

#experiment=['anterior_domain','sickle_domain']
experiment=['Domain5940_142','Domain5940_487','Domain5940_832','Domain5940_1177','Domain5940_1522','Domain5940_1867','sickle_domain','anterior_domain']
micron_start = np.array([142,487,832,1177,1522,1867])
#maxidx = [19,8]
maxidx = [8,8,8,8,8,8,8,19]
ntranche=10
naverage=7

# pixel to micron conversion (needed from Kees)
# Note that since everything is now dimensionless, we just need it for the labels of the plot along the streak
# everything is in pixels and frame units 1 pixel = 0.65 micron, 1 frame 3 minutes
# Find a circular region of x pixels, segment that
# First square, track over time. 2) edges flow in/out 3) inside circular region 150px. Compute all cells. Track those cells.
# Every two frames are anew from the same 150px region, we are not following them over time.

# V P, dUdt are all strain rates. Integrate them backwards to compute actual strains
# These should be dimensionless, and multiplied by the time difference between frames = 3 minutes
conv = 0.65
dt = 3 # minutes, until further notice

# After a lot of debate, we have decided to swap axes on everything: what is x becomes y, what is y becomes x
# and since it's nematic, we don't care about signs

plotraw = False
plotind = False

# We are keeping the following:
Mav_out = np.zeros((len(experiment),2,2))
polangle_out = np.zeros((len(experiment),))
polmag_out = np.zeros((len(experiment),))
flowangle_out = np.zeros((len(experiment),))
flowmag_out = np.zeros((len(experiment),))
VintR_out = np.zeros((len(experiment),ntranche*maxidx[0],2,2))
UrefR_out = np.zeros((len(experiment),ntranche*maxidx[0],2,2))
Vint_noR_out = np.zeros((len(experiment),ntranche*maxidx[0],2,2))
Uref_noR_out = np.zeros((len(experiment),ntranche*maxidx[0],2,2))
# first painstakingly disentangle the matlab structure / dictionary mess
for u in range(len(experiment)):
#for u in range(1,2):
    V = np.zeros((ntranche*maxidx[u],4))
    P = np.zeros((ntranche*maxidx[u],4))
    dUdt = np.zeros((ntranche*maxidx[u],4))
    M = np.zeros((ntranche*maxidx[u],4))
    nMlinks = np.zeros((ntranche*maxidx[u],))
    #T = np.zeros((ntranche*maxidx[u],4))
    for k in range(maxidx[u]):
        # because ... because
        if u>=6:
            filename = experiment[u] + '/tensor_store/tensor_store_' + f"{10*k+1:04}" + '_' + f"{10*(k+1)+1:04}" + ".mat"
        else:
            filename = experiment[u] + '/tensor_store/tensor_store_' + f"{10*(k+1)+2:04}" + '_' + f"{10*(k+2)+2:04}" + ".mat"
        print(filename)
        data = loadmat(filename)['tensor_store']
        V0 = data['V'][0]
        P0 = data['P'][0]
        M0  = data['M'][0]
        dUdt0 = data['dUdt'][0]
        #dMdt0 = data['dM_dt'][0]
        #T0 = data['T'][0]
        for l in range(ntranche):
            V[k*ntranche+l,:] = V0[l][0]
            P[k*ntranche+l,:] = P0[l][0]
            dUdt[k*ntranche+l,:] = dUdt0[l][0][:4]
            M[k*ntranche+l,:] = M0[l][0][:4]/M0[l][0][4] # possibly the right normalisation (number of links)
            nMlinks[k*ntranche+l]=M0[l][0][4]

    time = np.linspace(0,ntranche*maxidx[u]-1,ntranche*maxidx[u])*dt
    Vint = dt*np.cumsum(V,axis=0)
    Pint = dt*np.cumsum(P,axis=0)
    # Whatever this is, it is *not* U as we define it in the simulation. Disregard
    Uint = dt*np.cumsum(dUdt,axis=0) # Highly unclear what this really is
    
    if plotraw:
        plt.figure()
        plt.plot(time,Vint[:,0],'r-',label='exx')
        plt.plot(time,Vint[:,1],'g-',label='exy')
        #plt.plot(time,Vint[:,1],'y-',label='eyx')
        plt.plot(time,Vint[:,3],'b-',label='eyy')
        plt.xlabel('time')
        plt.legend()
        plt.ylim(-0.4,0.25)
        plt.ylabel('Integrated V stress ' + experiment[u])

        # plt.figure()
        # plt.plot(time,Pint[:,0],'r-',label='exx')
        # plt.plot(time,Pint[:,1],'g-',label='exy')
        # #plt.plot(time,Pint[:,1],'y-',label='eyx')
        # plt.plot(time,Pint[:,3],'b-',label='eyy')
        # plt.xlabel('time')
        # plt.legend()
        # plt.ylim(-0.4,0.4)
        # plt.ylabel('Integrated P stress '+ experiment[u])


        # plt.figure()
        # plt.plot(time,M[:,0],'r-',label='exx')
        # plt.plot(time,M[:,1],'g-',label='exy')
        # #plt.plot(time,M[:,1],'y-',label='eyx')
        # plt.plot(time,M[:,3],'b-',label='eyy')
        # plt.xlabel('time')
        # plt.legend()
        # plt.ylabel('Mean fabric tensor '+ experiment[u])

    # Meanwhile, M is the fabric tensor (with unsure normalisation -- probably a number of links in here
    # Cell elongation polarisation = its mean
    Mav = np.reshape(np.mean(M,axis=0),(2,2))
    Mav_out[u,:,:] = Mav

    # Find out the orienatin of the fabric tensor, later reorient with mean fabric tensor eigensystem
    eval, evec0 = np.linalg.eig(Mav)
    # 0 element needs to be the largest one (geometric polarisation direction), swap if necessary
    evec = np.zeros((2,2))
    # if it's the wrong way round (as it will be), swap
    if eval[0]<eval[1]:
        evec[:,0] = evec0[:,1]
        emax = eval[1]
        emin = eval[0]
    # otherwise, don't
    else:
        evec[:,0] = evec0[:,0]
        emax = eval[0]
        emin = eval[1]
    # ... and make sure that we actually do have a rotation matrix (may need to swap sign)
    evec[0,1] = -evec[1,0]
    evec[1,1] = evec[0,0]
    # rotation matrix: (cos theta, - sin theta, sin theta, cos theta)
    # polarisation angle in degrees
    polangle = np.arccos(evec[0,0])/(2*np.pi)*360
    # dimensionless polarisation magnitude
    polmag = (emax-emin)/(emax+emin)
    polangle_out[u] = polangle 
    polmag_out[u] = polmag 

    # eigenvalues on diagonal
    Mrot = np.dot(np.transpose(evec),np.dot(Mav,evec))

    # Variant: Diagonalise using final slice of V instead
    Vfinal = np.reshape(np.mean(Vint[70:80,:],axis=0),(2,2))
    # Reorient with mean fabric tensor eigensystem
    evalV, evec0V = np.linalg.eig(Vfinal)
    print(evalV)
    # 1 element needs to be the largest one (orthogonal to flow), swap if necessary
    # same procedure as with M
    evecV = np.zeros((2,2))
    if evalV[0]>evalV[1]:
        evecV[:,0] = evec0V[:,1]
        emax = evalV[1]
        emin = evalV[0]
    else:
        evecV[:,0] = evec0V[:,0]
        emax = evalV[0]
        emin = evalV[1]
    # ... and make sure that we actually do have a rotation matrix (may need to swap sign)
    evecV[0,1] = -evecV[1,0]
    evecV[1,1] = evecV[0,0]
    # rotation matrix: (cos theta, - sin theta, sin theta, cos theta)
    # flow angle is the orthogonal one in this definition
    flowangle= np.arccos(evecV[1,0])/(2*np.pi)*360
    # note misnomer here: max is smaller than min as we choose the flow to be along y
    # also, the normalised version would measue shear / compression, not what we want
    flowmag = (-emax+emin)# note sign; same as conv-ext measure in simulation
    flowangle_out[u] = flowangle
    flowmag_out[u] = flowmag


    # Calculate U tensor the same way as in the simulation (Graneraverages.py)
    # Definition of U = 1/2 ( log M - log M0), where log M = R^(-1) diag (log lambda) R
    # Here we define M0 to be the spherically symmetrised average M
    M0ref = np.diag([0.5*(eval[0]+eval[1]),0.5*(eval[0]+eval[1])])
    # Redundant in this case since diagonal
    M0vals, M0vecs = np.linalg.eigh(M0ref)
    logM0 = np.dot(np.transpose(M0vecs),np.dot(np.diag(np.log(M0vals)),M0vecs))

    # And now the full matrices
    Uref = np.zeros((ntranche*maxidx[u],2,2))
    for k in range(ntranche*maxidx[u]):
        Mvals, Mvecs = np.linalg.eigh(M[k,:].reshape(2,2))
        logM = np.dot(np.transpose(Mvecs),np.dot(np.diag(np.log(Mvals)),Mvecs))
        Uref[k,:,:] = 0.5*(logM - logM0)

    if plotraw:

        plt.figure()
        plt.plot(time,Uref[:,0,0],'r-',label='exx')
        plt.plot(time,Uref[:,0,1],'g-',label='exy')
        #plt.plot(time,Uint[:,1],'y-',label='eyx')
        plt.plot(time,Uref[:,1,1],'b-',label='eyy')
        plt.xlabel('time')
        plt.legend()
        plt.ylim(-0.35,0.3)
        plt.ylabel('Integrated U stress '+ experiment[u])

    # Slow and steady through the gemetric thicket ...
    VintR = np.zeros((ntranche*maxidx[u],2,2))
    PintR = np.zeros((ntranche*maxidx[u],2,2))
    UintR = np.zeros((ntranche*maxidx[u],2,2))
    MR = np.zeros((ntranche*maxidx[u],2,2))
    UrefR = np.zeros((ntranche*maxidx[u],2,2))
    rotate='flow' # flow is the alternative
    if rotate=='pol':
        for k in range(ntranche*maxidx[u]):
            VintR[k,:,:] = np.dot(np.transpose(evec),np.dot(Vint[k,:].reshape(2,2),evec))
            PintR[k,:,:] = np.dot(np.transpose(evec),np.dot(Pint[k,:].reshape(2,2),evec))
            UintR[k,:,:] = np.dot(np.transpose(evec),np.dot(Uint[k,:].reshape(2,2),evec))
            MR[k,:,:] = np.dot(np.transpose(evec),np.dot(M[k,:].reshape(2,2),evec))
            UrefR[k,:,:] = np.dot(np.transpose(evec),np.dot(Uref[k,:].reshape(2,2),evec))
    else:
        for k in range(ntranche*maxidx[u]):
            VintR[k,:,:] = np.dot(np.transpose(evecV),np.dot(Vint[k,:].reshape(2,2),evecV))
            PintR[k,:,:] = np.dot(np.transpose(evecV),np.dot(Pint[k,:].reshape(2,2),evecV))
            UintR[k,:,:] = np.dot(np.transpose(evecV),np.dot(Uint[k,:].reshape(2,2),evecV))
            MR[k,:,:] = np.dot(np.transpose(evecV),np.dot(M[k,:].reshape(2,2),evecV))
            UrefR[k,:,:] = np.dot(np.transpose(evecV),np.dot(Uref[k,:].reshape(2,2),evecV))

    VintR_out[u,:,:,:] = VintR[:ntranche*maxidx[0],:,:] # don't care about rotated anterior domain anyway
    UrefR_out[u,:,:,:] = UrefR[:ntranche*maxidx[0],:,:]
    # and here is our funny 'do nothing' guy
    Vint_noR_out[u,:,:,:] = Vint[:ntranche*maxidx[0],:].reshape(ntranche*maxidx[0],2,2)
    Uref_noR_out[u,:,:,:] = Uref[:ntranche*maxidx[0],:].reshape(ntranche*maxidx[0],2,2)
    if plotind:
        plt.figure()
        plt.plot(time,VintR[:,0,0],'r-',label='exx')
        plt.plot(time,VintR[:,0,1],'g-',label='exy')
        plt.plot(time,VintR[:,1,1],'b-',label='eyy')
        plt.xlabel('time')
        plt.legend()
        plt.ylabel('Integrated V stress ' + experiment[u])

        plt.figure()
        plt.plot(time,PintR[:,0,0],'r-',label='exx')
        plt.plot(time,PintR[:,0,1],'g-',label='exy')
        plt.plot(time,PintR[:,1,1],'b-',label='eyy')
        plt.xlabel('time')
        plt.legend()
        plt.ylabel('Integrated P stress '+ experiment[u])

        plt.figure()
        plt.plot(time,UintR[:,0,0],'r-',label='exx')
        plt.plot(time,UintR[:,0,1],'g-',label='exy')
        plt.plot(time,UintR[:,1,1],'b-',label='eyy')
        plt.xlabel('time')
        plt.legend()
        plt.ylabel('Integrated U stress '+ experiment[u])

        plt.figure()
        plt.plot(time,UrefR[:,0,0],'r-',label='exx')
        plt.plot(time,UrefR[:,0,1],'g-',label='exy')
        plt.plot(time,UrefR[:,1,1],'b-',label='eyy')
        plt.xlabel('time')
        plt.legend()
        plt.ylabel('U with respect to isotropic '+ experiment[u])

        plt.figure()
        plt.plot(time,MR[:,0,0],'r-',label='exx')
        plt.plot(time,MR[:,0,1],'g-',label='exy')
        plt.plot(time,MR[:,1,1],'b-',label='eyy')
        plt.xlabel('time')
        plt.legend()
        plt.ylabel('Mean fabric tensor '+ experiment[u])
        plt.title('Angle: ' + str(round(polangle))+ ' deg, polaris]ation: ' + str(round(polmag,3)))

# Final statistics of the whole thing
# focus on sickle simulations now only
time = np.linspace(0,ntranche*maxidx[0]-1,ntranche*maxidx[0])*dt
# plt.figure()
# cols = ['y','r','m','b','c','g','k','k']
# for u in range(len(experiment)):
#     if u==0:
#         plt.plot(time,VintR_out[u,:,0,0],cols[u],label=experiment[u])
#         #plt.plot(time,VintR_out[u,:,0,1],cols[u],label='exy')
#         #plt.plot(time,VintR_out[u,:,1,1],cols[u],label='eyy')
#     else:
#         plt.plot(time,VintR_out[u,:,0,0],cols[u],label=experiment[u])
#         #plt.plot(time,VintR_out[u,:,0,1],cols[u])
#         #plt.plot(time,VintR_out[u,:,1,1],cols[u])
# plt.xlabel('time')
# plt.legend()
# plt.ylabel('Integrated V stress')


Vintav = np.average(VintR_out[:naverage,:,:,:],axis=0)
dVintav = np.std(VintR_out[:naverage,:,:,:],axis=0)
plt.figure()
plt.plot(time,Vintav[:,0,0],'r-',label='exx')
plt.fill_between(time, Vintav[:,0,0]-dVintav[:,0,0], Vintav[:,0,0]+dVintav[:,0,0], color='r', alpha=.2)
plt.plot(time,Vintav[:,0,1],'g-',label='exy')
plt.fill_between(time, Vintav[:,0,1]-dVintav[:,0,1], Vintav[:,0,1]+dVintav[:,0,1], color='g', alpha=.2)
plt.plot(time,Vintav[:,1,1],'b-',label='eyy')
plt.fill_between(time, Vintav[:,1,1]-dVintav[:,1,1], Vintav[:,1,1]+dVintav[:,1,1], color='b', alpha=.2)
plt.xlabel('time')
plt.legend()
plt.ylabel('Integrated V stress')  

# plt.figure()
# for u in range(len(experiment)):
#     if u==0:
#         plt.plot(time,UrefR_out[u,:,0,0],cols[u],label=experiment[u])
#         #plt.plot(time,UrefR_out[u,:,0,1],cols[u],label='exy')
#         #plt.plot(time,UrefR_out[u,:,1,1],cols[u],label='eyy')
#     else:
#         plt.plot(time,UrefR_out[u,:,0,0],cols[u],label=experiment[u])
#         #plt.plot(time,UrefR_out[u,:,0,1],cols[u])
#         #plt.plot(time,UrefR_out[u,:,1,1],cols[u])
# plt.xlabel('time')
# plt.legend()
# plt.ylabel('U with respect to isotropic')

Urefav = np.average(UrefR_out[:naverage,:,:,:],axis=0)
dUrefav = np.std(UrefR_out[:naverage,:,:,:],axis=0)
plt.figure()
plt.plot(time,Urefav[:,0,0],'r-',label='exx')
plt.fill_between(time, Urefav[:,0,0]-dUrefav[:,0,0], Urefav[:,0,0]+dUrefav[:,0,0], color='r', alpha=.2)
plt.plot(time,Urefav[:,0,1],'g-',label='exy')
plt.fill_between(time, Urefav[:,0,1]-dUrefav[:,0,1], Urefav[:,0,1]+dUrefav[:,0,1], color='g', alpha=.2)
plt.plot(time,Urefav[:,1,1],'b-',label='eyy')
plt.fill_between(time, Urefav[:,1,1]-dUrefav[:,1,1], Urefav[:,1,1]+dUrefav[:,1,1], color='b', alpha=.2)
plt.xlabel('time')
plt.legend()
plt.ylabel('U with respect to isotropic')


# Try some spatial variation. Second bit of name corresponds to micron distance from ... not sure
# take these first 6 which are sequential downward and see if anything interesting is going on
# Swap x and y here with the 90 degree thing
fig, ax1 = plt.subplots()
npts= len(experiment)-2
ax1.plot(conv*micron_start[:npts],polangle_out[:npts],'o-r',label='pol angle')
ax1.plot(conv*micron_start[:npts],flowangle_out[:npts],'x-r',label='flow angle')
ax1.set_xlabel('relative streak y (micron)')
ax1.set_ylabel('angle',color='tab:red')
#ax1.set_ylim(0,130)
ax1.legend()
ax2=ax1.twinx()
ax2.plot(conv*micron_start[:npts],polmag_out[:npts],'o-k',label='pol mag')
ax2.plot(conv*micron_start[:npts],flowmag_out[:npts],'x-k',label='flow mag')
ax2.set_ylabel('magnitude')
ax2.set_ylim(0,1)
ax2.legend()

# Different layout
fig, axs = plt.subplots(2, sharex=True)
#fig = plt.figure()
#gs = fig.add_gridspec(2,1, hspace=0)
#axs = gs.subplots(sharex=True)
fig.suptitle('Sickle polarisation')
plt.subplots_adjust(hspace = 0.0)
# note rotation
axs[0].plot(conv*micron_start[:npts],polangle_out[:npts]-90,'o-g',lw=2,label='pol angle')
axs[0].plot(conv*micron_start[:npts],-(flowangle_out[:npts]-90),'x-r',lw=2,label='flow angle')
axs[0].plot(conv*micron_start[:npts],0.0*micron_start[:npts],'--k',lw=2)
axs[0].plot(conv*micron_start[:npts],90.0*micron_start[:npts]/micron_start[:npts],'--k',lw=2)
axs[0].set_ylabel('angle')
axs[0].set_ylim(-35,100)
axs[0].set_xlim(conv*micron_start[0]-20,conv*micron_start[npts-1]+20)
axs[0].legend()
axs[1].plot(conv*micron_start[:npts],polmag_out[:npts],'o-g',lw=2,label='pol mag')
axs[1].plot(conv*micron_start[:npts],flowmag_out[:npts],'x-r',lw=2,label='flow mag')
axs[1].set_ylabel('magnitude')
axs[1].set_xlabel('streak x (micron)')
axs[1].legend()
axs[1].set_ylim(0,0.75)
axs[1].set_xlim(conv*micron_start[0]-20,conv*micron_start[npts-1]+20)
# Hide x labels and tick labels for all but bottom plot.
for ax in axs:
    ax.label_outer()

# And what do we decide to focus on? Just the two middle simulations with no rotation whatsoever ...
# except in the labels for plotting :)
euse = [2,3]
mrkr = ['o','x']
plt.figure()
v=0
for u in euse:
    plt.plot(time,Vint_noR_out[u,:,0,0],'b-',marker=mrkr[v],label='eyy',markevery=7)
    plt.plot(time,Vint_noR_out[u,:,0,1],'g-',marker=mrkr[v],label='exy',markevery=7)
    plt.plot(time,Vint_noR_out[u,:,1,1],'r-',marker=mrkr[v],label='exx',markevery=7)
    v+=1
plt.xlabel('time')
plt.legend()
plt.xlim(-5,242)
plt.ylim(-0.45,0.25)
plt.ylabel('Integrated V stress')

plt.figure()
v=0
for u in euse:
    plt.plot(time,Uref_noR_out[u,:,0,0],'b-',marker=mrkr[v],label='eyy',markevery=7)
    plt.plot(time,Uref_noR_out[u,:,0,1],'g-',marker=mrkr[v],label='exy',markevery=7)
    plt.plot(time,Uref_noR_out[u,:,1,1],'r-',marker=mrkr[v],label='exx',markevery=7)
    v+=1
plt.xlabel('time')
plt.xlim(-5,242)
plt.legend()
plt.ylim(-0.32,0.25)
plt.ylabel('U stress')


plt.show()