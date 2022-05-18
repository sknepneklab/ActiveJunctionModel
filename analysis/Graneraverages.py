import pickle
import glob
from AJMAnalysis import *
import numpy as np
import matplotlib.pyplot as plt
import sys as system

class Graneraverages:
	
	def __init__(self, dt,filemax = 'all', basename = 'mesh_', directory = './', start = 0, step = 1, usetypes = ["active","buffer"]):
		
		# locate data and collect set of file names
		self.directory = directory
		self.allfiles = sorted(glob.glob(self.directory + basename+"*.json.bz2")) # read compressed by default
		if len(self.allfiles) == 0:
			self.allfiles = sorted(glob.glob(self.directory + basename+"*.json"))
		if len(self.allfiles) == 0:
			raise Exception("Graneraverages: No files.")
		self.start = start
		self.step = step
		self.filemax = filemax
		if (len(self.allfiles) < self.filemax) or (self.filemax == 'all'):
			self.filemax = len(allfiles)
			print("Graneraverages: Warning: reducing maximum file label to {:d}.".format(self.filemax))
		self.Nsnap = int((self.filemax-self.start)/self.step)
		print("Working on {:d} files.".format(self.Nsnap))
		# Figure out actual time difference between frames
		self.Ntot = np.zeros((self.Nsnap,))
		
			
		# Which data we want
		self.usetypes = usetypes
		
		# Check out some things: First of all, the labels of our cells won't change (?)
		tensor = MTensor(self.allfiles[self.start])
		# Time step
		currtime = tensor.mesh.system['time']
		print(currtime)
		try:
			currtime = tensor.mesh.system['time']
			self.tval = np.zeros((self.Nsnap,))
			print('Getting time from the json files')
		except: 
			self.tval = self.start*dt + dt*np.arange(self.Nsnap) # Get it from the one provided in the script
			print('Getting time manually')
			
		self.labels = []
		huh = []
		i = 0
		for fc in tensor.mesh.faces:
			if self.usetypes == 'all':
				if not fc.outer:
					self.labels.append(i)
					huh.append(fc.idx)
			else:
				if fc.type in self.usetypes:
					self.labels.append(i)
					huh.append(fc.idx)
			i+=1
		self.N = len(self.labels)
		print("Found " + str(self.N) + " cells of types " + str(self.usetypes))
		print("Their labels are " + str(self.labels))
		print("And their internal indices are " + str(huh))

		
		# This tells us what size data we should expect
		self.Neigh = np.zeros((self.Nsnap,self.N))
		
		
	# Run this first, mandatory
	def statMtensor(self,plotvtk=False,verbose=True):
		
		print("Processing M tensor")
		self.Mavg = np.empty((self.Nsnap,2,2))
		u=0
		for i in range(self.start,self.filemax,self.step):
			if i % 50 == 0:
				print("file number {:d}.".format(i))
			txtr = MTensor(self.allfiles[i])
			try: 
				currtime = txtr.mesh.system['time']
				self.tval[i]=currtime
			except:
				pass
			# Global: Compute number of neighbours of the cells we want for every frame
			# That info is only (easily) accessible in texture, so get it out of here
			self.Neigh[u,:] = txtr.Nneigh[self.labels]
			#'np.einsum('j,ik->ijk',np.ones((self.Ntrack,)),hmm)
			linktotal = np.sum(self.Neigh[u,:])
			self.Ntot[u] = linktotal/2.0
			mavg = np.einsum('ijk,i->jk',txtr.M.T[self.labels,:,:],self.Neigh[u,:])/(linktotal)
			self.Mavg[u,:,:] = mavg

			if plotvtk:
				txtrvtk = "mtensor_{:08d}.vtp".format(i)
				txtr.plot_vtk_tensor(txtrvtk)
			u+=1

		#self.Mvals, self.Mvecs = np.linalg.eigh(self.Mavg)
		# ??? decrypt later
		#strain =  Mvecs @ np.array([np.diag(np.log(v)) for v in Mvals]) @ np.transpose(Mvecs, (0,2,1))
		#strain -= strain[n0] # select the reference configuration to be the first one
		
	def statMyotensor(self, plotvtk = False, verbose = True):
		
		print("Processing Myosin tensor")
		self.Myoavg = np.empty((self.Nsnap,2,2))
		u=0
		for i in range(self.start,self.filemax,self.step):
			if i % 50 == 0:
				print("file number {:d}.".format(i))
			myosintr = MyoTensor(self.allfiles[i])
			linktotal = np.sum(self.Neigh[u,:])
			myoavg = np.einsum('ijk,i->jk',myosintr.Myo.T[self.labels,:,:],self.Neigh[u,:])/(linktotal)
			self.Myoavg[u,:,:] = myoavg

			if plotvtk:
				myovtk = "myosin_tensor_{:08d}.vtp".format(i)
				myosintr.plot_vtk_tensor(myovtk)
			u+=1

		#self.Myovals, self.Myovecs = np.linalg.eigh(self.Myoavg)
		
	def statTensiontensor(self, plotvtk = False, verbose = True):
		
		print("Processing Tension tensor")
		self.Tensavg = np.empty((self.Nsnap,2,2))
		u=0
		for i in range(self.start, self.filemax, self.step):
			if i % 50 == 0:
				print("file number {:d}".format(i))
			tensiontr = TensionTensor(self.allfiles[i])
			linktotal = np.sum(self.Neigh[u,:])
			tensavg = np.einsum('ijk,i->jk',tensiontr.Tension.T[self.labels,:,:],self.Neigh[u,:])/(linktotal)
			self.Tensavg[u,:,:] = tensavg

			if plotvtk:
				tenvtk = "tension_tensor_{:08d}.vtp".format(i)
				tensiontr.plot_vtk_tensor(tenvtk)
			u+=1
		#self.Tenvsals, self.Tensvecs = np.linalg.eigh(self.Tensavg)
		
	#Averaging this one: see notes
	def statBCtensor(self, plotvtk = False, verbose = True):
		print("Processing B tensor")
		self.Bavg = np.empty((self.Nsnap-1,2,2))
		self.Cavg = np.empty((self.Nsnap-1,2,2))
		u=0
		for i in range(self.start,self.filemax-self.step,self.step):
			if i % 50 == 0:
				print("files number {:d} and {:d}.".format(i, i + self.step))
			f1 = self.allfiles[i]
			f2 = self.allfiles[i + self.step]
			btensor = BTensor(f1,f2)
			# Careful averaging: See notes. The dealing with conserved links is already done inside the BTensor.
			# The prefactor is the time-symmetrised z of the cell. Undo that, and divide by 2*linktotal as other tensors. symmetrise that too.
			zav = 0.5*(self.Neigh[u,:]+self.Neigh[u+1,:])
			linktotal=np.sum(zav)
			bavg = np.einsum('ijk,i->jk',btensor.B.T[self.labels,:,:],zav)/(linktotal)
			# We will also need the C tensor later. Same averaging.
			cavg = np.einsum('ijk,i->jk',btensor.C.T[self.labels,:,:],zav)/(linktotal)
			self.Bavg[u,:,:] = bavg # time step is now taken care of inside the tensor
			self.Cavg[u,:,:] = cavg
			if plotvtk:
				bvtk = "btensor_{:08d}.vtp".format(i)
				btensor.plot_vtk_tensor(bvtk)
			u+=1
		#self.Tenvsals, self.Tensvecs = np.linalg.eigh(self.Tensavg)
		
	# Finally: T-tensor. This is so that we can coarse-grain the plastic contribution.
	# Also collect individual T1 events, of course.
	# Averaging this one: see notes
	def statTtensor(self, plotvtk = False, verbose = True):
		print("Processing T tensor")
		self.Tavg = np.empty((self.Nsnap-1, 2, 2))
		self.T1angle_pos = []
		self.T1angle_neg = []
		self.NT1 = np.zeros((self.Nsnap,))
		u=0
		for i in range(self.start,self.filemax-self.step,self.step):
			if i % 50 == 0:
				print("files number {:d} and {:d}.".format(i, i + self.step))
			f1 = self.allfiles[i]
			f2 = self.allfiles[i+self.step]
			ttensor = TTensor(f1,f2)
			# Careful averaging: See notes. The dealing with appearing and disappearing links already done inside the TTensor.
			# The prefactor is the time-symmetrised z of the cell. Undo that, and divide by 2*linktotal as other tensors. symmetrise that too.
			zav = 0.5*(self.Neigh[u,:] + self.Neigh[u+1,:])
			linktotal=np.sum(zav)
			tavg = np.einsum('ijk,i->jk',ttensor.T.T[self.labels,:,:], zav)/(linktotal)
			self.Tavg[u,:,:] = tavg # time scale now taken care of in T tensor

			if plotvtk:
				tvtk = "ttensor_{:08d}.vtp".format(i)
				ttensor.plot_vtk_tensor(tvtk)
			
			# Now do T1 statistics
			# First locate nonzero ones
			for j in self.labels:
				tt = ttensor.T.T[j]
				e, v = np.linalg.eigh(tt)
				tr = np.trace(tt)
				# diagonalise the T tensor, the sign of the trace tells us is there are production of new links (+) or a destruction term (-)
				# Only true if done per cell, else tr(T)=0 over a region that contains all of the cells involved in the T1!
				if tr  < 0:
					self.T1angle_neg.append(np.arctan2(v[1,0],v[0,0]))
					self.NT1[u] +=1
				elif tr > 0:
					self.T1angle_pos.append(np.arctan2(v[1,1],v[0,1]))
					self.NT1[u] +=1
			# This quadruple counts T1s: Every appearing / disappearing link shows up twice in the list. Also, every T1 has both an appearing and a disappearing link.
			# Note that fractional T1 numbers are possible if some of the 4 cells involved are part of the region and if some aren't
			if self.NT1[u] > 0:
				print("Found {:.4f} T1s between steps {:d} and {:d}.".format(self.NT1[u]/4, i, i + self.step))
			u+=1
			
			
##################### Continuum plastic and elastic strain functions ###################

	def getU(self, reference = 'default', plot = True):
		if reference == 'default':
			uref = 0 # element 0 in our matrix
		else:
			uref = int((reference-self.start)/self.step)
		print ("Starting from reference element {:d} corresponding to frame.".format(uref, uref*self.step + self.start))
		
		# Definition of U = 1/2 ( log M - log M0), where log M = R^(-1) diag (log lambda) R
		# 
		M0vals, M0vecs = np.linalg.eigh(self.Mavg[uref,:,:])
		logM0 = np.dot(np.transpose(M0vecs),np.dot(np.diag(M0vals),M0vecs))
		
		# And now the full matrices
		self.U = np.zeros((self.Nsnap,2,2))
		for u in range(uref,self.Nsnap):
			Mvals, Mvecs = np.linalg.eigh(self.Mavg[u,:,:])
			logM = np.dot(np.transpose(Mvecs),np.dot(np.diag(Mvals),Mvecs))
			self.U[u,:,:] = 0.5*(logM - logM0)
			
		if plot:
			plt.figure()
			# Plot the different components. Note this is a symmetric tensor
			plt.plot(self.tval,self.U[:,0,0],'-r',label='Uxx')
			plt.plot(self.tval,self.U[:,1,1],'-k',label='Uyy')
			plt.plot(self.tval,self.U[:,0,1],'-g',label='Uxy')
			plt.xlabel('time')
			plt.ylabel('Statistical internal strain U')
			plt.legend()
			
	def getP(self, plot = True):
		self.P = np.zeros((self.Nsnap-1,2,2))
		for u in range(self.Nsnap-1):
			minv = np.linalg.inv(self.Mavg[u,:,:])
			self.P[u,:,:] = 0.5*(np.dot(minv,self.Tavg[u,:,:]) + np.dot(self.Tavg[u,:,:],minv))
			
		if plot:
			plt.figure()
			# Plot the different components. Note this is a symmetric tensor
			plt.plot(self.tval[:(self.Nsnap-1)],self.P[:,0,0],'-r',label='Pxx')
			plt.plot(self.tval[:(self.Nsnap-1)],self.P[:,1,1],'-k',label='Pyy')
			plt.plot(self.tval[:(self.Nsnap-1)],self.P[:,0,1],'-g',label='Pxy')
			plt.xlabel('time')
			plt.ylabel('Topological rearrangement rate P')
			plt.legend()
			
			
	def getV(self, plot = True):
		self.V = np.zeros((self.Nsnap-1,2,2))
		for u in range(self.Nsnap-1):
			minv = np.linalg.inv(self.Mavg[u,:,:])
			self.V[u,:,:] = 0.5*(np.dot(minv,self.Cavg[u,:,:]) + np.dot(np.transpose(self.Cavg[u,:,:]),minv))
			
		if plot:
			plt.figure()
			# Plot the different components. Note this is a symmetric tensor
			plt.plot(self.tval[:(self.Nsnap-1)],self.V[:,0,0],'-r',label='Vxx')
			plt.plot(self.tval[:(self.Nsnap-1)],self.V[:,1,1],'-k',label='Vyy')
			plt.plot(self.tval[:(self.Nsnap-1)],self.V[:,0,1],'-g',label='Vxy')
			plt.xlabel('time')
			plt.ylabel('Symmetrised velocity gradient V')
			plt.legend()
			
	def getOmega(self, plot = True):
		self.Omega = np.zeros((self.Nsnap-1,2,2))
		for u in range(self.Nsnap-1):
			minv = np.linalg.inv(self.Mavg[u,:,:])
			self.Omega[u,:,:] = 0.5*(np.dot(minv,self.Cavg[u,:,:]) - np.dot(np.transpose(self.Cavg[u,:,:]),minv))
			
		if plot:
			plt.figure()
			# Plot the different components. Note this is a symmetric tensor
			plt.plot(self.tval[:(self.Nsnap-1)],self.Omega[:,0,0],'-r',label='Oxx')
			plt.plot(self.tval[:(self.Nsnap-1)],self.Omega[:,1,1],'-k',label='Oyy')
			plt.plot(self.tval[:(self.Nsnap-1)],self.Omega[:,0,1],'-g',label='Oxy')
			plt.xlabel('time')
			plt.ylabel('Rotational velocity gradient Omega')
			plt.legend()
			
################ Basic statistics ####### limit bins given data issues
	def getBasestats(self,nbin=40):
		# mean values
		avArea = np.zeros((self.Nsnap,))
		avPerim = np.zeros((self.Nsnap,))
		avp0 = np.zeros((self.Nsnap,))
		# Some histograms for quantities
		p0bin=np.linspace(3.6,4.6,nbin+1)
		dp0=p0bin[2]-p0bin[1]
		myobin=np.linspace(0.0,2.0,nbin+1)
		dmy=myobin[1]-myobin[0]
		tensbin=np.linspace(0.0,2.0,nbin+1)
		dten=tensbin[1]-tensbin[0]
		p0hist=np.zeros((self.Nsnap,nbin))
		myohist=np.zeros((self.Nsnap,nbin))
		tenshist=np.zeros((self.Nsnap,nbin))
		for i in range(self.start,self.filemax,self.step):
			if i % 50 == 0:
				print("file number {:d}.".format(i))
			# get values
			bst = BasicStatistics(self.allfiles[i])
			area = bst.getAreas()
			perimeter = bst.getPerimeters()
			tension = bst.getEdgeTension(self.labels)
			myosin = bst.getEdgeMyosin(self.labels)
			# produce statistics
			p0 = perimeter[self.labels]/np.sqrt(area[self.labels])
			avArea[i] = np.mean(area[self.labels])
			avPerim[i] = np.mean(perimeter[self.labels])
			avp0[i] = np.mean(p0)
			p0hist[i,:],edges=np.histogram(p0,bins=p0bin,normed=True)
			myohist[i,:],edges=np.histogram(myosin,bins=myobin,normed=True)
			tenshist[i,:],edges=np.histogram(tension,bins=tensbin,normed=True)
		# return the results
		return avArea, avPerim, avp0, p0bin, p0hist, myobin, myohist, tensbin, tenshist
		
			
		
				
	
			
	
	
	

		
