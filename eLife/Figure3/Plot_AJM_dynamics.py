import sys as system
import argparse
import pickle
import numpy as np
import matplotlib.pyplot as plt



Nsnap  =899

picklefile = 'junction_dynamics_timeseries.p'
print(picklefile)
data  = pickle.load(open(picklefile,'rb'))
# data.update({'isdata':isdata,'myosins':myosins,'tensions':tensions,'l0s':l0s,'lengths':lengths})
myosins = data['myosins']
tensions = data['tensions']
lengths = data['lengths']
l0s = data['l0s']

tval = np.linspace(0,Nsnap,Nsnap)-100

plt.figure()
plt.plot(tval,myosins[:,0],'-g',label='myosin',lw=2)
plt.plot(tval,myosins[:,1],'-g',lw=2)
plt.plot(tval,tensions[:,0],'-r',label='tension',lw=2)
plt.plot(tval,lengths[:,0],'-k',label='length',lw=2)
plt.plot(tval,l0s[:,0],'-b',label='l0',lw=2)
plt.xlabel('time')
plt.ylabel('m,T,l,l_0')
plt.legend()
plt.title('Central junction')
plt.xlim(0,206)
plt.ylim(0,1.75)

#plt.figure()
#plt.plot(tval,myosins[:,2],'-g',label='myosin',lw=2)
#for u in range(3,10):
	#plt.plot(tval,myosins[:,u],'-g',lw=2)
#plt.plot(tval,tensions[:,1],'-r',label='tension',lw=2)
#for u in range(2,5):
	#plt.plot(tval,tensions[:,u],'-r',lw=2)
#plt.plot(tval,lengths[:,1],'-k',label='length',lw=2)
#for u in range(2,5):
	#plt.plot(tval,lengths[:,u],'-k',lw=2)
#plt.plot(tval,l0s[:,1],'-b',label='l0',lw=2)
#for u in range(2,5):
	#plt.plot(tval,l0s[:,u],'-b',lw=2)
#plt.xlabel('time')
#plt.ylabel('value')
#plt.legend()
#plt.title('Shoulder junctions')
#plt.xlim(0,206)
#plt.ylim(0,2.3)

plt.figure()
shoulder1=[2,5,6,9]
myav = np.average(myosins[:,shoulder1],axis=1)
dmyo = np.std(myosins[:,shoulder1],axis=1)
plt.plot(tval,myav,'-g',label='myosin',lw=2)
plt.fill_between(tval, myav-dmyo, myav+dmyo, color='g', alpha=.2)

shoulder2=[3,4,7,8]
myav = np.average(myosins[:,shoulder2],axis=1)
dmyo = np.std(myosins[:,shoulder2],axis=1)
plt.plot(tval,myav,'-g',label='myosin',lw=2)
plt.fill_between(tval, myav-dmyo, myav+dmyo, color='g', alpha=.2)

tensav = np.average(tensions[:,1:5],axis=1)
dtens = np.std(tensions[:,1:5],axis=1)
plt.plot(tval,tensav,'-r',label='tension',lw=2)
plt.fill_between(tval, tensav-dtens, tensav+dtens, color='r', alpha=.2)

lenav = np.average(lengths[:,1:5],axis=1)
dlen = np.std(lengths[:,1:5],axis=1)
plt.plot(tval,lenav,'-k',label='length',lw=2)
plt.fill_between(tval, lenav-dlen, lenav+dlen, color='k', alpha=.2)

l0av = np.average(l0s[:,1:5],axis=1)
dl0 = np.std(l0s[:,1:5],axis=1)
plt.plot(tval,l0av,'-b',label='l0',lw=2)
plt.fill_between(tval, l0av-dl0, l0av+dl0, color='b', alpha=.2)

plt.xlabel('time')
plt.ylabel('m,T,l,l_0')
plt.legend()
plt.title('Shoulder junctions')
plt.xlim(0,206)
plt.ylim(0,2.3)




plt.show()

