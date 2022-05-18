# Script for plotting statistics of the processed AJM data from analysis.py
import numpy as np
import pickle
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm
from copy import deepcopy 
FIG_FOLDER = "./Figs/"
plt.rcParams["figure.figsize"] = [4.8,4.8]


def plot_slice_r(r,y,yerr,cond=True,c=None,ls="-",marker="o"):
    '''loops at fixed r value. Optional cond arguments allows for extra constraints''' 
    
    # convert to float
    tau_l = data['tau_l'] 
    # ignore nan values
    m1 = np.isfinite(y) 
    # take corresponding r values
    m2 = data['r'] == r 
    m = cond & m1 & m2
    y,yerr,tau_l = y[m],yerr[m],tau_l[m]

    # need to sort by ascending tau_l
    sortt = tau_l.argsort() #keep sorting order tau_l
    y,yerr,tau_l = y[sortt],yerr[sortt],tau_l[sortt]

    # then plot
    ln = plt.errorbar(tau_l,y,yerr,marker=marker,linestyle=ls,color=c)
    return ln


def plot_slice_r_fill(r,y,yerr,cond=True,c=None,ls="-"):
    '''loops at fixed r value. Optional cond arguments allows for extra constraints''' 
    
    # convert to float
    tau_l = data['tau_l'] 
    # ignore nan values
    m1 = np.isfinite(y) 
    # take corresponding r values
    m2 = data['r'] == r 
    m = cond & m1 & m2
    y,yerr,tau_l = y[m],yerr[m],tau_l[m]

    # need to sort by ascending tau_l
    sortt = tau_l.argsort() #keep sorting order tau_l
    y,yerr,tau_l = y[sortt],yerr[sortt],tau_l[sortt]

    # then plot
    plt.fill_between(tau_l, y-yerr, y+yerr,linestyle=ls,color=c,alpha=.2)

    
def make_regular_grid(tvs,rs):
    ''' Returns the maximal set of tau_l allowing to get a regular grid ''' 
    x = tvs #maximum set 
    y = rs 

    # Strip values which don't exist in one of the r=cst line
    for yc in y :
        idx =  data['r'] == yc
        tv = data['tau_l'][idx]
        x2 = np.unique(tv)
        x = np.intersect1d(x,x2)    
    x.sort()
    return x,y

        
def grid(x,y,f) :
    ''' Find corresponding points for array "f" on a grid based on x and y values '''
    c = np.full((len(x),len(y)),np.nan)
    for i,xc in enumerate(x) :
        for j,yc in enumerate(y) :
            c1 = data['tau_l'] == xc # get the value
            c2 = data['r'] == yc # get the value
            c3 = np.where( c1 & c2)
            val = f[c3]
            c[i][j] = val            
    return c

with open('data_ft1.p','rb') as f :
    # this pickle file contains:
    # # Parameter values: 
    #   tau_l,tau_r
    # - Information relative to the T1: 
    #    time to first T1(+stddev), first positive eigenvalue (+stddev)
    #   angle of first negative eigenvalue (+stddev)
    # - The extremum of the strain (integrated strain-rate)
    #    max(\int V_yy) - min(\int V_xx), min(\int V_xx) , max(\int V_yy)
    # - Information relative to the first *central* T1:
    # the type of T1, the collapse time(+stddev), the rate of collapse vs non-collapse
    # - Some raw data:
    # first T1 time(+sttdev), first central T1 time, positive and negative ev orientation 
    data = pickle.load(f)


# find unique r,tau_l values
rs  = np.unique(data['r'])
rs = rs[:-1] #discard some values
tvs = np.unique(data['tau_l'])
discard_no_ct1 = (data['n_fct1']/data['n_runs']) > 0.25 # take only those with 25% of T1
# Define colormap 
color = plt.cm.nipy_spectral(np.linspace(0,1,len(rs)+1)) 


    
# Plot the central t1 collapse time  (WITH errorbar)
with PdfPages(FIG_FOLDER+'fig_5b.pdf') as pdf :
    fig = plt.figure()
    plt.xlabel(r'$\tau_v$')  
    plt.ylabel(r'$\tau_c$')
    lns = []
    for r,c in zip(rs,color):
        ln = plot_slice_r(r,data['time_fct1'][0]
                          ,np.full_like(data['time_fct1'][0],np.NaN)
                          ,c=c,ls='-',cond=discard_no_ct1)
        plot_slice_r_fill(r,data['time_fct1'][0]
                     ,data['time_fct1'][1]/np.sqrt(data['n_fct1'])
                     ,c=c,ls='-',cond=discard_no_ct1)
        ln.set_label(r'${:d}$'.format(round(1/r)))
        
    # add 1/2 slope
    ln = plt.plot(np.linspace(min(tvs),max(tvs),10),np.sqrt(np.linspace(min(tvs),max(tvs),10))*10,ls='--',color='black',label=r'$\propto \sqrt{\tau_m}$')

    plt.legend(ncol=2,loc="upper left")
    
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(min(tvs),max(tvs))
    pdf.savefig()


# --------------------- Strain ----------------------     
with PdfPages(FIG_FOLDER+'fig_5c.pdf') as pdf :
    fig = plt.figure()
    #plt.title(r'$max(ϵ_{yy}-ϵ_{xx})$ vs $\tau_v$')
    plt.xlabel(r'$\tau_v$')
    plt.ylabel(r'$ϵ_{yy}-ϵ_{xx}$')
    plt.xlim(min(tvs),max(tvs))
    plt.ylim(0,0.5)
    for r,c in zip(rs,color):
        plot_slice_r(r,data['v_diffs'][0],np.full_like(data['v_diffs'][0],np.NaN),c=c,ls='-',marker='.',cond=discard_no_ct1)
        plot_slice_r_fill(r,data['v_diffs'][0],data['v_diffs'][1]/np.sqrt(data['n_fct1']),c=c,cond=discard_no_ct1)
    plt.legend([ r'${:d}$'.format(round(1/r)) for r in rs],loc='lower right',ncol=2)

    plt.xscale('log')
    plt.yscale('linear')
    pdf.savefig()
    plt.close()


# Probability of T1
# Find our gridpoints: Just in case we don't have a regular grid, 
x,y = make_regular_grid(tvs,rs) 

with PdfPages(FIG_FOLDER+'fig_5a.pdf') as pdf :
    c0 = grid(x,y,data['n_fct1']/data['n_runs'])
    c1 = grid(x,y,np.array(data['n_t1']/data['n_runs']))

    for c,name in zip([c0,c1],[r'$\tau^{central}_c$',r'$\tau_c$']):    
        cmap = cm.get_cmap()
        cmap.set_bad(color='red')

        plt.clf()
        fig,ax = plt.subplots(1,figsize=(5.8,4.8))

        ax.set_xlabel(r'$\tau_v$')
        ax.set_ylabel(r'$\tau_m$')
        ax.set_xscale('log')
        ax.set_yscale('log')

        pcm0 = ax.pcolormesh(x,1/y,c0.T,cmap=cmap)
        cbar0 = fig.colorbar(pcm0,ax=ax, orientation='vertical',label='probability central T1')
        plt.tight_layout()
        pdf.savefig()
        plt.close()
