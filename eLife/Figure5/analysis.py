# Script retrieving, refactorizing and computing statistics
# over simulations replica into an easily processible and compact format

#regular imports
import numpy as np
import glob as glob
import pickle as pickle
import sys as sys
import re
import glob
import matplotlib.pyplot as plt
import re
import pickle
import traceback 
from os import path

def moving_average(a, n=20) :
    '''A rolling average function which averages over n points'''
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    # we also need the first element 
    return ret[n - 1:] / n  

def int_graner_err(data=None,start=100,roll=False,mask=...):
    '''returns the standard deviation for the integrated strain rates'''
    #broadcast requires extra dimensions
    dt = np.diff(data['tval'][mask,start:,np.newaxis,np.newaxis],axis=1)
    intg = np.cumsum(data['V'][mask,start:,:,:]*dt,axis=1) #integral
    v_mean = np.nanmean(intg,axis=0) #sum V over all runs 
    v_std = np.nanstd(intg,axis=0) #std deviation of V over runs
    vdiff_std = np.nanstd(intg[...,1,1]-intg[...,0,0],axis=0) #std deviation of the difference

    # Rolling average version 
    if roll : 
        dims = np.asarray(v_mean.shape) #its a tuple 
        dims[0] -= n-1
        mean_roll = np.zeros(dims) # create array of right dimension 
        # loop over dimensions
        for i in range(2):
            for j in range(2):
                mean_roll[:,i,j] = moving_average(v_mean[:,i,j])
        txx = np.nanargmin(mean_roll[:,0,0])
        tyy = np.nanargmax(mean_roll[:,0,0])
        tdiff = np.nanargmax(mean_roll[:,1,1] - mean_roll[:,0,0])
        vxx_min =  mean_roll[txx]
        vyy_max =  mean_roll[tyy]
        vxx_min_std = v_std[txx]
        vyy_max_std = v_std[tyy]
        vdiff_max =  (mean_roll[:,1,1] - mean_roll[:,0,0])[tdiff]
        vdiff_max_std = vdiff_std[diff]
    else :
        txx = np.nanargmin(v_mean[:,0,0])
        vxx_min =  v_mean[txx,0,0]
        vxx_min_std = v_std[txx,0,0]        
        tyy = np.nanargmax(v_mean[:,0,0])
        vyy_max =  v_mean[tyy,1,1]
        vyy_max_std = v_std[tyy,1,1]        
        tdiff = np.nanargmax(v_mean[:,1,1] - v_mean[:,0,0])
        vdiff_max = (v_mean[:,1,1] - v_mean[:,0,0])[tdiff]
        vdiff_max_std = vdiff_std[tdiff]
    return vxx_min_std,vyy_max_std,vdiff_max_std

def int_graner(data=None,start=100,roll=False,mask=...):
    ''' returns the extremum of the integrated strain rate'''        
    dt = np.diff(data['tval'][mask,start:,np.newaxis,np.newaxis],axis=1) #broadcast requires extra dimensions
    intg = np.cumsum(data['V'][mask,start:,:,:]*dt,axis=1) #integral over time 
    v_mean = np.nanmean(intg,axis=0) # sum over runs
    if roll : 
        dims = np.asarray(v_mean.shape) #its a tuple 
        dims[0] -= n-1
        mean_roll = np.zeros(dims) # create array of right dimension 
        # loop over dimensions
        for i in range(2):
            for j in range(2):
                mean_roll[:,i,j] = moving_average(v_mean[:,i,j])
                vxx_min =  np.min(mean_roll[:,0,0])
                vyy_max =  np.max(mean_roll[:,1,1])
                v_diff =  np.nanmax(mean_roll[:,1,1] - mean_roll[:,0,0])
    else :
        vxx_min =  np.min(v_mean[:,0,0])
        vyy_max =  np.max(v_mean[:,1,1])
        v_diff =  np.nanmax(v_mean[:,1,1] - v_mean[:,0,0])
    return vxx_min,vyy_max,v_diff

def time_to_t1(data,start):
    ''' compute time to central T1 statistics for a set of data '''
    NT1 = data['NT1']
    ts = data['tval'] - data['tval'][:,start][...,np.newaxis]
    t1ts = np.empty(NT1.shape[:-1])
    indx = []
    for i,d in enumerate(NT1):
        t1i = np.where(d>0) # find T1 events, might be empty 
        t1t = ts[i][t1i]  #convert in real time 
        t1ts[i] = t1t[0] if t1t.size > 0 else np.NaN # store first T1 time 
        if t1t.size > 0 :
            indx.append(i)
    # then return statistics 
    return t1ts,np.nanmean(t1ts),np.nanstd(t1ts),indx 

def get_ft1_data(data=None,start=100,mask=...):
    ''' the optional mask argument allows to filter the runs we want to compute ''' 
    ft1_raw,ft1,ft1d,indx = time_to_t1(data,start)
    ft1,ft1d = ft1[mask],ft1d[mask]
    ft1op_raw = np.array([np.arcsin(np.sin(t1op[0])) if len(t1op) > 0 else np.NaN for t1op in data['T1angle_pos']])
    ft1op = np.nanmean(ft1op_raw[mask])
    ft1opd = np.nanstd(ft1op_raw[mask])
    ft1on_raw = np.array([np.arccos(np.cos(t1on[0])) if len(t1on) > 0 else np.NaN for t1on in data['T1angle_neg']])
    ft1on = np.nanmean(ft1on_raw[mask])
    ft1ond = np.nanstd(ft1on_raw[mask])
    return ft1_raw,ft1,ft1d,ft1op_raw[mask],ft1op,ft1opd,ft1on_raw[mask],ft1on,ft1ond,indx


def get_fct1_data(data=None,data2=None,start=100):
    ''' Retrieves first central and non central T1 
    T1 are filtered as follows :
    1) A first mask "mask" determinates if there are even T1 events
    2) A second mask "mask2" then attribute a value to T1 events :
       - -1 If no central T1 
       - 0 If central T1 is after the first T1 
       - +1 If the central T1 is the first T1 
    ''' 
    n = data['n_runs']
    data['t1_type'] = np.full(n,-1)
    mask =  np.full(n,False) # used to ignore non central "clean" T1s
    fct1 =  np.full(n,np.NaN)  # collapse time for central T1 (be it first or not)
    for i in range(n) :
        cj_ft1 = data['firstT1'][i]
        all_t1s = data2['NT1'][i]        
        mask[i] = False #mask value in case no T1 at all 
        t1_type = -1
        
        if np.argwhere(all_t1s>0).size == 0:
            mask[i] = False   # no T1 at all 
        else :
            mask[i] = True # there are T1s, don't mask 
            
            if cj_ft1 == -1  :
                data['t1_type'][i] = -1 #no central T1 
                
            elif cj_ft1 - 1 == np.argwhere(all_t1s>0)[0]: #compare first recorded T1s
                data['t1_type'][i] = +1 #central and first t1 match
                
            else :
                data['t1_type'][i] = 0 #central is after first T1                
    mask2 = data['t1_type'] >= 1 # only consider when the first t1 is central
    idxs = data['firstT1'][mask2] 
    fct1 = data2['tval'][mask2,idxs] - data2['tval'][mask2,start] # transform in real time
    return np.mean(data['t1_type'][mask]),mask2,fct1,np.mean(fct1),np.std(fct1),len(fct1) 


bdir = "/media/ilyas/Data/Sims/Unified/4Cells/Myosin_pool/Ref/PS_taul_r/"
beta = 1
regex = '{}/B{}/Tl*R*'.format(bdir,beta)
files = glob.glob(regex)
tau_l = np.array([ re.search('Tl(.*)_R(.*)', x).group(1) for x in files])
r = np.array([re.search('Tl(.*)_R(.*)', x).group(2) for x in files])

# Central T1 statistics
fct1s = np.NaN*np.ones(len(r)) #first central T1 time collapse time
fct1ds = np.NaN*np.ones(len(r)) #first central T1 collapse time (stddev)
fct1ns = np.NaN*np.ones(len(r)) #first central T1 number
mask_fct1s = [[] for _ in range(len(r))] # to filter central T1s 
t1_types = np.NaN*np.ones(len(r)) # "average" T1 type 

# First T1 statistics
ft1s = np.NaN*np.ones(len(r)) #first T1 collapse time
ft1ds = np.NaN*np.ones(len(r)) #first T1 collapse time (stddev)
ft1ops = np.NaN*np.ones(len(r)) #first T1 + EV orientation
ft1opds = np.NaN*np.ones(len(r)) # std dev
ft1ons = np.NaN*np.ones(len(r))  #first T1 - EV orientation
ft1onds = np.NaN*np.ones(len(r)) #first T1 orientation 
ft1ns = np.NaN*np.ones(len(r)) # number of T1s 

# Graner statistics
vxx_mins = np.NaN*np.ones(len(r)) 
vyy_maxs = np.NaN*np.ones(len(r))
v_diffs = np.NaN*np.ones(len(r))
vxx_mins_std = np.NaN*np.ones(len(r))
vyy_maxs_std = np.NaN*np.ones(len(r))
v_diffs_std = np.NaN*np.ones(len(r))


#Raw data
fct1s_raw =  [[] for _ in range(len(r))] # first central T1 time
ft1s_raw =  [[] for _ in range(len(r))] # first T1 time
ft1ops_raw =  [[] for _ in range(len(r))] # first T1 + EV orientation
ft1ons_raw =  [[] for _ in range(len(r))] # first T1 - EV orientation


for i,(R,TAUL) in enumerate(zip(r,tau_l)):
    # folders containing the central T1 data 
    folder = '{}/B{}/Tl{}_R{}/'.format(bdir,beta,TAUL,R) 
    regex = '{}/run*/junction_dynamics_timeseries.p'.format(folder)
    f_jcs = glob.glob(regex)
    n = len(f_jcs)
    
    if n == 0:
        print("no junction_dynamics_timeseries files in {}".format(folder))
        continue
    else: 
        print("Treating {} junction_dynamics_timeseries in folder {}".format(n,folder))

    data_jcs = [pickle.load(open(f,'rb')) for f in f_jcs] #take data
    data_jcs = {k:np.array([d[k] for d in data_jcs]) for k in data_jcs[0].keys()} #make data addressable by keys
    data_jcs['n_runs'] = n
    
    # take corresponding Graner tensor folders (ie in the same order)
    f_tsrs = [ f.replace('junction_dynamics_timeseries.p','tensor_timeseries.p') for f in f_jcs ]
    n2 = len([ True  for f in f_tsrs if path.isfile(f)])
    
    if n2 == 0:
        print("no tensor_timeseries files in {}".format(folder))
        continue

    elif n2  != n:
        print("Found {} tensor_timeseries  which differs from {} junction_dynamics_timeseries, skipping folder {}".format(n2,n,folder))
        continue
    else:
        print("Treating {} tensor_timeseries in folder {}".format(n2,folder))
        
    data_tsrs = [pickle.load(open(f,'rb')) for f in f_tsrs] #take data
    data_tsrs = {k:np.array([d[k] for d in data_tsrs]) for k in data_tsrs[0].keys()} #make data addressable by keys
    
    # Get T1 types    
    t1_type,mask_fct1,fct1_raw,fct1,fct1d,fct1n = get_fct1_data(data_jcs,data_tsrs)
    t1_types[i] = t1_type
    mask_fct1s[i] = mask_fct1
    fct1s_raw[i] =  fct1_raw 
    fct1s[i] = fct1
    fct1ds[i] = fct1d
    fct1ns[i] = fct1n
    
     
    # Get any first T1 collapse time
    ft1_raw,ft1,ft1d,ft1op_raw,ft1op,ft1opd,ft1on_raw,ft1on,ft1ond, indx = get_ft1_data(data_tsrs)
    ft1s_raw[i] =  ft1_raw
    ft1s[i]= ft1
    ft1ds[i]= ft1d
    ft1ops[i]= ft1op
    ft1opds[i]= ft1opd
    ft1ops_raw[i] = ft1op_raw
    ft1ons[i]= ft1on
    ft1onds[i]= ft1ond
    ft1ons_raw[i] = ft1on_raw
    ft1ns[i] = len(indx)     
    
    # Get Graner extremum
    vxx_min,vyy_max,v_diff = int_graner(data_tsrs,start=100,mask=mask_fct1)
    vxx_mins[i] = vxx_min
    vyy_maxs[i] = vyy_max
    v_diffs[i] = v_diff    

    
    # Get graner extremum std deviation
    vxx_min_std,vyy_max_std,v_diff_std = int_graner_err(data_tsrs,start=100,mask=mask_fct1)
    vxx_mins_std[i] = vxx_min_std
    vyy_maxs_std[i] = vyy_max_std
    v_diffs_std[i] = v_diff_std

    
r = np.asarray(r,float)
tau_l = np.asarray(tau_l,float)
n_runs = np.array([len(m) for m in mask_fct1s ])

with open('data_ft1.p','wb') as f :
    pickle.dump({'r':r,'tau_l':tau_l,'n_runs':n_runs, # parameters 
                 'time_ft1':(ft1s,ft1ds),'orientation_ft1_neg':(ft1ops,ft1opds),'orientation_ft1_pos':(ft1ons,ft1onds),'n_t1':ft1ns, # first T1 (central or not) stats 
                 'vxx_mins':(vxx_mins,vxx_mins_std),'vyy_maxs':(vyy_maxs,vyy_maxs_std),'v_diffs':(v_diffs,v_diffs_std), # Graner stats
                 't1_types':t1_types,'time_fct1':(fct1s,fct1ds),'n_fct1':fct1ns,'mask_fct1':mask_fct1s, # First central T1 stats
                 'raw_data':(ft1s_raw,fct1s_raw,ft1ops_raw,ft1ons_raw) #raw data
                 }
                ,f)
