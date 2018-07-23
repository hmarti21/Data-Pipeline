
# coding: utf-8

# In[1]:


import sys
sys.path.append('/global/homes/h/hmarti21/Pipeline')
import particleFileIO
import corFunc
import numpy as np
from matplotlib import pyplot as plt
from Corrfunc.theory.wp import wp
from astropy.table import Table
import h5py

import importlib
reload=importlib.reload
import corFunc
reload(corFunc)
import corFunc

filename='/global/u2/h/hmarti21/data_sukhdeep.hdf5'
col_names={'x':'halos.x','y':'halos.y','z':'halos.z','mass':'halos.mass','q':'shapesStar.q2d','ax':'shapesStar.a2d_x','ay':'shapesStar.a2d_y','bx':'shapesStar.b2d_x','by':'shapesStar.b2d_y','e1':None,'e2':None}
cuts='/global/u2/h/hmarti21/cuts.txt'
sightBins=5
rscale=1
nbins=10
min_sep=1
max_sep=10
min_box=0
max_box=100
rpar_min=0
rpar_max=50
key='hydro_full'
logfile='/global/u2/h/hmarti21/log.txt'
savefile='/global/u2/h/hmarti21/TreeCorrData.txt'

#data=corFunc.read_Data_hdf5(filename,col_names,cuts,key,logfile)
#rbins=np.logspace(np.log10(min_sep/1000),np.log10(max_sep/1000),nbins+1)
#data['x']=data['x']/1000
#data['y']=data['y']/1000
#data['z']=data['z']/1000
#nthreads=2
#pimax=30.
#x,y,z=np.array(data['x']),np.array(data['y']),np.array(data['z'])
#results_wp=wp(max_box/1000,pimax,nthreads,rbins,x,y,z)
# results_wp = wp(boxsize, pimax, nthreads, rbins, X, Y, Z)

filepath='/global/cscratch1/sd/sukhdeep/mb2_subsample/'
#col_def=[('Position', ('f8', 3), 'all'), ('Mass','auto',None)]
ptype=1
snapArray=[85]

#particleFileIO.particleFile(filepath,col_def,fraction,num_files,ptype)
particles=particleFileIO.readFITSFile('/global/cscratch1/sd/sukhdeep/mb2_subsample/',ptype,snapArray)
#NNPairs=np.zeros((sightBins,nbins))
#NRPairs=np.zeros((sightBins,nbins))
#RRPairs=np.zeros((sightBins,nbins))
#RNPairs=np.zeros((sightBins,nbins))
#print(filename,col_names,cuts,sightBins,rscale,nbins,min_sep,max_sep,rpar_step,min_box,max_box,logfile,savefile,key)
corFunc.corFunc(filename,col_names,cuts,sightBins,rscale,nbins,min_sep,max_sep,rpar_min,rpar_max,min_box,max_box,'NG',logfile,savefile,key=key,fname2='yes',dat2=particles)

