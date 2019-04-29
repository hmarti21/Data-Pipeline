
# coding: utf-8

# In[1]:


import sys
sys.path.append('/global/homes/h/hmarti21/Pipeline')
import particleFileIO
import corFunc
import numpy as np
from matplotlib import pyplot as plt
from Corrfunc.theory.wp import wp
from astropy import table
import h5py

import importlib
reload=importlib.reload
import corFunc
reload(corFunc)
import corFunc

print(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6])
filename='/global/u2/h/hmarti21/dataFrameTest'
col_names={'x':8,'y':9,'z':10,'mass':0,'q':11,'ax':12,'ay':13,'bx':14,'by':15,'e1':None,'e2':None}
cuts=sys.argv[1]  #'/global/u2/h/hmarti21/cuts.txt'
sightBins=5
rscale=np.int(sys.argv[2])
nbins=np.int(sys.argv[3])
min_sep=np.float(sys.argv[4])
max_sep=10
min_box=0
max_box=100
rpar_min=0
rpar_max=50

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
ptype=0
snapArray=[85]

#particleFileIO.particleFile(filepath,col_def,fraction,num_files,ptype)
#particles=particleFileIO.readFITSFile('/global/cscratch1/sd/sukhdeep/mb2_subsample/',ptype,snapArray)
#NNPairs=np.zeros((sightBins,nbins))
#NRPairs=np.zeros((sightBins,nbins))
#RRPairs=np.zeros((sightBins,nbins))
#RNPairs=np.zeros((sightBins,nbins))
#print(filename,col_names,cuts,sightBins,rscale,nbins,min_sep,max_sep,rpar_step,min_box,max_box,logfile,savefile,key)
#corFunc.corFunc(filename,col_names,cuts,sightBins,rscale,nbins,min_sep,max_sep,rpar_min,rpar_max,min_box,max_box,'NG',logfile,savefile,key=key,fname2='yes',dat2=particles)

ptype=[0,1,4]
combine=np.bool(sys.argv[5])
dir=sys.argv[6]

print(cuts,rscale,nbins,min_sep,combine,dir)
if combine==True:
    print('hello')
    for j in snapArray:
        particles=table.Table()
        particles['x']=[]
        particles['y']=[]
        particles['z']=[]
        particles ['mass']=[]
        key='snap_'+str(j)
        savefile='/global/homes/h/hmarti21/Results/'+dir+'/TreeCorrData_'+str(j)+'combine.txt'
        logfile='/global/homes/h/hmarti21/Logs/'+dir+'/log_'+str(j)+'combine.txt'
        print(savefile,logfile)
        for i in ptype:
            particle=particleFileIO.readFITSFile('/global/cscratch1/sd/sukhdeep/mb2_subsample/',i,j)
            particles=table.vstack([particles,particle])
        print (particles)
        corFunc.corFunc(filename,col_names,cuts,sightBins,rscale,nbins,min_sep,max_sep,rpar_min,rpar_max,min_box,max_box,'NN',logfile,savefile,key=key,fname2='yes',dat2=particles)

if combine==False:
    for j in ptype:
        for i in snapArray:
            key='snap_'+str(i)
            savefile='/global/homes/h/hmarti21/Results/'+dir+'/TreeCorrData_'+str(i)+'ptype_'+str(j)+'.txt'
            logfile='/global/homes/h/hmarti21/Logs/'+dir+'/log_'+str(i)+'ptype_'+str(j)+'.txt'
            particles=particleFileIO.readFITSFile('/global/cscratch1/sd/sukhdeep/mb2_subsample/',j,i)    
            corFunc.corFunc(filename,col_names,cuts,sightBins,rscale,nbins,min_sep,max_sep,rpar_min,rpar_max,min_box,max_box,'NG',logfile,savefile,key=key,fname2='yes',dat2=particles)
