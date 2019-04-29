
# coding: utf-8

# In[1]:

import time
import sys
sys.path.append('/global/homes/h/hmarti21/Pipeline')
import particleFileIO
import corFunc
import galaxyFileIO
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

filename='/global/u2/h/hmarti21/dataFrameTest'
col_names={'x':8,'y':9,'z':10,'mass':0,'q':11,'ax':12,'ay':13,'bx':14,'by':15,'e1':None,'e2':None}
cuts=sys.argv[1]  #'/global/u2/h/hmarti21/cuts.txt'
sightBins=5
rscale=np.int(sys.argv[2])
nbins=np.int(sys.argv[3])
min_sep=np.float(sys.argv[4])
max_sep=np.float(sys.argv[5])
min_box=0
max_box=100
rpar_min=0
rpar_max=50
filePath='/global/cscratch1/sd/sukhdeep/mb2_subsample/'
snapArray=[85]
ptype=[0,1,4]
combine=np.int(sys.argv[6])
dir=sys.argv[7]
func=sys.argv[8]
galFiles=np.int(sys.argv[9])
logfile=sys.argv[10]

if galFiles==0:
    func='NN'
    functype='MM'
    if combine==1:
        for j in snapArray:
            key='snap_'+str(j)
            savefile='/global/homes/h/hmarti21/Results/Current_Results/'+dir+'/DelMM_'+str(j)+'combine'
            #logfile='/global/homes/h/hmarti21/Results/Current_Results/'+dir+'/log_DelMM_'+str(j)+'combine'
            #corFunc.makeLog(logfile,func)
            particles=particleFileIO.combinedParticles(filePath,ptype,j)
            data=table.Table()
            data['x']=[]
            data['y']=[]
            data['z']=[]
            data['mass']=[]
            arr=particles
            fraction=0.1
            rand=np.random.uniform(0,1,size=len(arr['mass'][:]))
            rand_mask=rand<fraction
            arr_masked=arr[rand_mask]
            arr_masked['x']=arr_masked['x']/1000
            arr_masked['y']=arr_masked['y']/1000
            arr_masked['z']=arr_masked['z']/1000
            #print ('NumParticles:',len(particles['x']))
            corFunc.corFunc(arr_masked,sightBins,rscale,nbins,min_sep,max_sep,rpar_min,rpar_max,min_box,max_box,func,logfile,savefile,combine,functype,key=key)
    if combine==0:
        for j in ptype:
            for i in snapArray:
                key='snap_'+str(i)
                savefile='/global/homes/h/hmarti21/Results/Current_Results/'+dir+'/DelMM_'+str(i)+'ptype_'+str(j)
                #logfile='/global/homes/h/hmarti21/Results/Current_Results/'+dir+'/log_'+str(i)+'ptype_'+str(j)
                #corFunc.makeLog(logfile,func)
                particles=particleFileIO.readFITSFile('/global/cscratch1/sd/sukhdeep/mb2_subsample/',j,i)    
                corFunc.corFunc(particles,sightBins,rscale,nbins,min_sep,max_sep,rpar_min,rpar_max,min_box,max_box,func,logfile,savefile,combine,functype,key=key)

if galFiles==1:
    functype='GM'
    if combine==1:
        for j in snapArray:
            key='snap_'+str(j)
            savefile='/global/homes/h/hmarti21/Results/Current_Results/'+dir+'/DelGM_'+str(j)+'combine'
            #logfile='/global/homes/h/hmarti21/Results/Current_Results/'+dir+'/log_DelGM_'+str(j)+'combine'
            #corFunc.makeLog(logfile,func)
            particles=particleFileIO.combinedParticles(filePath,ptype,j)
            data=table.Table()
            data['x']=[]
            data['y']=[]
            data['z']=[]
            data['mass']=[]
            arr=particles
            fraction=0.1
            rand=np.random.uniform(0,1,size=len(arr['mass'][:]))
            rand_mask=rand<fraction
            arr_masked=arr[rand_mask]
            galaxies=galaxyFileIO.read_Data_MyHDF5(filename,col_names,cuts,key,logfile)
            corFunc.corFunc(galaxies,sightBins,rscale,nbins,min_sep,max_sep,rpar_min,rpar_max,min_box,max_box,func,logfile,savefile,combine,functype,key=key,fname2='yes',dat2=arr_masked)
    if combine==0:
        for j in ptype:
            for i in snapArray:
                key='snap_'+str(i)
                savefile='/global/homes/h/hmarti21/Results/Current_Results/'+dir+'/DelGM_'+str(i)+'ptype_'+str(j)
                #logfile='/global/homes/h/hmarti21/Results/Current_Results/'+dir+'/log_DelGM_'+str(i)+'ptype_'+str(j)
                #corFunc.makeLog(logfile,func)
                particles=particleFileIO.readFITSFile(filePath,j,i)
                galaxies=galaxyFileIO.read_Data_MyHDF5(filename,col_names,cuts,key,logfile)
                corFunc.corFunc(galaxies,sightBins,rscale,nbins,min_sep,max_sep,rpar_min,rpar_max,min_box,max_box,func,logfile,savefile,combine,functype,key=key,fname2='yes',dat2=particles)

if galFiles==2:
    #t1=time.time()
    #print(t1)
    combine=0
    functype='GG'
    for i in snapArray:
        key='snap_'+str(i)
        savefile='/global/homes/h/hmarti21/Results/Current_Results/'+dir+'/DelGG_'+str(i)
        #logfile='/global/homes/h/hmarti21/Results/Current_Results/'+dir+'/log_DelGG_'+str(i)
        #corFunc.makeLog(logfile,func)
        galaxies=galaxyFileIO.read_Data_MyHDF5(filename,col_names,cuts,key,logfile)
        corFunc.corFunc(galaxies,sightBins,rscale,nbins,min_sep,max_sep,rpar_min,rpar_max,min_box,max_box,func,logfile,savefile,combine,functype,key=key)
        #t2=time.time()
        #time=t2-t1
        #print(time)



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
