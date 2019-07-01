
# coding: utf-8

# In[1]:

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
filePath='/global/homes/h/hmarti21/mb2_subsample/'
snapArray=[85]
ptype=[0,1,4]
combine=np.int(sys.argv[6])
dir=sys.argv[7]
func=sys.argv[8]
galFiles=np.int(sys.argv[9])
logfile=sys.argv[10]
frac=np.float(sys.argv[11])

def subsample(frac,data):
    rand=np.random.uniform(0,1,size=len(data['Mass'][:]))
    rand_mask=rand<frac
    return data[rand_mask]

def kpctoMpc(data):
    data['x']=data['x']/1000
    data['y']=data['y']/1000
    data['z']=data['z']/1000
    return data

if galFiles==0:
    func='NN'
    functype='MM'
    if combine==1:
        for j in snapArray:
            key='snap_'+str(j)
            savefile='/global/homes/h/hmarti21/Results/Final_Results/'+dir+'/DelMM_'+str(j)+'combine'
            particles=particleFileIO.combinedParticles(filePath,ptype,j)
            particle_sample=subsample(frac,particles)
            particle_sample=kpctoMpc(particle_sample)
            corFunc.corFunc(particle_sample,sightBins,rscale,nbins,min_sep,max_sep,rpar_min,rpar_max,min_box,max_box,func,logfile,savefile,combine,functype,key=key)
    if combine==0:
        for j in ptype:
            for i in snapArray:
                key='snap_'+str(i)
                savefile='/global/homes/h/hmarti21/Results/Final_Results/'+dir+'/DelMM_'+str(i)+'ptype_'+str(j)
                particles=particleFileIO.readFITSFile('/global/cscratch1/sd/sukhdeep/mb2_subsample/',j,i)
                particles=kpctoMpc(particles)
                particles_sample=subsample(particles)
                corFunc.corFunc(particles_sample,sightBins,rscale,nbins,min_sep,max_sep,rpar_min,rpar_max,min_box,max_box,func,logfile,savefile,combine,functype,key=key)

if galFiles==1:
    functype='GM'
    if combine==1:
        for j in snapArray:
            key='snap_'+str(j)
            savefile='/global/homes/h/hmarti21/Results/Final_Results/'+dir+'/DelGM_'+str(j)+'combine'
            particles=particleFileIO.combinedParticles(filePath,ptype,j)
            particle_sample=subsample(frac,particles)
            particle_sample=kpctoMpc(particle_sample)
            galaxies=galaxyFileIO.read_Data_MyHDF5(filename,col_names,cuts,key,logfile)
            corFunc.corFunc(galaxies,sightBins,rscale,nbins,min_sep,max_sep,rpar_min,rpar_max,min_box,max_box,func,logfile,savefile,combine,functype,key=key,fname2='yes',data2=particle_sample)
    if combine==0:
        for j in ptype:
            for i in snapArray:
                key='snap_'+str(i)
                savefile='/global/homes/h/hmarti21/Results/Final_Results/'+dir+'/DelGM_'+str(i)+'ptype_'+str(j)
                particles=particleFileIO.readFITSFile(filePath,j,i)
                particles=kpctoMpc(particles)
                galaxies=galaxyFileIO.read_Data_MyHDF5(filename,col_names,cuts,key,logfile)
                corFunc.corFunc(galaxies,sightBins,rscale,nbins,min_sep,max_sep,rpar_min,rpar_max,min_box,max_box,func,logfile,savefile,combine,functype,key=key,fname2='yes',data2=particles)

if galFiles==2:
    combine=0
    functype='GG'
    for i in snapArray:
        key='snap_'+str(i)
        print('Testing pipeline')
        savefile='/global/homes/h/hmarti21/Results/Final_Results/'+dir+'/DelGG_'+str(i)
        galaxies=galaxyFileIO.read_Data_MyHDF5(filename,col_names,cuts,key,logfile)
        corFunc.corFunc(galaxies,sightBins,rscale,nbins,min_sep,max_sep,rpar_min,rpar_max,min_box,max_box,func,logfile,savefile,combine,functype,key=key)


