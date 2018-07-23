
# coding: utf-8

# In[3]:

#export HDF5_USE_FILE_LOCKING=FALSE
import sys
sys.path.append('/global/homes/h/hmarti21/Pipeline')
import particleFileIO
import corFunc
import numpy as np
from matplotlib import pyplot as plt

import importlib
reload=importlib.reload
import corFunc
reload(corFunc)
import corFunc


# In[4]:


filename='/global/homes/h/hmarti21/data_sukhdeep.hdf5'
col_names={'x':'halos.x','y':'halos.y','z':'halos.z','mass':'halos.mass','q':'shapesStar.q2d','ax':'shapesStar.a2d_x','ay':'shapesStar.a2d_y','bx':'shapesStar.b2d_x','by':'shapesStar.b2d_y'}    #,'e1':None,'e2':None}
cuts='/global/homes/h/hmarti21/cuts.txt'
sightBins=3
rscale=5
nbins=10
min_sep=1000
max_sep=10000
rpar_step=10000
min_box=0
max_box=100000
key='hydro_full'
logfile='/global/homes/h/hmarti21/log.txt'
savefile='/global/homes/h/hmarti21/data.txt'

filepath='/global/cscratch1/sd/sukhdeep/snapdir_194/snapshot_194.'
col_def=[('Position', ('f8', 3), 'all'), ('Mass','auto',None)]
fraction=0.05
num_files=1
ptype=1
particleFileIO.particleFile(filepath,col_def,fraction,num_files,ptype)
particles=particleFileIO.readFITSFile('/global/homes/h/hmarti21/hmarti21_sampledData.fits')
NNPairs=np.zeros((sightBins,nbins))
NRPairs=np.zeros((sightBins,nbins))
RRPairs=np.zeros((sightBins,nbins))
RNPairs=np.zeros((sightBins,nbins))

print(filename,col_names,cuts,sightBins,rscale,nbins,min_sep,max_sep,rpar_step,min_box,max_box,logfile,savefile,key)
NNPairs,NRPairs,RNPairs,RRPairs=corFunc.corFunc(filename,col_names,cuts,sightBins,rscale,nbins,min_sep,max_sep,rpar_step,min_box,max_box,'NN',logfile,savefile,key=key) #fname2=particles)

??corFunc.corFunc
# In[5]:


import numpy as np
data=np.loadtxt('/global/homes/h/hmarti21/data.txt')
print(data)


# In[38]:


from matplotlib import pyplot as plt
plt.plot(np.logspace(0,1,10),data[0])
plt.xscale('log')
plt.yscale('log')


# In[ ]:


plt.plot(np.logspace(0,1,10),data[0])
plt.xscale('log')


# In[ ]:


plt.plot(np.logspace(0,1,10),1/np.logspace(0,1,10))
plt.xscale('log')
plt.yscale('log')


# In[ ]:


plt.plot(np.arange(10),np.logspace(0,1,10))
plt.xscale('log')


# In[5]:


NNPairs


# In[27]:


for i in np.arange(sightBins):
    plt.plot(np.arange(nbins),NNPairs[i],label='n'+str(i))
    plt.plot(np.arange(nbins),RRPairs[i],label='r'+str(i))
    plt.legend()


# In[31]:


plt.plot(np.arange(nbins),NNPairs[0],label='n')
plt.plot(np.arange(nbins),RRPairs[0],label='r')
plt.legend()


# In[32]:


plt.plot(np.logspace(0,1,10),NNPairs[0],label='n')
plt.plot(np.logspace(0,1,10),RRPairs[0],label='r')
plt.plot(np.logspace(0,1,10),NNPairs[1],label='n')
plt.plot(np.logspace(0,1,10),RRPairs[1],label='r')
plt.legend()
plt.xscale('log')
plt.yscale('log')


# In[15]:


rpar_step/1000


# In[19]:


SD=np.sum(NNPairs/RRPairs*rscale**2-1,axis=0)*rpar_step/1000
SD0=(NNPairs[0]/RRPairs[0]*rscale**2-1)*rpar_step/1000
SR=NRPairs[0]/RRPairs[0]*rscale
DR=SR #true only for auto correlation
plt.plot(np.logspace(0,1,10),SD,'b-')
plt.plot(np.logspace(0,1,10),SD0,'b--')
#plt.plot(np.logspace(0,1,10),SR-1,'c-')
#plt.plot(np.logspace(0,1,10),SD-SR-DR+1,'r--')
#plt.plot(np.logspace(0,1,10),data[0],'g--')
plt.xscale('log')
plt.yscale('log')


# In[7]:


rb0=np.logspace(0,1,11)
rb=0.5*(rb0[1:]+rb0[:-1])
dr=rb0[1:]-rb0[:-1]


# In[8]:


get_ipython().run_line_magic('pylab', 'inline')


# In[9]:


plt.plot(rb,NRPairs[0],label='n')
plt.plot(rb,RRPairs[0]/rscale,'r--',label='r')
plt.plot(rb,rb*dr*1.5e5)
# plt.plot(np.logspace(0,1,10),NNPairs[1],label='n')
# plt.plot(np.logspace(0,1,10),RRPairs[1],label='r')
# plt.plot(np.logspace(0,1,10),NNPairs[2])

plt.legend()
plt.xscale('log')
plt.yscale('log')


# In[11]:


get_ipython().run_line_magic('pinfo2', 'corFunc.NN')


# In[22]:


np.sum(NNPairs,axis=0)

