
# coding: utf-8

# In[2]:


import numpy as np
from astropy import table
from nbodykit.io.gadget import Gadget1File
from astropy.io import fits


# In[1]:


def particleFile(path,col_def,fraction,num_files,ptype):
    data=table.Table()
    data['x']=[]
    data['y']=[]
    data['z']=[]
    data['mass']=[]
    for i in np.arange(num_files):
        file=Gadget1File(path+str(i),columndefs=col_def,ptype=ptype)
        arr=file[:]
        rand=np.random.uniform(0,1,size=len(arr['Mass'][:]))
        rand_mask=rand<fraction
        arr_masked=arr[rand_mask]
        dat_table=table.Table()
        dat_table['x']=arr_masked['Position'][:,0]
        dat_table['y']=arr_masked['Position'][:,1]
        dat_table['z']=arr_masked['Position'][:,2]
        dat_table['mass']=arr_masked['Mass']
        data=table.vstack([data,dat_table])
    data.write('hmarti21_sampledData.fits',format='fits',overwrite=True)

def readFITSFile(path,ptype,snapArray):
    particleTable=table.Table()
    for i in snapArray:
        read=fits.open(path+'snapshot_'+'0'+str(i)+'_type'+str(ptype)+'.fits')
        particleTable=table.vstack([particleTable,table.Table(read[1].data)])
    return particleTable


# In[4]:


filepath='/global/cscratch1/sd/sukhdeep/snapdir_194/snapshot_194.'
col_def=[('Position', ('f8', 3), 'all'), ('Mass','auto',None)]
fraction=0.05
num_files=2
ptype=1
particleFile(filepath,col_def,fraction,num_files,ptype)


# In[22]:


i=0
file=Gadget1File(filepath+str(i),columndefs=col_def,ptype=ptype)
arr=file[['Position','Mass']][:]


# In[24]:


len(arr)


# In[5]:


read=fits.open('hmarti21_sampledData.fits')


# In[6]:


dat=read[1].data


# In[9]:


print(dat)
dat.dtype


# In[10]:


print(table.Table(dat))
table.Table(dat).dtype


# In[11]:


len(table.Table(dat)['z'])


# In[12]:


len(table.Table(dat)['x'])


# In[13]:


len(table.Table(dat)['y'])


# In[14]:


len(table.Table(dat)['mass'])

