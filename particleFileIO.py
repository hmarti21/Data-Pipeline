
# coding: utf-8

# In[31]:

import numpy as np
from astropy import table
from nbodykit.io.gadget import Gadget1File
from astropy.io import fits


# In[45]:

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

def readFITSFile(path):
    read=fits.open(path)
    return read[1].data

