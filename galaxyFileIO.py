import h5py
import numpy as np
import treecorr
from astropy.table import Table

def read_Data_hdf5(filename,col_dict,cuts,key,logfile):
    full_dat=Table.read(filename,path=key)
    dat=Table()
    log=open(logfile,'r').read()
    logwrite=open(logfile,'w')
    logwrite.write(log)
    for i,v in col_dict.items():
        if v!=None:
            dat[i]=full_dat[v]
    cuts=open(cuts)
    for line in cuts:
        dat=dat[eval(line)]
        logwrite.write(line)
    dat['x']=dat['x']/1000
    dat['y']=dat['y']/1000
    dat['z']=dat['z']/1000
    logwrite.write(str(len(dat)))
    return dat

def read_Data_MyHDF5(filename,col_dict,cuts,key,logfile):
    full_dat=h5py.File('/global/homes/h/hmarti21/dataFrameTest','r')
    dat=Table()
    log=open(logfile,'r').read()
    logwrite=open(logfile,'w')
    logwrite.write(log)
    for i,v, in col_dict.items():
        if v!=None:
            dat[i]=full_dat[key]['block0_values'][:,v]
    dat['cent']=full_dat[key]['block1_values'][:,1]
    cuts=open(cuts)
    for line in cuts:
        dat=dat[eval(line)]
        logwrite.write(line+'\n')
    dat['x']=dat['x']/1000
    dat['y']=dat['y']/1000
    dat['z']=dat['z']/1000
    #print (dat)
    return dat
