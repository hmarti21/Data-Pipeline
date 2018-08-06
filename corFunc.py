
# coding: utf-8

# In[37]:

import numpy as np
from astropy.table import Table
import treecorr
import h5py
import matplotlib.pyplot as plt
from astropy.io import ascii

def div0( a, b ):
    """ ignore / 0, div0( [-1, 0, 1], 0 ) -> [0, 0, 0] """
    with np.errstate(divide='ignore', invalid='ignore'):
        c = np.true_divide( a, b )
        c[ ~ np.isfinite( c )] = 0  # -inf inf NaN
    return c

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
    cuts=open(cuts)
    for line in cuts:
        dat=dat[eval(line)]
        logwrite.write(line)
    logwrite.write(str(len(dat)))
    return dat
def toMPC(table):
    table['x']=table['x']/1000
    table['y']=table['y']/1000
    table['z']=table['z']/1000
    print(table)
    return table

def noPadding(table, primary_key):
    return table[table[primary_key]!=0]

def rParStepCal(rpar_min,rpar_max,sightBins):
    return (rpar_max-rpar_min)/sightBins

def projShapeCal(table):
    ax=table['ax']
    ay=table['ay']
    q=table['q']
    phi=np.arctan2(ay,ax)
    e1=np.divide(1-q,1+q)*np.cos(2*phi)
    e2=np.divide(1-q,1+q)*np.sin(2*phi)
    table['e1']=e1
    table['e2']=e2
    print(table['e1'],table['e2'])
    return table

def CorProcess(catN,catG,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,
		which_corr='NG'):
    corList=list()
    pairCounts=np.zeros((sightBins,nbins))
    if which_corr=='NG':
        for i in np.arange(sightBins):
            corList.append(treecorr.NGCorrelation(nbins=nbins,min_sep=min_sep,max_sep=max_sep,metric='Rperp',
                                             min_rpar=i*rpar_step+rpar_min,max_rpar=(i+1)*rpar_step+rpar_min))
            corList[i].process(catN,catG,metric='Rperp')
            pairCounts[i]=corList[i].npairs
        return corList,pairCounts
    if which_corr=='NN':
        for i in np.arange(sightBins):
            corList.append(treecorr.NNCorrelation(nbins=nbins,min_sep=min_sep,max_sep=max_sep,min_rpar=i*rpar_step+rpar_min,max_rpar=(i+1)*rpar_step+rpar_min))
            corList[i].process(catN,catG,metric='Rperp')
        return corList

def NG(table,rscale,rand_pcounts,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,fname,extraTable=False,tableN=None):
    table=projShapeCal(table)
    if extraTable==False:
        catG=treecorr.Catalog(x=table['x'],y=table['y'],z=table['z'],g1=table['e1'],g2=table['e2'])
        NGCor,NGCorCounts=CorProcess(catG,catG,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min)
    else:
        catG=treecorr.Catalog(x=table['x'],y=table['y'],z=table['z'],g1=table['e1'],g2=table['e2'])
        catN=treecorr.Catalog(x=tableN['x'],y=tableN['y'],z=tableN['z'])
        NGCor,NGCorCounts=CorProcess(catN,catG,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min)
    NG_ksi=np.zeros((sightBins,nbins))
    if extraTable==False:
        for i in np.arange(sightBins):
            NG_ksi[i]=div0(NGCorCounts[i],rand_pcounts[i])
    else:
        for i in np.arange(sightBins):
            NG_ksi[i]=div0(NGCorCounts[i],rand_pcounts[i])/(len(tableN['x'])/len(table['x']))
    NG_w=np.sum(NG_ksi,axis=0)
    np.savetxt(fname,NG_w)

def NN(table,rand_cat,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,fname,extraTable=False,table2=None):
    catN=treecorr.Catalog(x=table['x'],y=table['y'],z=table['z'])
    #NNCounts=np.zeros((sightBins,nbins))
    #NRCounts=np.zeros((sightBins,nbins))
    #RNCounts=np.zeros((sightBins,nbins))
    #RRCounts=np.zeros((sightBins,nbins))
    if extraTable==True:
        catN2=treecorr.Catalog(x=table2['x'],y=table2['y'],z=table2['z'])
    if extraTable==False:
        NNCor=CorProcess(catN,catN,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,which_corr='NN')
        NRCor=CorProcess(catN,rand_cat,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,which_corr='NN')
        RRCor=CorProcess(rand_cat,rand_cat,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,which_corr='NN')
    else:
        NNCor=CorProcess(catN,catN2,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,which_corr='NN')
        NRCor=CorProcess(catN,rand_cat,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,which_corr='NN')
        RNCor=CorProcess(catN2,rand_cat,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,which_corr='NN')
        RRCor=CorProcess(rand_cat,rand_cat,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,which_corr='NN')
    xi=np.zeros((nbins,sightBins))
    varxi=np.zeros((nbins,sightBins))
    for i in np.arange(sightBins):
        if extraTable==False:
            xi[:,i],varxi[:,i]=NNCor[i].calculateXi(RRCor[i],dr=NRCor[i])
        else:
            xi[:,i],varxi[:,i]=NNCor[i].calculateXi(RRCor[i],dr=NRCor[i],rd=RNCor[i])
    xi_tot=np.sum(xi,axis=1)*rpar_step*2
    varxi_tot=np.sum(varxi,axis=1)
    np.savetxt(fname, np.array([xi_tot,varxi_tot]))
    #return NNCounts,NRCounts,RNCounts,RRCounts

def RR(table,scale,min_box,max_box,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,cat_out=False):
    rand=np.random.uniform(min_box,max_box,size=len(table['x'])*scale*3)
    rand=np.random.permutation(rand).reshape(3,len(table['x'])*scale)
    catR=treecorr.Catalog(x=rand[0],y=rand[1],z=rand[2])
    if cat_out==True:
        return catR
    RRCor=list()
    rand_pcounts=np.zeros((sightBins,nbins))
    RRCor=CorProcess(catR,catR,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,which_corr='NN')
    for i in np.arange(sightBins):
        rand_pcounts[i]=RRCor[i].npairs
    return rand_pcounts

def corFunc(fname1,col_name,cuts,sightBins,rscale,nbins,min_sep,max_sep,rpar_min,rpar_max,min_box,max_box,func,logfile,savefile,fname2=None,col_name2=None,key=None,key2=None,cuts2=None,dat2=None):
    reset=open(logfile,'w')
    reset.write('------BEGIN------')
    reset.write('------FUNCTION------')
    reset.write(func)
    reset.write('------1ST FILE------')
    reset.close()
    data=read_Data_MyHDF5(fname1,col_name,cuts,key,logfile)
    data=toMPC(data)
    rpar_step=rParStepCal(rpar_min,rpar_max,sightBins)
    log_file=open(logfile,'r').read()
    log=open(logfile,'w')
    log.write(log_file)
    log.write(fname1)
    log.write('LoS Bin Step:'+str(rpar_step))
    #NNCounts=np.zeros((sightBins,nbins))
    #NRCounts=np.zeros((sightBins,nbins))
    #RNCounts=np.zeros((sightBins,nbins))
    #RRCounts=np.zeros((sightBins,nbins))
    if func=='NN' and fname2==None:
        catR=RR(data,rscale,min_box,max_box,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,cat_out=True)
        #NNCounts,NRCounts,RNCounts,RRCounts=
        NN(data,catR,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,savefile)
        log.close()
        #return NNCounts,NRCounts,RNCounts,RRCounts
    if func=='NG' and fname2==None:
        data=noPadding(data,'ax')
        rand_pairs=RR(data,rscale,min_box,max_box,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min)
        NG(data,rscale,rand_pairs,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,savefile)
        log.close()
    if fname2!=None:
        log.write('------2ND FILE------')
        data2=dat2 #read_Data_hdf5(fname2,col_name2,cuts2,key2)
        data2=toMPC(data2)
        if func=='NN':
            catR=RR(data,rscale,min_box,max_box,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,cat_out=True)
            NN(data,catR,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,savefile,extraTable=True,table2=data2)
            log.close()
            
        if func=='NG':
            data=noPadding(data,'ax')
            rand_pairs=RR(data,rscale,min_box,max_box,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min)
            NG(data,rscale,rand_pairs,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,savefile,extraTable=True,tableN=data2)
            log.close()


# In[ ]:



