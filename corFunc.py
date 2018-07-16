
# coding: utf-8

# In[37]:

import numpy as np
from astropy.table import Table
import treecorr
import h5py
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

def noPadding(table, primary_key):
    return table[table[primary_key]!=0]

def projShapeCal(table):
    ax=table['ax']
    ay=table['ay']
    q=table['q']
    phi=np.arctan2(ay,ax)
    e1=np.divide(1-q,1+q)*np.cos(2*phi)
    e2=np.divide(1-q,1+q)*np.sin(2*phi)
    table['e1']=e1
    table['e2']=e2
    return table

def NG(table,rand_pcounts,sightBins,nbins,min_sep,max_sep,rpar_step,fname,tableN=None):
    table=noPadding(table,'ax')
    table=projShapeCal(table)
    NGCor=list()
    if tableN==None:
        catG=treecorr.Catalog(x=table['x'],y=table['y'],z=table['z'],g1=table['e1'],g2=table['e2'])
        for i in np.arange(sightBins):
            NGCor.append(treecorr.NGCorrelation(nbins=nbins,min_sep=min_sep,max_sep=max_sep,metric="Rperp",
                                                min_rpar=i*rpar_step,max_rpar=(i+1)*rpar_step))
            NGCor[i].process(catG,catG,metric='Rperp')
    else:
        catG=treecorr.Catalog(g1=table['e1'],g2=table['e2'])
        catN=treecorr.Catalog(x=tableN['x'],y=tableN['y'],z=tableN['z'])
        for i in np.arange(sightBins):
            NGCor.append(treecorr.NGCorrelation(nbins=nbins,min_sep=min_sep,max_sep=max_sep,metric="Rperp",min_rpar=i*rpar_step,max_rpar=(i+1)*rpar_step))
            NGCor[i].process(catN,catG,metric='Rperp')
    NG_ksi=np.zeros((sightBins,nbins))
    for i in np.arange(sightBins):
        NG_ksi[i]=div0(NGCor[i].npairs,rand_pcounts[i])
    NG_w=np.sum(NG_ksi,axis=0)
    np.savetxt(fname,NG_w)

def NN(table,rand_cat,sightBins,nbins,min_sep,max_sep,rpar_step,fname,table2=None):
    catN=treecorr.Catalog(x=table['x'],y=table['y'],z=table['z'])
    if table2!=None:
        catN2=treecorr.Catalog(x=table2['x'],y=table2['y'],z=table2['z'])
    NNCor=list()
    NRCor=list()
    RRCor=list()
    RNCor=list()
    NNPairs=np.zeros((sightBins,nbins))
    NRPairs=np.zeros((sightBins,nbins))
    RNPairs=np.zeros((sightBins,nbins))
    RRPairs=np.zeros((sightBins,nbins))
    for i in np.arange(sightBins):
        NNCor.append(treecorr.NNCorrelation(nbins=nbins,min_sep=min_sep,max_sep=max_sep,metric="Rperp",
                                               min_rpar=i*rpar_step,max_rpar=(i+1)*rpar_step))
        NRCor.append(treecorr.NNCorrelation(nbins=nbins,min_sep=min_sep,max_sep=max_sep,metric="Rperp",
                                               min_rpar=i*rpar_step,max_rpar=(i+1)*rpar_step))
        RRCor.append(treecorr.NNCorrelation(nbins=nbins,min_sep=min_sep,max_sep=max_sep,metric="Rperp",
                                               min_rpar=i*rpar_step,max_rpar=(i+1)*rpar_step))
        if table2==None:
            NNCor[i].process(catN,metric='Rperp')
            NRCor[i].process(catN,rand_cat,metric='Rperp')
            RRCor[i].process(rand_cat,metric='Rperp')
            NNPairs[i]=NNCor[i].npairs
            NRPairs[i]=NRCor[i].npairs
            RRPairs[i]=RRCor[i].npairs
            print(NNCor[i].npairs)
            print(NRCor[i].npairs)
            print(RRCor[i].npairs)
        else:
            RNCor.append(treecorr.NNCorrelation(nbins=nbins,min_sep=min_sep,max_sep=max_sep,metric='Rperp',
                                              min_rpar=i*rpar_step,max_rpar=(i+1)*rpar_step))
            NNCor[i].process(catN,catN2,metric='Rperp')
            NRCor[i].process(catN,rand_cat,metric='Rperp')
            RNCor[i].process(catN2,rand_cat,metric='Rperp')
            RRCor[i].process(rand_cat,rand_cat,metric='Rperp')
            NNPairs[i]=NNCor[i].npairs
            NRPairs[i]=NRCor[i].npairs
            RNPairs[i]=RNCor[i].npairs
            RRPairs[i]=RRCor[i].npairs
    xi=np.zeros((nbins,sightBins))
    varxi=np.zeros((nbins,sightBins))
    for i in np.arange(sightBins):
        if table2==None:
            xi[:,i],varxi[:,i]=NNCor[i].calculateXi(RRCor[i],dr=NRCor[i])
        else:
            xi[:,i],varxi[:,i]=NNCor[i].calculateXi(RRCor[i],dr=NRCor[i],rd=RNCor[i])
    xi_tot=np.sum(xi,axis=1)
    varxi_tot=np.sum(varxi,axis=1)
    np.savetxt(fname, np.array([xi_tot,varxi_tot]))
    print(NNPairs,NRPairs,RNPairs,RRPairs)
    return NNPairs,NRPairs,RNPairs,RRPairs

def RR(table,scale,min_box,max_box,sightBins,nbins,min_sep,max_sep,rpar_step,cat_out=False):
    rand=np.random.uniform(min_box,max_box,size=len(table['x'])*scale*3)
    rand=np.random.permutation(rand).reshape(3,len(table['x'])*scale)
    print(min_box,max_box,np.amin(rand),np.amax(rand))
    catR=treecorr.Catalog(x=rand[0],y=rand[1],z=rand[2])
    if cat_out==True:
        return catR
    RRCor=list()
    rand_pcounts=np.zeros((sightBins,nbins))
    for i in np.arange(sightBins):
        RRCor.append(treecorr.NNCorrelation(nbins=nbins,min_sep=min_sep,max_sep=max_sep,metric="Rperp",
                                            min_rpar=i*rpar_step,max_rpar=(i+1)*rpar_step))
        RRCor[i].process(catR,catR,metric='Rperp')
        rand_pcounts[i]=RRCor[i].npairs
    return rand_pcounts

def corFunc(fname1,col_name,cuts,sightBins,rscale,nbins,min_sep,max_sep,rpar_step,min_box,max_box,func,logfile,savefile,
            fname2=None,col_name2=None,key=None,key2=None,cuts2=None):
    reset=open(logfile,'w')
    reset.write('------BEGIN------')
    reset.write('------FUNCTION------')
    reset.write(func)
    reset.write('------1ST FILE------')
    reset.close()
    data=read_Data_hdf5(fname1,col_name,cuts,key,logfile)
    log_file=open(logfile,'r').read()
    log=open(logfile,'w')
    log.write(log_file)
    log.write(fname1)
    NNPairs=np.zeros((sightBins,nbins))
    NRPairs=np.zeros((sightBins,nbins))
    RRPairs=np.zeros((sightBins,nbins))
    RNPairs=np.zeros((sightBins,nbins))

    if func=='NN' and fname2==None:
        catR=RR(data,rscale,min_box,max_box,sightBins,nbins,min_sep,max_sep,rpar_step,cat_out=True)
        NNPairs,NRPairs,RNPairs,RRPairs=NN(data,catR,sightBins,nbins,min_sep,max_sep,rpar_step,savefile)
        log.close()
        return NNPairs,NRPairs,RNPairs,RRPairs
    if func=='NG' and fname2==None:
        data=noPadding(data,'ax')
        rand_pairs=RR(data,rscale,min_box,max_box,sightBins,nbins,min_sep,max_sep,rpar_step)
        NG(data,rand_pairs,sightBins,nbins,min_sep,max_sep,rpar_step,savefile)
        log.close()
    if fname2!=None:
        log.write('------2ND FILE------')
        data2=fname2 #read_Data_hdf5(fname2,col_name2,cuts2,key2)
        if func=='NN':
            print(data2)
            catR=RR(data,rscale,min_box,max_box,sightBins,nbins,min_sep,max_sep,rpar_step,cat_out=True)
            NN(data,catR,sightBins,nbins,min_sep,max_sep,rpar_step,savefile,table2=data2)
            log.close()
        if func=='NG':
            data=noPadding(data,'ax')
            rand_pairs=RR(data,rscale,min_box,max_box,sightBins,nbins,min_sep,max_sep,rpar_step)
            NG(data,rand_pairs,sightBins,nbins,min_sep,max_sep,rpar_step,savefile,tableN=data2)
            log.close()

