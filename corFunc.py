
# coding: utf-8

# In[37]:

#Illustrous TNG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

import time
import numpy as np
from astropy.table import Table
import treecorr
import h5py
import matplotlib.pyplot as plt
from astropy.io import ascii

def div0( a, b ):
    """ ignore / 0, div0( [-1, 0, 1], 0 ) -> [0, 0, 0] """
    with np.errstate(divide='ignore', invalid='ignore'):
        #print(a,b)
        c = np.true_divide( a, b )
        #print(np.isfinite(c))
        #print(c)
        c[ ~ np.isfinite( c )] = 0  # -inf inf NaN
        #print(c)
    return c

def noPadding(table, primary_key):
    return table[table[primary_key]!=0]

def rParStepCal(rpar_min,rpar_max,sightBins):
    return (rpar_max-rpar_min)/sightBins

def projShapeCal(table,functype):
    if functype=='MM':
        return table
    else:
        ax=table['ax']
        ay=table['ay']
        q=table['q']
        phi=np.arctan2(ay,ax)
        e1=np.divide(1-q,1+q)*np.cos(2*phi)
        e2=np.divide(1-q,1+q)*np.sin(2*phi)
        table['e1']=e1
        table['e2']=e2
        return table

def CorProcess(catN,catG,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,logfile,which_corr='NG',RR=False):
    #print('test corProcess',file=logfile)
    corList=list()
    pairCounts=np.zeros((sightBins,nbins))
    wpairCounts=np.zeros((sightBins,nbins))
    for i in np.arange(sightBins):
        if which_corr=='NG' and not RR:
            corList.append(treecorr.NGCorrelation(nbins=nbins,min_sep=min_sep,max_sep=max_sep,min_rpar=i*rpar_step+rpar_min,max_rpar=(i+1)*rpar_step+rpar_min)) #metric='Rperp'
        if which_corr=='NN' or RR:
            #print('RR test', file=logfile)
            corList.append(treecorr.NNCorrelation(nbins=nbins,min_sep=min_sep,max_sep=max_sep,min_rpar=i*rpar_step+rpar_min,max_rpar=(i+1)*rpar_step+rpar_min)) #metric='Rperp'
        corList[i].process(catN,catG,metric='Rperp')
        pairCounts[i]=corList[i].npairs
        wpairCounts[i]=corList[i].weight
    return corList,pairCounts,wpairCounts

# def NG(table,rscale,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,min_box,max_box,fname,combine,extraTable=False,tableN=None):
#     table=projShapeCal(table)

#     rand=np.random.rand(3,len(table['x'])*rscale)*max_box+min_box
#     catR=treecorr.Catalog(x=rand[0],y=rand[1],z=rand[2])
#     RRCor=list()
#     rand_pcounts=np.zeros((sightBins,nbins))
#     RRCor=CorProcess(catR,catR,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,which_corr='NN')
#     for i in np.arange(sightBins):
#         rand_pcounts[i]=RRCor[i].npairs
# #get random pair counts here, to use in the formula, then need catR to build NRcor counts

#     catG=treecorr.Catalog(x=table['x'],y=table['y'],z=table['z'],g1=table['e1'],g2=table['e2'])
#     NRCor,NRCorCounts=CorProcess(catR,catG,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min)
#     if extraTable==False:
#         #galaxy position-shape autocorrelation
#         NGCor,NGCorCounts=CorProcess(catG,catG,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min)
#     else:
#         if combine==0:
#             #particle/galaxy cross correlation, per particle type
#             catN=treecorr.Catalog(x=tableN['x'],y=tableN['y'],z=tableN['z'])
#         else:
#             #particle/galaxy cross correlation, for all particles at once
#             catN=treecorr.Catalog(x=tableN['x'],y=tableN['y'],z=tableN['z'],w=tableN['Mass'])
#         NGCor,NGCorCounts=CorProcess(catN,catG,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min)
#     NG_ksi=np.zeros((sightBins,nbins))
#     if extraTable==False:
#         for i in np.arange(sightBins):
#             NG_ksi[i]=div0(NGCorCounts[i]/(len(table['x'])**2)-NRCorCounts[i]/(len(table['x'])*len(rand[0])),rand_pcounts[i]/len(rand[0])**2)
#     else:
#         for i in np.arange(sightBins):
#             NG_ksi[i]=div0(NGCorCounts[i]/(len(table['x'])*catN.sumw)-NRCorCounts[i]/(len(table['x'])*len(rand[0])),rand_pcounts[i]/len(rand[0])**2)
#             #print('NG:',NGCorCounts[i],'weight:',len(table['x'])*catN.sumw,'NR:',NRCorCounts[i],'weight:', len(rand[0])*len(table['x']),'RR:',rand_pcounts[i],'weight:',len(rand[0])**2)
#     NG_w=np.sum(NG_ksi,axis=0)
#     np.savetxt(fname+'.txt',NG_w)

def getCatR(table,rscale,min_box,max_box):
    rand=np.random.uniform(min_box,max_box,size=(3,len(table['x'])*rscale))
    catR=treecorr.Catalog(x=rand[0],y=rand[1],z=rand[2])
    return catR, rand

def genCat(table,rscale,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,min_box,max_box,fname,logfile,combine,which_corr,functype,table2=None,
	wt1_col=None,wt2_col=None):
    catR,rand=getCatR(table,rscale,min_box,max_box)
    #print('test genCat',file=logfile)
    """if functype=='MM':
        if wt1_col is None:
            wt1_col='Mass'
        if wt2_col is None:
            wt2_col='Mass
    
    wt1=None
    wt2=None
    if wt1_col is not None:
         wt1=table[wt1_col]
    if wt2_col is not None:
        if table2 is None:
            wt2=table[wt2_col]
        else:
            wt2=table2[wt2_col]
    catN=treecorr.Catalog(x=table['x'],y=table['y'],z=table['z'],w=wt1)"""
    
    table=projShapeCal(table,functype)
    if functype=='MM':
        catN=treecorr.Catalog(x=table['x'],y=table['y'],z=table['z'],w=table['Mass'])
        catN2=catN
        #print('test MM',file=logfile)
        print('NNWeight:'+str(np.sum(table['Mass'])**2))
        print('NRWeight:'+str(np.sum(table['Mass'])*len(rand[0])))
        print('RRWeight:'+str(len(rand[0])**2))
    if functype=='GM':
        catN=treecorr.Catalog(x=table['x'],y=table['y'],z=table['z'],g1=table['e1'],g2=table['e2'])
        catN2=treecorr.Catalog(x=table2['x'],y=table2['y'],z=table2['z'],w=table2['Mass'])
        #print('test GM',file=logfile)
        print('NNWeight:'+str(np.sum(table2['Mass'])*len(table['x'])))
        print('NRWeight:'+str(len(table['x'])*len(rand[0])))
        print('RRWeight:'+str(len(rand[0])**2))
        print('RNWeight:'+str(len(rand[0])*np.sum(table2['Mass'])))
    if functype=='GG':
        catN=treecorr.Catalog(x=table['x'],y=table['y'],z=table['z'],g1=table['e1'],g2=table['e2'])
        catN2=catN
        #print('test GG',file=logfile)
        #print(catN.g1,catN.g2,file=logfile)
        #print(catN2.g1,catN2.g2,file=logfile)
        print('NNWeight:'+str(len(table['x'])**2))
        print('NRWeight:'+str(len(table['x'])*len(rand[0])))
        print('RRWeight:'+str(len(rand[0])**2))
    return catN,catN2,catR

def pairCounts(cat1,cat2,catR,table,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,fname,logfile,combine,which_corr,functype,table2=None):
    #catalogs go N,N2 for NN and N,G for NG for cat1, cat2
#     if False:
#         if which_corr=='NN':
#             NNCor,NNCorCounts,wNNCorCounts=CorProcess(cat1,cat2,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,which_corr)
#             NRCor,NRCorCounts,wNRCorCounts=CorProcess(cat1,catR,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,which_corr)
#             RRCor,RRCorCounts,wRRCorCounts=CorProcess(catR,catR,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,which_corr)
#             RNCor,RNCorCounts,wRNCorCounts=CorProcess(cat2,catR,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,which_corr)
#             weightRNCounts=np.zeros((sightBins,nbins))
#             NN_xi=np.zeros((sightBins,nbins))
#             for i in np.arange(sightBins):
#                 weightNRCounts=wNRCorCounts[i]/(cat1.sumw*catR.sumw)
#                 weightRRCounts=wRRCorCounts[i]/(catR.sumw**2)
#                 if functype=='GG':
#                     weightNNCounts=wNNCorCounts[i]/(cat1.sumw**2)
#                     NN_xi[i]=div0(weightNNCounts-2*weightNRCounts+weightRRCounts,weightRRCounts)
#                 if functype == 'GM':
#                     weightNNCounts=wNNCorCounts[i]/(cat1.sumw*cat2.sumw)
#                     weightRNCounts=wRNCorCounts[i]/(cat2.sumw*catR.sumw)
#                     NN_xi[i]=div0(weightNNCounts-weightNRCounts-weightRNCounts+weightRRCounts,weightRRCounts)
#                 if functype == 'MM':
#                     weightNNCounts=wNNCorCounts[i]/(cat1.sumw**2)
#                     NN_xi[i]=div0(weightNNCounts-2*weightNRCounts+weightRRCounts, weightRRCounts)
#             NN_w=np.sum(NN_ksi,axis=0)
#             np.savetxt(fname+'bin.txt',NNCor[0].logr)
#             np.savetxt(fname+'NNCount.txt',NNCorCounts)
#             np.savetxt(fname+'NRCount.txt',NRCorCounts)
#             np.savetxt(fname+'RRCount.txt',RRCorCounts)
#             np.savetxt(fname+'RNCount.txt',RNCorCounts)
#             np.savetxt(fname+'normNNCount.txt',weightNNCounts)
#             np.savetxt(fname+'normNRCount.txt',weightNRCounts)
#             np.savetxt(fname+'normRRCount.txt',weightRRCounts)
#             np.savetxt(fname+'normRNCount.txt',weightRNCounts)
#             np.savetxt(fname+'Result.txt',NN_w)
#             np.savetxt(fname+'mass.txt',table['Mass'])

#         if which_corr=='NG':
#             NGCor,NGCorCounts,wNGCorCounts=CorProcess(cat1,cat2,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,which_corr)
#             NRCor,NRCorCounts,wNRCorCounts=CorProcess(cat1,catR,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,which_corr)
#             RRCor,RRCorCounts,wRRCorCounts=CorProcess(catR,catR,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,which_corr)
#             NG_xi=np.zeros((sightBins,nbins))
#             for i in np.arange(sightBins):
#                 weightNGCounts=wNGCorCounts/(cat1.sumw*cat2.sumw)
#                 weightNRCounts=wNRCorCounts/(cat2.sumw*catR.sumw)
#                 weightRRCounts=wRRCorCounts/(catR.sumw**2)
#                 NG_xi[i]=div0(weightNGCounts-weightNRCounts,weightRRCounts)
#             NG_w=np.sum(NG_xi,axis=0)
#             np.savetxt(fname+'Result.txt',NG_w)
    
#     else:
            #print('test pairCounts', file=logfile)
            NNCor,NNCorCounts,wNNCorCounts=CorProcess(cat1,cat2,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,logfile,which_corr)
            NRCor,NRCorCounts,wNRCorCounts=CorProcess(catR,cat2,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,logfile,which_corr)
            RRCor,RRCorCounts,wRRCorCounts=CorProcess(catR,catR,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,logfile,which_corr,RR=True)
            if which_corr=='NN' and functype=='GM':
                RNCor,RNCorCounts,wRNCorCounts=CorProcess(catR,cat1,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,which_corr)
            weightRNCounts=np.zeros((sightBins,nbins))
            weightNRCounts=np.zeros((sightBins,nbins))
            weightRRCounts=np.zeros((sightBins,nbins))
            weightNNCounts=np.zeros((sightBins,nbins))
            NN_xi=np.zeros((sightBins,nbins))
            #xi=np.zeros((nbins,sightBins))
            #var_xi=np.zeros((nbins,sightBins))
            for i in np.arange(sightBins):
                weightNRCounts[i]=wNRCorCounts[i]/(cat2.sumw*catR.sumw)
                print('weight sum N1: '+str(cat1.sumw),'weight sum R: '+str(catR.sumw), 'weight sum N2: '+str(cat2.sumw))
                weightRRCounts[i]=wRRCorCounts[i]/(catR.sumw**2)
                weightNNCounts[i]=wNNCorCounts[i]/(cat1.sumw*cat2.sumw)
                if which_corr=='NN' and functype=='GM':
                    weightRNCounts[i]=wRNCorCounts[i]/(cat1.sumw*catR.sumw)
                    NN_xi[i]=div0(weightNNCounts[i]-weightNRCounts[i]-weightRNCounts[i]+weightRRCounts[i],weightRRCounts[i])
                elif which_corr=='NG':
                    NN_xi[i]=div0(weightNNCounts[i]-weightNRCounts[i],weightRRCounts[i])
                #elif which_corr=='NN' and functype=='GG':
                    #xi[:,i],var_xi[:,i]=NNCor[i].calculateXi(RRCor[i],dr=NRCor[i])
                    #xi_tot=np.sum(xi,axis=1)*rpar_step*2
                    #varxi_tot=np.sum(var_xi,axis=1)
                    #np.savetxt(fname+'w_test',np.array([xi_tot,varxi_tot]))
                else:
                    NN_xi[i]=div0(weightNNCounts[i]-2*weightNRCounts[i]+weightRRCounts[i],weightRRCounts[i])
            NN_w=np.sum(NN_xi,axis=0)*rpar_step*2
            np.savetxt(fname+'bin.txt',NNCor[0].logr)
            np.savetxt(fname+'NNCount.txt',NNCorCounts)
            np.savetxt(fname+'NRCount.txt',NRCorCounts)
            np.savetxt(fname+'RRCount.txt',RRCorCounts)
            np.savetxt(fname+'weightNNCount.txt',wNNCorCounts)
            np.savetxt(fname+'weightNRCount.txt',wNRCorCounts)
            np.savetxt(fname+'weightRRCount.txt',wRRCorCounts)
            if which_corr=='NN' and functype=='GM':
                np.savetxt(fname+'RNCount.txt',RNCorCounts)
                np.savetxt(fname+'normRNCount.txt',weightRNCounts)
            np.savetxt(fname+'normNNCount.txt',weightNNCounts)
            np.savetxt(fname+'normNRCount.txt',weightNRCounts)
            np.savetxt(fname+'normRRCount.txt',weightRRCounts)
            np.savetxt(fname+'Result.txt',NN_w)
            #np.savetxt(fname+'mass.txt',table['Mass'])


def corFunc(data,sightBins,rscale,nbins,min_sep,max_sep,rpar_min,rpar_max,min_box,max_box,which_corr,logfile,savefile,combine,functype,fname2=None,col_name2=None,key=None,key2=None,cuts2=None,data2=None):
    rpar_step=rParStepCal(rpar_min,rpar_max,sightBins)
    log_file=open(logfile,'r').read()
    log=open(logfile,'w')
    log.write(log_file)
    log.write('LoS Bin Step:'+str(rpar_step)+'\n')
    np.random.seed(0)
#     if False:
#         if which_corr=='NN' and fname2==None:
#             log.write('Number of objects of first file: '+str(len(data['x']))+'\n')
#             catN,catN2,catR=genCat(data,rscale,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,min_box,max_box,savefile,logfile,combine,which_corr,functype)
#             pairCounts(catN,catN2,catR,table,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,savefile,logfile,combine,which_corr,functype)
#         if which_corr=='NG' and fname2==None:
#             data=noPadding(data,'ax')
#             log.write('Number of objects of first file after removing padding: '+str(len(data['x']))+'\n')
#             catG,catN,catR=genCat(data,rscale,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,min_box,max_box,savefile,logfile,combine,which_corr,functype)
#             pairCounts(catG,catN,catR,data,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,savefile,logfile,combine,which_corr,functype)
#         if fname2 !=None:
#             log.write('------2ND FILE------\n')
#             data2=dat2 #read_Data_hdf5(fname2,col_name2,cuts2,key2)
#             data2['x']=data2['x']/1000
#             data2['y']=data2['y']/1000
#             data2['z']=data2['z']/1000
#             if which_corr=='NN':
#                 log.write('Number of objects of first file: '+str(len(data['x']))+'\n')
#                 log.write('Number of objects of second file: '+str(len(data2['x']))+'\n')
#                 catN,catN2,catR=genCat(data,rscale,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,min_box,max_box,savefile,logfile,combine,which_corr,functype,table2=data2)
#                 pairCounts(catN,catN2,catR, data,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,savefile,logfile,combine,which_corr,functype,table2=data2)
#                 log.close()

#             if which_corr=='NG':
#                 data=noPadding(data,'ax')
#                 log.write('Number of objects of first file after removing padding: '+str(len(data['x']))+'\n')
#                 log.write('Number of objects of second file after removing padding: '+str(len(data2['x']))+'\n')
#                 log.write('Mass: '+ str(data2['Mass'])+'\n')
#                 catG,catN,catR=genCat(data,rscale,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,min_box,max_box,savefile,logfile,combine,which_corr,functype,table2=data2)
#                 pairCounts(catG,catN,catR,data,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,savefile,logfile,combine,which_corr,functype,table2=data2)
#                 log.close()
    
    if which_corr=='NG':
        data=noPadding(data,'ax')
    log.write('Number of objects of first file: '+str(len(data['x']))+'\n')
    if fname2==None:
        print('test corFunc no 2nd File',file=log)
        catG,catN2,catR=genCat(data,rscale,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,min_box,max_box,savefile,log,combine,which_corr,functype)
        pairCounts(catN2,catG,catR,data,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,savefile,log,combine,which_corr,functype)
    else:
        print('test corFunc 2nd File',file=log)
        log.write('------2ND FILE-------\n')
        log.write('Number of objects of second file: '+str(len(data2['x']))+'\n')
        catG,catN2,catR=genCat(data,rscale,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,min_box,max_box,savefile,log,combine,which_corr,functype,table2=data2)
        pairCounts(catN2,catG,catR,data,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,savefile,log,combine,which_corr,functype,table2=data2)

# In[ ]:



