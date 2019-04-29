
# coding: utf-8

# In[37]:
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

def projShapeCal(table):
    ax=table['ax']
    ay=table['ay']
    q=table['q']
    phi=np.arctan2(ay,ax)
    e1=np.divide(1-q,1+q)*np.cos(2*phi)
    e2=np.divide(1-q,1+q)*np.sin(2*phi)
    table['e1']=e1
    table['e2']=e2
    #print(table['e1'],table['e2'])
    return table

def CorProcess(catN,catG,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,
		which_corr='NG'):
    corList=list()
    pairCounts=np.zeros((sightBins,nbins))
    wpairCounts=np.zeros((sightBins,nbins))
    if which_corr=='NG':
        for i in np.arange(sightBins):
            corList.append(treecorr.NGCorrelation(nbins=nbins,min_sep=min_sep,max_sep=max_sep,metric='Rperp',
                                             min_rpar=i*rpar_step+rpar_min,max_rpar=(i+1)*rpar_step+rpar_min))
            corList[i].process(catN,catG,metric='Rperp')
            pairCounts[i]=corList[i].npairs
            wpairCounts[i]=corList[i].weight
            #print('count',i,pairCounts[i])
        #print('catG',catG.x,'end')
        #print('catN',catN.x,'end')
        return corList,pairCounts,wpairCounts
    if which_corr=='NN':
        for i in np.arange(sightBins):
            corList.append(treecorr.NNCorrelation(nbins=nbins,min_sep=min_sep,max_sep=max_sep,min_rpar=i*rpar_step+rpar_min,max_rpar=(i+1)*rpar_step+rpar_min))
            corList[i].process(catN,catG,metric='Rperp')
            pairCounts[i]=corList[i].npairs
            wpairCounts[i]=corList[i].weight
            #print('count',i,corList[i].npairs)
        print('cat1', catN.x, 'end')
        print('cat2', catG.x, 'end')
        return corList,pairCounts,wpairCounts

def NG(table,rscale,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,min_box,max_box,fname,combine,extraTable=False,tableN=None):
    table=projShapeCal(table)
    rand=np.random.rand(3,len(table['x'])*rscale)*max_box+min_box
    catR=treecorr.Catalog(x=rand[0],y=rand[1],z=rand[2])
    RRCor=list()
    rand_pcounts=np.zeros((sightBins,nbins))
    RRCor=CorProcess(catR,catR,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,which_corr='NN')
    for i in np.arange(sightBins):
        rand_pcounts[i]=RRCor[i].npairs
    catG=treecorr.Catalog(x=table['x'],y=table['y'],z=table['z'],g1=table['e1'],g2=table['e2'])
    NRCor,NRCorCounts=CorProcess(catR,catG,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min)
    if extraTable==False:
        NGCor,NGCorCounts=CorProcess(catG,catG,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min)
    else:
        if combine==0:
            catN=treecorr.Catalog(x=tableN['x'],y=tableN['y'],z=tableN['z'])
        else:
            catN=treecorr.Catalog(x=tableN['x'],y=tableN['y'],z=tableN['z'],w=tableN['Mass'])
        NGCor,NGCorCounts=CorProcess(catN,catG,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min)
    NG_ksi=np.zeros((sightBins,nbins))
    if extraTable==False:
        for i in np.arange(sightBins):
            NG_ksi[i]=div0(NGCorCounts[i]/(len(table['x'])**2)-NRCorCounts[i]/(len(table['x'])*len(rand[0])),rand_pcounts[i]/len(rand[0])**2)
    else:
        for i in np.arange(sightBins):
            NG_ksi[i]=div0(NGCorCounts[i]/(len(table['x'])*catN.sumw)-NRCorCounts[i]/(len(table['x'])*len(rand[0])),rand_pcounts[i]/len(rand[0])**2)
            #print('NG:',NGCorCounts[i],'weight:',len(table['x'])*catN.sumw,'NR:',NRCorCounts[i],'weight:', len(rand[0])*len(table['x']),'RR:',rand_pcounts[i],'weight:',len(rand[0])**2)
    NG_w=np.sum(NG_ksi,axis=0)
    np.savetxt(fname+'.txt',NG_w)

def NN(table,rscale,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,min_box,max_box,fname,logfile,combine,functype,table2=None):
    rand=np.random.uniform(min_box,max_box,size=(3,len(table['x'])*rscale))
    rand2=np.random.uniform(min_box,max_box,size=(3,len(table['x'])*rscale*100))
    print('Before:'+str(rand2.size))
    rand2=rand2[:,::100]
    print('After:'+str(rand2.size))
    rand_cat2=treecorr.Catalog(x=rand[0],y=rand[1],z=rand[2])
    rand_cat=treecorr.Catalog(x=rand2[0],y=rand2[1],z=rand2[2])
    if functype=='MM':
        catN=treecorr.Catalog(x=table['x'],y=table['y'],z=table['z'],w=table['Mass'])
        catN2=catN
        print('NNWeight:'+str(np.sum(table['Mass'])**2))
        print('NRWeight:'+str(np.sum(table['Mass'])*len(rand[0])))
        print('RRWeight:'+str(len(rand[0])**2))
    if functype=='GM':
        catN=treecorr.Catalog(x=table['x'],y=table['y'],z=table['z'])
        catN2=treecorr.Catalog(x=table2['x'],y=table2['y'],z=table2['z'],w=table2['Mass'])
        print('NNWeight:'+str(np.sum(table2['Mass'])*len(table['x'])))
        print('NRWeight:'+str(len(table['x'])*len(rand[0])))
        print('RRWeight:'+str(len(rand[0])**2))
        print('RNWeight:'+str(len(rand[0])*np.sum(table2['Mass'])))
    if functype=='GG':
       catN=treecorr.Catalog(x=table['x'],y=table['y'],z=table['z'])
       catN2=catN
       print('NNWeight:'+str(len(table['x'])**2))
       print('NRWeight:'+str(len(table['x'])*len(rand[0])))
       print('RRWeight:'+str(len(rand[0])**2))
    NNCor,NNCorCounts,wNNCorCounts=CorProcess(rand_cat2,rand_cat2,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,which_corr='NN') #catN,catN2
    NRCor,NRCorCounts,wNRCorCounts=CorProcess(rand_cat2,rand_cat,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,which_corr='NN') #catN,rand_cat
    RRCor,RRCorCounts,wRRCorCounts=CorProcess(rand_cat,rand_cat,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,which_corr='NN') #rand_cat,rand_cat
    RNCor,RNCorCounts,wRNCorCounts=CorProcess(catN2,rand_cat,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,which_corr='NN') #catN2,rand_cat
    weightRNCounts=np.zeros((sightBins,nbins))
    #print('NNCor: ', NNCor, 'NRCor: ', NRCor, 'RRCor: ', RRCor)
    NN_ksi=np.zeros((sightBins,nbins))
    for i in np.arange(sightBins):
        #print('test')
        #print(i)
        #print('NN: ',NNCorCounts[i]/(catN.sumw)**2,'NR: ',NRCorCounts[i]/(catN.sumw*len(rand[0])),'RR: ',RRCorCounts[i]/len(rand[0])**2,'sum/subtract: ', NNCorCounts[i]/(catN.sumw)**2-2*NR)
        if functype=='GG':
            weightNNCounts=wNNCorCounts[i]/(len(table['x'])**2)
            weightNRCounts=wNRCorCounts[i]/(len(table['x'])*len(rand[0]))
            weightRRcounts=wRRCorCounts[i]/(len(rand[0])**2)
            NN_ksi[i]=div0(weightNNCounts-2*weightNRCounts+weightRRCounts,weightRRCounts)
        if functype == 'GM':
            weightNNCounts=wNNCorCounts[i]/(len(table['x'])*catN2.sumw)
            weightNRCounts=wNRCorCounts[i]/(len(table['x'])*len(rand[0]))
            weightRRCounts=wRRCorCounts[i]/(len(rand[0])**2)
            weightRNCounts=wRNCorCounts[i]/(catN2.sumw*len(rand[0]))
            NN_ksi[i]=div0(weightNNCounts-weightNRCounts-weightRNCounts+weightRRCounts,weightRRCounts)
        if functype == 'MM':
            weightNNCounts=wNNCorCounts[i]/(rand_cat2.sumw**2)     #/(catN.sumw**2)
            weightNRCounts=wNRCorCounts[i]/(rand_cat.sumw*rand_cat2.sumw)     #/(catN.sumw*len(rand[0]))
            weightRRCounts=wRRCorCounts[i]/(rand_cat.sumw**2)     #/(len(rand[0])**2)
            #print('NN: ', NNCorCounts[i], 'NR: ', NRCorCounts[i], 'RR: ', RRCorCounts[i])
            #print('weighted NN: ', weightNNCounts, 'weighted NR: ', weightNRCounts, 'weighted RR: ', weightRRCounts, 'add/subtract: ', weightNNCounts-2*weightNRCounts+weightRRCounts)
            NN_ksi[i]=div0(weightNNCounts-2*weightNRCounts+weightRRCounts, weightRRCounts)
            #print(NN_ksi[i])
    NN_w=np.sum(NN_ksi,axis=0)
    np.savetxt(fname+'bin.txt',NNCor[0].logr)
    np.savetxt(fname+'NNCount.txt',NNCorCounts)
    np.savetxt(fname+'NRCount.txt',NRCorCounts)
    np.savetxt(fname+'RRCount.txt',RRCorCounts)
    np.savetxt(fname+'RNCount.txt',RNCorCounts)
    np.savetxt(fname+'normNNCount.txt',weightNNCounts)
    np.savetxt(fname+'normNRCount.txt',weightNRCounts)
    np.savetxt(fname+'normRRCount.txt',weightRRCounts)
    np.savetxt(fname+'normRNCount.txt',weightRNCounts)
    np.savetxt(fname+'Result.txt',NN_w)
    np.savetxt(fname+'mass.txt',table['Mass'])
    #return NNCounts,NRCounts,RNCounts,RRCounts

def makeLog(logfile,func):
    reset=open(logfile,'w')
    reset.write('------BEGIN------')
    reset.write('------FUNCTION-------')
    reset.write(func)
    reset.write('------WRITE------')
    reset.close()

def corFunc(data,sightBins,rscale,nbins,min_sep,max_sep,rpar_min,rpar_max,min_box,max_box,func,logfile,savefile,combine,functype,fname2=None,col_name2=None,key=None,key2=None,cuts2=None,dat2=None):
    rpar_step=rParStepCal(rpar_min,rpar_max,sightBins)
    log_file=open(logfile,'r').read()
    log=open(logfile,'w')
    log.write(log_file)
    log.write('LoS Bin Step:'+str(rpar_step)+'\n')
    np.random.seed(0)
    if func=='NN' and fname2==None:
        log.write('Number of objects of first file: '+str(len(data['x']))+'\n')
        #NNCounts,NRCounts,RNCounts,RRCounts=
        NN(data,rscale,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,min_box,max_box,savefile,logfile,combine,functype)
        log.close()
        #return NNCounts,NRCounts,RNCounts,RRCounts
    if func=='NG' and fname2==None:
        data=noPadding(data,'ax')
        log.write('Number of objects of first file after removing padding: '+str(len(data['x']))+'\n')
        NG(data,rscale,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,min_box,max_box,savefile,combine)
        log.close()
    if fname2!=None:
        log.write('------2ND FILE------\n')
        data2=dat2 #read_Data_hdf5(fname2,col_name2,cuts2,key2)
        data2['x']=data2['x']/1000
        data2['y']=data2['y']/1000
        data2['z']=data2['z']/1000
        if func=='NN':
            log.write('Number of objects of first file: '+str(len(data['x']))+'\n')
            log.write('Number of objects of second file: '+str(len(data2['x']))+'\n')
            NN(data,rscale,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,min_box,max_box,savefile,logfile,combine,functype,table2=data2)
            log.close()
            
        if func=='NG':
            data=noPadding(data,'ax')
            log.write('Number of objects of first file after removing padding: '+str(len(data['x']))+'\n')
            log.write('Number of objects of second file after removing padding: '+str(len(data2['x']))+'\n')
            NG(data,rscale,sightBins,nbins,min_sep,max_sep,rpar_step,rpar_min,min_box,max_box,savefile,combine,extraTable=True,tableN=data2)
            log.close()


# In[ ]:



