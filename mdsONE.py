# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 11:58:19 2020

@author: Wisam
"""

import sys
import Oligo.File
import Oligo
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr 
import pandas as pd


#%%
"""
METHODS:
    GET DATA, SPLIT IT
    read_chromo_from_genome
    region_get_seq
    bucLim
    split_region_seq
"""

#GETS CHROMO, RETURNS A CHROMOSOME OBJECT
def read_chromo_from_genome(species, seqfile, chromoname):
    genome = Oligo.File.read_genome(species, seqs_filename=Oligo.File.search(seqfile))
    for chromo in genome:
        if str(chromo) == species + " " + chromoname:
            return chromo 
#GETS SEQUENCE OF A CHROMO REGION

def region_get_seq(start, end, chromo):
    seq=chromo.get_seq()[start:end]
    return seq 

#DEFINES BUCKETS OF INTEREST DEPENDING ON THE REGION COORDINATES, RETURNS BUCKET RANGE OF A SUSPECTED REGION
def bucLim(start, end, bucsize):
    bucFirst=start//bucsize
    bucLast=(end-1)//bucsize + 1
    return int(bucFirst), int(bucLast) 

#SPLITS MY REGION SEQUENCE IN BUCKETS EACH OF THEM CONTAINING THE CORRESPONDING SEQUENCE
#RETURNS A LIST OF SEQUENCES OF BUCKETS
#IMPORTANT: BEGIN FROM 1st BASE of 1st BUC 
def split_region_seq(start, end, chromoOfRegion, bucsize): 
    bucFirst, bucLast = bucLim(start, end, bucsize)
    startRegion=bucFirst*bucsize
    endRegion= bucLast*bucsize
    seqRegion = chromoOfRegion.get_seq()[startRegion:endRegion]
    seqBuckets=[]
    for i in range(0,bucLast-bucFirst):
        seqBuckets.append([bucFirst+i,seqRegion[i*bucsize:(i+1)*bucsize]])
    return seqBuckets 
#%%
"""
METHODS:
    RAW KMER/HIC_DATA
    kmer_data_creator
    kmer_read_spectra
    hic_plot
    kmer_plot
"""
#LOOKS FOR ALL KMER DISTRIBUTIONS IN A GIVEN LIST CONTAINING SEQUENCES OF BUCKETS
def kmer_data_creator(species, optiseq, cID,start,end,bucsize,k):
    chromo=read_chromo_from_genome(species, optiseq, cID)
    listseq= split_region_seq(start,end,chromo,bucsize)
    for i in range(len(listseq)):
        Oligo.Search.search_kmer(data_seq=listseq[i][1], k=k, data_seq_name="Bucket" + str(listseq[i][0]), output_filename="results/kmer"+str(k)+"/"+cID+"_"+str(bucsize)+"_" + str(listseq[i][0])+".kmer")

#READS THE SEQUENCES OF THE LIST DESCRIBED ABOVE, RETURNS A LIST OF SPECTRA
def kmer_read_spectra(optiseq,species,start,end,cID,bucsize,k):
    chromo=read_chromo_from_genome(species, optiseq, cID)
    listseq= split_region_seq(start,end,chromo,bucsize)
    spectra=[]
    for i in range(len(listseq)):
        spectrum = Oligo.Kmer.KmerSpectrum.read("results/kmer"+str(k)+"/"+cID+"_"+str(bucsize)+"_" + str(listseq[i][0])+".kmer",name= str(listseq[i][0]))
        spectra.append(spectrum)
    return spectra

#PLOTS THE HIC DATA AFTER DEFINING THE START AND END POINTS OF THE ANALYZED BUCKETS
def hic_plot(dataset, cID, name,bucfirst,buclast):
    fig=plt.figure(figsize=(12,11))#36,33
    axes1 = fig.add_axes([0.1,0.1,0.9,0.8])
    matrix = Oligo.Matrix.HiCMatrix.read("data/HiC_Data/"+dataset+"_observed/HiCtool_"+cID+"_40kb_observed.txt", resolution=bucSize)
    #bucFirst, bucLast = bucLim(mdsStart,mdsEnd,bucSize)
    data=matrix.sub_matrix(bucfirst, buclast).complete_matrix()
    sns.heatmap(data,cmap="hot")
    plt.title(dataset+" "+ name +" HiC 40k resolution")
    plt.savefig("results/"+dataset+"_"+name+"_Hic.png")

#PLOTS THE KMER MATRIX AFTER GETTING THE DATA OF THE KMER DISTRIBUTION OF EACH BUCKET    
def kmer_plot(optiseq, name,species,start,end,cID,bucsize,k, figsize):
    chromo=read_chromo_from_genome(species, optiseq, cID)    
    spectra=kmer_read_spectra(optiseq,species,start,end,cID,bucsize,k)
    matrix = Oligo.Kmer.KmerSpectrum.correlate_spectra(spectra)
    matrix.save("results/"+cID+"_"+name+".matrix")
    #listseq=split_region_seq(mdsStart,mdsEnd,chromo,bucSize)
    drawer = Oligo.Plot.HeatmapDrawer(data=matrix.matrix, xticklabels=matrix.col_names, yticklabels=matrix.col_names, vmin=0.5, vmax=1.0)#, show_values=show_values)
    drawer.plot("results/"+cID+"_"+name+".png", figsize=figsize)

#%%
"""
METHODS:
    SPECTRA CALCULATION
    spectra_diffs
    spectra_diffs_list
"""

#SPECTRA DIFFERENCES CALCULATED AND DISPLAYED
def spectra_diffs(mainspectrum, mainlen, optiseq,species,start,end,cID,bucsize,k,limit,difflist,rule):
    chromo=read_chromo_from_genome(species, optiseq, cID)    
    spectra=kmer_read_spectra(optiseq,species,start,end,cID,bucsize,k)
    words=sorted(mainspectrum.data.keys())
    
    #print(len(words))
    diffs=difflist
    for i in range(len(spectra)):
        """
        diff=0.0
        for j in words:
            mainval=float(mainspectrum.data.get(j))/mainlen
            if i.data.get(j)==None:
                continue
            val=float(i.data.get(j))/bucsize
            diff=diff+np.abs(mainval-val)
        diffs.append(diff)"""
        diff=diffs[i]
        if diff>limit:
            #print i.name, diff
            wordoi=[]
            for w in words:
                mainval=float(mainspectrum.data.get(w))/mainlen
                if spectra[i].data.get(w)==None:
                    continue
                val=float(spectra[i].data.get(w))/bucsize
                #if (val/mainval -1)*100<-67 or (val/mainval -1)*100>200:
                wordoi.append([w,round((val/mainval -1)*100,2)]) 
                ###percentage relative to c5 ave
                    #print i.name, j, mainval, val
            #print i.name, wordoi 
            #print(rule)
            f=open("results/boi/"+cID+"_"+spectra[i].name+".kmerdiff","w")#_str(rule).kmerdiff gelöscht
            for word in  wordoi:
                #print word[1]
                mainspectrum.data
                mainval=float(mainspectrum.data.get(word[0]))/mainlen
                f.write(word[0]+": "+str(round(mainval,5))+"  " + str(word[1])+"% \n")
            f.close()
            
    h = open("results/kmerdiff_" + cID + "_" + str(start) + "_"+ str(end) + ".txt","w")
    for i in diffs:
        h.write(str(diffs.index(i)) + " " + str(round(i,4)) + "\n")
    h.close()
    return diffs

#SPECTRA DIFFERENCES CALCULATED AND DISPLAYED
def spectra_diffs_list(mainspectrum, mainlen, optiseq,species,start,end,cID,bucsize,k):
    spectra=kmer_read_spectra(optiseq,species,start,end,cID,bucsize,k)
    words=sorted(mainspectrum.data.keys())
    
    print(len(words))
    diffs=[]
    for i in spectra:
        diff=0.0
        for j in words:
            mainval=float(mainspectrum.data.get(j))/mainlen
            if i.data.get(j)==None:
                continue
            val=float(i.data.get(j))/bucsize
            diff=diff+np.abs(mainval-val)
        diffs.append(diff)
    return diffs
    
#%%
"""
METHODS:
    ANALYSIS OF RAW DATA
    bucChain
    gcCalc
    chainBucRegion
    chainBucsPlot
"""    
#LOOKS FOR A SPECIFIC CHAIN IN A GIVEN BUCKET
def bucChain(bucID, chromo, chainSeed,bucsize, cID, pur, n):
    bucSeq = region_get_seq(bucID*bucSize, (bucID+1)*bucSize, chromo)
    Oligo.Search.search_chains(data_seq=bucSeq, chain_seed_seq=chainSeed, output_filename="results/chain/"+str(bucID)+'chain_search_'+ str(bucsize)+chainSeed +cID+ '.loci', search_type=0, min_purity=pur, min_length=n)
  
#CALCS GC CONTENT IN A GIVEN BUCKET (SEQ)
def gcCalc(seq):
    count = 0
    for i in seq:
        if i == "G" or i == "C":
            count = count+1
    gc = float(count)/float(len(seq))
    return gc #RETURN GC CONTENT OF A SEQUENCE
            
#LOOKS FOR CHAINS IN A LIST OF BUCKETS
def chainBucRegion(bucs,seed,species, seqfile, chromoname,bucsize, pur,n):
    chromo = read_chromo_from_genome(species, seqfile, chromoname)
    for i in range(len(bucs)):
        bucChain(int(bucs[i]),chromo,seed,bucsize,chromoname, pur,n)
        
#READS THE DATA AND PLOTS IT
def chainBucsPlot(inputseed, bucs,bucsize,cID,cutoff):
    y = []
    for i in bucs:
        address = "results/chain/"+str(int(i))+'chain_search_'+ str(bucsize)+inputseed+cID + '.loci'
        loci = Oligo.Locus.read(address)
        y.append(len(loci))
    fig= plt.figure()
    plt.figure(figsize=(12,10))
    plt.plot(bucs,y,"-.o")
    plt.grid(True)
    plt.title("Number of "+ inputseed + "-Elements in Buckets " + str(bucs[0])+"-"+str(bucs[-1]))
    outliers=[inputseed]
    plt.savefig("results/figures/"+'chain_search_'+ str(bucsize)+inputseed+cID+"_"+str(bucs[0])+"-"+str(bucs[-1]) + '.png')
    for i in range(len(y)):
        if y[i]>cutoff:
            outliers.append(bucs[i])
    return y,outliers # LIST OF CHAIN COUNTS IN BUCKETS, LIST OF BUCKET IDS THAT ARE OUTLIERS

#%%
"""
METHODS:
    PLOTTING
    genomeElBucs
    gcRegion
"""
#READS THE REPEATMASKER DATA AND PLOTS IT    
#RETURN: LIST OF ELEMENT COUNTS IN BUCKETS, LIST OF BUCKET IDS THAT ARE OUTLIERS
def genomeElBucs(chromo,directory, el,start,end,bucsize,cutoff):
    s,e = bucLim(start,end,bucsize)
    bucs = np.arange(s,e)
    y = [0]*(e-s)
    loci = Oligo.Locus.read("data/"+directory+"/Homo sapiens "+chromo+"_"+el+"_repeatmasker.loci")
    for l in loci:
        start = l.start
        buc = start//bucSize
        for b in range(len(bucs)):
            if buc==bucs[b]:
                y[b]=y[b]+1
    fig= plt.figure()
    plt.plot(bucs,y,"-.o")
    plt.grid(True)
    plt.title("Number of "+ el + "-Elements in Buckets " + str(bucs[0])+"-"+str(bucs[-1]))
    outliers=[el]
    plt.savefig("results/figures/"+chromo+'_'+ el + "_" + str(bucsize)+ '.png')
    for i in range(len(y)):
        if y[i]>cutoff:
            outliers.append(bucs[i])
    return y, outliers 
    
#PLOTS GC CONTENT IN THE BUCKETS
#REETURNS LIST OF GC CONTENT VALUES OF BUCKETS
def gcRegion(start,end,bucsize,species, optiseq,cID):
    chromo=read_chromo_from_genome(species, optiseq, cID)
    bigseq = region_get_seq(int(start), int(end), chromo)
    s,e = bucLim(start,end,bucsize)
    bucs = np.arange(s,e)
    gcVals=[]
    for i in bucs:
        gcBuc = gcCalc(bigseq[bucsize*(i-s):bucsize*(i+1-s)])
        gcVals.append(gcBuc)
    #print gcVals
    plt.figure()
    plt.plot(bucs,gcVals,"-.o")
    plt.grid(True)
    plt.title("GC content in Buckets " + str(s)+"-"+str(e))
    plt.savefig("results/figures/GC_"+cID+"_"+str(bucsize)+"_"+str(s)+"-"+str(e)+ '.png')
    return gcVals  

#%%
"""
METHODS:
    ADVANCED DATA ANALYSIS
    averageCorr
    wordGrouping
    wordGroupIncrease
    getStdWord
    getKmerStd
    allwordGroups
    bucketFilter
"""
#CALCULATES THE AVERAGES CORRELATION VALUE OF A BUCKET BY AVERAGING THE VALUES OF THE CORRESPONDING MATRIX ROW AND COLUMN OF A BUCKET
#DATA STRUCTURE: [[bucketID, aveCorr],[bucketID2, aveCorr2],...]
#RETURNS: LIST OF [BUCKET ID, AVE CORR] AND OUTLIER BUCKET ID ABOVE A CERTAIN THRESHOLD
def averageCorr(matrix,cutoff):
    means=[]
    x=[]
    for i in range(matrix.len_cols()):
        a=np.array(matrix.row(ix=i))
        b=np.array(matrix.col(ix=2))
        vals=np.concatenate((a,b))
        mean=np.mean(vals)
        means.append(mean)
        x.append(int(matrix.col_names[i]))
    plt.figure(figsize=(15,8))
    plt.plot(x, np.array(means),".")
    plt.grid(True)
    """"DECORATION"""
    
    plt.title("Average Correlation in Buckets "+ str(matrix.col_names[0]) + "-" + str(matrix.col_names[-1]))
    res=[]
    for i in range(len(x)):
        res.append(np.array([x[i],means[i]]))
    outliers=["AveCorr:"]
    for i in range(len(means)):
        if means[i]<cutoff:
            outliers.append(x[i])
    
    return np.array(res),outliers,means 

#groups of A rich or AT rich words for example
def wordGrouping(words,lim,char,char2=None,anti=False):
    woi=[]
    for word in words:
        if char2 == None:
            c=word.count(char)
            if c>=lim and anti==False:
                woi.append(word)
            elif c<lim and anti==True:
                woi.append(word)
        else:
            c=word.count(char)+word.count(char2)
            if c>=lim and anti==False:
                woi.append(word)
            elif c<lim and anti==True:
                woi.append(word)
    return woi
      
#Calculates Difference of Wordgroups in Buckets to the main spectrum      
def wordgroupIncrease(mainspectrum, kmerdiff, lim, letter1, letter2=None, anti=False,sigma=False):    
    words=sorted(mainspectrum.data.keys())
    #print(kmerdiff)
    g= wordGrouping(sorted(words),lim,letter1,char2=letter2,anti=anti)
    f=open(kmerdiff)#"results/boi/c5_3394.kmerdiff")
    weight=0
    summa=0
    for line in f:
        items=line.split(" ")
        if items[0][:-1] in g:
            #print items
            if sigma==False:
                weight=weight+float(items[1])*float(items[3][:-1])
            elif sigma==True:
                weight=weight+float(items[1])*float(items[3][:-6])
            summa=summa+float(items[1])   
    ave=summa/len(g)
    if summa == 0:
        aveweight = 0
    else:
        aveweight=weight/summa    
    #print lim, letter1, letter2
    #print "anti", anti, "sigma", sigma
    ######print g, len(g)
    #print aveweight
    return aveweight

#Saves the Difference of Words in STD Units after having the known difference values
def getStdWord(kmerdiffdir,mainspectrum,bucketsize,savefile,cID):
    words=sorted(mainspectrum.data.keys())
    results=[]
    for w in words:
        vals=[]
        for fpath in sorted(os.listdir(kmerdiffdir)):
            if cID+"_" in fpath:
                splitpath= fpath.split("_")
                if int(splitpath[1])==bucketsize:
                    f=open(kmerdiffdir+"/"+fpath)
                    for line in f:
                        items=line.split("\t")
                        if items[0]==w:
                            item=items[1].split("\n")
                            val=float(item[0])/bucketsize
                            vals.append(val)
                        
        ave=np.mean(np.array(vals))
        std=round(np.abs(np.std(np.array(vals))/ave*100),2)
        #std in Percent
        res=[w,ave,std]
        print(w, std)
        results.append(res)
    
    g=open(savefile,"w")
    for r in results:
        for rr in r:
            g.write(str(rr))
            g.write(" ")
        g.write("\n")
    g.close()
    
#Gets difference of kmer densities in Sigma Units
def kmerdiffStd(kmerdiffdir, resfile,cID):
    g=open(resfile)
    gres=[]
    #os.mkdir(kmerdiffdir+"std")
    for line in g:
        items=line.split(" ")
        gres.append([items[0],float(items[2])])
    g.close()  
    for fpath in sorted(os.listdir(kmerdiffdir)):
        if cID+"_" in fpath:
            f=open(kmerdiffdir+"/"+fpath)
            h=open(kmerdiffdir+"std/"+fpath+"2","w")
            for line in f:
                items=line.split(" ")
                #print items
                pc=float(items[3][:-1])
                for i in gres:
                    #print items
                    #print gres[0]
                    if items[0][:-1]==i[0]:
                        res=pc/float(i[1])
                        h.write(items[0]+" "+items[1]+"  "+str(round(res,2))+"sigma\n")
            f.close()
            h.close()
            
#Defines Subset of relevant wordgroups
def allwordGroups(spectrum,kmerdifffile,anti=False,sigma=False):
    cg= wordgroupIncrease(spectrum,kmerdifffile,5,"C",letter2="G",anti=anti,sigma=sigma)   
    at = wordgroupIncrease(spectrum,kmerdifffile,5,"A",letter2="T",anti=anti,sigma=sigma) 
    ca = wordgroupIncrease(spectrum,kmerdifffile,5,"C",letter2="A",anti=anti,sigma=sigma) 
    ct = wordgroupIncrease(spectrum,kmerdifffile,5,"C",letter2="T",anti=anti,sigma=sigma) 
    gt = wordgroupIncrease(spectrum,kmerdifffile,5,"G",letter2="T",anti=anti,sigma=sigma)  
    ga = wordgroupIncrease(spectrum,kmerdifffile,5,"G",letter2="A",anti=anti,sigma=sigma)
    c = wordgroupIncrease(spectrum,kmerdifffile,4,"C",anti=anti,sigma=sigma)
    g= wordgroupIncrease(spectrum,kmerdifffile,4,"G",anti=anti,sigma=sigma)
    a = wordgroupIncrease(spectrum,kmerdifffile,4,"A",anti=anti,sigma=sigma)
    t = wordgroupIncrease(spectrum,kmerdifffile,4,"T",anti=anti,sigma=sigma)
    return a,c,g,t,cg,at,ca,ct,gt,ga
    
#filters Buckets out that are not that interesting in their deviations
#NOTE only buckets with n-k+1=39996 kmers are considered!! 
def bucketFilter(cID, cIDdir,spec, rule,anti=False,sigma=False):
    wordfiles = os.listdir(cIDdir)
    #print(wordfiles)
    listboi=[]
    nonincluded=[]
    #print(listboi)
    for w in wordfiles:
        if (cID + "_40000_") in w and w.count("_")==2: #and ("_" + str(rule) + ".") in w
            #valsBucket = allwordGroups(spec,cIDdir+"/"+w,anti=False,sigma=False)
            ww = w.split(".")
            www = ww[0].split("_")
            #print(www[2])
            f = open("results/kmer5/"+cID+"_40000_"+str(www[2])+".kmer")
            freq=0
            allowed_chars = set('ACGT')
            for line in f:
                l = line.split()
                if set(l[0]).issubset(allowed_chars):
                    freq = freq + int(l[1])
                else:
                    freq=freq
            if freq == 39996:
                listboi.append(int(www[2]))
                #print freq
            else:
                nonincluded.append(int(www[2]))
    return sorted(listboi), sorted(nonincluded)

#%%
"""
VISUALS
"""
#VISUAL REPRESENTATION OF BOIS AND DEVIATION OF WORDGROUPS  BY HAND             
def visualBoi(cID,spec,limit,selectionrulelimit,length,mask,zoom=False,start=0,end=1):
    listboi5 = np.array(bucketFilter(cID,"results/kmer5",spec,selectionrulelimit)[0])
    nonincluded=np.array(bucketFilter(cID,"results/kmer5",spec,selectionrulelimit)[1])
    listboinew=[]
    for i in mask:
        if np.any(listboi5==i):
            listboinew.append(i)
    domains=[]
    finishStart=True
    #print(len(nonincluded))
    for i in range(len(nonincluded)):
        if finishStart==True: 
            domainStart = nonincluded[i]-0.5
        if  i+2 > len(nonincluded):
            domainEnd = nonincluded[i]+0.5
            finishStart = True
            domains.append([domainStart,domainEnd])
            break
        if finishStart==False:
            pass
        if nonincluded[i+1]==nonincluded[i]+1:
            finishStart=False
            continue
        else: 
            domainEnd = nonincluded[i]+0.5
            finishStart = True
            domains.append([domainStart,domainEnd])
            
    listboi5=listboinew
    #print(listboi5)
    ###
    
    la=[]
    lc=[]
    lg=[]
    lt=[]
    lcg=[]
    lat=[]
    lca=[]
    lct=[]
    lgt=[]
    lga=[]
    
    for i in listboi5:
        a,c,g,t,cg,at,ca,ct,gt,ga = allwordGroups(spec,"results/boi/"+cID+"_"+str(i)+".kmerdiff",anti=False)#"_"+str(selectionrulelimit)+ gelöscht
        la.append(a)
        lc.append(c)
        lg.append(g)
        lt.append(t)
        lat.append(at)
        lcg.append(cg)
        lca.append(ca)
        lct.append(ct)
        lgt.append(gt)
        lga.append(ga)
        #vals2 = allwordGroups(spec,"results/boistd/"+cID+"_"+str(i)+"_"+str(selectionrulelimit)+".kmerdiff2",anti=False,sigma=True)
    fig= plt.figure(figsize=(24,15))
    axes= fig.add_axes([0.1,0.1,0.8,0.8])
    ###CUSTOMIZATION
    axes.set_title("BoI's in " + cID + ", 95% of buckets filtered out", fontsize=28) #
    #axes.text(130,600,"BAC probes region\n(Falk et al.)", fontsize=18, color="red")
    #THRESHOLD TEXT
    #axes.text(0.5,0.9,"Threshold for filter: " + str(round(limit*100,3)) +
    #          "% relative average deviation\n from average chromosome spectrum for every kmer", fontsize=18,
    #          horizontalalignment='center', verticalalignment='center', transform = axes.transAxes)
    #plt.title("visual Boi of " + cID + " with boi-limit of " + str(round(limit*100,3)) + "% relative average deviation for every kmer, " + str(selectionrulelimit) + "% screened")
    xvals=np.array(listboi5)*0.04
    axes.set_xlabel("Bucket start coordinate [Mbp]", fontsize=20)
    axes.set_ylabel("Deviation of word group in % (up- and downside)", fontsize=20)
    #print(xvals)
    ###CUSTOMIZATION
    #axes.axvspan(128.952312, 147.2582, alpha=0.3, color='red')
    #axes.axvspan(46.485901,50.059807, alpha=0.3, color='yellow')
    #axes.text(43,500,"Centromer", fontsize=18, color="orange")
    axes.plot(xvals, la,".",label="A-rich")
    axes.plot(xvals, lc,".",label="C-rich")
    axes.plot(xvals, lg,".",label="G-rich")
    axes.plot(xvals, lt,".",label="T-rich")
    axes.plot(xvals, lat,".",label="AT-rich")
    axes.plot(xvals, lcg,".",label="CG-rich")
    length=length/1e6
    if zoom==False:
        axes.set_xlim(0,length)
    elif zoom==True:
        axes.set_xlim(start,end)
    ###CUSTOMIZATION
    axes.legend(loc="upper left", bbox_to_anchor=(0.05,1), 
                title="Word group", title_fontsize=20, fontsize=18)
    axes.tick_params(axis='x', labelsize=16)
    axes.tick_params(axis='y', labelsize=16)
    for i in domains:
        axes.fill_betweenx([min([*la, *lc, *lg, *lt, *lat, *lcg]),max([*la, *lc, *lg, *lt, *lat, *lcg])],
                           i[0]*40000/1e6, i[1]*40000/1e6 , color='red', alpha=.1)
    plt.show()
    fig.savefig("results/visualboi/visualBoiV2_"+cID+"_"+str(selectionrulelimit)+".png")          
    return listboi5

#%%
"""
RANDOM
"""
def getTADbucs(taddir,tadfile,bucsize):
    f=open(taddir+"/"+tadfile)
    bucs=[]
    for l in f:
        coords=l.split()
        buc=int(coords[0])//bucSize
        bucs.append(buc)
    return bucs
        
def splitIUPAC(seq):
    new=[]
    c=seq.count("N")
    if c==0:
        return(seq)
    else:
        new.append(splitIUPAC(seq.replace("N","A",1)))
        new.append(splitIUPAC(seq.replace("N","C",1)))
        new.append(splitIUPAC(seq.replace("N","G",1)))
        new.append(splitIUPAC(seq.replace("N","T",1))) 
    return new

def flattenNestedList(nestedList):
    ''' Converts a nested list to a flat list '''
    flatList = []
    # Iterate over all the elements in given list
    for elem in nestedList:
        # Check if type of element is list
        if isinstance(elem, list):
            # Extend the flat list by adding contents of this element (list)
            flatList.extend(flattenNestedList(elem))
        else:
            # Append the elemengt to the list
            flatList.append(elem)    
    return flatList
                    


#%%VARIABLES     
mdsC = "c5"
mdsChromo = "chr5"
mdsStart=0
mdsEnd=181538259
species="Homo sapiens"
bucSize=40000
dataset="GSM1551599" #HiC dataset
optiseqs= "data/Seq_Files/Homo_sapiens_opti.seqs"
purity=0.9
chainlen=10

#BACPROBE DATA
i15=[128952312,129126486]
g6=[132400890,132585084]
j22=[137514367,137706209]
b8=[141718264,141894845]
c10=[147078981,147258200]

print(bucLim(i15[0],i15[1],40000))
print(bucLim(g6[0],g6[1],40000))
print(bucLim(j22[0],j22[1],40000))
print(bucLim(b8[0],b8[1],40000))
print(bucLim(c10[0],c10[1],40000))

#KMER Bucket data, here in kmer5 folder
#kmer_data_creator(species, optiseqs, mdsC,mdsStart,mdsEnd,bucSize,5)
"""
spectra=[]
for i in range(3000,4000):
    spectra.append(Oligo.Kmer.KmerSpectrum.read("results/kmer5/c5_40000_"+str(i)+".kmer", name=str(i)))
matrix = Oligo.Kmer.KmerSpectrum.correlate_spectra(spectra)
# Save Matrix
matrix.save("c5_3_4000.matrix")"""
#makes a plot of the kmer spectra data as a correlation matrix above, matrix is stored

#reads the above correlation matrix
"""
mat = Oligo.Matrix.Matrix.read("results/c5_3_4000.matrix")
#produces a list of average correlation of buckets [bucket, avecorr] with a list of outliers
res, oCorr,means = averageCorr(mat,0.7)
plt.savefig("avecorr_c53_4000.png")
plt.show()
f = open('results/c5_MDSfull.avecorr','w') 
for i in res:
    f.write(str(int(i[0]))+ " " + str(round(i[1],4)) + "\n")
f.close()    
  
print means

#bucket definition
#s,e = bucLim(mdsStart,mdsEnd,bucSize)
coord1=128e6
coord2=150e6
s,e = bucLim(coord1,coord2,bucSize)
bucs = np.arange(s,e)

#CHAIN DATA and OUTLIERS (also plotted)
#of all (N)n or (NN)n
#with seeds and thresholds
seeds=[["A",25],
       ["C",1],
       ["G",1],
       ["T",24],
       ["AC",9],
       ["AG",4],
       ["AT",15],
       ["CG",4],
       ["CT",4],
       ["GT",5]]
datasets=[]
for i in seeds:
    #chainBucRegion(bucs,i[0],species,"data/Seq_Files/Homo_sapiens_opti.seqs",mdsC,bucSize,purity,chainlen)
    data, outliers=chainBucsPlot(i[0],bucs,bucSize,mdsC,i[1])
    datasets.append([data, outliers])

#GENOME ELEMENT (Alu,L1) DATA PLOTTING, returns list of counts and outlier buckets
#dataAlu, oAlu = genomeElBucs(mdsC, "LOCI_alu_repeatmasker", "Alu",mdsStart,mdsEnd,bucSize,60)
#dataL1, oL1 = genomeElBucs(mdsC, "LOCI_L1_hg38","L1",mdsStart,mdsEnd,bucSize,30)

dataAlu, oAlu = genomeElBucs(mdsC, "LOCI_alu_repeatmasker", "Alu",coord1,coord2,bucSize,60)
dataL1, oL1 = genomeElBucs(mdsC, "LOCI_L1_hg38","L1",coord1,coord2,bucSize,30)
#same but with GC% without outliers
dataGCp = gcRegion(coord1,coord2,bucSize,species,optiseqs,mdsC)
#dataGCp = gcRegion(mdsStart,mdsEnd,bucSize,species,optiseqs,mdsC)

f=open("results/kmerdiff_c5_0_181538259.txt")
dataDev = []
for line in f:
    l = line.split(" ")
    if int(l[0]) >= 3200 and int(l[0]) < 3750:
        dataDev.append(float(l[1][:-2]))
#print dataDev
        

data = {'AveCorr': means,
        'Dev' : dataDev,
        'Alu': dataAlu,
        'L1': dataL1,
        'A': datasets[0][0],
        'C': datasets[1][0],
        'G': datasets[2][0],
        'T': datasets[3][0],
        'AC': datasets[4][0],
        'AG': datasets[5][0],
        'AT': datasets[6][0],
        'CG': datasets[7][0],
        'CT': datasets[8][0],
        'GT': datasets[9][0],
        'GC%': dataGCp
        }

#Pearson correlation of those datasets
df = pd.DataFrame(data,columns=['AveCorr', 'Dev', 'Alu','L1','A','C','G','T','AC','AG','AT','CG','CT','GT','GC%'])
corrMatrix = df.corr()
print (corrMatrix)
csv_data = corrMatrix.to_csv('results/corr_dataset_'+mdsC+"_"+str(s)+"_"+str(e)+'.csv', index = True) 
"""
#KMER SPECTRUM OF C5 AS A WHOLE
#Oligo.Search.search_kmer(data_seq=read_chromo_from_genome(species, optiseqs, "c9").get_seq(), k=5, data_seq_name="c9", output_filename="results/c9.kmer")
#Oligo.Search.search_kmer(data_seq=read_chromo_from_genome(species, optiseqs, "c8").get_seq(), k=5, data_seq_name="c8", output_filename="results/c8.kmer")
#Oligo.Search.search_kmer(data_seq=read_chromo_from_genome(species, optiseqs, "c14").get_seq(), k=5, data_seq_name="c14", output_filename="results/c14.kmer")
#Oligo.Search.search_kmer(data_seq=read_chromo_from_genome(species, optiseqs, "c2").get_seq(), k=5, data_seq_name="c2", output_filename="results/c2.kmer")
#Oligo.Search.search_kmer(data_seq=read_chromo_from_genome(species, optiseqs, "c19").get_seq(), k=5, data_seq_name="c19", output_filename="results/c19.kmer")
#Oligo.Search.search_kmer(data_seq=read_chromo_from_genome(species, optiseqs, "c7").get_seq(), k=5, data_seq_name="c7", output_filename="results/c7.kmer")
#Oligo.Search.search_kmer(data_seq=read_chromo_from_genome(species, optiseqs, "c15").get_seq(), k=5, data_seq_name="c15", output_filename="results/c15.kmer")
#Oligo.Search.search_kmer(data_seq=read_chromo_from_genome(species, optiseqs, "c17").get_seq(), k=5, data_seq_name="c17", output_filename="results/c17.kmer")
#for i in ["c1","c3","c4","c6","c10","c11","c12","c13","c16","c18","c20","c21","cX","cY"]:
#    Oligo.Search.search_kmer(data_seq=read_chromo_from_genome(species, optiseqs, i).get_seq(), k=5, data_seq_name=i, output_filename="results/"+i+".kmer")

spectrum_c5 = Oligo.Kmer.KmerSpectrum.read("results/c5.kmer",name= "c5")  #MDS
spectrum_c22 = Oligo.Kmer.KmerSpectrum.read("results/c22.kmer",name= "c22") #BCR
spectrum_c9 = Oligo.Kmer.KmerSpectrum.read("results/c9.kmer",name= "c9") #ABL
spectrum_c8 = Oligo.Kmer.KmerSpectrum.read("results/c8.kmer",name= "c8") #cMYC
spectrum_c14 = Oligo.Kmer.KmerSpectrum.read("results/c14.kmer",name= "c14") #IGH
spectrum_c2 = Oligo.Kmer.KmerSpectrum.read("results/c2.kmer",name= "c2") #IGK
spectrum_c19 = Oligo.Kmer.KmerSpectrum.read("results/c19.kmer",name= "c19") #ERCC
spectrum_c7 = Oligo.Kmer.KmerSpectrum.read("results/c7.kmer",name= "c7") #HOXa
spectrum_c15 = Oligo.Kmer.KmerSpectrum.read("results/c15.kmer",name= "c15") #PML
spectrum_c17 = Oligo.Kmer.KmerSpectrum.read("results/c17.kmer",name= "c17") #RARa

spectrum_c1 = Oligo.Kmer.KmerSpectrum.read("results/c1.kmer",name= "c1") #RARa
spectrum_c3 = Oligo.Kmer.KmerSpectrum.read("results/c3.kmer",name= "c3") #RARa
spectrum_c4 = Oligo.Kmer.KmerSpectrum.read("results/c4.kmer",name= "c4") #RARa
spectrum_c6 = Oligo.Kmer.KmerSpectrum.read("results/c6.kmer",name= "c6") #RARa
spectrum_c10 = Oligo.Kmer.KmerSpectrum.read("results/c10.kmer",name= "c10") #RARa
spectrum_c11 = Oligo.Kmer.KmerSpectrum.read("results/c11.kmer",name= "c11") #RARa
spectrum_c12 = Oligo.Kmer.KmerSpectrum.read("results/c12.kmer",name= "c12") #RARa
spectrum_c13 = Oligo.Kmer.KmerSpectrum.read("results/c13.kmer",name= "c13") #RARa
spectrum_c16 = Oligo.Kmer.KmerSpectrum.read("results/c16.kmer",name= "c16") #RARa
spectrum_c18 = Oligo.Kmer.KmerSpectrum.read("results/c18.kmer",name= "c18") #RARa
spectrum_c20 = Oligo.Kmer.KmerSpectrum.read("results/c20.kmer",name= "c20") #RARa
spectrum_c21 = Oligo.Kmer.KmerSpectrum.read("results/c21.kmer",name= "c21") #RARa
spectrum_cX = Oligo.Kmer.KmerSpectrum.read("results/cX.kmer",name= "cX") #RARa
spectrum_cY = Oligo.Kmer.KmerSpectrum.read("results/cY.kmer",name= "cY") #RARa


#kmer_data_creator(species, optiseqs, "c22",mdsStart,50818468,bucSize,5)
#kmer_data_creator(species, optiseqs, "c9",mdsStart,138394717,bucSize,5)
#kmer_data_creator(species, optiseqs, "c8",mdsStart,145138636,bucSize,5)
#kmer_data_creator(species, optiseqs, "c14",mdsStart,107043718,bucSize,5)
#kmer_data_creator(species, optiseqs, "c2",mdsStart,242193529,bucSize,5)
#kmer_data_creator(species, optiseqs, "c19",mdsStart,58617616,bucSize,5)
#kmer_data_creator(species, optiseqs, "c7",mdsStart,159345973,bucSize,5)
#kmer_data_creator(species, optiseqs, "c15",mdsStart,101991189,bucSize,5)
#kmer_data_creator(species, optiseqs, "c17",mdsStart,83257441,bucSize,5)
#for i in ["c1","c3","c4","c6","c10","c11","c12","c13","c16","c18","c20","c21","cX","cY"]:
#    l = read_chromo_from_genome(species, optiseqs, i).length
#    kmer_data_creator(species, optiseqs, i,mdsStart,l,bucSize,5)
    


#KMER HISTOGRAM OF C5
#kmers = sorted(spectrum.get_kmers())
#drawer = Oligo.Plot.HistoDrawer(bins=spectrum.get_bins(), label='k=5')
#drawer.plot('c5.kmer.png', xticklabels=kmers, xticks=xrange(len(kmers)), figsize=(10,3), xticklabels_rotation=90)

#KMER HISTOGRAM OF A BUCKET (3394)
#spectrum = Oligo.Kmer.KmerSpectrum.read("results/kmer5/c5_40000_3394.kmer",name= "c5 40000 3394")
#kmers = sorted(spectrum.get_kmers())
#drawer = Oligo.Plot.HistoDrawer(bins=spectrum.get_bins(), label='k=5')
#drawer.plot('c5_40000_3394.kmer.png', xticklabels=kmers, xticks=xrange(len(kmers)), figsize=(10,3), xticklabels_rotation=90)


#PRINTS A LIST OF KMERS THAT DEVIATE SIGNIFICANTLY IN BUCKETS OF INTEREST
"""
f = open("results/c5.kmer")
freq=0
allowed_chars = set('ACGT')
for line in f:
    l = line.split()
    if set(l[0]).issubset(allowed_chars):
        freq = freq + int(l[1])
g=open("results/difflists/c5.diffbucket")
difflist=[]      
for ll in g:
    difflist.append(float(ll))
diffsnp=np.array(difflist)
diffs = spectra_diffs(spectrum_c5,freq,optiseqs, species, mdsStart,mdsEnd,mdsC,bucSize,5,0.5,difflist,0)
plt.figure(figsize=(15,6))
plt.grid()
y=np.array(diffs)
plt.ylim(max(y)+0.1,min(y)-0.05)
plt.plot(np.arange(0,4539)*40000/1e6,y,".")
plt.xlabel("Starting coordinate in a bucket [Mbp]")
plt.title("Deviation of each bucket kmer spectrum from the chromosomal spectrum (chr 5)")
#data={"Avecorr": np.array(means),
#      "Deviation": np.array(y)}
#corrdata=pd.DataFrame(data,columns=['Avecorr','Deviation'])
#corrMatrix = corrdata.corr()
#print (corrMatrix)
plt.savefig("deviation_c5_3200_3750.png")
"""
#spectra_diffs(spectrum_c22,50818468,optiseqs, species, mdsStart,50818468,"c22",bucSize,5,0.5)
#spectra_diffs(spectrum_c9,138394717,optiseqs, species, mdsStart,138394717,"c9",bucSize,5,0.5)
#spectra_diffs(spectrum_c8,145138636,optiseqs, species, mdsStart,145138636,"c8",bucSize,5,0.5)
#spectra_diffs(spectrum_c14,107043718,optiseqs, species, mdsStart,107043718,"c14",bucSize,5,0.5)
#spectra_diffs(spectrum_c2,242193529,optiseqs, species, mdsStart,242193529,"c2",bucSize,5,0.5)
#spectra_diffs(spectrum_c19,58617616,optiseqs, species, mdsStart,58617616,"c19",bucSize,5,0.5)
#spectra_diffs(spectrum_c7,159345973,optiseqs, species, mdsStart,159345973,"c7",bucSize,5,0.5)
#spectra_diffs(spectrum_c15,101991189,optiseqs, species, mdsStart,101991189,"c15",bucSize,5,0.5)
#spectra_diffs(spectrum_c17,83257441,optiseqs, species, mdsStart,83257441,"c17",bucSize,5,0.5)

cs1=["c1",spectrum_c1,248956422]
cs2=["c2",spectrum_c2,242193529]
cs3=["c3",spectrum_c3,198295559]
cs4=["c4",spectrum_c4,190214555]
cs5=["c5",spectrum_c5,181538259]
cs6=["c6",spectrum_c6,170805979]
cs7=["c7",spectrum_c7,159345973]
cs8=["c8",spectrum_c8,145138636]
cs9=["c9",spectrum_c9,138394717]
cs10=["c10",spectrum_c10,133797422]
cs11=["c11",spectrum_c11,135086622]
cs12=["c12",spectrum_c12,133275309]
cs13=["c13",spectrum_c13,114364328]
cs14=["c14",spectrum_c14,107043718]
cs15=["c15",spectrum_c15,101991189]
cs16=["c16",spectrum_c16,90338345]
cs17=["c17",spectrum_c17,83257441]
cs18=["c18",spectrum_c18,80373285]
cs19=["c19",spectrum_c19,58617616]
cs20=["c20",spectrum_c20,64444167]
cs21=["c21",spectrum_c21,46709983]
cs22=["c22",spectrum_c22,50818468]
csX=["cX",spectrum_cX,156040895]
csY=["cY",spectrum_cY,57227415]
allcs=[cs1,cs2,cs3,cs4,cs5,cs6,cs7,cs8,cs9,cs10,cs11,cs12,cs13,cs14,cs15,cs16,cs17,cs18,cs19,cs20,cs21,cs22,csX,csY]
pcrule=95
#allcs=[cs5]
"""
for i in allcs:
    g=open("results/difflists/"+i[0]+".diffbucket","w") #saves diff of each bucket
    le = read_chromo_from_genome(species, optiseqs, i[0]).length
    f = open("results/"+i[0]+".kmer")
    freq=0
    allowed_chars = Set('ACGT')
        for line in f:
            l = line.split()
            if Set(l[0]).issubset(allowed_chars):
                freq = freq + int(l[1])
    difflist=spectra_diffs_list(i[1],freq,optiseqs, species, mdsStart,le,i[0],bucSize,5)
    for j in difflist:
        g.write(str(j)+"\n")
    diffsnp=np.array(difflist)
    limit95=np.percentile(diffsnp,pcrule)
    print(i[0],limit95)
    g.close()

"""
    

####
####
####
#IMPORTANT CODE
#GC content comparison
"""
for i in allcs:
    gcchromo= read_chromo_from_genome(species, optiseqs, i[0])
    buclims= bucLim(0, i[2], bucSize)
    bucs = split_region_seq(0, i[2], gcchromo, bucSize)
    fracs=[]
    for j in bucs:
        countc=j[1].count("C")
        countg=j[1].count("G")
        count=countc+countg
        fraction=count/bucSize
        if j[1].count("N")!=0:
            fraction=0
        fracs.append(fraction)
        #print(count)
    plt.figure(figsize=(15,10))
    plt.plot(np.arange(buclims[0],buclims[1])*bucSize/1e6,fracs,".")
    plt.xlabel("Start of the bucket [Mbp]")
    plt.ylabel("GC content fraction")
    plt.ylim(0.25,0.65)
    plt.title("GC content in chromosome 5, not fully-sequenced areas are omitted")
    plt.show()
"""
"""
data={"Deviation": np.array(y),
      "GC_content": np.array(fracs)}
corrdata=pd.DataFrame(data,columns=['Deviation',"GC_content"])
corrMatrix = corrdata.corr()
print (corrMatrix)
 """  
#VISUAL BOI PLOTTER
"""
boiBucs=[]
for i in allcs:
    g=open("results/difflists/"+i[0]+".diffbucket")
    #l = read_chromo_from_genome(species, optiseqs, i[0]).length
    difflist=[]      
    for ll in g:
        difflist.append(float(ll))
    diffsnp=np.array(difflist)
    #print(difflist)
    limit95=np.percentile(diffsnp,pcrule)
    maskedbuc=[]
    for j in range(len(diffsnp)):
        if diffsnp[j] > limit95:
            maskedbuc.append(j)
    maskedbuc=np.array(maskedbuc)
    f = open("results/"+i[0]+".kmer")
    freq=0
    allowed_chars = set('ACGT')
    for line in f:
        l = line.split()
        if set(l[0]).issubset(allowed_chars):
            freq = freq + int(l[1])
    print(freq)
    #spectra_diffs(i[1],freq,optiseqs, species, mdsStart,i[2],i[0],bucSize,5,limit95,difflist,pcrule)
    listboibuc=visualBoi(i[0],i[1],limit95,pcrule,i[2],maskedbuc)
    for j in listboibuc:
        boiBucs.append((i[0],j))
    #print(boiBucs)
#print(boiBucs)
n = open("results/boi_buckets.bucloc", "a")
for i in boiBucs:
    n.write (i[0] + " " + str(i[1]) + "\n")
n.close()
"""
"""   
g=open("results/difflists/"+allcs[4][0]+".diffbucket")
difflist=[]      
for ll in g:
    difflist.append(float(ll))
diffsnp=np.array(difflist)
#print(difflist)
limit95=np.percentile(diffsnp,pcrule)
maskedbuc=[]
for j in range(len(diffsnp)):
    if diffsnp[j] > limit95:
        maskedbuc.append(j)
maskedbuc=np.array(maskedbuc)
visualBoi(allcs[4][0],allcs[4][1],limit95,pcrule,allcs[4][2],maskedbuc,
          zoom=False,start=125,end=150)
"""
"""
new=splitIUPAC("CCACCAGGGGGC") 
allseqs=flattenNestedList(new)

seq1=read_chromo_from_genome(species, optiseqs, "c1")
listseq=split_region_seq(0,len(seq1),seq1,40000)
for i in range(len(listseq)):
    j="CCACCAGGGGGC"
    Oligo.Search.search_seq(data_seq=listseq[i][1], query_seq=j,output_filename="results/CTCF/"+str(i)+"_"+j+".loci",
                            allowed_missmatches=2,verbose=False)


tadbucs1=getTADbucs("data/GSM1909121_tad", "HiCtool_chr1_topological_domains.txt", 40000)
alltadvals=np.zeros(61)
for i in tadbucs1:
    
    seq=[j for j in range(i-30,i+31)]
    numbers=[]
    for j in seq:
        if j>0 and j<6223:
            number=0
            ii="CCACCAGGGGGC"
            loci = Oligo.Locus.read("results/CTCF/"+str(j)+"_"+ii+".loci")
            #f=open("results/CTCF/"+str(j)+"_"+ii+".loci")
            number=number+len(loci)
            numbers.append(number)
    for k in range(len(numbers)):
        alltadvals[k]=alltadvals[k]+np.array(numbers)[k]
#alltadvals=alltadvals/len(alltadvals)
plt.plot(alltadvals)

            found=False
            for l in f:
                if l.split()[0]=="CCGCG":
                    found=True
                    seqvals.append(float(l.split()[1]))
            if found==False:
                seqvals.append(0.)
    for k in range(len(seqvals)):
        alltadvals[k]=alltadvals[k]+np.array(seqvals)[k]
    #plt.figure()
    #plt.plot(seqvals)
alltadvals=alltadvals/(len(tadbucs1))
plt.plot(alltadvals)

            
"""


"""
for i in allcs:
    getStdWord("results/kmer5",i[1],40000,"results/stdWords"+i[0]+".kmerstd",i[0])#
for i in allcs:
    kmerdiffStd("results/boi", "results/stdWords"+i[0]+".kmerstd",i[0])



allwordGroups(spectrum_c5,"results/boi/c5_3394.kmerdiff",anti=False)
allwordGroups(spectrum_c5,"results/boistd/c5_3394.kmerdiff2",anti=False,sigma=True)

    
#print os.listdir("results/boistd")

for i in allcs:
    visualBoi(i[0],i[1],limit,pcrule)
"""

#TRASH
#TRASH
#TRASH
"""
plt.figure(figsize=(10,2))
plt.grid(True)
plt.plot(oA[1:],np.ones(len(oA[1:]))*0,".")
plt.plot(oC[1:],np.ones(len(oC[1:]))*1,".")
plt.plot(oG[1:],np.ones(len(oG[1:]))*2,".")
plt.plot(oT[1:],np.ones(len(oT[1:]))*3,".")
plt.plot(oAC[1:],np.ones(len(oAC[1:]))*4,".")
plt.plot(oAG[1:],np.ones(len(oAG[1:]))*5,".")
plt.plot(oAT[1:],np.ones(len(oAT[1:]))*6,".")
plt.plot(oCG[1:],np.ones(len(oCG[1:]))*7,".")
plt.plot(oCT[1:],np.ones(len(oCT[1:]))*8,".")
plt.plot(oGT[1:],np.ones(len(oGT[1:]))*9,".")
plt.plot(oCorr[1:],np.ones(len(oCorr[1:]))*10,".")
plt.plot(oAlu[1:],np.ones(len(oAlu[1:]))*11,".")
plt.plot(oL1[1:],np.ones(len(oL1[1:]))*12,".")
"""

#alternative representation of kmer corr matrix
#kmer_plot(optiseqs, "MDS_3200_3750", species, 128e6,150e6,mdsC,bucSize,5,(24,20))
#mat = Oligo.Matrix.Matrix.read("results/c5_MDS.matrix")
#plt.figure(figsize=(18,15))
#sns.heatmap(mat.matrix,cmap="hot",vmin=0.5, vmax=1.0)
#plt.savefig("results/c5_MDS_2.png")
"""
#hic_plot("GSM1551599", "chr5", "c5_3350_3450",3350,3450)
#hic_plot("GSM1551550", "chr5", "c5_3350_3450",3350,3450)
hic_plot("GSM1909121", "chr5", "c5_3350_3450",3350,3450)
#hic_plot("GSM1608505", "chr5", "c5_3350_3450",3350,3450)

"""
hic_plot("GSM1551599", "chr9", "c9_3250_3300",3250,3300)
hic_plot("GSM1909121", "chr9", "c9_3250_3300",3250,3300)

cenStart=np.array([122.0,92.1,90.7,49.7,46.4,58.5,58.1,44.0,43.2,39.6,51.0,34.7,
          16.0,16.0,17.0,36.3,22.8,15.4,24.4,26.4,10.8,12.9,58.6,10.3])
cenEnd=np.array([125.2,94.1,93.7,51.8,50.1,59.9,60.9,45.9,45.6,41.6,54.5,37.2,18.1,18.2,
        19.8,38.3,26.9,20.9,27.2,30.1,13.0,15.1,62.4,10.6])
cenStart=cenStart*1e6/40000
cenStart=cenStart.astype(int)
cenEnd=cenEnd*1e6/40000
cenEnd=cenEnd.astype(int)
chrs=["22"]#["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]

###BP FUSION GDB###
t=open("data/TCGA_ChiTaRS_combined_fusion_information_on_hg19.txt")

bucLoc=[]
for l in t:
    bp = l.split("\t") #chr:5, coord:6 chr2:9, coord2:10
    #print(bp)
    bucLoc.append(("c"+bp[5][3:],int(int(bp[6])/40000))) 
    bucLoc.append(("c"+bp[9][3:],int(int(bp[10])/40000))) 
t.close()
bpCorr=[]
#print(bucLoc)
n = open("results/bp_buckets.bucloc", "a")
for i in bucLoc:
    n.write (i[0] + " " + str(i[1]) + "\n")
n.close()


for chrn in range(len(chrs)) :
    matrix = Oligo.Matrix.HiCMatrix.read("data/HiC_Data/GSM1551599_observed/HiCtool_chr"+chrs[chrn]+"_40kb_observed.txt", resolution=bucSize)
    #bucFirst, bucLast = bucLim(mdsStart,mdsEnd,bucSize)
    data=np.array(matrix.complete_matrix())#.sub_matrix(3350,3450)
    print(data)
    means=[]
    x=[]
    for i in range(len(data)):
        a=np.array(data[i][i-5:i+6])
        #b=np.array(data[:][i])
        vals=a#np.concatenate((a,b))
        mean=np.mean(vals)#np.sum(vals[:10])-np.sum(vals[10:])#np.mean(vals)
        means.append(mean)
        x.append(int(i)*40000/1e6)
    """plt.figure(figsize=(12,8))
    plt.plot(x, np.array(means),"-.")
    plt.grid(True)
    plt.xlabel("Coordinate of the bucket [Mbp]")
    plt.ylabel("Average HiC-frequency in the neighborhood")
    plt.title("Healthy cell(GSM1551599) chr"+chrs[chrn]+": Average Freq in Buckets "+ str(x[0]) + "-" + str(x[-1]))
    plt.savefig("results/HIC_FREQ/1599_cancer_chr"+chrs[chrn]+".png")"""
    """"""
    #matrix = Oligo.Matrix.HiCMatrix.read("data/HiC_Data/GSM1909121_observed/HiCtool_chr"+chrs[chrn]+"_40kb_observed.txt", resolution=bucSize)
    #bucFirst, bucLast = bucLim(mdsStart,mdsEnd,bucSize)
    data=np.array(matrix.complete_matrix())#.sub_matrix(3350,3450)
    #print(data)
    means=[]
    x=[]
    for i in range(len(data)):
        a=np.array(data[i][i-5:i+6])
        #b=np.array(data[:][i])
        vals=a#np.concatenate((a,b))
        mean=np.mean(vals)#np.sum(vals[:10])-np.sum(vals[10:])#np.mean(vals)
        means.append(mean)
        x.append(int(i)*40000/1e6)
    """plt.figure(figsize=(12,8))
    plt.xlabel("Coordinate of the bucket [Mbp]")
    plt.ylabel("Average HiC-frequency in the neighborhood")
    plt.plot(x, np.array(means),"-.") #3350-3450
    plt.grid(True)
    plt.title("GSM1909121 cell chr"+chrs[chrn]+": Average Freq in Buckets "+ str(x[0]) + "-" + str(x[-1]))
    plt.savefig("results/HIC_FREQ/9121_healthy_chr"+chrs[chrn]+".png")"""
    g=open("results/difflists/c"+chrs[chrn]+".diffbucket")
    difflist=[]      
    for ll in g:
        difflist.append(float(ll))
    diffsnp=np.array(difflist)
    #diffs = spectra_diffs(spectrum_c5,freq,optiseqs, species, mdsStart,mdsEnd,mdsC,bucSize,5,0.5,difflist,0)
    y=np.array(diffsnp[:-1])
    
    ff = open("results/HiC_cancerous_specDev_pqArm.correlation", "a")
    ff.write("CHROMOSOME: " + chrs[chrn] + "\n")
    print("CHROMOSOME: " + chrs[chrn])
    print(len(diffsnp))
    print(len(means))
    data={"HiC_freq": np.array(means[550:650]),#np.array(means[cenEnd[chrn]:]),
          "Deviation": np.array(y[550:650])}##np.array(y[cenEnd[chrn]:])}
    corrdata=pd.DataFrame(data,columns=['HiC_freq','Deviation'])
    corrMatrix = corrdata.corr()
    print("pArm")
    print (corrMatrix)
    corrVal=round(corrMatrix.loc["HiC_freq","Deviation"],3)
    bpCorr.append(corrVal)
    """ff.write("pArm"  + "\n")
    dfAsString = corrMatrix.to_string(header=True, index=True)
    ff.write(dfAsString)
    ff.write("\n")
    data={"HiC_freq": np.array(means[:cenStart[chrn]]),
          "Deviation": np.array(y[:cenStart[chrn]])}
    corrdata=pd.DataFrame(data,columns=['HiC_freq','Deviation'])
    corrMatrix = corrdata.corr()
    print("qArm")
    print (corrMatrix)
    ff.write("qArm"  + "\n")
    dfAsString = corrMatrix.to_string(header=True, index=True)
    ff.write(dfAsString)
    ff.write("\n")"""
    ff.close()
    
print(bpCorr)
print(np.mean(np.array(bpCorr)))
"""
"""
"""
g=open("results/difflists/c5.diffbucket")
difflist=[]      
for ll in g:
    difflist.append(float(ll))
diffsnp=np.array(difflist)
plt.figure(figsize=(15,8))
plt.plot(np.array(range(len(difflist)))*40000/1e6,difflist,".")#SPECTRUM DEVIATION PLOT
plt.xlabel("Starting coordinate of a bucket [Mbp]")
plt.ylabel("Average absolute deviation (k=5)")
plt.title("Deviation of bucket kmer spectra in chromosome 5")
plt.grid()
plt.legend(title="Average deviation: " + str(np.round(np.mean(diffsnp),4)))
plt.figure()
#plt.plot(np.array(means),".")#*np.array(means)[:-1])


plt.figure(figsize=(15,8))
plt.plot(np.array(range(len(difflist[3000:4000])))*40000/1e6,difflist[3000:4000],".", label="Deviation")#SPECTRUM DEVIATION PLOT
plt.xlabel("Starting coordinate of a bucket [Mbp]")
plt.ylabel("Average absolute deviation (k=5)")
plt.title("Deviation of bucket kmer spectra in chromosome 5")
plt.grid()
plt.legend(title="Average deviation: " + str(np.round(np.mean(diffsnp[3000:4000]),4)))
#plt.figure()
#plt.plot(np.array(means),".")#*np.array(means)[:-1])
data={"Avecorr": np.array(means),
      "Deviation": np.array(difflist[3000:4000])}
corrdata=pd.DataFrame(data,columns=['Avecorr','Deviation'])
corrMatrix = corrdata.corr()
print (corrMatrix)

# create figure and axis objects with subplots()

fig,ax = plt.subplots(figsize=(12,10))
ax.set_title("Relationship between the correlation and deviation of bucket spectra")
# make a plot
ax.set_ylim(0.5,1.2)
ax.invert_yaxis()
ax.plot(np.array(range(len(difflist[3000:4000])))*40000/1e6+120, means, "ro")
# set x-axis label
ax.set_xlabel("Starting coordinate [Mbp]",fontsize=14)
# set y-axis label
ax.set_ylabel("Average correlation",color="red",fontsize=14)
# twin object for two different y-axis on the sample plot
ax2=ax.twinx()
ax2.set_ylim(0.1,1.5)
# make a plot with different y-axis using second axis object
ax2.plot(np.array(range(len(difflist[3000:4000])))*40000/1e6+120, diffsnp[3000:4000],"bo")
ax2.set_ylabel("Deviation",color="blue",fontsize=14)
ax2.text(0.2,0.9,"Correlation: " + str(round(corrMatrix.at["Deviation","Avecorr"],4)), fontsize=18,
              horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)
plt.show()
# save the plot as a file
fig.savefig('correlation.jpg',
            format='jpeg',
            dpi=100,
            bbox_inches='tight')

"""
"""
data={"Mean Contact freq": np.array(means)[3000:-500],
      "Dev of Gen Code": diffsnp[3000:-501]}
corrdata=pd.DataFrame(data,columns=['Mean Contact freq','Dev of Gen Code'])
corrMatrix = corrdata.corr()
print (corrMatrix)


data={"Mean Contact freq": np.array(means)[:],
      "Dev of Gen Code": diffsnp[:-1]}
corrdata=pd.DataFrame(data,columns=['Mean Contact freq','Dev of Gen Code'])
corrMatrix = corrdata.corr()
print (corrMatrix)
"""