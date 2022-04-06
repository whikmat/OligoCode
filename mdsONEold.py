# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 11:58:19 2020

@author: Wisam
"""

import Oligo
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr 
import pandas as pd



#GETS CHROMO
def read_chromo_from_genome(species, seqfile, chromoname):
    genome = Oligo.File.read_genome(species, seqs_filename=Oligo.File.search(seqfile))
    for chromo in genome:
        if str(chromo) == species + " " + chromoname:
            return chromo # RETURNS A CHROMOSOME OBJECT

#GETS SEQUENCE OF A CHROMO REGION
def region_get_seq(start, end, chromo):
    seq=chromo.get_seq()[start:end]
    return seq #GET SEQUENCE OF A REGION

#DEFINES BUCKETS OF INTEREST DEPENDING ON THE REGION COORDINATES
def bucLim(start, end, bucsize):
    bucFirst=start//bucsize
    bucLast=(end-1)//bucsize + 1
    return int(bucFirst), int(bucLast) #RETURNS BUCKET RANGE OF A SUSPECTED REGION

#SPLITS MY REGION SEQUENCE IN BUCKETS EACH OF THEM CONTAINING THE CORRESPONDING SEQUENCE
def split_region_seq(start, end, chromoOfRegion, bucsize): #IMPORTANT: BEGIN FROM 1st BASE of 1st BUC 
    bucFirst, bucLast = bucLim(start, end, bucsize)
    startRegion=bucFirst*bucsize
    endRegion= bucLast*bucsize
    seqRegion = chromoOfRegion.get_seq()[startRegion:endRegion]
    seqBuckets=[]
    for i in range(0,bucLast-bucFirst):
        seqBuckets.append([bucFirst+i,seqRegion[i*bucsize:(i+1)*bucsize]])
    return seqBuckets #RETURNS A LIST OF SEQUENCES OF BUCKETS

#LOOKS FOR ALL KMER DISTRIBUTIONS IN A GIVEN LIST CONTAINING SEQUENCES OF BUCKETS
def kmer_data_creator(species, optiseq, cID,start,end,bucsize,k):
    chromo=read_chromo_from_genome(species, optiseq, cID)
    listseq= split_region_seq(start,end,chromo,bucsize)
    for i in range(len(listseq)):
        Oligo.Search.search_kmer(data_seq=listseq[i][1], k=k, data_seq_name="Bucket" + str(listseq[i][0]), output_filename="results/kmer"+str(k)+"/"+cID+"_"+str(bucsize)+"_" + str(listseq[i][0])+".kmer")

#READS THE SEQUENCES OF THE LIST DESCRIBED ABOVE
def kmer_read_spectra(optiseq,species,start,end,cID,bucsize,k):
    chromo=read_chromo_from_genome(species, optiseq, cID)
    listseq= split_region_seq(start,end,chromo,bucsize)
    spectra=[]
    for i in range(len(listseq)):
        spectrum = Oligo.Kmer.KmerSpectrum.read("results/kmer"+str(k)+"/"+cID+"_"+str(bucsize)+"_" + str(listseq[i][0])+".kmer",name= str(listseq[i][0]))
        spectra.append(spectrum)
    return spectra #RETURNS A LIST OF SPECTRA

    
#PLOTS THE HIC DATA AFTER DEFINING THE START AND END POINTS OF THE ANALYZED BUCKETS
def hic_plot(dataset, cID, name):
    fig=plt.figure(figsize=(36,33))
    axes1 = fig.add_axes([0.1,0.1,0.9,0.8])
    matrix = Oligo.Matrix.HiCMatrix.read("data/HiC_Data/"+dataset+"_observed/HiCtool_"+cID+"_40kb_observed.txt", resolution=bucSize)
    bucFirst, bucLast = bucLim(mdsStart,mdsEnd,bucSize)
    data=matrix.sub_matrix(bucFirst, bucLast).complete_matrix()
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
  
#SPECTRA DIFFERENCES CALCULATED AND DISPLAYED
def spectra_diffs(mainspectrum, mainlen, optiseq,species,start,end,cID,bucsize,k):
    chromo=read_chromo_from_genome(species, optiseq, cID)    
    spectra=kmer_read_spectra(optiseq,species,start,end,cID,bucsize,k)
    words=sorted(mainspectrum.data.keys())
    
    print len(words)
    for i in spectra:
        diff=0.0
        for j in words:
            mainval=float(mainspectrum.data.get(j))/mainlen
            if i.data.get(j)==None:
                continue
            val=float(i.data.get(j))/bucsize
            diff=diff+np.abs(mainval-val)
        if diff>0.5:
            #print i.name, diff
            wordoi=[]
            for w in words:
                mainval=float(mainspectrum.data.get(w))/mainlen
                if i.data.get(w)==None:
                    continue
                val=float(i.data.get(w))/bucsize
                #if (val/mainval -1)*100<-67 or (val/mainval -1)*100>200:
                wordoi.append([w,round((val/mainval -1)*100,2)]) #percentage relative to c5 ave
                    #print i.name, j, mainval, val
            #print i.name, wordoi 
            f=open("results/boi/c5_"+i.name+".kmerdiff","w")
            for word in  wordoi:
                #print word[1]
                mainval=float(mainspectrum.data.get(word[0]))/mainlen
                f.write(word[0]+": "+str(round(mainval,5))+"  " + str(word[1])+"% \n")
            f.close()

    
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
    plt.plot(bucs,y,".")
    plt.grid(True)
    plt.title("Number of "+ inputseed + "-Elements in Buckets " + str(bucs[0])+"-"+str(bucs[-1]))
    outliers=[inputseed]
    plt.savefig("results/figures/"+'chain_search_'+ str(bucsize)+inputseed+cID+"_"+str(bucs[0])+"-"+str(bucs[-1]) + '.png')
    for i in range(len(y)):
        if y[i]>cutoff:
            outliers.append(bucs[i])
    return y,outliers # LIST OF CHAIN COUNTS IN BUCKETS, LIST OF BUCKET IDS THAT ARE OUTLIERS

#READS THE REPEATMASKER DATA AND PLOTS IT    
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
    plt.plot(bucs,y,".")
    plt.grid(True)
    plt.title("Number of "+ el + "-Elements in Buckets " + str(bucs[0])+"-"+str(bucs[-1]))
    outliers=[el]
    plt.savefig("results/figures/"+chromo+'_'+ el + "_" + str(bucsize)+ '.png')
    for i in range(len(y)):
        if y[i]>cutoff:
            outliers.append(bucs[i])
    return y, outliers # LIST OF ELEMENT COUNTS IN BUCKETS, LIST OF BUCKET IDS THAT ARE OUTLIERS
    
#PLOTS GC CONTENT IN THE BUCKETS
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
    plt.plot(bucs,gcVals,".")
    plt.grid(True)
    plt.title("GC content in Buckets " + str(s)+"-"+str(e))
    plt.savefig("results/figures/GC_"+cID+"_"+str(bucsize)+"_"+str(s)+"-"+str(e)+ '.png')
    return gcVals  #LIST OF GC CONTENT VALUES OF BUCKETS

#CALCULATES THE AVERAGES CORRELATION VALUE OF A BUCKET BY AVERAGING THE VALUES OF THE CORRESPONDING MATRIX ROW AND COLUMN OF A BUCKET
#DATA STRUCTURE: [[bucketID, aveCorr],[bucketID2, aveCorr2],...]
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
    plt.plot(x, np.array(means),".")
    plt.grid(True)
    plt.title("Average Correlation in Buckets "+ str(matrix.col_names[0]) + "-" + str(matrix.col_names[-1]))
    res=[]
    for i in range(len(x)):
        res.append(np.array([x[i],means[i]]))
    outliers=["AveCorr:"]
    for i in range(len(means)):
        if means[i]<cutoff:
            outliers.append(x[i])
    
    return np.array(res),outliers #LIST OF [BUCKET ID, AVE CORR] AND OUTLIER BUCKET ID ABOVE A CERTAIN THRESHOLD

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
            
def wordgroupIncrease(mainspectrum, kmerdiff, lim, letter1, letter2=None, anti=False):    
    words=sorted(mainspectrum.data.keys())

    g= wordGrouping(sorted(words),lim,letter1,char2=letter2,anti=anti)
    f=open(kmerdiff)#"results/boi/c5_3394.kmerdiff")
    weight=0
    summa=0
    for line in f:
        items=line.split(" ")
        if items[0][:-1] in g:
            #print items
            weight=weight+float(items[1])*float(items[3][:-1])
            summa=summa+float(items[1])   
    ave=summa/len(g)
    aveweight=weight/summa    
    print g, len(g)
    print aveweight

    
#VARIABLES     
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

#KMER Bucket data, here in kmer5 folder
#kmer_data_creator(species, optiseqs, mdsC,mdsStart,mdsEnd,bucSize,5)

#makes a plot of the kmer spectra data as a correlation matrix above, matrix is stored
#kmer_plot(optiseqs, "MDS", species, mdsStart,mdsEnd,mdsC,bucSize,5,(12,10))
"""
#reads the above correlation matrix
mat = Oligo.Matrix.Matrix.read("results/c5_MDS.matrix")
#produces a list of average correlation of buckets [bucket, avecorr] with a list of outliers
res, oCorr = averageCorr(mat,0.7)
f = open('results/c5_MDS.avecorr','w') 
for i in res:
    f.write(str(int(i[0]))+ " " + str(round(i[1],4)) + "\n")
f.close()    
    


#bucket definition
s,e = bucLim(mdsStart,mdsEnd,bucSize)
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
dataAlu, oAlu = genomeElBucs(mdsC, "LOCI_alu_repeatmasker", "Alu",mdsStart,mdsEnd,bucSize,60)
dataL1, oL1 = genomeElBucs(mdsC, "LOCI_L1_hg38","L1",mdsStart,mdsEnd,bucSize,30)

#same but with GC% without outliers
dataGCp = gcRegion(mdsStart,mdsEnd,bucSize,species,optiseqs,mdsC)

data = {'Alu': dataAlu,
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
df = pd.DataFrame(data,columns=['Alu','L1','A','C','G','T','AC','AG','AT','CG','CT','GT','GC%'])
corrMatrix = df.corr()
print (corrMatrix)
csv_data = corrMatrix.to_csv('results/corr_dataset_'+mdsC+"_"+str(s)+"_"+str(e)+'.csv', index = True) 
"""
#KMER SPECTRUM OF C5 AS A WHOLE
#Oligo.Search.search_kmer(data_seq=read_chromo_from_genome(species, optiseqs, "c5").get_seq(), k=5, data_seq_name="c5", output_filename="results/c5.kmer")
spectrum_c5 = Oligo.Kmer.KmerSpectrum.read("results/c5.kmer",name= "c5")

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
#spectra_diffs(spectrum_c5,181538259,optiseqs, species, mdsStart,mdsEnd,mdsC,bucSize,5)


wordgroupIncrease(spectrum_c5,"results/boi/c5_3394.kmerdiff",5,"C",letter2="G",anti=False)















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
#sns.heatmap(mat.matrix,cmap="hot",vmin=0.5, vmax=1.0)
#plt.savefig("results/c5_MDS_2.png")