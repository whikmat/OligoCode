# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 21:19:54 2021

@author: Wisam
"""

import sys
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


def get_chromoseq(species, optiseq, cID,start,end=None,excludemiddle=False, middlestart=None,middleend=None):    
    chromo=read_chromo_from_genome(species, optiseq, cID)
    if excludemiddle==False:
        seq=chromo.get_seq()[start:end]
    elif excludemiddle==True:
        seq=chromo.get_seq()[start:middlestart]+chromo.get_seq()[middleend:end]
    return seq

def merge_loci_files(file1,file2,resultfile):
    with open(file1) as f:
        content = f.readlines()
    i=0
    while "position" not in content[i]:
        i=i+1
    
    lst1=[]
    for l in content[i+1:]:
        pos=l.split()[0]
        #ind=content.index(l)
        lst1.append([int(pos),l])
        
    with open(file2) as f:
        content = f.readlines()
    i=0
    while "position" not in content[i]:
        i=i+1
    
    lst2=[]
    for l in content[i+1:]:
        pos=l.split()[0]
        #ind=content.index(l)
        lst2.append([int(pos),l])
    
    lst=sorted(lst1+lst2)
    #print(lst)
    
    with open(resultfile,"w") as f:
        for l in content[:i+1]:
            f.write(l)
        for j in lst:
            f.write(j[1])
        
        
def get_distances(locifile,length):
    with open(locifile) as f:
        content = f.readlines()
    i=0
    while "position" not in content[i]:
        i=i+1
    
    lst=[]
    for l in content[i+1:]:
        pos=l.split()[0]
        #ind=content.index(l)
        lst.append(int(pos))
    
    dist=[]
    for j in range(len(lst)-1):
        dist.append(lst[j+1]-lst[j]-length)
    plt.figure(figsize=(18,9))
    plt.title(locifile.split("/")[-1].split(".")[0]+" monomer motif distances on both strands")
    plt.hist(dist,1000,(0,2000))
    plt.savefig(locifile.split(".")[0]+"_all.png")
    """plt.figure()
    plt.title(locifile.split(".")[0]+" monomer motif distances on both strands")
    plt.hist(dist,50,(0,200))
    plt.figure()
    plt.title(locifile.split(".")[0]+" monomer motif distances on both strands")
    plt.hist(dist,50,(0,50))"""
    return dist
    

def get_distancesOneStrand(locifile,length,strand=0): #strand 0 or 1
    with open(locifile) as f:
        content = f.readlines()
    i=0
    while "position" not in content[i]:
        i=i+1
    
    lst=[]
    for l in content[i+1:]:
        pos=l.split()[0]
        if l.split()[1] == str(strand):
            #ind=content.index(l)
            lst.append(int(pos))
    
    dist=[]
    for j in range(len(lst)-1):
        dist.append(lst[j+1]-lst[j]-length)
    plt.figure(figsize=(18,9))
    plt.title(locifile.split("/")[-1].split(".")[0]+" monomer motif distances on strand " + str(strand))
    plt.hist(dist,1000,(0,2000))
    plt.savefig(locifile.split(".")[0]+"_strand"+str(strand)+".png")
    """plt.figure()
    plt.title(locifile.split(".")[0]+" monomer motif distances on strand " + str(strand))
    plt.hist(dist,50,(0,200))
    plt.figure()
    plt.title(locifile.split(".")[0]+" monomer motif distances on strand " + str(strand))
    plt.hist(dist,50,(0,50))"""
    return dist

def get_distancesAlternatingStrand(locifile,length): #strand 0 or 1
    with open(locifile) as f:
        content = f.readlines()
    i=0
    while "position" not in content[i]:
        i=i+1
    
    lst=[]
    prevStrand=None
    for l in content[i+1:]:
        pos=l.split()[0]
        if prevStrand == None:
            prevStrand=l.split()[1]
            lst.append(int(pos))
        if  prevStrand== str(0):
            if l.split()[1]=="1":
                lst.append(int(pos))
                prevStrand=l.split()[1]
        if prevStrand == str(1):
            if l.split()[1]=="0":
                lst.append(int(pos))
                prevStrand=l.split()[1] 
    #print(lst)
    
    dist=[]
    for j in range(len(lst)-1):
        dist.append(lst[j+1]-lst[j]-length)
    plt.figure(figsize=(18,9))
    plt.title(locifile.split("/")[-1].split(".")[0]+" monomer motif distances on alternating strand")
    plt.hist(dist,1000,(0,2000))
    plt.savefig(locifile.split(".")[0]+"_alternating.png")
    """plt.figure()
    plt.title(locifile.split(".")[0]+" monomer motif distances on alternating strand")
    plt.hist(dist,50,(0,200))
    plt.figure()
    plt.title(locifile.split(".")[0]+" monomer motif distances on alternating strand")
    plt.hist(dist,50,(0,50))"""
    return dist
        

def rare_multimerge(dic,cid): #cid is chromoName
    rarelst=[dic+"rare1_"+cid+"_"+regData[0]+".loci",dic+"rare2_"+cid+"_"+regData[0]+".loci",dic+"rare3_"+cid+"_"+regData[0]+".loci",dic+"rare4_"+cid+"_"+regData[0]+".loci"]
    rarenew=[dic+"rare1_"+cid+"_"+regData[0]+".loci"]
    for i in range(len(rarelst)-1):
        if i !=len(rarelst)-2:
            merge_loci_files(rarenew[-1],rarelst[i+1],dic+"rareTEMP"+str(i)+"_"+cid+"_"+regData[0]+".loci")
            rarenew.append(dic+"rareTEMP"+str(i)+"_"+cid+"_"+regData[0]+".loci")
        else:
            merge_loci_files(rarenew[-1],rarelst[i+1],dic+"rare_"+cid+"_"+regData[0]+".loci")
#VARIABLES     

#mdsChromo = "chr5"
mdsStart=128000000
mdsEnd=150000000
species="Homo sapiens"
bucSize=40000
dataset="GSM1551599" #HiC dataset
optiseqs= "data/Seq_Files/Homo_sapiens_opti.seqs"
purity=0.9
chainlen=10
mdsC = "c5"
allc=[]
for i in range(1,22):
    allc.append("c"+str(i))
allc.append("cX")
allc.append("cY")
for i in allc:
    #regData=("MDS",int(128e6),int(150e6),False,None,None) 
    regData=("",0,None,False,None,None) 
    #regData=("exclCEN",0,None,True,int(46845e3),int(50059e3))          
    #seq=get_chromoseq(species, optiseqs, mdsC, mdsStart,mdsEnd)
    seq=get_chromoseq(species, optiseqs, i, regData[1],regData[2],regData[3],regData[4],regData[5])
    #RARE
    dic="results/RARE/"
    #reg="MDS"
    
    Oligo.Search.search_seq(data_seq=seq, query_seq="AGGTCA",output_filename=dic+"rare1_"+i+"_"+regData[0]+".loci") #then gapn then same seq again
    Oligo.Search.search_seq(data_seq=seq, query_seq="AGTTCA",output_filename=dic+"rare2_"+i+"_"+regData[0]+".loci")
    Oligo.Search.search_seq(data_seq=seq, query_seq="GGGTCA",output_filename=dic+"rare3_"+i+"_"+regData[0]+".loci")
    Oligo.Search.search_seq(data_seq=seq, query_seq="GGTTCA",output_filename=dic+"rare4_"+i+"_"+regData[0]+".loci")
    
      
    rare_multimerge(dic,i)
    d1=get_distances(dic+"rare_"+i+"_"+regData[0]+".loci",6)
    d2=get_distancesOneStrand(dic+"rare_"+i+"_"+regData[0]+".loci",6,0)
    d3=get_distancesOneStrand(dic+"rare_"+i+"_"+regData[0]+".loci",6,1)
    d4=get_distancesAlternatingStrand(dic+"rare_"+i+"_"+regData[0]+".loci",6)
    
    #2nd part of RARE if not the same again:
    """
    Oligo.Search.search_seq(data_seq=seq, query_seq="TGAACC",output_filename="rare5_c5.loci")
    Oligo.Search.search_seq(data_seq=seq, query_seq="TGACCC",output_filename="rare6_c5.loci")
    Oligo.Search.search_seq(data_seq=seq, query_seq="TGAACT",output_filename="rare7_c5.loci")
    Oligo.Search.search_seq(data_seq=seq, query_seq="TGACCT",output_filename="rare8_c5.loci")"""
    
