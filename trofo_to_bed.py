# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 22:53:01 2022

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

"""
chrs=["c1","c2","c3","c4","c5","c6","c7","c8","c9","c10","c11","c12","c13",
      "c14","c15","c16","c17","c18","c19","c20","c21","c22","cX","cY"]
for j in chrs:
    n = open("data/bp_loci/bp_"+j+".loci", "r")
    f = open("data/TCGA_ChiTaRS_combined_fusion_information_on_hg19.bed", "a")
    for i in n:
        if i[0:5]!= "#date" and i[0:5]!= "start":
            pos = i.split()[0]
            f.write("chr"+j[1:]+"\t"+pos+"\t"+str(int(pos)+1)+"\n")   
    n.close()
    f.close()
"""
chrs=["c1","c2","c3","c4","c5","c6","c7","c8","c9","c10","c11","c12","c13",
      "c14","c15","c16","c17","c18","c19","c20","c21","c22","cX","cY"]
for i in chrs:
    g=open("results/difflists/"+i+".diffbucket")
    difflist=[]      
    for ll in g:
        difflist.append(float(ll))
    diffsnp=np.array(difflist)
    plt.figure(figsize=(15,8))
    plt.plot(np.array(range(len(difflist)))*40000/1e6,difflist,".")#SPECTRUM DEVIATION PLOT
    plt.xlabel("Starting coordinate of a bucket [Mbp]")
    plt.ylabel("Average absolute deviation (k=5)")
    plt.title("Deviation of bucket kmer spectra in chromosome " + i[0:])
    plt.grid()
    plt.legend(title="Average deviation: " + str(np.round(np.mean(diffsnp),4)))
    plt.savefig("results/deviation_plot/deviationplot"+i+".png")