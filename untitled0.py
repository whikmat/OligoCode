# -*- coding: utf-8 -*-
"""
Created on Mon Jan 24 12:54:42 2022

@author: whikm
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

cenStart=np.array([122.0,92.1,90.7,49.7,46.4,58.5,58.1,44.0,43.2,39.6,51.0,34.7,
          16.0,16.0,17.0,36.3,22.8,15.4,24.4,26.4,10.8,12.9,58.6,10.3])
cenEnd=np.array([125.2,94.1,93.7,51.8,50.1,59.9,60.9,45.9,45.6,41.6,54.5,37.2,18.1,18.2,
        19.8,38.3,26.9,20.9,27.2,30.1,13.0,15.1,62.4,10.6])
cenStart=cenStart*1e6/40000
cenStart=cenStart.astype(int)
cenEnd=cenEnd*1e6/40000
cenEnd=cenEnd.astype(int)

###BP FUSION GDB###
bucLoc=[]
"""
t=open("data/TCGA_ChiTaRS_combined_fusion_information_on_hg38_CONVERTED.bed")


for l in t:
    bp = l.split("\t") #chr:5, coord:6 chr2:9, coord2:10
    #print(bp)
    bucLoc.append(("c"+bp[0][3:],int(bp[1]))) 
t.close()
"""
#OLD DATA
"""
bucLoc=[]
for l in t:
    bp = l.split("\t") #chr:5, coord:6 chr2:9, coord2:10
    #print(bp)
    bucLoc.append(("c"+bp[5][3:],int(bp[6]))) 
    bucLoc.append(("c"+bp[9][3:],int(bp[10]))) 
t.close()
"""
"""
chrs=["c1","c2","c3","c4","c5","c6","c7","c8","c9","c10","c11","c12","c13",
      "c14","c15","c16","c17","c18","c19","c20","c21","c22","cX","cY"]
for j in range(len(chrs)):
    n = open("data/bp_loci/bp_"+chrs[j]+".loci", "a")
    n.write("#date:	2021-10-24 20:38:59.263811\n")
    n.write("start\tlength\n")
    for i in bucLoc:
        if i[0]==chrs[j]:
            n.write(str(i[1])+ "\t1\n")
    n.close()

"""
cs1=["c1",248956422]
cs2=["c2",242193529]
cs3=["c3",198295559]
cs4=["c4",190214555]
cs5=["c5",181538259]
cs6=["c6",170805979]
cs7=["c7",159345973]
cs8=["c8",145138636]
cs9=["c9",138394717]
cs10=["c10",133797422]
cs11=["c11",135086622]
cs12=["c12",133275309]
cs13=["c13",114364328]
cs14=["c14",107043718]
cs15=["c15",101991189]
cs16=["c16",90338345]
cs17=["c17",83257441]
cs18=["c18",80373285]
cs19=["c19",58617616]
cs20=["c20",64444167]
cs21=["c21",46709983]
cs22=["c22",50818468]
csX=["cX",156040895]
csY=["cY",57227415]
chrs=[cs1,cs2,cs3,cs4,cs5,cs6,cs7,cs8,cs9,cs10,cs11,cs12,cs13,cs14,cs15,cs16,cs17,cs18,cs19,cs20,cs21,cs22,csX,csY]

f = open("results/boi_buckets.bucloc", "r")
buclocs=[]
for i in f:
    cID=i.split(" ")[0][1:]
    bucID=int(i.split(" ")[1])
    buclocs.append([cID,bucID])
f.close()
print(len(buclocs))
for i in buclocs:
    if i[0] == "X":
        cNo=22
    elif i[0] =="Y":
        cNo=23
    else:
        cNo=int(i[0])-1
        
    if i[1] < cenEnd[cNo] and i[1] > cenStart[cNo]:
        buclocs.remove(i)
   
       
print(len(buclocs))
"""
for j in range(len(chrs)):
    n = open("data/boi_loci/boi_"+chrs[j][0]+".loci", "a")
    n.write("#date:	2021-10-24 20:38:59.263811\n")
    n.write("start\tlength\n")
    for i in buclocs:
        #print(buclocs)
        if "c"+i[0]==chrs[j][0]:
            n.write(str(i[1]*40000)+ "\t40000\n")
    n.close()
"""


doc = open("results/summary.results","a")
for j in range(len(chrs)):
    print(chrs[j][0])
    doc.write(chrs[j][0]+"\n")
    loci1 = Oligo.Locus.read('data/bp_loci/bp_'+chrs[j][0]+'.loci')
    ###loci2 = Oligo.Locus.read('data/Lukas_loci/Final_pc/Loci/Final_pc Homo sapiens '+chrs[j][0]+'.loci')
    loci2 = Oligo.Locus.read('data/LOCI_L1_hg38/Homo sapiens '+chrs[j][0]+'_L1_repeatmasker.loci')
    #loci2 = Oligo.Locus.read("data/boi_loci/"+"boi_"+chrs[j][0]+".loci")
    #loci2 = Oligo.File.read_genes('E:/Daten/genbank/Homo sapiens c1.gb')
    results=Oligo.Loci.loci_in_loci(loci1, loci2, target_length=chrs[j][1])
    doc.write("bp in L1: \n")
    contained_perc1=int(round(100*float(results['k'])/results['n']))
    contained_perc2=int(round(100*float(results['c_sum'])/results['target_length']))
    doc.write(str(contained_perc1) + "% of bp in L1 which covers " + str(contained_perc2) + "%\n")
doc.close()



#print(buclocs)
"""
m = open("results/bp_buckets.bucloc")
bpbuclocs=[]
for i in m:
    cID=i.split(" ")[0][1:]
    bucID=int(i.split(" ")[1])
    bpbuclocs.append([cID,bucID])
m.close()
count=0
for i in bpbuclocs:
    if i in buclocs:
        count=count+1
print("BP")
print(count)
print(count/len(bpbuclocs))

m = open("results/b_alu_buckets.bucloc")
alubuclocs=[]
for i in m:
    cID=i.split(" ")[0][1:]
    bucID=int(i.split(" ")[1])
    alubuclocs.append([cID,bucID])
m.close()
count=0
for i in alubuclocs:
    if i in buclocs:
        count=count+1
print("ALU")
print(count)
print(count/len(alubuclocs))


m = open("results/b_l1_buckets.bucloc")
l1buclocs=[]
for i in m:
    cID=i.split(" ")[0][1:]
    bucID=int(i.split(" ")[1])
    l1buclocs.append([cID,bucID])
m.close()
count=0
for i in l1buclocs:
    if i in buclocs:
        count=count+1
print("L1")
print(count)
print(count/len(l1buclocs))
"""

###L1 HG38###
"""
chrs=["c1","c2","c3","c4","c5","c6","c7","c8","c9","c10","c11","c12","c13",
      "c14","c15","c16","c17","c18","c19","c20","c21","c22","cX","cY"]
for c in range(len(chrs)):
    t=open("data/LOCI_L1_hg38/Homo sapiens "+chrs[c]+"_L1_repeatmasker.loci")

    bucLoc=[]
    for l in t:
        loci = l.split("\t") #chr:5, coord:6 chr2:9, coord2:10
        #print(bp)
        #print(loci[0])
        if loci[0][0] in "0123456789":
            if int(int(loci[0])/40000) > cenEnd[c] or int(int(loci[0])/40000) < cenStart[c]:
                bucLoc.append(str(int(int(loci[0])/40000))+"\n")
                #print(int(int(loci[0])/40000))
    t.close()
    bpCorr=[]
    #print(bucLoc)
    n = open("results/b_l1_buckets.bucloc", "a")
    for i in bucLoc:
        n.write (chrs[c] + " " + str(i))
    n.close()
""" 