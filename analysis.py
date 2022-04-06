# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 11:32:32 2020

@author: Wisam
"""


import Oligo 

class Breakpoint(object):

    def __init__(self, name, species, chromosome, start, length=None, end=None, strand=None, band=None, partner=None):
        self.name = name
        self.species = species
        self.chromosome = chromosome
        if end is not None and end < start:
           start, end = end, start
        self.start = start
        self.strand = strand
        self.band = band
        self.seq = None
        if length is not None:
            self.length = length
        else:
            self.length = end-self.start
        if end is not None:
            self.end = end
        else: 
            self.end = self.start+length
    
    def get_name(self):
        return self.name

    def set_partner(self, partnerbp):
        self.partner = partnerbp.get_name()
        partnerbp.partner = self.get_name()
        
    def bp_bucket(self, bucketsize):
        bucstart = self.start/bucketsize
        bucend = self.end/bucketsize + 1
        return bucstart, bucend
    

    def bp_wrap(self, wrapsize, seqfile):
        bpwrap = Breakpoint(self.name, self.species, self.chromosome, self.start, self.length,  self.end, self.strand, self.band, self.partner)
        if self.start - wrapsize < 0:
            bpwrap.start = 0
        else:
            bpwrap.start = self.start - wrapsize
            
        genome = Oligo.File.read_genome(self.species, seqs_filename=Oligo.File.search(seqfile))
        for chromo in genome:
            if str(chromo) == self.species + " " + self.chromosome:
                l = chromo.length
        if self.end + wrapsize > l:
            bpwrap.end = l
        else:
            bpwrap.end = self.end + wrapsize
        return bpwrap
    
    def bp_seq(self, seqfile):
        if self.seq == None:
            print "Assigning sequence."
            genome = Oligo.File.read_genome(self.species, seqs_filename=Oligo.File.search(seqfile))
            for chromo in genome:
                if str(chromo) == self.species + " " + self.chromosome:
                    self.seq=chromo.get_seq()[self.start:self.end]
                    #print len(self.seq)
                    #print self.end-self.start
                    #print self.start
                    #print self.end
        else: 
            print "Breakpoint contains sequence."
        
    def bp_save(self, seqfile, filename=None):
        self.bp_seq(seqfile)
        if filename == None:
            filename= self.name+".fasta"
        with open(filename, "w+") as f:
            f.write(">" + self.name + "\n")
            f.write(self.seq)

def bucketlim(bucsize, buc, bucend=None, offset=0):
        if bucend==None:
            s = buc*bucsize + offset
            e = (buc+1)*bucsize + offset
        else:
            s=buc*bucsize + offset
            e=(bucend+1)*bucsize +offset
        return s,e
        
    
    
    
    ##
    ##
    ##
    
    
    
    
    
import Oligo
from source import Breakpoint as Bp

def read_chromo_from_genome(species, seqfile, chromoname):
    genome = Oligo.File.read_genome(species, seqs_filename=Oligo.File.search(seqfile))
    for chromo in genome:
        if str(chromo) == species + " " + chromoname:
            return chromo
    
def kmer_breakpoint(breakpoint, seqfile, outputfile, k):
    if breakpoint.seq == None:
        breakpoint.bp_seq(seqfile)
    seq = breakpoint.seq
    Oligo.Search.search_kmer(data_seq=seq, k=k, data_seq_name=str(breakpoint.name), output_filename=outputfile)
    
def kmer_bucket(chromo, outputfile, k, bucsize, buc, bucend = None,offset=0,GC=False):
    if GC==False:
        c = chromo
        s, e = Bp.bucketlim(bucsize, buc, bucend,offset)
        stop=False
        seq = c.get_seq()[s:e]
        if len(seq)==0:
            stop=True
        Oligo.Search.search_kmer(data_seq=seq, k=k, data_seq_name="Bucket" + str(buc), output_filename=outputfile)
    if GC==True:
        c = chromo
        s, e = Bp.bucketlim(bucsize, buc, bucend,offset)
        stop=False
        seq = c.get_seq()[s:e]
        stop = Oligo.Kmer.gc_content(seq)
    return stop

def chain_breakpoint(breakpoint, seqfile, outputfile, chain_seed, minp=1.0,minl=12):
    chromo = read_chromo_from_genome("Homo sapiens", seqfile, breakpoint.chromosome)
    seq = chromo.get_seq()
    #if breakpoint.seq == None:
    #    breakpoint.bp_seq(seqfile)
    #seq = breakpoint.seq
    Oligo.Search.search_chains(data_seq=seq, chain_seed_seq=chain_seed, data_seq_name=str(seq), chain_seed_seq_name=chain_seed, output_filename=outputfile , search_type=0, min_purity=minp, min_length=minl)
  
def chromo_kmer(k, species, seqfile, chromo, outputfile):
    c = read_chromo_from_genome(species, seqfile, chromo)
    c.load_seq()
    Oligo.Search.search_kmer(data_seq=c.get_seq(), k=k, data_seq_name=str(c), output_filename=outputfile % str(c))
    c.unload_seq()
    

##
    ##
    ##
    
def kmer_bucket_hm(name,start,end,species, chromoname, seqfile, kmerbucdir, specdir, matrixfile, drawerfile,k,bucsize,vmin=0.7,vmax=1., add_spectra=[], show_values=False, comb=False,
                   gbfile=None, exonstart=0, exonend=0,bpstart=0,CDS=False, bpbuc1=[], bpbuc2=[]):
    if comb==False:
        spectra=[] 
        gcs=[]
        gcx=[]
        cov1=[]
        cov2=[]
        spectra = spectra + add_spectra
        genome = Oligo.File.read_genome(species, seqs_filename=Oligo.File.search(seqfile))
        for chromo in genome:
            if str(chromo) == species + " " + chromoname:
                c = chromo
        exons = exonlist(gbfile, exonstart, exonend,CDS)
        for n in range(start,end):
            #print Bp.bucketlim(40000,n)
            cov1.append((bucketexon_coverage(n,exons,bucsize,bpstart),n*bucsize))
            stop = Km.kmer_bucket(c, kmerbucdir+str(bucsize)+name+ str(k)+ " " +str(n) +  ".kmer", k, bucsize, n)
            gc = Km.kmer_bucket(c, kmerbucdir+str(bucsize)+name+ str(k)+ " " +str(n) +  ".kmer", k, bucsize, n,GC=True)
            if stop==False:
                spectra.append(Oligo.Kmer.KmerSpectrum.read(kmerbucdir + str(bucsize)+name+str(k) + " " +str(n) +".kmer", name= chromoname[1]+" " +str(n*bucsize)+ " E%: " + str(round(cov1[-1][0]*100,1))))
                gc = Km.kmer_bucket(c, kmerbucdir+str(bucsize)+name+ str(k)+ " " +str(n) +  ".kmer", k, bucsize, n,GC=True)
                gcs.append(gc)
                gcx.append(n*bucsize)
        """for spectrum in spectra: 
            spectrum.basic_stats()
            kmers = sorted(spectrum.get_kmers())
            drawer = Oligo.Plot.HistoDrawer(bins=spectrum.get_bins(), label=str(spectrum))
            drawer.plot(specdir +  str(k)+ " " + str(spectrum)[19:-6] + '_kmer.png', xticklabels=kmers, xticks=xrange(len(kmers)), figsize=(20,7), xticklabels_rotation=90)"""
        # Correlate
        matrix = Oligo.Kmer.KmerSpectrum.correlate_spectra(spectra)
        # Save Matrix
        matrix.save(matrixfile)
        # Display Matrix
        drawer = Oligo.Plot.HeatmapDrawer(data=matrix.matrix, xticklabels=matrix.col_names, yticklabels=matrix.row_names, vmin=vmin, vmax=vmax, show_values=show_values)#, cmap='rainbow'
        if len(matrix.col_names)>150:
            drawer.plot(drawerfile, figsize=(72,60))
        elif len(matrix.col_names)>70:
            drawer.plot(drawerfile, figsize=(48,40))
        elif len(matrix.col_names)>30:
            drawer.plot(drawerfile, figsize=(36,30))
        else:
            drawer.plot(drawerfile, figsize=(18,15))
        
        drawer2 = plt.figure()
        drawerfileGC=drawerfile.split(".")
        dGC=drawerfileGC[0].split("/")
        dGCdir=""
        for i in dGC[:-2]:
            dGCdir=dGCdir+i+"/"
        plt.plot(gcx,gcs,ls="-.")
        plt.title(name)
        plt.show()
        drawer2.savefig(dGCdir+"GC/"+name+str(bucsize)+"_"+str(start)+"GC.png")
        
        with open(dGCdir+"GC/"+name+str(chromoname)+"_"+str(bucsize)+"_"+str(start)+"GC.txt","w+") as f:
            for item in gcx:
                f.write(str(item) + " ")
            f.write("\n")
            for item in gcs:
                f.write(str(item) + " ")
        
        
    if comb==True:
        spectra=[]
        gcs1=[]
        gcx1=[]
        gcs2=[]
        gcx2=[]
        spectra = spectra + add_spectra
        genome = Oligo.File.read_genome(species, seqs_filename=Oligo.File.search(seqfile))
        for chromo in genome:
            if str(chromo) == species + " " + chromoname[0]:
                c0 = chromo
            if str(chromo) == species + " " + chromoname[1]:
                c1 = chromo
        cov1= []
        exons = exonlist(gbfile[0], exonstart[0], exonend[0],CDS[0])
        for n in range(start[0],end[0]):
            #print Bp.bucketlim(40000,n)
            cov1.append((bucketexon_coverage(n,exons,bucsize,bpstart[0]),n*bucsize))
            stop = Km.kmer_bucket(c0, kmerbucdir+str(bucsize)+name+chromoname[0]+ " "+str(k)+ " " +str(n) +  ".kmer", k, bucsize, n)
            if stop==False:
                bpcheck=False
                if len(bpbuc1)>0:
                    for i in bpbuc1:
                        if i==n*bucsize:
                            spectra.append(Oligo.Kmer.KmerSpectrum.read(kmerbucdir +str(bucsize)+ name+chromoname[0]+ " "+str(k) + " " +str(n) +".kmer",
                                                                        name= "*" + chromoname[0]+" " +str(n*bucsize)+ " E%: " + str(round(cov1[-1][0]*100,1))))
                            gc = Km.kmer_bucket(c0, kmerbucdir+str(bucsize)+name+ str(k)+ " " +str(n) +  ".kmer", k, bucsize, n,GC=True)
                            gcs1.append(gc)
                            gcx1.append(n*bucsize)
                            bpcheck=True
                if bpcheck==False:
                    spectra.append(Oligo.Kmer.KmerSpectrum.read(kmerbucdir +str(bucsize)+ name+chromoname[0]+ " "+str(k) + " " +str(n) +".kmer",
                                                                name= chromoname[0]+" " +str(n*bucsize)+ " E%: " + str(round(cov1[-1][0]*100,1))))
                    gc = Km.kmer_bucket(c0, kmerbucdir+str(bucsize)+name+ str(k)+ " " +str(n) +  ".kmer", k, bucsize, n,GC=True)
                    gcs1.append(gc)
                    gcx1.append(n*bucsize)
        exons = exonlist(gbfile[1], exonstart[1], exonend[1],CDS[1])
        cov2=[]
        for n in range(start[1],end[1]):
            #print Bp.bucketlim(40000,n)
            cov2.append((bucketexon_coverage(n,exons,bucsize,bpstart[1]),n*bucsize))
            stop = Km.kmer_bucket(c1, kmerbucdir+str(bucsize)+name+chromoname[1]+ " "+str(k)+ " " +str(n) +  ".kmer", k, bucsize, n)
            if stop==False:
                bpcheck=False
                if len(bpbuc2)>0:
                    for i in bpbuc2:
                        if i==n*bucsize:
                            spectra.append(Oligo.Kmer.KmerSpectrum.read(kmerbucdir +str(bucsize)+ name+chromoname[1]+ " "+str(k) + " " +str(n) +".kmer",
                                                                        name= "*" + chromoname[1]+" " +str(n*bucsize) + " E%: " + str(round(cov2[-1][0]*100,1))))
                            bpcheck=True
                            gc = Km.kmer_bucket(c1, kmerbucdir+str(bucsize)+name+ str(k)+ " " +str(n) +  ".kmer", k, bucsize, n,GC=True)
                            gcs2.append(gc)
                            gcx2.append(n*bucsize)
                if bpcheck==False:
                    spectra.append(Oligo.Kmer.KmerSpectrum.read(kmerbucdir +str(bucsize)+ name+chromoname[1]+ " "+str(k) + " " +str(n) +".kmer",
                                                            name= chromoname[1]+" " +str(n*bucsize) + " E%: " + str(round(cov2[-1][0]*100,1))))
                    gc = Km.kmer_bucket(c1, kmerbucdir+str(bucsize)+name+ str(k)+ " " +str(n) +  ".kmer", k, bucsize, n,GC=True)
                    gcs2.append(gc)
                    gcx2.append(n*bucsize)
        """for spectrum in spectra: 
            spectrum.basic_stats()
            kmers = sorted(spectrum.get_kmers())
            drawer = Oligo.Plot.HistoDrawer(bins=spectrum.get_bins(), label=str(spectrum))
            drawer.plot(specdir +  str(k)+ " " + str(spectrum)[19:-6] + '_kmer.png', xticklabels=kmers, xticks=xrange(len(kmers)), figsize=(20,7), xticklabels_rotation=90)"""
        # Correlate
        matrix = Oligo.Kmer.KmerSpectrum.correlate_spectra(spectra)
        # Save Matrix
        matrix.save(matrixfile)
        # Display Matrix
        drawer = Oligo.Plot.HeatmapDrawer(data=matrix.matrix, xticklabels=matrix.col_names, yticklabels=matrix.row_names, vmin=vmin, vmax=vmax, show_values=show_values)#, cmap='rainbow'
        
        if len(matrix.col_names)>400:
            drawer.plot(drawerfile, figsize=(144,120))
        elif len(matrix.col_names)>150:
            drawer.plot(drawerfile, figsize=(84,70))
        elif len(matrix.col_names)>70:
            drawer.plot(drawerfile, figsize=(48,40))
        elif len(matrix.col_names)>30:
            drawer.plot(drawerfile, figsize=(36,30))
        else:
            drawer.plot(drawerfile, figsize=(24,20))
        
        drawer2 = plt.figure(figsize=(15,8))
        drawerfileGC=drawerfile.split(".")
        dGC=drawerfileGC[0].split("/")
        dGCdir=""
        for i in dGC[:-2]:
            dGCdir=dGCdir+i+"/"
        ax = drawer2.add_subplot(211)
        plt.plot(gcx1,gcs1,ls="-.",marker=".")
        ax.set_title(str(chromoname[0]))
        ax2 = drawer2.add_subplot(212)
        plt.plot(gcx2,gcs2,ls="-.",marker=".")
        ax2.set_title(str(chromoname[1]))
        plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
        plt.show()
        drawer2.savefig(dGCdir+"GC/"+name+"_"+str(bucsize)+"_"+str(start[0])+"_"+str(start[1])+"GC.png")

                
        with open(dGCdir+"GC/"+name+str(chromoname[0])+"_"+str(bucsize)+"_"+str(start[0])+"GC.txt","w+") as f:
            for item in gcx1:
                f.write(str(item) + " ")
            f.write("\n")
            for item in gcs1:
                f.write(str(item) + " ")
        with open(dGCdir+"GC/"+name+str(chromoname[1])+"_"+str(bucsize)+"_"+str(start[1])+"GC.txt","w+") as f:
            for item in gcx2:
                f.write(str(item) + " ")
            f.write("\n")
            for item in gcs2:
                f.write(str(item) + " ")
                

        
    return cov1,cov2

def read_breakpointGB(gbfile, bpname):
    GB = Oligo.File.read_sequence_file(gbfile)
    lociG = (Oligo.File.read_genes(gbfile), "lociG " + bpname)
    lociE = (Oligo.File.read_exons(gbfile), "lociE "+ bpname)
    lociI = (Oligo.File.read_introns(gbfile), "lociI "+ bpname)
    lociC = (Oligo.File.read_cds(gbfile), "lociC "+ bpname)
    reg=[lociG, lociE, lociI, lociC]
    regions=[]
    for i in reg:
        if len(i[0])>0:
            regions.append(i)
    return GB, regions

def bp_HiC(rawfile, outputdir, s,e, bp, wrapbuc):
    Hic.draw_HiC_Matrix(rawfile, outputdir + "/" + bp.name + "_HiC.png", substart= s, subend=e)
    Hic.draw_HiC_Matrix(rawfile, outputdir + "/" + bp.name + "_Wrap_HiC.png", substart= s-wrapbuc, subend=e+wrapbuc)    
    

#ablbcr=[abl,bcr]
def bp_kmer_hm(bps,bucsize,vmin=0.7,vmax=1.0,offsetmin=0.0,offsetmax=0.0,estart=0,eend=-1,cds=True):
    for b in bps:
        s,e=b.bp_bucket(bucsize)
        #bp_HiC("data/HiC_Data/HiCtool_chr"+str(b.chromosome[1:])+"_40kb_observed.txt","results/ABL_BCR/ABL", s,e,b,10)
        #b.bp_save('data/Seq_files/Homo_sapiens.seqs')
        #b.bp_seq('data/Seq_files/Homo_sapiens.seqs')
        bGB, bRegions = read_breakpointGB("data/Breakpoints/"+b.name+".gb",b.name)  
        for k in [3,4,5]: #k5 für vergleichbar
            spec=[]
            for a in bRegions:
                l = [""]
                for i in a[0]:
                    l[0]=str(l[0])+str(bGB[1][i.start:i.start+i.length])
                print len(l[0])
                Oligo.Search.search_kmer(data_seq=l[0], k=k, data_seq_name=str(a[1]), output_filename="results/ABL_BCR/ABL/"+str(bucsize)+b.name+str(k)+str(a[1])+".kmer")
                spectrum = Oligo.Kmer.KmerSpectrum.read("results/ABL_BCR/ABL/"+str(bucsize)+b.name+str(k)+str(a[1])+".kmer",name= str(a[1])+" " + str(len(l[0])))
                spec.append(spectrum)
                bucmp = 15
                if bucsize==5000:
                    bucmp = 120
                if k==5:
                    cov1,cov2 = kmer_bucket_hm(b.name,s-bucmp,e+bucmp,"Homo sapiens", "c"+str(b.chromosome[1:]), 'data/Seq_files/Homo_sapiens.seqs', "results/ABL_BCR/ABL/KMER_BUCKET/", 'results/ABL_BCR/ABL/KMER_BUCKET_SPECTRUM/Homo sapiens ',
                                               'results/ABL_BCR/ABL/'+str(bucsize)+b.name+str(k)+'_UMG.matrix', 'results/ABL_BCR/ABL/'+str(bucsize)+b.name+'_kmer'+str(k)+'_hm.png',k,bucsize,vmin-offsetmin,vmax-offsetmax, add_spectra=spec,comb=False,
                                               gbfile="data/Breakpoints/"+b.name+".gb",exonstart=estart,exonend=eend,bpstart=b.start,CDS=cds,
                                               bpbuc1=[],bpbuc2=[])
                    cov3,cov4 = kmer_bucket_hm(b.name,s,e,"Homo sapiens", "c"+str(b.chromosome[1:]), 'data/Seq_files/Homo_sapiens.seqs', "results/ABL_BCR/ABL/KMER_BUCKET/", 'results/ABL_BCR/ABL/KMER_BUCKET_SPECTRUM/Homo sapiens ',
                                               'results/ABL_BCR/ABL/'+str(bucsize)+b.name+str(k)+'.matrix', 'results/ABL_BCR/ABL/'+str(bucsize)+b.name+'_kmer'+str(k)+'_hmEXACT.png',k,bucsize,vmin-offsetmin,vmax-offsetmax, add_spectra=spec,comb=False,
                                               gbfile="data/Breakpoints/"+b.name+".gb",exonstart=estart,exonend=eend,bpstart=b.start,CDS=cds,
                                               bpbuc1=[],bpbuc2=[])
                else:
                    cov1,cov2 = kmer_bucket_hm(b.name,s-bucmp,e+bucmp,"Homo sapiens", "c"+str(b.chromosome[1:]), 'data/Seq_files/Homo_sapiens.seqs', "results/ABL_BCR/ABL/KMER_BUCKET/", 'results/ABL_BCR/ABL/KMER_BUCKET_SPECTRUM/Homo sapiens ',
                                               'results/ABL_BCR/ABL/'+str(bucsize)+b.name+str(k)+'_UMG.matrix', 'results/ABL_BCR/ABL/'+str(bucsize)+b.name+'_kmer'+str(k)+'_hm.png',k,bucsize,vmin,vmax, add_spectra=spec,comb=False,
                                               gbfile="data/Breakpoints/"+b.name+".gb",exonstart=estart,exonend=eend,bpstart=b.start,CDS=cds,
                                               bpbuc1=[],bpbuc2=[])
                    cov3,cov4 = kmer_bucket_hm(b.name,s,e,"Homo sapiens", "c"+str(b.chromosome[1:]), 'data/Seq_files/Homo_sapiens.seqs', "results/ABL_BCR/ABL/KMER_BUCKET/", 'results/ABL_BCR/ABL/KMER_BUCKET_SPECTRUM/Homo sapiens ',
                                               'results/ABL_BCR/ABL/'+str(bucsize)+b.name+str(k)+'.matrix', 'results/ABL_BCR/ABL/'+str(bucsize)+b.name+'_kmer'+str(k)+'_hmEXACT.png',k,bucsize,vmin,vmax, add_spectra=spec,comb=False,
                                               gbfile="data/Breakpoints/"+b.name+".gb",exonstart=estart,exonend=eend,bpstart=b.start,CDS=cds,
                                               bpbuc1=[],bpbuc2=[])

def comb_kmer_hm(ablbcr,bucsize,vmin=0.7,vmax=1.0,offsetmin=0.0,offsetmax=0.0,exonstart=0,exonend=0,CDS=False,bpbuc1=[],bpbuc2=[]):
    for k in [3,4,5]:  
        spec=[]
        s=[]
        e=[]
        for b in ablbcr:
            sb,eb=b.bp_bucket(bucsize)
            s.append(sb)
            e.append(eb)
            #b.bp_save('data/Seq_files/Homo_sapiens.seqs')
            #b.bp_seq('data/Seq_files/Homo_sapiens.seqs')
            bGB, bRegions = read_breakpointGB("data/Breakpoints/"+b.name+".gb", b.name)  
            for a in bRegions:
                l = [""]
                for i in a[0]:
                    l[0]=str(l[0])+str(bGB[1][i.start:i.start+i.length])
                print len(l[0])
                Oligo.Search.search_kmer(data_seq=l[0], k=k, data_seq_name=str(a[1]), output_filename="results/ABL_BCR/ABL/"+str(k)+str(a[1])+".kmer")
                spectrum = Oligo.Kmer.KmerSpectrum.read("results/ABL_BCR/ABL/"+str(k)+str(a[1])+".kmer",name= str(a[1]))#+" " + str(len(l[0]))
                spec.append(spectrum)
        bucmp = 15
        if bucsize==5000:
            bucmp = 120
        if k==5:
            cov1,cov2 = kmer_bucket_hm(ablbcr[0].name+ablbcr[1].name,[s[0]-bucmp,s[1]-bucmp],[e[0]+bucmp,e[1]+bucmp],"Homo sapiens", ["c"+str(ablbcr[0].chromosome[1:]),"c"+str(ablbcr[1].chromosome[1:])], 'data/Seq_files/Homo_sapiens.seqs', "results/ABL_BCR/ABL/KMER_BUCKET/", 'results/ABL_BCR/ABL/KMER_BUCKET_SPECTRUM/Homo sapiens ',
                           'results/ABL_BCR/ABL/'+str(bucsize)+ablbcr[0].name+ablbcr[1].name+str(k)+'.matrix', 'results/ABL_BCR/ABL/'+str(bucsize)+ablbcr[0].name+ablbcr[1].name+'_kmer'+str(k)+'_hm.png',k,bucsize,vmin-offsetmin,vmax-offsetmax, add_spectra=spec,comb=True,
                           gbfile=["data/Breakpoints/"+ablbcr[0].name+".gb","data/Breakpoints/"+ablbcr[1].name+".gb"],exonstart=exonstart,exonend=exonend,bpstart=[ablbcr[0].start,ablbcr[1].start],CDS=CDS,
                           bpbuc1=bpbuc1,bpbuc2=bpbuc2)
            cov3,cov4 = kmer_bucket_hm(ablbcr[0].name+ablbcr[1].name,s,e,"Homo sapiens", ["c"+str(ablbcr[0].chromosome[1:]),"c"+str(ablbcr[1].chromosome[1:])], 'data/Seq_files/Homo_sapiens.seqs', "results/ABL_BCR/ABL/KMER_BUCKET/", 'results/ABL_BCR/ABL/KMER_BUCKET_SPECTRUM/Homo sapiens ',
                           'results/ABL_BCR/ABL/'+str(bucsize)+ablbcr[0].name+ablbcr[1].name+str(k)+'EXACT.matrix', 'results/ABL_BCR/ABL/'+str(bucsize)+ablbcr[0].name+ablbcr[1].name+'_kmer'+str(k)+'_hmEXACT.png',k,bucsize,vmin-offsetmin,vmax-offsetmax, add_spectra=spec,show_values=False,comb=True,
                           gbfile=["data/Breakpoints/"+ablbcr[0].name+".gb","data/Breakpoints/"+ablbcr[1].name+".gb"],exonstart=exonstart,exonend=exonend,bpstart=[ablbcr[0].start,ablbcr[1].start],CDS=CDS,
                           bpbuc1=bpbuc1,bpbuc2=bpbuc2)
        else:
            cov1,cov2 = kmer_bucket_hm(ablbcr[0].name+ablbcr[1].name,[s[0]-bucmp,s[1]-bucmp],[e[0]+bucmp,e[1]+bucmp],"Homo sapiens", ["c"+str(ablbcr[0].chromosome[1:]),"c"+str(ablbcr[1].chromosome[1:])], 'data/Seq_files/Homo_sapiens.seqs', "results/ABL_BCR/ABL/KMER_BUCKET/", 'results/ABL_BCR/ABL/KMER_BUCKET_SPECTRUM/Homo sapiens ',
                           'results/ABL_BCR/ABL/'+str(bucsize)+ablbcr[0].name+ablbcr[1].name+str(k)+'.matrix', 'results/ABL_BCR/ABL/'+str(bucsize)+ablbcr[0].name+ablbcr[1].name+'_kmer'+str(k)+'_hm.png',k, bucsize, vmin, vmax,add_spectra=spec,comb=True,
                           gbfile=["data/Breakpoints/"+ablbcr[0].name+".gb","data/Breakpoints/"+ablbcr[1].name+".gb"],exonstart=exonstart,exonend=exonend,bpstart=[ablbcr[0].start,ablbcr[1].start],CDS=CDS,
                           bpbuc1=bpbuc1,bpbuc2=bpbuc2)
            cov3,cov4 = kmer_bucket_hm(ablbcr[0].name+ablbcr[1].name,s,e,"Homo sapiens", ["c"+str(ablbcr[0].chromosome[1:]),"c"+str(ablbcr[1].chromosome[1:])], 'data/Seq_files/Homo_sapiens.seqs', "results/ABL_BCR/ABL/KMER_BUCKET/", 'results/ABL_BCR/ABL/KMER_BUCKET_SPECTRUM/Homo sapiens ',
                           'results/ABL_BCR/ABL/'+str(bucsize)+ablbcr[0].name+ablbcr[1].name+str(k)+'EXACT.matrix', 'results/ABL_BCR/ABL/'+str(bucsize)+ablbcr[0].name+ablbcr[1].name+'_kmer'+str(k)+'_hmEXACT.png',k, bucsize, vmin,vmax,add_spectra=spec,show_values=False,comb=True,
                           gbfile=["data/Breakpoints/"+ablbcr[0].name+".gb","data/Breakpoints/"+ablbcr[1].name+".gb"],exonstart=exonstart,exonend=exonend,bpstart=[ablbcr[0].start,ablbcr[1].start],CDS=CDS,
                           bpbuc1=bpbuc1,bpbuc2=bpbuc2)
        
        with open('results/ABL_BCR/ABL/exoncov/'+str(k)+ablbcr[0].name+ablbcr[1].name+str(bucsize)+".exoncov", "w+") as f:
            f.write(ablbcr[0].name + ":\n")
            for i in cov3:
                f.write("Bucket at " +str(i[1]) + ": " + str(i[0]) + "%\n")
                        
            f.write(" \n")
            f.write(ablbcr[1].name + ":\n")
            for i in cov4:
                f.write("Bucket at " +str(i[1]) + ": " + str(i[0]) + "%\n")
                        
            f.write(" \n")
            f.write(ablbcr[0].name + " and surroundings:\n")
            for i in cov1:
                f.write("Bucket at " +str(i[1]) + ": " + str(i[0]) + "%\n")
                        
            f.write(" \n")
            f.write(ablbcr[1].name + " and surroundings:\n")
            for i in cov2:
                f.write("Bucket at " +str(i[1]) + ": " + str(i[0]) + "%\n")
                        
            f.write(" \n")
            
            
def exonlist(gbfile, start, end,CDS):
    if CDS==False:
        return Oligo.File.read_mRNA(gbfile)[start:end]
    if CDS==True:
        return Oligo.File.read_cds(gbfile)[start:]
def bucketexon_coverage(bucket, exonlist,bucsize,bpstart):
    cov =0.0
    for e in exonlist:
        end=  Oligo.Locus.get_end(e) + bpstart
        start = e.start + bpstart
        bucs= bucket*bucsize
        buce= (bucket+1)*bucsize
        if bucs<=start<=buce or bucs<=end<=buce:
            cov= cov+(min(float(buce),float(end))-max(float(bucs),float(start)))/bucsize
    return cov

abl = Bp.Breakpoint(name="ABL", species="Homo sapiens", chromosome="c9", start=130713881, length=None, end=130887675, strand=None, band = "q34,12", partner=None)
bcr = Bp.Breakpoint(name="BCR", species="Homo sapiens", chromosome="c22", start=23180509, length=None, end=23318037, strand=None, band = "q11.23", partner=None)
myc = Bp.Breakpoint(name="cMYC", species="Homo sapiens", chromosome="c8", start=127735434, length=None, end=127742951, strand=None, partner=None)
igh = Bp.Breakpoint(name="IGH", species="Homo sapiens", chromosome="c14", start=105586437, length=None, end=106879844, strand=None, partner=None)
igk = Bp.Breakpoint(name="IGK", species="Homo sapiens", chromosome="c2", start=88857361, length=None, end=90235368, strand=None, partner=None)
pml = Bp.Breakpoint(name="PML", species="Homo sapiens", chromosome="c15", start=73994673, length=None, end=74047827, strand=None, partner=None)
rara = Bp.Breakpoint(name="RARa", species="Homo sapiens", chromosome="c17", start=40309180, length=None, end=40357643, strand=None, partner=None)
ercc1 = Bp.Breakpoint(name="ERCC1",species="Homo sapiens", chromosome="c19", start=45407334, length=None, end=45451547, strand=None, partner=None) #als reguläre Region
hoxa = Bp.Breakpoint(name="HOXa",species="Homo sapiens", chromosome="c7", start=27092993, length=None, end=27200091, strand=None, partner=None)
bcr.set_partner(abl)


def chains():

    abl = Bp.Breakpoint(name="ABL", species="Homo sapiens", chromosome="c9", start=130713881, length=None, end=130887675, strand=None, band = "q34,12", partner=None)
    bcr = Bp.Breakpoint(name="BCR", species="Homo sapiens", chromosome="c22", start=23180509, length=None, end=23318037, strand=None, band = "q11.23", partner=None)
    myc = Bp.Breakpoint(name="cMYC", species="Homo sapiens", chromosome="c8", start=127735434, length=None, end=127742951, strand=None, partner=None)
    igh = Bp.Breakpoint(name="IGH", species="Homo sapiens", chromosome="c14", start=105586437, length=None, end=106879844, strand=None, partner=None)
    igk = Bp.Breakpoint(name="IGK", species="Homo sapiens", chromosome="c2", start=88857361, length=None, end=90235368, strand=None, partner=None)
    pml = Bp.Breakpoint(name="PML", species="Homo sapiens", chromosome="c15", start=73994673, length=None, end=74047827, strand=None, partner=None)
    rara = Bp.Breakpoint(name="RARa", species="Homo sapiens", chromosome="c17", start=40309180, length=None, end=40357643, strand=None, partner=None)
    ercc1 = Bp.Breakpoint(name="ERCC1",species="Homo sapiens", chromosome="c19", start=45407334, length=None, end=45451547, strand=None, partner=None) #als reguläre Region
    hoxa = Bp.Breakpoint(name="HOXa",species="Homo sapiens", chromosome="c7", start=27092993, length=None, end=27200091, strand=None, partner=None)
    bplist=[abl,bcr,myc,igh,igk,pml,rara,ercc1,hoxa]
    chainlist=["A","C","G","T","AC","AT","AG","GC","TG","TC"]
    for b in bplist:
        for c in chainlist:
            if os.path.isfile("results/ABL_BCR/chains/"+b.name+"_("+c+")n.loci") == False:
                Km.chain_breakpoint(b, 'data/Seq_files/Homo_sapiens.seqs', "results/ABL_BCR/chains/"+b.name+"_("+c+")n.loci", c, minp=0.9,minl=10)
            l1=Oligo.Locus.read("results/ABL_BCR/chains/"+b.name+"_("+c+")n.loci")
            
            ldist=[]
            for i in range(len(l1[1:-1])):
                if l1[i].start>b.start and l1[i].start<b.end:
                    lminus=abs(l1[i].start-l1[i-1].start+l1[i].length)
                    lplus=abs(l1[i].start+l1[i].length-l1[i+1].start)
                    ldist.append(min(lminus,lplus))
                    #print min(lminus,lplus)
            plt.figure() 
            plt.hist(ldist,bins=50,range=[0,(b.end-b.start)/5])
            plt.title(b.name+c+ " binsize: "+str((b.end-b.start)/250))
            plt.savefig("results/ABL_BCR/chains/distances/"+b.name+c+"_distances.png")
            #plt.figure()
            aa,bb,f1 = La.loci_binning("results/ABL_BCR/chains/"+b.name+"_("+c+")n.loci", 5000, b.start, b.end,c + "purity 0.9, > 10",b.name)
            with open("results/ABL_BCR/chains/"+b.name+str(5000)+c+".txt","w+") as f:
                for item in bb:
                    f.write(str(item) + " ")
                f.write("\n")
                for item in aa:
                    f.write(str(item) + " ")
            plt.savefig("results/ABL_BCR/chains/"+b.name+str(5000)+c+".png")
            aa,bb,f2 = La.loci_binning("results/ABL_BCR/chains/"+b.name+"_("+c+")n.loci", 40000, b.start, b.end,c+ "purity 0.9, > 10",b.name)
            plt.savefig("results/ABL_BCR/chains/"+b.name+str(40000)+c+".png")
            with open("results/ABL_BCR/chains/"+b.name+str(40000)+c+".txt","w+") as f:
                for item in bb:
                    f.write(str(item) + " ")
                f.write("\n")
                for item in aa:
                    f.write(str(item) + " ")
            plt.show()
            
 
def exoncovumg(bp, buc1, bucm1):
    lociABL = Oligo.File.read_mRNA("data/Genome/Homo sapiens_"+bp.chromosome+".gb")
    #lociABL = Oligo.File.read_cds("data/Genome/Homo sapiens_"+bp.chromosome+".gb")
    lociABL600 = []
    for l in lociABL:
        if buc1 < l.start < bucm1:
            lociABL600.append(l)

    buckets= buc1/5000
    buckete=bucm1/5000
    #Oligo.Loci.coverage_sum(lociABL600)
    lociABL600f = Oligo.Loci.remove_overlap(lociABL600)
    #print lociABL600f
    covumg1=[]
    for i in range(buckets,buckete+1):
        cov =0.0
        for e in lociABL600f:
            end=  Oligo.Locus.get_end(e) 
            start = e.start
            bucs= i*5000
            buce= (i+1)*5000
            #print end,start, bucs, buce
            if bucs<=start<=buce or bucs<=end<=buce:
                cov= cov+(min(float(buce),float(end))-max(float(bucs),float(start)))/5000
        covumg1.append(cov*100)
    #print covumg1
    return covumg1
