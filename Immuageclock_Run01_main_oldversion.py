## usage:call index:python ~/mao/Codeman/Immuagepipeline/0309testrun.py -f test -g   genome.fa  -s Star_CE --gtf /home/maolp/mao/Ref/Caenorhabditis_elegans/UCSC/ce10/Annotation/Genes/genes.gtf  -m ID
## usage: python ~/mao/Codeman/Immuagepipeline/0309testrun.py -f A02.Clean -a bowtie2 -m A -s /home/maolp/mao/Ref/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome -e 40

## the packages that you need to install  is :bowtie2,star,htseq-count,featureCounts,resm,samtools,parafly,multiqc,R,fastqc,trim_galore,python=3.7 above
## index you can download from:https://support.illumina.com/sequencing/sequencing_software/igenome.html




from doctest import script_from_examples
from email.policy import default
from ensurepip import version

from optparse import OptionParser
import os
import subprocess
import re
from socketserver import ThreadingUnixStreamServer
import sys
# from tkinter import E
# from typing_extensions import Self
# from unittest.util import three_way_cmp
# import pandas as pd
# import numpy as np
import logging
import datetime
import glob
import time
import configparser
import Immuageclock_Run05_TCRandHLA as IHT
from multiprocessing import Pool
import subprocess

# from tensorboard import notebook
#import matplotlib.pyplot as plt
#import seaborn as sns

version="1.0"

datename=time.strftime('%Y%m_%d%H%M',time.localtime(time.time()))

parser = OptionParser(usage="usage: %prog -f readFolde [-p <int>] [-e <int>] ", version="%prog 1")
parser.add_option("-f",
        action="store", 
        dest="readFolder",
        help="Folder with RNA-seq data (gzipped fastq files)")

parser.add_option("-e", "--threads",
                    action="store",
                    dest="threads",
                    default="8",
                    help="Parafly used. Default : 8")

parser.add_option("-p", "--Parallel",
                    action="store",
                    dest="Parallel",
                    default="True",
                    help="Parafly for parallel processing. Default: True")
parser.add_option("-s","--starindex",  
                    action="store",
                    dest="starindex",
                    default="/home/maolp/mao/Ref/hg38/genome",
                    help="STAR index.Default: ' /home/maolp/mao/Ref/hg38/genome|/home/maolp/mao/Ref/Homo_sapiens/UCSC/hg38/Annotation/star_hg38'")
parser.add_option("-m","--Module",
                    action="store",
                    dest="Module",
                    default="RA",
                    help="Many Module for you can choose about: \
                    S1(fastQC),\
                        \
                    S2(TrimQC),\
                        \
                    S3_1(star index),\
                        \
                    S3(align),S14(from QC to Count). \
                        \
                    S34(align to count),\
                    S5_1(rmarkdown),\
                        \
                    S5(Normalization),\
                        \
                    S6(Regression),\
                        \
                    S78(cirHLATCR)\
                        \
                    S7(TCRHLA),\
                        \
                    S8cir(circleRNA)\
                        \
                    RAR(runallreal)\
                        \
                    Default: RA(runall)")                                      
parser.add_option("-g","--genome",
                    action="store",
                    dest="genome",
                    default="/home/maolp/mao/Ref/",
                    help="Genome fasta file. Default: '/home/maolp/mao/Ref/'")

parser.add_option("-u","--gtf",
                    action="store",
                    dest="gtf",
                    default="/home/maolp/mao/Ref/Homo_sapiens/UCSC/hg38/Annotation/Genes.gencode/genes.gtf",
                    help="GTF file. Default: '/home/maolp/mao/Ref/Homo_sapiens/UCSC/hg38/Annotation/Genes.gencode/genes.gtf'")

parser.add_option("-c","--count",
                    action="store",
                    dest="count",
                    default="featureCounts",
                    help="featureCounts or rsem htseq-count. Default: featureCounts")

parser.add_option("-a","--alignmethods",
                    action="store",
                    dest="alignmethods",
                    default="hisat2",
                    help="Align methods:hisat2, bowtie2 or STAR. Default: hisat2")
parser.add_option("-t","--alignthreads",
                    action="store",         
                    dest="alignthreads",
                    default="40",
                    help="Bowtie option: Launch <int> parallel search threads. Default : 40")
parser.add_option("-x","--rsemindex",
                    action="store",
                    dest="rsemindex",
                    default="/home/maolp/mao/Ref/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/HG38RSEM/HG38RSEM", 
                    help="starindex that use by rsem. Default: '/home/maolp/mao/Ref/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/HG38RSEM'")

parser.add_option("-l","--less",
                    action="store",
                    dest="less",    
                    default="True",
                    help="only need a few file ,others can be move . Default: True")

parser.add_option("-n","--globname",
                    action="store",
                    dest="globname",
                    default="NEW",
                    help="globname. Default: ''")

parser.add_option("-w","--rmarkdown",
                    action="store",
                    dest="rmarkdown",
                    default="/home/maolp/mao/Codeman/Immuagepipeline/Immuageclock_Run02_PCA.Rmd",
                    help="ramrkdown file . Default: '/home/maolp/mao/Codeman/Immuagepipeline/Immuclock_trandata.Rmd'")

parser.add_option("-r","--rcut",
                    action="store",
                    dest="rcut",
                    default="/home/maolp/mao/Codeman/Immuagepipeline/Immuage_Rmain.R",
                    help="Rscript file. Default: '/home/maolp/mao/Codeman/Immuagepipeline/Immuage_Rmain.R'")

parser.add_option("--coldata",
                    action="store",
                    dest="coldata",
                    help="all.count.txt.No Default")

parser.add_option("--add",
                    action="store",
                    dest="add",
                    default="ap",
                    help="add all or add part or add not Default: 'ap'")

parser.add_option("--health",
                    action="store",
                    dest="helath",
                    default="False",    
                    help="update the old program. Default: False") 

parser.add_option("--rootdir",
                    action="store", 
                    dest="rootdir",
                    default="/data4/maolp/All_myself_data4/Immu_age/Immuageroot",
                    help="rootdir. Default: '/data4/maolp/All_myself_data4/Immu_age/Immuageroot'")

parser.add_option("--pipepath",
                    action="store",
                    dest="pipepath",
                    default="/home/maolp/mao/Codeman/All_Archived_Project/Immuagepipeline",
                    help="pipeline dir.Default: /home/maolp/mao/Codeman/All_Archived_Project/Immuagepipeline" )
parser.add_option("--resmode",
                    action="store", 
                    dest="resmode",
                    default="train_pre",
                    help="resmode. Default: train,test,train_pre")
parser.add_option("--config",
                    action="store",
                    dest="config",
                    help="config file. Default: ''")






def Fastqc(fold,fqcthereds=8):
    print("------------------------STEP1:QCreport-----------------------------")
    import datetime
    starttime = datetime.datetime.now()
    print("Start time: %s" % starttime)
    path=os.getcwd()
    os.system("mkdir -p "+path+"/"+globname+"_A01.raw.fastqc_report")
    os.system("ls "+fold+"/*gz |xargs -I [] echo ' fastqc -t 15 [] -o "+path+"/"+globname+"_A01.raw.fastqc_report'> "+globname+"_A01.fastqc.sh")
    os.system("ParaFly -c "+globname+"_A01.fastqc.sh -CPU "+str(fqcthereds)+"-failed_cmds "+globname+"_A01.fastqc.failed -v;rm -rf "+globname+"_A01.fastqc.sh")

    os.system("multiqc ./"+globname+"_A01.raw.fastqc_report")
    endtime = datetime.datetime.now()
    print("End time: %s" % endtime)
    print("Total time: %s" % (endtime - starttime))
    print("------------------------STEP1:QCreport END-----------------------------")

class cirRNA():
    def __init__(self,starfold,globname,threads=40):
        self.starfold=starfold
        self.globname=globname
        self.theread=threads
        self.aloutdirname=self.globname+"_A08.CiR_align"


    # def hello(self):
    #    print("hellow")     

    def circexplore(self,samplename):
        print("------------------------STEP2:circexplore parse -----------------------------")
        outdir_cicr_out=self.aloutdirname+"_circ"
        outdir=self.aloutdirname
        os.system("mkdir -p "+  outdir_cicr_out)
        os.system("CIRCexplorer2 parse -t STAR "+outdir+"/"+ samplename+"_Chimeric.out.junction  -b"+   outdir_cicr_out+"/"+ samplename+"_back_spliced_junction.bed")
        print(ref_cir)
        # print(outdir+"/"+ samplename+"_back_spliced_junction.bed" )
        # print(outdir_cicr_out+"/"+ samplename+"_back_spliced_junction.bed")
        os.system("CIRCexplorer2 annotate -r "+ ref_cir +" -g "+genomefa+ " -b "+  outdir_cicr_out+"/"+ samplename+"_back_spliced_junction.bed"  +" -o " + outdir_cicr_out+"/"+ samplename+"_cirRNA_known.txt" )
        os.system("rm -rf "+outdir+"/"+"*sam")
        
    def ParrunAlign_circ(self):
        print("------------------------STEP1:Alignment_circ START-----------------------------")
        therads=self.theread
        starfold=self.starfold

        
        for line in glob.glob(starfold+"/*_1.f*"):
            samplename=line.split("/")[-1].split("_1")[0]
            sample_realname=line.split("/")[-1]
            sample2=sample_realname.replace("_1","_2")
            print(samplename)
            
            print(starfold+"/"+sample2)
            if os.path.exists(starfold+"/"+sample2):
                print("Found paired-end reads")
                read1=os.path.join(starfold,sample_realname)
                read2=os.path.join(starfold,sample2)
                print(read1)
                print(read2)
                outdir= self.aloutdirname
                
                Align(readFile1=read1,readFile2=read2,threads=therads,samename=samplename,index=statrindex2,alignoutname=outdir)
                # self.hello()
                self.circexplore(samplename=samplename)

            else:
                pass
        
            # circexplore(self.starfold,self.globname)
        print("------------------------STEP2:circexplore END-----------------------------")
            
def thellow():
    print("hellow")

class TrimQCclass:
    def __init__(self, readFolder, readFile1,readFile2, fqcthereds=8):
        self.fq1=os.path.join(readFolder,readFile1)
        self.fq2=os.path.join(readFolder,readFile2)
        self.fqcthereds=fqcthereds

    def TrimQC(self):
        FQ1=self.fq1
        FQ2=self.fq2
        fqcthereds=self.fqcthereds
        
        print("------------------------STEP2:TrimQC-----------------------------")
        
        path=os.getcwd()
        

        os.system("echo \"trim_galore -q 25  --phred33 --length 36 -e 0.1 --stringency 3 --paired  "+FQ1+" "+FQ2+"  -o "+path+"/"+globname+"_A02.Clean \" >> "+globname+"_A02.trim_galore.sh")
    
    def ParrunTrimQC(self):
        print("------------------------STEP2:TrimQC_Parafly START-----------------------------")
        os.system("ParaFly -c "+globname+"_A02.trim_galore.sh -CPU "+str(self.fqcthereds)+"-failed_cmds "+globname+"_A02.trim_galore.failed -v;")
    
  
        Fastqc(fold=globname+""+globname+"_A02.Clean")
        os.system("mv "+globname+"_A02.Clean/*txt"+globname+"_A02.Clean.fastqc_report")
        print("------------------------STEP2:TrimQC_Parafly END-----------------------------")
        
        os.system("mkdir "+globname+"_A02.Clean.fastqc_report")
        
        os.system("mv "+globname+"_A02.Clean/*txt "+globname+"_A02.Clean.fastqc_report")
        os.system("rm -rf mulit*")
        os.system("multiqc "+globname+"_A02.Clean.fastqc_report/.")
        os.system("mv *html "+globname+"_A02.Clean.fastqc_report")



def Align(readFile1,readFile2,threads=40,Parallel="False",samename="Atest",alignoutname="",index=""):
    print("------------------------STEP3:Aliging-----------------------------")
    alignstarttime = datetime.datetime.now()
    print("Start time: %s" % alignstarttime)

    Al_fq1=readFile1
    Al_fq2=readFile2

    os.system("mkdir -p  "+globname+"_A03.Alignment")

    if options.alignmethods=="bowtie2":
        os.system("bowtie2  -p "+str(threads)+" -x "+options.starindex+" -1 "+Al_fq1+" -2 "+Al_fq2+" -S  "+globname+"_A03.Alignment/"+samename+".sam")
        os.system("samtools view -bS  "+globname+"_A03.Alignment/"+samename+".sam >  "+globname+"_A03.Alignment/"+samename+".bam")
        os.system("samtools sort -@ 40  "+globname+"_A03.Alignment/"+samename+".bam  -o  "+globname+"_A03.Alignment/"+samename+".sorted.bam")
        os.system("samtools index  "+globname+"_A03.Alignment/"+samename+".sorted.bam")
        os.system("rm  "+globname+"_A03.Alignment/"+samename+".sam")
        os.system("rm  "+globname+"_A03.Alignment/"+samename+".bam")


    elif options.alignmethods=="STAR":
        os.system("STAR  --runThreadN "+str(threads)+" \
        --runMode alignReads \
        --readFilesCommand zcat --quantMode TranscriptomeSAM GeneCounts \
        --twopassMode Basic \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped None \
        --genomeDir "+ starindex + " \
        --readFilesIn  " + Al_fq1 + " "+ Al_fq2 +" --outFileNamePrefix  "+globname+"_A03.Alignment/"+samename+"_")
        #os.system("samtools sort -@ 40  "+globname+"_A03.Alignment/"+samename+"_Aligned.out.bam  -o  "+globname+"_A03.Alignment/"+samename+"_Aligned.out.sorted.bam")
        
        # os.system("samtools sort -@ "+str(threads)+" -o  "+globname+"_A03.Alignment/"+samename+"_Aligned.sorted.bam  "+globname+"_A03.Alignment/"+samename+"_Aligned.out.bam")
        # os.system("samtools index  "+globname+"_A03.Alignment/"+samename+"_Aligned.out.sorted.bam")
        # os.system("rm -rf   "+globname+"_A03.Alignment/"+samename+"_Aligned.out.bam")


    elif options.alignmethods=="STARfast":
        os.system("STAR  --runThreadN "+str(threads)+" \
        --runMode alignReads \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --genomeDir "+ starindex + " \
        --readFilesIn  " + Al_fq1 + " "+ Al_fq2 +" --outFileNamePrefix  "+globname+"_A03.Alignment/"+samename+"_")
        #os.system("samtools sort -@ 40  "+globname+"_A03.Alignment/"+samename+"_Aligned.out.bam  -o  "+globname+"_A03.Alignment/"+samename+"_Aligned.out.sorted.bam")
        
        os.system("samtools sort -@ "+str(threads)+" -o  "+globname+"_A03.Alignment/"+samename+"_Aligned.sorted.bam  "+globname+"_A03.Alignment/"+samename+"_Aligned.out.bam")
        os.system("samtools index  "+globname+"_A03.Alignment/"+samename+"_Aligned.out.sorted.bam")
        os.system("rm -rf   "+globname+"_A03.Alignment/"+samename+"_Aligned.out.bam")

    elif options.alignmethods=="hisat2":
        #datawrite=open(globname+"_A03.Alignment/"+samename+".sh","a+")
        #seq="hisat2 -p "+str(options.alignthreads)+" -x "+options.starindex+" -1 "+Al_fq1+" -2 "+Al_fq2+" -S  "+globname+"_A03.Alignment/"+samename+".sam"
        #datawrite.write(seq)
        os.system("hisat2 -p "+str(options.alignthreads)+" -x "+options.starindex+" -1 "+Al_fq1+" -2 "+Al_fq2+" -S  "+globname+"_A03.Alignment/"+samename+".sam")
       # os.system("samtools view -bS  "+globname+"_A03.Alignment/"+samename+".sam >  "+globname+"_A03.Alignment/"+samename+".bam")
        os.system("/home/maolp/mao/Biosoft/samtools/samtools-1.9/samtools sort -@  "+str(options.alignthreads)+" "+globname+"_A03.Alignment/"+samename+".sam  -o  "+globname+"_A03.Alignment/"+samename+".sorted.bam")
        os.system("/home/maolp/mao/Biosoft/samtools/samtools-1.9/samtools index  "+globname+"_A03.Alignment/"+samename+".sorted.bam")
        os.system("rm  "+globname+"_A03.Alignment/"+samename+".sam")
        # os.system("rm  "+globname+"_A03.Alignment/"+samename+".bam")
    
    elif options.alignmethods=="salmon":
        os.system("/home/maolp/mao/Biosoft/salmon-1.8.0_linux_x86_64/bin/salmon quant -p "+str(options.alignthreads)+" -i "+options.starindex+" -l A -1 "+Al_fq1+" -2 "+Al_fq2+" -o  "+globname+"_A03.Alignment/"+samename)
        #os.system("hisat2 -p "+str(options.alignthreads)+" -x "+options.starindex+" -1 "+Al_fq1+" -2 "+Al_fq2+" -S  "+globname+"_A03.Alignment/"+samename+".sam")

    elif options.alignmethods=="STARcir":
        os.system("STAR  --runThreadN "+str(threads)+" \
        --readFilesCommand zcat \
        --chimSegmentMin 10 \
        --chimOutType Junctions \
        --genomeDir "+ index + " \
        --readFilesIn  " + Al_fq1 + " "+ Al_fq2 +" --outFileNamePrefix  "+alignoutname+"/"+samename+"_")
        

    print("------------------------STEP3:Alignment END-----------------------------")
    alignendtime = datetime.datetime.now()
    print("End time: %s" % alignendtime)
    print("Alignment time:",alignendtime-alignstarttime)
    
    print("\n")


def Count(readFolder="A03.Alignment",threads=40,Parallel="False",samename="All"):
    datetime.datetime.now()
    os.system("mkdir -p  "+globname+"_A04.Count")
    print("------------------------STEP4:Count-----------------------------")
    if options.count=="htseq-count":

        os.system("htseq-count -f bam -r pos -s no -t exon -i gene_id -a  "+globname+"_A04.Count/"+options.globname+samename+"_count.txt "+readFolder+"/"+samename+"_Aligned.sorted.bam "+gtf)
    elif options.count=="featureCounts":
        
        os.system("featureCounts -t exon  -p -g gene_id -T "+str(options.alignthreads)+" -a "+gtf+" -o  "+globname+"_A04.Count/"+str(options.globname)+samename+"_count.txt " +readFolder+"/*.bam")
    elif options.count=="rsem":
        os.system("mkdir -p  "+globname+"_A04.Count")
        for samp in glob.glob(readFolder+"/Run02.AlltoTranscriptomeBAM/*.bam"):
            newname=samp.split("/")[-1].split(".")[0]
            #os.system("resm  "+globname+"04.Count/"+samename+"_count.txt "+i+" "+gtf)
            os.system("rsem-calculate-expression --paired-end -no-bam-output --alignments -p "+str(threads)+" " + samp+" "+options.rsemindex+"  "+globname+"_A04.Count/"+newname+"_RSEM")
    
    print("------------------------STEP4:Count END-----------------------------")


def ParrunAlign(starfold,threads=40,Parallel="False"):
    print("------------------------STEP3:Alignment_Parafly START-----------------------------")

    
    for line in glob.glob(starfold+"/*_1.fq.gz"):
        samplename=line.split("/")[-1].split("_1_val_1.fq.gz")[0]
        print(samplename)
        print(starfold+"/"+samplename+"_2_val_2.fq.gz")
        if os.path.exists(starfold+"/"+samplename+"_2_val_2.fq.gz"):
            print("Found paired-end reads")
            read1=os.path.join(starfold,samplename+"_1_val_1.fq.gz")
            read2=os.path.join(starfold,samplename+"_2_val_2.fq.gz")
            print(read1)
            print(read2)
            
            Align(readFile1=read1,readFile2=read2,threads=options.alignthreads,Parallel=Parallel,samename=samplename)
    if options.alignmethods=="STAR":
        os.system("mkdir -p  "+globname+"_A03.Alignment/Aun01.AllalignedBAM")
        os.system("mv  "+globname+"_A03.Alignment/*Aligned.out.sorted.bam*  "+globname+"_A03.Alignment/Aun01.AllalignedBAM")
        os.system("mkdir -p  "+globname+"_A03.Alignment/Run02.AlltoTranscriptomeBAM")
        os.system("mv  "+globname+"_A03.Alignment/*toTranscriptome.out.bam  "+globname+"_A03.Alignment/Run02.AlltoTranscriptomeBAM")
        os.system("mkdir -p  "+globname+"_A03.Alignment/Run03.All_log.out")
        os.system("mv  "+globname+"_A03.Alignment/*Log.out  "+globname+"_A03.Alignment/Run03.All_log.out")
        os.system("mkdir -p  "+globname+"_A03.Alignment/Run04.All_log.progress.out")
        os.system("mv  "+globname+"_A03.Alignment/*Log.progress.out  "+globname+"_A03.Alignment/Run04.All_log.progress.out")
        os.system("mkdir -p  "+globname+"_A03.Alignment/Run05.All_sj.out.tab")
        os.system("mv  "+globname+"_A03.Alignment/*SJ.out.tab  "+globname+"_A03.Alignment/Run05.All_sj.out.tab")
        os.system("mkdir -p  "+globname+"_A03.Alignment/Run06.All_ReadsPerGene.final.out")
        os.system("mv  "+globname+"_A03.Alignment/*ReadsPerGene.out  "+globname+"_A03.Alignment/Run06.All_ReadsPerGene.final.out")
        os.system("mkdir -p  "+globname+"_A03.Alignment/Run07.All_genome.out")
        os.system("mv  "+globname+"_A03.Alignment/*STARgenome  "+globname+"_A03.Alignment/Run07.All_genome.out")
        os.system("mkdir -p  "+globname+"_A03.Alignment/Run08.All_STARpass1")
        os.system("mv  "+globname+"_A03.Alignment/*STARpass1  "+globname+"_A03.Alignment/Run08.All_STARpass1")
        os.system("mkdir -p  "+globname+"_A03.Alignment/Run09.All_ReadsPerGene.out.tab")
        os.system("mv  "+globname+"_A03.Alignment/*ReadsPerGene.out.tab  "+globname+"_A03.Alignment/Run09.All_ReadsPerGene.out.tab")
        os.system("mkdir -p  "+globname+"_A03.Alignment/Run10.All_Log.final.out")
        os.system("mv  "+globname+"_A03.Alignment/*Log.final.out  "+globname+"_A03.Alignment/Run10.All_Log.final.out")


        if options.less=="True":
            os.system("rm -rf Run*")
         


    elif options.alignmethods=="bowtie2":
        print("bowtie2")
        # os.system("mkdir -p  "+globname+"_A03.Alignment/Aun01.AllalignedBAM")
        # os.system("mv  "+globname+"_A03.Alignment/*.bam  "+globname+"_A03.Alignment/Aun01.AllalignedBAM")

    

    
    print("------------------------STEP3:Alignment_Parafly END-----------------------------")


def Star_genomeindex():
    print("------------------------STEP4:STAR_genomeindex-----------------------------")
    os.system("STAR --runMode genomeGenerate --genomeDir "+starindex+" --genomeFastaFiles "+genome+" --runThreadN "+str(threads)+" --sjdbGTFfile "+gtf+" --sjdbOverhang 149")
    print("------------------------STEP4:STAR_genomeindex END-----------------------------")
    print("\n")
                

class Pepline:
    def __init__(self, readFolder, threads):                                                                                                                                                                                              
        self.readFolder=readFolder
        self.threads=threads

    
    
    def Fastqrun(self,readFolder,module):
        os.system("mkdir -p "+globname+"_A02.Clean")
        os.system("rm -rf "+globname+"_A02.trim_galore.sh") 
        for file in glob.glob(readFolder+"/*_1.f*"):
            
            samplename=file.split("/")[-1].split("_1.f")[0]
            print(samplename)
            
            read1=file.split("/")[-1]
            print(samplename)
            #print(file.split("/")[-1].split("_")[0])
   
            read2=file.split("/")[-1].split("_1.f")[0]+"_2.f"+file.split("/")[-1].split("_1.f")[1]
            read2=read2.replace("_1_","_2_")
            if os.path.exists(readFolder+"/"+read2) and read2.endswith(".gz") and read1.endswith(".gz"):
                print("Found paired-end reads")
                read2=read2
                print(read2)

                run=TrimQCclass(readFolder,read1,read2, threads)   
                run.TrimQC() 
            
        if options.Parallel=="True":
            
            TrimQCclass(readFolder,read1,read2, threads).ParrunTrimQC()
        else:
            os.system("bash "+globname+"_A02.trim_galore.sh") 
        



def fastpfun(readFolder,threads=40):
    for file in glob.glob(readFolder+"/*_1.f*"):
    
        samplename=file.split("/")[-1].split("_1.f")[0]
        print(samplename)
        
        read1=file.split("/")[-1]
        print(samplename)
        #print(file.split("/")[-1].split("_")[0])

        read2=file.split("/")[-1].split("_1.f")[0]+"_2.f"+file.split("/")[-1].split("_1.f")[1]
        read2=read2.replace("_1_","_2_")
        if os.path.exists(readFolder+"/"+read2) and read2.endswith(".gz") and read1.endswith(".gz"):
            print("Found paired-end reads")
            print(read1)
            read1=readFolder+"/"+read1
            read2=readFolder+"/"+read2
            print(read2)
            outname_dir=globname+"_A02.Clean_fastp"
            outname_sample=outname_dir+"/"+samplename

            os.system("mkdir -p "+outname_dir)
            os.system("fastp -i "+ read1+ " -I "+read2+" -o "+outname_sample+"R1.fq.gz"+" -O "+outname_sample+"R2.fq.gz"+" -z 4 -q 20 -u 30 -n 0 -A  -w "+str(threads)+" -j "+outname_sample+".json"+" -h "+outname_sample+".html")

         


def rmats_pepline(readFolder,threads):
    os.system("mkdir -p A05.rmats")


def Alignbreak():
    print("------------------------STEP3:Alignment_breakpoint-----------------------------")
    if os.path.exists(globname+"_A03.Alignment"):
        lookbam=glob.glob(globname+"_A03.Alignment/*.bam")
        lookfq=glob.glob(globname+"_A02.Clean/*.gz")
        
        if len(lookbam)>0:
            for line in lookbam:
                foldline=line.split("/")[-1].split(".sorted.bam")[0]
                nameline=globname+"_A02.Clean"+"/"+line.split("/")[-1].split(".sorted.bam")[0]+"_1_val_1.fq.gz"
                print(nameline)
                if nameline in lookfq:
                    print(line)
                    os.system("mkdir -p "+globname+"_A02.Clean/RunDone")
                    os.system("mv "+globname+"_A02.Clean/"+foldline+"* "+globname+"_A02.Clean/RunDone")
        print("------------------------THis fastq don't need to run again-----------------------------")
        print("------------------------STEP3:Alignment_breakpoint END-----------------------------")


def Normalization(rootdir):
    coldata_newarg=options.coldata
    if not coldata_newarg and not options.pipepath and not rootdir :
        parser.error('Soryy,new coldata Rscript file is not given.')
    else:
        rootdir=rootdir
        os.system("mkdir -p "+globname+"_A05.Count_table")
        Newfeaturedata=globname+"_A04.Count/"+globname+"All_count.txt"
        Oldfeaturedata=rootdir+"/Root_A04.Count/Root_Allraw.count.txt"
        
        coldata_oldarg=rootdir+"/Root_A04.Count/Root_Allfeature.txt"
        data_newloaddir=globname+"_A05.Count_table"
        runmode="single"
        runname=globname
        norwayname="tpm"
        oldsave=rootdir
        rscriptmode=options.pipepath+"/Immuageclock_Run02_add.R"
        runseq=" ".join(("Rscript",rscriptmode,Newfeaturedata,Oldfeaturedata,coldata_newarg,coldata_oldarg,data_newloaddir,runmode,runname,norwayname,oldsave))
        print("\n")
        print(runseq)
        print("\n")
        os.system(runseq)


def Regression(resmode="train_pre"):
    Regression_script=options.pipepath+"/Immuageclock_Run03_regression.py"
    os.system("mkdir -p "+globname+"_A06.Regression")
    resfile= globname+"_A05.Count_table/"+globname+"_A05.combat_eata.txt "
    rescolfile=globname+"_A05.Count_table/"+globname+"_mergecoldata.txt"
    resnewcolfile=options.coldata
    if not resnewcolfile:
        parser.error('Soryy,new coldata Rscript file is not given.')
    else:
        os.system("python "+Regression_script+" -w "+resmode+" -f "+resfile+" -c "+rescolfile+" --newcol "+resnewcolfile+" --rootdir "+options.rootdir)
       





if __name__ == '__main__':
    (options, args) = parser.parse_args()

    try:
        if not options.readFolder and  options.Module=="RA":   
            parser.error('Soryy,Folder with RNA-seq data (uncompressed or gzipped fastq files) is not given.')
        else:
            nunm=glob.glob(str(options.readFolder)+"/*.gz")
            if len(nunm)==0 and options.Module=="RA":
                parser.error('Soryy,there is no uncompressed or gzipped fastq files in this fold.')
            else:
                print("\n")
                starttime = datetime.datetime.now()
                print("Start time: %s" % starttime)
                threads=str(options.threads)
                readFolder=options.readFolder
                globname=options.globname
                starindex=options.starindex
                Module=options.Module
                genome=options.genome
                gtf=options.gtf
                alignmethods=options.alignmethods
                lessmode=options.less
                countmode=options.count
                alignthreadsmod=options.alignthreads
                rscriptmode=os.path.join(options.pipepath,"Immuageclock_Run02_main.R")
                rmarkdownmode=options.rmarkdown
                coldatamode=options.coldata
                rootdir=options.rootdir

                True_satrt_init=input("Do you want to make init config?(y/n)")
                if True_satrt_init=="y":
                    generate_config=options.pipepath+"/Immuageclock_Run00_geneconfig.py"
                    os.system("python "+ generate_config)
                

                True_init=input("Do you want to run the init config?(y/n)")
                if True_init=="y":
                    con = configparser.ConfigParser()
                    init_config=os.path.join(options.pipepath,"Immuagerootconfig.ini")
                    con.read(init_config)
                    sections = con.sections()

                    items = con.items('dir') # 返回结果为元组
                    # print(items)

                    items = dict(items)
                    starindex=items["starindex"]
                    gtf=items['gtf']
                    rootdir=items['immuageroot']
                    genomefa=items['genomefa']
                    ref_cir=items['ref_cir']
                    statrindex2=items['starindex2']
                else:
                    pass



                print("\n")
                logging.basicConfig(filename=f"{globname}.log", level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')


                console = logging.StreamHandler()
                console.setLevel(logging.INFO)
                formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
                console.setFormatter(formatter)
                logging.getLogger('').addHandler(console)
                
                # 记录脚本开始运行
                logging.info(f"Start running the script with globname: {globname}")



                print("Note:Welcome this is the RNA-seq pipeline")
                print("Note:This is the version "+version)
                print("Note:The folder with RNA-seq data is : "+str(readFolder))
                print("Note:The genome you choose is: "+genome)
                print("Note:The gtf file you choose is: "+gtf)
                print("Note:The star index you choose is: "+starindex)
                print("Note:The parafly threads you choose is: "+threads)
                print("Note:The globname you choose is: "+globname)
                print("Note:The alignmethod you choose is: "+alignmethods)
                print("Note:The lessmodule you choose is: "+lessmode)
                print("Note:The countmode you choose is: "+countmode)
                print("Note:The Module you choose is: "+Module)
                print("Note:The alignthreadsmod you choose is: "+alignthreadsmod)
                print("Note:The rscriptmode you choose is: "+rscriptmode)
                # print("Note:The rmarkdownmode you choose is: "+rmarkdownmode)
                
                
                
                print("\n")
                print("------------------------STEP:RUN MODULE-----------------------------")
                print("The Module you want to run  is: "+Module)
                print("\n")

                keyopen=input("Do you sure want to run?(y/n)")
                emailget=input("Do you want to get the email?(y/n)")
                if keyopen=="y":
                    os.system("rm -rf  multiqc_data*")

                    if Module == "S1":
                            ##step1:Fastqc
                        Fastqc(readFolder)  
                    elif Module == "S2":
                        ##step2:TrimQC
                        Poline=Pepline(readFolder, threads)
                        Poline.Fastqrun(readFolder,module="T")
                    elif Module == "S2P":
                        print(readFolder)
                        fastpfun(readFolder,threads)

                    elif Module == "S3":
                        ##step3:Alignmen
                        ParrunAlign(readFolder,threads=threads,Parallel=options.Parallel)
                    elif Module == "S24":
                        # Fastqc(readFolder,fqcthereds=threads)
                        Poline=Pepline(readFolder, threads)
                        Poline.Fastqrun(readFolder,module="T")
                        Alignbreak()
                        ParrunAlign(globname+"_A02.Clean",threads=threads,Parallel=options.Parallel)
                        Count(readFolder=globname+"_A03.Alignment",threads=threads,Parallel=options.Parallel)
                    elif Module == "S14":
                        Fastqc(readFolder,fqcthereds=threads)
                        Poline=Pepline(readFolder, threads)
                        Poline.Fastqrun(readFolder,module="T")
                        Alignbreak()
                        ParrunAlign(globname+"_A02.Clean",threads=threads,Parallel=options.Parallel)
                        Count(readFolder=globname+"_A03.Alignment",threads=threads,Parallel=options.Parallel)
                        # Normalization()
                        # Regression(resmode="train_pre")
                    elif  Module == "S123":
                        # Fastqc(readFolder,fqcthereds=threads)
                        Poline=Pepline(readFolder, threads)
                        Poline.Fastqrun(readFolder,module="T")
                        Alignbreak()
                        ParrunAlign(globname+"_A02.Clean",threads=threads,Parallel=options.Parallel)


                    elif Module == "RA":
                    
                        Fastqc(readFolder,fqcthereds=threads)
                        Poline=Pepline(readFolder, threads)
                        Poline.Fastqrun(readFolder,module="T")
                        Alignbreak()
                      
                        ParrunAlign(globname+"_A02.Clean",threads=threads,Parallel=options.Parallel)
                        Count(readFolder=globname+"_A03.Alignment",threads=threads,Parallel=options.Parallel)
                        Normalization()
                        Regression(resmode="train_pre")
                    
                    elif Module == "RAR":
                        outsaveHT=globname+"_A07"
                        Fastqc(readFolder,fqcthereds=threads)
                        Poline=Pepline(readFolder, threads)
                        Poline.Fastqrun(readFolder,module="T")
                        # CirRNA=cirRNA(globname+"_A02.Clean",threads=options.alignthreads,globname=globname). ParrunAlign_circ()
                        IHT.HLAtype(readFolder=globname+"_A02.Clean",outdirname=outsaveHT)
                        IHT.TCR_get(readFolder=globname+"_A02.Clean",outdirname=outsaveHT)
                        IHT.MergeHLA(indirname=outsaveHT)
                        Alignbreak()
                        ParrunAlign(globname+"_A02.Clean",threads=threads,Parallel=options.Parallel)
                        Count(readFolder=globname+"_A03.Alignment",threads=threads,Parallel=options.Parallel)
                        Normalization()
                        Regression(resmode="train_pre")

                    
                    elif Module == "S3_1":
                        Star_genomeindex()
                    elif Module == "S4":
                        Count(readFolder=options.readFolder,threads=threads,Parallel=options.Parallel)

                    elif Module == "S34":
                        print("------------------------STEP:Alignment_and_COUNT-----------------------------")
                        ParrunAlign(readFolder,threads=threads,Parallel=options.Parallel) 
                        Count(readFolder=globname+"_A03.Alignment",threads=threads,Parallel=options.Parallel)
                        print("------------------------STEP:Alignment_and_COUNT END-----------------------------")
                    elif Module =="S5_1" :
                        # if not options.rcut and not options.rmarkdown and not options.coldata:
                        #     parser.error('Soryy,Rscript file is not given.')
                        # else:
                        os.system("mkdir -p "+globname+"_A05.Count_table")
                        data_load2=globname+"_A05.Count_table"
                        rmarkdownmode=os.path.join(options.pipepath,"Immuageclock_Run02_main.R")
                        print(rscriptmode)
                        print(globname)
                        print(data_load2)
                        print(rmarkdownmode)
                        os.system("Rscript  "+rscriptmode+" "+globname+"_A04.Count/"+globname+"All_count.txt"+" "+coldatamode+" "+data_load2+" "+rmarkdownmode)   
                    elif options.add=="ap" and Module =="S5":
                        Normalization(rootdir=rootdir)
                    elif Module =="S6":
                        Regression(resmode=options.resmode)
                    
                    elif Module =="S8cir" : 
                        if  options.alignmethods=="STARcir":
                            print("------------------------STEP:Circle_Alignment-----------------------------")
                            # print("------------------------STEP:AlignmentCIR-----------------------------")
                            CirRNA=cirRNA(readFolder,threads=options.alignthreads,globname=globname). ParrunAlign_circ()
                            # ParrunAlign(readFolder,threads=threads,Parallel=options.Parallel)
                        else: 
                            print("\n")
                            print("Error:if paraments Module is circle and paraments alignmethods must be  STARcir")

                        print("\n")                        
                        print("------------------------STEP:Circle_Alignment END-----------------------------")
                    elif Module =="S7":
                        #  HLAtype(readFolder=args.input[0],outdirname=outsave)
                        print("------------------------STEP:HLA_TCR_Type-----------------------------")
                        outsaveHT=globname+"_A07"                      
                        IHT.HLAtype(readFolder=readFolder,outdirname=outsaveHT)
                        IHT.TCR_get(readFolder=readFolder,outdirname=outsaveHT)
                        IHT.MergeHLA(indirname=outsaveHT)
                        print("------------------------STEP:HLA_TCR_Type END-----------------------------")
                    
                    elif Module =="S7HLA":
                        print("------------------------STEP:HLA_Type-----------------------------")
                        outsaveHT=globname+"_A07"
                        IHT.HLAtype(readFolder=readFolder,outdirname=outsaveHT)
                        print("------------------------STEP:HLA_Type END-----------------------------")
                    elif Module =="S7TCR":
                        print("------------------------STEP:TCR_Type-----------------------------")                 
                        outsaveHT=globname+"_A07"
                        os.system("mkdir -p "+globname+"_A07")
                        IHT.TCR_get(readFolder=readFolder,outdirname=outsaveHT)
                        print("------------------------STEP:TCR_Type END-----------------------------")
                    elif Module =="S78":
                        if  options.alignmethods=="STARcir":
                            print("------------------------STEP:Cir_HLA_TCR_Alignment-----------------------------")
                            # print("------------------------STEP:AlignmentCIR-----------------------------")
                            CirRNA=cirRNA(readFolder,threads=options.alignthreads,globname=globname). ParrunAlign_circ()
                            
                            outsaveHT=globname+"_A07"                       
                            IHT.HLAtype(readFolder=readFolder,outdirname=outsaveHT)
                            IHT.TCR_get(readFolder=readFolder,outdirname=outsaveHT)
                            IHT.MergeHLA(indirname=outsaveHT)
                            print("------------------------STEP:CIR_HLA_TCR_Type END-----------------------------") 

                            # ParrunAlign(readFolder,threads=threads,Parallel=options.Parallel)
                        else: 
                            print("\n")
                            print("Error:Module is circle and alignmethods is  STARcir")

                    


                    endtime = datetime.datetime.now()
                    print("End time: %s" % endtime)     

                    
                    print("Total time: %s" % (endtime - starttime))
                    
                    logging.info("Script finished")
                    if emailget=="y":
                        keeptime = str(endtime - starttime)
                        messageattach = globname + "05.Count_table/Immuclock_trandata.html"
                        emailscript = "/home/maolp/mao/Codeman/All_Archived_Project/SentEmail.py"

                        # Build the command with safe argument passing
                        command = ['python', emailscript, f"{globname}_{keeptime}", messageattach]

                        # Run the command safely with subprocess.run
                        subprocess.run(command, check=True)

                    else:
                        print("\n")
                        print("Note:You do not want to get the email")
                        print("\n")
                
                else:
                    print("You choose not to run")
    except KeyboardInterrupt:
        print("\n")
        os.system("rm -rf  multiqc_data_*")
        print("Note:Interrupt..")
        print("\n")
                        


