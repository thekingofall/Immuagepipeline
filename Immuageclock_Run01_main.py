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

parser = OptionParser(usage="""
Usage: %prog -f READ_FOLDER -m MODULE [options]

Description:
  RNA-seq analysis pipeline for immune age prediction, supporting multiple analysis modules and species.

Examples:
  # Run complete analysis on human samples
  %prog -f /path/to/fastq_files -m RA --species human
  
  # Run only QC and trimming on mouse samples
  %prog -f /path/to/fastq_files -m S2 --species mouse
  
  # Run with custom genome and GTF
  %prog -f /path/to/fastq_files -m RA -g /path/to/genome.fa -u /path/to/annotation.gtf
""", version="%prog "+version)

# Required parameters
parser.add_option("-f",
        action="store", 
        dest="readFolder",
        help="[REQUIRED] Path to folder containing RNA-seq data (gzipped FASTQ files). Required for most analysis modules.")

# Analysis configuration
parser.add_option("-m","--Module",
                    action="store",
                    dest="Module",
                    default="S14",
                    help="""Analysis module to run. Options include:
  S1: FastQC quality control only
  S2: TrimQC (adapter/quality trimming)
  S2P: FastP trimming (alternative to S2)
  S3: Alignment only
  S3_1: Build STAR index
  S4: Count only (from existing alignments)
  S14: QC through count (S1+S2+S3+S4) (default)
  S24: Trim through count (S2+S3+S4)
  S123: QC through alignment (S1+S2+S3)
  S5: Normalization
  S6: Regression analysis
  S78: TCR/HLA analysis
  S7: TCR/HLA typing only
  S8cir: circRNA analysis
  RA: Complete analysis 
  RAR: Complete analysis with TCR/HLA typing""")

# Reference files
parser.add_option("-o", "--species",
                  action="store",
                  dest="species",
                  default="human",
                  help="Species to analyze. Options: 'human' (default) or 'mouse'. Sets appropriate reference genome automatically.")

parser.add_option("-g","--genome",
                    action="store",
                    dest="genome",
                    default="/home/maolp/mao/Ref/",
                    help="Path to reference genome FASTA file. Default: '/home/maolp/mao/Ref/'")

parser.add_option("-u","--gtf",
                    action="store",
                    dest="gtf",
                    default="", 
                    help="Path to gene annotation GTF file. If empty, will be set based on --species.")

parser.add_option("-s","--starindex",  
                    action="store",
                    dest="starindex",
                    default="", 
                    help="Path to genome index for alignment. If empty, will be set based on --species.")

parser.add_option("-x","--rsemindex",
                    action="store",
                    dest="rsemindex",
                    default="/home/maolp/mao/Ref/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/HG38RSEM/HG38RSEM", 
                    help="Path to RSEM index for transcript quantification. Default: '/home/maolp/mao/Ref/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/HG38RSEM/HG38RSEM'")

# Alignment and counting options
parser.add_option("-a","--alignmethods",
                    action="store",
                    dest="alignmethods",
                    default="hisat2",
                    help="Alignment method to use. Options: 'hisat2' (default), 'bowtie2', 'STAR', 'STARfast', 'STARcir', or 'salmon'.")

parser.add_option("-t","--alignthreads",
                    action="store",         
                    dest="alignthreads",
                    default="40",
                    help="Number of threads for alignment tools. Default: 40")

parser.add_option("-c","--count",
                    action="store",
                    dest="count",
                    default="featureCounts",
                    help="Method for gene counting. Options: 'featureCounts' (default), 'htseq-count', or 'rsem'.")

# Parallelization options
parser.add_option("-e", "--threads",
                    action="store",
                    dest="threads",
                    default="8",
                    help="Number of threads for overall pipeline parallelization with ParaFly. Default: 8")

parser.add_option("-p", "--Parallel",
                    action="store",
                    dest="Parallel",
                    default="True",
                    help="Enable parallel processing with ParaFly. Options: 'True' or 'False'. Default: True")

# Output and environment options
parser.add_option("-n","--globname",
                    action="store",
                    dest="globname",
                    default="RNAc",
                    help="Prefix for output files and directories. Default: 'RNAc'")

parser.add_option("--rootdir",
                    action="store", 
                    dest="rootdir",
                    default="/data4/maolp/All_myself_data4/Immu_age/Immuageroot",
                    help="Root directory for results. Default: '/data4/maolp/All_myself_data4/Immu_age/Immuageroot'")

parser.add_option("--pipepath",
                    action="store",
                    dest="pipepath",
                    default="/home/maolp/mao/Codeman/All_Archived_Project/Immuagepipeline",
                    help="Path to pipeline directory with supporting scripts. Default: '/home/maolp/mao/Codeman/All_Archived_Project/Immuagepipeline'")

parser.add_option("-l","--less",
                    action="store",
                    dest="less",    
                    default="True",
                    help="Clean up intermediate files to save disk space. Options: 'True' (default) or 'False'.")

parser.add_option("--conda-env",
                  action="store",
                  dest="conda_env",
                  default="TCRRNA",
                  help="Conda environment to use for analysis. Default: 'TCRRNA'. Will be created if it doesn't exist.")

# R analysis options
parser.add_option("-w","--rmarkdown",
                    action="store",
                    dest="rmarkdown",
                    default="/home/maolp/mao/Codeman/Immuagepipeline/Immuageclock_Run02_PCA.Rmd",
                    help="Path to RMarkdown template for results visualization. Default: '/home/maolp/mao/Codeman/Immuagepipeline/Immuageclock_Run02_PCA.Rmd'")

parser.add_option("-r","--rcut",
                    action="store",
                    dest="rcut",
                    default="/home/maolp/mao/Codeman/Immuagepipeline/Immuage_Rmain.R",
                    help="Path to R script for analysis. Default: '/home/maolp/mao/Codeman/Immuagepipeline/Immuage_Rmain.R'")

# Sample metadata
parser.add_option("--coldata",
                    action="store",
                    dest="coldata",
                    help="Path to sample metadata file (required for regression analysis modules). No default.")

parser.add_option("--add",
                    action="store",
                    dest="add",
                    default="ap",
                    help="Sample addition mode. Options: 'aa' (add all), 'ap' (add part, default), or 'an' (add none).")

# Advanced options
parser.add_option("--resmode",
                    action="store", 
                    dest="resmode",
                    default="train_pre",
                    help="Regression model mode. Options: 'train', 'test', or 'train_pre' (default).")

parser.add_option("--config",
                    action="store",
                    dest="config",
                    help="Path to configuration file for pipeline settings. Default: None")

parser.add_option("--health",
                    action="store",
                    dest="helath",
                    default="False",    
                    help="Update/upgrade pipeline components. Options: 'True' or 'False' (default).")




def Fastqc(fold, fqcthereds=8):
    print("------------------------STEP1:QCreport-----------------------------")
    starttime = datetime.datetime.now()
    print("Start time: %s" % starttime)
    path = os.getcwd()
    
    # Create fastqc output directory
    fastqc_dir = f"{path}/{globname}_A01.raw.fastqc_report"
    os.system(f"mkdir -p {fastqc_dir}")
    
    # Run fastqc on all files
    os.system(f"ls {fold}/*gz | xargs -I [] echo 'fastqc -t 15 [] -o {fastqc_dir}' > {globname}_A01.fastqc.sh")
    os.system(f"ParaFly -c {globname}_A01.fastqc.sh -CPU {str(fqcthereds)} -failed_cmds {globname}_A01.fastqc.failed -v")
    # os.system(f"rm -rf {globname}_A01.fastqc.sh")

    # Run multiqc in the fastqc directory
    os.system(f"cd {fastqc_dir} && multiqc .")
    
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
        os.system(f"ParaFly -c {globname}_A02.trim_galore.sh -CPU {str(self.fqcthereds)} -failed_cmds {globname}_A02.trim_galore.failed -v")
        
        # Create fastqc report directory
        fastqc_dir = f"{globname}_A02.Clean.fastqc_report"
        os.system(f"mkdir -p {fastqc_dir}")
        
        # Run fastqc on cleaned files
        Fastqc(fold=f"{globname}_A02.Clean")
        
        # Move txt files and run multiqc
        os.system(f"mv {globname}_A02.Clean/*txt {fastqc_dir}")
        os.system(f"cd {fastqc_dir} && multiqc .")
        
        print("------------------------STEP2:TrimQC_Parafly END-----------------------------")



def Align(readFile1,readFile2,threads=40,Parallel="False",samename="Atest",alignoutname="",index=""):
    print("------------------------STEP3:Aliging-----------------------------")
    alignstarttime = datetime.datetime.now()
    print("Start time: %s" % alignstarttime)

    Al_fq1=readFile1
    Al_fq2=readFile2

    os.system("mkdir -p  "+globname+"_A03.Alignment")

    if options.alignmethods=="bowtie2":
        cmd = f"bowtie2 -p {threads} -x {options.starindex} -1 {Al_fq1} -2 {Al_fq2} -S {globname}_A03.Alignment/{samename}.sam 2> {globname}_A03.Alignment/{samename}.align.log"
        logging.info(f"Running command: {cmd}")
        os.system(cmd)
        
        # Convert to BAM and sort
        os.system(f"samtools view -bS {globname}_A03.Alignment/{samename}.sam > {globname}_A03.Alignment/{samename}.bam 2>> {globname}_A03.Alignment/{samename}.align.log")
        os.system(f"samtools sort -@ {threads} {globname}_A03.Alignment/{samename}.bam -o {globname}_A03.Alignment/{samename}.sorted.bam 2>> {globname}_A03.Alignment/{samename}.align.log")
        os.system(f"samtools index {globname}_A03.Alignment/{samename}.sorted.bam 2>> {globname}_A03.Alignment/{samename}.align.log")
        os.system(f"rm {globname}_A03.Alignment/{samename}.sam")
        os.system(f"rm {globname}_A03.Alignment/{samename}.bam")

    elif options.alignmethods=="STAR":
        cmd = f"STAR --runThreadN {threads} \
        --runMode alignReads \
        --readFilesCommand zcat --quantMode TranscriptomeSAM GeneCounts \
        --twopassMode Basic \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped None \
        --genomeDir {options.starindex} \
        --readFilesIn {Al_fq1} {Al_fq2} \
        --outFileNamePrefix {globname}_A03.Alignment/{samename}_ 2> {globname}_A03.Alignment/{samename}.align.log"
        
        logging.info(f"Running STAR alignment for {samename}")
        logging.info(f"Command: {cmd}")
        os.system(cmd)

    elif options.alignmethods=="STARfast":
        os.system("STAR  --runThreadN "+str(threads)+" \
        --runMode alignReads \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --genomeDir "+ options.starindex + " \
        --readFilesIn  " + Al_fq1 + " "+ Al_fq2 +" --outFileNamePrefix  "+globname+"_A03.Alignment/"+samename+"_")
        #os.system("samtools sort -@ 40  "+globname+"_A03.Alignment/"+samename+"_Aligned.out.bam  -o  "+globname+"_A03.Alignment/"+samename+"_Aligned.out.sorted.bam")
        
        os.system("samtools sort -@ "+str(threads)+" -o  "+globname+"_A03.Alignment/"+samename+"_Aligned.sorted.bam  "+globname+"_A03.Alignment/"+samename+"_Aligned.out.bam")
        os.system("samtools index  "+globname+"_A03.Alignment/"+samename+"_Aligned.out.sorted.bam")
        os.system("rm -rf   "+globname+"_A03.Alignment/"+samename+"_Aligned.out.bam")

    elif options.alignmethods=="hisat2":
        cmd = f"hisat2 -p {options.alignthreads} -x {options.starindex} -1 {Al_fq1} -2 {Al_fq2} -S {globname}_A03.Alignment/{samename}.sam 2> {globname}_A03.Alignment/{samename}.align.log"
        logging.info(f"Running hisat2 alignment for {samename}")
        logging.info(f"Command: {cmd}")
        os.system(cmd)

        # Convert to BAM and sort
        os.system(f"/home/maolp/mao/Biosoft/samtools/samtools-1.9/samtools sort -@ {options.alignthreads} {globname}_A03.Alignment/{samename}.sam -o {globname}_A03.Alignment/{samename}.sorted.bam 2>> {globname}_A03.Alignment/{samename}.align.log")
        os.system(f"/home/maolp/mao/Biosoft/samtools/samtools index {globname}_A03.Alignment/{samename}.sorted.bam 2>> {globname}_A03.Alignment/{samename}.align.log")
        os.system(f"rm {globname}_A03.Alignment/{samename}.sam")

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


def process_rnaseq_logs(align_dir, count_log, output_dir, prefix):
    """处理RNA-seq比对和定量日志文件"""
    # ... copy content from process_rnaseq_logs.py ...
    os.makedirs(output_dir, exist_ok=True)
    
    # 读取比对日志
    align_stats = []
    align_logs = glob.glob(os.path.join(align_dir, "*.log"))
    
    if not align_logs:
        print(f"Warning: No log files found in {align_dir}")
        return
        
    print(f"Found {len(align_logs)} alignment log files")
    
    # ... rest of the log processing code ...

def run_rnaseq_analysis(input_file, name="New", out=None, min_count=0, 
                       width=10, height=6, dpi=300):
    """运行RNA-seq分析"""
    # 构建R脚本命令
    r_script = os.path.join(os.path.dirname(__file__), "rna_seq_analysis.R")
    cmd = ["Rscript", r_script,
           "-i", input_file,
           "-n", name,
           "-m", str(min_count),
           "-w", str(width),
           "--height", str(height),
           "-d", str(dpi)]
    
    if out:
        cmd.extend(["-o", out])
        
    # 执行R脚本
    try:
        subprocess.run(cmd, check=True)
        print(f"RNA-seq analysis completed successfully. Results in: {name}_A05.Rcount")
    except subprocess.CalledProcessError as e:
        print(f"Error running RNA-seq analysis: {e}")
        logging.error(f"RNA-seq analysis failed: {e}")

def Count(readFolder="A03.Alignment", threads=40, Parallel="False", samename="All"):
    """计数和分析RNA-seq数据"""
    datetime.datetime.now()
    os.system("mkdir -p "+globname+"_A04.Count")
    print("------------------------STEP4:Count-----------------------------")
    
    if options.count == "htseq-count":
        cmd = f"htseq-count -f bam -r pos -s no -t exon -i gene_id -a {globname}_A04.Count/{options.globname}{samename}_count.txt {readFolder}/{samename}_Aligned.sorted.bam {options.gtf} 2> {globname}_A04.Count/{samename}.count.log"
        logging.info(f"Running htseq-count for {samename}")
        logging.info(f"Command: {cmd}")
        os.system(cmd)

    elif options.count == "featureCounts":
        # 构建正确的BAM文件路径
        if readFolder == "A03.Alignment":  # 使用默认值的情况
            bam_path = f"{globname}_{readFolder}/*.sorted.bam"
        else:  # 使用自定义路径的情况
            bam_path = f"{readFolder}/*.sorted.bam"
            
        cmd = f"featureCounts -t exon -p -g gene_id -T {str(options.alignthreads)} \
            -a {options.gtf} \
            -o {globname}_A04.Count/{globname}{samename}_count.txt \
            {bam_path} 2> {globname}_A04.Count/{samename}.count.log"
            
        logging.info(f"Running featureCounts with command:")
        logging.info(cmd)
        print(f"Looking for BAM files in: {bam_path}")
        os.system(cmd)

    elif options.count == "rsem":
        os.system("mkdir -p  "+globname+"_A04.Count")
        for samp in glob.glob(readFolder+"/Run02.AlltoTranscriptomeBAM/*.bam"):
            newname=samp.split("/")[-1].split(".")[0]
            cmd = f"rsem-calculate-expression --paired-end -no-bam-output --alignments -p {str(threads)} {samp} {options.rsemindex} {globname}_A04.Count/{newname}_RSEM 2> {globname}_A04.Count/{newname}.rsem.log"
            logging.info(f"Running RSEM for {newname}")
            logging.info(f"Command: {cmd}")
            os.system(cmd)
    
    # 添加RNA-seq分析
    count_file = f"{globname}_A04.Count/{globname}{samename}_count.txt"
    if os.path.exists(count_file):
        print("\nRunning RNA-seq analysis...")
        run_rnaseq_analysis(
            input_file=count_file,
            name=globname,
            out=f"{globname}_A05.Rcount"
        )
        
        # 处理日志文件
        print("\nProcessing alignment and count logs...")
        process_rnaseq_logs(
            align_dir=f"{globname}_A03.Alignment",
            count_log=f"{globname}_A04.Count/{samename}.count.log",
            output_dir=f"{globname}_A04.Count/logs",
            prefix=globname
        )
    else:
        print(f"Warning: Count file not found: {count_file}")
        logging.warning(f"Count file not found: {count_file}")

    print("------------------------STEP4:Count END-----------------------------")


def ParrunAlign(starfold, threads=40, Parallel="False"):
    print("------------------------STEP3:Alignment_Parafly START-----------------------------")
    
    # 同时匹配 _1.f* 和 _R1.f* 模式
    for line in glob.glob(starfold+"/*[_R]1.f*"):
        # 分别处理两种命名模式
        if "_R1.f" in line:
            samplename = line.split("/")[-1].split("_R1")[0]
            read1 = line
            read2 = line.replace("_R1", "_R2")
        else:
            samplename = line.split("/")[-1].split("_1.f")[0]
            read1 = line
            read2 = line.replace("_1_val_1", "_2_val_2")
            
        print(f"Processing sample: {samplename}")
        print(f"Read1: {read1}")
        print(f"Read2: {read2}")
        
        if os.path.exists(read2):
            print("Found paired-end reads")
            Align(readFile1=read1, readFile2=read2, threads=options.alignthreads, 
                  Parallel=Parallel, samename=samplename)
        else:
            print(f"Warning: Could not find matching read2 file: {read2}")
            logging.warning(f"Missing read2 file for {read1}")

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
    os.system("STAR --runMode genomeGenerate --genomeDir "+options.starindex+" --genomeFastaFiles "+options.genome+" --runThreadN "+str(options.threads)+" --sjdbGTFfile "+options.gtf+" --sjdbOverhang 149")
    print("------------------------STEP4:STAR_genomeindex END-----------------------------")
    print("\n")
                

class Pepline:
    def __init__(self, readFolder, threads):                                                                                                                                                                                              
        self.readFolder=readFolder
        self.threads=threads

    
    
    def Fastqrun(self,readFolder,module):
        os.system("mkdir -p "+globname+"_A02.Clean")
        os.system("rm -rf "+globname+"_A02.trim_galore.sh") 
        
        # 标记是否找到了有效的文件对
        found_valid_pairs = False
        
        # 同时匹配 _1.f* 和 _R1.f* 模式的文件
        for file in glob.glob(readFolder+"/*[_R]1.f*"):
            # 根据文件名模式分别处理
            if "_R1.f" in file:
                samplename = file.split("/")[-1].split("_R1.f")[0]
                read1 = file.split("/")[-1]
                read2 = file.split("/")[-1].split("_R1.f")[0]+"_R2.f"+file.split("/")[-1].split("_R1.f")[1]
            else:
                samplename = file.split("/")[-1].split("_1.f")[0]
                read1 = file.split("/")[-1]
                read2 = file.split("/")[-1].split("_1.f")[0]+"_2.f"+file.split("/")[-1].split("_1.f")[1]
                read2 = read2.replace("_1_","_2_")
            
            print(f"Processing sample: {samplename}")
            print(f"Read1: {read1}")
            print(f"Read2: {read2}")
            
            if os.path.exists(readFolder+"/"+read2) and read2.endswith(".gz") and read1.endswith(".gz"):
                print("Found paired-end reads")
                found_valid_pairs = True
                
                run = TrimQCclass(readFolder,read1,read2, self.threads)   
                run.TrimQC() 
        
        if found_valid_pairs:
            if options.Parallel=="True":
                run.ParrunTrimQC()
            else:
                os.system("bash "+globname+"_A02.trim_galore.sh")
        else:
            print("No valid paired-end reads found in the input folder")
            logging.warning("No valid paired-end reads found in folder: %s", readFolder)



def fastpfun(readFolder, threads=40):
    # 同时匹配 _1.f* 和 _R1.f* 模式
    for file in glob.glob(readFolder+"/*[_R]1.f*"):
        # 分别处理两种命名模式
        if "_R1.f" in file:
            samplename = file.split("/")[-1].split("_R1.f")[0]
            read1 = file.split("/")[-1]
            read2 = file.split("/")[-1].split("_R1.f")[0]+"_R2.f"+file.split("/")[-1].split("_R1.f")[1]
        else:
            samplename = file.split("/")[-1].split("_1.f")[0]
            read1 = file.split("/")[-1]
            read2 = file.split("/")[-1].split("_1.f")[0]+"_2.f"+file.split("/")[-1].split("_1.f")[1]
            read2 = read2.replace("_1_","_2_")
            
        if os.path.exists(readFolder+"/"+read2) and read2.endswith(".gz") and read1.endswith(".gz"):
            print("Found paired-end reads")
            print(read1)
            read1 = readFolder+"/"+read1
            read2 = readFolder+"/"+read2
            print(read2)
            outname_dir = globname+"_A02.Clean_fastp"
            outname_sample = outname_dir+"/"+samplename

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
       

def setup_conda_environment(env_name="TCRRNA"):
    """
    Check if specified conda environment exists, create it if it doesn't,
    and activate it for the current script execution.
    """
    print(f"Checking conda environment '{env_name}'...")
    logging.info(f"Checking conda environment '{env_name}'")
    
    # Check if conda is available
    try:
        subprocess.run(["conda", "--version"], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("Error: conda is not installed or not available in PATH")
        logging.error("conda is not installed or not available in PATH")
        return False
        
    # Check if environment exists
    env_check = subprocess.run(
        ["conda", "env", "list"], 
        check=True, 
        stdout=subprocess.PIPE,
        universal_newlines=True
    )
    
    env_exists = False
    for line in env_check.stdout.split('\n'):
        if line.startswith(env_name + ' ') or f"/{env_name}" in line:
            env_exists = True
            break
    
    if not env_exists:
        create_env = input(f"Conda environment '{env_name}' not found. Create it? (y/n): ")
        if create_env.lower() == 'y':
            print(f"Creating conda environment '{env_name}'...")
            logging.info(f"Creating conda environment '{env_name}'")
            
            # Create environment with necessary packages
            create_cmd = [
                "conda", "create", "-y", "-n", env_name,
                "python=3.7", "bowtie2", "star", "htseq-count", "bioconductor-rsubread", 
                "rsem", "samtools", "parafly", "multiqc", "r-base", "fastqc", "trim-galore"
            ]
            
            try:
                subprocess.run(create_cmd, check=True)
                print(f"Successfully created conda environment '{env_name}'")
                logging.info(f"Successfully created conda environment '{env_name}'")
            except subprocess.CalledProcessError as e:
                print(f"Error creating conda environment: {e}")
                logging.error(f"Error creating conda environment: {e}")
                return False
        else:
            print(f"Environment creation skipped. Using current environment.")
            return True
    
    # Check if the environment is already active
    current_env = os.environ.get('CONDA_DEFAULT_ENV')
    if current_env == env_name:
        print(f"Conda environment '{env_name}' is already active")
        return True
        
    # Try to activate the environment
    print(f"Activating conda environment '{env_name}'...")
    try:
        # Get the conda base path
        conda_path = subprocess.run(
            ["conda", "info", "--base"],
            check=True,
            stdout=subprocess.PIPE,
            universal_newlines=True
        ).stdout.strip()
        
        # Source the conda.sh script to get conda functions
        activate_cmd = f"source {conda_path}/etc/profile.d/conda.sh && conda activate {env_name}"
        
        # Use bash to execute the activate command and modify the current environment
        activate_env = subprocess.run(
            ["bash", "-c", f"{activate_cmd} && exec python -c \"import os; print('CONDA_DEFAULT_ENV=' + os.environ.get('CONDA_DEFAULT_ENV', ''))\""],
            check=True,
            stdout=subprocess.PIPE,
            universal_newlines=True
        )
        
        if f"CONDA_DEFAULT_ENV={env_name}" in activate_env.stdout:
            print(f"Successfully activated conda environment '{env_name}'")
            logging.info(f"Successfully activated conda environment '{env_name}'")
            
            # Update PATH and other environment variables
            env_vars = subprocess.run(
                ["bash", "-c", f"{activate_cmd} && env"],
                check=True,
                stdout=subprocess.PIPE,
                universal_newlines=True
            ).stdout.strip().split('\n')
            
            for line in env_vars:
                if '=' in line:
                    key, value = line.split('=', 1)
                    os.environ[key] = value
                    
            return True
        else:
            print(f"Warning: Failed to activate conda environment '{env_name}'")
            print(f"Please manually activate with 'conda activate {env_name}' before running this script")
            logging.warning(f"Failed to activate conda environment '{env_name}'")
            return False
            
    except subprocess.CalledProcessError as e:
        print(f"Error activating conda environment: {e}")
        logging.error(f"Error activating conda environment: {e}")
        return False


if __name__ == '__main__':
    (options, args) = parser.parse_args()

    try:
        # Setup conda environment
        conda_env_ready = setup_conda_environment(options.conda_env)
        if not conda_env_ready:
            print(f"Warning: Could not set up conda environment '{options.conda_env}'.")
            use_anyway = input("Continue with current environment anyway? (y/n): ")
            if use_anyway.lower() != 'y':
                print("Exiting.")
                sys.exit(1)
        
        # Set defaults based on species
        if options.species.lower() == "mouse" or options.species.lower() == "mm10":
            print("Setting up reference files for mouse (mm10)...")
            if not options.starindex:
                options.starindex = "/data1/Ref/Histat2/fast/mm10/genome"
            if not options.gtf:
                options.gtf = "/data1/Ref/Histat2/fast/mm10.ncbiRefSeq.gtf"
            print(f"Using mouse genome index: {options.starindex}")
            print(f"Using mouse GTF file: {options.gtf}")
        else:
            # Default to human
            if not options.starindex:
                options.starindex = "/home/maolp/mao/Ref/hg38/genome"
            if not options.gtf:
                options.gtf = "/home/maolp/mao/Ref/Homo_sapiens/UCSC/hg38/Annotation/Genes.gencode/genes.gtf"
            print(f"Using human genome index: {options.starindex}")
            print(f"Using human GTF file: {options.gtf}")
        
        if not options.readFolder and options.Module=="RA":   
            parser.error('Sorry, folder with RNA-seq data (uncompressed or gzipped fastq files) is not given.')
        else:
            nunm=glob.glob(str(options.readFolder)+"/*.gz")
            if len(nunm)==0 and options.Module=="RA":
                parser.error('Sorry,there is no uncompressed or gzipped fastq files in this fold.')
            else:
                print("\n")
                starttime = datetime.datetime.now()
                print("Start time: %s" % starttime)
                
                # Add current working directory logging
                current_dir = os.getcwd()
                print(f"Working directory: {current_dir}")
                logging.info(f"Working directory: {current_dir}")
                
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
                print("Note:The species you selected is: "+options.species)
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
                    # Remove any existing multiqc data at start
                    os.system("rm -rf multiqc_data")

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
                        


