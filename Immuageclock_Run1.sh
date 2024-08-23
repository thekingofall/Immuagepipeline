#!/bin/bash


while read samename; do


    FQ1="${samename}_R1.fastq"
    FQ2="${samename}_R2.fastq"


    trim_galore -q 25  --phred33 --length 36 -e 0.1 --stringency 3 --paired  "$FQ1" "$FQ2"  -o "${path}/${globname}_A02.Clean"


    hisat2 -p "$options_alignthreads" -x "$options_starindex" -1 "$Al_fq1" -2 "$Al_fq2" -S  "${globname}_A03.Alignment/${samename}.sam"

    
    featureCounts -t exon  -g gene_id -T "$options_alignthreads" -a "$gtf" -o  "${globname}_A04.Count/${options_globname}${samename}_count.txt" "$readFolder/*.bam"


done < samples.txt