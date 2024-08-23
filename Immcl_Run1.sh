#!/bin/bash

# 读取 samples.txt 文件的每一行
while read samename; do

    # 假设对于每个样本，你有一对特定的 FASTQ 文件
    FQ1="${samename}_R1.fastq"
    FQ2="${samename}_R2.fastq"

    # 其他的变量也可以在这里定义，例如 path, globname, options_alignthreads, options_starindex, Al_fq1, Al_fq2, gtf, options_globname, readFolder

    # 调用 trim_galore
    trim_galore -q 25  --phred33 --length 36 -e 0.1 --stringency 3 --paired  "$FQ1" "$FQ2"  -o "${path}/${globname}_A02.Clean"

    # 调用 hisat2
    hisat2 -p "$options_alignthreads" -x "$options_starindex" -1 "$Al_fq1" -2 "$Al_fq2" -S  "${globname}_A03.Alignment/${samename}.sam"

    # 调用 featureCounts
    featureCounts -t exon  -g gene_id -T "$options_alignthreads" -a "$gtf" -o  "${globname}_A04.Count/${options_globname}${samename}_count.txt" "$readFolder/*.bam"

# 从 samples.txt 文件中读取输入
done < samples.txt