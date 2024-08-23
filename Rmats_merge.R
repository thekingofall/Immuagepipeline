
# args <- commandArgs(trailingOnly = TRUE)
# print(args)
library(tidyverse)
    name=paste("*","JC.txt",sep="")
    
for (i in list.files(pattern = name) ){
    
    data_name=substr(i,1,nchar(i)-8)

    print(data_name)
    assign(data_name,read.table(i,header=T,sep = "\t"))

    
    assign(paste(data_name,"_FDR",sep=""),filter(get(data_name),FDR<0.05 & abs(IncLevelDifference) > 0.05))

    data_type=substr(i,1,nchar(i)-12)
    print(data_type)
    type=rep(data_type,nrow(get(paste(data_name,"_FDR",sep=""))))
    print(type)
    # assign(get(paste(data_name,"_FDR",sep="")),cbind(type,get(paste(data_name,"_FDR",sep=""))))
    assign(paste(data_name,"_FDR2",sep=""),cbind(type,get(paste(data_name,"_FDR",sep=""))))
    assign(paste(data_name,"_FDR2",sep=""),get(paste(data_name,"_FDR2",sep=""))[,c("type","ID","GeneID", "geneSymbol","FDR","IncLevelDifference" )])

}

    

data_list=list(A3SS.MAT_FDR2,A5SS.MAT_FDR2,MXE.MAT_FDR2,RI.MAT_FDR2,SE.MAT_FDR2)
do.call(rbind,data_list)-> data_last
sname=paste("eff.txt",sep="")
write.table(data_last,file=sname,sep="\t",quote = F,row.names = F,col.names = T)


# data_last=data_last[data_last$FDR<0.05 & abs(data_last$IncLevelDifference) > 0.05,]

# setwd("/home/maolp/mao/Alice/xianchong_data/A03.Alignment/Run01.AllalignedBAM/outbam/OUTCOME/CK_NR_D1_vs_CK_D1")
# data=read.table("A3SS.MATS.JCEC.txt",header=T,sep = "\t")
# head(data)
# data2=filter(data,FDR<0.05 & abs(IncLevelDifference) > 0.05)
# head(data2)
# data2[,"GeneID"]
# # dir()
