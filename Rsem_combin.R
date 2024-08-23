rm(list = ls())
# setwd("/home/maolp/mao/Alice/xianchong_data/A03.Alignment/Run01.AllalignedBAM/tranbam")
setwd("/data/maolp/Alice/mouse_data/mouse_raw/Mouse_starnew_A03.Alignment/Run02.AlltoTranscriptomeBAM/ISOresultes")
suppressMessages(library(tidyverse))
list.files(pattern = "*isoforms.results")



# data1=read.table(isoforms_results[1],sep="\t",header=T)
# data2=read.table(isoforms_results[3],sep="\t",header=T)
# all(data1$transcript_id==data2$transcript_id)
# help(splice)
# ??str_split
Rsem_get<-function(isoforms_results,Types){
if (Types=="count"){
    outtab=data.frame()
for(i in 1:length(isoforms_results)){
    iso_name=str_split(isoforms_results[i], "\\.")[[1]][1]
    print(iso_name)
    if(i==1){
    
    data=read.table(isoforms_results[i],sep="\t",header=T)
    data_count=data[,c(1,2,3,5)]
    colnames(data_count)[4]<-paste(iso_name,"_count",sep="")
    outtab=rbind(outtab,data_count)
    }else{
    data=read.table(isoforms_results[i],sep="\t",header=T)
    data_count=data.frame(data[,c(5)])
    names(data_count)=paste(iso_name,"_count",sep="")
    outtab=cbind(outtab,data_count[,1])
    colnames(outtab)[length(colnames(outtab))]<-paste(iso_name,"_count",sep="")
    }
    }
}

else if(Types=="TPM"){
    outtab=data.frame()
 for(i in 1:length(isoforms_results)){
    iso_name=str_split(isoforms_results[i], "\\.")[[1]][1]
    print(iso_name)
    if(i==1){
    
    data=read.table(isoforms_results[i],sep="\t",header=T)
    data_count=data[,c(1,2,3,6)]
    colnames(data_count)[4]<-paste(iso_name,"_TPM",sep="")
    outtab=rbind(outtab,data_count)
    }else{
    data=read.table(isoforms_results[i],sep="\t",header=T)
    data_count=data.frame(data[,c(6)])
    names(data_count)=paste(iso_name,"_TPM",sep="")
    outtab=cbind(outtab,data_count[,1])
    colnames(outtab)[length(colnames(outtab))]<-paste(iso_name,"_TPM",sep="")
    }
    }
}   




else if (Types=="FPKM"){
    outtab=data.frame()
 for(i in 1:length(isoforms_results)){
    iso_name=str_split(isoforms_results[i], "\\.")[[1]][1]
    print(iso_name)
    if(i==1){
    
    data=read.table(isoforms_results[i],sep="\t",header=T)
    data_count=data[,c(1,2,3,7)]
    colnames(data_count)[4]<-paste(iso_name,"_FPKM",sep="")
    outtab=rbind(outtab,data_count)
    }else{
    data=read.table(isoforms_results[i],sep="\t",header=T)
    data_count=data.frame(data[,c(7)])
    names(data_count)=paste(iso_name,"_FPKM",sep="")
    outtab=cbind(outtab,data_count[,1])
    colnames(outtab)[length(colnames(outtab))]<-paste(iso_name,"_FPKM",sep="")
    }
    }

}



else if (Types=="ISO"){
    outtab=data.frame()
 for(i in 1:length(isoforms_results)){
    iso_name=str_split(isoforms_results[i], "\\.")[[1]][1]
    print(iso_name)
    if(i==1){
    
    data=read.table(isoforms_results[i],sep="\t",header=T)
    data_count=data[,c(1,2,3,8)]
    colnames(data_count)[4]<-paste(iso_name,"_IsoPct",sep="")
    outtab=rbind(outtab,data_count)
    }else{
    data=read.table(isoforms_results[i],sep="\t",header=T)
    data_count=data.frame(data[,c(8)])
    names(data_count)=paste(iso_name,"_IsoPct",sep="")
    outtab=cbind(outtab,data_count[,1])
    colnames(outtab)[length(colnames(outtab))]<-paste(iso_name,"_IsoPct",sep="")
    }
    }

}

return(outtab)
}

list.files(pattern = "*isoforms.results")->isoforms_results

data_tpm<-Rsem_get(isoforms_results,"TPM")





dir() %>% grep("isoforms.results",.) %>% dir()[.]-> isoforms_results
data_count=Rsem_get(isoforms_results,"count")
tail(data_count)

exprset_count=Rsem_get(isoforms_results,"count")
rownames(exprset_count)=exprset_count[,1]
head(exprset_count)

for(i in 1:length(colnames(exprset_count))){
    print(i)
colnames(exprset_count)[i]<-str_split(colnames(exprset_count)[i], "\\.")[[1]][1]
}
colnames(exprset_count)<-gsub("_Aligned","",colnames(exprset_count))
# exprset_count<-exprset_count[,-c(1,2,3)]
group_list=c("TAU_NR","TAU","WT_NR","WT") #就是分成多少组
exprset_count<-exprset_count[,-1]

head(exprset_count)

num=c(4,4,4,4) #这几组每组分别有多少个
num

# rm(group_list)
group_list

num
head(exprset_count)
exprset_count<-ceiling(exprset_count)

library(DESeq2)
Batch_Deseq_differnece<-function(exprSet,group,num,save_dir="Alldiffenece",save_dir2="NEW_MA"){
  ##create a folder 
  save_dir<-paste0(save_dir,"/")
  dir.create(save_dir)
  ## creat a group
  group_list= factor(rep(group_list,times=num))
  group_list
  colData=data.frame(row.names = colnames(exprSet),
                     group=group_list)
  
  #dat<-data.frame()
  ## use the Deseq2 to have Diffence analyse
  for (i in 1:length(group)){
    name=unique(group)[i]
    print(name)
    colData$group<-relevel(colData$group,ref=name)
    dds=DESeq2::DESeqDataSetFromMatrix(countData = exprSet,
                             colData = colData,
                             design = ~group) 
    dds <- dds[ rowSums(DESeq2::counts(dds)) > 10, ]
    dds <- DESeq2::DESeq(dds)
    for (j in 2:length(DESeq2::resultsNames(dds))){
      
    resname=DESeq2::resultsNames(dds)[j]
    
    res=DESeq2::results(dds, name=resname)
    
    res_lfc <- lfcShrink(dds, coef=j, res=res, type="apeglm")
    res_lfc
    #res=res_lfc
    
    summary(res_lfc)
    summary(res)
    
    dir.create(save_dir)
    write.csv(res,paste0(save_dir,resname,".csv"))
    save_dir2=paste0(save_dir2,"/")
    dir.create(save_dir2)
   
    
    
    save_dir_MA=paste0(save_dir2,"/",resname)
    dir.create(save_dir_MA)
    write.csv(res,paste0(save_dir_MA,"/",resname,"_res.csv"))
    write.csv(res_lfc,paste0(save_dir_MA,"/",resname,"_reslfc.csv"))
    png(paste0(save_dir_MA,"/",resname,"_MA.png"),width=600*3,height=3*600,res=72*3) 
    plotMA(res, ylim=c(-3,3),main=paste0(resname," MA"))

    dev.off()
    png(paste0(save_dir_MA,"/",resname,"_MAlfc.png"),width=600*3,height=3*600,res=72*3) 
    xlim <- c(1,1e5); ylim<-c(-3,3)
    plotMA( res_lfc, xlim=xlim, ylim=ylim, main=paste0(resname," apeglm"))

    dev.off()
    
    }
   
  }
  
}
Batch_Deseq_differnece(exprset_count,group=group_list,num,save_dir = "20220414New",save_dir2="202204141536NEW_MA")

