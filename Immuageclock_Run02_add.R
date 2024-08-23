
## 整合新样本和旧的样本的脚本
rm(list = ls())
options (warn = -1)





args = commandArgs(TRUE)

Newfeaturedata=args[1] #新生成的原始表达量文件 # nolint
Oldfeaturedata=args[2] #以前的原始表达量文件
coldata_newarg=args[3] #新样本信息文件
coldata_oldarg=args[4] #新样本信息文件 # nolint
data_newloaddir=args[5] #输出文件目录
runmode=args[6] #运行模式,1为全部，2为更新
runname=args[7] #运行名称
norwayname=args[8] #tpm
oldsave=args[9] #旧的输出文件目录

print("___________start___________")



# exp_raw_merge=cbind(exp_raw_old,exp_raw_new)

if (length(args)<9) {
 
  stop("Usage: <Newfeaturedata> <Oldfeaturedata> <coldata_newarg> <coldata_oldarg> 
  <data_newloaddir> <runmode>['cover','single'] <runname> <norwayname>['deseq','tpm'] <oldsave>['Immuageroot']")
}else{


suppressMessages(library(data.table))
suppressMessages(library(tidyverse,quietly =TRUE))


exp_raw_old= read.table(Oldfeaturedata,header = T)


exp_raw_new=read.table(Newfeaturedata,header = T)

# print(head(exp_raw_old)[,1])

# print(head(exp_raw_new)[,1])

# if (exp_raw_new$Geneid == exp_raw_old$Geneid){
# exp_raw_new=exp_raw_new[,c(7:ncol(exp_raw_new))]
# exp_raw_merge=cbind(exp_raw_old,exp_raw_new)
# }else{
#   print("Geneid not match")
# }

print("___________load___________")
  
data_ten=suppressMessages(fread(Newfeaturedata,header = T) %>% as.data.frame() %>%
    separate(col = Geneid,
             into = c("Geneid"),
             sep = "\\."))


             
exp_new=suppressMessages(fread(Newfeaturedata,header = T) %>% as.data.frame() %>%
    separate(col = Geneid,
             into = c("Geneid"),
             sep = "\\."))
# exp_new
print("___________load olddata___________")

exp_new=exp_new[, c(1, 7:ncol(exp_new))]
exp_old=suppressMessages(fread(Oldfeaturedata) %>% as.data.frame()  %>%
    separate(col = Geneid,
             into = c("Geneid"),
             sep = "\\."))


# print(all(exp_new$Geneid==exp_old$Geneid)


# exp_raw_old= read.table(Oldfeaturedata,header = T)


# exp_raw_new=read.table(Newfeaturedata,header = T)

# print(head(exp_raw_old)[,1])

# print(head(exp_raw_new)[,1])

# if (exp_new$Geneid == exp_old$Geneid){
#   print("Geneid not match")
# }else {
#    exp_raw_new=exp_new[,c(7:ncol(exp_new))]
#   #  exp_raw_merge=cbind(exp_raw_old,exp_raw_new)

# ···exp_raw_merge=cbind(exp_raw_old,exp_raw_new)
# }

# help(stop)

print("__________part1:deal with data___________")
# exp_old
exp_old=exp_old[, c(1, 7:ncol(exp_old))]

exp_last=merge(exp_new,exp_old,by="Geneid",all=T)


print("__________part2.1: merge the data___________")
# print(coldata_newarg)
coldata_new=read.table(coldata_newarg,sep = "\t",header=T)
print(coldata_oldarg)

# head(coldata_new)

# coldatatest=read.table()

coldata_old=read.table(coldata_oldarg,sep = "\t",header=T)
coldata_old$real_name=coldata_old$name

print(coldata_old)
# coldata_old$batch=rep(1,nrow(coldata_old))
# head(coldata_old)
bacthnum=as.numeric(length(unique(coldata_old$batch)))+1

print(coldata_new)
coldata_new$batch=rep(bacthnum,nrow(coldata_new))

coldata_last=rbind(coldata_new,coldata_old)
# coldata_last
# dim(coldata_last)
# dim(exp_last)
print("__________part2.2: change the name___________")

if (length(colnames(exp_last)[2:length(colnames(exp_last))])==nrow(coldata_last)){
    print("the colnames is right")
    
for(i in 1:nrow(coldata_last)){

    temname=coldata_last$name[i]
    #print(temname)
   
    colnames(exp_last)[grepl(temname,colnames(exp_last))]=temname
    
}
}else{
    stop("The coldata file is not the same length as the rawcount file")
}
data3=exp_last

# colnames(data3)[2]<-"OT4"
### 转换ID
suppressMessages(library(clusterProfiler,quietly =TRUE))
suppressMessages(library(org.Hs.eg.db,quietly =TRUE))
id =suppressMessages( bitr(data3$Geneid,
          fromType = "ENSEMBL",
          toType = "SYMBOL",
          OrgDb = "org.Hs.eg.db"))

exprset <- merge(id, data3, by.x = "ENSEMBL", by.y = "Geneid")
exprset <- exprset[which(!duplicated(exprset$SYMBOL)), ]

rownames(exprset) <- exprset$SYMBOL


exprset_name = exprset[, c(1, 2)]

lastlen = data_ten[, c(1, 6)]
exprset_len <-merge(exprset_name, lastlen, by.x = "ENSEMBL", by.y = "Geneid")
# head(exprset_len)
# data$length
exprset_last <- exprset[, -c(1, 2)] %>% .[rowSums(.) > 10, ]# 去掉第一列和第二列
exprset_len <-exprset_len[exprset_len$SYMBOL %in% rownames(exprset_last), ]
exprset <- exprset[exprset_len$SYMBOL, ] 
# print(coldata2$name)
# print(colnames(exprset_last))
exprset_last <- exprset_last[,coldata_last$name ]

print("__________part2:Normalize the data__________")
suppressMessages(library(sva,quietly = T))
suppressMessages(library(bladderbatch,quietly = T))
suppressMessages(library(DESeq2,quietly = T))
countToTpm <- function(counts, effLen){
      rate <- log(counts) - log(effLen)
      denom <- log(sum(exp(rate)))
      exp(rate - denom + log(1e6))
    }

coldata2=coldata_last

suppressMessages(library(tidyverse,quietly = T))

normlizefunc <-
  function(data_exp=exprset_last, coldata2=coldata,exp_len=exprset_len ,normway = "tpm") {
    if (normway == "deseq") {
      coldata2$batch<-factor(coldata2$batch)
      dds <- DESeqDataSetFromMatrix(countData = data_exp,
                                    colData = coldata2,
                                    design = ~ batch)
      dds <- DESeq(dds) #标准化
      rld <-rlogTransformation(dds)  
      exprSet_new = assay(rld) %>% as.data.frame()
      return(exprSet_new)
    }else if(normway == "tpm") {
      if (all(rownames(data_exp) == exp_len$SYMBOL)) {
  
        # print("B")
        data_TPM <-
          apply(data_exp, 2, function(x) {
            countToTpm(x,exp_len$Length)
          }) %>% as.data.frame()
        
  
            
   return(data_TPM)

      }
      
    }


    
    
    
  }





data_norm=normlizefunc(normway = norwayname)
# data_norm=log2(data_norm+1)
if (all(colnames(data_norm)==coldata2$name )){
    print("the colnames is the same")

 combat_edata <-ComBat(dat = data_norm, batch = coldata_last$batch) %>% as.data.frame()
}


print("__________part3:save the data_________")

# savename=unlist(strsplit(data_newloaddir,"_"))[[1]]

# print("test")
# dir.create(data_newloaddir)
savdir=paste0(data_newloaddir,"/",runname,"_A05.rundata.Rdata")

if (length(coldata_new$name)==1){

  singeldata=combat_edata[,colnames(combat_edata) %in% coldata_new$name] %>% as.data.frame()
  colnames(singeldata)[1]<-coldata_new$name
  rownames(singeldata)=rownames(combat_edata)
  
}else {
   singeldata=combat_edata[,c("YM3","YF1" )] %>% as.data.frame()
}

if (runmode=="single"){
# head(singeldata)

save(data,exprset_len ,exprset,exprset_last,norwayname,coldata_last,data_norm,combat_edata,file=savdir)
# help(render)
write.table(coldata_new,paste0(data_newloaddir,"/",runname,"_humancoldata.txt"),sep="\t",quote=F,row.names=F,col.names=T)
write.table(coldata_last,paste0(data_newloaddir,"/",runname,"_mergecoldata.txt"),sep="\t",quote=F,row.names=F,col.names=T)
write.table(singeldata,paste0(data_newloaddir,"/",runname,"_PBMCsingledata.txt"),sep="\t",quote=F,row.names=T,col.names=T)
write.table(combat_edata,paste0(data_newloaddir,"/",runname,"_A05.combat_eata.txt"),sep="\t",quote=F,row.names=T,col.names=T)



}else if(runmode=="cover"){

feature=rownames(combat_edata)
write.table(coldata_new,paste0(data_newloaddir,"/",runname,"_humancoldata.txt"),sep="\t",quote=F,row.names=F,col.names=T)
write.table(coldata_last,paste0(data_newloaddir,"/",runname,"_mergecoldata.txt"),sep="\t",quote=F,row.names=F,col.names=T)
write.table(singeldata,paste0(data_newloaddir,"/",runname,"_PBMCsingledata.txt"),sep="\t",quote=F,row.names=T,col.names=T)
write.table(combat_edata,paste0(data_newloaddir,"/",runname,"_A05.combat_eata.txt"),sep="\t",quote=F,row.names=T,col.names=T)

write.table(coldata_last,paste0(oldsave,"/Root_A04.Count/","Root_humancoldata.txt"),sep="\t",quote=F,row.names=F,col.names=T)
write.table(feature,paste0(oldsave,"/Root_A05.Count_table/","Root_humanfeature.txt"),sep="\t",quote=F,row.names=F,col.names=T)
#write.table(singeldata,paste0(data_newloaddir,"/""PBMCsingeldata.txt"),sep="\t",quote=F,row.names=F,col.names=T)
write.table(combat_edata,paste0(oldsave,"/Root_A05.Count_table/","Root_A05.combat_edata.txt"),sep="\t",quote=F,row.names=T,col.names=T)

write.table(exp_raw_merge,paste0(data_newloaddir,"/Root_A04.Count/"," Root_Allraw.count.txt"),sep="\t",quote=F,row.names=F,col.names=T)


}

# data=read.table("/home/maolp/data2/Immu_age/test/Root_A05.Count_table/Root_PBMCsingledata.txt",header=T)


}


print("__________add:end_________")

# Newfeaturedata="Root_A04.Count/RootAll_count.txt"
# Oldfeaturedata="../Immuageroot/Root_A04.Count/NewAll_count.txt"
# coldata_newarg="Coldata/coldata.txt"
# coldata_oldarg="../Immuageroot/Root_A05.Count_table/PBMCcoldata.txt"
# data_newloaddir="Root_A05.Count_table"

# str_splitfun=function(x,y){
#   return(unlist(strsplit(x,split=y)))
# }
# str_splitfun(Oldfeaturedata,"/")[length(str_splitfun(Oldfeaturedata,"/"))]
# oldsave=paste(str_splitfun(Oldfeaturedata,"/")[-length(str_splitfun(Oldfeaturedata,"/"))], sep = "/")
