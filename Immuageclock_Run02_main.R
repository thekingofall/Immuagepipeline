library(data.table)
library(tidyverse)
library(rmarkdown)

args=commandArgs(TRUE)


featurecountfile=args[1] #生成的原始表达量文件
coldata=args[2] #样本信息文件
data_load=args[3] #输出文件目录
rmdfile=args[4] #运行rmd文件

getwd()
dir()
# View(featurecountfile)
# featurecountfile=fread("/data2/users/maolp/PBMC/20220317data/New04.Count/NewAll_count.txt")
# coldata=fread("/data2/users/maolp/PBMC/20220317data/NewA05.Count_table/PBMCpartone.txt")
# print(featurecountfile)
# print(coldata)
# print(data_load)
# print(rmdfile)

data = fread(featurecountfile) %>% as.data.frame()
datacol=fread(coldata)
if (length(colnames(data)[7:length(colnames(data))])==nrow(datacol)){
    print("yes")
    
for(i in 1:nrow(datacol)){

    temname=datacol$name[i]
    print(temname)
   
    colnames(data)[grepl(temname,colnames(data))]=temname
    
}
}else{
    print("The coldata file is not the same length as the rawcount file")
}

data %>% separate(col = Geneid,
                  into = c("Geneid"),
                  sep = "\\.") -> data
data_last = data[, c(1, 7:ncol(data))]
# View(data_last)

coldata=datacol
norwayname="tpm"

data3 = data_last


head(data)




### 转换ID
library(clusterProfiler)
library(org.Hs.eg.db)
id = bitr(data3$Geneid,
          fromType = "ENSEMBL",
          toType = "SYMBOL",
          OrgDb = "org.Hs.eg.db")

exprset <- merge(id, data3, by.x = "ENSEMBL", by.y = "Geneid")
exprset <- exprset[which(!duplicated(exprset$SYMBOL)), ]
# head(data_last)
rownames(exprset) <- exprset$SYMBOL
# exprset
# head(exprset)

exprset_name = exprset[, c(1, 2)]
# exprset_name
lastlen = data[, c(1, 6)]
exprset_len <-merge(exprset_name, lastlen, by.x = "ENSEMBL", by.y = "Geneid")
exprset_last <- exprset[, -c(1, 2)] %>% .[rowSums(.) > 10, ]
exprset_len <-exprset_len[exprset_len$SYMBOL %in% rownames(exprset_last), ]
exprset <- exprset[exprset_len$SYMBOL, ]
exprset_last <- exprset_last[,coldata$name ]


library(sva)
library(bladderbatch)
library(DESeq2)
countToTpm <- function(counts, effLen)
    {
      rate <- log(counts) - log(effLen)
      denom <- log(sum(exp(rate)))
      exp(rate - denom + log(1e6))
    }

coldata2=coldata

normlizefunc <-
  function(data_exp=exprset_last, coldata2=coldata,exp_len=exprset,normway = "tpm") {
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
      # print("a")
    
      
      if (all(rownames(data_exp) == exp_len$SYMBOL)) {
        # print("B")
        data_TPM <-
          apply(data_exp, 2, function(x) {
            countToTpm(x, exprset_len$Length)
          }) %>% as.data.frame()
        
  
            
   return(data_TPM)

      }
      
    }

    
    
    
  }

data_norm=normlizefunc(normway = norwayname)

if (all(colnames(data_norm)==coldata$name )){
    print("yes")

 combat_edata <-ComBat(dat = data_norm, batch = coldata$batch) %>% as.data.frame()
}


# dir.create(data_load)
savdir=paste0(data_load,"/","A05.rawCount.Rdata")


save(data,data_last,exprset_len ,exprset,exprset_last,norwayname,coldata,data_norm,combat_edata,file=savdir)
# help(render)
write.table(coldata,paste0(data_load,"/","PBMCcoldata.txt"),sep="\t",quote=F,row.names=F,col.names=T)
write.table(combat_edata,paste0(data_load,"/","A05.combat_eata.txt"),sep="\t",quote=F,row.names=T,col.names=T)
print("hello")
print(paste0("cp -r ",rmdfile," ",data_load))
system(paste0("rm -rf ",data_load,"/","*.Rmd"))
system(paste0("cp -r ",rmdfile," ",data_load))



rmdname=unlist(strsplit(rmdfile, split = "/"))[length(unlist(strsplit(rmdfile, split = "/")))]
rmdname
rmdnow=paste0(data_load,"/",rmdname)
print(rmdnow)


# rmarkdown::render(rmdnow, output_format = "html_document", output_dir = data_load)

# Rscript ~/mao/Codeman/0308PBMCpepline/Immuage_Rmain.R A04.Count/combat_data.csv  
#Rscript ~/mao/Codeman/0308PBMCpepline/Immuage_Rmain.R A04.Count/combat_data.csv  A04.Count/coldata.txt A05test /home/maolp/mao/Codeman/0308PBMCpepline/Immucolck_trandata.Rm