setwd("I:\\HNSC\\")##设置工作路径
rm(list=ls())
library(dplyr)
library(tidyr)
lnc<-read.csv("hg19lnc.txt",sep="\t",header=T)#####读取lncRNA基因列表
lnc2<-lnc %>%
       separate(gene_id,c("ENSEMBL","id"),"[.]")
lnc22<-as.data.frame(lnc2[,5])
colnames(lnc22)<-c("ENSEMBL")

fpkm<-read.csv("TCGA-HNSC.htseq_fpkm.tsv.gz",sep="\t",header=T)#####读取TCGA-HNSC-RNA seq-FPKM数据
fpkm2<-fpkm %>%
       separate(Ensembl_ID,c("ENSEMBL","id"),"[.]")      ####删除ensemble ID的版本号
fpkm2<-fpkm2[,-2]
lnc_fpkm<-merge(lnc22,fpkm2,,by="ENSEMBL")    ####提取RNA-seq表达谱中lncRNA的表达 根据基因ID匹配
write.table(lnc_fpkm,"HNSC_lnc_fpkm.csv",sep=",",col.names=T,row.names=F,quote=F) ###将lncRNA的RNA-seq表达谱写出备用




#######survival####
file<-read.csv("HNSC_lnc_fpkm.csv",sep=",",header=T)
file2<-aggregate(.~Ensembl,data=file,mean)###将同一个lnc的表达值求平均
rownames(file2)<-file2[,1] ###将第一列设置为行名
file22<-file2[,-1]
file222<-as.data.frame(t(file22)) ###将表达谱转置
col<-as.data.frame(row.names(file222))
file2222<-cbind(col,file222)
file2222[1:4,1:4]
colnames(file2222)[1]<-c("sample")



sur<-read.csv("TCGA-HNSC.survival.tsv",sep="\t",header=T)  ####读取TCGA-HNSC-SURVIVAL 生存信息的表格
sur2<-as.data.frame(gsub(sur[,1], pattern = '-', replacement = '.'))  ###将样本名中的“-”转换成“."
sur22<-cbind(sur2,sur)
colnames(sur22)[1]<-c("sample")
sur22<-sur22[,-2]
file3<-merge(sur22,file2222,by="sample") ###根据ensemble id 将lncRNA RNA-seq表达和生存数据合并
library(dplyr)
library(tidyr)
library(stringr)
file33<-file3[which(str_sub(file3[,1],14,15)<10),]  ##根据TCGA样本命名规则，样本ID第14,15位小于10的为cancer样本，提取出cancer样本
write.table(file33,"HNSC_lnc_survival_cancer.csv",sep=",",col.names=T,row.names=F,quote=F)#写出备用




