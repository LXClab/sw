#####lasso8个基因的表达热图 
diff<-read.csv("E:\\汕大\\TCGA-expression\\diff\\TCGA-GTEX-stomach-count-diff.csv",sep=",",header=T)
se<-read.csv("lnc_km_cox_lasso_cancer.csv",sep=",",header=T)
re<-merge(se,diff,by.x="geneid",by.y="Ensembl")
write.table(re,"lnc_km_cox_lasso_diff.csv",sep=",",col.names=T,row.names=F,quote=F)
ree<-re[which(re[,9]=="up"|re[,9]=="down"),]
write.table(ree,"lnc_km_cox_lasso_diff_updown.csv",sep=",",col.names=T,row.names=F,quote=F)

###
count<-read.csv("E:\\汕大\\TCGA-expression\\tcga_gtex_stomach_count.csv",sep=",",header=T)
row<-as.data.frame(rownames(count))
count2<-cbind(row,count)
colnames(count2)[1]<-c("geneid")

laex<-merge(se,count2,by="geneid")
write.table(laex,"lnc_km_cox_lasso_tcga_gtex_count.csv",sep=",",col.names=T,row.names=F,quote=F)
laex<-read.csv("lnc_km_cox_lasso_tcga_gtex_count.csv",sep=",",header=T)
laex2<-laex[,-2]
rownames(laex2)<-laex2[,1]
laex2<-laex2[,-1]
library(dplyr)
library(tidyr)
library(stringr)
#tumor<-laex2[,which(str_sub(colnames(laex2),14,15)<10)]
#normal<-laex2[,which(str_sub(colnames(laex2),14,15)>10)]
red<-laex2
annotation_col = data.frame(
  group = rep(c("Normal","Tumor"),c(206,375)),
  row.names = colnames(red))
head(annotation_col)

library(pheatmap)
pheatmap(red,
         scale = "column",cluster_cols =FALSE,show_colnames=F,
         color = colorRampPalette(colors = c("blue","white","red"))(100),
         annotation_col = annotation_col,
         cellwidth =1, cellheight = 5)
##up down hetamap
laex<-merge(ree,count2,by="geneid")
write.table(laex[,-c(2:9)],"lnc_km_cox_lasso_tcga_gtex_count_updown.csv",sep=",",col.names=T,row.names=F,quote=F)
laex<-read.csv("lnc_km_cox_lasso_tcga_gtex_count_updown.csv",sep=",",header=T)
laex2<-laex
rownames(laex2)<-laex2[,1]
laex2<-laex2[,-1]
library(dplyr)
library(tidyr)
library(stringr)
#tumor<-laex2[,which(str_sub(colnames(laex2),14,15)<10)]
#normal<-laex2[,which(str_sub(colnames(laex2),14,15)>10)]
red<-laex2
annotation_col = data.frame(
  group = rep(c("Normal","Tumor"),c(206,375)),
  row.names = colnames(red))
head(annotation_col)

library(pheatmap)
pheatmap(red,
         scale = "column",cluster_cols =FALSE,show_colnames=F,
         color = colorRampPalette(colors = c("blue","white","red"))(100),
         annotation_col = annotation_col,
         cellwidth =1, cellheight = 10)