#########################################
library("rtracklayer")
gtf_data = import('./ref/gencode.v41.long_noncoding_RNAs.gtf') #gtf的路径
#这里使用import导入gtf文件， 生成一个GRangs对象
gtf_data = as.data.frame(gtf_data)
gtf<-gtf_data[,c(1,2,3,7,10,11,12)]
gtf2<-gtf[which(gtf$type=="gene"),]
write.table(gtf2,"./hg19lnc.txt",sep="\t",col.names=T,row.names=F,quote=F)