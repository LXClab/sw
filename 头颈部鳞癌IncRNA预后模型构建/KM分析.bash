setwd("I:\\HNSC\\")##设置工作路径
##############logrank
file3<-read.csv("HNSC_lnc_survival_cancer.csv",sep=",",header=T)  ###读取预处理得到的cancer样本中lnc的RNA-seq表达和生存汇总数据
library(survival)
library(pROC)
#给基因表达值分高低
Y2 <- Surv(file3$OS.time,file3$OS)
result_KMAUC <- data.frame(colnames(file3)[-c(1:4)],NA,NA) #新生成一个表，以基因名为列名
colnames(result_KMAUC ) <- c('geneID','KM_P','KM_AUC')
##循环对基因做KM分析
for(i in 5:ncol(file3)){ 

  lisan=ifelse(file3[,i]>median(file3[,i]),1,0) ###根据每个基因的表达分成两组
  if((length(lisan[lisan==0])>1)&(length(lisan[lisan==1])>1)){
    sdf <- survdiff(Y2~lisan,rho=0)     ###KM_P
    result_KMAUC [i,2]=p2=pchisq(sdf$chisq, df=1, lower=FALSE)
    rocR=roc(file3$OS~lisan)
    result_KMAUC [i,3]=rocR$auc
  }
}
table(result_KMAUC$KM_P<0.01) 
table(result_KMAUC$KM_P<0.05) 
result_KMAUC0.05<-result_KMAUC[which(result_KMAUC$KM_P<0.05),] ####提取KM分析结果中p<0.05的
write.table(result_KMAUC0.05,"HNSC_lnc_KM0.05_cancer.csv",sep=",",col.names=T,row.names=F,quote=F)