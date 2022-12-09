#########用KM分析中p<0.05的基因做单因素cox回归分析
setwd("I:\\HNSC\\km+cox+lasso\\")
km<-read.csv("HNSC_lnc_KM0.05_cancer.csv",sep=",",header=T)########读取KM分析p<0.05的基因列表
km<-km[-c(1513,1514),]
file3<-read.csv("HNSC_lnc_survival_cancer.csv",sep=",",header=T)##########读取所有的lnc RNA-seq表达谱和生存数据
###匹配出相应的基因的表达谱和生存信息
file33<-file3[,-c(1:4)]
re2<-file33[,match(km[,1],colnames(file33))]
re3<-cbind(file3[,1:4],re2)
write.table(re3,"lnc_KM_cancer_survival.csv",sep=",",col.names=T,row.names=F,quote=F)##读出备用
file3<-read.csv("lnc_KM_cancer_survival.csv",sep=",",header=T)

library(survival)
######cox##3
outTab=data.frame() ###构建新的空表，存储数据
for(i in colnames(file3[,5:ncol(file3)])){ ####循环对基因进行单因素cox回归
    #cox分析
    cox <- coxph(Surv(file3$OS.time, file3$OS) ~ file3[,i], data = file3) ###利用生存状态和生存时间
    coxSummary = summary(cox)
   outTab=rbind(outTab,
                   cbind(id=i,
                         coef=coxSummary$coefficients[,"coef"],
                         HR=coxSummary$conf.int[,"exp(coef)"],
                         HR.95L=coxSummary$conf.int[,"lower .95"],
                         HR.95H=coxSummary$conf.int[,"upper .95"],
                         pvalue=coxSummary$coefficients[,"Pr(>|z|)"]))
    }
colnames(outTab)<-c("geneid","cox_coef","cox_HR","cox_HR95L","cox_HR95H","cox_P")
exout<-outTab[which(outTab$cox_P<0.05),]####挑选p<0.05的基因
dim(exout)
#write.table(outTab,"lnc_km_cox_cancer.csv",sep=",",col.names=T,row.names=F,quote=F)
write.table(exout,"lnc_km_coxp0.05_cancer.csv",sep=",",col.names=T,row.names=F,quote=F)###写出结果



