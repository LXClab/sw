##########时间独立性ROC检验
###########3个特征基因预测3 年 5年
rm(list=ls())
setwd("I:\\HNSC\\km+cox+lasso\\")
data<-read.csv("./ROC2/Training拟合系数校正-3.csv",sep=",",header=T)
re<-data
colnames(re)[3]<-c("risk_score")
op<-data.frame()
for(i in 1:nrow(re)){
if(re[i,3]>median(re[,3])){
op[i,1]=c("high")
}else {
op[i,1]=c("low")
}
}
colnames(op)<-c("group")
re2<-cbind(re,op)
library(survminer) # 加载包
library(survival) 
library(timeROC)
head(re2)
###time independent
ROC <- timeROC(T=re2$time,   
               delta=re2$status,   
               marker=re2$risk_score,   
               cause=1,                #阳性结局指标数值
               weighting="marginal",   #计算方法，默认为marginal
               times=c(365, 1095, 1825),       #时间点，选取1年，3年和5年的生存率
               iid=TRUE)
pdf("3个基因predict-year.pdf")
plot(ROC, 
     time=365, col="red", lwd=2, title = "")   #time是时间点，col是线条颜色
plot(ROC,
     time=1095, col="blue", add=TRUE, lwd=2)    #add指是否添加在上一张图中
plot(ROC,
     time=1825, col="orange", add=TRUE, lwd=2)
legend("bottomright",
       c(paste0("AUC at 1 year: ",round(ROC[["AUC"]][1],2)), 
         paste0("AUC at 2 year: ",round(ROC[["AUC"]][2],2)), 
         paste0("AUC at 3 year: ",round(ROC[["AUC"]][3],2))),
       col=c("red", "blue", "orange"),
       lty=1, lwd=2,bty = "n") 

dev.off()
