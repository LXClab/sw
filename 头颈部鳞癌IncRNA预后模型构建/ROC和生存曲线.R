#########三种组合的ROC曲线
##7个基因
rm(list=ls())
setwd("I:\\HNSC\\km+cox+lasso\\")
data<-read.csv("./ROC/Training拟合系数校正-7.csv",sep=",",header=T)##读取系数*表达的模型表
re<-data
##根据risk_score值的中位数分high risk和low risk
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
head(re2)
###绘制ROC曲线
library(pROC)
roc1=roc(re2$status~re2[,3])
plot(roc1,print.auc=TRUE,col="red")
##绘制KM生存曲线
library(survminer) # 加载包
library(survival)
fit <- survfit(Surv(time,status) ~ group,data = re2)
pdf(file="./ROC2/7-survival.pdf",width=5,height=5)
ggsurvplot(fit,
          pval = TRUE, conf.int = TRUE,
          risk.table = TRUE, # 添加风险表
          risk.table.col = "strata", # 根据分层更改风险表颜色
          linetype = "strata", # 根据分层更改线型
          surv.median.line = "hv", # 同时显示垂直和水平参考线
          ggtheme = theme_bw(), # 更改ggplot2的主题
          palette = c("red", "blue"),#定义颜色
          xlab="Time (day)")
dev.off()
###6个基因
data7<-read.csv("./ROC/Training拟合系数校正-6-1.csv",sep=",",header=T)
re7<-data7
op7<-data.frame()
for(i in 1:nrow(re7)){
if(re7[i,3]>median(re7[,3])){
op7[i,1]=c("high")
}else {
op7[i,1]=c("low")
}
}
colnames(op7)<-c("group")
re27<-cbind(re7,op7)
##ROC
library(pROC)
roc7=roc(re27$status~re27[,3])
plot(roc7,add=TRUE,print.auc=TRUE,col="blue")
fit <- survfit(Surv(time,status) ~ group,data = re27)
pdf(file="./ROC2/6-survival.pdf",width=5,height=5)
ggsurvplot(fit,
          pval = TRUE, conf.int = TRUE,
          risk.table = TRUE, # 添加风险表
          risk.table.col = "strata", # 根据分层更改风险表颜色
          linetype = "strata", # 根据分层更改线型
          surv.median.line = "hv", # 同时显示垂直和水平参考线
          ggtheme = theme_bw(), # 更改ggplot2的主题
          palette = c("red", "blue"),#定义颜色
          xlab="Time (day)")
dev.off()

###4个基因
data5<-read.csv("./ROC/Training拟合系数校正-4-1.csv",sep=",",header=T)
re5<-data5
op5<-data.frame()
for(i in 1:nrow(re5)){
if(re5[i,3]>median(re5[,3])){
op5[i,1]=c("high")
}else {
op5[i,1]=c("low")
}
}
colnames(op5)<-c("group")
re25<-cbind(re5,op5)

fit <- survfit(Surv(time,status) ~ group,data = re25)
pdf(file="./ROC2/4-survival.pdf",width=5,height=5)
ggsurvplot(fit,
          pval = TRUE, conf.int = TRUE,
          risk.table = TRUE, # 添加风险表
          risk.table.col = "strata", # 根据分层更改风险表颜色
          linetype = "strata", # 根据分层更改线型
          surv.median.line = "hv", # 同时显示垂直和水平参考线
          ggtheme = theme_bw(), # 更改ggplot2的主题
          palette = c("red", "blue"),#定义颜色
          xlab="Time (day)")
dev.off()


library(pROC)
roc5=roc(data5$status~data5[,3])
plot(roc5,add=TRUE,print.auc=TRUE,col="green")

round(auc(roc1),3)
round(ci(roc1),3)
round(ci(roc7),3)
round(ci(roc5),3)
pdf("ALL_roc.pdf")
plot(roc1,print.auc=TRUE,col="red")
plot(roc7,add=TRUE,print.auc=TRUE,col="blue")
plot(roc5,add=TRUE,print.auc=TRUE,col="green")
dev.off()


