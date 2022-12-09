##3个基因的ROC曲线
rm(list=ls())
setwd("I:\\HNSC\\km+cox+lasso\\")
data<-read.csv("./ROC2/Training拟合系数校正-3.csv",sep=",",header=T) ###读取AUC最大时的基因模型
re<-data
op<-data.frame()
###根据中位表达值区分高低风险
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
library(pROC)
roc1=roc(re2$status~re2[,3])
##生存曲线
library(survminer) # 加载包
library(survival)
fit <- survfit(Surv(time,status) ~ group,data = re2)
pdf(file="./ROC2/3-survival.pdf",width=5,height=5)
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
########3ROC曲线
pdf("ALL_roc.pdf")
plot(roc1,print.auc=TRUE,col="red")

dev.off()