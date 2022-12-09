#########km cox lasso ROC
rm(list=ls())
setwd("I:\\HNSC\\km+cox+lasso\\")
library(gtools)

####得到对应基因系数### 
xishu<-read.csv("lnc_km_cox_lasso_gene_updown.csv",sep=",",header=T)
dim(xishu)
head(xishu)
colnames(xishu)=c('geneID','coef')
###特征基因的表达及生存
clini_dataTraining<-read.csv("HNSC_lnc_survival_cancer.csv",sep=",",header=T)
clini_dataTraining[1:3,1:5]
dim(clini_dataTraining)
phe<-read.csv("phe.csv",sep=",")####读取癌症临床数据集
clini_dataTraining2<-merge(phe,clini_dataTraining,by="sample")###合并数据
clini_data<-da2[,c(1:7,match(xishu[,1],colnames(da2)))]
sample<-colnames(clini_data)[8:dim(clini_data)[2]]##取基因名  特征基因名
name_location<-match(xishu[,1],colnames(clini_data))###拿到对应系数位置
name_location

for(q in 1:length(sample))
  #q=1###q=2
  #colnames(clini_data)[name_location[q]]
{clini_data[,name_location[q]]<-as.numeric(xishu[q,2])*clini_data[,name_location[q]]} ##筛选出的基因乘以你和系数的表达量
colnames(clini_data)[5]<-c("status")
colnames(clini_data)[7]<-c("time")
clini_data[1:2,1:8]
#####以中位表达值区分高低风险
risk<-data.frame()
for(i in 1:nrow(clini_data)){
risk[i,1]=sum(clini_data[i,8:dim(clini_data)[2]])
}
colnames(risk)<-c("Risk")
group<-data.frame()
for(i in 1:nrow(risk)){
if(risk[i,1]>median(risk[,1])){
group[i,1]=c("Low")
}else {
group[i,1]=c("High")
}
}
colnames(group)<-c("Group")

rda<-cbind(clini_data,risk,group)
library(survival)
library(survminer)
fd<-rda[which(rda$TNM!="Not reported"),] ###剔除没有TNM分期的
fd[1:2,1:dim(fd)[2]]
model<- coxph(Surv(fd$time, fd$status) ~ Age+Gender+TNM+Group, data = fd)


summary(model)
options(scipen=1)
pdf("3个基因基因集和年龄多因素.pdf")
ggforest(model, data =fd, 
         main = "Hazard ratio", 
         cpositions = c(0.10, 0.22, 0.4), 
         fontsize = 1.0, 
         refLabel = "1", noDigits = 4)
dev.off()