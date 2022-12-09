########测试集检验
rm(list=ls())
setwd("I:\\HNSC\\km+cox+lasso\\")
library(gtools)
####得到3个基因对应基因系数###
xishu<-read.csv("lnc_km_cox_lasso_gene_updown.csv",sep=",",header=T)
xishu<-xishu[c(2,4,5),]
dim(xishu)
head(xishu)
colnames(xishu)=c('geneID','coef')
###测试集基因集
clini_dataTraining<-read.csv("HNSC_lnc_survival_cancer.csv",sep=",",header=T)
da1<-clini_dataTraining[,1:4]
da2<-clini_dataTraining[,-c(1:4)]
da22<-as.data.frame(sample(da2,10000))####随机提取10000个样本作为测试集
da3<-da22[,match(xishu[,1],colnames(da22))]
da<-cbind(da1,da3)####eight genes clinical expression(only tumor)
clini_data<-da
sample<-colnames(clini_data)[5:dim(clini_data)[2]]##取基因名  特征基因名
name_location<-match(xishu[,1],colnames(clini_data))###拿到对应系数位置
name_location
for(q in 1:length(sample))
  #q=1###q=2
  #colnames(clini_data)[name_location[q]]
{clini_data[,name_location[q]]<-as.numeric(xishu[q,2])*clini_data[,name_location[q]]}##筛选出的基因乘以你和系数的表达量
colnames(clini_data)[2]<-c("status")
colnames(clini_data)[4]<-c("time")
#write.csv(clini_data,"./9-GSE84437-survival.txt",sep="\t",row.names=F,quote=F)
##
data<-clini_data
vec<-data.frame()
for(i in 1:nrow(data)){
vec[i,1]= sum(data[i,5:dim(data)[2]])
}
colnames(vec)<-c("risk_score")
re<-cbind(data,vec)
#write.table(re,"./7-GSE84437-survival.txt",sep="\t",row.names=F,quote=F)

op<-data.frame()
for(i in 1:nrow(re)){
if(re[i,dim(re)[2]]>median(re[,dim(re)[2]])){
op[i,1]=c("high")
}else {
op[i,1]=c("low")
}
}
table(op)
colnames(op)<-c("group")
re2<-cbind(re,op)
table(re2$group)
library(survminer) # 加载包
library(survival) 
fit <- survfit(Surv(time,status) ~ group,data = re2)
pdf(file="./test/test-survival.pdf",width=5,height=5)
ggsurvplot(fit,
          pval = TRUE, conf.int = TRUE,
          risk.table = TRUE, # 添加风险表
          risk.table.col = "strata", # 根据分层更改风险表颜色
          linetype = "strata", # 根据分层更改线型
          surv.median.line = "hv", # 同时显示垂直和水平参考线
          ggtheme = theme_bw(), # 更改ggplot2的主题
          palette = c("red", "blue"),#定义颜色
          xlab="Time (Day)")
dev.off()

