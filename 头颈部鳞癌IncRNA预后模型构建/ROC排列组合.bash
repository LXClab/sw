rm(list=ls())
setwd("I:\\HNSC\\km+cox+lasso\\")
library(gtools)
######排列组合

xishu<-read.csv("lnc_km_cox_lasso_gene_updown.csv",sep=",",header=T)####得到对应基因系数###
dim(xishu)
head(xishu)
colnames(xishu)=c('geneID','coef')
clini_dataTraining<-read.csv("HNSC_lnc_survival_cancer.csv",sep=",",header=T)###特征基因的表达及生存
clini_data<-clini_dataTraining[,c(1:4,match(xishu[,1],colnames(clini_dataTraining)))] ###匹配基因的表达和生存
colnames(clini_data)
sample<-colnames(clini_data)[5:dim(clini_data)[2]] ##取基因名  特征基因名
name_location<-match(xishu[,1],colnames(clini_data))###拿到对应系数位置
name_location
#####计算单独基因的risk——score
for(q in 1:length(sample))
  #q=1###q=2
  #colnames(clini_data)[name_location[q]]
{clini_data[,name_location[q]]<-as.numeric(xishu[q,2])*clini_data[,name_location[q]]} ##筛选出的基因乘以你和系数的表达量
clini_data[1:3,1:5]
colnames(clini_data)[2]<-c("status")
colnames(clini_data)[4]<-c("time")
clini_data[1:4,1:5]
write.csv(clini_data[,-c(1,3)],"./ROC/Training拟合系数校正-1.csv",row.names=F)


##将基因集随机组合，生成risk score表
n=length(sample)###从几个里选
v=sample
for(p in 2:n)
{
  y=p
  com<-combinations(n,y,v)
  ####choose(10,3)
  cox_names<-function(com,y)
  {
    multi_names1<-""
    for(u in 1:y){multi_names1<-paste0(multi_names1,"_",com[u]) }
    multi_names1<-substr(multi_names1,2,nchar(multi_names1))
    (multi_names1)
  }
  multi_names<-apply(com,1,cox_names,y)####对应的基因名字连一块##length(multi_names)
  ###组合找出多基因序列名####combination(n=length(sample),y=3,v=sample)###排列
  newcox<-clini_data[,c(4,2)] ###提取生存时间和状态
  for(i in 1:choose(length(sample),y))
  {
    location<-match(com[i,],colnames(clini_data))
    newcox1<-apply(clini_data[,location],1,sum)####按行求和
    newcox<-cbind(newcox,newcox1)####按行求和
  }
  colnames(newcox)<-c("time","status",multi_names)
  ###newcox[1:3,1:4]###dim(newcox)
  write.csv(newcox,paste0("./ROC/Training拟合系数校正-",p,".csv"),row.names=F)
print(p)
}



###############Train组筛选
#######要与上setwd()路径一致#######


#########计算每种组合的的AUC值
library(pROC)
library(survival)
setwd('./ROC')
mfile<-dir()
resultTraining<-data.frame()
AUC_MAX_every<-data.frame()
gene_number<-length(colnames(clini_data)[-c(1:4)])
gene_number
for(f in 1:gene_number)
{##f=1
  ###############ROC
  mayo<-read.csv(mfile[f])
  endvalue <- length(mayo[1,])
  auc<-c()
  for(j in 3:endvalue)
  {
    #lisan=ifelse(mayo[,j]>median(mayo[,j]),1,0)
    rocR=roc(mayo$status~mayo[,j])
    auc<-c(auc,rocR$auc)
  }
  a<-unlist(strsplit(mfile[f],'-'))[2]
  b<-as.numeric(strsplit(a,'.csv'))
  AUC_MAX_every<-rbind(AUC_MAX_every,t(data.frame(c(mfile[f],b,max(auc)))))
  result_AUC<-cbind(colnames(mayo[3:endvalue]),auc)
  colnames(result_AUC)<-c("gene","AUC")
  #########################################################生存分析KM####
  stad<-mayo
  endvalue_KM<-length(stad[1,])
  KM_P<-c()
  COX_xishu<-c()
  COX_P<-c()
  for(i in 3:endvalue_KM)
  {
    Y2 <- Surv(stad$time,stad$status)
    gene1<-as.double(stad[,i])
    gene2<-ifelse(gene1>=median(gene1),1,0)
    if((length(gene2[gene2==0])!=0)&(length(gene2[gene2==1])!=0))
    {
      survfit(Y2~gene2)
      sdf <- survdiff(Y2~gene2,rho=0)
      p2 <- pchisq(sdf$chisq, df=1, lower=FALSE)
      res <- coxph(Surv(stad$time,stad$status)~gene2)##多基因就+as.numeric(stad[,i])+a+b逐步回归
      COX_xishu<- c(COX_xishu,res$coef)  #取回归系数
      COX_P<- c(COX_P,summary(res)$coefficients[5])  #取回归系数
      KM_P<-c(KM_P,p2)
    }
    else
    {
      COX_xishu<-c(COX_xishu,2)
      KM_P<-c(KM_P,2)}
  }
  result_KM_P<-data.frame(cbind(colnames(mayo[3:endvalue]),KM_P,COX_xishu,COX_P))
  colnames(result_KM_P)<-c("gene","KM_P","COX_xishu","COX_P")
  result_KM_P$Num=b
  rownames(AUC_MAX_every)<-NULL
  #####################################合并AUC与KM写出##################
  result1<-cbind(result_AUC,result_KM_P)
  #write.table(result1,paste0(f,".csv",sep=""),sep=",",row.names=F,quote=F)
  resultTraining<-rbind(resultTraining,result1)
  print(i)
print(f)}

projectname=20221124
write.csv(resultTraining,paste(projectname,'-',"Training-组合合并-ROC-KM.csv"),row.names=F)
AUC_MAX_every[,2]<-as.numeric(as.character(AUC_MAX_every[,2]))###提取组合数
AUC_MAX_every[,3]<-as.numeric(as.character(AUC_MAX_every[,3]))###提取每种组合的最大auc值
AUC_MAX_every<-AUC_MAX_every[order(AUC_MAX_every[,2],decreasing = F),]
write.csv(AUC_MAX_every,paste(projectname,'-',"Train-组合合并-AUC_MAX_every.csv"),row.names=F)
pdf(paste(projectname,'-',"AUC_MAX-图.pdf"))
plot(AUC_MAX_every[,2],AUC_MAX_every[,3],type="b")
dev.off()


max(AUC_MAX_every[,3])

AUC_MAX_everyTrain=AUC_MAX_every

