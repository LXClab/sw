rm(list=ls())
setwd("I:\\HNSC\\km+cox+lasso\\")
library(gtools)
######�������

xishu<-read.csv("lnc_km_cox_lasso_gene_updown.csv",sep=",",header=T)####�õ���Ӧ����ϵ��###
dim(xishu)
head(xishu)
colnames(xishu)=c('geneID','coef')
clini_dataTraining<-read.csv("HNSC_lnc_survival_cancer.csv",sep=",",header=T)###��������ı�Ｐ����
clini_data<-clini_dataTraining[,c(1:4,match(xishu[,1],colnames(clini_dataTraining)))] ###ƥ�����ı�������
colnames(clini_data)
sample<-colnames(clini_data)[5:dim(clini_data)[2]] ##ȡ������  ����������
name_location<-match(xishu[,1],colnames(clini_data))###�õ���Ӧϵ��λ��
name_location
#####���㵥�������risk����score
for(q in 1:length(sample))
  #q=1###q=2
  #colnames(clini_data)[name_location[q]]
{clini_data[,name_location[q]]<-as.numeric(xishu[q,2])*clini_data[,name_location[q]]} ##ɸѡ���Ļ���������ϵ���ı����
clini_data[1:3,1:5]
colnames(clini_data)[2]<-c("status")
colnames(clini_data)[4]<-c("time")
clini_data[1:4,1:5]
write.csv(clini_data[,-c(1,3)],"./ROC/Training���ϵ��У��-1.csv",row.names=F)


##�����������ϣ�����risk score��
n=length(sample)###�Ӽ�����ѡ
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
  multi_names<-apply(com,1,cox_names,y)####��Ӧ�Ļ���������һ��##length(multi_names)
  ###����ҳ������������####combination(n=length(sample),y=3,v=sample)###����
  newcox<-clini_data[,c(4,2)] ###��ȡ����ʱ���״̬
  for(i in 1:choose(length(sample),y))
  {
    location<-match(com[i,],colnames(clini_data))
    newcox1<-apply(clini_data[,location],1,sum)####�������
    newcox<-cbind(newcox,newcox1)####�������
  }
  colnames(newcox)<-c("time","status",multi_names)
  ###newcox[1:3,1:4]###dim(newcox)
  write.csv(newcox,paste0("./ROC/Training���ϵ��У��-",p,".csv"),row.names=F)
print(p)
}



###############Train��ɸѡ
#######Ҫ����setwd()·��һ��#######


#########����ÿ����ϵĵ�AUCֵ
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
  #########################################################�������KM####
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
      res <- coxph(Surv(stad$time,stad$status)~gene2)##������+as.numeric(stad[,i])+a+b�𲽻ع�
      COX_xishu<- c(COX_xishu,res$coef)  #ȡ�ع�ϵ��
      COX_P<- c(COX_P,summary(res)$coefficients[5])  #ȡ�ع�ϵ��
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
  #####################################�ϲ�AUC��KMд��##################
  result1<-cbind(result_AUC,result_KM_P)
  #write.table(result1,paste0(f,".csv",sep=""),sep=",",row.names=F,quote=F)
  resultTraining<-rbind(resultTraining,result1)
  print(i)
print(f)}

projectname=20221124
write.csv(resultTraining,paste(projectname,'-',"Training-��Ϻϲ�-ROC-KM.csv"),row.names=F)
AUC_MAX_every[,2]<-as.numeric(as.character(AUC_MAX_every[,2]))###��ȡ�����
AUC_MAX_every[,3]<-as.numeric(as.character(AUC_MAX_every[,3]))###��ȡÿ����ϵ����aucֵ
AUC_MAX_every<-AUC_MAX_every[order(AUC_MAX_every[,2],decreasing = F),]
write.csv(AUC_MAX_every,paste(projectname,'-',"Train-��Ϻϲ�-AUC_MAX_every.csv"),row.names=F)
pdf(paste(projectname,'-',"AUC_MAX-ͼ.pdf"))
plot(AUC_MAX_every[,2],AUC_MAX_every[,3],type="b")
dev.off()


max(AUC_MAX_every[,3])

AUC_MAX_everyTrain=AUC_MAX_every

