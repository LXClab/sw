#########��KM������p<0.05�Ļ�����������cox�ع����
setwd("I:\\HNSC\\km+cox+lasso\\")
km<-read.csv("HNSC_lnc_KM0.05_cancer.csv",sep=",",header=T)########��ȡKM����p<0.05�Ļ����б�
km<-km[-c(1513,1514),]
file3<-read.csv("HNSC_lnc_survival_cancer.csv",sep=",",header=T)##########��ȡ���е�lnc RNA-seq����׺���������
###ƥ�����Ӧ�Ļ���ı���׺�������Ϣ
file33<-file3[,-c(1:4)]
re2<-file33[,match(km[,1],colnames(file33))]
re3<-cbind(file3[,1:4],re2)
write.table(re3,"lnc_KM_cancer_survival.csv",sep=",",col.names=T,row.names=F,quote=F)##��������
file3<-read.csv("lnc_KM_cancer_survival.csv",sep=",",header=T)

library(survival)
######cox##3
outTab=data.frame() ###�����µĿձ��洢����
for(i in colnames(file3[,5:ncol(file3)])){ ####ѭ���Ի�����е�����cox�ع�
    #cox����
    cox <- coxph(Surv(file3$OS.time, file3$OS) ~ file3[,i], data = file3) ###��������״̬������ʱ��
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
exout<-outTab[which(outTab$cox_P<0.05),]####��ѡp<0.05�Ļ���
dim(exout)
#write.table(outTab,"lnc_km_cox_cancer.csv",sep=",",col.names=T,row.names=F,quote=F)
write.table(exout,"lnc_km_coxp0.05_cancer.csv",sep=",",col.names=T,row.names=F,quote=F)###д�����



