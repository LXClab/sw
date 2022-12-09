setwd("I:\\HNSC\\km+cox+lasso\\")
rm(list=ls())
file3<-read.csv("HNSC_lnc_survival_cancer.csv",sep=",",header=T)###��ȡ��������
file33<-file3[,-c(1:4)]
re1<-read.csv("lnc_km_coxp0.05_cancer.csv",sep=",",header=T)##��ȡcox�ع�p<0.05�Ļ����б�
re2<-file33[,match(re1[,1],colnames(file33))] ###ƥ�������
re3<-cbind(file3[,1:4],re2)
write.table(re3,"lnc_KM_cox_cancer_survival.csv",sep=",",col.names=T,row.names=F,quote=F)##д������
##
rm(list=ls())
file3<-read.csv("lnc_KM_cox_cancer_survival.csv",sep=",",header=T)
library(glmnet)
library(survival)
####
survival_cancer<-file3
x <-as.matrix(survival_cancer[,5:ncol(survival_cancer)]) ##��ȡ������
y <-survival_cancer[,c('OS.time', 'OS')] ##��ȡ����״̬������ʱ��
names(y) <-c('time', 'status')
y$time <-as.double(y$time)
y$status <-as.double(y$status)
y <-as.matrix(survival::Surv(y$time, y$status))
set.seed(13098)##�������������Ե�ģ�Ϳɸ���

fit<- glmnet(x, y, family = "cox", maxit = 10000)
pdf("lnc_km_cox_lasso_cancer24_lambda.pdf")
plot(fit, xvar = "lambda", label = TRUE)
dev.off()

cvfit <- cv.glmnet(x, y, family="cox", maxit = 10000)
pdf("lnc_km_cox_lasso_cancer24_cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()

lasso_fit<-cvfit
lasso_fit$lambda.1se
lasso_fit$lambda.min
coefficient <-coef(lasso_fit, s=lasso_fit$lambda.min)
#���ʣ�µı�����ָ��������λ�ã�
Active.Index <-which(as.numeric(coefficient) != 0)
Active.Index
#�ͱ�����ϵ��
active.coefficients <-as.numeric(coefficient)[Active.Index]
##��Ӧ�ı���
sig_gene_multi_cox <-rownames(coefficient)[Active.Index]
lasso_re<-cbind(gene= sig_gene_multi_cox,coef=active.coefficients)
dim(lasso_re)
head(lasso_re)
colnames(lasso_re)<-c("geneid","lasso_coef")
lasso_re
write.csv(lasso_re,file="lnc_km_cox_lasso_cancer.csv",row.names=F)
##
####�Լ�Ԥ��
###�Լ�Ԥ���Լ���ʽͼ
########lasso�������Ԥ��
lasso.prob<-predict(cvfit,new=x,s=c(cvfit$lambda.min,cvfit$lambda.1se))
re=cbind(y,lasso.prob)
#write.table(re,"./SE/SE_LASSO31_probre.csv",sep=",",col.names=T,row.names=F,quote=F)
###
re=as.data.frame(re)
colnames(re)<-c("time","event","prob_min","prob_1se")
re$event=as.factor(re$event)
library(ggpubr)
p1=ggboxplot(re,x="event",y="prob_min",color="event",
             palette="jco",add="jitter")+stat_compare_means()

p2=ggboxplot(re,x="event",y="prob_1se",color="event",
             palette="jco",add="jitter")+stat_compare_means()

library(patchwork)
pdf("lnc_km_cox_pro_min_boxplot.pdf")
p1
dev.off()
pdf("lnc_km_cox_pro_1se_boxplot.pdf")
p2
dev.off()

###ROC
library(ROCR)
library(pROC)
roc_min=roc(re$event~re[,3])
roc_1se=roc(re$event~re[,4])
pdf("lnc_km_cox_LASSO_ROC.pdf")
plot(roc_min,colorize=FALSE,col="red")
plot(roc_1se,colorize=FALSE,col="blue",add=T)
lines(c(0,1),c(0,1),col="gray",lty=4)
text(0.8,0.3,col="red",
     labels=paste0("pred_min_AUC= ",round(auc(roc_min),3)))####���AUCֵ
text(0.8,0.2,col="blue",
     labels=paste0("pred_1se_AUC= ",round(auc(roc_1se),3)))
dev.off()






