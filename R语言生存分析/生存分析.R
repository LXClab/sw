#===============
#1.生存分析
#===============
# R包准备
#survival 用于计算生存分析
#survminer 用于统计和可视化
# 安装R包
#install.packages(c("survival", "survminer"))
# 加载R包
library("survival")
library("survminer")
#示例数据集survival包内置的lung数据
#inst: 机构代码
#time: 生存时间(天)
#status: 审查状态，1=截尾，2=死亡
#age: 年龄
#sex: Male=1 Female=2
#ph.ecog: ECOG 评分 (0=good 5=dead)
#ph.karno: 医生进行的Karnofsky 评分 (bad=0-good=100)
#pat.karno: 病人自己进行的Karnofsky 评分
#meal.cal: 用餐时摄入的卡路里
#wt.loss: 过去六个月体重减轻
head(lung)
help(lung)
## 计算生存曲线:survfit() 
#可以使用survival 包中的survfit() 计算kaplan-Meier生存估计
fit <- survfit(Surv(time, status) ~ sex, data = lung)
fit
summary(fit)
summary(fit)$table
d <- data.frame(time = fit$time,
                n.risk = fit$n.risk,
                n.event = fit$n.event,
                n.censor = fit$n.censor,
                surv = fit$surv,
                upper = fit$upper,
                lower = fit$lower
)
head(d)
#结果解读
#n：每条曲线中的对象数。
#time：曲线上的时间点。
#n.riks：在时间t处受试者人数
#n.event：在时间t处发生的事件数
#n.censor：在时间t处退出的删失者的数量
#lower,upper：曲线的置信度上限和下限。
#strata：表示曲线估计的分层。如果strata不为NULL，则结果中有多条曲线。strata的水平（一个因子）是曲线的标签。
#可视化
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,#展示置信区间和pvalue
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF")) #颜色
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,#展示置信区间和pvalue
           risk.table = TRUE, # 画出风险表
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF")) #颜色
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  pval = TRUE,conf.int = TRUE,         
  # point estimaes of survival curves.
  conf.int.style = "step",  # customize style of confidence intervals
  xlab = "Time in days",   # customize X axis label.
  break.time.by = 200,     # break X axis in time intervals by 200.
  ggtheme = theme_light(), # customize plot and risk table with a theme.
  risk.table = "abs_pct",  # absolute number and percentage at risk.
  risk.table.y.text.col = T,# colour risk table text annotations.
  risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  surv.median.line = "hv",  # add the median survival pointer.
  legend.labs = 
    c("Male", "Female"),    # change legend labels.
  palette = 
    c("#E7B800", "#2E9FDF") # custom color palettes.
)
#事件发生曲线图
ggsurvplot(fit,
           conf.int = TRUE,
           risk.table.col = "strata", # Change risk table color by groups
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"),
           fun = "event")
## 累计风险概率（cumulative hazard）图
ggsurvplot(fit,
           conf.int = TRUE,
           risk.table.col = "strata", # Change risk table color by groups
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"),
           fun = "cumhaz")
### 使用多个因素的组合来计算生存曲线。
#使用结肠数据集拟合生存曲线
require("survival")
fit2 <- survfit( Surv(time, status) ~ sex + rx + adhere,
                 data = colon )
#基于survminer可视化生存曲线
library(survminer)
# Plot survival curves by sex and facet by rx and adhere
ggsurv <- ggsurvplot(fit2, conf.int = TRUE,
                     ggtheme = theme_bw())
ggsurv$plot +theme_bw() + 
  theme (legend.position = "right")+
  facet_grid(rx ~ adhere)
#===============
#2.log-rank检验
#===============
surv_diff<-survdiff(Surv(time,status)~sex,data=lung)
surv_diff
#============================
#3.cox构建Cox比例风险回归模型
#============================
#单因素cox回归
res.cox <- coxph(Surv(time, status) ~ sex, data = lung)
res.cox
summary(res.cox)
#森林图
ggforest(res.cox, data = lung, 
         main = "Hazard ratio",
         cpositions = c(0.10, 0.22, 0.4),
         fontsize = 1.0) 
#多因素cox回归
res.cox.plus <- coxph(Surv(time, status) ~ age + sex + ph.ecog, data =  lung)
summary(res.cox.plus)
#森林图
ggforest(res.cox.plus, data = lung, 
         main = "Hazard ratio",
         cpositions = c(0.10, 0.22, 0.4),
         fontsize = 1.0) 

#============================
#4.复现SCI论文中的生存分析
#============================
lung_new<-lung[lung$ph.ecog<2,]
fit <- survfit(Surv(time, status) ~ sex + ph.ecog, data = lung_new)
fit
logrank<-survdiff(Surv(time, status) ~ sex + ph.ecog, data = lung_new)
logrank
library(ggplot2)
ggsurv<-ggsurvplot(fit,palette = c("#2E8B57","#800080","#4169E1","#EEE8AA"))
ggsurv
ggsurv$plot+
  labs(x="Time after surgery(month)",y="Overall survival probability")+
  theme(legend.position = "right")
  
#============================
#5.TCGA数据生存分析
#============================
#载入必要的包
#BiocManager::install("cgdsr")
#install.packages("cgdsr")
library("cgdsr")
library("tidyverse")
#获取乳腺癌临床数据
mycgds <- CGDS("http://www.cbioportal.org/")
mycancerstudy <- "brca_tcga"
mycaselist <- getCaseLists(mycgds,mycancerstudy)[1,1] 
myclinicaldata <- getClinicalData(mycgds,mycaselist) 
choose_columns <- c('AJCC_METASTASIS_PATHOLOGIC_PM','AJCC_NODES_PATHOLOGIC_PN','AJCC_PATHOLOGIC_TUMOR_STAGE','AJCC_TUMOR_PATHOLOGIC_PT','AGE','SEX','OS_STATUS','OS_MONTHS',"DFS_MONTHS","DFS_STATUS")
choose_clinicaldata <- myclinicaldata[,choose_columns] #只选择其中的10列用于后续分析
colnames(choose_clinicaldata)[1:4] <- c('M','N','STAGE','T')
###简单生存分析
library(survival)
library(survminer)
library(Cairo)
#choose_clinicaldata<-read.table("TCGA_breast_clinicaldata.txt",header=T,sep="\t")
table(choose_clinicaldata$OS_STATUS)
kmfit1 <- survfit(Surv(OS_MONTHS,OS_STATUS=='1:DECEASED')~1,data=choose_clinicaldata) 
CairoPDF('tcga_surv2.pdf')
ggsurvplot(kmfit1,data = choose_clinicaldata)
dev.off()

