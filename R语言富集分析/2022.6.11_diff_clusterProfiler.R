#差异分析
rm(list=ls())
setwd("C:/Users/34750/Desktop/R教程-富集分析")
library(edgeR)
exprSet_1 <- read.csv("normal_expression.csv", header=T, stringsAsFactors = FALSE)
exprSet_2 <- read.csv("GH_expression.csv", header=T, stringsAsFactors = FALSE)
exprSet <- merge(exprSet_1,exprSet_2,by="Gene")
group_list <- factor(c(rep("normal",7), rep("GH",11)))
##设置分组信息，去除低表达量的gene以及做TMM标准化
rownames(exprSet)<-exprSet[,1]
exprSet<-exprSet[,-1]
exprSet <- exprSet[rowSums(cpm(exprSet) > 1) >= 2,]
exprSet <- DGEList(counts = exprSet, group = group_list)
exprSet <- calcNormFactors(exprSet)
##使用qCML（quantile-adjusted conditional maximum likelihood）估计离散度（只针对单因素实验设计）
exprSet <- estimateCommonDisp(exprSet)
exprSet <- estimateTagwiseDisp(exprSet)
##寻找差异gene(这里的exactTest函数还是基于qCML并且只针对单因素实验设计)，然后按照阈值进行筛选即可
et <- exactTest(exprSet)
tTag <- topTags(et, n=nrow(exprSet))
tTag <- as.data.frame(tTag)

tTag[which(tTag$PValue < 0.05 & tTag$logFC <= -1  %in% NA),'sig'] <- 'Down'
tTag[which(tTag$PValue < 0.05 & tTag$logFC >= 1   %in% NA),'sig'] <- 'Up'
tTag[which(tTag$PValue >= 0.05| abs(tTag$logFC) < 1  %in% NA),'sig'] <- 'None'
head(tTag)
table(tTag$sig)
write.csv(tTag,file = "normal-GH.csv")
tTag$PValue<--log(tTag$PValue,base=10)
colnames(tTag)[3]<-c("log10PValue")

#绘制火山图
library(ggplot2)
ggplot(tTag, aes(x = logFC, y = log10PValue, color = sig)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_colour_manual(values  = c('red2', 'green2', 'gray'), 
                      limits = c('Up', 'Down', 'None')) + #更改颜色
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', 
                                        fill = 'transparent'), 
        plot.title = element_text(hjust = 0.5)) + #更改背景框
  theme(legend.key = element_rect(fill = 'transparent'), 
        legend.background = element_rect(fill = 'transparent'), 
        legend.position = c(0.88, 0.9)) + #更改图例位置
  xlim(-6, 6) + ylim(0, 6) +
  labs(x = '\nLog2 Fold Change', 
       y = '-log10(p_value)\n', 
       color = '') #更改坐标轴范围和横纵坐标的标题


#富集分析
#BiocManager::install("clusterProfiler")
#BiocManager::install("GO.db")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("topGO")
#BiocManager::install("pathview")
###############
# libraries
library(ggplot2)
library(gridExtra)
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(DO.db)
library(R.utils)
R.utils::setOption("clusterProfiler.download.method","auto")

#通过keytypes函数可查看org.Hs.eg.db所有支持及可转化类型
keytypes(org.Hs.eg.db) 
setwd("C:/Users/34750/Desktop/R教程-富集分析")
data<-read.csv("genelist.csv",header = F)
x <- unique(data[,1])
eg <- bitr(x, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")

genelist<-eg[,2]
genelist<-unique(genelist)
genelist[duplicated(genelist)]##检验是否还有重复值
length(genelist)
##go
go <- enrichGO(genelist, OrgDb = org.Hs.eg.db, ont='ALL',
               pvalueCutoff = 0.05,qvalueCutoff = 0.2,
               keyType = 'ENTREZID')
dim(go)
dim(go[go$ONTOLOGY=='BP',])
dim(go[go$ONTOLOGY=='CC',])
dim(go[go$ONTOLOGY=='MF',])
# 看来这些差异基因主要富集到BP中了
go.res <- data.frame(go)
goBP <- subset(go.res,subset = (ONTOLOGY == "BP"))[1:5,]
goCC <- subset(go.res,subset = (ONTOLOGY == "CC"))[1:5,]
goMF <- subset(go.res,subset = (ONTOLOGY == "MF"))[1:5,]
go.df <- rbind(goBP,goCC,goMF)
go.df$Description <- factor(go.df$Description,levels = rev(go.df$Description))

G<-ggplot(data = go.df, # 绘图使用的数据
        aes(x = Description, y = Count,fill = ONTOLOGY))+ # 横轴坐标及颜色分类填充
        geom_bar(stat = "identity",width = 0.9)+ # 绘制条形图及宽度设置
        coord_flip()+theme_bw()+ # 横纵坐标反转及去除背景色
        labs(x = "GO terms",y = "GeneNumber",title = "Barplot of Enriched GO Terms")+ # 设置坐标轴标题及标题
        theme(axis.title = element_text(size = 13), # 坐标轴标题大小
        axis.text = element_text(size = 11), # 坐标轴标签大小
        plot.title = element_text(size = 14,hjust = 0.5,face = "bold"), # 标题设置
        legend.title = element_text(size = 13), # 图例标题大小
        legend.text = element_text(size = 11), # 图例标签大小
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")) # 图边距
G

# 可视化
g1<-barplot(go)
g1
g2<-dotplot(go)
g2
## 还可以绘制GO的网络关系图，但是值得注意的是这里的数据只能是富集一个GO通路（BP、CC或MF）的数据
go.BP <- enrichGO(genelist, 
                  OrgDb = org.Hs.eg.db, 
                  ont='BP',pAdjustMethod = 'BH',
                  pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                  keyType = 'ENTREZID')
plotGOgraph(go.BP)

# Plots
grid.arrange(G, arrangeGrob(g1,g2 ,ncol=2), nrow = 2)

# kegg
kegg <- enrichKEGG(genelist, organism = "hsa",keyType = 'kegg', 
                   pvalueCutoff = 0.05,qvalueCutoff = 0.2,
                   use_internal_data = F) #借助在线数据库进行分析
#install.packages("R.utils")
library(R.utils)
R.utils::setOption("clusterProfiler.download.method","auto")

dim(kegg)
# 可视化
barplot(kegg)
dotplot(kegg)

