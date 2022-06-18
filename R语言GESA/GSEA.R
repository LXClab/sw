rm(list = ls())  
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(tidyverse)
library(ggstatsplot)
setwd("C:/Users/34750/Desktop/R教程（四）GSEA/")
data<-read.csv("normal-GH.csv")
head(data)
eg <- bitr(data$X, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
head(eg)
genelist<-data$logFC
names(genelist)<-eg$ENTREZID
geneList<-sort(genelist,decreasing = T)

library(R.utils)
R.utils::setOption("clusterProfiler.download.method","auto")
GO_entrez <- gseGO(geneList     = geneList,
                      ont          = "ALL",  # "BP"、"MF"和"CC"或"ALL"
                      OrgDb        = "org.Hs.eg.db",
                      keyType      = "ENTREZID",
                      pvalueCutoff = 0.35)   #实际为padj阈值可调整
GO <- DOSE::setReadable(GO_entrez, 
                        OrgDb="org.Hs.eg.db",
                        keyType='ENTREZID')#转化id 

KEGG_entrez <- gseKEGG(geneList     = geneList,
                        organism     = "hsa", 
                        pvalueCutoff = 0.35)  #实际为padj阈值可调整
KEGG <- DOSE::setReadable(KEGG_entrez, 
                             OrgDb="org.Hs.eg.db",
                             keyType='ENTREZID')#转化id

##选取富集结果
kk_gse <- KEGG
kk_gse_entrez <- KEGG_entrez

###条件筛选 
#一般认为|NES|>1，NOM pvalue<0.05，FDR（padj）<0.25的通路是显著富集的
kk_gse_cut <- kk_gse[kk_gse$pvalue<0.05 & kk_gse$p.adjust<0.25 & abs(kk_gse$NES)>1]
kk_gse_cut_down <- kk_gse_cut[kk_gse_cut$NES < 0,]
kk_gse_cut_up <- kk_gse_cut[kk_gse_cut$NES > 0,]

#选择展现NES前几个通路 
down_gsea <- kk_gse_cut_down[tail(order(kk_gse_cut_down$NES,decreasing = T),10),]
up_gsea <- kk_gse_cut_up[head(order(kk_gse_cut_up$NES,decreasing = T),10),]
diff_gsea <- kk_gse_cut[head(order(abs(kk_gse_cut$NES),decreasing = T),10),]


#### 经典的GSEA图 
up_gsea$Description
i=2
gseap1 <- gseaplot2(kk_gse,
                    up_gsea$ID[i],#富集的ID编号
                    title = up_gsea$Description[i],#标题
                    color = "red", #GSEA线条颜色
                    base_size = 20,#基础字体大小
                    rel_heights = c(1.5, 0.5, 1),#副图的相对高度
                    subplots = 1:3,   #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图
                    ES_geom = "line", #enrichment score用线还是用点"dot"
                    pvalue_table = T) #显示pvalue等信息
gseap1
#### 合并 GSEA通路 
gseap2 <- gseaplot2(kk_gse,
                    up_gsea$ID,#富集的ID编号
                    title = "UP_GSEA_all",#标题
                    color = "red",#GSEA线条颜色
                    base_size = 20,#基础字体大小
                    rel_heights = c(1.5, 0.5, 1),#副图的相对高度
                    subplots = 1:3, #要显示哪些副图 如subplots=c(1,3) #只要第一和第三个图
                    ES_geom = "line",#enrichment score用线还是用点"dot"
                    pvalue_table = T) #显示pvalue等信息
gseap2

#山脊图
ridgep <- ridgeplot(kk_gse_entrez,
                    showCategory = 15,
                    fill = "p.adjust",
                    core_enrichment = TRUE,
                    label_format = 30, #设置轴标签文字的每行字符数长度，过长则会自动换行。
                    orderBy = "NES",
                    decreasing = F) 
ridgep

#复现
gseaplot2(kk_gse,# gseaResult object
          up_gsea$ID[1], # gene set ID
          title = up_gsea$Description[1],# plot title
          color = "green",# color of running enrichment score line
          base_size = 11, # base font size
          rel_heights = c(1.5, 0.5, 1),# relative heights of subplots
          subplots = 1:3,# which subplots to be displayed
          pvalue_table = FALSE, # whether add pvalue table
          ES_geom = "line") # geom for plotting running enrichment score,
#AI添加网格线，注释信息

