#=======================
#	1.散点图
#=======================
#install.packages("ggplot2")
library(ggplot2)
head(mtcars)
#mpg：Miles Per Gallon 油耗(每加仑英里(美国))：功能更强大, 更重的汽车往往消耗更多的燃油。
ggplot(data = mtcars , mapping = aes(x = wt , y = mpg))
ggplot(data = mtcars , mapping = aes(x = wt , y = mpg)) + geom_point()
ggplot(data = mtcars , mapping = aes(x = wt , y = mpg , size = mpg)) + geom_point()
ggplot(data = mtcars , mapping = aes(x = wt , y = mpg , size = mpg , shape = as.factor(cyl) )) + geom_point()
ggplot(data = mtcars , mapping = aes(x = wt , y = mpg , size = mpg , color = mpg))+geom_point()
ggplot(mtcars , aes(x = wt , y = mpg , size = mpg , color = mpg))+geom_point()
ggplot(mtcars , aes(x = wt , y = mpg , size = mpg , color = mpg))+geom_point()+
  annotate("text",x=4,y=30,label="point_test",size=15)+
  labs(title = "mtcars_mpg vs mtcars_wt",x = "mtcars_wt", y = "mtcars_mpg")

#=======================
#	2.直方图（长度分布图）
#=======================
library(ggplot2)
setwd("C:\\Users\\34750\\Desktop\\")
data<-as.data.frame(read.csv("ecc_length_example.csv"))
ggplot(data,aes(x=Length,fill="red"))+geom_histogram(bins = 10)
ggplot(data,aes(x=Length,fill="red"))+geom_histogram(bins = 10)+
  labs(title = "Length of eccDNA" ,x="Length(bp)",y="Number")
ggplot(data,aes(x=Length,fill="red"))+geom_histogram(bins = 10)+
  labs(title = "Length of eccDNA" ,x="Length(bp)",y="Number")+
  theme_bw()
ggplot(data,aes(x=Length,fill="red"))+geom_histogram(bins = 10)+
  labs(title = "Length of eccDNA" ,x="Length(bp)",y="Number")+
  theme_bw() + 
  theme(panel.grid=element_blank())
ggplot(data,aes(x=Length,fill="red"))+geom_histogram(bins = 10)+
  labs(title = "Length of eccDNA" ,x="Length(bp)",y="Number")+
  theme_classic() + 
  theme(panel.grid=element_blank())+
  guides(fill=F)
#theme_bw()使用黑白主题
#theme_classic使用经典主题
#theme(panel.grid=element_blank())去掉网格线

#=======================
#	3.箱线图
#=======================
head(ToothGrowth)
#ToothGrowth为R内置数据集。它包含一项评估维生素C对豚鼠牙齿生长的影响的研究数据。实验在60只豚鼠上进行，其中每只豚鼠通过两种递送方法（橙汁，OJ，或抗坏血酸，VC）分别接受三种剂量水平的维生素C量（0.5、1和2 mg /天， VC）。
ToothGrowth$dose <- as.factor(ToothGrowth$dose)
ggplot(data = ToothGrowth,aes(x=supp,y=len))+geom_boxplot()
ggplot(data = ToothGrowth,aes(x=supp,y=len,fill=supp))+geom_boxplot()
ggplot(data = ToothGrowth,aes(x=supp,y=len,color=supp))+geom_boxplot()
#分组箱线图
ggplot(data = ToothGrowth,aes(x=supp,y=len, group=supp:dose))+geom_boxplot()
ggplot(data = ToothGrowth,aes(x=supp,y=len, group=supp:dose,fill=supp:dose))+geom_boxplot()

#=======================
#	4.小提琴图
#=======================
# Libraries
#install.packages("tidyverse")
#install.packages("hrbrthemes")
#install.packages("viridis")
library(tidyverse)
library(hrbrthemes)
library(viridis)
# create a dataset
data <- data.frame(
  name=c( rep("A",500), rep("B",500), rep("B",500), rep("C",20), rep('D', 100)  ),
  value=c( rnorm(500, 10, 5), rnorm(500, 13, 1), rnorm(500, 18, 1), rnorm(20, 25, 4), rnorm(100, 12, 1) )
)
# Plot
ggplot(data, aes(x=name, y=value, fill=name)) +
  geom_violin()
ggplot(data, aes(x=name, y=value, fill=name)) +
  geom_violin() +
  theme_ipsum() +
  ggtitle("Violin chart") +
  xlab("")
ggplot(data, aes(x=name, y=value, fill=name)) +
  geom_violin() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  theme_ipsum() +
  ggtitle("Violin chart") +
  xlab("") +
  ylab("") +
  guides(fill=F)
#scale_fill_viridis()：更改区域的填充颜色

#=======================
#	5.motif_logo图
#=======================
#install.packages("ggseqlogo")
library(ggseqlogo)
library(ggplot2)

data<-read.table("MA0268.1.txt",header=F,sep="\t")
seq<-data[,9]
ggseqlogo(seq)

#=======================
#	6.科研文献配色ggsci
#=======================
#install.packages("ggsci")
library(ggsci)
mtcars$cyl <- as.factor(mtcars$cyl)
p <- ggplot(mtcars, aes(x=cyl,y=mpg,fill=cyl))+geom_boxplot()
p+scale_fill_manual(values = c("red","green","blue"))
p+scale_fill_lancet() #柳叶刀杂志
p+scale_fill_aaas() #science
p+scale_fill_npg() #nature

help(package="ggsci")

#=======================
#	7.可交互式绘图plotly
#=======================
#install.packages("plotly")
library(plotly)
mtcars$cyl <- as.factor(mtcars$cyl)
p <- ggplot(mtcars, aes(x=cyl,y=mpg,fill=cyl))+geom_boxplot()
ggplotly(p)


#================================
#--------------------------------
#     ggplot2细节调整      
#--------------------------------
#================================

# install CRAN packages
#install.packages(c("tidyverse", "colorspace", "corrr",  "cowplot",
#                   "ggdark", "ggforce", "ggrepel", "ggridges", "ggsci",
#                  "ggtext", "ggthemes", "grid", "gridExtra", "patchwork",
#                   "rcartocolor", "scico", "showtext", "shiny",
#                   "plotly", "highcharter", "echarts4r"))
#
# install from GitHub since not on CRAN
#install.packages("devtools")
#devtools::install_github("JohnCoene/charter")

library(ggplot2)
head(mtcars)
#mpg：Miles Per Gallon 油耗(每加仑英里(美国))：功能更强大, 更重的汽车往往消耗更多的燃油。
g<-ggplot(data = mtcars , mapping = aes(x = wt , y = mpg)) 
g
g + geom_point() #散点图
g + geom_point() +geom_line() #点线图 
g + geom_point(color = "firebrick", shape = "diamond", size = 5) #设定颜色火红色，形状砖石状
theme_set(theme_bw()) #改变主题

#=======================
# 1.调整坐标轴
#=======================
#散点图
g + geom_point(color = "firebrick", shape = "diamond", size = 5) + labs(x = "Weight", y = "Miles Per Gallon")
g + geom_point(color = "firebrick", shape = "diamond", size = 5) + xlab("Weight") +  ylab("Miles Per Gallon")
p<-g + geom_point(color = "firebrick", shape = "diamond", size = 5) + xlab("Weight") +  ylab("Miles Per Gallon")
p
#Theme()是修改特定主题元素(文本、标题、框、符号、背景等)的必要命令。element_text()来改变所有或特定文本元素的属性(这里是axis标题)
#vjust指的是垂直对齐，它的范围通常在0和1之间
p+theme(axis.title.x = element_text(vjust = 1, size = 15),axis.title.y = element_text(vjust = 1, size = 15))
#在element_text()中,更改大小、颜色和外观
p+theme(axis.title = element_text(size = 15, color = "firebrick", face = "italic"))
#更改坐标轴文本
p+theme(axis.text = element_text(color = "dodgerblue", size = 12),axis.text.x = element_text(face = "italic"))
#旋转坐标轴文本
p+theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1, size = 12))
#删除坐标轴标题
p+labs(x = "", y = "")
#限制坐标轴范围(不在范围内的数据点删除，会报warning missing values)
p+ylim(c(0, 30))
p+scale_y_continuous(limits = c(0,30))
#限制坐标轴范围(保留不在范围内的数据点)
p+coord_cartesian(ylim = c(0,30))
#使图表从坐标原点开始
p+expand_limits(x = 0, y = 0)

#=======================
# 2.调整标题
#=======================
#添加标题
p+ggtitle("mtcars_mpg vs mtcars_wt")
p+labs(title = "mtcars_mpg vs mtcars_wt", x = "mtcars_weight", y = "mtcars_mpg")
#加粗标题
p+ggtitle("mtcars_mpg vs mtcars_wt")+theme(plot.title = element_text(face = "bold"))
#调整标题的位置、字体、大小
p+ggtitle("mtcars_mpg vs mtcars_wt")+theme(plot.title = element_text(face = "bold"))+theme(plot.title = element_text(hjust = 1, size = 16, face = "bold.italic"))

#=======================
# 3.调整图例
#=======================
#直方图
library(ggplot2)
setwd("C:\\Users\\34750\\Desktop\\")
data<-as.data.frame(read.csv("ecc_length_example.csv"))
theme_set(theme_classic()) #改变主题
l<-ggplot(data,aes(x=Length,fill="red"))+geom_histogram(bins = 10)
l
#去掉图例
l+ theme(legend.position = "none")
#去掉图例的标题
l+ theme(legend.title = element_blank())
#改变图例的位置，可以是top,bottom,left,right
l+ theme(legend.position = "top")

#=======================
# 4.调整背景和网格线
#=======================
#散点图
p
#更改背景颜色
p+theme(panel.background = element_rect(fill = "#64D2AA", color = "#64D2AA"))
#添加网格线
p+theme(panel.grid.major = element_line(color = "black", size = 0.5))
#dashed虚线、dotted斑点
p+theme(panel.grid.major = element_line(size = 0.5, linetype = "dashed"),
        panel.grid.minor = element_line(size = 0.25, linetype = "dotted"),
        panel.grid.major.x = element_line(color = "red1"),
        panel.grid.major.y = element_line(color = "blue1"),
        panel.grid.minor.x = element_line(color = "red4"),
        panel.grid.minor.y = element_line(color = "blue4"))
#去除网格线
p+theme(panel.grid = element_blank())

#=======================
# 5.调整边距、主题
#=======================
#调整边距margin
p+theme(plot.background = element_rect(fill = "gray60"),
      plot.margin = margin(t = 1, r = 3, b = 1, l = 2, unit = "cm"))
#改变主题,常用的有8种(详见PPT)
theme_set(theme_bw()) 
#改变所有文本的字体、大小
p + theme_bw(base_family = "", base_size = 20)

#=======================
# 6.调整线
#=======================
#添加垂直线
p+geom_vline(xintercept = c(2.5,4))
#添加水平线
p+geom_hline(yintercept = c(16, 28))
#添加虚线
p+geom_vline(aes(xintercept = median(wt)), size = 1, linetype = "dashed") +
  geom_hline(aes(yintercept = median(mpg)), size = 1,linetype = "dashed")

#=============================
# 7.调整文本、坐标，添加平滑线
#=============================
#添加文本注释
p+geom_text(aes(x = 4, y = 30,label = "This is a useful annotation"),size = 7, color = "darkcyan")

#翻转图表
p+coord_flip()
#反转轴（180度颠倒）
p+scale_y_reverse()
p+scale_x_reverse()

#添加平滑线
p+ stat_smooth()
#添加线性拟合
p+stat_smooth(method = "lm", se = FALSE, size = 1.3)


#=============================
# 8.facet分面_多面板图
#=============================
library(ggplot2)
head(mtcars)
#绘制单独的面板 – facet_null ()
ggplot(mtcars, aes(mpg, wt)) + geom_point()+facet_null()
#facet_wrap 单一变量：var ~ .
ggplot(mtcars, aes(mpg, wt)) + geom_point()+facet_wrap(facets = carb ~ .)
#facet_wrap 多变量：var1~var2
ggplot(mtcars, aes(mpg, wt)) + geom_point()+facet_wrap(facets = carb ~ cyl)

#facet_grid 一行多列或者是多行一列：. ~ var 或 var ~ .
ggplot(mtcars, aes(mpg, wt)) + geom_point()+facet_grid(facets = carb ~ .)
ggplot(mtcars, aes(mpg, wt)) + geom_point()+facet_grid(facets = . ~ carb)
#facet_grid 多行多列：var1~var2
ggplot(mtcars, aes(mpg, wt)) + geom_point()+facet_grid(facets = carb ~ cyl)

#facet_wrap 是一维，facet_grid 是二维。

###################
#构建数据集
set.seed(20220219)
#set.seed()用于设定随机数种子，一个特定的种子可以产生一个特定的伪随机序列。
#这个函数的主要目的，是让你的模拟能够可重复出现。
#因为很多时候我们需要取随机数，但这段代码再跑一次的时候，结果就不一样了。
#如果需要重复出现同样的模拟结果的话，就可以用set.seed()。
create_df <- rbind(data.frame(group="a", x = runif(100), y = rnorm(100, mean = 5)),
                   data.frame(group="b", x = runif(100), y = rnorm(100, mean = 5, sd = 3)+20),
                   data.frame(group="c", x = runif(100), y = rnorm(100, mean = 5, sd = 5)+30))
#runif()生成随机偏差
str(create_df)
head(create_df)
tail(create_df)

p1<-ggplot(data = create_df, aes(x = x, y = y, colour = group)) + geom_point(size = 3)
p1
p1+facet_wrap( ~ group) 
# 使用scale ="free_y"允许y轴在保持x轴不变的情况下变化
p1+facet_wrap( ~ group, scales = "free_y")
p1+facet_wrap( ~ group, scales = "free_y") + theme_bw()

###############
# libraries
library(ggplot2)
library(gridExtra)

# Make 3 simple graphics:
g1 <- ggplot(mtcars, aes(x=qsec)) + geom_density(fill="slateblue") #密度图
g1
g2 <- ggplot(mtcars, aes(x=drat, y=qsec, color=cyl)) + geom_point(size=5) + theme(legend.position="none")#散点图 去掉图例
g2
g3 <- ggplot(mtcars, aes(x=factor(cyl), y=qsec, fill=cyl)) + geom_boxplot() + theme(legend.position="none")#箱线图
g3
g4 <- ggplot(mtcars , aes(x=factor(cyl), fill=factor(cyl))) +  geom_bar()#柱状图
g4

# Plots
grid.arrange(g2, arrangeGrob(g3, g4, ncol=2), nrow = 2)
grid.arrange(g1, g2, g3, nrow = 3)
grid.arrange(g2, arrangeGrob(g3, g4, nrow=2), nrow = 1)


#=======================================
# 9.ggplot2实例1：nature箱线图+小提琴图
#=======================================
library(ggplot2)
setwd("C:/Users/34750/Desktop/")
data <- read.csv("violin-plot.CSV",header = T)
head(data)

#查看数据可知是一种宽数据，我们先转换成长数据格式，即把多个列名转换成一个变量列中的多个因子。

library(tidyverse) #包含了ggplot2包和数据处理用到的dplyr包
#像这种数据直接用gather()函数合并多列，然后修改合并后数据列名：
data <- gather(data)
colnames(data)<-c("Ancestry","Allele frequency")
head(data)

ggplot(data, aes(Ancestry,`Allele frequency`))+ #注意这种文本需要``保护，不是''
  geom_violin()+              #小提琴图
  geom_boxplot()              #箱线图

ggplot(data,aes(Ancestry,`Allele frequency`))+ 
  geom_violin(cex=1.2)+                 #cex设置边框粗细 width设置图的宽度
  geom_boxplot(width=0.1,cex=1.2)+      #设置箱线图宽度     
  theme_classic(base_size = 20)+        #可设置不同的主题 base_size设置主题基础字体大小
  theme(axis.text = element_text(color = 'black'), #默认坐标轴刻度颜色其实是灰色，这里我们修改成黑色
        legend.position = 'none')       #去掉图例

#使用ggplot2中的fill自动填充颜色
ggplot(data,aes(Ancestry,`Allele frequency`))+
  geom_violin(aes(fill=Ancestry),cex=1.2)+  
  geom_boxplot(width=0.1,cex=1.2)+
  theme_classic(base_size = 20)+
  theme(axis.text = element_text(color = 'black'),
        legend.position = 'none')

#颜色微调，复现文章中的图
ggplot(data,aes(Ancestry,`Allele frequency`))+
  geom_violin(aes(fill=Ancestry),cex=1.2)+  
  scale_fill_manual(values = c('#FB5554','#868B31','#42F203','#579ABB','#B978AE'))+
  geom_boxplot(width=0.1,cex=1.2)+
  theme_classic(base_size = 20)+
  theme(axis.text = element_text(color = 'black'),
        legend.position = 'none')

#===============================================
# 10.ggplot2实例2：nc箱线图、小提琴图、点图混合
#===============================================
library(ggplot2)
library(tidyverse)
library(gghalves)

df <- data.frame(subtills=c(rnorm(100,12,1.5),rnorm(20,16,5.5),rnorm(30,10,6)),
                 cereus=c(rnorm(60,13,2),rnorm(40,10,1.5),rnorm(50,6,1.5)),
                 megaterium=c(rnorm(130,6,1.5),rnorm(10,10,4.5),rnorm(10,5,1)),
                 circulans=c(rnorm(100,5,1.5),rnorm(40,3,1),rnorm(10,1,0.5))) %>% #tidyverse包的管道符
  gather() %>% #合并列，转换成长数据
  rename("Species"=key)  #修改列名

head(df)

ggplot(df)+
  geom_violin(aes(x=Species,y=value,fill=Species),cex=0.8)+
  geom_boxplot(aes(Species,value),width=0.2,cex=0.8)+
  scale_fill_manual(values = c('#84a354','#b15053','#cf9a2c','#7d7c7f'))+
  theme_test(base_size = 15)+
  labs(x=NULL,y='NO.of BGCs/genome')+
  scale_y_continuous(limits = c(0,25))+
  theme(legend.title = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(color = 'black'),
        axis.text.x = element_text(face = 'italic'))

#当然，一般情况画图到这个层面基本也就差不多了，但我们还是尽量还原文章中的图，以下几点我们需要解决：
#1.在箱线图、小提琴图旁边画出抖点图，像是雨滴；2.小提琴图只显示一半，像是一片云朵；3.几个变量顺序的调整。

#解决方案：

#1.如果直接加上抖点图，则会重叠到箱线图上，要绘出并排的效果需要一个小操作，直接用aes(x=Species)映射x轴，
#Species是一个因子向量，我们需要将其转换成一个数字型向量：as.numeric(factor(Species)，
#然后就可以对其进行±一定范围，使得抖点和箱线图分开并排。

#2.绘出一半的小提琴图，需要用到一个R包：gghalves，它给我们提供了一个geom_half_violin()的函数，
#当然也可以绘制半边的箱线图。

#3.变量顺序的自定义在我之前的文章中示例过几次了，有多种方法，这里我们可以先创建一个排好序的向量my_sort，
#将aes(x=Species)中替换成factor(Species,levels = my_sort)就好了。

variable <- c('subtills','cereus','megaterium','circulans')
my_sort <-factor(variable,levels = variable) #排好序

ggplot(df)+
  geom_half_violin(aes(as.numeric(factor(Species,levels = my_sort))+0.1,
                       value,fill=factor(Species,levels = my_sort)),
                   side = 'r',cex=0.8)+  #右侧半边显示 +0.1为了出的图形分布情况更合理
  geom_boxplot(aes(as.numeric(factor(Species,levels = my_sort))+0.1,
                   value,fill=factor(Species,levels = my_sort)),
               outlier.colour="black",width=0.1,cex=0.8)+
  geom_jitter(aes(as.numeric(factor(Species,levels = my_sort))-0.2,
                  value,color=factor(Species,levels = my_sort)), 
              width = 0.1,size=1.5) +
  scale_fill_manual(values = c('#84a354','#b15053','#cf9a2c','#7d7c7f'))+
  scale_color_manual(values = c('#84a354','#b15053','#cf9a2c','#7d7c7f'),guide='none')+
  scale_x_continuous(breaks = unique(as.numeric(factor(df$Species,levels = my_sort))), 
                     labels = unique(factor(df$Species,levels = my_sort)))+
  theme_test(base_size = 15)+
  labs(x=NULL,y='NO.of BGCs/genome')+
  scale_y_continuous(limits = c(0,25))+
  theme(legend.title = element_blank(),
        axis.title = element_text(size = 18),
        axis.text = element_text(color = 'black'),
        axis.text.x = element_text(face = 'italic'))
