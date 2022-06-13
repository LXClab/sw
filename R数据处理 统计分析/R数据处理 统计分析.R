#=======================
#	1.R判断
#=======================
#if判断
x=120
if(x>100){print("x is good")}

#if..else判断
x <- c("google","runoob","taobao")
class(x)
if("runoob" %in% x)
{
	print("包含 runoob")
}else
  {
	print("不包含 runoob")
	}
#=======================
#	2.R循环
#=======================
#repeat循环
v <- c("Google","Runoob")
cnt <- 2

repeat{
	print(v)
	cnt <- cnt+1
	if(cnt > 5){
		break
	}
}

#while循环
v <- c("Google","Runoob")
cnt <- 2
while(cnt<7){
	print(v)
	cnt = cnt +1
}

#for循环
v <- LETTERS[1:4]
for (i in v){
	print(i)
}

#=======================
#	3.R函数
#=======================
#内置函数
print(seq(32,44))
print(mean(25:82))
print(sum(41:68))

#自定义函数(传参)
new.function<-function(a)
{
	for(i in 1:a)
	{
		b<-i^2
		print(b)
	}
}
new.function(6)

#自定义函数（不传参）
new.function<-function()
{
	for(i in 1:5)
	{
		print(i^2)
	}
}
new.function()
#=======================
#	4.R数据处理
#=======================
#paste(连接字符串)
a <- "Google"
b <- "Runoob"
d <- "Taobao"
print(paste(a,b,d))
print(paste(a,b,d, sep = ""))
print(paste(letters[1:6],1:6,sep = ""))
print(paste(letters[1:6],1:6,sep = "",collapse = "="))
paste(letters[1:6],1:6,collapse = ".")

#substring(提取字符串)
result<-substring("Runoob",2,5)
print(result)

#merge(合并数据框)
g<-c("gene1","gene2","gene3","human","mouse","human","high","high","1ow")
gene<-matrix(g,3,3)
gene<-as.data.frame (gene)
colnames(gene)<-c ("Gene" , "Species", "Expression") 
gene
p<-c("gene2","gene3","gene1","0.0723","0.0143","0.0053")
pvalue<-matrix(p,3,2)
pvalue<-data.frame(pvalue)
colnames(pvalue)<-c("Gene","P-value")
pvalue
result<-merge(gene,pvalue,by="Gene") 
result	

#=======================
#	1.独立性检验
#=======================
##独立性检验示例：通过Arthritis关节炎病人的数据资料，分析疾病改善水平与哪些因素相关(不相关，即为独立)
library(vcd) ##加载数据包
mytable<-table(Arthritis$Treatment,Arthritis$Improved) #Arthritis关节炎 table()函数计算频数 选择Treatment和Improved作为实验数据
mytable
chisq.test(mytable) #使用chisq.test()函数对二维表的行和列进行卡方独立性检验
# 由于p<0.01 ，患者接受的治疗和改善的水平一定程度上存在某种关系

mytable<-table(Arthritis$Sex,Arthritis$Improved) #Arthritis关节炎 table()函数计算频数 选择Sex和Improved作为实验数据
mytable
chisq.test(mytable) #使用chisq.test()函数对二维表的行和列进行卡方独立性检验
# 由于p>0.01 ，患者的性别和改善的水平不存在某种关系，说明患者的性别和疾病改善是相互独立的

#=======================
#	2、相关性分析
#=======================
library(ggcorrplot) ##ggcorrplot 基因相关性矩阵可视化的R包
setwd("C:/Users/34750/Desktop/")
exp=read.csv("expression.csv") #读取数据
rownames(exp)=exp[,1] #更改行名
exp=exp[,-1]  #去掉第一列
data=t(exp) #转置
co<-data[,1:10]
co
corr.p<-ggcorrplot::cor_pmat(co) #再计算P值
corr.p
test=cor(co)
test
ggcorrplot(test)
ggcorrplot(test,type = 'full',p.mat = corr.p,#P-Value
  	sig.level = 0.05 #P-Value大于0.05的在图中标记出来
)

#=======================
#	3.回归分析
#=======================
##简单线性回归
women   #为身高与体重的数据表
fit <- lm(weight ~ height, data=women)  
#拟合简单线性回归，通过身高（x）去预测体重（y）
summary(fit)   #可以得到很多信息

women$weight   #实际值
fitted(fit)     #模型对应的预测值
residuals(fit)   #模型误差
confint(fit)   #提供模型参数的置信区间
plot(women$height,women$weight,  #实际数据的点图
     main="Women Age 30-39", 
     xlab="Height (in inches)", 
     ylab="Weight (in pounds)")
abline(fit)         #模型拟合曲线

##多元线性回归
states <- as.data.frame(state.x77[,c("Murder", "Population", 
                                     "Illiteracy", "Income", "Frost")])
head(states)
fit2 <- lm(Murder ~ Population + Illiteracy + Income + Frost, data=states)
summary(fit2)

##logistics回归
library(DAAG)
fm<-glm(nomove~conc, family=binomial(link="logit"),data=anesthetic)
# nomove为上表中的y, conc为x conc表示麻醉剂的用量,move则表示手术病人是否有所移动，
#而我们用nomove做为因变量，因为研究的重点在于conc的增加是否会使nomove的概率增加。
summary(fm)

#=======================
#	4.方差分析
#=======================
data<- PlantGrowth ##示例数据
head(data)
summary(aov(weight~group, data=data)) ##aov方差分析
pvalue <- summary(aov(weight~group, data=data))[[1]][,"Pr(>F)"][1]
##Df为自由度，Sum Sq为总方差和，Mean Sq为平均方差和，F value为F统计量，Pr(>F)为P值
##Residuals 为残差部分的结果

#=======================
#	5.ROC曲线绘制
#=======================
library(ROCR)
data(ROCR.simple)
df <- data.frame(ROCR.simple)
head(df)
pred<-prediction(df$predictions, df$labels) ##使用ROCR的prediction函数进行模型预测
perf <-performance(pred,"tpr","fpr") ##计算模型的真阳率（True positive rate）和假阳率（False positive rate）
plot(perf , colorize=TRUE) ##绘制ROC曲线 

#=======================
#	6.主成分分析
#=======================
install.packages("psych") ##下载R包
library("psych") ##加载R包
Harman23.cor ##示例数据 对305名年龄在7至17岁之间的女孩进行的8项身体测量的相关矩阵。
pc <- principal(Harman23.cor$cov,nfactors = 2,rotate = 'none') 
##主成分分析,使用Harman23.cor数据集,nfactors代表主成分的个数,rotate代表是否进行主成分旋转
pc ##查看分析结果


