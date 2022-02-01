###-------------------------------------------------------
### GSE1297 Bayesian Network External Validation by ROSMAP
### 此脚本用ROSMAP真实世界数据对module_11外部验证
### 研究查看M11表达量高低与认知下降/临床分期/病理级别间的关系
### 代码行数较多，每段都可直接载入数据运行，不需要从头运行
### 20200828 针对审稿人意见进行p值多重检验
### lastupdate:20200828
### copyright@Yixuan Zhang
###-------------------------------------------------------

rm(list=ls())
library(ggplot2)
library(rayshader)
library(pheatmap)
library(RColorBrewer)
library(Rtsne)
library(ggpubr)

##--------------------------------------------------------
## 数据导入与清洗
load("GSE1297 25%geneinput all Module expression.Rdata") #GSE1297Module的基因名单

ROSMAPclic <- read.csv("G:\\AD reserve omics\\Rproject\\GSE1297_BayesianNetwork\\External_validation\\ROSMAP_clinical.csv",
                       header=TRUE, stringsAsFactors = FALSE)
ROSSMAPid <- read.csv("G:\\AD reserve omics\\Rproject\\GSE1297_BayesianNetwork\\External_validation\\ROSMAP_IDkey.csv",
                      header=TRUE, stringsAsFactors = FALSE)
ROSMAPexpn <- read.csv("G:\\AD reserve omics\\Rproject\\GSE1297_BayesianNetwork\\External_validation\\ROSMAP_arrayExpression_normalized.tsv",
                       header = TRUE,sep ="\t", stringsAsFactors = FALSE)
#提取关键临床信息
clinical <- data.frame(ID=character(1788),
                       projid=ROSMAPclic$projid,study=ROSMAPclic$study,
                       pmi=ROSMAPclic$pmi,
                       age=ROSMAPclic$age_at_visit_max,
                       sex=ROSMAPclic$msex,apoe=ROSMAPclic$apoe_genotype,
                       mmse_dx=ROSMAPclic$cts_mmse30_first_ad_dx,
                       mmse_lv=ROSMAPclic$cts_mmse30_lv,
                       braak=ROSMAPclic$braaksc,cerad=ROSMAPclic$ceradsc,
                       dx=ROSMAPclic$cogdx)

clinical$projid <- as.character(clinical$projid)
clinical$study <- as.character(clinical$study)
clinical$ID <- as.character(clinical$ID)

#projid补齐到8位数字
for (i in 1:nrow(clinical)) {
  abc <- nchar(clinical$projid[i])
  if (abc == 6) {
    clinical$projid[i]=paste("00",clinical$projid[i], sep = "")
  }
  if (abc == 7) {
    clinical$projid[i]=paste("0",clinical$projid[i], sep = "")
  }
}

#创建合并study和proj的正式ID
for (i in 1:nrow(clinical)) {
  clinical$ID[i] <- paste(clinical$study[i],clinical$projid[i],sep = "_")
}
#挑选有MMSE随访的受试者，命名为follow
follow <-clinical[!is.na(clinical$mmse_dx),]; head(follow)
#无MMSE_lv，Dx的受试者剔除
clinicalabc<-clinical[!is.na(clinical$mmse_lv) & !is.na(clinical$dx),]
clinical <- clinicalabc

#ROSMAP原始表达数据(normalized)的去重
ROSMAPexpn[1:5,1:5]
dat <- ROSMAPexpn[,-2:-3]#新建dat用于处理数据
rownames(dat) <- dat$ProbeID
ids <- ROSMAPexpn[,1:3];head(ids) #新建ids用于建立探针ID与基因名的联系

ids$median <- apply(dat,1,median)#求中位数
ids$median <- ids[order(ids$Symbol,ids$median,decreasing = TRUE),]#中位数降序
ids <- ids[!duplicated(ids$Symbol),]#IDS降序后去重
dat <- dat[rownames(dat) %in% ids$ProbeID,]#去重后的ProbeID再提取表达矩阵
rownames(dat) <- ids$Symbol#去重后的表达矩阵加上基因名称
ROSMAPexpn <- dat[,-1]#重建去重后的EXPN
ROSMAPexpn[1:5,1:5]

##--------------------------------------------------------
## 受试者分组(subre,sublre,subna,submci,subad)及表达矩阵提取

#AD reserve组受试者筛选，命名：subre
subre <- data.frame(ID=character(0),
                    projid=numeric(0),study=character(0),pmi=numeric(0),
                    age=character(0),sex=numeric(0),apoe=numeric(0),
                    mmse_dx=numeric(0),mmse_lv=numeric(0),
                    braak=numeric(0),cerad=numeric(0),dx=numeric(0))
for (i in 1:nrow(clinical)) {
  abc1 <- clinical$mmse_lv[i]  #MMSE>=27,即无认知障碍
  abc2 <- clinical$cerad[i] #改用CERAD评分，=为Definite AD
  abc4 <- clinical$dx[i] #不要Dx=1的受试者
  if (abc1 >= 27 & abc2<=2 & abc4!=1) {
    subre <- rbind(subre,clinical[i,])
  }
}

#AD loss reserve组受试者筛选，命名：sublre
#MMSE<10为重度AD
#此处根据Dx结论进行分组，Dx=4为AD dementia
sublre <- data.frame(ID=character(0),
                     projid=numeric(0),study=character(0),pmi=numeric(0),
                     age=character(0),sex=numeric(0),apoe=numeric(0),
                     mmse_dx=numeric(0),mmse_lv=numeric(0),
                     braak=numeric(0),cerad=numeric(0),dx=numeric(0))
for (i in 1:nrow(clinical)) {
  abc1 <- clinical$mmse_lv[i]
  abc2 <- clinical$dx[i]
  abc4 <- clinical$cerad[i]
  if (abc1 < 27 & abc2 == 4 & abc4 == 1) {
    sublre <- rbind(sublre,clinical[i,])
  }
}

#Normal aging组受试者筛选，命名：subna
subna <- data.frame(ID=character(0),
                    projid=numeric(0),study=character(0),pmi=numeric(0),
                    age=character(0),sex=numeric(0),apoe=numeric(0),
                    mmse_dx=numeric(0),mmse_lv=numeric(0),
                    braak=numeric(0),cerad=numeric(0),dx=numeric(0))
for (i in 1:nrow(clinical)) {
  abc1 <- clinical$mmse_lv[i]
  abc2 <- clinical$dx[i]
  abc3 <- clinical$cerad[i]
  if (abc1 >= 27 & abc2 ==1 & abc3 == 4) {
    subna <- rbind(subna,clinical[i,])
  }
}

#Dx MCI组受试者筛选，命名：submci
submci <- data.frame(ID=character(0),
                    projid=numeric(0),study=character(0),pmi=numeric(0),
                    age=character(0),sex=numeric(0),apoe=numeric(0),
                    mmse_dx=numeric(0),mmse_lv=numeric(0),
                    braak=numeric(0),cerad=numeric(0),dx=numeric(0))
for (i in 1:nrow(clinical)) {
  abc1 <- clinical$dx[i]
  abc2 <- clinical$mmse_lv[i]
  if (abc1 == 2 & abc2 <=26) {
    submci <- rbind(submci,clinical[i,])
  }
}

#Dx AD组受试者筛选，命名：subad
subad <- data.frame(ID=character(0),
                     projid=numeric(0),study=character(0),pmi=numeric(0),
                     age=character(0),sex=numeric(0),apoe=numeric(0),
                     mmse_dx=numeric(0),mmse_lv=numeric(0),
                     braak=numeric(0),cerad=numeric(0),dx=numeric(0))
for (i in 1:nrow(clinical)) {
  abc1 <- clinical$dx[i]
  abc2 <- clinical$mmse_lv[i]
  if (abc1 == 4 & abc2 <=26) {
    subad <- rbind(subad,clinical[i,])
  }
}

#AD reserve 受试者表达数据，命名expnre
expnre <- ROSMAPexpn[,colnames(ROSMAPexpn) %in% subre$ID]
expnre[1:5,1:4]
#AD loss reserve 受试者表达数据，命名expnlre
expnlre <- ROSMAPexpn[,colnames(ROSMAPexpn) %in% sublre$ID]
expnlre[1:5,1:4]
#Normal aging 受试者表达数据，命名expnna
expnna <- ROSMAPexpn[,colnames(ROSMAPexpn) %in% subna$ID]
expnna[1:5,1:4]
#Dx AD 受试者表达数据，命名expnad
expnad <- ROSMAPexpn[,colnames(ROSMAPexpn) %in% subad$ID]
expnad[1:5,1:4]
#Dx MCI 受试者表达数据，命名expnmci
expnmci <- ROSMAPexpn[,colnames(ROSMAPexpn) %in% submci$ID]
expnmci[1:5,1:4]
#follow 受试者表达数据，命名expnfl
expnfl <- ROSMAPexpn[,colnames(ROSMAPexpn) %in% follow$ID]
expnfl[1:5,1:4]

##--------------------------------------------------------
## 整理好的数据保存
save(ROSMAPexpn,ROSMAPclic,clinical,
     follow,expnfl,
     subre,sublre,subna,submci,subad,
     expnre,expnlre,expnna,expnmci,expnad,
     exp_yellow,exp_black,
     file = "G:\\AD reserve omics\\Rproject\\GSE1297_BayesianNetwork\\External_validation\\ROSMAP_EVdata.Rdata")
rm(list = ls())

##--------------------------------------------------------
## 各个分组Module11基因表达矩阵
load("G:\\AD reserve omics\\Rproject\\GSE1297_BayesianNetwork\\External_validation\\ROSMAP_EVdata.Rdata")
# yellow模块即为M11
# 函数z-score：表达矩阵的归一化
zscore <- function(x) {
  dat <- x
  pheatmap(dat,show_colnames =FALSE, show_rownames = FALSE)
  n=t(scale(t(dat)))
  n[n>2]=2 #限定上限，使表达量大于2的等于2
  n[n< -2]= -2 #限定下限，使表达量小于-2的等于-2
  n[1:4,1:4]
  dat <- as.data.frame(n)
  # pheatmap(dat,show_colnames =FALSE, show_rownames = FALSE)
  return(dat)
}

m11re <- expnre[rownames(expnre) %in% rownames(exp_yellow),];m11re <- zscore(m11re)
m11lre <- expnlre[rownames(expnlre) %in% rownames(exp_yellow),];m11lre <- zscore(m11lre)
m11na <- expnna[rownames(expnna) %in% rownames(exp_yellow),];m11na <- zscore(m11na)
m11mci <- expnmci[rownames(expnmci) %in% rownames(exp_yellow),];m11mci <- zscore(m11mci)
m11ad <- expnad[rownames(expnad) %in% rownames(exp_yellow),];m11ad <- zscore(m11ad)
m11fl <- expnfl[rownames(expnfl) %in% rownames(exp_yellow),];m11fl <- zscore(m11fl)
m11ROSMAP <- ROSMAPexpn[rownames(ROSMAPexpn) %in% rownames(exp_yellow),]
m11ROSMAP <- zscore(m11ROSMAP)

save(m11re,m11lre,m11na,m11mci,m11ad,m11fl,m11ROSMAP,
     file = "G:\\AD reserve omics\\Rproject\\GSE1297_BayesianNetwork\\External_validation\\ROSMAP M11 gene expression matrix.Rdata")
# 20200320发现：原本GSE1297的M11（yellow）模块中有262个基因；
# 但是上面在ROSMAP数据集中挑选M11基因的时候，发现只有225个基因入选
# GPL96与GPL17518芯片差异所致，不影响整体结果

##--------------------------------------------------------
## 模块表达高低的求解：非参数检验-核密度估计KDE
## 注：此方法由于性能不佳，不采用，详见20200317 log

# 计算模块总体分布
all_gene <- as.numeric(as.matrix(m11fl))
plot(density(all_gene, bw = "SJ"))
abc1 <- density(all_gene, bw = "SJ"); abc1
kde <- as.matrix(summary(abc1$x)) # 模块总体四分位数信息

# 计算具体患者的分布信息
expn <- m11fl
hl <- as.data.frame(colnames(expn));rownames(hl)<-hl[,1];
colnames(hl)<-"HighorLow";hl$HighorLow<-as.character(hl$HighorLow)
hl["median"] <- "NA"
for (i in 1:ncol(expn)) {
  abc1 <- expn[,i]
  abc2 <- density(abc1, bw = "SJ")
  abc3 <- as.matrix(summary(abc2$x))
  if ( kde[3,1] - abc3[3,1] > 0.0051) {
    hl[i,1] <- "Low_expression"; hl[i,2] <- abc3[3,1]
  } else if (abc3[3,1] - kde[3,1] > 0.005 ) {
    hl[i,1] <- "High_expression"; hl[i,2] <- abc3[3,1]
  } else {
    hl[i,1] <- "No_difference"; hl[i,2] <- abc3[3,1]
  }
}

# M11高低表达下随访认知功能的下降数据提取
rm(abc1,abc2)
abc1 <- data.frame("Highorlow"=character(0),"median"=character(0))
abc2 <- data.frame("Highorlow"=character(0),"median"=character(0))
for (i in 1:nrow(hl)) {
  if (hl$HighorLow[i] %in% "High_expression" ) {
    abc1 <- rbind(abc1,hl[i,]) 
  } else if (hl$HighorLow[i] %in% "Low_expression" ) {
    abc2 <- rbind(abc2,hl[i,])
  }
}
hm11fl <- follow[follow$ID %in% rownames(abc1),]
lm11fl <- follow[follow$ID %in% rownames(abc2),]

# M11高表达的MMSE下降状况出图
abc1 <- as.data.frame(hm11fl[,8]); abc2 <- as.data.frame(hm11fl[,9]) 
colnames(abc1) <- "MMSE"; colnames(abc2) <- "MMSE"
rownames(abc2) <- c(14:26); abc3 <- rbind(abc1,abc2)
abc3["X"] <- c(rep("Diagnostic MMSE",times=13),rep("Last vaild MMSE",times=13))
abc3["gp"] <- c(1:13,1:13); phm11fl <- abc3

display.brewer.all(type = "all")
ggplot(phm11fl,aes(x=X,y=MMSE,group=gp)) + 
  geom_line(color=brewer.pal(9,"YlOrRd")[6]) + 
  geom_point() +
  geom_hline(aes(yintercept=26),colour=brewer.pal(9,"Set1")[9],linetype="dashed") +
  expand_limits(y=c(0,30)) +
  theme_bw() + theme(panel.grid=element_blank())

##--------------------------------------------------------
## M11模块表达高低的求解：t-SNE 2D 降维考察主成分
# 对各组受试者的M11基因降维，取差异大的主成分进行作图
# 为了各组间降维尺度一致，每次降维必须运行种子设置

# 载入M11表达矩阵数据
load("G:\\AD reserve omics\\Rproject\\GSE1297_BayesianNetwork\\External_validation\\ROSMAP M11 gene expression matrix.Rdata")

set.seed(1024) # AD 
m11ad <- as.data.frame(t(m11ad))
tsad <- Rtsne(m11ad,dims=2,pca=TRUE,perplexity=10,theta=0.0)
plot(tsad$Y,asp=1,xlab = "t-SNE out [1]",ylab = "t-SNE out [2]",col="red",
     main = "AD Module11 t-SNE dimensionality reduction")
tsad2 <- tsad$Y[,2]; summary(tsad2)

set.seed(1024) # MCI
m11mci <- as.data.frame(t(m11mci))
tsmci <- Rtsne(m11mci,dims=2,pca=TRUE,perplexity=10,theta=0.0)
plot(tsmci$Y,asp=1,xlab = "t-SNE out [1]",ylab = "t-SNE out [2]",col="green3",
     main = "MCI Module11 t-SNE dimensionality reduction")
tsmci2 <- tsmci$Y[,2]; summary(tsmci2)

set.seed(1024) # NA
m11na <- as.data.frame(t(m11na))
tsna <- Rtsne(m11na,dims=2,pca=TRUE,perplexity=10,theta=0.0)
plot(tsna$Y,asp=1,xlab = "t-SNE out [1]",ylab = "t-SNE out [2]",col="blue",
     main = "NA Module11 t-SNE dimensionality reduction")
tsna2 <- tsna$Y[,2]; summary(tsna2)
  
set.seed(1024) # reserve
m11re <- as.data.frame(t(m11re))
tsre <- Rtsne(m11re,dims=2,pca=TRUE,perplexity=10,theta=0.0)
plot(tsre$Y,asp=1,xlab = "t-SNE out [1]",ylab = "t-SNE out [2]",col="turquoise3",
     main = "Rserve Module11 t-SNE dimensionality reduction")
tsre1 <- tsre$Y[,1]; summary(tsre1)

set.seed(1024) # loss reserve
m11lre <- as.data.frame(t(m11lre))
tslre <- Rtsne(m11lre,dims=2,pca=TRUE,perplexity=10,theta=0.0)
plot(tslre$Y,asp=1,xlab = "t-SNE out [1]",ylab = "t-SNE out [2]",col="pink3",
     main = "Loss rserve Module11 t-SNE dimensionality reduction")
tslre1 <- tslre$Y[,1]; summary(tslre1)

set.seed(1024) # Follow-up
tsfl <- Rtsne(m11fl,dims=2,pca=TRUE,perplexity=30,theta=0.0)
plot(tsfl$Y,asp=1,xlab = "t-SNE out [1]",ylab = "t-SNE out [2]",
     main = "Follow-up Module11 t-SNE dimensionality reduction")
tsfl1 <- tsfl$Y[,1]; summary(tsfl1) 
tsfl2 <- tsfl$Y[,2]; summary(tsfl2) 

save(tsfl1,tsfl2,tsad2,tsmci2,tsna2,tsre1,tslre1,
     file = "G:\\AD reserve omics\\Rproject\\GSE1297_BayesianNetwork\\External_validation\\ROSMAP M11 tSNE.Rdata")
##--------------------------------------------------------
## ROSMAP M11的t-SNE结果散点图可视化:2D 

# t-SNE降维图的2D作图(AD/MCI/NA三个的降维图合一)
dtsad <- as.data.frame(tsad$Y)
dtsmci <- as.data.frame(tsmci$Y)
dtsna <- as.data.frame(tsna$Y) 
dtsad$group <- character(151); dtsad$group <- "AD"
dtsmci$group <- character(57); dtsmci$group <- "MCI"
dtsna$group <- character(67); dtsna$group <- "NCI" # "NA"易混淆，改称NCI
dotts <- rbind(dtsad,dtsmci,dtsna)
# 普通2D散点图
ggplot(dotts,aes(x=V1,y=V2,color=group)) + 
  geom_point(size=2) + theme_bw() + 
  theme(panel.grid=element_blank(),plot.title = element_text(hjust = 0.5)) +
  labs(title="ROSMAP Module11 t-SNE dimensionality reduction",
       x="dim 1", y = "dim 2")
# 绘制密度等高线图
gg <- ggplot(dotts, aes(V1, V2)) +
        stat_density_2d(aes(fill = stat(nlevel)), 
                        geom = "polygon",
                        n = 100,bins = 10, contour = TRUE) +
        facet_wrap(group~.) +    # 按group分类
        scale_fill_viridis_c(option = "D"); gg

plot_gg(gg,multicore=TRUE,width=5,height=5,scale=250)

##--------------------------------------------------------
## ROSMAP M11的t-SNE结果散点图可视化:3D 
# 3D作图的数据需要重新降维，故重新载入数据计算
rm(list = ls())
load("G:\\AD reserve omics\\Rproject\\GSE1297_BayesianNetwork\\External_validation\\ROSMAP M11 gene expression matrix.Rdata")
set.seed(1024) # AD 
m11ad <- as.data.frame(t(m11ad))
tsad3 <- Rtsne(m11ad,dims=3,pca=TRUE,perplexity=10,theta=0.0)
dtsad3 <- as.data.frame(tsad3[["Y"]])
set.seed(1024) # MCI 
m11mci <- as.data.frame(t(m11mci))
tsmci3 <- Rtsne(m11mci,dims=3,pca=TRUE,perplexity=10,theta=0.0)
dtsmci3 <- as.data.frame(tsmci3[["Y"]])
set.seed(1024) # NA 
m11na <- as.data.frame(t(m11na))
tsna3 <- Rtsne(m11na,dims=3,pca=TRUE,perplexity=10,theta=0.0)
dtsna3 <- as.data.frame(tsna3[["Y"]])

dtsad3$group <- character(151); dtsad3$group <- "AD"
dtsmci3$group <- character(57); dtsmci3$group <- "MCI"
dtsna3$group <- character(67); dtsna3$group <- "NCI" # "NA"易混淆，改称NCI
dotts3 <- rbind(dtsad3,dtsmci3,dtsna3)

library(plot3D)
attach(dotts3)
mycolor <- c(rep("#F000007F",151),rep("#00FF007F",57),rep("#004EFF7F",67))

scatter3D(V1, V2, V3, col= mycolor, phi = 0, theta = 45,
          main = "ROSMAP Module11 t-SNE dimensionality reduction in 3D",
          xlab = "Dim 1", ylab = "Dim 2", zlab = "Dim 3",
          bty = "g", cex = 0.75, colkey = FALSE)  

save(dotts3,
     file = "G:\\AD reserve omics\\Rproject\\GSE1297_BayesianNetwork\\External_validation\\ROSMAP M11 tSNE 2D&3D data.Rdata")


##--------------------------------------------------------
## ROSMAP M11 AD/MCI/NA，reserve/lossreserve组M11模块表达的箱线图

rm(list = ls())
library(reshape2) #用于配合ggplot2进行统计学检验
load("G:\\AD reserve omics\\Rproject\\GSE1297_BayesianNetwork\\External_validation\\ROSMAP M11 tSNE.Rdata")

abc1 <- matrix(numeric(94)[NA])
tsmci2 <- rbind(as.matrix(tsmci2),abc1)
abc2 <- matrix(numeric(84)[NA])
tsna2 <- rbind(as.matrix(tsna2),abc2)
# AD/MCI/NA
amn <- as.data.frame(cbind(tsad2,tsmci2,tsna2))
colnames(amn) <- c("AD","MCI","NCI")
amn_melt <- melt(amn, value.name = "counts")
amn_melt$counts = log(amn_melt$counts) # normalized
# 出图
p1 <- ggplot(amn_melt, aes(x=variable,y=counts)) +
  geom_boxplot(aes(fill=variable)) +
  geom_jitter(shape=21, position=position_jitter(0.2)) +
  labs(title="Expression Level of Module11 ",
       x="Subjects Type", y = "Module11 expression levels", 
       fill = "Cell Type") +
  theme_bw() + guides(fill=FALSE) +
  theme(panel.grid=element_blank(),plot.title = element_text(hjust = 0.5),
        axis.text.x=element_text(size = 13),axis.text.y=element_text(size = 13));p1
# Multi test of p-value
p_adj <- compare_means(counts~variable,amn_melt,method = "anova",
                       p.adjust.method = "bonferroni");p_adj
write.csv(p_adj,file = "G:\\AD reserve omics\\Rproject\\GSE1297_BayesianNetwork\\External_validation\\p.adjust of M11_ADMCINCI exp.csv")
# 显著性检验：注意此方法存在bug，原图位置会被压缩，故只用数据不用图
my_cp <- list(c("AD","MCI"),c("AD","NCI"),c("MCI","NCI"))
p1 + stat_compare_means(comparisons = my_cp) +
  stat_compare_means(method = "anova",label.y = 6)

# reserve/lossreserve
abc1 <- matrix(numeric(33)[NA])
tsre1 <- rbind(as.matrix(tsre1),abc1)

rl <- as.data.frame(cbind(tsre1,tslre1))
colnames(rl) <- c("loss reserve","reserve")
rl_melt <- melt(rl, value.name = "counts")
rl_melt$counts = log(rl_melt$counts)

p2 <- ggplot(rl_melt, aes(x=variable,y=counts)) +
  geom_boxplot(aes(fill=variable)) +
  geom_jitter(shape=21, position=position_jitter(0.2)) +
  labs(title="Expression Level of Module11 ",
       x="Subjects Type", y = "Module11 expression levals", 
       fill = "Cell Type") +
  theme_bw() + guides(fill=FALSE) +
  theme(panel.grid=element_blank(),plot.title = element_text(hjust = 0.5),
        axis.text.x=element_text(size = 13),axis.text.y=element_text(size = 13));p2
# Multi test of p-value
p_adj <- compare_means(counts~variable,rl_melt,method = "t.test",
                       p.adjust.method = "bonferroni");p_adj
write.csv(p_adj,file = "G:\\AD reserve omics\\Rproject\\GSE1297_BayesianNetwork\\External_validation\\p.adjust of M11_reservelossreserve exp.csv")
# 显著性检验：注意此方法存在bug，原图位置会被压缩，故只用数据不用图
my_cp <- list(c("loss reserve","reserve"))
p2 + stat_compare_means(comparisons = my_cp) +
  stat_compare_means(method = "t.test",label.y = 2)



##--------------------------------------------------------
## ROSMAP M11 Follow-up组受试者按照M11表达高/低划分的认知下降趋势图

rm(list = ls())
library(Rtsne)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(reshape2)

load("G:\\AD reserve omics\\Rproject\\GSE1297_BayesianNetwork\\External_validation\\ROSMAP_EVdata.Rdata")
load("G:\\AD reserve omics\\Rproject\\GSE1297_BayesianNetwork\\External_validation\\ROSMAP M11 gene expression matrix.Rdata")
# 转置M11fl表达矩阵，将受试者作为因变量进行降维，获得145个主成分用来判断高低:
set.seed(1024)
abc1 <- as.data.frame(t(m11fl))
flts <- Rtsne(abc1,dims=2,pca=TRUE,perplexity=5,theta=0.0)
plot(flts$Y,asp=1,xlab = "t-SNE out [1]",ylab = "t-SNE out [2]",col="red",
     main = "Follow-up subjects t-SNE dimensionality reduction")
flts1 <- flts$Y[,2]; summary(flts1)
LOW <- as.numeric( summary(flts1)[2]); HIGH <- as.numeric(summary(flts1)[5])  
# 提取fl组各个受试者降维的结果
tsneout <-as.data.frame(flts$Y)  
rownames(tsneout) <- rownames(abc1)
tsneout$exp <- character(145)
# 判断每一个受试者主成分的高低
for (i in 1:nrow(tsneout)) {
  if (tsneout$V2[i] <= LOW) {
    tsneout$exp[i] <- "low"
  } else if (tsneout$V2[i] >= HIGH) {
    tsneout$exp[i] <- "high"
  } else {
    tsneout$exp[i] <- "NA"
  }
}
# 分拣高低表达fl受试者的临床表型
lowfl <- follow[follow$ID %in% rownames(tsneout[tsneout$exp %in% "low",]),]
highfl <- follow[follow$ID %in% rownames(tsneout[tsneout$exp %in% "high",]),]
# 剔除MMSE_lv小于MMSE_dx的样本
lowfl <- lowfl[lowfl$mmse_dx >= lowfl$mmse_lv,]
highfl <- highfl[highfl$mmse_dx >= highfl$mmse_lv,]

# M11 cog decline 具有显著性检验的普通箱线图
hmmse_dcl <- as.data.frame(highfl$mmse_dx - highfl$mmse_lv) 
colnames(hmmse_dcl) <- "high"; summary(hmmse_dcl)
lmmse_dcl <- as.data.frame(lowfl$mmse_dx - lowfl$mmse_lv)
colnames(lmmse_dcl) <- "low"; summary(lmmse_dcl)
abc <- data.frame(matrix(NA,3,1)); colnames(abc) <- "high" # 长度补齐
hmmse_dcl <- rbind(hmmse_dcl,abc)
mmse_dcl <- cbind(lmmse_dcl,hmmse_dcl)

mmse_dcl <- melt(mmse_dcl, value.name = "counts") # 出图部分程序
set.seed(1024); p1 <- ggplot(mmse_dcl,aes(x=variable,y=counts)) +
  geom_boxplot(aes(fill=variable)) +
  geom_jitter(shape=21, position=position_jitter(0.2)) +
  labs(title="MMSE decline during the follow-up",
       x="Expression levels of Module 11", y = "MMSE decline") +
  theme_bw() + guides(fill=FALSE) +
  theme(panel.grid=element_blank(),plot.title = element_text(hjust = 0.5),
        axis.text.x=element_text(size = 13),axis.text.y=element_text(size = 13));p1
# Multi test of p-value
p_adj <- compare_means(counts~variable,mmse_dcl,method = "t.test",
                       p.adjust.method = "bonferroni");p_adj
write.csv(p_adj,file = "G:\\AD reserve omics\\Rproject\\GSE1297_BayesianNetwork\\External_validation\\p.adjust of MMSEdecline exp.csv")
# 统计学检验部分
my_cp <- list(c("high","low")) 
p1 + stat_compare_means(comparisons = my_cp) +
  stat_compare_means(method = "t.test",label.y = -2)

# M11高表达的MMSE下降状况出图，注意受试者数量要手动修改,
# 注意需要载入柱状图代码前面的数据！
abc1 <- as.data.frame(highfl[,8]); abc2 <- as.data.frame(highfl[,9]) 
colnames(abc1) <- "MMSE"; colnames(abc2) <- "MMSE"
rownames(abc2) <- c(33:64); abc3 <- rbind(abc1,abc2)
abc3["X"] <- c(rep("Frist diagnostic",times=32),rep("Last vaild",times=32))
abc3["gp"] <- c(1:32,1:32); ghighfl <- abc3
#display.brewer.all(type = "all") #显示绘图色块
ggplot(ghighfl,aes(x=X,y=MMSE,group=gp)) + 
  geom_line(color=brewer.pal(9,"Set3")[5],size=0.75) + 
  xlab("Module11 High expression") + 
  geom_hline(aes(yintercept=26),colour=brewer.pal(9,"Set1")[9],linetype="dashed") +
  expand_limits(y=c(0,30)) +
  theme_bw() + theme(panel.grid=element_blank()) +
  theme(axis.text.x=element_text(size = 13),axis.text.y=element_text(size = 13))

# M11低表达的MMSE下降状况出图，注意受试者数量要手动修改
abc1 <- as.data.frame(lowfl[,8]); abc2 <- as.data.frame(lowfl[,9]) 
colnames(abc1) <- "MMSE"; colnames(abc2) <- "MMSE"
rownames(abc2) <- c(36:70); abc3 <- rbind(abc1,abc2) #手动检查数量
abc3["X"] <- c(rep("Frist diagnostic",times=35),rep("Last vaild",times=35))
abc3["gp"] <- c(1:35,1:35); glowfl <- abc3
# 传统出图
ggplot(glowfl,aes(x=X,y=MMSE,group=gp)) + 
  geom_line(color=brewer.pal(9,"Set3")[4],size=0.75) + 
  xlab("Module11 Low expression") +
  geom_hline(aes(yintercept=26),colour=brewer.pal(9,"Set1")[9],linetype="dashed") +
  expand_limits(y=c(0,30)) +
  theme_bw() + theme(panel.grid=element_blank()) +
  theme(axis.text.x=element_text(size = 13),axis.text.y=element_text(size = 13))

# save(tsneout,highfl,lowfl,ghighfl,glowfl, # 保存
#      file = "G:\\AD reserve omics\\Rproject\\GSE1297_BayesianNetwork\\External_validation\\ROSMAP M11 highlow exp MMSE decline.Rdata")
