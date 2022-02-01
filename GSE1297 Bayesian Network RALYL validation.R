###-------------------------------------------------------
### GSE1297 Bayesian Network RALYL Validation 
### 此脚本用于验证RALYL在多个AD数据集中的表达情况
### 验证RALYL表达与疾病进程、cog decline的关系
### lastupdate:20200619
### copyright@Yixuan Zhang
###-------------------------------------------------------

rm(list = ls())
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(reshape2)

##--------------------------------------------------------
## RALYL 在多个数据集中的表达情况验证
load("G:\\AD reserve omics\\Rproject\\GSE1297_BayesianNetwork\\RALYL validation\\GSE1297 original expression matrix.Rdata")
load("G:\\AD reserve omics\\Rproject\\GSE1297_BayesianNetwork\\RALYL validation\\GSE12685 and GSE28146 original expression matrix.Rdata")

# GSE12685 RALYL各组表达的分拣，但没有mo和se数据，不采用
Rna12 <- con12685[con12685$symbol %in% "RALYL",]
rownames(Rna12) <- "RALYL"; Rna12 <- Rna12[,c(-1,-2)]
Rin12 <- in12685[in12685$symbol %in% "RALYL",]
rownames(Rin12) <- "RALYL"; Rin12 <- Rin12[,c(-1,-2)]

## GSE28146 RALYL各组表达
Rna28 <- con28146[con28146$symbol %in% "RALYL",]
rownames(Rna28) <- "RALYL"; Rna28 <- as.data.frame(t(log(Rna28[,c(-1,-2)]))) 
Rin28 <- in28146[in28146$symbol %in% "RALYL",]
rownames(Rin28) <- "RALYL"; Rin28 <- as.data.frame(t(log(Rin28[,c(-1,-2)])))   
Rmo28 <- mo28146[mo28146$symbol %in% "RALYL",]
rownames(Rmo28) <- "RALYL"; Rmo28 <- as.data.frame(t(log(Rmo28[,c(-1,-2)]))) 
Rse28 <- se28146[se28146$symbol %in% "RALYL",]
rownames(Rse28) <- "RALYL"; Rse28 <- as.data.frame(t(log(Rse28[,c(-1,-2)]))) 
abc1 <- as.data.frame(matrix(NA,1,1)); colnames(abc1) <- "RALYL"
Rin28 <- rbind(Rin28,abc1); Rse28 <- rbind(Rse28,abc1)
ims28 <- cbind(Rna28,Rin28,Rmo28,Rse28)
colnames(ims28) <- c("na","inAD","moAD","seAD")
ims28_melt <- melt(ims28, value.name = "counts")
# 箱线图
p1 <- ggplot(ims28_melt, aes(x=variable,y=counts)) +
  geom_boxplot(aes(fill=variable)) +
  geom_jitter(shape=21, position=position_jitter(0.2)) +
  labs(title="Expression Levels of RALYL in GSE28146 ",
       x="Subjects Type", y = "RALYL gene expression levels", 
       fill = "Cell Type") +
  theme_bw() + guides(fill=FALSE) +
  theme(panel.grid=element_blank(),plot.title = element_text(hjust = 0.5),
        axis.text.x=element_text(size = 13),axis.text.y=element_text(size = 13));p1
# Multi test of p-value
p_adj <- compare_means(counts~variable,ims28_melt,method = "anova",
                       p.adjust.method = "bonferroni");p_adj
write.csv(p_adj,file = "G:\\AD reserve omics\\Rproject\\GSE1297_BayesianNetwork\\RALYL validation\\p.adjust of GSE28146 exp.csv")
# 显著性检验
my_cp <- list(c("na","inAD"),c("na","moAD"),c("na","seAD")) 
p1 + stat_compare_means(comparisons = my_cp) +
  stat_compare_means(method = "anova",label.y = 4)

## GSE1297 RALYL各组表达
R129 <- as.data.frame(t(exp1297[rownames(exp1297) %in% "RALYL",])) 
R129$group <- character(31)
for (i in 1:nrow(R129)) {
  for (j in 1:nrow(trait)) {
    if (rownames(R129)[i] %in% rownames(trait)[j]) {
      R129$group[i] <- trait$group[j]
        if (R129$group[i] %in% "0") {
          R129$group[i] <- "na"
        } else if (R129$group[i] %in% "1") {
          R129$group[i] <- "inAD"
        } else if (R129$group[i] %in% "2") {
          R129$group[i] <- "moAD"
        } else if (R129$group[i] %in% "3") {
          R129$group[i] <- "seAD"
        }
    }
  }
}
# 箱线图
R129$group = factor(R129$group, levels=c('na','inAD','moAD','seAD')) #柱条顺序
p1 <- ggplot(R129, aes(x=group,y=RALYL)) +
  geom_boxplot(aes(fill=group)) +
  geom_jitter(shape=21, position=position_jitter(0.2)) +
  labs(title="Expression Levels of RALYL in GSE1297 ",
       x="Subjects Type", y = "RALYL gene expression levels", 
       fill = "Cell Type") +
  theme_bw() + guides(fill=FALSE) +
  theme(panel.grid=element_blank(),plot.title = element_text(hjust = 0.5),
        axis.text.x=element_text(size = 13),axis.text.y=element_text(size = 13));p1
# Multi test of p-value
R129$group <- as.character(R129$group)# R129中的colname "group"在此处有歧义，替换
p_adj <- compare_means(counts~variable,R129,method = "anova",
                       p.adjust.method = "bonferroni");p_adj
write.csv(p_adj,file = "G:\\AD reserve omics\\Rproject\\GSE1297_BayesianNetwork\\RALYL validation\\p.adjust of GSE1297 exp.csv")
# 显著性检验
my_cp <- list(c("na","inAD"),c("na","moAD"),c("na","seAD")) 
p1 + stat_compare_means(comparisons = my_cp) +
  stat_compare_means(method = "anova",label.y = -2)

### 20200828注：以下代码出图无统计学差异，不用


##--------------------------------------------------------
## RALYL 在NCI/MCI/AD中的表达；在reserve/loss reserve的表达
# 载入M11表达矩阵数据，从中挑选RALYL表达数据
# load("G:\\AD reserve omics\\Rproject\\GSE1297_BayesianNetwork\\RALYL validation\\ROSMAP data sorting for RALYL validation.Rdata")
load("G:\\AD reserve omics\\Rproject\\GSE1297_BayesianNetwork\\External_validation\\ROSMAP_EVdata.Rdata")
# load("G:\\AD reserve omics\\Rproject\\GSE1297_BayesianNetwork\\External_validation\\ROSMAP M11 gene expression matrix.Rdata")

## RALYL高低表达的患者的cog decline情况
Xfl <- as.data.frame(t(expnfl[rownames(expnfl) %in% "RALYL",])); summary(Xfl)  
LOW = 10.875; HIGH = 11.149

for (i in 1:nrow(Xfl)) {
  if (Xfl$RALYL[i] <= LOW ) {
    Xfl$exp[i] <- "low"
  } else if (Xfl$RALYL[i] > HIGH) {
    Xfl$exp[i] <- "high"
  } else {
    Xfl$exp[i] <- NA
  }
}
# 分拣高低表达fl受试者的临床表型
lowfl <- follow[follow$ID %in% rownames(Xfl[Xfl$exp %in% "low",]),]
highfl <- follow[follow$ID %in% rownames(Xfl[Xfl$exp %in% "high",]),]
# 剔除MMSE_lv小于MMSE_dx的样本
lowfl <- lowfl[lowfl$mmse_dx >= lowfl$mmse_lv,]
highfl <- highfl[highfl$mmse_dx >= highfl$mmse_lv,]

# RALYL cog decline 具有显著性检验的普通箱线图
hmmse_dcl <- as.data.frame(highfl$mmse_dx - highfl$mmse_lv) 
colnames(hmmse_dcl) <- "high"; summary(hmmse_dcl)
lmmse_dcl <- as.data.frame(lowfl$mmse_dx - lowfl$mmse_lv)
colnames(lmmse_dcl) <- "low"; summary(lmmse_dcl)
abc <- data.frame(matrix(NA,1,1)); colnames(abc) <- "high" # 长度补齐
hmmse_dcl <- rbind(hmmse_dcl,abc)
mmse_dcl <- cbind(lmmse_dcl,hmmse_dcl)
mmse_dcl <- melt(mmse_dcl, value.name = "counts") # 出图部分程序
set.seed(1024); p1 <- ggplot(mmse_dcl,aes(x=variable,y=counts)) +
  geom_boxplot(aes(fill=variable)) +
  geom_jitter(shape=21, position=position_jitter(0.2)) +
  labs(title="MMSE decline during the follow-up",
       x="Expression levels of RALYL", y = "MMSE decline") +
  theme_bw() + guides(fill=FALSE) +
  theme(panel.grid=element_blank(),plot.title = element_text(hjust = 0.5),
        axis.text.x=element_text(size = 13),axis.text.y=element_text(size = 13));p1

my_cp <- list(c("high","low")) # 统计学检验部分
p1 + stat_compare_means(comparisons = my_cp) +
  stat_compare_means(method = "t.test",label.y = -2)

## 按照Braak评分分级后的RALYL表达情况
bk <- clinical[,c(1,10)]
Xall <- ROSMAPexpn[rownames(ROSMAPexpn) %in% "MAPK4" ,]
for (i in 1:nrow(bk)) {
  for (j in 1:ncol(Xall)) {
    if (bk$ID[i] %in% colnames(Xall)[j] ) {
      bk$RALYL[i] <- Xall[1,j]
    }
  }
}
NROW <- 1000000
for (i in 1:NROW) {
  if (is.na(bk$braak[i]) ) {
    bk <- bk[-i,]
  }
}
for (i in 1:NROW) {
  if (is.na(bk$RALYL[i])) {
    bk <- bk[-i,]
  }
}
bk$braak = factor(bk$braak, levels=c('1','2','3','4','5','6')) #柱条顺序
p1 <- ggplot(bk,aes(x=braak,y=RALYL)) + 
  geom_boxplot(aes(fill=braak)) +
  geom_jitter(shape=21, position=position_jitter(0.2)) +
  labs(title="MMSE decline during the follow-up",
       x="Braak stage", y = "RALYL expression levels") +
  theme_bw() + guides(fill=FALSE) +
  theme(panel.grid=element_blank(),plot.title = element_text(hjust = 0.5),
        axis.text.x=element_text(size = 13),axis.text.y=element_text(size = 13));p1
my_cp <- list(c("1","3")) # 统计学检验部分
p1 + stat_compare_means(comparisons = my_cp) +
  stat_compare_means(method = "t.test",label.y = -2)

## 各患者分组中RALYL的表达
Xna <- expnna[rownames(expnna) %in% "RALYL",]
Xmci <- expnmci[rownames(expnmci) %in% "RALYL",] 
Xad <- expnad[rownames(expnad) %in% "RALYL",]  
Xre <- expnre[rownames(expnre) %in% "RALYL",] 
Xlre <- expnlre[rownames(expnlre) %in% "RALYL",] 


abc1 <- matrix(numeric(46)[NA])
Xmci <- as.data.frame(rbind(t(Xmci),abc1))
abc2 <- matrix(numeric(2)[NA])
Xna <- as.data.frame(rbind(t(Xna),abc2)) 
Xad <- as.data.frame(t(Xad))

## AD/MCI/NA
library(reshape2) #用于配合ggplot2进行统计学检验
amn <- as.data.frame(cbind(Xad,Xmci,Xna))
colnames(amn) <- c("AD","MCI","NCI")
amn_melt <- melt(amn, value.name = "counts")

amn_melt$counts = log(amn_melt$counts)  # normalized

amn_melt$variable = factor(amn_melt$variable, levels=c('NCI','MCI','AD'))
p1 <- ggplot(amn_melt, aes(x=variable,y=counts)) +
  geom_boxplot(aes(fill=variable)) +
  geom_jitter(shape=21, position=position_jitter(0.2)) +
  labs(title="Expression Level of Module11 ",
       x="Subjects Type", y = "Module11 expression levals", 
       fill = "Cell Type") +
  theme_bw() + guides(fill=FALSE) +
  theme(panel.grid=element_blank(),plot.title = element_text(hjust = 0.5),
        axis.text.x=element_text(size = 13),axis.text.y=element_text(size = 13));p1
# 显著性检验：注意此方法存在bug，原图位置会被压缩，故只用数据不用图
my_cp <- list(c("AD","MCI"),c("AD","NCI"),c("MCI","NCI"))
p1 + stat_compare_means(comparisons = my_cp) +
  stat_compare_means(method = "anova",label.y =8)

## reserve/lossreserve
abc1 <- matrix(numeric(21)[NA])
Xre <- rbind(t(Xre),abc1)
Xlre <- as.data.frame(t(Xlre)) 

rl <- as.data.frame(cbind(Xre,Xlre))
colnames(rl) <- c("loss reserve","reserve")
rl_melt <- melt(rl, value.name = "counts")
rl_melt$counts = log(rl_melt$counts)

p2 <- ggplot(rl_melt, aes(x=variable,y=counts)) +
  geom_boxplot(aes(fill=variable)) +
  geom_jitter(shape=21, position=position_jitter(0.2)) +
  labs(title="Expression Level of RALYL ",
       x="Subjects Type", y = "RALYL gene expression levels", 
       fill = "Cell Type") +
  theme_bw() + guides(fill=FALSE) +
  theme(panel.grid=element_blank(),plot.title = element_text(hjust = 0.5),
        axis.text.x=element_text(size = 13),axis.text.y=element_text(size = 13));p2
# 显著性检验：注意此方法存在bug，原图位置会被压缩，故只用数据不用图
my_cp <- list(c("loss reserve","reserve"))
p2 + stat_compare_means(comparisons = my_cp) +
  stat_compare_means(method = "t.test",label.y =8)

##--------------------------------------------------------
## ROSMAP Follow-up组受试者按照RALYL表达高/低划分的认知下降趋势
Xfl <- as.numeric(t(expnfl[rownames(expnfl) %in% "RALYL",])) 
summary(Xfl)
LOW <- as.numeric(summary(Xfl)[2]); HIGH <- as.numeric(summary(Xfl)[5])
Xfl <- as.data.frame(t(expnfl[rownames(expnfl) %in% "RALYL",]))
Xfl$exp <- character(145)
# 判断每一个受试者主成分的高低
for (i in 1:nrow(Xfl)) {
  if (Xfl$RALYL[i] <= LOW) {
    Xfl$exp[i] <- "low"
  } else if (Xfl$RALYL[i] >= HIGH) {
    Xfl$exp[i] <- "high"
  } else {
    Xfl$exp[i] <- "NA"
  }
}
# 分拣高低表达fl受试者的临床表型
lowfl <- follow[follow$ID %in% rownames(Xfl[Xfl$exp %in% "low",]),]
highfl <- follow[follow$ID %in% rownames(Xfl[Xfl$exp %in% "high",]),]
# 剔除MMSE_lv小于MMSE_dx的样本
highfl <- highfl[highfl$mmse_dx >= highfl$mmse_lv,]
lowfl <- lowfl[lowfl$mmse_dx >= lowfl$mmse_lv,]
hmmse_dcl <- as.data.frame(highfl$mmse_dx - highfl$mmse_lv) 
colnames(hmmse_dcl) <- "high"; summary(hmmse_dcl)
lmmse_dcl <- as.data.frame(lowfl$mmse_dx - lowfl$mmse_lv)
colnames(lmmse_dcl) <- "low"; summary(lmmse_dcl)
abc <- data.frame(matrix(NA,2,1)); colnames(abc) <- "low"
lmmse_dcl <- rbind(lmmse_dcl,abc)
 
## 出图
# 箱线图
mmse_dcl <- cbind(hmmse_dcl,lmmse_dcl)
library(reshape2)
mmse_dcl <- melt(mmse_dcl, value.name = "counts")
p1 <- ggplot(mmse_dcl,aes(x=variable,y=counts)) +
  geom_boxplot(aes(fill=variable)) +
  geom_jitter(shape=21, position=position_jitter(0.2)) +
  labs(title="MMSE decline during the follow-up",
       x="Expression levels of RALYL", y = "MMSE decline ") +
  theme_bw() + guides(fill=FALSE) +
  theme(panel.grid=element_blank(),plot.title = element_text(hjust = 0.5),
        axis.text.x=element_text(size = 13),axis.text.y=element_text(size = 13));p1

my_cp <- list(c("high","low"))
p1 + stat_compare_means(comparisons = my_cp) +
  stat_compare_means(method = "t.test",label.y = 2.5)
 
## RALYL高表达的MMSE下降状况出图，注意受试者数量要手动修改
abc1 <- as.data.frame(highfl[,8]); abc2 <- as.data.frame(highfl[,9]) 
colnames(abc1) <- "MMSE"; colnames(abc2) <- "MMSE"
rownames(abc2) <- c(36:70); abc3 <- rbind(abc1,abc2)
abc3["X"] <- c(rep("Frist diagnostic",times=35),rep("Last vaild",times=35))
abc3["gp"] <- c(1:35,1:35); ghighfl <- abc3
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
rownames(abc2) <- c(32:62); abc3 <- rbind(abc1,abc2) #手动检查数量
abc3["X"] <- c(rep("Frist diagnostic",times=31),rep("Last vaild",times=31))
abc3["gp"] <- c(1:31,1:31); glowfl <- abc3
# 传统出图
ggplot(glowfl,aes(x=X,y=MMSE,group=gp)) + 
  geom_line(color=brewer.pal(9,"Set3")[4],size=0.75) + 
  xlab("Module11 Low expression") +
  geom_hline(aes(yintercept=26),colour=brewer.pal(9,"Set1")[9],linetype="dashed") +
  expand_limits(y=c(0,30)) +
  theme_bw() + theme(panel.grid=element_blank()) +
  theme(axis.text.x=element_text(size = 13),axis.text.y=element_text(size = 13))

