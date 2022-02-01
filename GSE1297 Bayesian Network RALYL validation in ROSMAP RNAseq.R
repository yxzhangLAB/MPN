###-------------------------------------------------------
### GSE1297 Bayesian Network RALYL Validation - ROSMAP RNAseq
### 此脚本用ROSMAP RNAseq，并挑选临床表型
### 验证RALYL在多个AD数据集中的表达情况
### lastupdate:20200629
### copyright@Yixuan Zhang
###-------------------------------------------------------

rm(list = ls())
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(reshape2)

##--------------------------------------------------------
## 五组受试者微阵列数据的提取

load("G:\\AD reserve omics\\Rproject\\GSE1297_BayesianNetwork\\RALYL validation\\ROSMAP RNAseq expression matrix")
load("G:\\AD reserve omics\\Rproject\\GSE1297_BayesianNetwork\\External_validation\\ROSMAP_EVdata.Rdata")
#AD reserve 受试者表达数据，命名expnre
expnre <- ralyl[,colnames(ralyl) %in% subre$ID]
expnre <- cbind(ralyl[,c(1,2,3)],expnre)
rownames(expnre) <- expnre$SYMBOL; expnre <- expnre[,-1:-3]
expnre[1:2,1:6]
#AD loss reserve 受试者表达数据，命名expnlre
expnlre <- ralyl[,colnames(ralyl) %in% sublre$ID]
expnlre <- cbind(ralyl[,c(1,2,3)],expnlre)
rownames(expnlre) <- expnlre$SYMBOL; expnlre <- expnlre[,-1:-3]
expnlre[1:2,1:6]
#Normal aging 受试者表达数据，命名expnna
expnna <- ralyl[,colnames(ralyl) %in% subna$ID]
expnna <- cbind(ralyl[,c(1,2,3)],expnna)
rownames(expnna) <- expnna$SYMBOL; expnna <- expnna[,-1:-3]
expnna[1:2,1:6]
#MCI 受试者表达数据，命名expmci
expnmci <- ralyl[,colnames(ralyl) %in% submci$ID]
expnmci <- cbind(ralyl[,c(1,2,3)],expnmci)
rownames(expnmci) <- expnmci$SYMBOL; expnmci <- expnmci[,-1:-3]
expnmci[1:2,1:6]
#AD dementia 受试者表达数据，命名expnna
expnad <- ralyl[,colnames(ralyl) %in% subad$ID]
expnad <- cbind(ralyl[,c(1,2,3)],expnad)
rownames(expnad) <- expnad$SYMBOL; expnad <- expnad[,-1:-3]
expnad[1:2,1:6]

save(ralyl,chd8,subre,sublre,subna,submci,subad,
     expnre,expnlre,expnna,expnmci,expnad,clinical,
     follow,rexp3,idkey,
     file = "G:\\AD reserve omics\\Rproject\\GSE1297_BayesianNetwork\\RALYL validation\\ROSMAP RNAseq data sorting for RALYL validation.Rdata")

##--------------------------------------------------------
## RALYL 在NCI/MCI/AD中的表达；在reserve/loss reserve的表达
load("G:\\AD reserve omics\\Rproject\\GSE1297_BayesianNetwork\\RALYL validation\\ROSMAP RNAseq data sorting for RALYL validation.Rdata")

## RALYL高低表达的患者的cog decline情况
Xfl <- as.data.frame(t(ralyl)); summary(Xfl)  
LOW = 14.119; HIGH = 21.607 # 手动设置常量

# 分拣高低表达fl受试者的临床表型
for (i in 1:nrow(Xfl)) {
  if (Xfl$RALYL[i] <= LOW ) {
    Xfl$exp[i] <- "low"
  } else if (Xfl$RALYL[i] > HIGH) {
    Xfl$exp[i] <- "high"
  } else {
    Xfl$exp[i] <- NA
  }
}
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
mmse_dcl <- cbind(lmmse_dcl,hmmse_dcl)
mmse_dcl <- melt(mmse_dcl, value.name = "counts") # 出图部分程序
set.seed(1024); p1 <- ggplot(mmse_dcl,aes(x=variable,y=counts)) +
  geom_boxplot(aes(fill=variable)) +
  geom_jitter(shape=21, position=position_jitter(0.2)) +
  labs(title="MMSE decline during the follow-up",
       x="Expression levels of RALYL in ROSMAP RNAseq", y = "MMSE decline") +
  theme_bw() + guides(fill=FALSE) +
  theme(panel.grid=element_blank(),plot.title = element_text(hjust = 0.5),
        axis.text.x=element_text(size = 13),axis.text.y=element_text(size = 13));p1
my_cp <- list(c("high","low")) # 统计学检验部分
p1 + stat_compare_means(comparisons = my_cp) +
  stat_compare_means(method = "t.test",label.y = -2)

## 按照Braak评分分级后的RALYL表达情况
bk <- clinical[,c(1,10)]
for (i in 1:nrow(bk)) {
  for (j in 1:ncol(ralyl)) {
    if (bk$ID[i] %in% colnames(ralyl)[j] ) {
      bk$RALYL[i] <- ralyl[1,j]
    }
  }
}
for (i in 1:nrow(bk)) {
  if (is.na(bk$braak[i]) | is.na(bk$RALYL[i]) ) {
    bk <- bk[-i,]
  }
}

bk$braak = factor(bk$braak, levels=c('1','2','3','4','5','6')) #柱条顺序
set.seed(1024)
p1 <- ggplot(bk,aes(x=braak,y=RALYL)) + 
  geom_boxplot(aes(fill=braak)) +
  geom_jitter(shape=21, position=position_jitter(0.2)) +
  labs(title="RALYL expression levels and braak stage",
       x="Braak stage", y = "RALYL expression levels") +
  theme_bw() + guides(fill=FALSE) +
  theme(panel.grid=element_blank(),plot.title = element_text(hjust = 0.5),
        axis.text.x=element_text(size = 13),axis.text.y=element_text(size = 13));p1
my_cp <- list(c("1","6"),c("1","5"),c("1","4"),c("1","3"),c("1","2")) # 统计学检验
  stat_compare_means(method = "anova",label.x = 4, label.y = 40)

## AD/MCI/NA差异表达
abc1 <- matrix(numeric(71)[NA])
expnmci <- as.data.frame(rbind(t(expnmci),abc1))
abc2 <- matrix(numeric(72)[NA])
expnna <- as.data.frame(rbind(t(expnna),abc2)) 
expnad <- as.data.frame(t(expnad))
amn <- as.data.frame(cbind(expnad,expnmci,expnna))
colnames(amn) <- c("AD","MCI","NCI")
amn_melt <- melt(amn, value.name = "counts")

amn_melt$variable = factor(amn_melt$variable, levels=c('NCI','MCI','AD'))
set.seed(1024)
p1 <- ggplot(amn_melt, aes(x=variable,y=counts)) +
  geom_boxplot(aes(fill=variable)) +
  geom_jitter(shape=21, position=position_jitter(0.2)) +
  labs(title="Expression Level of RALYL in ROSMAP RNAseq ",
       x="Subjects Type", y = "RALYL gene expression levels", 
       fill = "Cell Type") +
  theme_bw() + guides(fill=FALSE) +
  theme(panel.grid=element_blank(),plot.title = element_text(hjust = 0.5),
        axis.text.x=element_text(size = 13),axis.text.y=element_text(size = 13));p1
my_cp <- list(c("AD","MCI"),c("AD","NCI"),c("MCI","NCI")) #统计学检验
p1 + stat_compare_means(comparisons = my_cp) +
  stat_compare_means(method = "anova",label.y =30)

## reserve/lossreserve差异表达

abc1 <- matrix(numeric(23)[NA])
expnre <- as.data.frame(rbind(t(expnre),abc1)) 
expnlre <- as.data.frame(t(expnlre)) 

rl <- as.data.frame(cbind(expnre,expnlre))
colnames(rl) <- c("reserve","loss reserve")
rl_melt <- melt(rl, value.name = "counts")
# rl_melt$counts = log(rl_melt$counts)
set.seed(1024)
p2 <- ggplot(rl_melt, aes(x=variable,y=counts)) +
  geom_boxplot(aes(fill=variable)) +
  geom_jitter(shape=21, position=position_jitter(0.2)) +
  labs(title="Expression Level of RALYL in reserve & loss_reserve",
       x="Subjects Type", y = "RALYL gene expression levels", 
       fill = "Cell Type") +
  theme_bw() + guides(fill=FALSE) +
  theme(panel.grid=element_blank(),plot.title = element_text(hjust = 0.5),
        axis.text.x=element_text(size = 13),axis.text.y=element_text(size = 13));p2
my_cp <- list(c("reserve","loss reserve"))
p2 + stat_compare_means(comparisons = my_cp) +
  stat_compare_means(method = "t.test",label.y =8)
