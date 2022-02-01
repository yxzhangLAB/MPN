###-------------------------------------------------------
### GSE1297 Bayesian Network Target genes Validation 
### 此脚本用真实世界数据对GSE197 m11的target genes外部验证
### 验证数据集包括：ROSMAP，GSE12685，GSE28146
### lastupdate:20200409
### copyright@Yixuan Zhang
###-------------------------------------------------------

library(reshape2)
library(ggplot2)
library(ggpubr)
library(pheatmap)
rm(list = ls())

##--------------------------------------------------------
##  ROSMAP Validation

# Data cleaning
load("G:\\AD reserve omics\\Rproject\\GSE1297_BayesianNetwork\\External_validation\\ROSMAP_EVdata.Rdata")
# 函数z-score：表达矩阵的归一化
zscore <- function(x) {
  dat <- x
  n=t(scale(t(dat)))
  n[n>2]=2 #限定上限，使表达量大于2的等于2
  n[n< -2]= -2 #限定下限，使表达量小于-2的等于-2
  n[1:4,1:4]
  dat <- as.data.frame(n)
  return(dat)
}
expnad <- zscore(expnad);expnmci <- zscore(expnmci);expnna <- zscore(expnna)
expnre <- zscore(expnre);expnlre <- zscore(expnlre)

# 靶基因导入，原CEP112在ROSMAP中symbol为"CCDC46"
tgs <- as.data.frame(c("RALYL","DYNLT3","EBI3","TSPYL5","CCDC46","SLC9A6"))
colnames(tgs) <- "hubG"
tgs$hubG %in% rownames(ROSMAPexpn)

# ROSMAP AD/MCI/NCI - boxplot
box_plot1 <- function(x) {
  dat <- as.character(x) 
  tpad <- as.data.frame(t(expnad[rownames(expnad) %in% dat,]))
  tpmci <- as.data.frame(t(expnmci[rownames(expnmci) %in% dat,]))
  tpna <- as.data.frame(t(expnna[rownames(expnna) %in% dat,])) 
  
  abc1 <- matrix(numeric(94)[NA])
  tpmci <- rbind(as.matrix(tpmci),abc1)
  abc2 <- matrix(numeric(84)[NA])
  tpna <- rbind(as.matrix(tpna),abc2)
  # AD/MCI/NA
  amn <- as.data.frame(cbind(tpad,tpmci,tpna))
  colnames(amn) <- c("AD","MCI","NCI")
  amn_melt <- melt(amn, value.name = "counts",na.rm = TRUE)
  # amn_melt$counts = log(amn_melt$counts)
  # 出图
  p1 <-ggplot(amn_melt, aes(x=variable,y=counts)) +
    geom_boxplot(aes(fill=variable)) +
    geom_jitter(shape=21, position=position_jitter(0.2)) +
    labs(title=paste("Expression Level of",dat),
         x="Subjects Type", y = "Expression levels", 
         fill = "Cell Type") +
    theme_bw() + guides(fill=FALSE) +
    theme(panel.grid=element_blank(),plot.title = element_text(hjust = 0.5),
          axis.text.x=element_text(size = 13),axis.text.y=element_text(size = 13))
  p1
  # 显著性检验：注意此方法存在bug，原图位置会被压缩，故只用数据不用图
  my_cp <- list(c("AD","MCI"),c("AD","NCI"),c("MCI","NCI"))
  p2 <- p1 + stat_compare_means(comparisons = my_cp) +
    stat_compare_means(method = "anova",label.y = 2.5)
  p2
  abc <-list(p1,p2)
  return(abc)
} # box plot drawing function
box_plot1("PSEN1")


##--------------------------------------------------------
## GSE12685，GSE28146 Validation

# Data Cleaning
# GSE12685,GSE28146数据导入，数据来源于2019年6月处理
load("G:\\AD reserve omics\\Rproject\\GSE1297_BayesianNetwork\\targetgene_validation\\step1_output.Rdata")
# 注：zscore后数据趋势被掩盖，故不再归一化
rownames(con12685) <- con12685$symbol;con12685 <- con12685[,-1:-2];
rownames(in12685) <- in12685$symbol;in12685 <- in12685[,-1:-2];
rownames(con28146) <- con28146$symbol;con28146 <- con28146[,-1:-2];
rownames(in28146) <- in28146$symbol;in28146 <- in28146[,-1:-2];
rownames(mo28146) <- mo28146$symbol;mo28146 <- mo28146[,-1:-2];
rownames(se28146) <- se28146$symbol;se28146 <- se28146[,-1:-2];
save(con12685,in12685,
     con28146,in28146,mo28146,se28146,
     Dx12685,Dx28146,exp12685,exp28146,
     file = "G:\\AD reserve omics\\Rproject\\GSE1297_BayesianNetwork\\targetgene_validation\\GSE12685_28146 all dataclean.Rdata")
rm(list = ls())
load("G:\\AD reserve omics\\Rproject\\GSE1297_BayesianNetwork\\targetgene_validation\\GSE12685_28146 all dataclean.Rdata")

# GSE12685 box plot
box_plot2 <- function(x) {
  dat <- as.character(x) 
  tpcon <- as.data.frame(t(con12685[rownames(con12685) %in% dat,]))
  tpin <- as.data.frame(t(in12685[rownames(in12685) %in% dat,]))
  
  abc1 <- matrix(numeric(2)[NA])
  tpin <- rbind(as.matrix(tpin),abc1)
  
  # AD/MCI/NA
  amn <- as.data.frame(cbind(tpcon,tpin))
  colnames(amn) <- c("Control","Incipient AD")
  amn_melt <- melt(amn, value.name = "counts",na.rm = TRUE)
  # amn_melt$counts = log(amn_melt$counts)
  # 出图
  p1 <-ggplot(amn_melt, aes(x=variable,y=counts)) +
    geom_boxplot(aes(fill=variable)) +
    geom_jitter(shape=21, position=position_jitter(0.2)) +
    labs(title=paste("Expression Level of",dat),
         x="Subjects Type", y = "Expression levels", 
         fill = "Cell Type") +
    theme_bw() + guides(fill=FALSE) +
    theme(panel.grid=element_blank(),plot.title = element_text(hjust = 0.5),
          axis.text.x=element_text(size = 13),axis.text.y=element_text(size = 13))
  p1
  # 显著性检验：注意此方法存在bug，原图位置会被压缩，故只用数据不用图
  my_cp <- list(c("Control","Incipient AD"))
  p2 <- p1 + stat_compare_means(comparisons = my_cp) +
    stat_compare_means(method = "t.test",label.y = 2.5)
  p2
  abc <-list(p1,p2)
  return(abc)
}
box_plot2("PSEN2")

# GSE28146 box plot
box_plot3 <- function(x) {
  dat <- as.character(x) 
  tpcon <- as.data.frame(t(con28146[rownames(con28146) %in% dat,]))
  tpin <- as.data.frame(t(in28146[rownames(in28146) %in% dat,]))
  tpmo <- as.data.frame(t(mo28146[rownames(mo28146) %in% dat,]))
  tpse <- as.data.frame(t(se28146[rownames(se28146) %in% dat,]))
  
  abc1 <- matrix(numeric(1)[NA])
  tpin <- rbind(as.matrix(tpin),abc1)
  tpse <- rbind(as.matrix(tpse),abc1)
  
  amn <- as.data.frame(cbind(tpcon,tpin,tpmo,tpse))
  colnames(amn) <- c("Control","Incipient AD","Moderate AD","Severe AD")
  amn_melt <- melt(amn, value.name = "counts",na.rm = TRUE)
  # amn_melt$counts = log(amn_melt$counts)
  # 出图
  p1 <-ggplot(amn_melt, aes(x=variable,y=counts)) +
    geom_boxplot(aes(fill=variable)) +
    geom_jitter(shape=21, position=position_jitter(0.2)) +
    labs(title=paste("Expression Level of",dat),
         x="Subjects Type", y = "Expression levels", 
         fill = "Cell Type") +
    theme_bw() + guides(fill=FALSE) +
    theme(panel.grid=element_blank(),plot.title = element_text(hjust = 0.5),
          axis.text.x=element_text(size = 13),axis.text.y=element_text(size = 13))
  p1
  # 显著性检验：注意此方法存在bug，原图位置会被压缩，故只用数据不用图
  my_cp <- list(c("Control","Incipient AD"),c("Control", "Moderate AD"),
                c("Control","Severe AD"))
  p2 <- p1 + stat_compare_means(comparisons = my_cp) +
    stat_compare_means(method = "anova",label.y = 2.5) +
    stat_compare_means(symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1),symbols = c("****", "***", "**", "*", "ns")))
  p2
  abc <-list(p1,p2)
  return(abc)
}
box_plot3("EBI3")
