###-------------------------------------------------------
### GSE1297 BN data cleaning
### 此脚本用于BN建模数据整理,挑选各module的表达数据
### 构建包含病理、认知等的BN用数据
### Lastupdate:20200202
### copyright@Yixuan Zhang
###-------------------------------------------------------

rm(list = ls())
library(Rtsne)
library(bnlearn)
library(Rgraphviz)

##--------------------------------------------------------
## 数据导入
exp1297 <- read.csv("Exp matrix GSE1297 Normalized.csv",header = TRUE)
trait1297 <- read.csv("Trait GSE1297.csv",header = TRUE)

##--------------------------------------------------------
## 导入module信息
if (FALSE) 
  {
module <- as.data.frame(c("BLACK","BLUE","BROWN","GREEN","GREY","MAGENTA",
                         "PINK","PURPLE","RED","TURQUOISE","YELLOW"))
colnames(module) <- "MODULE"
module_l <- list("BLACK","BLUE","BROWN","GREEN","GREY","MAGENTA",
                 "PINK","PURPLE","RED","TURQUOISE","YELLOW")
for (i in 1:11){
  abc1=paste('25% SD geneinput', as.character(module[i,1]), 
          'module gene name.csv',sep=' ' )
  module_l[[i]]=read.csv(as.character(abc1))
}
# module加入表达信息
exp_black <- exp1297[exp1297$X %in% as.character(module_l[[1]]$modProbes),]
rownames(exp_black) <- exp_black$X
exp_black <- exp_black[,-1] #black
exp_blue <- exp1297[exp1297$X %in% as.character(module_l[[2]]$modProbes),]
rownames(exp_blue) <- exp_blue$X
exp_blue <- exp_blue[,-1] #blue
exp_brown <- exp1297[exp1297$X %in% as.character(module_l[[3]]$modProbes),]
rownames(exp_brown) <- exp_brown$X
exp_brown <- exp_brown[,-1] #brown
exp_green <- exp1297[exp1297$X %in% as.character(module_l[[4]]$modProbes),]
rownames(exp_green) <- exp_green$X
exp_green <- exp_green[,-1] #green
exp_grey <- exp1297[exp1297$X %in% as.character(module_l[[5]]$modProbes),]
rownames(exp_grey) <- exp_grey$X
exp_grey <- exp_grey[,-1] #grey
exp_magenta <- exp1297[exp1297$X %in% as.character(module_l[[6]]$modProbes),]
rownames(exp_magenta) <- exp_magenta$X
exp_magenta <- exp_magenta[,-1] #magenta
exp_pink <- exp1297[exp1297$X %in% as.character(module_l[[7]]$modProbes),]
rownames(exp_pink) <- exp_pink$X
exp_pink <- exp_pink[,-1] #pink
exp_purple <- exp1297[exp1297$X %in% as.character(module_l[[8]]$modProbes),]
rownames(exp_purple) <- exp_purple$X
exp_purple <- exp_purple[,-1] #purple
exp_red <- exp1297[exp1297$X %in% as.character(module_l[[9]]$modProbes),]
rownames(exp_red) <- exp_red$X
exp_red <- exp_red[,-1] #red
exp_turquoise <- exp1297[exp1297$X %in% as.character(module_l[[10]]$modProbes),]
rownames(exp_turquoise) <- exp_turquoise$X
exp_turquoise <- exp_turquoise[,-1] #turquoise
exp_yellow <- exp1297[exp1297$X %in% as.character(module_l[[11]]$modProbes),]
rownames(exp_yellow) <- exp_yellow$X
exp_yellow <- exp_yellow[,-1] #yellow

save(exp_black,exp_blue,exp_brown,exp_green,exp_grey,exp_magenta,
     exp_pink,exp_purple,exp_red,exp_turquoise,exp_yellow,
     file="GSE1297 25%geneinput all Module expression.Rdata")
 }

##--------------------------------------------------------
## t-SNE降维,获取每个受试者的1维表达矩阵
if (FALSE) 
  {
load("GSE1297 25%geneinput all Module expression.Rdata") 
set.seed(1024)
  
exp_black <- t(as.matrix(exp_black))
tsne_out <- Rtsne(exp_black,dims=1,pca=FALSE,perplexity=5,theta=0.0)
plot(tsne_out$Y)
dr_black <- tsne_out[["Y"]] #black
exp_blue <- t(as.matrix(exp_blue))
tsne_out <- Rtsne(exp_blue,dims=1,pca=FALSE,perplexity=5,theta=0.0)
plot(tsne_out$Y)
dr_blue <- tsne_out[["Y"]] #blue 
exp_brown <- t(as.matrix(exp_brown))
tsne_out <- Rtsne(exp_brown,dims=1,pca=FALSE,perplexity=5,theta=0.0)
plot(tsne_out$Y)
dr_brown <- tsne_out[["Y"]] #brown  
exp_green <- t(as.matrix(exp_green))
tsne_out <- Rtsne(exp_green,dims=1,pca=FALSE,perplexity=5,theta=0.0)
plot(tsne_out$Y)
dr_green <- tsne_out[["Y"]] #green   
exp_grey <- t(as.matrix(exp_grey))
tsne_out <- Rtsne(exp_grey,dims=1,pca=FALSE,perplexity=5,theta=0.0)
plot(tsne_out$Y)
dr_grey <- tsne_out[["Y"]] #grey   
exp_magenta <- t(as.matrix(exp_magenta))
tsne_out <- Rtsne(exp_magenta,dims=1,pca=FALSE,perplexity=5,theta=0.0)
plot(tsne_out$Y)
dr_magenta <- tsne_out[["Y"]] #magenta 
exp_pink <- t(as.matrix(exp_pink))
tsne_out <- Rtsne(exp_pink,dims=1,pca=FALSE,perplexity=5,theta=0.0)
plot(tsne_out$Y)
dr_pink <- tsne_out[["Y"]] #pink 
exp_purple <- t(as.matrix(exp_purple))
tsne_out <- Rtsne(exp_purple,dims=1,pca=FALSE,perplexity=5,theta=0.0)
plot(tsne_out$Y)
dr_purple <- tsne_out[["Y"]] #purple 
exp_red <- t(as.matrix(exp_red))
tsne_out <- Rtsne(exp_red,dims=1,pca=FALSE,perplexity=5,theta=0.0)
plot(tsne_out$Y)
dr_red <- tsne_out[["Y"]] #red 
exp_turquoise <- t(as.matrix(exp_turquoise))
tsne_out <- Rtsne(exp_turquoise,dims=1,pca=FALSE,perplexity=5,theta=0.0)
plot(tsne_out$Y)
dr_turquoise <- tsne_out[["Y"]] #turquoise 
exp_yellow <- t(as.matrix(exp_yellow))
tsne_out <- Rtsne(exp_yellow,dims=1,pca=FALSE,perplexity=5,theta=0.0)
plot(tsne_out$Y)
dr_yellow <- tsne_out[["Y"]] #yellow 

save(dr_black,dr_blue,dr_brown,dr_green,dr_grey,dr_magenta,
     dr_pink,dr_purple,dr_red,dr_turquoise,dr_yellow,
     file="GSE1297 25%geneinput module Dimensionality Reduction.Rdata")

}

##--------------------------------------------------------
## 临床表型数据的离散与清洗
MMSE <- as.matrix(as.numeric(trait1297$MMSE)); colnames(MMSE) <- "MMSE"
NFTs <- as.matrix(trait1297$NFT); colnames(NFTs) <- "NFTs"
braak <- as.matrix(trait1297$braak);colnames(braak) <- "braak"
reserve <- as.matrix(trait1297$reserve);colnames(reserve) <- "reserve"
# 根据轻、中、重度离散MMSE
for (i in 1:nrow(MMSE)) {
  if (MMSE[i,1] >= 26) {
    MMSE[i,1] <- "normal"
  } else if (MMSE[i,1] < 26 & MMSE[i,1] >= 21) {
    MMSE[i,1] <- "incipient"
  } else if (MMSE[i,1] < 21 & MMSE[i,1] >= 10) {
    MMSE[i,1] <- "moderate" 
  } else if (MMSE[i,1] <10) {
    MMSE[i,1] <- "severe"
  }
}
# 基于CERAD离散NFT
for (i in 1:nrow(NFTs)) {
  if (NFTs[i,1] <= 2) {
    NFTs[i,1] <- "stageA"
  } else if (NFTs[i,1] > 2 & NFTs[i,1] <=6 ) {
    NFTs[i,1] <- "stageB"
  } else if (NFTs[i,1] > 6) {
    NFTs[i,1] <- "stageC"
  }
}
# 离散braak
for (i in 1:nrow(braak)) {
  if (braak[i,1] == 1) {
    braak[i,1] <- "stage_Ⅰ"
  } else if (braak[i,1] == 2) {
    braak[i,1] <- "stage_Ⅱ"
  } else if (braak[i,1] == 3) {
    braak[i,1] <- "stage_Ⅲ"
  } else if (braak[i,1] == 4) {
    braak[i,1] <- "stage_Ⅳ"
  } else if (braak[i,1] == 5) {
    braak[i,1] <- "stage_Ⅴ"
  } else if (braak[i,1] == 6) {
    braak[i,1] <- "stage_Ⅵ"
  }
# 离散reserve
reserve <- as.matrix(as.factor(reserve)); colnames(reserve) <- "reserve"

}
save(MMSE,NFTs,braak,reserve,
     file = "GSE1297 MMSE & NFTs & braak & reserve discreted.Rdata")
##--------------------------------------------------------
## BN用matrix整理，命名：bn_1297
load("GSE1297 25%geneinput module Dimensionality Reduction.Rdata")
load("GSE1297 MMSE & NFTs & braak & reserve discreted.Rdata")

# 定义module与数字间的对应关系
abc1 <- c("dr_black","dr_blue","dr_brown","dr_green","dr_grey",
          "dr_magenta","dr_pink","dr_purple","dr_red","dr_turquoise","dr_yellow")
abc2 <- c("M1","M2","M3","M4","M5","M6","M7","M8","M9","M10","M11")
num_m <- cbind(abc1,abc2);colnames(num_m) <- c("module","number")

bn_1297 <- cbind(dr_black,dr_blue,dr_brown,dr_green,dr_grey,dr_magenta,
                 dr_pink,dr_purple,dr_red,dr_turquoise,dr_yellow)
colnames(bn_1297) <- c("M1","M2","M3","M4","M5","M6",
                       "M7","M8","M9","M10","M11")
save(bn_1297,num_m,file = "GSE1297 BN inputdata.Rdata")
