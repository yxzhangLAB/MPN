###-------------------------------------------------------
### GSE1297 Bayesian Network Module11 target gene
### 此脚本用于构建M11内部网络，并找到target gene
### lastupdate:20200325
### copyright@Yixuan Zhang
###-------------------------------------------------------

rm(list = ls())
library(bnlearn)
library(Rgraphviz)
library(beepr)

##--------------------------------------------------------
## 并行计算
library(parallel)
clnum <- detectCores()
cl <- makeCluster(getOption("cl.cores", clnum)) #用所有的6核12线程计算

##--------------------------------------------------------
## 载入Module11的symbol及表达矩阵
m11in <- read.csv("G:\\AD reserve omics\\Rproject\\GSE1297_BayesianNetwork\\targetgene\\Module11 gene symbol.csv",
                  header = FALSE)
m11in <- m11in[m11in$V3 ==1,]
m11exp <- read.csv("G:\\AD reserve omics\\Rproject\\GSE1297_BayesianNetwork\\targetgene\\GSE1297 Module11 gene expression matrix.csv")
m11exp <- m11exp[m11exp$X %in% m11in$V2,]
rownames(m11exp) <- m11exp$X

##--------------------------------------------------------
## M11基因网络构建
set.seed(1024)
mbn <- as.data.frame(t(m11exp[,-1]))
dbn <- discretize(mbn, method = "hartemink", breaks = 6,
                  idisc = "quantile",debug = TRUE)
# 禁忌搜索法
set.seed(1024)
bnabc <- tabu(dbn, start = NULL, whitelist = NULL, score ="k2", 
               debug = TRUE, tabu = 100, max.tabu = 100, max.iter = 1000000, maxp = Inf,
               optimized = TRUE);beep()
# 绘图
bn.gr <- bnabc
hl2 <- list(arcs = vstructs(bn.gr, arcs = TRUE),lwd = 2, col = "black")
gr <- graphviz.plot(bn.gr, highlight = hl2, layout = "dot",
                    shape = "ellipse", 
                    main = "Module_11 inner Gene Bayesian Network",
                    sub ="Algorithm: tabu search with k2 score")
m11bn <- bnabc[["arcs"]]
# 数据导出
write.csv(m11bn,file="G:\\AD reserve omics\\Rproject\\GSE1297_BayesianNetwork\\targetgene\\GSE1297 M11 inner Bayesian Network.csv")
save(dbn,bn.gr,file="G:\\AD reserve omics\\Rproject\\GSE1297_BayesianNetwork\\targetgene\\GSE1297 Module 11 inner bayeian Network tabu k2.Rdata")

