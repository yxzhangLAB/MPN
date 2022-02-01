###-------------------------------------------------------
### GSE1297 Bayesian Network in Module_11
### 此脚本用于M11内部基因网络的构建
### Lastupdate:20200224
### copyright@Yixuan Zhang
###-------------------------------------------------------

rm(list = ls())
library(bnlearn)
library(Rgraphviz)

##--------------------------------------------------------
## 并行计算
library(parallel)
clnum <- detectCores()
cl <- makeCluster(getOption("cl.cores", clnum)) #用所有的6核12线程计算

##--------------------------------------------------------
## 载入所有module的基因数据
load("GSE1297 25%geneinput all Module expression.Rdata")

##--------------------------------------------------------
## M11基因网络构建

mbn <- as.data.frame(t(exp_yellow))
dbn <- discretize(mbn, method = "hartemink", breaks = 4,
                       idisc = "quantile",debug = TRUE)
# 禁忌搜索法
set.seed(1024)
bn1297 <- tabu(dbn, start = NULL, whitelist = NULL, score ="k2", 
     debug = TRUE, tabu = 100, max.tabu = 100, max.iter = 1000000, maxp = Inf,
     optimized = TRUE)
# 绘图
bn.gr <- bn1297
hl2 <- list(arcs = vstructs(bn.gr, arcs = TRUE),lwd = 2, col = "black")
gr <- graphviz.plot(bn.gr, highlight = hl2, layout = "dot",
                    shape = "ellipse", 
                    main = "Module_11 inner Gene Bayesian Network",
                    sub ="Algorithm: tabu search with k2 score")
m11bn <- bn1297[["arcs"]]

save(mbn,gr,m11bn,file = "Module 11 inner bayeian Network tabu k2.Rdata")
