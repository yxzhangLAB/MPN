###-------------------------------------------------------
### GSE1297 Bayesian Network develop
### 此脚本用于BN建模，优化，评价
### Lastupdate:202002015
### copyright@Yixuan Zhang
###-------------------------------------------------------

rm(list = ls())
set.seed(1024)
library(bnlearn)
library(Rgraphviz)
library(pcalg)
library(deal)

##--------------------------------------------------------
## 并行计算
library(parallel)
clnum <- detectCores()
cl <- makeCluster(getOption("cl.cores", clnum)) #用所有的6核12线程计算

##--------------------------------------------------------
## 数据导入
load("GSE1297 BN inputdata.Rdata")
load("GSE1297 MMSE & NFTs & braak & reserve discreted.Rdata")

##--------------------------------------------------------
## bootstraping抽样与白名单设置
if (FALSE) 
  
  {
bn_1297 <- as.data.frame(bn_1297)
dbn_1297 <- discretize(bn_1297, method = "hartemink", breaks = 6,
                       idisc = "quantile",debug = FALSE)
dbn_1297 <- cbind(dbn_1297,MMSE,NFTs,braak,reserve)
tdbn_1297 <- as.data.frame(t(dbn_1297))

btsp_1297 <- sample(tdbn_1297, size = 10, replace = TRUE) #抽样10个，循环抽30次
for (i in 1:30) {
  abc1 <- sample(tdbn_1297, size = 10, replace = TRUE)
  btsp_1297 <- cbind(btsp_1297,abc1)
}
btsp_1297 <- as.data.frame(t(btsp_1297))
summary(btsp_1297);class(btsp_1297)

#设置白名单，NFTs & braak to MMSE  
wl <- data.frame(from=character(1), to=character(1)) #设置白名单，NFTs2MMSE
wl$from <- as.character(wl$from); wl$to <- as.character(wl$to)
wl[1,1] <- "M1"; wl[1,2] <- "MMSE"
wl[2,1] <- "M11"; wl[2,2] <- "reserve"
wl[3,1] <- "M11"; wl[3,2] <- "MMSE"
wl[4,1] <- "braak"; wl[4,2] <- "MMSE"
save(btsp_1297,wl,file = "GSE1297 bootstraping data for BN use.Rdata")
  }
##--------------------------------------------------------
## BN建模与调参： hc 方法
load("GSE1297 bootstraping data for BN use.Rdata")

set.seed(1024)
bn1297 <- hc(btsp_1297, whitelist = wl, score = "k2",debug = TRUE, #爬山法
              restart = 100, perturb = 1,maxp = 5, max.iter = Inf,optimized = TRUE) 
# 注：v1.3版本 restart = 1000
# 注：20200824 paper001修回过程中发现，一直以来使用的经典的BN结构来自于v1.1版本；
# 原始Rdata名称为：GSE1297 Bayesian network hc-k2 v1_1 .Rdata

bn.gr <- bn1297
hl2 = list(arcs = vstructs(bn.gr, arcs = TRUE),lwd = 3, col = "black")
gr <- graphviz.plot(bn.gr, highlight = hl2, layout = "dot",
                    main = "GSE1297 Bayesian Network by_hc-k2 v1.3",
                    shape = "rectangle")
mb(bn1297, "M11")
save(bn1297,gr,file = "GSE1297 Bayesian network hc-k2 v1_4 .Rdata")

load("GSE1297 Bayesian network hc-k2 v1_1 .Rdata")
# 获取网络nodes相关参数
nd <- "M11"
mb(bn1297, nd)
nbr(bn1297, nd)
parents(bn1297, nd)
children(bn1297, nd)
in.degree(bn1297, nd)
out.degree(bn1297, nd)
root.nodes(bn1297)
leaf.nodes(bn1297)
nnodes(bn1297)

##--------------------------------------------------------
## BN建模- 其他几种未采用方法：tabu、mmhc、boot

# bn1297 <- tabu(btsp_1297, start = NULL, whitelist = NULL, score ="k2", #禁忌搜索法
#      debug = TRUE, tabu = 100, max.tabu = 100, max.iter = 1000000, maxp = Inf,
#      optimized = TRUE)

# avg.boot <- mmhc(btsp_1297,whitelist = wl,debug = TRUE, #大小爬山法
#                  restrict.args = list(test="sp-mi",alpha=0.05),
#                  maximize.args =list(restart = 1000,perturb = 1,score = "k2",
#                                     maxp=Inf,optimized = TRUE)) 

# boot = boot.strength(data = btsp_1297, R = 100, algorithm = "hc",cluster = cl, 
#                      algorithm.args = list(score = "k2",iss = 100),debug = TRUE)
# 
# boot[(boot$strength > 0.8) & (boot$direction >= 0.5), ]
# plot(boot)
# avg.boot = averaged.network(boot, threshold = 0.5)
# 
# bn.gr <- avg.boot
# hl2 <- list(arcs = vstructs(bn.gr, arcs = TRUE),lwd = 4, col = "black")
# gr <- graphviz.plot(bn.gr, highlight = hl2, layout = "dot", 
#                     shape = "ellipse", main = "***",sub ="****")





