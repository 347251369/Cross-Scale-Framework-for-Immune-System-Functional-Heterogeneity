# 考虑到Bulk高相关性免疫基因数目过高，需要划分成若干个模块进行贝叶斯分析
# PBMC
bayesGeneSplitPBMC = data.frame(ID = NULL, Description = NULL, geneSet = NULL, geneCount = NULL)
goIDs = c()
description = c()
geneSet = c()
geneCount = c()
for(i in 1:length(immuneGOInfo$ID)){
  goIDs = c(goIDs, immuneGOInfo$ID[i])
  description = c(description, immuneGOInfo$Description[i])
  geneAll = unlist(strsplit(immuneGOInfo$geneSet[i], "/"))
  geneAll = geneAll[which(geneAll %in% highImmuneGenePBMC)]
  geneSet = c(geneSet, paste(geneAll, collapse = "/"))
  geneCount = c(geneCount, length(geneAll))
}
bayesGeneSplitPBMC = data.frame(ID = goIDs, Description = description, geneSet = geneSet, geneCount = geneCount)
# 修改那些大通路（geneCount > 400）对应的基因集
smallPathGene = unique(unlist(lapply(strsplit(as.character(bayesGeneSplitPBMC$geneSet[which(bayesGeneSplitPBMC$geneCount <= 2)]), "/"), trimws)))
goodPathGene = unique(unlist(lapply(strsplit(as.character(bayesGeneSplitPBMC$geneSet[which(bayesGeneSplitPBMC$geneCount <= 400)]), "/"), trimws)))
bigPathID = "GO:0002443"
bigPathGene = unique(unlist(lapply(strsplit(as.character(bayesGeneSplitPBMC$geneSet[which(bayesGeneSplitPBMC$ID == bigPathID)]), "/"), trimws)))
bayesGeneSplitPBMC[which(bayesGeneSplitPBMC$ID == bigPathID),]$geneSet = paste(unique(c(bigPathGene[!(bigPathGene %in% goodPathGene)])), collapse = "/")
bayesGeneSplitPBMC[which(bayesGeneSplitPBMC$ID == bigPathID),]$geneCount = length(unique(c(bigPathGene[!(bigPathGene %in% goodPathGene)])))
# 调整了10条通路，0002376,0006955,0002520,0045321,0048534,0045087,0046649,0002250,0042110,0002443
# 保留geneCount>=3的通路对应的基因集用于贝叶斯分析
bayesGeneSplitPBMC = bayesGeneSplitPBMC[which(bayesGeneSplitPBMC$geneCount >= 3),]
save(bayesGeneSplitPBMC, file = "BayesianAnalysis/bayesGeneSplitPBMC.RData")

#LUAD
bayesGeneSplitLUAD = data.frame(ID = NULL, Description = NULL, geneSet = NULL, geneCount = NULL)
goIDs = c()
description = c()
geneSet = c()
geneCount = c()
for(i in 1:length(immuneGOInfo$ID)){
  goIDs = c(goIDs, immuneGOInfo$ID[i])
  description = c(description, immuneGOInfo$Description[i])
  geneAll = unlist(strsplit(immuneGOInfo$geneSet[i], "/"))
  geneAll = geneAll[which(geneAll %in% highImmuneGeneLUAD)]
  geneSet = c(geneSet, paste(geneAll, collapse = "/"))
  geneCount = c(geneCount, length(geneAll))
}
bayesGeneSplitLUAD = data.frame(ID = goIDs, Description = description, geneSet = geneSet, geneCount = geneCount)
# 修改那些大通路（geneCount > 400）对应的基因集
smallPathGene = unique(unlist(lapply(strsplit(as.character(bayesGeneSplitLUAD$geneSet[which(bayesGeneSplitLUAD$geneCount <= 2)]), "/"), trimws)))
goodPathGene = unique(unlist(lapply(strsplit(as.character(bayesGeneSplitLUAD$geneSet[which(bayesGeneSplitLUAD$geneCount <= 600)]), "/"), trimws)))
bigPathID = "GO:0006955"
bigPathGene = unique(unlist(lapply(strsplit(as.character(bayesGeneSplitLUAD$geneSet[which(bayesGeneSplitLUAD$ID == bigPathID)]), "/"), trimws)))
bayesGeneSplitLUAD[which(bayesGeneSplitLUAD$ID == bigPathID),]$geneSet = paste(unique(c(intersect(bigPathGene, smallPathGene), setdiff(bigPathGene, goodPathGene))), collapse = "/")
bayesGeneSplitLUAD[which(bayesGeneSplitLUAD$ID == bigPathID),]$geneCount = length(unique(c(intersect(bigPathGene, smallPathGene), setdiff(bigPathGene, goodPathGene))))
# 调整了两个通路0002376和0006955，最后保留geneCount>2的基因集进行贝叶斯分析
bayesGeneSplitLUAD = bayesGeneSplitLUAD[which(bayesGeneSplitLUAD$geneCount >= 3),]
save(bayesGeneSplitLUAD, file = "bayesianAnalysis/bayesGeneSplitLUAD.RData")

# 开始进行贝叶斯分析
library(bnlearn)
library(parallel)
# 加载Bulk高相关性免疫基因表达矩阵和划分后的结果
corMatrixPBMC = read.csv(file = "correlationAnalysis/corMatrixPBMC.csv", check.names = F, row.names = 1)
corMatrixLUAD = read.csv(file = "correlationAnalysis/corMatrixLUAD.csv", check.names = F, row.names = 1)
bayesGeneSplitPBMC = get(load(file = "bayesianAnalysis/bayesGeneSplitPBMC.RData"))
bayesGeneSplitLUAD = get(load(file = "bayesianAnalysis/bayesGeneSplitLUAD.RData"))
# PBMC
cl = makeCluster(15)
bayesResPBMC = list()
for(i in 1:length(bayesGeneSplitPBMC$geneSet)){
  datExpr = t(corMatrixPBMC[,unlist(strsplit(bayesGeneSplitPBMC$geneSet[i], "/"))])
  res = boot.strength(datExpr, algorithm = "hc", R=100, cluster = cl)
  bayesResPBMC[[i]] = res
}
save(bayesResPBMC, file = "bayesianAnalysis/bayesResPBMC.RData")
# LUAD
cl = makeCluster(15)
bayesResLUAD = list()
for(i in 1:length(bayesGeneSplitLUAD$geneSet)){
  datExpr = as.data.frame(t(corMatrixLUAD[unlist(strsplit(bayesGeneSplitLUAD$geneSet[i], "/")),]))
  res = boot.strength(datExpr, algorithm = "hc", R=100, cluster = cl)
  bayesResLUAD[[i]] = res
}
save(bayesResLUAD, file = "bayesianAnalysis/bayesResLUAD.RData")

# 处理bayesRes
bayesResProcessing = function(bayesRes){
  bayesGenePair = do.call("rbind",bayesRes)
  # 设定阈值，过滤掉低置信度的gene pair
  bayesGenePair = bayesGenePair[bayesGenePair$direction >= 0.5,]
  bayesGenePair = bayesGenePair[bayesGenePair$strength >= 0.7,]
  # 在确定gene pair都符合条件后，删除重复的gene pair
  # 根据from和to列判断重复行
  duplicated_rows = duplicated(subset(bayesGenePair, select = c("from", "to")))
  # 保留不重复的行
  bayesGenePair = bayesGenePair[!duplicated_rows, ]
  
  return(bayesGenePair)
}

# PBMC
bayesGenePairPBMC = bayesResProcessing(bayesResPBMC)
# 保存处理结果
write.csv(bayesGenePairPBMC, file = "bayesianAnalysis/bayesGenePairPBMC.csv", row.names = F)
# LUAD
bayesGenePairLUAD = bayesResProcessing(bayesResLUAD)
# 保存处理结果
write.csv(bayesGenePairLUAD, file = "bayesianAnalysis/bayesGenePairLUAD.csv", row.names = F)
