# 加载Bulk免疫基因表达数据
immuneGenePBMC = read.csv("loadBulkData/immuneGenePBMC.csv", row.names = 1, check.names = F)
immuneGeneLUAD = read.csv("loadBulkData/immuneGeneLUAD.csv", row.names = 1, check.names = F)

# 利用WGCNA包进行相关性分析
library(WGCNA)

corAnalysis = function(immuneGeneMatrix){
  corMatrix = cor(t(immuneGeneMatrix))
  # 可视化相关性网络
  cyt = WGCNA::exportNetworkToCytoscape(corMatrix, weighted = TRUE, threshold = 0.7,
                                        nodeNames = colnames(t(immuneGeneMatrix)))
  # 从网络中抽取出gene pairs信息
  genePairs = cyt[["edgeData"]][,c(1:3)]
  
  return(genePairs)
}



# PBMC相关性分析
genePairsPBMC = corAnalysis(immuneGenePBMC)
write.csv(genePairsPBMC, "correlationAnalysis/corGenePairsPBMC.csv",row.names = F)
# 从genePairsPBMC中抽取高相关性的免疫基因
highImmuneGenePBMC = unique(c(genePairsPBMC$fromNode, genePairsPBMC$toNode))
corMatrixPBMC = immuneGenePBMC[highImmuneGenePBMC,]
write.csv(highImmuneGenePBMC, "correlationAnalysis/highImmuneGenePBMC.csv")
write.csv(corMatrixPBMC, "correlationAnalysis/corMatrixPBMC.csv")

# LUAD相关性分析
genePairsLUAD = corAnalysis(immuneGeneLUAD)
write.csv(genePairsLUAD, "correlationAnalysis/corGenePairsLUAD.csv")
# 从genePairsLUAD中抽取高相关性的免疫基因
highImmuneGeneLUAD = unique(c(genePairsLUAD$fromNode, genePairsLUAD$toNode))
corMatrixLUAD = immuneGeneLUAD[highImmuneGeneLUAD,]
write.csv(highImmuneGeneLUAD, "correlationAnalysis/highImmuneGeneLUAD.csv")
write.csv(corMatrixLUAD, "correlationAnalysis/corMatrixLUAD.csv")

