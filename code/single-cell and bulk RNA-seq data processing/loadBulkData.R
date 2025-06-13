
# 过滤低表达基因 TPM<0.01
# 基因表达数据预处理，包括下载公共数据集、Ensembl转换、"_PBMC"样本提取、删除LUAD中的健康样本和非冰冻组织样本
# 检测过滤离群样本


# 详细的bulk数据处理代码未给出，主要列出原始数据和处理后的数据
# 直接加载预处理好的Bulk数据
rawBulkPBMC = read.csv("bulkDataProcessing/rawBulkPBMC_TPM.csv", row.names = 1, check.names = F)
bulkPBMC = read.csv("bulkDataProcessing/bulkPBMC.csv", row.names = 1, check.names = F)
# rawBulkLUAD
load("bulkDataProcessing/rawBulkLUAD_TPM.RData")
# bulkLUAD
load("bulkDataProcessing/bulkLUAD.RData")

# 获取Bulk免疫基因表达矩阵
# 加载免疫基因集immuneGeneSet

immuneGeneSet = read.csv(file = "immuneGeneIdentification/immuneGeneSet.csv")
immuneGeneSet = immuneGeneSet$immuneGene
# 获得bulk normal PBMC 免疫基因表达矩阵
immuneGenePBMC = bulkPBMC[which(row.names(bulkPBMC) %in% immuneGeneSet),]
write.csv(immuneGenePBMC, file = "loadBulkData/immuneGenePBMC.csv", row.names = T)
# 获得bulk LUAD 免疫基因表达矩阵 
immuneGeneLUAD = bulkLUAD[which(row.names(bulkLUAD) %in% immuneGeneSet),]
write.csv(immuneGeneLUAD, file = "loadBulkDatan/immuneGeneLUAD.csv", row.names = T)

