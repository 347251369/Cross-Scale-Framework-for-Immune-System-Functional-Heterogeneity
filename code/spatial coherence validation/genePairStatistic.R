# 加载genePair的ratio值-------------------------------
resTableLUAD = read.csv(file = "spatialDataProcess/genePairStatistic/resTableLUAD.csv",
                        check.names = FALSE, row.names = 1)
goodRowsLUAD = read.csv(file = "spatialDataProcess/genePairStatistic/goodRowsLUAD.csv",
                        check.names = FALSE, row.names = 1)
# 筛选出resTableLUAD和goodRowsLUAD中两个节点的genePairs--------------------
targetRow = c()
for(row in row.names(resTableLUAD)){
  if(length(unlist(str_split(row, "_"))) == 2){
    targetRow = c(targetRow, row)
  }
}
resTableLUAD = resTableLUAD[targetRow,]

targetRow = c()
for(row in row.names(goodRowsLUAD)){
  if(length(unlist(str_split(row, "_"))) == 2){
    targetRow = c(targetRow, row)
  }
}
goodRowsLUAD = goodRowsLUAD[targetRow,]
table(row.names(resTableLUAD) %in% row.names(goodRowsLUAD))
tmp = resTableLUAD[!(row.names(resTableLUAD) %in% row.names(goodRowsLUAD)),]

# 整理genePairStatisticResLUAD1020文件-----------------------------
library(readxl)
sheet_names = excel_sheets("spatialDataProcess/genePairStatistic/genePairStatisticResLUAD1020.xlsx")
# 读取所有sheet
sheets = lapply(sheet_names, function(sheet) {
  read_excel("spatialDataProcess/genePairStatistic/genePairStatisticResLUAD1020.xlsx", sheet = sheet)
})
allSheets = do.call(rbind, sheets)
allSheets_ls = list()
for(i in 1:nrow(allSheets)){
  allSheets_ls[[i]] = c(allSheets$from[i], allSheets$to[i])
}
allSheets_ls = lapply(allSheets_ls, function(x){x = x[order(x)]})
allSheets_ls = unique(allSheets_ls)
allSheets_ls = unlist(lapply(allSheets_ls, function(x){str_c(x, collapse = "_")}))
allSheets_ls.ratio = resTableLUAD[allSheets_ls,]
allSheets_ls.ratio = goodRowsLUAD
allSheets_ls.maxRatio = as.data.frame(matrix(data = 0, nrow = nrow(allSheets_ls.ratio), ncol = 2))
names(allSheets_ls.maxRatio) = c("genePair", "maxRatio")
for(i in 1:nrow(allSheets_ls.ratio)){
  allSheets_ls.maxRatio$genePair[i] = row.names(allSheets_ls.ratio)[i]
  allSheets_ls.maxRatio$maxRatio[i] = max(allSheets_ls.ratio[i,])
}

# 构造clusterTOcluster---------------------------------------------
library(stringr)
library(dplyr)
library(tidyr)
library(matrixStats) # rowMaxs
# 找到两个biomarker之间的bayesGenePair
# 输入两个networkBiomarker
# 输出networkBiomarker之间的bayesGenePair
getBayesGenePair = function(biomarker1, biomarker2, bayesGenePair){
  biomarker1_gene = unlist(str_split(biomarker1, "_"))
  biomarker2_gene = unlist(str_split(biomarker2, "_"))
  genePair = c()
  for(gene1 in biomarker1_gene){
    for(gene2 in biomarker2_gene){
      gene = c(gene1, gene2)
      gene = gene[order(gene)]
      genePair = c(genePair, str_c(gene, collapse = "_"))
    }
  }
  genePair = genePair[genePair %in% bayesGenePair]
  return(genePair)
}
# 设置Ratio阈值用于筛选不共表达genePair
ratioThreshold = 1
# 找到符合阈值条件的genePair
# 输入两个biomarker之间的bayesGenePair
# 输出符合阈值条件的bayesGenePair
getTargetGenePair = function(bayesGenePair, ratioThreshold, resTableLUAD){
  resTableLUADI = resTableLUAD[bayesGenePair,]
  # 找出满足ratioThreshold的bayesGenePair
  rowMaxs = rowMaxs(as.matrix(resTableLUADI), useNames = FALSE)
  if(length(which(rowMaxs <= ratioThreshold) != 0)){
    return(row.names(resTableLUADI)[which(rowMaxs <= ratioThreshold)])
  }else{
    return(NULL)
  }
}
networkBiomarkerSummary = read.csv(file = "spatialDataProcess/genePairStatistic/networkBiomarkerSummary.csv",
                                   check.names = FALSE)
# 建立clusterTOcluster矩阵
clusterTOcluster = matrix(data = 0, nrow = nrow(networkBiomarkerSummary), ncol = nrow(networkBiomarkerSummary))
for(i in 1:nrow(networkBiomarkerSummary)){
  for(j in 1:nrow(networkBiomarkerSummary)){
    # 相同cluster之间的关系定义为2
    if(i == j){
      clusterTOcluster[i,j] = 2
      next
    }
    # 定义不同cluster之间的关系
    biomarker1 = networkBiomarkerSummary$networkBiomarker[i]
    biomarker2 = networkBiomarkerSummary$networkBiomarker[j]
    # 找出两个biomarker之间存在的bayesGenePair
    bayesGenePairI = getBayesGenePair(biomarker1, biomarker2, bayesGenePair)
    if(length(bayesGenePairI) == 0){
      clusterTOcluster[i,j] = 0
    }else{
      # 找出两个biomarker之间满足阈值的bayesGenePair
      targetGenePairI = getTargetGenePair(bayesGenePairI, ratioThreshold, resTableLUAD)
      # 通过targetGenePairI判断cluster之间的关系
      if(length(targetGenePairI) != 0){
        # 不同cluster之间有因果关系定义为1
        clusterTOcluster[i,j] = 1
      }else{
        # 不同cluster之间无因果关系定义为0
        clusterTOcluster[i,j] = 0
      }
    }
  }
}
colnames(clusterTOcluster) = paste(networkBiomarkerSummary$cellType, networkBiomarkerSummary$cluster, sep = "_")
row.names(clusterTOcluster) = paste(networkBiomarkerSummary$cellType, networkBiomarkerSummary$cluster, sep = "_")
write.csv(clusterTOcluster, file = paste("spatialDataProcess/genePairStatistic/clusterTOcluster(threshold=", ratioThreshold, ").csv", sep = ""))

# 统计networkBiomarker之间的关系------------------------------------------------
# 数据框，包括genePair、cluster1、cluster2
networkBiomarkerSummary = read.csv(file = "spatialDataProcess/genePairStatistic/networkBiomarkerSummary.csv",
                                   check.names = FALSE)
genePair_cluster_info = list()
for(i in 1:nrow(networkBiomarkerSummary)){
  for(j in 1:nrow(networkBiomarkerSummary)){
    # 相同cluster之间的关系定义为2
    if(i == j){
      next
    }
    # 定义不同cluster之间的关系
    biomarker1 = networkBiomarkerSummary$networkBiomarker[i]
    cluster1 = paste(networkBiomarkerSummary$cellType[i], networkBiomarkerSummary$cluster[i], sep = "_")
    biomarker2 = networkBiomarkerSummary$networkBiomarker[j]
    cluster2 = paste(networkBiomarkerSummary$cellType[j], networkBiomarkerSummary$cluster[j], sep = "_")
    # 找出两个biomarker之间存在的bayesGenePair
    bayesGenePairI = getBayesGenePair(biomarker1, biomarker2, bayesGenePair)
    if(length(bayesGenePairI) == 0){
      clusterTOcluster[i,j] = 0
    }else{
      # 找出两个biomarker之间满足阈值的bayesGenePair
      targetGenePairI = bayesGenePairI[which(bayesGenePairI %in% row.names(resTableLUAD))]
      # 通过targetGenePairI判断cluster之间的关系
      if(length(targetGenePairI) != 0){
        # 不同cluster之间有因果关系定义为1
        genePair = targetGenePairI
        clusterName1 = rep(cluster1, length(genePair))
        clusterName2 = rep(cluster2, length(genePair))
        loc = paste(i, j, sep = "_")
        genePair_cluster_info[[loc]] = data.frame(genePair = genePair, clusterName1 = clusterName1,
                                                  clusterName2 = clusterName2)
      }
    }
  }
}
genePair_cluster_info = do.call(rbind, genePair_cluster_info)
# 输入genePair，计算在cluster中的ratio
# 问题转换成如何在cluster中计算ratio
getGenePairRatio_cluster = function(genePair, cellType, cluster, expCounts){
  cellClusterRes = read.csv(file = paste("spatialDataProcess/cellClusterResL/cellClusterRes_", str_replace_all(cellType, " ", "_"), ".csv", sep = ""),
                           check.names = FALSE, row.names = 1)
  cell = row.names(cellClusterRes)[cellClusterRes$clusterRes == cluster]
  gene = unlist(str_split(genePair, "_"))
  expCountsI = expCounts[cell, gene]
  cellCount = nrow(expCountsI)
  targetCell = c()
  for(row in 1:nrow(expCountsI)){
    cellExp = expCountsI[row,]
    if(length(which(cellExp > 0)) == length(gene)){
      targetCell = c(targetCell, row)
    }
  }
  ratio = length(targetCell) / cellCount * 100
  return(ratio)
}
# 计算ratio
counts = as.matrix(luad@assays[["RNA"]]@counts)
counts = t(counts)
allGene = unique(c(bayesGenePairLUAD$from, bayesGenePairLUAD$to)) # 855
allGene = allGene[which(allGene %in% colnames(counts))] # 811
expCounts = counts[,allGene]
rm(counts)
# cluster数为1且细胞数前后变化不大的亚群
cellTypeSet.1 = c("NK", "cDC", "pre-B cells", "Plasma cells", "Neutrophils", "Macrophages")
# 循环遍历genePair_cluster_info
genePair_cluster_info$ratio1 = 0
genePair_cluster_info$ratio2 = 0
for(i in 1:nrow(genePair_cluster_info)){
  genePair = genePair_cluster_info$genePair[i]
  cellType1 = unlist(str_split(genePair_cluster_info$clusterName1[i],"_"))[1]
  cluster1 = unlist(str_split(genePair_cluster_info$clusterName1[i],"_"))[2]
  cellType2 = unlist(str_split(genePair_cluster_info$clusterName2[i],"_"))[1]
  cluster2 = unlist(str_split(genePair_cluster_info$clusterName2[i],"_"))[2]
  if(cellType1 %in% cellTypeSet.1){
    genePair_cluster_info$ratio1[i] = resTableLUAD[genePair, cellType1]
  }else{
    genePair_cluster_info$ratio1[i] = getGenePairRatio_cluster(genePair, cellType1, cluster1, expCounts)
  }
  if(cellType2 %in% cellTypeSet.1){
    genePair_cluster_info$ratio2[i] = resTableLUAD[genePair, cellType2]
  }else{
    genePair_cluster_info$ratio2[i] = getGenePairRatio_cluster(genePair, cellType2, cluster2, expCounts)
  }
  print(i)
}
genePair_cluster_info$strength = 0
# 利用bayesGenePairLUAD为genePair_cluster_info中的strength赋值
for(i in 1:nrow(bayesGenePairLUAD)){
  gene1 = bayesGenePairLUAD$from[i]
  gene2 = bayesGenePairLUAD$to[i]
  strength = bayesGenePairLUAD$strength[i]
  genePair = c(gene1, gene2)
  genePair = paste(genePair[order(genePair)], collapse = "_")
  loc = which(genePair_cluster_info$genePair %in% genePair)
  if(length(loc) != 0){ 
    genePair_cluster_info$strength[loc] = strength
  }
}
write.csv(genePair_cluster_info, file = "spatialDataProcess/genePairStatistic/genePair_cluster_info.csv")

# 设置共表达阈值，利用genePair_cluster_info获得clusterTOcluster--------------------
genePair_cluster_info = read.csv(file = "spatialDataProcess/genePairStatistic/genePair_cluster_info.csv",
                                 check.names = FALSE, row.names = 1)
# 利用ratioThreshold和strengthThreshold筛选genePair_cluster_info---------------
ratioThreshold = 5
strengthThreshold = 0.7
targetGenePair_cluster_info = genePair_cluster_info[genePair_cluster_info$ratio1 <= ratioThreshold,]
targetGenePair_cluster_info = targetGenePair_cluster_info[targetGenePair_cluster_info$ratio2 <= ratioThreshold,] # 1330
targetGenePair_cluster_info = targetGenePair_cluster_info[targetGenePair_cluster_info$strength >= strengthThreshold,] # 1330
# 基于targetGenePair_cluster_info建立clusterTOcluster矩阵----------------------
clusterTOcluster = matrix(data = 0, nrow = nrow(networkBiomarkerSummary), ncol = nrow(networkBiomarkerSummary))
for(i in 1:nrow(networkBiomarkerSummary)){
  for(j in 1:nrow(networkBiomarkerSummary)){
    # 相同cluster之间的关系定义为2
    if(i == j){
      clusterTOcluster[i,j] = 2
      next
    }
    # 定义不同cluster之间的关系
    clusterName1 = paste(networkBiomarkerSummary$cellType[i], networkBiomarkerSummary$cluster[i], sep = "_")
    clusterName2 = paste(networkBiomarkerSummary$cellType[j], networkBiomarkerSummary$cluster[j], sep = "_")
    loc1 = which(targetGenePair_cluster_info$clusterName1 %in% clusterName1)
    loc2 = which(targetGenePair_cluster_info$clusterName2 %in% clusterName2)
    if(length(intersect(loc1, loc2)) == 0){
      clusterTOcluster[i,j] = 0
    }else{
      clusterTOcluster[i,j] = 1
    }
  }
}

colnames(clusterTOcluster) = paste(networkBiomarkerSummary$cellType, networkBiomarkerSummary$cluster, sep = "_")
row.names(clusterTOcluster) = paste(networkBiomarkerSummary$cellType, networkBiomarkerSummary$cluster, sep = "_")
write.csv(clusterTOcluster, file = paste("spatialDataProcess/genePairStatistic/clusterTOcluster(threshold=",
                                         ratioThreshold, "_", strengthThreshold, ").csv", sep = ""))

# 基于targetGenePair_cluster_info构建celltypeTOcelltype---------------------
cellTypeSet = c("memory CD4+ T cells", "memory CD8+ T cells", "naive CD4+ T cells", "naive CD8+ T cells",
                "T helper cells", "T effector cells", "NKT", "MAIT", "Tregs",
                "NK", "Plasma cells", "pre-B cells", "cDC", "pDC", "Monocytes",
                "Basophils", "Neutrophils", "Macrophages")
# 将targetGenePair_cluster_info转换成targetGenePair_celltype_info
targetGenePair_celltype_info = targetGenePair_cluster_info
celltype_cluster1 = unlist(str_split(targetGenePair_celltype_info$clusterName1, "_"))
celltype_cluster2 = unlist(str_split(targetGenePair_celltype_info$clusterName2, "_"))
targetGenePair_celltype_info$clusterName1 =  celltype_cluster1[celltype_cluster1 %in% cellTypeSet]
targetGenePair_celltype_info$clusterName2 =  celltype_cluster2[celltype_cluster2 %in% cellTypeSet]
colnames(targetGenePair_celltype_info)[2:3] = c("celltype1", "celltype2")
# 开始构建celltypeTOcelltype
celltypeTOcelltype = matrix(data = 0, nrow = length(cellTypeSet), ncol = length(cellTypeSet))
colnames(celltypeTOcelltype) = cellTypeSet
row.names(celltypeTOcelltype) = cellTypeSet
for(i in 1:length(cellTypeSet)){
  for(j in 1:length(cellTypeSet)){
    # 相同cluster之间的关系定义为2
    if(i == j){
      celltypeTOcelltype[i,j] = 2
      next
    }
    # 定义不同cluster之间的关系
    celltype1 = cellTypeSet[i]
    celltype2 = cellTypeSet[j]
    loc1 = which(targetGenePair_celltype_info$celltype1 %in% celltype1)
    loc2 = which(targetGenePair_celltype_info$celltype2 %in% celltype2)
    if(length(intersect(loc1, loc2)) == 0){
      celltypeTOcelltype[i,j] = 0
    }else{
      celltypeTOcelltype[i,j] = 1
    }
  }
}
write.csv(celltypeTOcelltype, file = paste("spatialDataProcess/genePairStatistic/celltypeTOcelltype(threshold=",
                                           ratioThreshold, "_", strengthThreshold, ").csv", sep = ""))





