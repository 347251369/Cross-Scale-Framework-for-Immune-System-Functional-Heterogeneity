library(stringr)

# 基于stereoSeurat注释结果实现networkBiomarker投射------------------------------
# 导入LUAD的networkBiomarker
networkBiomarkerSummary = read.csv(file = "spatialDataProcess/genePairStatistic/networkBiomarkerSummary.csv",
                                   check.names = FALSE)
cellTypeSet = c("memory CD4+ T cells","memory CD8+ T cells","naive CD4+ T cells",
                "naive CD8+ T cells","T helper cells","T effector cells","NKT","MAIT","Tregs",
                "NK","Plasma cells","pre-B cells","cDC","pDC","Monocytes",
                "Basophils","Neutrophils","Macrophages")
# 获得组成networkBiomarker的基因
bioGene = unique(unlist(str_split(networkBiomarkerSummary$networkBiomarker, "_"))) # 127
sampleSet = c("C04143A3", "C04140A5", "C04145D6", "C04143E4")
sampleID = sampleSet[2]
# 加载空间转录组信息
metaData = stereoSeurat@meta.data[,c("x", "y", "cellType_subset")]
expressionRes = LayerData(stereoSeurat, assay = "RNA", layer = "counts")
expressionRes = t(expressionRes)
targetBioGene = bioGene[which(bioGene %in% colnames(expressionRes))] # bioGene中的空转基因
expressionRes = expressionRes[, targetBioGene]
# 目标细胞亚群
targetCelltypeSet = unique(metaData$cellType_subset[metaData$cellType_subset %in% cellTypeSet])
metaData$bins = row.names(metaData)
targetBin = metaData$bins[metaData$cellType_subset %in% targetCelltypeSet]
metaData = metaData[targetBin,]
positionRes = metaData[,c("x", "y")]
expressionRes = expressionRes[targetBin,]
networkBioExp = list()
for(cellType in unique(metaData$cellType_subset)){
  targetBinI = metaData$bins[metaData$cellType_subset %in% cellType]
  expressionResI = expressionRes[targetBinI,]
  networkBiomarkerI = networkBiomarkerSummary$networkBiomarker[networkBiomarkerSummary$cellType %in% cellType]
  networkBioExpI = matrix(data=0, nrow = length(targetBinI), ncol = length(networkBiomarkerI))
  row.names(networkBioExpI) = targetBinI
  colnames(networkBioExpI) = networkBiomarkerI
  for(bin in targetBinI){
    #exp = expressionResI[bin,]
    for(networkBiomarker in networkBiomarkerI){
      networkBiomarker_gene = unlist(str_split(networkBiomarker, "_"))
      loc = which(colnames(expressionResI) %in% networkBiomarker_gene)
      if(length(loc) == length(networkBiomarker_gene)){
        exp = expressionResI[bin, loc]
        if(all(exp > 0)){
          networkBioExpI[bin, networkBiomarker] = 1
        }
      }
    }
  }
  networkBioExp[[cellType]] = as.data.frame(networkBioExpI)
  print(cellType)
}
# 利用networkBioExp中的统计结果对细胞进行分类
binClassificationRes = list()
for(cellType in names(networkBioExp)){
  networkBioExpI = networkBioExp[[cellType]]
  loc = which(colSums(networkBioExpI) > 0)
  if(length(loc) > 0){
    for(loci in loc){
      targetBinI = row.names(networkBioExpI)[which(networkBioExpI[,loci] > 0)]
      networkBiomarkerI = colnames(networkBioExpI)[loci]
      binClassificationRes[[cellType]][[networkBiomarkerI]] = targetBinI
    }
  }
}
# 优化binClassificationRes，如果同一亚群中networkBiomarker之间存在子集关系，需要做进一步处理
binClassificationRes.opt = list()
for(cellType in names(binClassificationRes)){
  binClassificationResI = binClassificationRes[[cellType]]
  for(networkBiomarkerI in names(binClassificationResI)){
    for(networkBiomarkerJ in names(binClassificationResI)){
      networkBiomarkerI_gene = unlist(str_split(networkBiomarkerI, "_"))
      networkBiomarkerJ_gene = unlist(str_split(networkBiomarkerJ, "_"))
      # I是J的子集
      if(all(networkBiomarkerI_gene %in% networkBiomarkerJ_gene) && !all(networkBiomarkerJ_gene %in% networkBiomarkerI_gene)){
        binClassificationResI[[networkBiomarkerI]] = setdiff(binClassificationResI[[networkBiomarkerI]], binClassificationResI[[networkBiomarkerJ]])
      }
    }
    if(length(binClassificationResI[[networkBiomarkerI]]) > 0){
      binClassificationRes.opt[[cellType]][[networkBiomarkerI]] = binClassificationResI[[networkBiomarkerI]]
    }
  }
  binClassificationRes[[cellType]] = binClassificationResI
}
binClassificationRes = binClassificationRes.opt
# 删除C04145D6中的BTLA_CCR7_CD247_CXCL13_ICOS_ITK对应的细胞
# binClassificationRes[["naive CD4+ T cells"]][["BTLA_CCR7_CD247_CXCL13_ICOS_ITK"]] = NULL
# 获得binPositionRes
binPositionRes = list()
for(cellType in names(binClassificationRes)){
  binClassificationResI = binClassificationRes[[cellType]]
  for(networkBiomarkerI in names(binClassificationResI)){
    targetBinI = binClassificationResI[[networkBiomarkerI]]
    binPositionRes[[cellType]][[networkBiomarkerI]] = positionRes[targetBinI,]
  }
}
# 保存三项信息：networkBioExp、binClassificationRes、binPositionRes
save(networkBioExp, binClassificationRes, binPositionRes, 
     file = paste("spatialDataProcess/stTargetData_tissue/",sampleID,"_annoNetworkBiomarker.RData",sep = ""))
# 汇总注释信息
expNetworkBiomarkerInfo = data.frame(cellType = as.character(), networkBiomarker = as.character(), cellNum = as.integer())
for(cellTypeI in names(binClassificationRes)){
  binClassificationResI = binClassificationRes[[cellTypeI]]
  for(networkBiomarkerI in names(binClassificationResI)){
    cellNum = length(binClassificationResI[[networkBiomarkerI]])
    expNetworkBiomarkerInfoI = data.frame(cellType = cellTypeI, networkBiomarker = networkBiomarkerI, cellNum = cellNum)
    expNetworkBiomarkerInfo = rbind(expNetworkBiomarkerInfo, expNetworkBiomarkerInfoI)
  }
}
write.csv(expNetworkBiomarkerInfo, file = paste("spatialDataProcess/expNetworkBiomarkerInfo_", sampleID, "_avg.csv", sep = ""))
# expNetworkBiomarkerInfo = expNetworkBiomarkerInfo[expNetworkBiomarkerInfo$cellNum > 1,]
# 统计所有networkBiomarker之间的关系-----------------------------------------
# 加载clusterTOcluster
clusterTOcluster = read.csv(file = "spatialDataProcess/genePairStatistic/clusterTOcluster(threshold=10_0.7).csv",
                            check.names = FALSE, row.names = 1)
rownames(clusterTOcluster) = networkBiomarkerSummary$networkBiomarker
colnames(clusterTOcluster) = networkBiomarkerSummary$networkBiomarker
# 从clusterTOcluster中整理networkBiomarker相互作用结果
iA.networkBiomarker = list()
noIA.networkBiomarker = list()
for(i in 1:nrow(clusterTOcluster)){
  networkBiomarkerI = rownames(clusterTOcluster)[i]
  for(j in 1:ncol(clusterTOcluster)){
    networkBiomarkerJ = colnames(clusterTOcluster)[j]
    if(clusterTOcluster[i,j] == 0){
      noIA.networkBiomarker = c(noIA.networkBiomarker, list(c(networkBiomarkerI, networkBiomarkerJ)))
    }else if(clusterTOcluster[i,j] == 1){
      iA.networkBiomarker = c(iA.networkBiomarker, list(c(networkBiomarkerI, networkBiomarkerJ)))
    }
  }
}
iA.networkBiomarker = lapply(iA.networkBiomarker, function(x){x=x[order(x)]})
iA.networkBiomarker = unique(iA.networkBiomarker)
noIA.networkBiomarker = lapply(noIA.networkBiomarker, function(x){x=x[order(x)]})
noIA.networkBiomarker = unique(noIA.networkBiomarker)
# 给定cellType1和cellType2的networkBiomarker，找出强因果关系和没有强因果关系的pairs--------------------
networkBiomarkerPairsInfo = data.frame(cellType1 = as.character(), networkBiomarker1 = as.character(),
                                       cellType2 = as.character(), networkBiomarker2 = as.character(),
                                       strength = as.integer(), distance = as.double())
for(i in 1:(length(unique(expNetworkBiomarkerInfo$cellType))-1)){
  for(j in (i+1):length(unique(expNetworkBiomarkerInfo$cellType))){
    cellType1 = unique(expNetworkBiomarkerInfo$cellType)[i]
    cellType2 = unique(expNetworkBiomarkerInfo$cellType)[j]
    networkBiomarkerSet1 = expNetworkBiomarkerInfo$networkBiomarker[expNetworkBiomarkerInfo$cellType == cellType1]
    networkBiomarkerSet2 = expNetworkBiomarkerInfo$networkBiomarker[expNetworkBiomarkerInfo$cellType == cellType2]
    for(networkBiomarker1 in networkBiomarkerSet1){
      for(networkBiomarker2 in networkBiomarkerSet2){
        strength = clusterTOcluster[networkBiomarker1, networkBiomarker2]
        binPositionRes1 = binPositionRes[[cellType1]][[networkBiomarker1]]
        binPositionRes2 = binPositionRes[[cellType2]][[networkBiomarker2]]
        distance = calDistance(binPositionRes1, binPositionRes2)
        networkBiomarkerPairsInfoI = data.frame(cellType1 = cellType1, networkBiomarker1 = networkBiomarker1,
                                                cellType2 = cellType2, networkBiomarker2 = networkBiomarker2,
                                                strength = strength, distance = distance)
        networkBiomarkerPairsInfo = rbind(networkBiomarkerPairsInfo, networkBiomarkerPairsInfoI)
      }
    }
  }
}
write.csv(networkBiomarkerPairsInfo, file = paste("spatialDataProcess/networkBiomarkerPairsInfo_", sampleID, "_avg.csv", sep = ""))
# 计算距离----------------------------------------------------
# 给定两个binPositionRes
# 参考配受体分析文献提供的距离计算方法
# 以细胞数少的cluster作为标准：nrow(binPositionRes1) <= nrow(binPositionRes2)
calDistance = function(binPositionRes1, binPositionRes2){
  distanceI = c()
  if(nrow(binPositionRes1) <= nrow(binPositionRes2)){
    for(bin1 in rownames(binPositionRes1)){
      distI = Inf
      for(bin2 in rownames(binPositionRes2)){
        if(!(bin1 %in% bin2)){
          currDist = dist(rbind(binPositionRes1[bin1,], binPositionRes2[bin2,]), method = "euclidean")
          distI = min(distI, currDist)
        }
      }
      if(distI != Inf){
        distanceI = c(distanceI, distI)
      }
    }
  }else{
    for(bin2 in rownames(binPositionRes2)){
      distI = Inf
      for(bin1 in rownames(binPositionRes1)){
        if(!(bin1 %in% bin2)){
          currDist = dist(rbind(binPositionRes2[bin2,], binPositionRes1[bin1,]), method = "euclidean")
          distI = min(distI, currDist)
        }
      }
      if(distI != Inf){
        distanceI = c(distanceI, distI)
      }
    }
  }
  distance = mean(distanceI)
  return(distance)
}
# 平均距离法
calDistance = function(binPositionRes1, binPositionRes2){
  distance = c()
  for(i in 1:nrow(binPositionRes1)){
    for(j in 1:nrow(binPositionRes2)){
      distanceI = dist(rbind(as.vector(binPositionRes1[i,]), as.vector(binPositionRes2[j,])), method = "euclidean")
      distance = c(distance, as.vector(distanceI))
    }
  }
  return(mean(distance))
}
# 质心法
calculate_centroid = function(cluster){
  return(colMeans(cluster))
}
calDistance = function(binPositionRes1, binPositionRes2){
  centroid1 = calculate_centroid(binPositionRes1)
  centroid2 = calculate_centroid(binPositionRes2)
  distance = dist(rbind(centroid1, centroid2), method = "euclidean")
  return(as.numeric(distance))
}


binPositionRes1 = binPositionRes[[networkBiomarker1]]
binPositionRes2 = binPositionRes[[networkBiomarker2]]
group0I = c()
if(nrow(binPositionRes1) <= nrow(binPositionRes2)){
  for(bin1 in rownames(binPositionRes1)){
    distI = Inf
    for(bin2 in rownames(binPositionRes2)){
      if(!(bin1 %in% bin2)){
        currDist = dist(rbind(binPositionRes1[bin1,], binPositionRes2[bin2,]), method = "euclidean")
        distI = min(distI, currDist)
      }
    }
    if(distI != Inf){
      group0I = c(group0I, distI)
    }
  }
}else{
  for(bin2 in rownames(binPositionRes2)){
    distI = Inf
    for(bin1 in rownames(binPositionRes1)){
      if(!(bin1 %in% bin2)){
        currDist = dist(rbind(binPositionRes2[bin2,], binPositionRes1[bin1,]), method = "euclidean")
        distI = min(distI, currDist)
      }
    }
    if(distI != Inf){
      group0I = c(group0I, distI)
    }
  }
}


