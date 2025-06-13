
# 函数设置
library(stringr)
library(ggplot2)
library(Seurat)
library(monocle3)

# getNetSetNet01
getNetSetNet01 = function(cellNetwork01){
  # 先去除重复的数据行
  cellNetwork01$netSetRes = 0
  for(i in 1:nrow(cellNetwork01)){
    netSet = names(cellNetwork01)[1:ncol(cellNetwork01)-1][which(cellNetwork01[i,1:ncol(cellNetwork01)-1] == 1)]
    # 将网络合并
    netSet = unlist(strsplit(netSet, split = "_"))
    #netSet = netSet[!duplicated(netSet)]
    netSet = paste(netSet, collapse = "_")
    cellNetwork01$netSetRes[i] = netSet
  }
  # 删除netSet相同的行
  cellNetwork01 = cellNetwork01[!duplicated(cellNetwork01$netSetRes),]
  # 将行名改为netSet
  row.names(cellNetwork01) = cellNetwork01$netSetRes
  # 删除netSet列
  netSetNet01 = cellNetwork01[,-ncol(cellNetwork01)]
  return(netSetNet01)
}

# 统计cellNetwork01P中每个细胞的netSet
getNetSetRes = function(cellNetwork01){
  cellNetwork01$netSetRes = 0
  for(i in 1:nrow(cellNetwork01)){
    netSet = names(cellNetwork01)[1:ncol(cellNetwork01)-1][which(cellNetwork01[i,1:ncol(cellNetwork01)-1] == 1)]
    # 将网络合并
    netSet = unlist(strsplit(netSet, split = "_"))
    #netSet = netSet[!duplicated(netSet)]
    netSet = paste(netSet, collapse = "_")
    cellNetwork01$netSetRes[i] = netSet
  }
  cellNetwork01_netSet = cellNetwork01
  return(cellNetwork01_netSet)
}

# 计算单个cluster中networked marker的ratio，以便于进一步判断cluster是否合格
getNetworkMarkerRatio_cluster = function(cellNetwork01_cluster){
  
  # cellNetwork01_cluster，某个cluster的cellNetwork01结果
  
  # 开始计算networkResTableCluster结果
  # 建表
  resTable = as.data.frame(matrix(ncol = 1, nrow = length(names(cellNetwork01_cluster))))
  colnames(resTable) = "networkMarkerRatio"
  # 行名命名
  rownames(resTable) = names(cellNetwork01_cluster)
  resTable[is.na(resTable)] = 0
  
  
  rowFlag = 1
  for (i in 1:length(names(cellNetwork01_cluster))){
    
    # 提取表达该networkMarker的细胞, length(cellId)代表在该cluster中共表达该networkMarker的细胞数
    cellId = row.names(cellNetwork01_cluster)[which(cellNetwork01_cluster[,i] == 1)]
    
    # 若没有细胞表达该networkMarker
    if(length(cellId) == 0){
      rowFlag = rowFlag + 1
      next
    }else{
      resTable[rowFlag, 1] = length(cellId)
    }
    
    rowFlag = rowFlag + 1
  }
  
  #将结果转换成百分比
  for(i in 1:dim(resTable)[2]){
    resTable[,i] = round(resTable[,i]/nrow(cellNetwork01_cluster)*100, 2)
  }
  
  return(resTable)
}

# 输入cluster的resTable结果，通过计算高表达networkMarker的数量来判断cluster的有效性
getEffectiveCluster = function(resTable){
  # resTable中的列名是networkMarkerRatio
  # 高表达networkMarker是指networkMarkerRatio不低于highNetworkMarkerRatioLowest的networkMarker
  highNetworkMarkerCount = length(which(resTable$networkMarkerRatio >= highNetworkMarkerRatioLowest))
  # 有效cluster要求高表达networkMarker数目不低于highNetworkMarkerCountLowest
  if(highNetworkMarkerCount >= highNetworkMarkerCountLowest){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

# 计算networkMarker在各个cluster中的ratio值
getNetworkMarkerRatio = function(cellNetwork01, cellClusterRes){
  
  # cellClusterRes
  # cellNetwork01
  
  # 统计每个cluster的数量
  staticForCluster = as.data.frame(table(cellClusterRes$clusterRes)) 
  staticForCluster$Var1 = as.character(staticForCluster$Var1)
  
  # 开始计算networkResTableCluster结果
  # 建表
  resTable = as.data.frame(matrix(ncol = length(unique(cellClusterRes$clusterRes)), nrow = length(names(cellNetwork01))))
  colnames(resTable) = staticForCluster$Var1
  # 行名命名
  rownames(resTable) = names(cellNetwork01)
  resTable[is.na(resTable)] = 0
  
  
  rowFlag = 1
  for (i in 1:length(names(cellNetwork01))){
    
    # 提取表达该networkMarker的细胞
    cellId = row.names(cellNetwork01)[which(cellNetwork01[,i] == 1)]
    
    # 若没有细胞表达该networkMarker
    if(length(cellId) == 0){
      rowFlag = rowFlag + 1
      next
    }
    
    # 统计这些细胞的cluster
    # 根据细胞id提取cluster
    tmpcellClusterRes = cellClusterRes[row.names(cellClusterRes) %in% cellId,]
    # 统计细胞id对应的细胞类型
    staticsFortmpCluster = as.data.frame(table(tmpcellClusterRes$clusterRes))
    staticsFortmpCluster$Var1 = as.character(staticsFortmpCluster$Var1)
    
    for(j in 1:dim(staticsFortmpCluster)[1]){
      resTable[rowFlag, staticsFortmpCluster$Var1[j]] = staticsFortmpCluster$Freq[j]
    }
    rowFlag = rowFlag + 1
  }
  
  #将结果转换成百分比
  for(i in 1:dim(resTable)[2]){
    resTable[,i] = round(resTable[,i]/staticForCluster$Freq[i]*100, 2)
  }
  
  return(resTable)
}

# 统计每个cluster对应的highNetworkMarker
# 处理聚类结果，从networkResTableCluster中提取每个cluster的高表达networkMarker，最后作为该cluster的bioMarker
getHighNetworkMarker = function(networkResTableCluster){
  res = as.data.frame(matrix(nrow = nrow(networkResTableCluster),ncol = ncol(networkResTableCluster)))
  colnames(res) = colnames(networkResTableCluster)
  maxRow = 0
  
  #一些行的最大值相同，提取出来看看怎么回事
  confusedRow = c()
  
  for(i in 1:ncol(networkResTableCluster)){
    networkResTableCluster = networkResTableCluster[order(networkResTableCluster[,i],decreasing = T),]
    highNetworkMarkerID = which(networkResTableCluster[,i] >= highNetworkMarkerRatioLowest)
    if(length(highNetworkMarkerID) == 0){
      next
    }
    res[1:length(highNetworkMarkerID),i] = rownames(networkResTableCluster)[highNetworkMarkerID]
    maxRow = max(maxRow, length(highNetworkMarkerID))
  }
  res[is.na(res)] = ""
  res = res[1:maxRow,]
  highNetworkMarker = as.data.frame(res)
  colnames(highNetworkMarker) = colnames(networkResTableCluster)
  return(highNetworkMarker)
}

# 对聚类结果进行统计，包括clusterRes、cellCount_cluster、genePairs、cellCount_coexpress、Ratio
getResultStatistic = function(networkResTableCluster, cellClusterRes){
  clusterRes = c()
  genePairs = c()
  Ratio = c()
  for(i in c(1:max(as.numeric(names(networkResTableCluster))))){
    networkResTableCluster = networkResTableCluster[order(networkResTableCluster[,as.character(i)], decreasing = T),]
    targetGenePairs = row.names(networkResTableCluster)[which(networkResTableCluster[,as.character(i)] >= highNetworkMarkerRatioLowest)]
    targetRatio = networkResTableCluster[,as.character(i)][which(networkResTableCluster[,as.character(i)] >= highNetworkMarkerRatioLowest)]
    targetClusterRes = rep(i, length(targetGenePairs))
    clusterRes = c(clusterRes, targetClusterRes)
    genePairs = c(genePairs, targetGenePairs)
    Ratio = c(Ratio, targetRatio)
  }
  resultStatistic = data.frame(clusterRes = clusterRes, genePairs = genePairs, Ratio = Ratio)
  # 在resultStatistic中添加cellCount_cluster和cellCount_express两列，这两项信息主要从cellClusterRes中提取
  cellCount_cluster = c()
  cellCount_express = c()
  for(i in unique(resultStatistic$clusterRes)){
    cellClusterRes_cluster = cellClusterRes[which(cellClusterRes$clusterRes == i),]
    targetCellCount_cluster = rep(nrow(cellClusterRes_cluster), length(which(resultStatistic$clusterRes == i)))
    targetGenePairs = resultStatistic$genePairs[which(resultStatistic$clusterRes == i)]
    targetCellCount_express = c()
    for(j in targetGenePairs){
      targetCellCount_express = c(targetCellCount_express, length(which(cellClusterRes_cluster[,j] == 1)))
    }
    cellCount_cluster = c(cellCount_cluster, targetCellCount_cluster)
    cellCount_express = c(cellCount_express, targetCellCount_express)
  }
  resultStatistic$cellCount_cluster = cellCount_cluster
  resultStatistic$cellCount_express = cellCount_express
  
  return(resultStatistic)
}

# 计算networkBiomarker在各个cluster中的ratio值
getNetworkBiomarkerRatio = function(networkMarkerSet, cellClusterRes){
  
  # cellClusterRes
  # cellNetwork01
  
  # 统计每个cluster的数量
  staticForCluster = as.data.frame(table(cellClusterRes$clusterRes)) 
  staticForCluster$Var1 = as.character(staticForCluster$Var1)
  
  # 开始计算networkResTableCluster结果
  # 建表
  resTable = as.data.frame(matrix(ncol = length(unique(cellClusterRes$clusterRes)), nrow = 1))
  colnames(resTable) = staticForCluster$Var1
  # 合并networkMarkerSet中的基因
  networkMarkerSetGene = unlist(strsplit(networkMarkerSet, split = "_"))
  networkMarkerSetGene = networkMarkerSetGene[!duplicated(networkMarkerSetGene)]
  networkMarkerSetGene = networkMarkerSetGene[order(networkMarkerSetGene, decreasing = F)]
  networkBiomarker = paste(networkMarkerSetGene, collapse = "_")
  # 行名命名
  rownames(resTable) = networkBiomarker
  resTable[is.na(resTable)] = 0
  
  
  # 提取表达该networkBiomarker的细胞
  if(length(networkMarkerSet) == 1){
    cellId = row.names(cellClusterRes)[which(cellClusterRes[,networkMarkerSet] == length(networkMarkerSet))]
  }else{
    cellId = row.names(cellClusterRes)[which(rowSums(cellClusterRes[,networkMarkerSet]) == length(networkMarkerSet))]
  }
  
  # 若没有细胞表达该networkBiomarker
  if(length(cellId) == 0){
    return(resTable)
  }
  
  # 统计这些细胞的cluster
  # 根据细胞id提取cluster
  tmpcellClusterRes = cellClusterRes[row.names(cellClusterRes) %in% cellId,]
  # 统计细胞id对应的细胞类型
  staticsFortmpCluster = as.data.frame(table(tmpcellClusterRes$clusterRes))
  staticsFortmpCluster$Var1 = as.character(staticsFortmpCluster$Var1)
  
  for(j in 1:dim(staticsFortmpCluster)[1]){
    resTable[1, staticsFortmpCluster$Var1[j]] = staticsFortmpCluster$Freq[j]
  }
  
  #将结果转换成百分比
  for(i in 1:dim(resTable)[2]){
    resTable[,i] = round(resTable[,i]/staticForCluster$Freq[i]*100, 2)
  }
  
  return(resTable)
}


# PBMC
cellTypeSet = c("memory CD4+ T cells", "memory CD8+ T cells", "naive CD4+ T cells",
                "naive CD8+ T cells", "NKT", "MAIT", "Tregs",
                "NK",
                "Plasma cells", "pre-B cells",
                "cDC", "pDC",
                "Basophils", "Neutrophils")

for(i in cellTypeSet){
  # 准备工作
  # 载入输入数据和变量设置
  cellType = i
  cellNetwork01 = get(load(paste("clusterDivision/pbmcClusterData/cellNetwork01P/", str_replace_all(cellType,pattern = " ",replacement = "_"), ".RData", sep = "")))
  allCellCount = nrow(cellNetwork01)
  highNetworkMarkerRatioLowest = 80
  cellCountRatioLowest = 0.8
  highNetworkMarkerCountLowest = 1
  
  # 输入数据处理
  netSetNet01 = getNetSetNet01(cellNetwork01)
  
  # 获得cell对应的netSet，方便后面利用netSet获取cell
  cellNetwork01_netSet = getNetSetRes(cellNetwork01)
  
  # 定义一个列表，该列表中用于存储聚类出来的cluster，它们等待被判断是否有效
  # clusterRes_ls初始化，首先将整个处理后的输入数据看作一个大类
  # clusterRes_ls中cluster的形式是netSetNet01形式
  clusterRes_ls = list(netSetNet01)
  
  # 定义一个列表，该列表存储满足条件的cluster
  # 列表中的cluster形式为cell——network+clusterRes
  cellClusterRes_ls = list()
  netSetClusterRes_ls = list()
  # 记录已有满足条件的cluster覆盖的细胞数
  satisfiedCellCount = 0
  
  # 在循环迭代过程中会有新的聚类出来的cluster被赋值给clusterRes_ls
  while(TRUE){
    cluster = clusterRes_ls[[1]]
    # 判断当前cluster是否满足条件，getNetworkMarkerRatio_cluster、getEffectiveCluster()
    # 获得cellNetworkMarker01_cluster
    cellName = row.names(cellNetwork01_netSet)[which(cellNetwork01_netSet$netSetRes %in% row.names(cluster))]
    if(length(cellName) <= allCellCount*0.01){
      print("此时cluster细胞数过少")
      break
    }
    cellNetwork01_cluster = cellNetwork01[cellName,]
    networkResTable_cluster = getNetworkMarkerRatio_cluster(cellNetwork01_cluster)
    if (getEffectiveCluster(networkResTable_cluster)) {
      print(paste0("有效",nrow(cluster)))
      # 如果满足条件，将当前cluster添加到cellClusterRes_ls中
      cellNetwork01_cluster$clusterRes = length(cellClusterRes_ls)+1
      cellClusterRes_ls = append(cellClusterRes_ls, list(cellNetwork01_cluster))
      cluster$clusterRes = length(netSetClusterRes_ls)+1
      netSetClusterRes_ls = append(netSetClusterRes_ls, list(cluster))
      satisfiedCellCount = satisfiedCellCount + nrow(cellNetwork01_cluster)
      clusterRes_ls = clusterRes_ls[-1]
      # 设置迭代终止条件
      # 计算已有满足条件的cluster中的细胞数占总细胞数的比例
      cellCountRatio = satisfiedCellCount/allCellCount
      
      # 如果满足条件的类别数目达到或超过指定比例，停止迭代
      if (cellCountRatio >= cellCountRatioLowest){
        print("细胞数占比满足条件而终止")
        break
      }
    }else{
      print(paste0("无效",nrow(cluster)))
      # 如果不满足条件，将当前cluster划分为两个cluster，并将新得到的cluster添加到newClusterRes_ls中
      # 再聚类
      set.seed(1000)
      newClusterRes = kmeans(cluster, centers = 2, nstart = 1, iter.max = 10)
      newClusterRes_df = as.data.frame(newClusterRes$cluster)
      names(newClusterRes_df) = c("clusterRes")
      #table(newClusterRes_df$clusterRes)
      newClusterRes_df$clusterRes = as.numeric(newClusterRes_df$clusterRes)
      if(length(unique(newClusterRes_df$clusterRes)) == 2){
        netSet1 = row.names(newClusterRes_df)[which(newClusterRes_df$clusterRes == 1)]
        netSet2 = row.names(newClusterRes_df)[which(newClusterRes_df$clusterRes == 2)]
        cluster1 = netSetNet01[netSet1,]
        cluster2 = netSetNet01[netSet2,]
        print(paste0("再聚类",nrow(cluster1)))
        print(paste0("再聚类",nrow(cluster2)))
        #newClusterRes_ls = append(newClusterRes_ls, list(cluster1))
        #newClusterRes_ls = append(newClusterRes_ls, list(cluster2))
        clusterRes_ls = append(clusterRes_ls, list(cluster1))
        clusterRes_ls = append(clusterRes_ls, list(cluster2))
        clusterRes_ls = clusterRes_ls[-1]
      }else{
        clusterRes_ls = clusterRes_ls[-1]
        clusterRes_ls = clusterRes_ls[order(sapply(clusterRes_ls, nrow), decreasing = TRUE)]
        next
      }
      
      
    }
    clusterRes_ls = clusterRes_ls[order(sapply(clusterRes_ls, nrow), decreasing = TRUE)]
    
  }
  
  # 迭代终止，得到cellClusterRes_ls
  # 合并cellClusterRes_ls，得到cellClusterRes (df)
  cellClusterRes = do.call(rbind, cellClusterRes_ls)
  # 没有被聚类到的细胞被分配到cluster0
  otherCellClusterRes = cellNetwork01[-which(row.names(cellNetwork01) %in% row.names(cellClusterRes)),]
  if(nrow(otherCellClusterRes) > 0){
    otherCellClusterRes$clusterRes = 0
    cellClusterRes = rbind(cellClusterRes, otherCellClusterRes)
  }
  # 调整cellClusterRes和cellNetwork01行顺序
  cellClusterRes = cellClusterRes[match(row.names(cellNetwork01), row.names(cellClusterRes)),]
  write.csv(cellClusterRes, paste("clusterDivision/pbmcClusterRes/networkClusterResPv2/cellClusterRes_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""))
  
  
  # 合并netSetClusterRes_ls，得到netSetClusterRes (df)
  netSetClusterRes = do.call(rbind, netSetClusterRes_ls)
  # 没有被聚类到的netSet被分配到cluster0
  otherNetSetClusterRes = netSetNet01[-which(row.names(netSetNet01) %in% row.names(netSetClusterRes)),]
  if(nrow(otherNetSetClusterRes) > 0){
    otherNetSetClusterRes$clusterRes = 0
    netSetClusterRes = rbind(netSetClusterRes, otherNetSetClusterRes)
  }
  # 调整cellClusterRes和cellNetwork01行顺序
  netSetClusterRes = netSetClusterRes[match(row.names(netSetNet01), row.names(netSetClusterRes)),]
  write.csv(netSetClusterRes, paste("clusterDivision/pbmcClusterRes/networkClusterResPv2/netSetClusterRes_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""))
  
  
  # 计算networkMarker在每个cluster中的ratio
  networkResTableCluster = getNetworkMarkerRatio(cellNetwork01, cellClusterRes)
  write.csv(networkResTableCluster, paste("clusterDivision/resTableClusterP/networkResTableClusterPv2/networkResTableCluster_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""))
  
  
  highNetworkMarker = getHighNetworkMarker(networkResTableCluster)
  write.csv(highNetworkMarker, paste("clusterDivision/highNetworkMarkerTablePv2/highNetworkMarkerTable_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""))
  
  # 对聚类结果进行统计，包括clusterRes、cellCount_cluster、genePairs、cellCount_coexpress、Ratio
  resultStatistic = getResultStatistic(networkResTableCluster, cellClusterRes)
  write.csv(resultStatistic, paste("clusterDivision/resultStatisticPv2/resultStatistic_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""))
  print(cellType)
}




# LUAD
cellTypeSet = c("cDC", "Monocytes", "Tregs")

for(i in cellTypeSet){
  # 准备工作
  # 载入输入数据和变量设置
  cellType = i
  cellNetwork01 = get(load(paste("clusterDivision/luadClusterData/cellNetwork01L/", str_replace_all(cellType,pattern = " ",replacement = "_"), ".RData", sep = "")))
  allCellCount = nrow(cellNetwork01)
  highNetworkMarkerRatioLowest = 80
  cellCountRatioLowest = 0.8
  highNetworkMarkerCountLowest = 1
  
  # 输入数据处理
  netSetNet01 = getNetSetNet01(cellNetwork01)
  
  # 获得cell对应的netSet，方便后面利用netSet获取cell
  cellNetwork01_netSet = getNetSetRes(cellNetwork01)
  
  # 定义一个列表，该列表中用于存储聚类出来的cluster，它们等待被判断是否有效
  # clusterRes_ls初始化，首先将整个处理后的输入数据看作一个大类
  # clusterRes_ls中cluster的形式是netSetNet01形式
  clusterRes_ls = list(netSetNet01)
  
  # 定义一个列表，该列表存储满足条件的cluster
  # 列表中的cluster形式为cell——network+clusterRes
  cellClusterRes_ls = list()
  netSetClusterRes_ls = list()
  # 记录已有满足条件的cluster覆盖的细胞数
  satisfiedCellCount = 0
  
  # 在循环迭代过程中会有新的聚类出来的cluster被赋值给clusterRes_ls
  while(TRUE){
    cluster = clusterRes_ls[[1]]
    if(nrow(cluster) <= 2){
      break
    }
    # 判断当前cluster是否满足条件，getNetworkMarkerRatio_cluster、getEffectiveCluster()
    # 获得cellNetworkMarker01_cluster
    cellName = row.names(cellNetwork01_netSet)[which(cellNetwork01_netSet$netSetRes %in% row.names(cluster))]
    cellNetwork01_cluster = cellNetwork01[cellName,]
    networkResTable_cluster = getNetworkMarkerRatio_cluster(cellNetwork01_cluster)
    if (getEffectiveCluster(networkResTable_cluster)) {
      print(paste0("有效",nrow(cluster)))
      # 如果满足条件，将当前cluster添加到cellClusterRes_ls中
      cellNetwork01_cluster$clusterRes = length(cellClusterRes_ls)+1
      cellClusterRes_ls = append(cellClusterRes_ls, list(cellNetwork01_cluster))
      cluster$clusterRes = length(netSetClusterRes_ls)+1
      netSetClusterRes_ls = append(netSetClusterRes_ls, list(cluster))
      satisfiedCellCount = satisfiedCellCount + nrow(cellNetwork01_cluster)
      clusterRes_ls = clusterRes_ls[-1]
      # 设置迭代终止条件
      # 计算已有满足条件的cluster中的细胞数占总细胞数的比例
      cellCountRatio = satisfiedCellCount/allCellCount
      
      # 如果满足条件的类别数目达到或超过指定比例，停止迭代
      if (cellCountRatio >= cellCountRatioLowest){
        print("细胞数占比满足条件而终止")
        break
      }
    }else{
      print(paste0("无效",nrow(cluster)))
      # 如果不满足条件，将当前cluster划分为两个cluster，并将新得到的cluster添加到newClusterRes_ls中
      # 再聚类
      set.seed(1000)
      newClusterRes = kmeans(cluster, centers = 2, nstart = 1, iter.max = 10)
      newClusterRes_df = as.data.frame(newClusterRes$cluster)
      names(newClusterRes_df) = c("clusterRes")
      #table(newClusterRes_df$clusterRes)
      newClusterRes_df$clusterRes = as.numeric(newClusterRes_df$clusterRes)
      if(length(unique(newClusterRes_df$clusterRes)) == 2){
        netSet1 = row.names(newClusterRes_df)[which(newClusterRes_df$clusterRes == 1)]
        netSet2 = row.names(newClusterRes_df)[which(newClusterRes_df$clusterRes == 2)]
        cluster1 = netSetNet01[netSet1,]
        cluster2 = netSetNet01[netSet2,]
        print(paste0("再聚类",nrow(cluster1)))
        print(paste0("再聚类",nrow(cluster2)))
        #newClusterRes_ls = append(newClusterRes_ls, list(cluster1))
        #newClusterRes_ls = append(newClusterRes_ls, list(cluster2))
        clusterRes_ls = append(clusterRes_ls, list(cluster1))
        clusterRes_ls = append(clusterRes_ls, list(cluster2))
        clusterRes_ls = clusterRes_ls[-1]
      }else{
        clusterRes_ls = clusterRes_ls[-1]
        clusterRes_ls = clusterRes_ls[order(sapply(clusterRes_ls, nrow), decreasing = TRUE)]
        next
      }
      
      
    }
    clusterRes_ls = clusterRes_ls[order(sapply(clusterRes_ls, nrow), decreasing = TRUE)]
    
  }
  
  # 迭代终止，得到cellClusterRes_ls
  # 合并cellClusterRes_ls，得到cellClusterRes (df)
  cellClusterRes = do.call(rbind, cellClusterRes_ls)
  # 没有被聚类到的细胞被分配到cluster0
  otherCellClusterRes = cellNetwork01[-which(row.names(cellNetwork01) %in% row.names(cellClusterRes)),]
  if(nrow(otherCellClusterRes) > 0){
    otherCellClusterRes$clusterRes = 0
    cellClusterRes = rbind(cellClusterRes, otherCellClusterRes)
  }
  # 调整cellClusterRes和cellNetwork01行顺序
  cellClusterRes = cellClusterRes[match(row.names(cellNetwork01), row.names(cellClusterRes)),]
  write.csv(cellClusterRes, paste("clusterDivision/luadClusterRes/networkClusterResL0730/cellClusterRes_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""))
  
  
  # 合并netSetClusterRes_ls，得到netSetClusterRes (df)
  netSetClusterRes = do.call(rbind, netSetClusterRes_ls)
  # 没有被聚类到的netSet被分配到cluster0
  otherNetSetClusterRes = netSetNet01[-which(row.names(netSetNet01) %in% row.names(netSetClusterRes)),]
  if(nrow(otherNetSetClusterRes) > 0){
    otherNetSetClusterRes$clusterRes = 0
    netSetClusterRes = rbind(netSetClusterRes, otherNetSetClusterRes)
  }
  # 调整cellClusterRes和cellNetwork01行顺序
  netSetClusterRes = netSetClusterRes[match(row.names(netSetNet01), row.names(netSetClusterRes)),]
  write.csv(netSetClusterRes, paste("clusterDivision/luadClusterRes/networkClusterResL0730/netSetClusterRes_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""))
  
  
  # 计算networkMarker在每个cluster中的ratio
  networkResTableCluster = getNetworkMarkerRatio(cellNetwork01, cellClusterRes)
  write.csv(networkResTableCluster, paste("clusterDivision/resTableClusterL/networkResTableClusterL0730/networkResTableCluster_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""))
  
  
  highNetworkMarker = getHighNetworkMarker(networkResTableCluster)
  write.csv(highNetworkMarker, paste("clusterDivision/highNetworkMarkerTableLv2/highNetworkMarkerTable_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""))
  
  # 对聚类结果进行统计，包括clusterRes、cellCount_cluster、genePairs、cellCount_coexpress、Ratio
  resultStatistic = getResultStatistic(networkResTableCluster, cellClusterRes)
  write.csv(resultStatistic, paste("clusterDivision/resultStatisticLv2/resultStatistic_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""))
  print(cellType)
}

# 利用findGoodNetworkMarker对resTableCluster进行处理
findGoodNetworkMarker = function(resTableCluster, nsd){
  goodNetworkMarker = c()
  for(i in 1:nrow(resTableCluster)){
    meantmp = mean(as.numeric(resTableCluster[i,]))
    sdtmp = sd(as.numeric(resTableCluster[i,]))
    if(max(as.numeric(resTableCluster[i,]))>0 & max(as.numeric(resTableCluster[i,])) >= meantmp+nsd*sdtmp){
      goodNetworkMarker = c(goodNetworkMarker, i)
    }
  }
  return(resTableCluster[goodNetworkMarker,])
}

findGoodNetworkMarkerV2 = function(resTableCluster, threshold){
  goodRow=c()
  for(i in 1:nrow(resTableCluster)){
    max = order(as.numeric(resTableCluster[i,]), decreasing = T)[1]
    second = order(as.numeric(resTableCluster[i,]), decreasing = T)[2]
    if(max(as.numeric(resTableCluster[i,]))>0 & resTableCluster[i, second] <= resTableCluster[i, max]*threshold){
      goodRow=c(goodRow,i)
    }
  }
  return(resTableCluster[goodRow,])
}

cellTypeSet = c("memory CD4+ T cells","memory CD8+ T cells","naive CD4+ T cells",
                "naive CD8+ T cells","NKT","MAIT","Tregs",
                "NK","Plasma cells","pre-B cells","cDC","pDC","Monocytes",
                "Basophils","Neutrophils")
for(i in cellTypeSet){
  cellType = i
  resTableCluster = read.csv(paste("clusterDivision/resTableClusterP/networkResTableClusterPv2/networkResTableCluster_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""),
                             check.names = F, row.names = 1)
  if("0" %in% names(resTableCluster)){
    resTableCluster = resTableCluster[,-1]
  }
  if(ncol(resTableCluster) < 2){
    next
  }
  goodNetworkMarker = findGoodNetworkMarkerV2(resTableCluster, 1/2)
  if(length(goodNetworkMarker) == 0){
    next
  }
  write.csv(goodNetworkMarker, file = paste("clusterDivision/resTableClusterP/goodNetworkMarkerPv2/goodNetworkMarker_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""))
  print(cellType)
  
}


# 因为luad中一些细胞类型只有一个cluster，且network marker数量较多，不宜计算每个network marker的ratio，所以选择几个单独计算
# cDC
specifyNetworkMarker = c("CD80_CD86", "CD74_CD80_CD86",
                         "CD80_CD86_HLA-DRA", "CD80_PLEK",
                         "CD80_CD86_PTPRC", "CD80_CD86_PLEK",
                         "CD80_CD86_HLA-DMA", "CD80_PLEK_PTPRC",
                         "HLA-DQA1_HLA-DRB1_NECTIN2", "CD80_LAPTM5")
cellNetwork01 = get(load(paste("clusterDivision/luadClusterData/cellNetwork01L/", str_replace_all(cellType,pattern = " ",replacement = "_"), ".RData", sep = "")))
cellNetwork01 = cellNetwork01[,specifyNetworkMarker]
cellClusterRes = cellNetwork01
cellClusterRes$clusterRes = 1
networkResTableCluster = getNetworkMarkerRatio(cellNetwork01, cellClusterRes)
write.csv(networkResTableCluster, paste("clusterDivision/resTableClusterL/networkResTableClusterLv2/networkResTableCluster_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""))

# pre-B cells
specifyNetworkMarker = c("BANK1_MS4A1", "MS4A1_PTPRC",
                         "CD79A_MS4A1", "BANK1_MS4A1_PTPRC",
                         "HLA-A_MS4A1_PTPRC", "BANK1_CD79A_MS4A1",
                         "HLA-E_MS4A1_PTPRC", "HLA-DQA1_MS4A1_PTPRC")
cellNetwork01 = get(load(paste("clusterDivision/luadClusterData/cellNetwork01L/", str_replace_all(cellType,pattern = " ",replacement = "_"), ".RData", sep = "")))
cellNetwork01 = cellNetwork01[,specifyNetworkMarker]
cellClusterRes = cellNetwork01
cellClusterRes$clusterRes = 1
networkResTableCluster = getNetworkMarkerRatio(cellNetwork01, cellClusterRes)
write.csv(networkResTableCluster, paste("clusterDivision/resTableClusterL/networkResTableClusterLv2/networkResTableCluster_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""))

# Macrophages
specifyNetworkMarker = c("C1QA_TYROBP", "C1QA_HLA-DRA_TYROBP",
                         "C1QA_CD74_TYROBP", "C1QA_FCER1G_TYROBP",
                         "C1QA_C1QB", "C1QA_LAPTM5_TYROBP",
                         "C1QA_C1QB_TYROBP", "C1QC_TYROBP",
                         "C1QC_HLA-DRA_TYROBP", "C1QC_CD74_TYROBP",
                         "C1QA_C1QC", "AIF1_C1QA",
                         "AIF1_C1QA_CD74", "C1QC_FCER1G",
                         "AIF1_C1QA_TYROBP", "AIF1_C1QA_HLA-DPA1",
                         "C1QA_C1QC_TYROBP", "C1QB_C1QC",
                         "C1QC_LAPTM5", "C1QC_FCER1G_TYROBP",
                         "AIF1_C1QA_HLA-DPB1", "C1QC_FCER1G_HLA-B",
                         "C1QC_HLA-DRB1_LAPTM5", "C1QC_FCER1G_HLA-A",
                         "C1QC_LAPTM5_TYROBP", "C1QC_HLA-DPB1_LAPTM5",
                         "C1QB_C1QC_TYROBP", "AIF1_C1QA_FCER1G",
                         "C1QA_C1QB_C1QC", "C1QA_C1QC_FCER1G")
cellNetwork01 = get(load(paste("clusterDivision/luadClusterData/cellNetwork01L/", str_replace_all(cellType,pattern = " ",replacement = "_"), ".RData", sep = "")))
cellNetwork01 = cellNetwork01[,specifyNetworkMarker]
cellClusterRes = cellNetwork01
cellClusterRes$clusterRes = 1
networkResTableCluster = getNetworkMarkerRatio(cellNetwork01, cellClusterRes)
write.csv(networkResTableCluster, paste("clusterDivision/resTableClusterL/networkResTableClusterLv2/networkResTableCluster_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""))














# PBMC 从networkResRableClusterV3中获得highMarker数据
cellTypeSet = c("pre-B cells",
                "cDC","pDC",
                "Monocytes",
                "Basophils","Neutrophils")
for(i in cellTypeSet){
  cellType = i
  cellNetwork01 = get(load(paste("clusterDivision/luadClusterData/cellNetwork01L/", str_replace_all(cellType,pattern = " ",replacement = "_"), ".RData", sep = "")))
  cellClusterRes = read.csv(paste("clusterDivision/luadClusterRes/networkClusterResLv2/cellClusterResv2/cellClusterRes_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""),
                            check.names = F, row.names = 1)
  # 计算networkMarker在每个cluster中的ratio
  networkResTableCluster = read.csv(paste("clusterDivision/resTableClusterL/networkResTableClusterLv3/networkResTableCluster_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""),
                                    check.names = F, row.names = 1)
  
  
  highNetworkMarker = getHighNetworkMarker(networkResTableCluster)
  if("0" %in% colnames(highNetworkMarker)){
    highNetworkMarker = highNetworkMarker[,-1]
  }
  write.csv(highNetworkMarker, paste("clusterDivision/highNetworkMarkerTablePv3/highNetworkMarkerTable_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""))
  networkBiomarkerRatio_ls = list()
  for(i in 1:ncol(highNetworkMarker)){
    networkMarkerSet = highNetworkMarker[,i]
    networkMarkerSet = networkMarkerSet[1]
    # 删除networkMarkerSet中的空格
    if(is.element("",networkMarkerSet)){
      networkMarkerSet = networkMarkerSet[-which(networkMarkerSet == "")] #排除空字符串
    }
    networkBiomarkerRatio_ls[[i]] = getNetworkBiomarkerRatio(networkMarkerSet, cellClusterRes)
  }
  networkBiomarkerRatio = do.call(rbind, networkBiomarkerRatio_ls)
  write.csv(networkBiomarkerRatio, paste("clusterDivision/networkBiomarkerRatioP/networkBiomarkerRatio_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""))
  
  
}


# LUAD 从networkResRableClusterV3中获得highMarker数据
cellTypeSet = c("pre-B cells",
                "cDC","pDC",
                "Monocytes",
                "Basophils","Neutrophils")
for(i in cellTypeSet){
  cellType = i
  cellNetwork01 = get(load(paste("clusterDivision/luadClusterData/cellNetwork01L/", str_replace_all(cellType,pattern = " ",replacement = "_"), ".RData", sep = "")))
  cellClusterRes = read.csv(paste("clusterDivision/luadClusterRes/networkClusterResLv2/cellClusterRes_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""),
                            check.names = F, row.names = 1)
  cellClusterRes = cellNetwork01
  cellClusterRes$clusterRes = 1 
  # 计算networkMarker在每个cluster中的ratio
  networkResTableCluster = read.csv(paste("clusterDivision/resTableClusterL/networkResTableClusterLv2/networkResTableCluster_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""),
                                    check.names = F, row.names = 1)
  
  
  highNetworkMarker = getHighNetworkMarker(networkResTableCluster)
  
  if("0" %in% colnames(highNetworkMarker)){
    colname = colnames(highNetworkMarker)[-1]
    highNetworkMarker = as.data.frame(highNetworkMarker[,-1])
    colnames(highNetworkMarker) = colname
  }
  write.csv(highNetworkMarker, paste("clusterDivision/highNetworkMarkerTableLv3/highNetworkMarkerTable_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""))
  networkBiomarkerRatio_ls = list()
  for(i in 1:ncol(highNetworkMarker)){
    networkMarkerSet = highNetworkMarker[,i]
    networkMarkerSet = networkMarkerSet[c(1:26)]
    # 删除networkMarkerSet中的空格
    if(is.element("",networkMarkerSet)){
      networkMarkerSet = networkMarkerSet[-which(networkMarkerSet == "")] #排除空字符串
    }
    networkBiomarkerRatio_ls[[i]] = getNetworkBiomarkerRatio(networkMarkerSet, cellClusterRes)
  }
  networkBiomarkerRatio = do.call(rbind, networkBiomarkerRatio_ls)
  write.csv(networkBiomarkerRatio, paste("clusterDivision/networkBiomarkerRatioL0730/networkBiomarkerRatio_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""))
  
  
}
