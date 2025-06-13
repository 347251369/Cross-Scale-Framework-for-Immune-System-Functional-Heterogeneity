
# 计算单个cluster中networked marker的ratio，以便于进一步判断cluster是否合格
getNetworkMarkerRatio_netSetCluster = function(netSetNet01_cluster){
  
  # cellNetwork01_cluster，某个cluster的cellNetwork01结果
  
  # 开始计算networkResTableCluster结果
  # 建表
  resTable = as.data.frame(matrix(ncol = 1, nrow = length(names(netSetNet01_cluster))))
  colnames(resTable) = "networkMarkerRatio"
  # 行名命名
  rownames(resTable) = names(netSetNet01_cluster)
  resTable[is.na(resTable)] = 0
  
  
  rowFlag = 1
  for (i in 1:length(names(netSetNet01_cluster))){
    
    # 提取表达该networkMarker的细胞, length(cellId)代表在该cluster中共表达该networkMarker的细胞数
    cellId = row.names(netSetNet01_cluster)[which(netSetNet01_cluster[,i] == 1)]
    
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
    resTable[,i] = round(resTable[,i]/nrow(netSetNet01_cluster)*100, 2)
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
getNetworkMarkerRatio_netSet = function(netSetNet01, netSetClusterRes){
  
  # netSetClusterRes
  # netSetNet01
  
  # 统计每个cluster的数量
  staticForCluster = as.data.frame(table(netSetClusterRes$clusterRes)) 
  staticForCluster$Var1 = as.character(staticForCluster$Var1)
  
  # 开始计算networkResTableCluster结果
  # 建表
  resTable = as.data.frame(matrix(ncol = length(unique(netSetClusterRes$clusterRes)), nrow = length(names(netSetNet01))))
  colnames(resTable) = staticForCluster$Var1
  # 行名命名
  rownames(resTable) = names(netSetNet01)
  resTable[is.na(resTable)] = 0
  
  
  rowFlag = 1
  for (i in 1:length(names(netSetNet01))){
    
    # 提取表达该networkMarker的细胞
    netSetId = row.names(netSetNet01)[which(netSetNet01[,i] == 1)]
    
    # 若没有细胞表达该networkMarker
    if(length(netSetId) == 0){
      rowFlag = rowFlag + 1
      next
    }
    
    # 统计这些细胞的cluster
    # 根据细胞id提取cluster
    tmpnetSetClusterRes = netSetClusterRes[row.names(netSetClusterRes) %in% netSetId,]
    # 统计细胞id对应的细胞类型
    staticsFortmpCluster = as.data.frame(table(tmpnetSetClusterRes$clusterRes))
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

# 通过netSetClusterRes直接更新cellClusterRes
getCellClusterRes = function(cellNetwork01, netSetClusterRes){
  # 在cellNetwork01中添加clusterRes列
  cellNetwork01$clusterRes = 0
  
  # 利用cellNetwork01的netSet列和netSetClusterRes确定clusterRes
  cellNetwork01$netSet = 0
  for(i in 1:nrow(cellNetwork01)){
    netSet = names(cellNetwork01)[1:ncol(cellNetwork01)-1][which(cellNetwork01[i,1:ncol(cellNetwork01)-1] == 1)]
    # 将网络合并
    netSet = unlist(strsplit(netSet, split = "_"))
    #netSet = netSet[!duplicated(netSet)]
    netSet = paste(netSet, collapse = "_")
    cellNetwork01$netSet[i] = netSet
  }
  
  for(i in 1:nrow(cellNetwork01)){
    if(length(which(rownames(netSetClusterRes) %in% cellNetwork01$netSet[i])) > 0)(
      cellNetwork01$clusterRes[i] = netSetClusterRes[which(rownames(netSetClusterRes) %in% cellNetwork01$netSet[i]),]$clusterRes
    )
  }
  
  cellClusterRes = cellNetwork01[,-ncol(cellNetwork01)]
  return(cellClusterRes)
}

for(i in cellTypeSet){
  # 准备工作
  # 载入输入数据和变量设置
  cellType = i
  cellNetwork01 = get(load(paste("clusterDivision/luadClusterData/cellNetwork01L/", str_replace_all(cellType,pattern = " ",replacement = "_"), ".RData", sep = "")))
  allNetSetCount = nrow(netSetNet01)
  highNetworkMarkerRatioLowest = 80
  netSetCountRatioLowest = 0.8
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
  satisfiedNetSetCount = 0
  
  # 在循环迭代过程中会有新的聚类出来的cluster被赋值给clusterRes_ls
  while(TRUE){
    cluster = clusterRes_ls[[1]]
    if(nrow(cluster) <= allNetSetCount * 0.01){
      break
    }
    # 判断当前cluster是否满足条件，getNetworkMarkerRatio_cluster、getEffectiveCluster()
    # 获得cellNetworkMarker01_cluster
    netSetName = row.names(cluster)
    netSetNet01_cluster = netSetNet01[netSetName,]
    networkResTable_cluster = getNetworkMarkerRatio_netSetCluster(netSetNet01_cluster)
    if (getEffectiveCluster(networkResTable_cluster)) {
      print(paste0("有效",nrow(cluster)))
      # 如果满足条件，将当前cluster添加到cellClusterRes_ls中
      netSetNet01_cluster$clusterRes = length(netSetClusterRes_ls)+1
      netSetClusterRes_ls = append(netSetClusterRes_ls, list(netSetNet01_cluster))
      
      satisfiedNetSetCount = satisfiedNetSetCount + nrow(netSetNet01_cluster)
      clusterRes_ls = clusterRes_ls[-1]
      # 设置迭代终止条件
      # 计算已有满足条件的cluster中的细胞数占总细胞数的比例
      netSetCountRatio = satisfiedNetSetCount/allNetSetCount
      
      # 如果满足条件的类别数目达到或超过指定比例，停止迭代
      if (netSetCountRatio >= netSetCountRatioLowest){
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
  
  
  cellClusterRes = getCellClusterRes(cellNetwork01, netSetClusterRes)
  
  # 计算networkMarker在每个cluster中的ratio
  networkResTableCluster = getNetworkMarkerRatio_netSet(cellNetwork01, cellClusterRes)
  write.csv(networkResTableCluster, paste("clusterDivision/resTableClusterL/networkResTableClusterL0730/networkResTableCluster_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""))
  
  
  highNetworkMarker = getHighNetworkMarker(networkResTableCluster)
  write.csv(highNetworkMarker, paste("clusterDivision/highNetworkMarkerTableLv2/highNetworkMarkerTable_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""))
  
  # 对聚类结果进行统计，包括clusterRes、cellCount_cluster、genePairs、cellCount_coexpress、Ratio
  resultStatistic = getResultStatistic(networkResTableCluster, cellClusterRes)
  write.csv(resultStatistic, paste("clusterDivision/resultStatisticLv2/resultStatistic_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""))
  print(cellType)
}