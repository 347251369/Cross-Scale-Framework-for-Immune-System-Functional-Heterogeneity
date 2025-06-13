# 准备工作
# 载入输入数据和变量设置
cellType = "Basophils"
cellNetwork01 = get(load(paste("clusterDivision/luadClusterData/cellNetwork01L/", str_replace_all(cellType,pattern = " ",replacement = "_"), ".RData", sep = "")))
allCellCount = nrow(cellNetwork01)


# 输入数据处理
netSetNet01 = getNetSetNet01(cellNetwork01)

# 获得cell对应的netSet，方便后面利用netSet获取cell
cellNetwork01_netSet = getNetSetRes(cellNetwork01)

# 修改cellClusterRes
cellClusterRes = read.csv(paste("clusterDivision/luadClusterRes/networkClusterResLv2/cellClusterRes_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""),
                          check.names = F, row.names = 1)
cellClusterRes$clusterRes[which(cellClusterRes$clusterRes %in% c(1))] = 1
cellClusterRes$clusterRes[which(cellClusterRes$clusterRes %in% c(2))] = 2
cellClusterRes$clusterRes[which(cellClusterRes$clusterRes %in% c(3))] = 3
cellClusterRes$clusterRes[which(cellClusterRes$clusterRes %in% c(4))] = 4
cellClusterRes$clusterRes[which(cellClusterRes$clusterRes %in% c(5))] = 0
cellClusterRes$clusterRes[which(cellClusterRes$clusterRes %in% c(6))] = 5
cellClusterRes$clusterRes[which(cellClusterRes$clusterRes %in% c(7))] = 6

cellClusterRes$clusterRes[which(cellClusterRes$clusterRes %in% c(8))] = 7
cellClusterRes$clusterRes[which(cellClusterRes$clusterRes %in% c(9))] = 8
cellClusterRes$clusterRes[which(cellClusterRes$clusterRes %in% c(11))] = 10
cellClusterRes$clusterRes[which(cellClusterRes$clusterRes %in% c(15))] = 0
cellClusterRes$clusterRes[which(cellClusterRes$clusterRes %in% c(16))] = 12
cellClusterRes$clusterRes[which(cellClusterRes$clusterRes %in% c(17))] = 13
cellClusterRes$clusterRes[which(cellClusterRes$clusterRes %in% c(13))] = 0

write.csv(cellClusterRes, paste("clusterDivision/luadAllDataProcessed/cellClusterRes_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""))

# 修改netSetClusterRes
netSetClusterRes = read.csv(paste("clusterDivision/luadClusterRes/networkClusterResLv2/netSetClusterRes_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""),
                            check.names = F, row.names = 1)
netSetClusterRes = read.csv(paste("clusterDivision/luadAllDataProcessed/netSetClusterRes_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""),
                            check.names = F, row.names = 1)
# cluster1聚类
cluster1 = netSetClusterRes[which(netSetClusterRes$clusterRes == 1),]
cluster1 = cluster1[,-ncol(cluster1)]
set.seed(1000)
newClusterRes = kmeans(cluster1, centers = 2, nstart = 1, iter.max = 10)
newClusterRes_df = as.data.frame(newClusterRes$cluster)
names(newClusterRes_df) = c("clusterRes")
#table(newClusterRes_df$clusterRes)
newClusterRes_df$clusterRes = as.numeric(newClusterRes_df$clusterRes)
netSet1 = row.names(newClusterRes_df)[which(newClusterRes_df$clusterRes == 1)]
netSet2 = row.names(newClusterRes_df)[which(newClusterRes_df$clusterRes == 2)]
cluster1A = netSetNet01[netSet1,]
cluster1B = netSetNet01[netSet2,]
ratio1A = getNetworkMarkerRatio_cluster(cluster1A)
ratio1B = getNetworkMarkerRatio_cluster(cluster1B)
# 修改cluster7结果
cluster7 = netSetClusterRes[which(netSetClusterRes$clusterRes == 7),]
cluster7$clusterRes = 0
cluster7$clusterRes[which(row.names(cluster7) %in% row.names(cluster7A))] = 7
cluster7$clusterRes[which(row.names(cluster7) %in% row.names(cluster7B))] = 18
# 更新netSetClusterRes中cluster7的结果
netSetClusterRes$clusterRes[which(netSetClusterRes$clusterRes == 7)] = cluster7$clusterRes

biomarker = row.names(networkBiomarkerRatio)

bioGene1 = unique(unlist(strsplit(biomarker[1], split = "_")))
bioGene2 = unique(unlist(strsplit(biomarker[2], split = "_")))

bioGene3 = unique(unlist(strsplit(biomarker[3], split = "_")))
bioGene4 = unique(unlist(strsplit(biomarker[4], split = "_")))
bioGene5 = unique(unlist(strsplit(biomarker[5], split = "_")))
bioGene6 = unique(unlist(strsplit(biomarker[6], split = "_")))
bioGene7 = unique(unlist(strsplit(biomarker[7], split = "_")))
bioGene8 = unique(unlist(strsplit(biomarker[8], split = "_")))
bioGene9 = unique(unlist(strsplit(biomarker[9], split = "_")))
bioGene10 = unique(unlist(strsplit(biomarker[10], split = "_")))
bioGene11 = unique(unlist(strsplit(biomarker[11], split = "_")))
bioGene12 = unique(unlist(strsplit(biomarker[12], split = "_")))

# bioGene4 = unique(unlist(strsplit(biomarker[4], split = "_")))


bioGene = unique(c(bioGene3,bioGene4))
loc1 = which(netSetClusterRes$clusterRes %in% c(8))
loc2 = which(rowSums(netSetGene01[,bioGene]) %in% length(bioGene))
loc = intersect(loc1, loc2)
netSetClusterRes$clusterRes[loc] = 8
netSetClusterRes$clusterRes[loc] = 6

bioGene = unique(c(bioGene1))
loc1 = which(cellClusterRes$clusterRes %in% c(5))
loc2 = which(rowSums(cellGene01[,bioGene]) %in% length(bioGene))
loc = intersect(loc1, loc2)
cellClusterRes$clusterRes[loc] = 9



netSetClusterRes$clusterRes[which(netSetClusterRes$clusterRes %in% c(1,2))] = 1
netSetClusterRes$clusterRes[which(netSetClusterRes$clusterRes %in% c(3))] = 2
netSetClusterRes$clusterRes[which(netSetClusterRes$clusterRes %in% c(4))] = 3
netSetClusterRes$clusterRes[which(netSetClusterRes$clusterRes %in% c(5,10))] = 4
netSetClusterRes$clusterRes[which(netSetClusterRes$clusterRes %in% c(6,11))] = 5

netSetClusterRes$clusterRes[which(netSetClusterRes$clusterRes %in% c(7))] = 6
netSetClusterRes$clusterRes[which(netSetClusterRes$clusterRes %in% c(8))] = 7
netSetClusterRes$clusterRes[which(netSetClusterRes$clusterRes %in% c(9))] = 8
netSetClusterRes$clusterRes[which(netSetClusterRes$clusterRes %in% c(13))] = 10
netSetClusterRes$clusterRes[which(netSetClusterRes$clusterRes %in% c(15))] = 0
netSetClusterRes$clusterRes[which(netSetClusterRes$clusterRes %in% c(16))] = 12
netSetClusterRes$clusterRes[which(netSetClusterRes$clusterRes %in% c(17))] = 13
netSetClusterRes$clusterRes[which(netSetClusterRes$clusterRes %in% c(18))] = 14

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
cellClusterRes = getCellClusterRes(cellNetwork01, netSetClusterRes)

write.csv(cellClusterRes, paste("clusterDivision/luadAllDataProcessed/cellClusterRes_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""))

write.csv(netSetClusterRes, paste("clusterDivision/luadAllDataProcessed/netSetClusterRes_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""))

cellClusterRes = read.csv(paste("clusterDivision/luadAllDataProcessed/cellClusterRes_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""),
                          check.names = F, row.names = 1)
# 计算networkMarker在每个cluster中的ratio
networkResTableCluster = getNetworkMarkerRatio(cellNetwork01, cellClusterRes)
write.csv(networkResTableCluster, paste("clusterDivision/luadAllDataProcessed/networkResTableCluster_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""))


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
  networkMarkerSet = networkMarkerSet[c(1:10)]
  # 删除networkMarkerSet中的空格
  if(is.element("",networkMarkerSet)){
    networkMarkerSet = networkMarkerSet[-which(networkMarkerSet == "")] #排除空字符串
  }
  networkBiomarkerRatio_ls[[i]] = getNetworkBiomarkerRatio(networkMarkerSet, cellClusterRes)
}
networkBiomarkerRatio = do.call(rbind, networkBiomarkerRatio_ls)
write.csv(networkBiomarkerRatio, paste("clusterDivision/networkBiomarkerRatioL/networkBiomarkerRatio_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""))


highNetworkMarker = getHighNetworkMarker(networkResTableCluster)
write.csv(highNetworkMarker, paste("clusterDivision/pbmcAllDataProcessed/highNetworkMarkerTable_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""))

# 对聚类结果进行统计，包括clusterRes、cellCount_cluster、genePairs、cellCount_coexpress、Ratio
resultStatistic = getResultStatistic(networkResTableCluster, cellClusterRes)
write.csv(resultStatistic, paste("clusterDivision/pbmcAllDataProcessed/resultStatistic_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""))




if("0" %in% names(networkResTableCluster)){
  networkResTableCluster = networkResTableCluster[,-1]
}
goodNetworkMarker = findGoodNetworkMarkerV2(networkResTableCluster, 1/2)
write.csv(goodNetworkMarker, file = paste("clusterDivision/pbmcAllDataProcessed/goodNetworkMarker_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""))
goodNetworkMarkerv2 = findGoodNetworkMarker(networkResTableCluster, 3)
print(cellType)

networkResTableCluster = read.csv(paste("clusterDivision/resTableClusterL/networkResTableClusterLv2/networkResTableCluster_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""),
                                  check.names = F, row.names = 1)
networkResTableCluster$"2B" = ratio1B$networkMarkerRatio
# cluster1聚类
cluster1 = netSetClusterRes[which(netSetClusterRes$clusterRes == 2),]
cluster1 = cluster1[,-ncol(cluster1)]
set.seed(1000)
newClusterRes = kmeans(cluster1, centers = 2, nstart = 1, iter.max = 10)
newClusterRes_df = as.data.frame(newClusterRes$cluster)
names(newClusterRes_df) = c("clusterRes")
#table(newClusterRes_df$clusterRes)
newClusterRes_df$clusterRes = as.numeric(newClusterRes_df$clusterRes)
netSet1 = row.names(newClusterRes_df)[which(newClusterRes_df$clusterRes == 1)]
netSet2 = row.names(newClusterRes_df)[which(newClusterRes_df$clusterRes == 2)]
netSet3 = row.names(newClusterRes_df)[which(newClusterRes_df$clusterRes == 3)]
netSet4 = row.names(newClusterRes_df)[which(newClusterRes_df$clusterRes == 4)]
netSet5 = row.names(newClusterRes_df)[which(newClusterRes_df$clusterRes == 5)]

cluster1A = netSetNet01[netSet1,]
cluster1B = netSetNet01[netSet2,]
cluster1C = netSetNet01[netSet3,]
cluster1D = netSetNet01[netSet4,]
cluster1E = netSetNet01[netSet5,]

ratio1A = getNetworkMarkerRatio_cluster(cluster1A)
ratio1B = getNetworkMarkerRatio_cluster(cluster1B)
ratio1C = getNetworkMarkerRatio_cluster(cluster1C)
ratio1D = getNetworkMarkerRatio_cluster(cluster1D)
ratio1E = getNetworkMarkerRatio_cluster(cluster1E)

