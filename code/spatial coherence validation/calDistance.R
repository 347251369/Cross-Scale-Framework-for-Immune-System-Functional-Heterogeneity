# 该代码用于计算cluster之间的距离
# 输入注释好的空间转录组数据
## 输入数据格式为数据框，包括细胞、细胞位置坐标、细胞分类label
## 输出cluster之间的距离
clusterInfo
clusterLabel = paste(networkBiomarkerSummary$cellType, networkBiomarkerSummary$cluster, sep = "_")
# 计算两个cluster之间的距离----------------------
## 最小距离法----------------
cluster_distance_min = function(cluster1, cluster2, method = "euclidean") {
  min_distance = Inf
  for(i in 1:nrow(cluster1)){
    for(j in 1:nrow(cluster2)){
      distance = dist(rbind(as.vector(cluster1[i,]), as.vector(cluster2[j,])), method = method)
      min_distance = min(min_distance, distance)
    }
  }
  return(min_distance)
}
## 平均距离法----------------
cluster_distance_mean = function(cluster1, cluster2, method = "euclidean") {
  all_distance = c()
  for(i in 1:nrow(cluster1)){
    for(j in 1:nrow(cluster2)){
      distance = dist(rbind(as.vector(cluster1[i,]), as.vector(cluster2[j,])), method = method)
      all_distance = c(all_distance, as.vector(distance))
    }
  }
  return(mean(all_distance))
}
## 最大距离法----------------
cluster_distance_max = function(cluster1, cluster2, method = "euclidean") {
  max_distance = -Inf
  for(i in 1:nrow(cluster1)){
    for(j in 1:nrow(cluster2)){
      distance = dist(rbind(as.vector(cluster1[i,]), as.vector(cluster2[j,])), method = method)
      max_distance = max(max_distance, distance)
    }
  }
  return(max_distance)
}
## 类中心距离法(重心/质心法)--------------
### 计算簇的重心
calculate_centroid = function(cluster) {
  return(colMeans(cluster))
}
#### 计算两个重心之间的距离
cluster_distance_centroid = function(cluster1, cluster2, method = "euclidean") {
  centroid1 = calculate_centroid(cluster1)
  centroid2 = calculate_centroid(cluster2)
  distance = dist(rbind(centroid1, centroid2), method = method)
  return(as.numeric(distance))
}

# 加载annoNetworkBiomarker中的数据
sampleSet = c("B03609G3", "B03607G1", "B03612D3", "B03609G4", "B03613D6", "B03613C5")
sampleID = "B03613C5"
load(file = paste("spatialDataProcess/stTargetData_tissue/",sampleID,"_annoNetworkBiomarker.RData",sep = ""))
clusterName = c()
for(majorCelltype in names(binPositionRes)){
  clusterName = c(clusterName, names(binPositionRes[[majorCelltype]]))
}
distMatrix = matrix(data = 0, nrow = length(clusterName), ncol = length(clusterName))
row.names(distMatrix) = clusterName
colnames(distMatrix) = clusterName
# 计算cluster之间的距离，获得distMatrix
for(majorCelltype1 in names(binPositionRes)){
  for(networkBiomarker1 in names(binPositionRes[[majorCelltype1]])){
    cluster1 = binPositionRes[[majorCelltype1]][[networkBiomarker1]]
    for(majorCelltype2 in names(binPositionRes)){
      for(networkBiomarker2 in names(binPositionRes[[majorCelltype2]])){
        cluster2 = binPositionRes[[majorCelltype2]][[networkBiomarker2]]
        distMatrix[networkBiomarker1, networkBiomarker2] = cluster_distance_centroid(cluster1, cluster2)
      }
    }
  }
}
save(distMatrix, file = paste("spatialDataProcess/stTargetData_tissue/",sampleID,"_distMatrix.csv",sep = ""))

# 加载已知的clusterTOcluster，划分distMatrix中的距离对（0/1）
clusterTOcluster = read.csv(file = "spatialDataProcess/genePairStatistic/clusterTOcluster(threshold=10_0.7).csv",
                            row.names = 1, check.names = FALSE)
# 利用networkBiomarkerSummary将clusterTOcluster的行列名更改为对应的networkBiomarker
for(i in 1:nrow(clusterTOcluster)){
  celltype = unlist(str_split(row.names(clusterTOcluster)[i], "_"))[1]
  loc = unlist(str_split(row.names(clusterTOcluster)[i], "_"))[2]
  networkBiomarkerSet = networkBiomarkerSummary$networkBiomarker[networkBiomarkerSummary$cellType %in% celltype]
  networkBiomarker = networkBiomarkerSet[as.numeric(loc)]
  row.names(clusterTOcluster)[i] = networkBiomarker
}
colnames(clusterTOcluster) = row.names(clusterTOcluster)

group0 = c()
group1 = c()
for(i in row.names(distMatrix)){
  for(j in colnames(distMatrix)){
    if(clusterTOcluster[i,j] == 0){
      group0 = c(group0, distMatrix[i,j])
    }
    if(clusterTOcluster[i,j] == 1){
      print(c(i,j,distMatrix[i,j]))
      group1 = c(group1, distMatrix[i,j])
    }
  }
}

save(group0, group1, file = paste("spatialDataProcess/stTargetData_tissue/",sampleID,"_group.RData",sep = ""))


cluster1 = binPositionRes[["Neutrophils"]][["S100A8_S100A9"]]
cluster2 = binPositionRes[["Neutrophils"]][["S100A8_S100A9"]]

# 计算两个簇之间的距离（平均距离法）
distance = cluster_distance_mean(cluster1, cluster2)
distance = cluster_distance_centroid(cluster1, cluster2)
print(distance)

# 计算三个及三个以上cluster之间的距离------------
## 平均距离法----------------
### 计算两簇之间的平均距离
cluster_distance_mean = function(cluster1, cluster2, method = "euclidean") {
  all_distance = c()
  for(i in 1:nrow(cluster1)){
    for(j in 1:nrow(cluster2)){
      distance = dist(rbind(as.vector(cluster1[i,]), as.vector(cluster2[j,])), method = method)
      all_distance = c(all_distance, as.vector(distance))
    }
  }
  return(mean(all_distance))
}
### 计算多个簇之间的平均距离
mul_cluster_distance_mean = function(clusters) {
  num_clusters = length(clusters)
  distances = c()
  for (i in 1:num_clusters) {
    for (j in 1:num_clusters) {
      if(i != j){
        distances = c(distances, cluster_distance_mean(clusters[[i]], clusters[[j]]))
      }
    }
  }
  return(mean(distances))  # 返回距离的平均值
}



## 类中心距离法--------------
### 计算簇的重心
calculate_centroid = function(cluster) {
  return(colMeans(cluster))
}
### 计算重心之间的距离
mul_cluster_distance_centroid = function(clusters) {
  num_clusters = length(clusters)
  centroids = lapply(clusters, calculate_centroid)
  distances = c()
  # 计算每两个重心之间的距离
  for (i in 1:num_clusters) {
    for (j in 1:num_clusters) {
      if(i != j){
        distances = c(distances, sqrt(sum((centroids[[i]] - centroids[[j]])^2)))
      }
      
    }
  }
  return(mean(distances))  # 返回距离的平均值
}
