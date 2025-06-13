# 获得网络中三个基因的组合
# geneList为一个网络中有两个基因的全部geneList
getAllGeneList3 = function(geneList){
  k = 2
  gene_listTmp = c()
  for (i in geneList) {
    for(j in geneList){
      if(length(intersect(i,j)) == (k-1)){
        listTmp = unique(union(i,j))
        gene_listTmp = c(gene_listTmp,list(listTmp[order(listTmp)])) 
      }
    }
  }
  return(gene_listTmp)
}

# 寻找关联
# getCoExpressCondition:该函数为了提取出因果网络中基因共表达的情况
# geneListAll:gene_list_all是组合完成的所有基因组合
# staticForCellTypeDF：对所有细胞的类型进行聚合统计后的DF（dataframe）
# expressionCount：单细胞基因表达值，只包含有细胞注释的细胞和可用的基因
# allCellType：每个细胞对应的细胞类型
getCoExpressCondition = function(geneListAll, staticForCellTypeDF, expressionCount, allCellType, targetCellType){
  # 建表。
  resTable = as.data.frame(matrix(ncol = dim(staticForCellTypeDF)[1],nrow = length(geneListAll)))
  colnames(resTable) = staticForCellTypeDF$Var1
  # 行名命名
  rownames(resTable) = unlist(lapply(geneListAll,function(x){str_c(x,collapse = "_")})) 
  resTable[is.na(resTable)] = 0
  
  rowFlag = 1
  for (i in 1:length(geneListAll)) {
    # print(paste0("****",rowFlag,"****"))
    tmpNet = geneListAll[[i]]
    
    # 把网络中的基因的表达值提取出来
    tmpCount = expressionCount[tmpNet,]
    
    # 提取共表达网络的细胞
    cellId = colnames(tmpCount)
    for(j in 1:dim(tmpCount)[1]){
      cellId = intersect(cellId,colnames(tmpCount)[tmpCount[j,]>0])
      if(length(cellId) == 0){
        break
      }
    }
    
    # 没有细胞共表达该网络
    if(length(cellId) == 0){
      rowFlag = rowFlag + 1
      next
    }
    
    # 统计这些细胞的细胞类型
    # 根据细胞id提取细胞类型
    tmpCelltype = allCellType[cellId,]
    # 统计细胞id对应的细胞类型
    staticsFortmpCelltype = as.data.frame(table(tmpCelltype[,targetCellType]))
    staticsFortmpCelltype$Var1 = as.character(staticsFortmpCelltype$Var1)
    
    for(j in 1:dim(staticsFortmpCelltype)[1]){
      resTable[rowFlag,staticsFortmpCelltype$Var1[j]] = staticsFortmpCelltype$Freq[j]
    }
    rowFlag = rowFlag + 1
  }
  return(resTable)
}

library(tidyr)
# 贝叶斯网络拆分和子网络ratio计算
getSubnetworkRatio = function(targetObject, targetCellType, bayesGenePair){
  # 加载targetObject的所有meta.data
  allCellType = targetObject@meta.data
  # 排除unknown的细胞（即不关心的细胞类型）
  allCellTypeNoOther = allCellType[which(allCellType[,c(targetCellType)]!="unknown"),]
  # 挑选meta.data中感兴趣的列，因为单细胞数据之间的meta.data不统一，故此处不操作
  # allCellTypeNoOther = allCellTypeNoOther[,c(9,10,11,13,14)]
  # 统计每种细胞的数量
  staticForCellType = as.data.frame(table(allCellTypeNoOther[,c(targetCellType)])) 
  staticForCellType$Var1 = as.character(staticForCellType$Var1)
  
  # 处理counts矩阵,保留非unknown的细胞
  expression_count = targetObject@assays$RNA@counts
  expressionCountNoOther = expression_count[,rownames(allCellTypeNoOther)]
  
  # 处理筛选后的贝叶斯分析结果，文件包含列：from,to,strength,direction
  # 合并from，to两列，用"%%"隔开，方便下面split隔开得到list
  bayesGenePair = unite(bayesGenePair,"fromto", from, to, sep = "%%", remove = F)
  
  bayesGenePair = str_split(bayesGenePair$fromto, pattern = "%%")
  # 将list中的每一个排序
  bayesGenePair = lapply(bayesGenePair, function(x){x=x[order(x)]})
  bayesGenePair = unique(bayesGenePair)
  
  # 把单细胞测序数据中没有的网络去掉
  removeNet = c()
  for(i in 1:length(bayesGenePair)){
    if(length(intersect(rownames(expressionCountNoOther),bayesGenePair[[i]]))<2){
      removeNet=c(removeNet,i)
    }
  }
  bayesGenePair2 = bayesGenePair[-removeNet]
  geneAvailable = unique(unlist(bayesGenePair2)) 
  bayesGenePair_all = bayesGenePair2
  # 获得网络中三个基因的组合
  # geneList为一个网络中有两个基因的全部geneList
  bayesGenePair3 = getAllGeneList3(bayesGenePair2)
  bayesGenePair3 = unique(bayesGenePair3)
  #2个基因的网络和3个基因的网络合在一起
  bayesGenePair_all = c(bayesGenePair2, bayesGenePair3)
  expressionCountNoOtherAvailable = expressionCountNoOther[geneAvailable,]
  # 开始寻找关联 getCoExpressCondition
  restable = getCoExpressCondition(bayesGenePair_all, staticForCellType, expressionCountNoOtherAvailable, allCellTypeNoOther, targetCellType)
  #将结果转换成百分比
  resTable = restable
  for(i in 1:dim(resTable)[2]){
    resTable[,i] = round(resTable[,i]/staticForCellType$Freq[i]*100, 2)
  }
  return(resTable)
}

# 筛选resTable
#如果一行中最大值大于等于均值+nsd倍方差，那么这一行被称作goodRows
#输入：得到的resTable
#输出：筛选后的resTable
findGoodRows = function(table, nsd){
  goodRow = c()
  for(i in 1:nrow(table)){
    meantmp = mean(as.numeric(table[i,]))
    sdtmp = sd(as.numeric(table[i,]))
    if(max(as.numeric(table[i,]))>0 & max(as.numeric(table[i,])) >= meantmp+nsd*sdtmp){
      goodRow = c(goodRow,i)
    }
  }
  return(table[goodRow,])
}

# networkMarker分类
#使用findGoodRows函数后得到想要的网络，但是还没归类
#输入：使用findGoodRows函数后得到的goodRows
#输出：归类好的网络
getHighestCellTypeAndNet = function(goodRows){
  res = as.data.frame(matrix(nrow = nrow(goodRows),ncol = ncol(goodRows)))
  colnames(res) = colnames(goodRows)
  maxRow = 0
  #一些行的最大值相同，提取出来看看怎么回事
  confusedRow = c()
  for(i in 1:nrow(goodRows)){
    colIdx = which(goodRows[i,] == max(goodRows[i,]))
    if(length(colIdx)>1){
      confusedRow = c(confusedRow,i)
      next
    }
    for(j in 1:nrow(res)){
      if(is.na(res[,colIdx][j])){
        res[,colIdx][j] = rownames(goodRows)[i]
        maxRow = max(maxRow,j)
        break
      }
    }
  }
  res[is.na(res)] = ""
  res = res[1:maxRow,]
  confusedMatrix<<-goodRows[confusedRow,]
  return(res)
}


# 加载贝叶斯分析结果
bayesGenePairPBMC = read.csv("bayesianAnalysis/bayesGenePairPBMC.csv")
resTablePBMC = getSubnetworkRatio(pbmc, "label1", bayesGenePairPBMC)
# 保存resTable
write.csv(resTablePBMC, "networkMarkerIdentification/resTablePBMC.csv")
goodRowsPBMC = findGoodRows(resTablePBMC, 3)
# 保存goodRows
write.csv(goodRowsPBMC, "networkMarkerIdentification/goodRowsPBMC.csv")
networkMarkePBMC = getHighestCellTypeAndNet(goodRowsPBMC)
# 保存networkMarker
write.csv(networkMarkerPBMC, "networkMarkerIdentification/networkMarkerPBMC.csv")

# 整理network marker数据
cellTypeSet = c("memory CD4+ T cells","memory CD8+ T cells","naive CD4+ T cells",
                "naive CD8+ T cells","T helper cells","T effector cells","NKT","MAIT","Tregs",
                "NK","Plasma cells","pre-B cells","cDC","pDC","Monocytes",
                "Basophils","Neutrophils","Macrophages")
ratio_ls = list()
# 创建工作簿
wb = createWorkbook()
for(i in c(cellTypeSet)){
  cellType = i
  #取对应细胞类型的网络列表
  networkMarker = networkMarkeLUAD[,cellType]
  if(is.element("",networkMarker)){
    networkMarker = networkMarker[-which(networkMarker == "")] #排除空字符串
  }
  # 提取网络的所有得分
  ratio = goodRowsLUAD[networkMarker,]
  # 排序
  order = order(ratio[,cellType], decreasing = T)
  ratio = ratio[order,]
  # 添加到list中
  ratio_ls[[cellType]] = ratio
  # 创建工作表
  addWorksheet(wb, sheetName = cellType)
  # 向工作表写入数据
  writeData(wb, sheet = cellType, ratio, rowNames = TRUE)
}
# 保存工作簿
saveWorkbook(wb, "networkMarkerIdentification/networkMarkerRatioLUAD.xlsx",
             overwrite = TRUE)


# 并行计算resTable
targetObject = pbmc
targetCellType = "label1" # LUAD是label0525
bayesGenePair = bayesGenePairPBMC
# 加载targetObject的所有meta.data
allCellType = targetObject@meta.data
# 排除unknown的细胞（即不关心的细胞类型）
allCellTypeNoOther = allCellType[which(allCellType[,c(targetCellType)]!="unknown"),]
# 统计每种细胞的数量
staticForCellType = as.data.frame(table(allCellTypeNoOther[,c(targetCellType)])) 
staticForCellType$Var1 = as.character(staticForCellType$Var1)
# 处理counts矩阵,保留非unknown的细胞
expression_count = targetObject@assays$RNA@counts
expressionCountNoOther = expression_count[,rownames(allCellTypeNoOther)]
# 处理筛选后的贝叶斯分析结果，文件包含列：from,to,strength,direction
# 合并from，to两列，用"%%"隔开，方便下面split隔开得到list
# 读取因果网络文件，文件包含列：from,to,strength,direction
gene_list_sheet = bayesGenePair
# 筛选网络中的基因pairs
gene_list_sheet = gene_list_sheet[gene_list_sheet$direction >= 0.5,]
geneListSheet = gene_list_sheet
gene_list_sheet = geneListSheet[geneListSheet$strength >= 0.7,]
#合并from，to两列，用"%%"隔开，方便下面split隔开得到list
gene_list_sheet = unite(gene_list_sheet,"fromto",from,to,sep = "%%",remove = F)
gene_list = str_split(gene_list_sheet$fromto,pattern = "%%")
#将list中的每一个排序
gene_list = lapply(gene_list,function(x){x=x[order(x)]})
gene_list = unique(gene_list)
#把单细胞测序数据中没有的网络去掉
removeNet = c()
for(i in 1:length(gene_list)){
  if(length(intersect(rownames(expressionCountNoOther), gene_list[[i]]))<2){
    removeNet=c(removeNet,i)
  }
}
gene_list = gene_list[-removeNet]
geneAvailable = unique(unlist(gene_list))
gene_list_all = gene_list
# 获得网络中三个基因的组合
# geneList为一个网络中有两个基因的全部geneList
gene_list3 = getAllGeneList3(gene_list)
gene_list3 = unique(gene_list3)
#2个基因的网络和3个基因的网络合在一起
gene_list_all = c(gene_list, gene_list3)
expressionCountNoOtherAvailable = expressionCountNoOther[geneAvailable,]
# 将gene_list划分成20份，并行计算resTable
library(parallel)
cl = makeCluster(20)
geneListALL = gene_list_all
# 共享数据
clusterExport(cl, c("staticForCellType", "expressionCountNoOtherAvailable", "allCellTypeNoOther", "targetCellType"))
clusterEvalQ(cl, {
  getCoExpressCondition = function(geneListAll, staticForCellTypeDF, expressionCount, allCellType, targetCellType){
    # 建表。
    resTable = as.data.frame(matrix(ncol = dim(staticForCellTypeDF)[1],nrow = length(geneListAll)))
    colnames(resTable) = staticForCellTypeDF$Var1
    # 行名命名
    rownames(resTable) = unlist(lapply(geneListAll,function(x){str_c(x,collapse = "_")})) 
    resTable[is.na(resTable)] = 0
    
    rowFlag = 1
    for (i in 1:length(geneListAll)) {
      # print(paste0("****",rowFlag,"****"))
      tmpNet = geneListAll[[i]]
      
      # 把网络中的基因的表达值提取出来
      tmpCount = expressionCount[tmpNet,]
      
      # 提取共表达网络的细胞
      cellId = colnames(tmpCount)
      for(j in 1:dim(tmpCount)[1]){
        cellId = intersect(cellId,colnames(tmpCount)[tmpCount[j,]>0])
        if(length(cellId) == 0){
          break
        }
      }
      
      # 没有细胞共表达该网络
      if(length(cellId) == 0){
        rowFlag = rowFlag + 1
        next
      }
      
      # 统计这些细胞的细胞类型
      # 根据细胞id提取细胞类型
      tmpCelltype = allCellType[cellId,]
      # 统计细胞id对应的细胞类型
      staticsFortmpCelltype = as.data.frame(table(tmpCelltype[,targetCellType]))
      staticsFortmpCelltype$Var1 = as.character(staticsFortmpCelltype$Var1)
      
      for(j in 1:dim(staticsFortmpCelltype)[1]){
        resTable[rowFlag,staticsFortmpCelltype$Var1[j]] = staticsFortmpCelltype$Freq[j]
      }
      rowFlag = rowFlag + 1
    }
    return(resTable)
  }
  library(stringr)
})
result = parSapply(cl, split(head(geneListALL, 84680), rep(1:20)), function(x) {
  restableX = getCoExpressCondition(x, staticForCellType, expressionCountNoOtherAvailable, allCellTypeNoOther, targetCellType)
})
# 整理result, 里面包括20个restableX
restableX = as.data.frame(matrix(ncol = dim(staticForCellType)[1], nrow = length(head(geneListALL, 84680))))
colnames(restableX) = staticForCellType$Var1
row.names(restableX) = unlist(lapply(head(geneListALL, 84680), function(x){str_c(x,collapse = "_")}))
for(i in row.names(result)){
  ratioi = c()
  for(j in colnames(result)){
    ratioi = c(ratioi, unlist(result[i,j]))
  }
  restableX[,i] = ratioi
}
#将结果转换成百分比
resTablePBMC = restableX
for(i in 1:dim(resTablePBMC)[2]){
  resTablePBMC[,i] = round(resTablePBMC[,i]/staticForCellType$Freq[i]*100, 2)
}
# 保存resTable
write.csv(resTablePBMC, "networkMarkerIdentification/resTablePBMC.csv")
goodRowsLUADv2 = findGoodRows(resTableLUAD, 3)
# 保存goodRows
write.csv(goodRowsLUADv2, "networkMarkerIdentification/goodRowsLUAD.csv")

# 并行计算networkMarker--------------------------------------------------------
library(parallel)
cl = makeCluster(20)
# 共享数据
clusterExport(cl, c("goodRowsPBMC"))
clusterEvalQ(cl, {
  getHighestCellTypeAndNet = function(goodRows){
    res = as.data.frame(matrix(nrow = nrow(goodRows),ncol = ncol(goodRows)))
    colnames(res) = colnames(goodRows)
    maxRow = 0
    
    #一些行的最大值相同，提取出来看看怎么回事
    confusedRow = c()
    
    for(i in 1:nrow(goodRows)){
      colIdx = which(goodRows[i,] == max(goodRows[i,]))
      if(length(colIdx)>1){
        confusedRow = c(confusedRow,i)
        next
      }
      for(j in 1:nrow(res)){
        if(is.na(res[,colIdx][j])){
          res[,colIdx][j] = rownames(goodRows)[i]
          maxRow = max(maxRow,j)
          break
        }
      }
    }
    res[is.na(res)] = ""
    res = res[1:maxRow,]
    
    confusedMatrix<<-goodRows[confusedRow,]
    return(res)
  }
})
resultNetworkMarker = parSapply(cl, split(goodRowsPBMC, rep(1:20)), function(x) {
  restableNetworkMarker = getHighestCellTypeAndNet(x)
})
# 整理result, 里面包括20个restableX
restableNetworkMarker = as.data.frame(matrix(nrow = nrow(goodRowsPBMC), ncol = ncol(goodRowsPBMC)))
colnames(restableNetworkMarker) = colnames(goodRowsPBMC)
maxRow = 0
for(i in row.names(resultNetworkMarker)){
  networkMarkeri = c()
  for(j in colnames(resultNetworkMarker)){
    networkMarkeri = c(networkMarkeri, unlist(resultNetworkMarker[i,j]))
  }
  restableNetworkMarker[,i] = c(networkMarkeri, rep("", nrow(goodRowsPBMC)-length(networkMarkeri)))
  if(length(networkMarkeri) > maxRow){
    maxRow = length(networkMarkeri)
  }
}
restableNetworkMarker = restableNetworkMarker[1:maxRow,]
networkMarkerPBMC = restableNetworkMarker
# 保存networkMarker
write.csv(networkMarkerPBMC, "networkMarkerIdentification/networkMarkerPBMC.csv")

# 整理network marker数据
cellTypeSet = c("memory CD4+ T cells","memory CD8+ T cells","naive CD4+ T cells",
                "naive CD8+ T cells","NKT","MAIT","Tregs",
                "NK","Plasma cells","pre-B cells","cDC","pDC","Monocytes",
                "Basophils","Neutrophils")
ratio_ls = list()
# 创建工作簿
wb = createWorkbook()
for(i in c(cellTypeSet)){
  cellType = i
  #取对应细胞类型的网络列表
  networkMarker = networkMarkerPBMC[,cellType]
  if(is.element("",networkMarker)){
    networkMarker = networkMarker[-which(networkMarker == "")] #排除空字符串
  }
  # 提取网络的所有得分
  ratio = goodRowsPBMC[networkMarker,]
  # 排序
  order = order(ratio[,cellType], decreasing = T)
  ratio = ratio[order,]
  # 添加到list中
  ratio_ls[[cellType]] = ratio
  # 创建工作表
  addWorksheet(wb, sheetName = cellType)
  # 向工作表写入数据
  writeData(wb, sheet = cellType, ratio, rowNames = TRUE)
}
# 保存工作簿
saveWorkbook(wb, "networkMarkerIdentification/networkMarkerRatioPBMC.xlsx",
             overwrite = TRUE)
