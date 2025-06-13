
# PBMC-networkCell01P
library(Seurat)
library(ggplot2)
library(stringr)
setwd("/home/tzy/project/clusterDivision")
allCelltypeName = c("memory CD4+ T cells","memory CD8+ T cells","naive CD4+ T cells",
                    "naive CD8+ T cells","NKT","Tregs",
                    "NK","Plasma cells","cDC","pDC","Monocytes",
                    "Basophils","Neutrophils")
networkMarker = networkMarkerPBMC
#循环celltype
for(a in allCelltypeName){
  cellType = a
  targetObject = subset(pbmc, label1 %in% a)
  #取对应细胞类型的网络列表
  cellTypeNetTmp = networkMarker[,cellType]
  if(is.element("",cellTypeNetTmp)){
    cellTypeNetTmp = cellTypeNetTmp[-which(cellTypeNetTmp == "")] #排除空字符串
  }
  if(length(cellTypeNetTmp) == 0){#如果排除空字符串后空了就进入下一个celltype进行循环
    next
  }
  
  
  
  #给每个网络共表达的细胞计数
  cellCountForcellTypeNetTmp = as.data.frame(matrix(nrow = length(cellTypeNetTmp),ncol = 2))
  colnames(cellCountForcellTypeNetTmp) = c("netIdx","cellCnt")
  #共表达umap图在一个list中
  gList=c()
  
  #提取细胞，网络的01矩阵，用于聚类
  dataForCluster = as.data.frame(matrix(
    nrow = nrow(targetObject@meta.data),
    ncol = length(cellTypeNetTmp)
  ))
  colnames(dataForCluster) = cellTypeNetTmp
  base::rownames(dataForCluster) = base::rownames(targetObject@meta.data)
  
  
  for(i in cellTypeNetTmp){
    if(i == ""){
      break
    }
    #net的行索引
    netIdx = which(cellTypeNetTmp == i)
    netTmp = str_split(i,pattern = "_")
    featurePlotData = c()
    #将网络的基因表达数据提取出来
    for(j in netTmp[[1]]){
      featurePlotDataTmp = FeaturePlot(targetObject,j)$data
      if(is.null(featurePlotData)){
        featurePlotData = featurePlotDataTmp
        colnames(featurePlotData)[ncol(featurePlotData)] = j
      }else{
        featurePlotData = cbind(featurePlotData,featurePlotDataTmp[,4])
        colnames(featurePlotData)[ncol(featurePlotData)] = j
      }
    }
    #将网络中基因的共表达细胞提取出来
    for(j in 4:ncol(featurePlotData)){
      if(is.null(featurePlotData)){
        gList = c(gList,list(NULL))
        break
      }else{
        featurePlotData = featurePlotData[featurePlotData[,j]>0,] 
      }
    }
    
    if(!is.null(featurePlotData)){
      umap_dat = DimPlot(targetObject)
      umap = umap_dat$data
      
      g = ggplot(umap) + 
        geom_point(aes(x =UMAP_1 , y = UMAP_2),data = umap,color="grey",size = 0.2, alpha = 1) + 
        geom_point(aes(x =UMAP_1 , y = UMAP_2,color=zscore),data = na.omit(featurePlotData),size = 0.2, alpha = 1)+
        scale_color_gradient(low="lightblue", high="darkblue")+
        ggtitle(i)+theme(plot.title = element_text(hjust = 0.5,size = 8))   
      #g=drawPicWithCluster(colors = colorSet,clusterFile = clusterFile,umapData = umap,featurePlotData = featurePlotData,
      #                     highExpCells = highExpCells,title = i)
      
      cellCountForcellTypeNetTmp[netIdx,1] = netIdx
      cellCountForcellTypeNetTmp[netIdx,2] = nrow(featurePlotData)
      if(is.null(gList)){
        dimplot = DimPlot(targetObject)
        dimplot = dimplot$data
        d = ggplot(dimplot) + 
          geom_point(aes(x = UMAP_1 , y = UMAP_2),data = dimplot,color="grey",size = 0.2, alpha = 1) + 
          geom_point(aes(x = UMAP_1 , y = UMAP_2),
                     data = dimplot[rownames(targetObject@meta.data[which(targetObject@meta.data$label0525 == cellType),]),],
                     color='blue',size = 0.2, alpha = 1)
        #d
        dimplot = DimPlot(targetObject, label = T)
        #dimplot
        #ggsave(filename = paste0(str_replace_all(cellType,pattern = " ",replacement = "_"),"/","dimplot",".png"),width = 40,height = 20,units = "cm")
        gList = list(g)
      }else{
        gList = c(gList,list(g))
      }
      
      #dimplot=DimPlot(targetObject,label = T)
      #dimplot+g+d
      #ggsave(filename = paste0(str_replace_all(cellType,pattern = " ",replacement = "_"),"/",i,".png"),width = 40,height = 20,units = "cm")
    }
  }
  
  cellCountForcellTypeNetTmp = na.omit(cellCountForcellTypeNetTmp)
  cellCountForcellTypeNetTmp = cellCountForcellTypeNetTmp[order(cellCountForcellTypeNetTmp$cellCnt, decreasing = T),]
  #挑选细胞数排序前3/4的网络图
  #netSelect=cellCountForcellTypeNetTmp$netIdx[c(1:ceiling(nrow(cellCountForcellTypeNetTmp)*3/4))] 
  #挑选细胞数不为0的网络图
  netSelect = cellCountForcellTypeNetTmp[which(cellCountForcellTypeNetTmp$cellCnt > 0),]$netIdx
  #修改dataForCluster
  #提取细胞，网络的01矩阵，用于聚类
  dataForCluster = as.data.frame(matrix(0,
    nrow = nrow(targetObject@meta.data),
    ncol = length(netSelect)
  ))
  colnames(dataForCluster) = cellTypeNetTmp[netSelect]
  base::rownames(dataForCluster) = base::rownames(targetObject@meta.data)
  
  #提取用于聚类的数据
  for(netSelecti in netSelect){
    netSelectiData = gList[[netSelecti]][["layers"]][[2]][["data"]]
    netSelectiNet = gList[[netSelecti]][["labels"]][["title"]]
    dataForCluster[base::rownames(netSelectiData),netSelectiNet] = 1
  }
  dataForCluster[is.na(dataForCluster)] = 0
  dataForCluster = dataForCluster[rowSums(dataForCluster)>0,]
  save(dataForCluster,file = paste0("pbmcClusterData/cellNetwork01P/", str_replace_all(cellType,pattern = " ",replacement = "_"),".RData"))
  print(cellType)
}


setwd("/home/tzy/project/clusterDivision")
allCelltypeName = c("memory CD4+ T cells","memory CD8+ T cells","naive CD4+ T cells",
                    "naive CD8+ T cells", "T helper cells","T effector cells", "NKT","Tregs",
                    "NK","Plasma cells","cDC","pDC","Monocytes",
                    "Basophils","Neutrophils","Macrophages")
networkMarker = networkMarkerLUAD
#循环celltype
for(a in allCelltypeName){
  cellType = a
  targetObject = subset(luad, label0525 %in% a)
  #取对应细胞类型的网络列表
  cellTypeNetTmp = networkMarker[,cellType]
  if(is.element("",cellTypeNetTmp)){
    cellTypeNetTmp = cellTypeNetTmp[-which(cellTypeNetTmp == "")] #排除空字符串
  }
  if(length(cellTypeNetTmp) == 0){#如果排除空字符串后空了就进入下一个celltype进行循环
    next
  }
  
  
  
  #给每个网络共表达的细胞计数
  cellCountForcellTypeNetTmp = as.data.frame(matrix(nrow = length(cellTypeNetTmp),ncol = 2))
  colnames(cellCountForcellTypeNetTmp) = c("netIdx","cellCnt")
  #共表达umap图在一个list中
  gList=c()
  
  #提取细胞，网络的01矩阵，用于聚类
  dataForCluster = as.data.frame(matrix(
    nrow = nrow(targetObject@meta.data),
    ncol = length(cellTypeNetTmp)
  ))
  colnames(dataForCluster) = cellTypeNetTmp
  base::rownames(dataForCluster) = base::rownames(targetObject@meta.data)
  
  
  for(i in cellTypeNetTmp){
    if(i == ""){
      break
    }
    #net的行索引
    netIdx = which(cellTypeNetTmp == i)
    netTmp = str_split(i,pattern = "_")
    featurePlotData = c()
    #将网络的基因表达数据提取出来
    for(j in netTmp[[1]]){
      featurePlotDataTmp = FeaturePlot(targetObject,j)$data
      if(is.null(featurePlotData)){
        featurePlotData = featurePlotDataTmp
        colnames(featurePlotData)[ncol(featurePlotData)] = j
      }else{
        featurePlotData = cbind(featurePlotData,featurePlotDataTmp[,4])
        colnames(featurePlotData)[ncol(featurePlotData)] = j
      }
    }
    #将网络中基因的共表达细胞提取出来
    for(j in 4:ncol(featurePlotData)){
      if(is.null(featurePlotData)){
        gList = c(gList,list(NULL))
        break
      }else{
        featurePlotData = featurePlotData[featurePlotData[,j]>0,] 
      }
    }
    
    if(!is.null(featurePlotData)){
      umap_dat = DimPlot(targetObject)
      umap = umap_dat$data
      
      g = ggplot(umap) + 
        geom_point(aes(x =UMAP_1 , y = UMAP_2),data = umap,color="grey",size = 0.2, alpha = 1) + 
        geom_point(aes(x =UMAP_1 , y = UMAP_2,color=zscore),data = na.omit(featurePlotData),size = 0.2, alpha = 1)+
        scale_color_gradient(low="lightblue", high="darkblue")+
        ggtitle(i)+theme(plot.title = element_text(hjust = 0.5,size = 8))   
      #g=drawPicWithCluster(colors = colorSet,clusterFile = clusterFile,umapData = umap,featurePlotData = featurePlotData,
      #                     highExpCells = highExpCells,title = i)
      
      cellCountForcellTypeNetTmp[netIdx,1] = netIdx
      cellCountForcellTypeNetTmp[netIdx,2] = nrow(featurePlotData)
      if(is.null(gList)){
        dimplot = DimPlot(targetObject)
        dimplot = dimplot$data
        d = ggplot(dimplot) + 
          geom_point(aes(x = UMAP_1 , y = UMAP_2),data = dimplot,color="grey",size = 0.2, alpha = 1) + 
          geom_point(aes(x = UMAP_1 , y = UMAP_2),
                     data = dimplot[rownames(targetObject@meta.data[which(targetObject@meta.data$label0525 == cellType),]),],
                     color='blue',size = 0.2, alpha = 1)
        #d
        dimplot = DimPlot(targetObject, label = T)
        #dimplot
        #ggsave(filename = paste0(str_replace_all(cellType,pattern = " ",replacement = "_"),"/","dimplot",".png"),width = 40,height = 20,units = "cm")
        gList = list(g)
      }else{
        gList = c(gList,list(g))
      }
      
      #dimplot=DimPlot(targetObject,label = T)
      #dimplot+g+d
      #ggsave(filename = paste0(str_replace_all(cellType,pattern = " ",replacement = "_"),"/",i,".png"),width = 40,height = 20,units = "cm")
    }
  }
  
  cellCountForcellTypeNetTmp = na.omit(cellCountForcellTypeNetTmp)
  cellCountForcellTypeNetTmp = cellCountForcellTypeNetTmp[order(cellCountForcellTypeNetTmp$cellCnt, decreasing = T),]
  #挑选细胞数排序前3/4的网络图
  #netSelect=cellCountForcellTypeNetTmp$netIdx[c(1:ceiling(nrow(cellCountForcellTypeNetTmp)*3/4))] 
  #挑选细胞数不为0的网络图
  netSelect = cellCountForcellTypeNetTmp[which(cellCountForcellTypeNetTmp$cellCnt > 0),]$netIdx
  #修改dataForCluster
  #提取细胞，网络的01矩阵，用于聚类
  dataForCluster = as.data.frame(matrix(0,
                                        nrow = nrow(targetObject@meta.data),
                                        ncol = length(netSelect)
  ))
  colnames(dataForCluster) = cellTypeNetTmp[netSelect]
  base::rownames(dataForCluster) = base::rownames(targetObject@meta.data)
  
  #提取用于聚类的数据
  for(netSelecti in netSelect){
    netSelectiData = gList[[netSelecti]][["layers"]][[2]][["data"]]
    netSelectiNet = gList[[netSelecti]][["labels"]][["title"]]
    dataForCluster[base::rownames(netSelectiData),netSelectiNet] = 1
  }
  dataForCluster[is.na(dataForCluster)] = 0
  dataForCluster = dataForCluster[rowSums(dataForCluster)>0,]
  save(dataForCluster,file = paste0("luadClusterData/cellNetwork01L/", str_replace_all(cellType,pattern = " ",replacement = "_"),".RData"))
  print(cellType)
  
}


