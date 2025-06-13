BiocParallel::register(BiocParallel::SerialParam())
library(GSEABase)
library(DESeq2) 
library(clusterProfiler)
# library(enrichplot)
library(org.Hs.eg.db)
library(stringr)
library(Seurat)


# 确定细胞类型
cellTypeSet = c("memory CD8+ T cells", "naive CD4+ T cells", "naive CD8+ T cells", "pre-B cells", "Basophils", "Neutrophils")
cellType = cellTypeSet[i]
cellClusterRes = read.csv(paste("developmentalTrajectory/pbmcClusterRes/cellClusterRes_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""),
                          check.names = F, row.names = 1)
# 删除cluster为0的细胞
if(length(which(cellClusterRes$clusterRes %in% 0)) > 0){
  cellClusterRes = cellClusterRes[-which(cellClusterRes$clusterRes %in% 0),]
}
targetCell = row.names(cellClusterRes)

# 筛选subPBMC中的细胞
subPBMC = subset(pbmc, label1 %in% cellType)
subPBMC$cellName = row.names(subPBMC@meta.data) 
subPBMC = subset(subPBMC, cellName %in% targetCell)

# 获得targetcell的count矩阵
cellCountM = as.matrix(subPBMC@assays$RNA@counts)

cellCountM = cellCountM[rowSums(cellCountM)>0,]



# DESeq2 差异分析，计算每个基因的log2FoldChange
getLOGFC = function(targetCluster, cellClusterRes, cellCountM){
  cellClusterRes$label = "no"
  cellClusterRes[which(cellClusterRes$clusterRes %in% targetCluster),]$label = "yes"
  # 降序处理
  cellClusterRes = cellClusterRes[order(cellClusterRes$label ,decreasing = T),]
  # 调整cellClusterRes行名和cellCountM列名一致
  # 判断顺序是否一致
  # all(rownames(cellClusterRes) == colnames(cellCountM))
  # 对cellCountM的列按照cellClusterRes的行顺序进行重排序，并过滤掉当初被分到pbmc亚群的细胞
  cellCountM = cellCountM[,rownames(cellClusterRes)]
  # 获得每一个细胞（样本）对应的属性信息
  group = as.factor(cellClusterRes$label)
  
  colData = data.frame(row.names = colnames(cellCountM), condition = group)
  # 使用count matrix(luad_cellCountM)和样本信息(coldata)表格, 构建DESeqDataSet(dds)对象；
  dds = DESeqDataSetFromMatrix(countData = cellCountM,
                               colData = colData,
                               design = ~condition)
  # 使用DESeq()函数进行差异分析流程
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  geoMeans = apply(counts(dds), 1, gm_mean)
  dds = estimateSizeFactors(dds, geoMeans = geoMeans)
  #dds = DESeq(dds)
  #, fitType="local"
  dds = DESeq(dds, parallel = TRUE)
  # 使用results()函数提取分析结果
  # yes vs no表示log2(yes/no)，no为对照
  res = results(dds, contrast = c("condition", "yes", "no"))
  resOrdered = res[order(res$pvalue),]
  deseq = as.data.frame(resOrdered)
  # deseq = na.omit(deseq)
  return(deseq)
}

# 使用DESeq包得到的差异分析结果  deseq_ls
# 使用GO数据库富集
getGOPath = function(targetCluster, deseq_ls){
  deseq = deseq_ls[[targetCluster]]
  # 准备进行GSEA富集分析
  # 从deseq中提取细胞名和log2FoldChange
  GSEA_df = data.frame(SYMBOL = row.names(deseq), logFC = deseq$log2FoldChange)
  # 删除logFC为NA的基因
  GSEA_df = GSEA_df[!(GSEA_df$logFC == "NA"),]
  # 转换基因ID，该步骤完成后会丢失一些基因
  GSEA_dfi = bitr(GSEA_df$SYMBOL, 
                  fromType = "SYMBOL",
                  toType =  "ENTREZID",
                  OrgDb = org.Hs.eg.db) #人数据库org.Hs.eg.db 小鼠org.Mm.eg.db
  # 合并两个数据框
  GSEA_all = merge(GSEA_df, GSEA_dfi, by = "SYMBOL", all = F)
  # 创建gsea分析的geneList（包含从大到小排列的log2FoldChange和ENTREZID信息）
  geneList = GSEA_all$logFC
  names(geneList) = GSEA_all$ENTREZID
  geneList = sort(geneList, decreasing = T)   #从大到小排序
  # 利用clusterProfiler包进行GSEA富集
  GO_kk_entrez = gseGO(geneList = geneList,
                       ont = "ALL",  # "BP"、"MF"和"CC"或"ALL"
                       OrgDb = org.Hs.eg.db,#人类org.Hs.eg.db 鼠org.Mm.eg.db
                       keyType = "ENTREZID",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2)   #实际为padj阈值可调整
  GO_kk = DOSE::setReadable(GO_kk_entrez, 
                            OrgDb = org.Hs.eg.db,
                            keyType='ENTREZID') #转化id
  return(GO_kk)
}

cellTypeSet = c("memory CD8+ T cells", "naive CD4+ T cells", "naive CD8+ T cells", "pre-B cells", "Basophils", "Neutrophils")
for(i in c(1:6)){
  # 确定细胞类型
  cellType = cellTypeSet[i]
  cellClusterRes = read.csv(paste("developmentalTrajectory/pbmcClusterRes/cellClusterRes_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""),
                            check.names = F, row.names = 1)
  # 删除cluster为0的细胞
  if(length(which(cellClusterRes$clusterRes %in% 0)) > 0){
    cellClusterRes = cellClusterRes[-which(cellClusterRes$clusterRes %in% 0),]
  }
  targetCell = row.names(cellClusterRes)
  
  # 筛选subPBMC中的细胞
  subPBMC = subset(pbmc, label1 %in% cellType)
  subPBMC$cellName = row.names(subPBMC@meta.data) 
  subPBMC = subset(subPBMC, cellName %in% targetCell)
  
  # 获得targetcell的count矩阵
  cellCountM = as.matrix(subPBMC@assays$RNA@counts)
  
  cellCountM = cellCountM[rowSums(cellCountM)>0,]
  
  
  # library("BiocParallel")
  # register(MulticoreParam(3))
  deseq_ls = list()
  for(i in unique(cellClusterRes$clusterRes)){
    deseq_ls[[i]] = getLOGFC(i, cellClusterRes, cellCountM)
    print(i)
  }
  save(deseq_ls, file = paste("GSEA/deseqPBMC/deseqRes_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".RData", sep = ""))
  
  
  gseaGO_ls = list()
  for(i in unique(cellClusterRes$clusterRes)){
    gseaGO_ls[[i]] = getGOPath(i, deseq_ls)
    print(i)
  }
  save(gseaGO_ls, file = paste("GSEA/gseaPBMC/gseaRes_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".RData", sep = ""))
  
  print(cellType)
}
# 通路筛选
cellTypeSet = c("memory CD8+ T cells", "naive CD4+ T cells", "naive CD8+ T cells", "pre-B cells", "Basophils", "Neutrophils")
for(i in 1:length(cellTypeSet)){
  cellType = cellTypeSet[i]
  cellClusterRes = read.csv(paste("developmentalTrajectory/pbmcClusterRes/cellClusterRes_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""),
                            check.names = F, row.names = 1)
  # 删除cluster为0的细胞
  if(length(which(cellClusterRes$clusterRes %in% 0)) > 0){
    cellClusterRes = cellClusterRes[-which(cellClusterRes$clusterRes %in% 0),]
  }
  gseaGO_ls = get(load(file = paste("GSEA/gseaPBMC/gseaRes_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".RData", sep = "")))
  
  # 筛选出每个亚群的上调通路——all
  all_up_pathway = data.frame(goID = NULL, ontoLogy = NULL, pathway = NULL, setSize= NULL, core_enrichment = NULL, leading_edge = NULL, clusterRes = NULL)
  for(i in 1:length(gseaGO_ls)){
    gseaGO = gseaGO_ls[[i]]
    condition = gseaGO$NES>1 & gseaGO$pvalue<0.05 & gseaGO$qvalue<0.25
    up_pathway = gseaGO[condition,]$Description
    up_GO = gseaGO[condition,]$ID
    up_ONTOLOGY = gseaGO[condition,]$ONTOLOGY
    setSize = gseaGO[condition,]$setSize
    coreEnrichment = gseaGO[condition,]$core_enrichment
    leading_edge = gseaGO[condition,]$leading_edge
    all_up_pathwayi = data.frame(goID = up_GO, ontoLogy = up_ONTOLOGY, pathway = up_pathway, setSize = setSize, core_enrichment = coreEnrichment, leading_edge = leading_edge, clusterRes = rep(i, length(up_pathway)))
    all_up_pathway = rbind(all_up_pathway, all_up_pathwayi)
  }
  all_up_pathway = all_up_pathway[which(all_up_pathway$ontoLogy %in% "BP"),]
  # 按照clusterRes排序
  all_up_pathway = all_up_pathway[order(all_up_pathway$clusterRes, decreasing = F),]
  all_up_pathway = leadingEdgeProcess(all_up_pathway)
  
  # 保存富集结果
  write.csv(all_up_pathway, paste("GSEA/gseaPBMC/allUpPathway/allUpPathway_clusterRes_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""),
            row.names = F)
  
  # 筛选出每个亚群的特异上调通路——only
  only_up_pathway = data.frame(goID = NULL, ontoLogy = NULL, pathway = NULL, setSize = NULL, core_enrichment = NULL, leading_edge = NULL, clusterRes = NULL)
  for(i in 1:length(gseaGO_ls)){
    gseaGO = gseaGO_ls[[i]]
    condition = gseaGO$NES>1 & gseaGO$pvalue<0.05 & gseaGO$qvalue<0.25
    up_pathway = gseaGO[condition,]$Description
    up_GO = gseaGO[condition,]$ID
    up_ONTOLOGY = gseaGO[condition,]$ONTOLOGY
    setSize = gseaGO[condition,]$setSize
    coreEnrichment = gseaGO[condition,]$core_enrichment
    leading_edge = gseaGO[condition,]$leading_edge
    only_up_pathwayi = data.frame(goID = up_GO, ontoLogy = up_ONTOLOGY, pathway = up_pathway, setSize = setSize, core_enrichment = coreEnrichment, leading_edge = leading_edge, clusterRes = rep(i, length(up_pathway)))
    only_up_pathway = rbind(only_up_pathway, only_up_pathwayi)
  }
  only_up_pathway = leadingEdgeProcess(only_up_pathway)
  only_up_pathway = only_up_pathway[which(only_up_pathway$tag >= 40),]
  
  only_up_pathway = only_up_pathway[which(only_up_pathway$pathway %in% names(which(table(only_up_pathway$pathway) == 1))),]
  only_up_pathway = only_up_pathway[which(only_up_pathway$ontoLogy %in% "BP"),]
  # 按照clusterRes排序
  only_up_pathway = only_up_pathway[order(only_up_pathway$clusterRes, decreasing = F),]
  
  
  #temp_up_pathway = only_up_pathway[which(only_up_pathway$tag >= 30),]
  print(cellType)
  print(table(only_up_pathway$clusterRes))
  #print(table(temp_up_pathway$clusterRes))
  # 保存富集结果
  write.csv(only_up_pathway, paste("GSEA/gseaPBMC/onlyUpPathway/onlyUpPathway_clusterRes_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""),
            row.names = F)
}

leadingEdgeProcess = function(up_pathway_df){
  up_pathway_df$tag = 0
  up_pathway_df$list = 0
  up_pathway_df$signal = 0
  for(i in 1:nrow(up_pathway_df)){
    leadingEdge = up_pathway_df$leading_edge[i]
    # 使用正则表达式提取所有数字
    numbers = gregexpr("[0-9]+", leadingEdge, perl=TRUE)
    numbers = regmatches(leadingEdge, numbers)
    numbers = unlist(numbers)
    up_pathway_df$tag[i] = as.numeric(numbers[1])
    up_pathway_df$list[i] = as.numeric(numbers[2])
    up_pathway_df$signal[i] = as.numeric(numbers[3])
  }
  return(up_pathway_df)
}


# 导出热图输入数据
all_upPathway_cluster_mat = matrix(data = 0, nrow = length(all_up_pathway$pathway), ncol = length(unique(cellClusterRes$clusterRes)))
colnames(all_upPathway_cluster_mat) = c(paste("Cluster", c(1:length(unique(cellClusterRes$clusterRes))), sep = ""))
rownames(all_upPathway_cluster_mat) = all_up_pathway$pathway
# 将矩阵元素值赋值为NES
# 先列遍历再行遍历
for(i in 1:ncol(all_upPathway_cluster_mat)){
  gseaGO = gseaGO_ls[[i]]
  condition = gseaGO$NES>1 & gseaGO$pvalue<0.05 & gseaGO$qvalue<0.25
  gseaGO = gseaGO[condition,]
  for(j in 1:nrow(all_upPathway_cluster_mat)){
    if(length(which(gseaGO$Description == rownames(all_upPathway_cluster_mat)[j])) > 0){
      all_upPathway_cluster_mat[j,i] = gseaGO$NES[which(gseaGO$Description == rownames(all_upPathway_cluster_mat)[j])]
    }else{
      next
    }
    
  }
}



only_upPathway_cluster_mat = matrix(data = 0, nrow = length(only_up_pathway$pathway), ncol = length(unique(cellClusterRes$clusterRes)))
colnames(only_upPathway_cluster_mat) = c(paste("Cluster", c(1:length(unique(cellClusterRes$clusterRes))), sep = ""))
rownames(only_upPathway_cluster_mat) = only_up_pathway$pathway
# 将矩阵元素值赋值为NES
# 先列遍历再行遍历
for(i in 1:ncol(only_upPathway_cluster_mat)){
  gseaGO = gseaGO_ls[[i]]
  condition = gseaGO$NES>1 & gseaGO$pvalue<0.05 & gseaGO$qvalue<0.25
  gseaGO = gseaGO[condition,]
  for(j in 1:nrow(only_upPathway_cluster_mat)){
    if(length(which(gseaGO$Description == rownames(only_upPathway_cluster_mat)[j])) > 0){
      only_upPathway_cluster_mat[j,i] = gseaGO$NES[which(gseaGO$Description == rownames(only_upPathway_cluster_mat)[j])]
    }else{
      next
    }
  }
}


# 筛选出每个亚群的指定上调通路——special
special_up_pathway = c()
special_up_pathway = only_up_pathway[which(only_up_pathway$pathway %in% special_up_pathway),]

# 按照clusterRes排序
special_up_pathway = special_up_pathway[order(special_up_pathway$clusterRes, decreasing = F),]

# 保存富集结果
write.csv(special_up_pathway, paste("GSEA/gseaPBMC/specialUpPathway/specialUpPathway_clusterRes_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""))

# 导出热图输入数据
special_upPathway_cluster_mat = only_upPathway_cluster_mat[special_up_pathway$pathway,]


# 绘制热图
library(ComplexHeatmap)
getHT = function(upPathway_cluster_mat, rowFont){
  # 列名注释条
  annotation_col1 = data.frame(clusterName = factor(colnames(upPathway_cluster_mat)))
  
  column_ann1 = HeatmapAnnotation(df = annotation_col1,
                                  show_legend = F, show_annotation_name = F,
                                  simple_anno_size = unit(2,"mm"))
  # 行名注释条
  annotation_row1 = data.frame(pathway = factor(rownames(upPathway_cluster_mat)), fontSize = 5)
  rowmark_ann1 = HeatmapAnnotation(df = annotation_row1,
                                   show_legend = F, show_annotation_name = T,
                                   which = "row")
  ht = Heatmap(upPathway_cluster_mat,
               name = "NES", 
               #col = col_fun, 
               cluster_rows = T, cluster_columns = F,
               show_row_dend = F, show_column_dend = F,
               show_row_names = T, show_column_names = T,
               row_names_side = "left",
               row_names_gp = gpar(fontsize = rowFont),
               column_names_side = "top", column_names_rot = 45,
               column_names_gp = gpar(fontsize = 10),
               top_annotation = column_ann1,
               width = unit(8, "cm"),
               height = unit(16, "cm"))
}
draw(getHT(all_upPathway_cluster_mat, 1))
draw(getHT(only_upPathway_cluster_mat, 1))
draw(getHT(special_upPathway_cluster_mat, 5))