
library(monocle)
# LUAD
cellTypeSet = c("memory CD8+ T cells", "naive CD4+ T cells", "Tregs", "Basophils")
cellType = cellTypeSet[i]
cellClusterRes = read.csv(paste("developmentalTrajectory/luadClusterRes/cellClusterRes_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""),
                          check.names = F, row.names = 1)
# 删除cluster为0的细胞
if(length(which(cellClusterRes$clusterRes %in% 0)) > 0){
  cellClusterRes = cellClusterRes[-which(cellClusterRes$clusterRes %in% 0),]
}
targetCell = row.names(cellClusterRes)

subLUAD = subset(luad, label0525 %in% cellType)
# 筛选subLUAD中的细胞
subLUAD$cellName = row.names(subLUAD@meta.data) 
subLUAD = subset(subLUAD, cellName %in% targetCell)
# 拟时序分析
##提取表型信息--细胞信息(建议载入细胞的聚类或者细胞类型鉴定信息、实验条件等信息)
expr_matrix = as.sparse(subLUAD@assays$RNA@counts)
##提取表型信息到p_ data(phenotype_data) 里面
p_data = subLUAD@meta.data
# p_data$celltype = subLUAD@meta.data ##整合每个细胞的细胞鉴定信息到p_data里面。如果已经)
#提取基因信息如生物类型、gc含量等
f_data = data.frame(gene_short_name = row.names(subLUAD) , row.names = row.names(subLUAD))
##expr_matrix的行数与f_data的行数相同(gene number), expr.matrix的列数与p_data的行数相同(cell number)
#构建CDS对象
pd = new('AnnotatedDataFrame', data = p_data)
fd = new('AnnotatedDataFrame', data = f_data)
#将p_data和f_data从data.frame转换AnnotatedDataFrame对象。
cds = newCellDataSet(expr_matrix,
                     phenoData = pd,
                     featureData = fd,
                     lowerDetectionLimit = 0.5,
                     expressionFamily = negbinomial.size())

# 标准化
cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)

##使用monocle选择的高变基因
disp_table = dispersionTable(cds)
disp.genes = subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
cds = setOrderingFilter(cds, row.names(subLUAD))
plot_ordering_genes(cds)
#出的图黑色的点表示用来构建轨迹的差异基因，灰色表示背景基因
#降维
cds = reduceDimension(cds, max_components = 2, method = "DDRTree")
#拟时间轴轨迹构建和在拟时间内排列细胞
cds = orderCells(cds)

pData(cds)$clusterRes = as.character(pData(cds)$clusterRes)

plot_cell_trajectory(cds, color_by = "clusterRes",  size = 1, show_bakbone = TRUE)+
  facet_wrap(~clusterRes, nrow = 3)