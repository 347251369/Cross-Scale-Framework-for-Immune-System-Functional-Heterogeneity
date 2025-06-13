
# PBMC -------------------------------------------------
# 获得seuratObject
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

# 重新处理Seurat Object，覆盖原始数据
# 归一化
subPBMC = NormalizeData(subPBMC, normalization.method = "LogNormalize", scale.factor = 10000)
# 鉴定高度变化网络
subPBMC = FindVariableFeatures(subPBMC, selection.method = "vst", nfeatures = 2000)
# 中心化，为下一步PCA做准备
all_gene = rownames(subPBMC)
subPBMC = ScaleData(subPBMC, features = all_gene)

subPBMC$cellClusterRes = as.factor(paste("cluster", cellClusterRes$clusterRes, sep = ""))


# 利用通路确定featureGene
# 加载通路信息
goInfo = read.csv(paste("GSEA/goAnalysisFinal/PBMC/goResFinal_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""),
                  check.names = F, row.names = 1)

featureGene_ls = list()
for (i in unique(goInfo$cluster)) {
  goInfoi = goInfo[goInfo$cluster == i,]
  # interested通路数
  goNum = nrow(goInfoi)
  allCoreGene = unlist(str_split(goInfoi$geneID, "/"))
  # 找出通路共有的coreGene作为featureGene
  statistic = as.data.frame(table(allCoreGene))
  featureGene = as.character(statistic[statistic$Freq == goNum,]$allCoreGene)
  featureGene_ls[[i]] = featureGene
}
allFeatureGene = unlist(featureGene_ls)
allFeatureGene = unique(allFeatureGene)
allFeatureGene = c("IFIT1", "ISG15", "LYN", "CXCL8", "EGR3", "FFAR2")
# 通过点图可视化
DotPlot(object = subPBMC, features = allFeatureGene, group.by = "cellClusterRes")


# LUAD ---------------------------------------------------------------
# 获得seuratObject
# 确定细胞类型
cellTypeSet = c("memory CD8+ T cells", "naive CD4+ T cells", "Tregs", "Basophils")
cellType = cellTypeSet[i]
cellClusterRes = read.csv(paste("developmentalTrajectory/luadClusterRes/cellClusterRes_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""),
                          check.names = F, row.names = 1)
# 删除cluster为0的细胞
if(length(which(cellClusterRes$clusterRes %in% 0)) > 0){
  cellClusterRes = cellClusterRes[-which(cellClusterRes$clusterRes %in% 0),]
}
targetCell = row.names(cellClusterRes)

# 筛选subLUAD中的细胞
subLUAD = subset(luad, label0525 %in% cellType)
subLUAD$cellName = row.names(subLUAD@meta.data)
subLUAD = subset(subLUAD, cellName %in% targetCell)

# 重新处理Seurat Object，覆盖原始数据
# 归一化
subLUAD = NormalizeData(subLUAD, normalization.method = "LogNormalize", scale.factor = 10000)
# 鉴定高度变化网络
subLUAD = FindVariableFeatures(subLUAD, selection.method = "vst", nfeatures = 2000)
# 中心化，为下一步PCA做准备
all_gene = rownames(subLUAD)
subLUAD = ScaleData(subLUAD, features = all_gene)

subLUAD$cellClusterRes = as.factor(paste("细胞簇", cellClusterRes$clusterRes, sep = ""))


# 利用通路确定featureGene
# 加载通路信息
goInfo = read.csv(paste("GSEA/goAnalysisFinal/LUAD/goResFinal_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".csv", sep = ""),
                  check.names = F, row.names = 1)

featureGene_ls = list()
for (i in unique(goInfo$cluster)) {
  goInfoi = goInfo[goInfo$cluster == i,]
  # interested通路数
  goNum = nrow(goInfoi)
  allCoreGene = unlist(str_split(goInfoi$geneID, "/"))
  # 找出通路共有的coreGene作为featureGene
  statistic = as.data.frame(table(allCoreGene))
  featureGene = as.character(statistic[statistic$Freq == goNum,]$allCoreGene)
  featureGene_ls[[i]] = featureGene
}
allFeatureGene = unlist(featureGene_ls)
allFeatureGene = unique(allFeatureGene)
allFeatureGene = allFeatureGene[!(allFeatureGene %in% c("ICOS", "IL2RG", "CCR7", "PTK2", "CD55"))]
allFeatureGene = c("CTSG", "LGALS3", "BTK", "LAPTM5")
allFeatureGene = c("FOXP3", "GATA3", "NR3C1", "ZEB2")
# 通过点图可视化
DotPlot(object = subLUAD, features = allFeatureGene, group.by = "cellClusterRes", cex = 1.2)

# 绘制差异基因火山图
deseq = get(load(paste("GSEA/deseqLUAD/deseqRes_", str_replace_all(cellType,pattern = " ",replacement = "_"), ".RData", sep = "")))
resData = deseq[[2]]
resData = na.omit(resData)
# 数据分类
# 将padj值＜0.05，FDR>=2设定UP;padj值＜0.05,FDR<=-2设定为Down；其他为NS
resData$group <- as.factor(ifelse(resData$padj<0.05&abs(resData$log2FoldChange)>=1,ifelse(resData$log2FoldChange>=1,"Up","Down"),"NS"))
#查看上调和下调基因数目
table(resData$group)
# 将padj换算成-log10（padj）
resData$padj[resData$padj == 0] = .Machine$double.eps
resData$logadj = -log10(resData$padj)
#绘图，以及 label 列，您想要用标签来标识点
p <- ggscatter(resData, x = "log2FoldChange", y = "logadj",
               color = "group",
               palette = c("#2f5688", "#BBBBBB", "#CC0000"),
               size = 2,
               font.label = 8,
               repel = TRUE,
               xlab = "log2FoldChange",
               ylab = "-log10(padj)")
#添加分界线,foldchange为2，padj为0.05（-log10（0.05）为1.3），-log10（0.01）为2
p <- p +  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = 1.30, linetype = "dashed")
# 设置图例位置，添加主题
p <- p + theme(legend.position = "right") +  theme_classic()
# 打印图形
print(p) 


