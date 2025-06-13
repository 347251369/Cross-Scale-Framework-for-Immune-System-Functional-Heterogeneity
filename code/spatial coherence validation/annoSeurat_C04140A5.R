library(Seurat)
sampleID = "C04140A5"
fileDirTarget = "spatialDataProcess/stTargetData_tissue/"
stereoSeurat_all = readRDS(file = paste(fileDirTarget, sampleID, "_stereoSeurat.rds", sep = ""))
# 注释immune和noImmune
stereoSeurat_all
stereoSeurat = stereoSeurat_all[["bin200"]]
stereoSeurat = NormalizeData(stereoSeurat, normalization.method = "LogNormalize", scale.factor = 10000)
stereoSeurat = FindVariableFeatures(stereoSeurat, selection.method = "vst", nfeatures = 2000)
all.genes = rownames(stereoSeurat)
stereoSeurat = ScaleData(stereoSeurat, features = all.genes)
stereoSeurat = RunPCA(stereoSeurat, features = VariableFeatures(object = stereoSeurat))
ElbowPlot(stereoSeurat)
stereoSeurat = FindNeighbors(stereoSeurat, dims = 1:15)
stereoSeurat = FindClusters(stereoSeurat, resolution = 5)
stereoSeurat = RunUMAP(stereoSeurat, dims = 1:15)
# 导入来自stereopy的cluster结果
classificationRes = read.csv(file = paste(fileDirTarget, sampleID, "_classificationRes_bin200.csv", sep = ""))
stereoSeurat$classificationRes = classificationRes$group
DimPlot(stereoSeurat, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "classificationRes")
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
gene.markers = unique(c(habermann_epi, habermann_imm, habermann_oth))
gene.features = gene.markers[gene.markers %in% row.names(stereoSeurat)]
stereoSeurat.markers = FindAllMarkers(stereoSeurat, only.pos = TRUE,
                                      features = gene.features)
stereoSeurat.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.25)

habermann_epi = c("ABCA3", "SFTPB", "SFTPC", "AGER", "PDPN",  "KRT5", "TRP63", "NGFR", "SCGB1A1", "MUC5B", "KRT17", "FOXJ1", "TMEM190", "CAPS", "CHGA", "CALCA", "ASCL1", "PTPRC", "EPCAM")
habermann_imm = c("CD3E", "CD4", "FOXP3", "IL7R", "IL2RA", "CD40LG", "CD8A", "CCL5", "NCR1", "KLRB1", "NKG7", "LYZ", "CD68", "ITGAX", "MARCO", "FCGR1A", "FCGR3A", "C1QA", "APOC1", "S100A12", "FCN1", "S100A9", "CD14", "FCER1A", "CD1C", "CD16", "CLEC9A", "LILRA4", "CLEC4C", "JCHAIN", "IGHG1", "IGLL5", "MS4A1", "CD19", "CD79A", "CPA3", "KIT", "MKI67", "CDK1", "EPCAM")
habermann_oth = c("VWF", "PECAM1", "CCL21", "PROX1", "ACTA2", "MYH11", "PDGFRB", "WT1", "UPK3B", "LUM", "PDGFRA", "MYLK", "HAS1", "PLIN2", "FAP", "PTPRC", "EPCAM")
immune.markers = c("CD3E", "CD4", "IL7R", "IL2RA", "LYZ", "CD68", "ITGAX", "MARCO", "FCGR3A", "C1QA", "JCHAIN", "IGHG1", "IGLL5", "MS4A1", "PTPRC", "CD14")
noImmune.markers = c("ABCA3", "SFTPB", "SFTPC", "PDPN",  "KRT5", "SCGB1A1", "MUC5B", "EPCAM",
                     "VWF", "PECAM1", "LUM", "DCN", "COL1A1", "COL1A2")
immune.markers = c("CD3E", "CD4", "CD3D", "MS4A1", "CD79A", "JCHAIN", "IGHG1", "LYZ", "CD68", "MARCO", "FCGR3A", "CD14")
noImmune.markers = c("EPCAM", "CDH1", "KRT19", "PECAM1", "VWF",
                     "DCN", "COL1A1", "COL1A2")
DotPlot(stereoSeurat, features = c(noImmune.markers, immune.markers), group.by = "cellType")
new.cluster.ids = c("noImmune", "noImmune", "noImmune", "immune", "noImmune",
                    "immune", "noImmune", "noImmune", "noImmune", "noImmune",
                    "noImmune", "noImmune", "immune", "noImmune", "noImmune",
                    "noImmune", "noImmune", "noImmune", "noImmune", "noImmune",
                    "noImmune", "immune", "noImmune", "immune", "immune",
                    "noImmune", "noImmune", "immune", "noImmune", "noImmune",
                    "noImmune", "immune", "noImmune", "noImmune", "noImmune",
                    "immune", "noImmune", "immune", "immune", "noImmune",
                    "noImmune", "immune", "noImmune", "noImmune", "immune",
                    "noImmune", "noImmune", "noImmune", "noImmune", "noImmune",
                    "immune", "noImmune")
new.cluster.ids = c("noImmune", "noImmune", "immune", "immune", "immune",
                    "immune", "noImmune", "noImmune", "noImmune", "immune",
                    "noImmune", "noImmune", "immune", "immune", "immune",
                    "immune", "immune", "noImmune", "noImmune", "noImmune",
                    "noImmune", "immune", "immune", "immune", "immune",
                    "noImmune", "noImmune", "immune", "noImmune", "noImmune",
                    "noImmune", "immune", "immune", "immune", "noImmune",
                    "immune", "immune", "immune", "immune", "noImmune",
                    "immune", "immune", "immune", "noImmune", "immune",
                    "immune", "noImmune", "noImmune", "noImmune", "noImmune",
                    "immune", "noImmune") # end
new.cluster.ids = c("noImmune", "noImmune", "immune", "immune", "immune",
                    "noImmune", "noImmune", "immune", "immune", "immune",
                    "noImmune", "noImmune", "noImmune", "noImmune", "noImmune",
                    "noImmune", "noImmune", "noImmune", "immune", "noImmune",
                    "noImmune", "immune", "immune", "immune", "immune",
                    "immune", "immune", "noImmune", "immune", "noImmune",
                    "noImmune", "immune", "immune", "noImmune", "noImmune",
                    "noImmune", "noImmune", "noImmune", "immune", "noImmune",
                    "noImmune", "noImmune", "immune", "noImmune", "immune",
                    "noImmune", "noImmune", "immune", "noImmune", "noImmune",
                    "noImmune", "immune", "noImmune", "noImmune", "noImmune",
                    "noImmune", "noImmune", "immune", "noImmune", "immune",
                    "noImmune", "noImmune", "noImmune", "noImmune", "immune",
                    "noImmune", "noImmune")
names(new.cluster.ids) = levels(stereoSeurat)
stereoSeurat = RenameIdents(stereoSeurat, new.cluster.ids)
stereoSeurat$cellType = Idents(stereoSeurat)
DimPlot(stereoSeurat, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
immune.markers = c("CD3E", "CD4", "CD3D", "MS4A1", "CD79A", "JCHAIN", "IGHG1", "LYZ", "CD68", "MARCO", "FCGR3A", "CD14", "S100A8", "S100A9")
noImmune.markers = c("EPCAM", "CDH1", "KRT19", "VWF",
                     "DCN", "COL1A1", "COL1A2")
DotPlot(stereoSeurat, features = c(noImmune.markers, immune.markers))
# stereoSeurat.immune----------------------------------------------------
stereoSeurat.immune = subset(stereoSeurat, cellType == "immune")
stereoSeurat.immune
stereoSeurat.immune = NormalizeData(stereoSeurat.immune, normalization.method = "LogNormalize", scale.factor = 10000)
stereoSeurat.immune = FindVariableFeatures(stereoSeurat.immune, selection.method = "vst", nfeatures = 2000)
all.genes = rownames(stereoSeurat.immune)
stereoSeurat.immune = ScaleData(stereoSeurat.immune, features = all.genes)
stereoSeurat.immune = RunPCA(stereoSeurat.immune, features = VariableFeatures(object = stereoSeurat.immune))
ElbowPlot(stereoSeurat.immune)
stereoSeurat.immune = FindNeighbors(stereoSeurat.immune, dims = 1:15)
stereoSeurat.immune = FindClusters(stereoSeurat.immune, resolution = 5)
stereoSeurat.immune = RunUMAP(stereoSeurat.immune, dims = 1:15)
DimPlot(stereoSeurat.immune, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "seurat_clusters")
immune.markers = c("IL7R", "PTPRC", "CD3E", "CD4", "CD3D", "NKG7", "GZMA", "MS4A1", "CD79A", "JCHAIN", "IGHG1", "LYZ", "CD68", "MARCO", "FCGR3A", "CD14")
DotPlot(stereoSeurat.immune, features = immune.markers, group.by = "seurat_clusters")
new.cluster.ids = c("U", "B", "T", "B", "M",
                    "B", "B", "B", "T", "M",
                    "T", "T", "T", "M", "M",
                    "T", "B", "M", "T", "B",
                    "B", "B", "B", "B", "M",
                    "U", "M", "M", "B", "B")
new.cluster.ids = c("B", "B", "U", "U", "B",
                    "B", "B", "B", "U", "T",
                    "B", "M", "B", "T", "T",
                    "T", "B", "T", "B", "B",
                    "B", "M", "T", "B", "B",
                    "T", "T", "T", "M", "M",
                    "M", "T", "T", "T", "B",
                    "M", "T", "T", "T", "U",
                    "B")
new.cluster.ids = c("B", "B", "U", "U", "B",
                    "B", "B", "B", "U", "T",
                    "B", "M", "B", "T", "T",
                    "T", "B", "T", "B", "B",
                    "B", "M", "T", "B", "B",
                    "T", "T", "T", "M", "M",
                    "M", "T", "T", "T", "B",
                    "M", "T", "T", "T", "U",
                    "B") # 0905
names(new.cluster.ids) = levels(stereoSeurat.immune)
stereoSeurat.immune = RenameIdents(stereoSeurat.immune, new.cluster.ids)
stereoSeurat.immune$cellType_immune = Idents(stereoSeurat.immune)
DimPlot(stereoSeurat.immune, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "cellType_immune")
immune.markers = c("IL7R", "CD3E", "CD4", "CD3D", "NKG7", "GZMA", "MS4A1", "CD79A", "JCHAIN", "IGHG1", "LYZ", "CD68", "MARCO", "FCGR3A", "CD14")
DotPlot(stereoSeurat.immune, features = immune.markers, group.by = "cellType_immune")

stereoSeurat.B = subset(stereoSeurat.immune, cellType_immune == "B")
stereoSeurat.T = subset(stereoSeurat.immune, cellType_immune == "T")
stereoSeurat.M = subset(stereoSeurat.immune, cellType_immune == "M")
# stereoSeurat.B----------------------------------------------------
stereoSeurat.B = NormalizeData(stereoSeurat.B, normalization.method = "LogNormalize", scale.factor = 10000)
stereoSeurat.B = FindVariableFeatures(stereoSeurat.B, selection.method = "vst", nfeatures = 2000)
all.genes = rownames(stereoSeurat.B)
stereoSeurat.B = ScaleData(stereoSeurat.B, features = all.genes)
stereoSeurat.B = RunPCA(stereoSeurat.B, features = VariableFeatures(object = stereoSeurat.B))
ElbowPlot(stereoSeurat.B)
stereoSeurat.B = FindNeighbors(stereoSeurat.B, dims = 1:10)
stereoSeurat.B = FindClusters(stereoSeurat.B, resolution = 5)
stereoSeurat.B = RunUMAP(stereoSeurat.B, dims = 1:10)
DimPlot(stereoSeurat.B, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "seurat_clusters")
B.markers = c("MS4A1", "CD79A", "JCHAIN", "IGHG1")
DotPlot(stereoSeurat.B, features = B.markers, group.by = "seurat_clusters")
new.cluster.ids = c("Plasma cells", "Plasma cells", "Plasma cells", "Plasma cells", "Plasma cells",
                    "pre-B cells", "Plasma cells", "Plasma cells", "Plasma cells", "Plasma cells",
                    "Plasma cells", "Plasma cells", "Plasma cells", "Plasma cells", "pre-B cells",
                    "Plasma cells", "Plasma cells", "Plasma cells", "Plasma cells", "Plasma cells")

new.cluster.ids = c("Plasma cells", "Plasma cells", "pre-B cells", "Plasma cells", "Plasma cells",
                    "Plasma cells", "Plasma cells", "Plasma cells", "Plasma cells", "Plasma cells",
                    "Plasma cells", "Plasma cells", "pre-B cells", "Plasma cells", "Plasma cells",
                    "Plasma cells", "Plasma cells", "Plasma cells", "Plasma cells", "Plasma cells",
                    "Plasma cells", "Plasma cells", "Plasma cells", "Plasma cells", "Plasma cells",
                    "Plasma cells", "Plasma cells", "Plasma cells", "Plasma cells", "Plasma cells") # cluster26 pre-B换plasma
names(new.cluster.ids) = levels(stereoSeurat.B)
stereoSeurat.B = RenameIdents(stereoSeurat.B, new.cluster.ids)
stereoSeurat.B$cellType_subset = Idents(stereoSeurat.B)
DimPlot(stereoSeurat.B, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "cellType_subset") + NoLegend()
B.markers = c("MS4A1", "CD79A", "JCHAIN", "IGHG1")
DotPlot(stereoSeurat.B, features = B.markers, group.by = "cellType_subset")
# stereoSeurat.M----------------------------------------
stereoSeurat.M = NormalizeData(stereoSeurat.M, normalization.method = "LogNormalize", scale.factor = 10000)
stereoSeurat.M = FindVariableFeatures(stereoSeurat.M, selection.method = "vst", nfeatures = 2000)
all.genes = rownames(stereoSeurat.M)
stereoSeurat.M = ScaleData(stereoSeurat.M, features = all.genes)
stereoSeurat.M = RunPCA(stereoSeurat.M, features = VariableFeatures(object = stereoSeurat.M))
ElbowPlot(stereoSeurat.M)
stereoSeurat.M = FindNeighbors(stereoSeurat.M, dims = 1:10)
stereoSeurat.M = FindClusters(stereoSeurat.M, resolution = 5)
stereoSeurat.M = RunUMAP(stereoSeurat.M, dims = 1:10)
DimPlot(stereoSeurat.M, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "seurat_clusters")
M.markers = c("SPP1", "LYZ", "CD1C", "S100A8", "S100A9", "CPA3", "MS4A2", "C1QA", "CD74", "CD1A", "XCR1", "CD83", "CD68", "CD14",
              "CD195", "TPSAB1", "BSP1")
M.markers = c("C1QA", "CD14", "CD68", "LYZ", "CD74", "CD1C", "S100A8", "S100A9", "CPA3", "MS4A2")
DotPlot(stereoSeurat.M, features = M.markers, group.by = "seurat_clusters")
new.cluster.ids = c("Basophils", "Macrophages", "cDC", "Macrophages", "Monocytes",
                    "cDC", "Macrophages", "Macrophages", "Macrophages", "Macrophages",
                    "Macrophages", "Monocytes", "cDC", "cDC", "Macrophages",
                    "Macrophages", "Macrophages", "Monocytes")
new.cluster.ids = c("Basophils", "Macrophages", "cDC", "Macrophages", "Monocytes",
                    "Macrophages", "Macrophages", "Macrophages", "Macrophages", "Macrophages",
                    "Macrophages", "Monocytes", "cDC", "Macrophages", "cDC",
                    "Macrophages", "Macrophages", "Monocytes")
new.cluster.ids = c("cDC", "Macrophages", "cDC", "Macrophages", "Monocytes",
                    "Macrophages", "Macrophages", "Macrophages", "Macrophages", "Macrophages",
                    "Macrophages", "Monocytes", "cDC", "Macrophages", "cDC",
                    "Macrophages", "Macrophages", "Monocytes")
new.cluster.ids = c("Macrophages", "Macrophages", "Macrophages", "Macrophages", "Monocytes",
                    "Basophils", "cDC", "cDC", "cDC", "Macrophages",
                    "Monocytes", "Monocytes", "Macrophages", "Monocytes", "Macrophages",
                    "Macrophages", "Macrophages", "Monocytes") # 0905
names(new.cluster.ids) = levels(stereoSeurat.M)
stereoSeurat.M = RenameIdents(stereoSeurat.M, new.cluster.ids)
stereoSeurat.M$cellType_subset = Idents(stereoSeurat.M)
DimPlot(stereoSeurat.M, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "cellType_subset")
M.markers = c("CD68", "LYZ", "CD74", "CPA3", "MS4A2", "CD83", "CD14")
DotPlot(stereoSeurat.M, features = M.markers, group.by = "cellType_subset")

# stereoSeurat.T-----------------------------------------------
stereoSeurat.T = NormalizeData(stereoSeurat.T, normalization.method = "LogNormalize", scale.factor = 10000)
stereoSeurat.T = FindVariableFeatures(stereoSeurat.T, selection.method = "vst", nfeatures = 2000)
all.genes = rownames(stereoSeurat.T)
stereoSeurat.T = ScaleData(stereoSeurat.T, features = all.genes)
stereoSeurat.T = RunPCA(stereoSeurat.T, features = VariableFeatures(object = stereoSeurat.T))
ElbowPlot(stereoSeurat.T)
stereoSeurat.T = FindNeighbors(stereoSeurat.T, dims = 1:10) # 1:5
stereoSeurat.T = FindClusters(stereoSeurat.T, resolution = 5)
stereoSeurat.T = RunUMAP(stereoSeurat.T, dims = 1:10)
DimPlot(stereoSeurat.T, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "seurat_clusters")
T.markers = c("CCR7", "IL7R", "CD4", "LEF1", "SELL", "GZMA", "FCER1G", "TYROBP", "CTLA4", "IL2RA", "CD3E")
DotPlot(stereoSeurat.T, features = T.markers, group.by = "seurat_clusters")
new.cluster.ids = c("memory CD4+ T cells", "NKT", "NKT", "naive CD4+ T cells", "naive CD4+ T cells",
                    "naive CD4+ T cells", "naive CD4+ T cells", "memory CD4+ T cells", "NKT", "NKT",
                    "U", "memory CD4+ T cells", "naive CD4+ T cells", "naive CD4+ T cells", "naive CD4+ T cells",
                    "Tregs", "Tregs", "memory CD4+ T cells", "memory CD8+ T cells", "Tregs")
new.cluster.ids = c("U", "U", "Tregs", "U", "memory CD4+ T cells",
                    "NKT", "naive CD4+ T cells", "U", "memory CD4+ T cells", "U",
                    "U", "NKT", "naive CD4+ T cells", "naive CD8+ T cells", "U",
                    "naive CD4+ T cells", "naive CD4+ T cells", "naive CD4+ T cells", "naive CD8+ T cells", "U",
                    "naive CD8+ T cells", "naive CD4+ T cells", "memory CD8+ T cells")

new.cluster.ids = c("memory CD4+ T cells", "memory CD4+ T cells", "Tregs", "U", "memory CD4+ T cells",
                    "NKT", "naive CD4+ T cells", "U", "naive CD4+ T cells", "U",
                    "U", "NKT", "naive CD4+ T cells", "naive CD8+ T cells", "U",
                    "naive CD4+ T cells", "naive CD4+ T cells", "naive CD4+ T cells", "naive CD8+ T cells", "U",
                    "naive CD8+ T cells", "naive CD4+ T cells", "memory CD8+ T cells") # 0902
new.cluster.ids = c("U", "U", "memory CD4+ T cells", "U", "memory CD4+ T cells",
                    "NKT", "naive CD4+ T cells", "U", "naive CD4+ T cells", "U",
                    "U", "NKT", "naive CD4+ T cells", "naive CD8+ T cells", "U",
                    "naive CD4+ T cells", "U", "naive CD4+ T cells", "naive CD8+ T cells", "U",
                    "naive CD8+ T cells", "naive CD4+ T cells", "memory CD8+ T cells") # 0902
new.cluster.ids = c("U", "U", "naive CD8+ T cells", "U", "naive CD8+ T cells",
                    "NKT", "naive CD4+ T cells", "U", "naive CD4+ T cells", "naive CD4+ T cells",
                    "naive CD8+ T cells", "NKT", "U", "memory CD4+ T cells", "U",
                    "naive CD4+ T cells", "memory CD4+ T cells", "naive CD4+ T cells", "naive CD4+ T cells", "U",
                    "memory CD4+ T cells", "naive CD4+ T cells", "memory CD8+ T cells") # 0905
names(new.cluster.ids) = levels(stereoSeurat.T)
stereoSeurat.T = RenameIdents(stereoSeurat.T, new.cluster.ids)
stereoSeurat.T$cellType_subset = Idents(stereoSeurat.T)
DimPlot(stereoSeurat.T, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "cellType_subset")
T.markers = c("CCR7", "IL7R", "CD4", "LEF1", "SELL", "GZMA", "FCER1G", "TYROBP", "CTLA4", "IL2RA")
DotPlot(stereoSeurat.T, features = T.markers)

# 查看每个cluster的networkBiomarker表达情况-------------------------------------
stereoSeurat_tmp = stereoSeurat.T
count = LayerData(stereoSeurat_tmp, assay = "RNA", layer = "counts")
count = t(count)
targetBioGene = bioGene[which(bioGene %in% colnames(count))] # bioGene中的空转基因
count = count[, targetBioGene]
count = as.data.frame(count)
# 输入bin-bioGene矩阵（count.M），输出bin-networkBiomarker矩阵
binNetworkBiomarker = matrix(data = 0, nrow = nrow(count), ncol = 75) # 75表示networkBiomarker数目
rownames(binNetworkBiomarker) = rownames(count)
colnames(binNetworkBiomarker) = networkBiomarkerSummary$networkBiomarker[1:75]
for(i in 1:nrow(count)){
  bin = rownames(count)[i]
  binExpr = count[i,]
  for(j in 1:ncol(binNetworkBiomarker)){
    networkBiomarker = colnames(binNetworkBiomarker)[j]
    networkBiomarker_gene = unlist(str_split(networkBiomarker, "_"))
    if(all(networkBiomarker_gene %in% colnames(binExpr))){
      if(all(as.numeric(binExpr[,networkBiomarker_gene]))){
        binNetworkBiomarker[i,j] = 1
      }
    }
  }
  print(nrow(count) - i)
}
binNetworkBiomarker = as.data.frame(binNetworkBiomarker)
binNetworkBiomarker$cluster = as.numeric(stereoSeurat_tmp$seurat_clusters)
binNetworkBiomarker = as.matrix(binNetworkBiomarker)
# 利用networkBiomarker表达情况反向注释
expNetworkBiomarker = data.frame(cellType = networkBiomarkerSummary$cellType[1:75], networkBiomarker = networkBiomarkerSummary$networkBiomarker[1:75])
for(cluster in 1:max(as.numeric(stereoSeurat_tmp$seurat_clusters))){
  expNetworkBiomarker[,cluster+2] = colSums(binNetworkBiomarker[binNetworkBiomarker[,ncol(binNetworkBiomarker)] == cluster,])[1:75]
}
colnames(expNetworkBiomarker)[3:ncol(expNetworkBiomarker)] = paste("cluster", 1:max(as.numeric(stereoSeurat_tmp$seurat_clusters)), sep = "")
write.csv(expNetworkBiomarker, file = paste("spatialDataProcess/expNetworkBiomarker_", sampleID, "_T_bin200.csv", sep = ""))

# metaBiomarker
metaBiomarkerInfo = stereoSeurat.T@meta.data[,c('x','y', 'cellType_subset', 'cellType', 'cellType_immune', 'seurat_clusters', 'cellType_biomarker')]
metaBiomarkerInfo$biomarker = "unknown"
biomarker = "A2M_EOMES_PTPRC"
metaBiomarkerInfo$biomarker[binNetworkBiomarker[,biomarker] > 0] = biomarker
metaBiomarkerInfo = metaBiomarkerInfo[metaBiomarkerInfo$biomarker != "unknown",]


# 整合stereoSeurat.immune、stereoSeurat.B、stereoSeurat.M和stereoSeurat中的注释结果----------------
# stereoSeurat中的cellType_immune
stereoSeurat$cellType_immune = "noImmune"
stereoSeurat$cellType_immune[match(row.names(stereoSeurat.immune@meta.data), row.names(stereoSeurat@meta.data))] = as.character(stereoSeurat.immune$cellType_immune)
# stereoSeurat中的cellType_subset
stereoSeurat$cellType_subset = stereoSeurat$cellType_immune
stereoSeurat$cellType_subset[match(row.names(stereoSeurat.B@meta.data), row.names(stereoSeurat@meta.data))] = as.character(stereoSeurat.B$cellType_subset)
stereoSeurat$cellType_subset[match(row.names(stereoSeurat.M@meta.data), row.names(stereoSeurat@meta.data))] = as.character(stereoSeurat.M$cellType_subset)
stereoSeurat$cellType_subset[match(row.names(stereoSeurat.T@meta.data), row.names(stereoSeurat@meta.data))] = as.character(stereoSeurat.T$cellType_subset)

# 保存注释结果
saveRDS(stereoSeurat, file = paste(fileDirTarget, sampleID, "_stereoSeurat_anno_avg.rds", sep = ""))
