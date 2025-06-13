# 细胞注释结果主要是在华西生信分析人员给出的注释结果基础上进行微调
# 具体的代码略过
# 数据说明
# rawScDataPBMC:原始PBMC单细胞数据，共计41825个细胞
# scDataPBMC:注释好的PBMC单细胞数据，mate.data中，major1列中T cell没有被划分，label1列中包含细分后的T cell亚群
# scDataPBMC_T:注释好的PBMC T细胞数据，mate.data中，label1列中包含细分后的T cell亚群，不含unknown细胞
# 重点关注label1

# LUAD
# rawScDataLUAD:无单独存在的原始LUAD数据
# scDataLUAD:注释好的LUAD单细胞数据，mate.data中，major列中T cell没有被划分，
# label0521列中包含细分后的NK cell亚群和不同的myeloid亚群划分，label0525列中NK cell亚群未划分
# scDataLUAD_T:注释好的LUAD T细胞数据，mate.data中，label0521不含Th和Te亚群，label0525列中包含更细的T cell亚群
# 重点关注label0525


# PBMC------------------------------------------------------------
# 可视化细胞注释结果
cols = c("#00CED1",
         "#FFC0CB", "#8B0000",
         "#9ACD32", "#32CD32", "#9932CC",
         "#00BFFF", "#6A5ACD",
         "#FF1493") # 指定每个组别的颜色
names(cols) = c("NK",
                "Plasma cells","pre-B cells",
                "cDC","pDC","Monocytes",
                "Basophils","Neutrophils",
                "T cells")
DimPlot(scDataPBMC, 
        group.by = "major1",
        cols = cols, # 指定每个组别的颜色
        label = TRUE,
        cells = rownames(scDataPBMC@meta.data)[which(scDataPBMC$major1 %in% c("NK",
                                                                              "Plasma cells","pre-B cells",
                                                                              "cDC","pDC","Monocytes",
                                                                              "Basophils","Neutrophils",
                                                                              "T cells"))])

# 保存
saveRDS(scDataPBMC, file = "cellAnnotation/scDataPBMC.rds")

scDataPBMC = readRDS(file = "cellAnnotation/scDataPBMC.rds")

# PBMC T细胞亚群重新降维和可视化
scDataPBMC_T = subset(scDataPBMC, label1 %in% c("memory CD4+ T cells","memory CD8+ T cells","naive CD4+ T cells",
                                                "naive CD8+ T cells", "MAIT", "Tregs", "NKT"))

sce = scDataPBMC_T
# 归一化、鉴定高度变化基因、中心化，为下一步PCA做准备
sce = NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
sce = FindVariableFeatures(sce, selection.method = "vst", nfeatures = 4000)
sce = ScaleData(sce)
# PCA、判断最终选取的主成分数
sce = RunPCA(sce, features = VariableFeatures(object = sce))
ElbowPlot(sce, ndims = 50)
sce = RunUMAP(sce, dims = 1:8)
scDataPBMC_T = sce
# 保存
saveRDS(scDataPBMC_T, file = "cellAnnotation/scDataPBMC_T.rds")

scDataPBMC_T = readRDS(file = "cellAnnotation/scDataPBMC_T.rds")

# 可视化
cols = c("#696969", "#DAA520", "#FF4500", "#698B22",
         "#FF8C00", "#FFD700", "#FF69B4") # 指定每个组别的颜色
names(cols) = c("memory CD4+ T cells","memory CD8+ T cells","naive CD4+ T cells",
                "naive CD8+ T cells", "NKT", "MAIT","Tregs")
DimPlot(scDataPBMC_T, 
        group.by = "label1",
        label = TRUE,
        cols = cols, # 指定每个组别的颜色
        cells = rownames(scDataPBMC_T@meta.data)[which(scDataPBMC_T$label1 %in% c("memory CD4+ T cells","memory CD8+ T cells","naive CD4+ T cells",
                                                                                  "naive CD8+ T cells", "NKT", "MAIT","Tregs"))])

# LUAD----------------------------------------------------------
# 可视化细胞注释结果
cols = c("#00CED1", 
         #"#B22222", "#87CEEB",
         "#FFC0CB", "#8B0000",
         "#9ACD32", "#32CD32", "#9932CC",
         "#00BFFF", "#6A5ACD", "#1E90FF",
         "#FF1493") # 指定每个组别的颜色
names(cols) = c("NK",
                "Plasma cells","pre-B cells",
                "cDC","pDC","Monocytes",
                "Basophils","Neutrophils",
                "Macrophages","T cells")
DimPlot(scDataLUAD,
        label.size = 3,
        group.by = "major",
        label = TRUE,
        cols = cols, # 指定每个组别的颜色
        cells = rownames(scDataLUAD@meta.data)[which(scDataLUAD$major %in% c("NK",
                                                                             "Plasma cells","pre-B cells",
                                                                             "cDC","pDC","Monocytes",
                                                                             "Basophils","Neutrophils",
                                                                             "Macrophages","T cells"))])
# 保存
saveRDS(scDataLUAD, file = "cellAnnotation/scDataLUAD.rds")

scDataLUAD = readRDS(file = "cellAnnotation/scDataLUAD.rds")

# 具体的细胞注释有点复杂，不再赘述，可参考给出的marker gene
# 可视化T细胞亚群
cols = c("#696969", "#DAA520", "#FF4500", "#698B22",
         "#800080", "#00FF7F",
         "#FF8C00", "#FFD700", "#FF69B4") # 指定每个组别的颜色
names(cols) = c("memory CD4+ T cells","memory CD8+ T cells","naive CD4+ T cells",
                "naive CD8+ T cells", "T helper cells", "T effector cells",
                "NKT", "MAIT","Tregs")
DimPlot(scDataLUAD_T, group.by = "label0525", label = TRUE, cols = cols,
        cells = rownames(scDataLUAD_T@meta.data)[which(scDataLUAD_T$label0525 %in% c("memory CD4+ T cells", "memory CD8+ T cells",
                                                                                     "naive CD4+ T cells","naive CD8+ T cells",
                                                                                     "T helper cells","T effector cells",
                                                                                     "NKT", "MAIT", "Tregs"))])
# 保存
saveRDS(scDataLUAD_T, file = "cellAnnotation/scDataLUAD_T.rds")

scDataLUAD_T = readRDS(file = "cellAnnotation/scDataLUAD_T.rds")

count = scDataLUAD@assays[["RNA"]]@counts
saveRDS(count, file = "cellAnnotation/rawLUAD.rds")
count = scDataPBMC@assays[["RNA"]]@counts
saveRDS(count, file = "cellAnnotation/rawPBMC.rds")
tmp = scDataLUAD@meta.data
write.csv(tmp, file = "cellAnnotation/annLUAD.csv")
tmp = scDataPBMC@meta.data
write.csv(tmp, file = "cellAnnotation/annPBMC.csv")

