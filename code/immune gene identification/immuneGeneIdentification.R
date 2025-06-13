
# 从AmiGO网站上导出的所有免疫相关的通路
immunePathways = read.csv("immuneGeneIdentification/immuneSystemProcess.csv")
immuneGO = unique(immunePathways$ID)
write.csv(immuneGO, file = "immuneGeneIdentification/immuneGO.csv", row.names = F)

# 从GO_DATA中抽取免疫相关通路的基因信息
immuneGO = names(GO_DATA$PATHID2EXTID)[which(names(GO_DATA$PATHID2EXTID) %in% immuneGO)] # 筛选出GO_DATA中的immuneGO
immuneGOInfo = data.frame(ID = NULL, Description = NULL, geneSet = NULL, geneCount = NULL)
goIDs = c()
description = c()
geneSet = c()
geneCount = c()

for(goID in immuneGO){
  goIDs = c(goIDs, goID)
  description = c(description, as.data.frame(GO_DATA$PATHID2NAME[goID])[,1])
  geneAll =  unique(unlist(GO_DATA$PATHID2EXTID[goID]))
  geneSet = c(geneSet, paste(geneAll, collapse = "/"))
  geneCount = c(geneCount, length(geneAll))
}
immuneGOInfo = data.frame(ID = goIDs, Description = description, geneSet = geneSet, geneCount = geneCount)
write.csv(immuneGOInfo, file = "immuneGeneIdentification/immuneGOInfo.csv", row.names = F)

# 免疫基因集immuneGeneSet
immuneGeneSet = unique(unlist(lapply(strsplit(as.character(immuneGOInfo$geneSet), "/"), trimws)))
immuneGeneSet = data.frame(immuneGene = immuneGeneSet)
write.csv(immuneGeneSet, file = "immuneGeneIdentification/immuneGeneSet.csv", row.names = F)

