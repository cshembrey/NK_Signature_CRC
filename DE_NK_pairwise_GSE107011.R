
# Xu et al
# download Xu et al data and assign the dataPath to where you include the data
# dataPath <- "/Users/foroutanm/Documents/data/blood/GEO/GSE107011_RNAseq_bloodCells/"
# tpm <- read.table(paste0(dataPath, "GSE107011_Processed_data_TPM.txt"), header = T, sep = "\t", stringsAsFactors = F, check.names = F)
# head(tpm)
# 
# ##------------------- Map IDs and generate logTPM
# tpm$ENSG <- sapply(tpm[, 1], function(x) unlist(strsplit(x, "[.]"))[[1]])
# ### remove 45 duplicated names
# tpm <- tpm[!duplicated(tpm$ENSG), ]
# 
# library(org.Hs.eg.db)               
# tpm$GeneIDs = mapIds(
#   org.Hs.eg.db,                     
#   keys = tpm$ENSG,          
#   keytype = "ENSEMBL",              
#   column = "ENTREZID",              
#   multiVals = "asNA"                
# )                                 
# 
# ## remove rows that do not map to a gene ID
# tpm <- tpm[complete.cases(tpm$GeneIDs), ]
# tpm <- tpm[ ! duplicated(tpm$GeneIDs), ]
# 
# dd <- as.matrix(tpm[, 2:(ncol(tpm)-2)])
# row.names(dd) <- tpm$GeneIDs
# 
# ## Remove PBMC samples
# dd <- dd[ , ! grepl("PBMC", colnames(dd))]
# 
# logTPM <- log2(dd+1)

# write.table(logTPM, "logTPM_GSE107011_bloodCells_2019.txt", row.names = T, sep = "\t")
logTPM <- read.table( "./data/logTPM_GSE107011_bloodCells_2019.txt", header = T, sep = "\t", check.names = F)

##---------------- filter 
kp <- rowSums(logTPM > 2) >= 4
summary(kp)

logTPMFilt <- logTPM[kp, ]

##---------------- annotation data:
ann <- read.csv("./data/annot_GSE107011.csv", stringsAsFactors = F)
ann <- ann[ ! grepl("PBMC", ann$X.Sample_description.2), ]

ann$cellType <- sapply(ann$X.Sample_description.2, function(x)
  unlist(strsplit(x, "[_]"))[[2]])

ann$cellType[ann$cellType == "mDC" | ann$cellType == "pDC"] <- "DC"
ann$cellType[ann$cellType == "C" | ann$cellType == "I"] <- "Mono"
ann$cellType[ann$cellType == "NC"] <- "Mono"

ann$cellType[grepl("^Th", ann$cellType)] <- "Th"
## Follicular helper T cells
ann$cellType[ann$cellType == "TFH"] <- "Th"
ann$cellType[ann$cellType == "VD2+"] <- "VD2p_gd"
ann$cellType[ann$cellType == "VD2-"] <- "VD2n_gd"


ann$Group[grepl("NK", ann$cellType)] <- "NK"
ann$Group[! grepl("NK", ann$cellType)] <- "nonNK"

all(ann$X.Sample_description.2 == colnames(logTPM))

##=============================== DE NK vs Non-NK
##----------------- Do DE analysis through limma-trend
## eBayes with trend=TRUE and by using arrayWeights() to try to partially recover the library sizes
library(limma)

design <- model.matrix(~ 0 + ann$Group)
colnames(design) <- c("NK", "nonNK")
design

contr.matrix <- makeContrasts(
  NK_nonNK = NK - nonNK,
  levels = colnames(design))

fit <- lmFit(logTPMFilt, design)
fit <- contrasts.fit(fit, contr.matrix)
fit <- eBayes(fit, trend = TRUE)
dt <- decideTests(fit)

DEstats   <- topTable(fit, coef = 1, n = Inf)

##----------------- Add Symbols to the results of DE analysis and export DE stat and DEGs
DEstats$Symbols <-
  mapIds(
    org.Hs.eg.db,
    keys = row.names(DEstats),
    keytype = "ENTREZID",
    column = c("SYMBOL"),
  ) 
DEstats$GeneIDs <- row.names(DEstats)

write.table(DEstats, 
  "./output/DEStat_NK_nonNK_GSE107011.txt", row.names = F, sep = "\t")

DEGs <- DEstats[DEstats$GeneIDs %in% row.names(dt)[! dt[, 1] == 0], ]

write.table(DEGs, 
  "./output/DEGs_NK_nonNK_GSE107011.txt", row.names = F, sep = "\t")

##--- only genes with positve log FC

DEGs_pos <- DEGs[DEGs$logFC > 0, ]

write.table(DEGs_pos, 
  "./output/DEGs_positiveLogFC_NK_nonNK_GSE107011.txt", row.names = F, sep = "\t")

##---------- Generate heatmap
d <- logTPMFilt[DEGs_pos$GeneIDs, ]
row.names(d) <- DEGs_pos$Symbols
d <- t(scale(t(d)))
minExpr <- min(d)
maxExpr <- max(d)


textSize <- 8

library(ComplexHeatmap)

p <- ComplexHeatmap::Heatmap(
  d,
  col = circlize::colorRamp2(c(minExpr, maxExpr), c("yellow3" , "navy")),
  name = "logExpr",
  column_title = "",
  cluster_rows = T,
  cluster_columns = T,
  show_row_names = F,
  show_column_names = T,
  row_names_gp = gpar(fontsize = textSize),
  column_names_gp = gpar(fontsize = textSize),
  # column_title_gp = gpar(fontsize = textSize), 
  heatmap_legend_param = list(
    title_position = "topcenter",
    title = "Scaled log(TPM)",
    color_bar = "continuous", 
    legend_direction = "horizontal",
    legend_width = unit(5, "cm")
    # legend_height = unit(3, "cm")
  ),
  na_col = "gray90"
)

pdf("./figure/DEGs_positiveLogFC_GSE107011_noNames.pdf", height = 10, width = 12)
ComplexHeatmap::draw(
  p,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "right",
  legend_title_gp = gpar(fontsize = textSize, fontface = "bold")
)
dev.off()


##======================== DE NK vs cell types


design <- model.matrix(~ 0 + ann$cellType)
colnames(design) <-
  c(
    "B",
    "Bas",
    "CD4",
    "CD8",
    "DC",
    "MAIT",
    "Mono",
    "Neut",
    "NK",
    "Plasm",
    "Prog",
    "Th",
    "Treg",
    "VD2n_gd",
    "VD2p_gd"
  )
design

contr.matrix <- makeContrasts(
  NK_B = NK - B,
  NK_Bas = NK - Bas,
  NK_CD4 = NK - CD4,
  NK_CD8 = NK - CD8,
  NK_DC = NK - DC,
  NK_MAIT = NK - MAIT,
  NK_Mono = NK - Mono,
  NK_Neut = NK - Neut,
  NK_Th = NK - Th,
  NK_Treg = NK - Treg,
  NK_Tgdn = NK - VD2n_gd,
  NK_Tgdp = NK - VD2p_gd,
  NK_Plasm = NK - Plasm,
  NK_Prog = NK - Prog,
  levels = colnames(design))

fit <- lmFit(logTPMFilt, design)
fit <- contrasts.fit(fit, contr.matrix)

efit <- eBayes(fit, trend = TRUE)
dt0 <- decideTests(efit)
summary(dt0)

tfit <- treat(fit, lfc=1)
dt <- decideTests(tfit)
summary(dt)


##----------------- Add Symbols to the results of DE analysis and export DE stat and DEGs

geneid <- rownames(logTPMFilt) 

genes <-
  select(
    Homo.sapiens,
    keys = geneid,
    columns = c("SYMBOL", "TXCHROM"),
    keytype = "ENTREZID"
  ) 
head(genes)

genes <- genes[!duplicated(genes$ENTREZID), ]


getDEGs <- function(fitObject, comparison, dt, geneIDs, DEtype = "treat"){
  DEStat   <- topTreat(fitObject, coef = comparison, n = Inf)
  DEStat$GeneID <- rownames(DEStat)
  DEStat <- merge(geneIDs, DEStat, by.x = "ENTREZID", by.y = "GeneID")
  DEGs <- DEStat[DEStat$ENTREZID %in% row.names(dt)[ ! dt[, comparison] == 0], ]
  ## remove genes without Symbols:
  DEGs <- DEGs[complete.cases(DEGs$SYMBOL) & complete.cases(DEGs$TXCHROM),]
  DEGs$Comparison <- comparison
  return(DEGs)
}

DEGsList <- list()
for (i in colnames(contr.matrix)){
  DEGsList[[i]] <- getDEGs(fitObject = tfit, comparison = i, dt = dt, geneIDs = genes, DEtype = "treat")
}

DEGsList[["DEGs_NK_nonNK"]] <- DEGs

GeneList_GSE107011 <- DEGsList

comDEGs <- intersect(GeneList_GSE107011$DEGs_NK_nonNK$Symbols[GeneList_GSE107011$DEGs_NK_nonNK$logFC > 0], 
  GeneList_GSE60424$DEGs_NK_nonNK$SYMBOL[GeneList_GSE60424$DEGs_NK_nonNK$DEG_type =="Treat_Str" & GeneList_GSE60424$DEGs_NK_nonNK$logFC > 0])


symbolList <- lapply(GeneList_GSE107011[1:(length(GeneList_GSE107011)-1)], function(x) x$SYMBOL)
Reduce(intersect, symbolList )

DEGsComAcrossAll <- GeneList_GSE107011$NK_Tgdp[GeneList_GSE107011$NK_Tgdp$SYMBOL  %in%   Reduce(intersect, symbolList ),]

comCD4 <- intersect(GeneList_GSE60424$DEGs_NK_CD4$SYMBOL[GeneList_GSE60424$DEGs_NK_CD4$logFC > 0], GeneList_GSE107011$NK_CD4$SYMBOL[GeneList_GSE107011$NK_CD4$logFC > 0])

comCD8 <- intersect(GeneList_GSE60424$DEGs_NK_CD8$SYMBOL[GeneList_GSE60424$DEGs_NK_CD8$logFC > 0], GeneList_GSE107011$NK_CD8$SYMBOL[GeneList_GSE107011$NK_CD8$logFC > 0])

## CAMK4 is common across all T cells
NKdownSet2 <- rbind(
  GeneList_GSE107011$NK_CD4[GeneList_GSE107011$NK_CD4$logFC < -5,],
  GeneList_GSE107011$NK_CD8[GeneList_GSE107011$NK_CD8$logFC < -5,],
  GeneList_GSE107011$NK_Th[GeneList_GSE107011$NK_Th$logFC < -5,],  ## -6
  GeneList_GSE107011$NK_Tgdn[GeneList_GSE107011$NK_Tgdn$logFC < -3,],  ## -6  -4: TRAT1 and CAMK4
  GeneList_GSE107011$NK_Tgdp[GeneList_GSE107011$NK_Tgdp$logFC < -3,],  ## -6  -4: TRAT1 and CAMK4
  GeneList_GSE107011$NK_Th[GeneList_GSE107011$NK_Th$logFC < -5,],  ## -6
  GeneList_GSE107011$NK_Treg[GeneList_GSE107011$NK_Treg$logFC < -5,],  ## -6
  GeneList_GSE107011$NK_MAIT[GeneList_GSE107011$NK_MAIT$logFC < -6,],  ## -5
 
  GeneList_GSE107011$NK_B[GeneList_GSE107011$NK_B$logFC < -7,],
  GeneList_GSE107011$NK_Bas[GeneList_GSE107011$NK_Bas$logFC < -9,],
  GeneList_GSE107011$NK_DC[GeneList_GSE107011$NK_DC$logFC < -6,],
  GeneList_GSE107011$NK_Mono[GeneList_GSE107011$NK_Mono$logFC < -7,],  ## -6
  GeneList_GSE107011$NK_Neut[GeneList_GSE107011$NK_Neut$logFC < -9,],  ## -8
  GeneList_GSE107011$NK_Plasm[GeneList_GSE107011$NK_Plasm$logFC < -9,],  ## IGLL5 and JCHAIN
  GeneList_GSE107011$NK_Prog[GeneList_GSE107011$NK_Prog$logFC < -9,]  
)


NKdownSet2 <- NKdownSet2 %>% group_by(SYMBOL) %>%  add_tally() %>%  data.frame()
NKdownSetGenes2 <- NKdownSet2$SYMBOL
names(NKdownSetGenes2) <- NKdownSet2$Comparison
NKdownSetGenes2[duplicated(NKdownSetGenes2)]
## 80 genes are duplicated.
NKdownSetGenes2 <- NKdownSetGenes2[!duplicated(NKdownSetGenes2)]



NKdownSet4 <-
  rbind(
    GeneList_GSE107011$NK_CD4[GeneList_GSE107011$NK_CD4$logFC < -5,],
    GeneList_GSE107011$NK_CD8[GeneList_GSE107011$NK_CD8$logFC < -5,], 
    GeneList_GSE107011$NK_B[GeneList_GSE107011$NK_B$logFC < -7,],
    GeneList_GSE107011$NK_Mono[GeneList_GSE107011$NK_Mono$logFC < -7,],  ## -6
    GeneList_GSE107011$NK_Neut[GeneList_GSE107011$NK_Neut$logFC < -9,])


NKdownSet4 <- NKdownSet4 %>% group_by(SYMBOL) %>%  add_tally() %>%  data.frame()
NKdownSetGenes4 <- NKdownSet4$SYMBOL
names(NKdownSetGenes4) <- NKdownSet4$Comparison
NKdownSetGenes4[!duplicated(NKdownSetGenes4)]

## NKdownSetGenes is from teh otehr study
comCells_NKDown <- NKdownSetGenes4[NKdownSetGenes4 %in% NKdownSetGenes]



## One data set has only CD4 and CD8 for T cells, so we add a few genes from other data set

## "SLAMF1" "CAMK4"  are up regulated in all the other T cells listed below
NKdownSet3 <-
  rbind(
    GeneList_GSE107011$NK_Th[GeneList_GSE107011$NK_Th$logFC < -6, ], ## -6
    GeneList_GSE107011$NK_Tgdn[GeneList_GSE107011$NK_Tgdn$logFC < -3, ],## -6  -4: TRAT1 and CAMK4
    GeneList_GSE107011$NK_Tgdp[GeneList_GSE107011$NK_Tgdp$logFC < -3, ], ## -6  -4: TRAT1 and CAMK4
    GeneList_GSE107011$NK_Th[GeneList_GSE107011$NK_Th$logFC < -6, ], ## -6
    GeneList_GSE107011$NK_Treg[GeneList_GSE107011$NK_Treg$logFC < -6, ],## -6
    GeneList_GSE107011$NK_MAIT[GeneList_GSE107011$NK_MAIT$logFC < -6, ],
    GeneList_GSE107011$NK_Bas[GeneList_GSE107011$NK_Bas$logFC < -10, ],
    GeneList_GSE107011$NK_DC[GeneList_GSE107011$NK_DC$logFC < -7, ]
  )

NKdownSet3 <- NKdownSet3 %>% group_by(SYMBOL) %>%  add_tally() %>%  data.frame()
NKdownSetGenes3 <- NKdownSet3$SYMBOL
names(NKdownSetGenes3) <- NKdownSet3$Comparison
NKdownSetGenes3 <- NKdownSetGenes3[!duplicated(NKdownSetGenes3)]

otherCells_NKDown <- NKdownSetGenes3[! NKdownSetGenes3  %in% comCells_NKDown]

## Also look at genes down regulated in NK vs all in two data sets: 
## "SPINT2" "SUSD3"  "GABBR1" are common
NKdownSet5 <- GeneList_GSE107011$DEGs_NK_nonNK$Symbols[GeneList_GSE107011$DEGs_NK_nonNK$logFC < -2]
# names(NKdownSet5) <- "DEGsDown_GSE107011"
NKdownSet6 <- GeneList_GSE60424$DEGs_NK_nonNK$SYMBOL[GeneList_GSE60424$DEGs_NK_nonNK$logFC < -3]
# names(NKdownSet6) <- "DEGsDown_GSE60424"

NKdownSetOverallDEG <- c(NKdownSet5, NKdownSet6)
names(NKdownSetOverallDEG) <- c(rep("DEGsDown_GSE107011", length(NKdownSet5)), rep("DEGsDown_GSE60424", length(NKdownSet6)))
names(NKdownSetOverallDEG)[NKdownSetOverallDEG  %in% c( "SPINT2", "SUSD3", "GABBR1")] <- "DEGsDown_Common"

NKdownSetOverallDEG <- NKdownSetOverallDEG[!duplicated(NKdownSetOverallDEG)]
NKdownSetOverallDEG[ NKdownSetOverallDEG  %in% NKDown]
# DEGsDown_Common DEGsDown_GSE107011 DEGsDown_GSE107011 
# "SUSD3"            "SOCS3"             "SIT1" 

## remove the overlapping ones:
NKdownSetOverallDEG <- NKdownSetOverallDEG[ ! NKdownSetOverallDEG  %in% NKDown]

NKDown <- c(comCells_NKDown, otherCells_NKDown, NKdownSetOverallDEG)

NKDownData <- data.frame(NKDown)
NKDownData$Comparison <- names(NKDown)
NKDownData$Data <- c(rep("PairwiseDEGs_BothBlood", length(comCells_NKDown)), rep("PairwiseDEGs_GSE107011", length(otherCells_NKDown)), NKDownData$Comparison[grepl("DEGsDown", NKDownData$Comparison)])

NKDownData <- NKDownData[!duplicated(NKDownData$NKDown), ]
rownames(NKDownData) <- NKDownData$NKDown
NKDownData$NKDown <- as.character(NKDownData$NKDown)


GeneList_GSE107011[["NK_downSet"]] <- NKDownData
GeneList_GSE107011[["CommonDEGs_NK_nonNK"]] <- comDEGs

saveRDS(GeneList_GSE107011, "./output/DEGs_List_NK_GSE107011.RDS")



gg <- NKDownData$NKDown[NKDownData$NKDown  %in% rownames(dd)]
aa <- NKDownData[gg, ]

dd <- merge(genes, logTPMFilt, by.x = "ENTREZID", by.y = "row.names")
dd <- dd[complete.cases(dd$SYMBOL), ]
dd <- dd[!duplicated(dd$SYMBOL), ]
rownames(dd) <- dd$SYMBOL
dd <- dd[, -c(1, 2, 3)]

d <- dd[gg, ]

dataCols <- brewer.pal(5, "Set1")
names(dataCols) <- unique(aa$Data)

gAnnot <-
  ComplexHeatmap::rowAnnotation(df = aa[, "Data", drop = F], col =  list(Data = dataCols))

sampleCol <- c("Green3", "Orange", brewer.pal(9, "Purples")[c(3, 4)], "Orangered", brewer.pal(9, "Purples")[c(5)], "blue", "lightblue", "black", "Gray50", "Gray70", brewer.pal(9, "Purples")[6:9])
names(sampleCol) <- names(table(ann$cellType))
sAnnot <- 
  ComplexHeatmap::columnAnnotation(df = ann[, "cellType", drop = F], col =  list(cellType = sampleCol))

pdf("./figure/Heatmap_DEGs_Down_NK.pdf", height = 18, width = 18)
p <- ComplexHeatmap::Heatmap(d, left_annotation = gAnnot, top_annotation = sAnnot)
draw(p, annotation_legend_side = "bottom", heatmap_legend_side = "bottom")
dev.off()



pdf("./figure/Heatmap_DEGsCom_Up_NK.pdf", height = 30, width = 20)
p <- ComplexHeatmap::Heatmap(dd[comDEGs, ], top_annotation = sAnnot)
draw(p, annotation_legend_side = "right", heatmap_legend_side = "right")
dev.off()




## ---- get genes that have positive logFC in Nk vs every other cell in both data:
##======== Data 1
nk_up_g1 <- lapply(GeneList_GSE60424[2:6], function(x) x$SYMBOL[x$logFC > 0])
## 56 UP genes overlapping between all the NK vs X Cell comparisons
nk_up_d1 <- Reduce(intersect, nk_up_g1)

nk_dn_g1 <- lapply(GeneList_GSE60424[2:6], function(x) x$SYMBOL[x$logFC < 0])
## 3 Down genes overlapping between all the NK vs X Cell comparisons
nk_dn_d1 <- Reduce(intersect, nk_dn_g1) ## "ID3"   "KCNQ1" "SNN"


##======= Data 2
nk_up_g2 <- lapply(GeneList_GSE107011[1:12], function(x) x$SYMBOL[x$logFC > 0])
## 9 UP genes are overlapping between NK and X cell type in data 2
nk_up_d2 <- Reduce(intersect, nk_up_g2)

## 5 genes are overlapping across both data:
nk_up_d2[nk_up_d2  %in%  nk_up_d1]
# [1] "SH2D1B" "LINGO2" "RNF165" "ZMAT4"  "LDB2" 

nk_dn_g2 <- lapply(GeneList_GSE107011[1:12], function(x) x$SYMBOL[x$logFC < 0])
## 3 Down genes overlapping between all the NK vs X Cell comparisons
nk_dn_d2 <- Reduce(intersect, nk_dn_g2) ## None!



##------ merge the two up_NK_data together and see how many of them overlap with 
## those from NK vs non_NK:

## 60 unique genes
nk_up <- unique(c(nk_up_d1, nk_up_d2))

## 10 genes do not present in the NK vs non_NK obtained in two blood data, so we add them to the list
add_nk_up <- nk_up[! nk_up  %in% GeneList_GSE107011$CommonDEGs_NK_nonNK]
# [1] "USP43"      "C19orf84"   "GRIK4"      "KIFC3"      "KRT81"      "MEIS1"     
# [7] "FAT4"       "BAALC"      "CABLES1"    "MIR181A2HG"

## 195 genes
upGenes <- c(GeneList_GSE107011$CommonDEGs_NK_nonNK, add_nk_up)




Tcolumns <- c(3, 4, 6, 9, 10, 11, 12)
nkT_up_g2 <- lapply(GeneList_GSE107011[Tcolumns], function(x) x$SYMBOL[x$logFC > 0])
## 9 UP genes are overlapping between NK and X cell type in data 2
nkT_up_d2 <- Reduce(intersect, nkT_up_g2)
# [1] "MIR181A2HG" "ADAM28"     "SH2D1B"     "COL13A1"    "LINGO2"     "MPEG1"     
# [7] "FCER1G"     "FES"        "LGALS9B"    "GSN"        "HBA1"       "HBA2"      
# [13] "HBB"        "HCK"        "IGFBP7"     "MEF2C"      "MEIS1"      "RNF165"    
# [19] "LGALS9C"    "SYK"        "CCNJL"      "ZMAT4"      "LDB2"      

nkT_up_d2[!nkT_up_d2 %in% upGenes]
# [1] "ADAM28" "MPEG1"  "FCER1G" "FES"    "GSN"    "HBA1"   "HBA2"   "HBB"    "HCK"   
# [10] "MEF2C"  "SYK" 


NK_genelist <- list(
  NK_upSet = upGenes,
  NK_downSet = GeneList_GSE107011$NK_downSet, 
  NK_panT_up = nkT_up_d2)

saveRDS(NK_genelist, "./output/NK_genelist_2bld.RDS")

