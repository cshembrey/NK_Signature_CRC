library(edgeR)
library(Homo.sapiens) 
dataPath <- "./data/"

##---- Read combined NK signature genes
# nkSig <- read.csv(paste0(dataPath, "Cursons_Guimaraes_NKsignature_CIR_2019.csv"),
#   stringsAsFactors = F)

# nk_signature <- nkSig$HGNC.Symbol

## Make sure that we have updated gene symbols
# nkIDs <- alias2SymbolUsingNCBI(nk_signature, 
#   paste0(dataPath, "Homo_sapiens.gene_info"))
# 
# nkIDs <- nkIDs[complete.cases(nkIDs$GeneID), ]

##----- Read DGEList object for the expression
load(paste0(dataPath, "bldCells_DGEList.RData"))
dge <- bldCells_DGEList

##---- Divide samples into Nk and non_NK samples
dge$samples$celltype_s <- as.character(dge$samples$celltype_s)
dge$samples$groups[dge$samples$celltype_s == "NK"] <- "NK"
dge$samples$groups[ ! dge$samples$celltype_s == "NK"] <- "non-NK"

##---- remove the whole blood samples
dge <- dge[, ! dge$samples$celltype_s == "Whole Blood", ]

##---- Filter genes
kp <- rowSums(cpm(dge$counts) > 2) >= 4 
summary(kp)
dgeFilt <- dge[kp, ]

##---- Perform TMM normalisation
dgeFiltNorm <- calcNormFactors(dgeFilt)

##================================== NK vs Others
##----- Perform DE analysis (we use voom-limma + TREAT)
design <- model.matrix(~ 0 + dge$samples$groups)
colnames(design) <- c("NK", "nonNK")

contr.matrix <- makeContrasts(
  NK_nonNK = NK - nonNK,
  levels = colnames(design))

v <- voom(dgeFiltNorm, design, plot=TRUE)

## Two genes from the combined NK signature gene list do not present in the expression data:
# nkIDs$GeneID[ ! nkIDs$GeneID %in% row.names(dgeFiltNorm$counts)]
## "3803"  "28639"

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)

## First we use eBayes but it gives many DEGs
efit <- eBayes(vfit)
plotSA(efit)
summary(decideTests(efit))

dt0 <- decideTests(efit)
summary(decideTests(efit))


##----- Map gene IDs to symbols
geneid <- rownames(dgeFiltNorm) 

genes <-
  select(
    Homo.sapiens,
    keys = geneid,
    columns = c("SYMBOL", "TXCHROM"),
    keytype = "ENTREZID"
  ) 
head(genes)

genes <- genes[!duplicated(genes$ENTREZID), ]

##------

DEstats_eBayes   <- topTable(efit, coef = 1, n = Inf)
DEstatsID_eBayes <- merge(genes, DEstats_eBayes, by.x = "ENTREZID", by.y = "GeneID")


DEGs_pos_eBayes <- DEstatsID_eBayes[DEstatsID_eBayes$ENTREZID %in% row.names(dt0)[ dt0[, 1] > 0], 
  c(1, 2, 3, 8:13) ]

DEGs_eBayes <- DEstatsID_eBayes[DEstatsID_eBayes$ENTREZID %in% row.names(dt0)[ ! dt0[, 1] == 0], 
  c(1, 2, 3, 8:13) ]

## remove genes without Symbols:
DEGs_pos_eBayes <- DEGs_pos_eBayes[complete.cases(DEGs_pos_eBayes$SYMBOL),]

##-------------------- Use TREAT to lower the number of DEGs
tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)


DEstats_treat   <- topTreat(tfit, coef = 1, n = Inf)
DEstatsID_treat <- merge(genes, DEstats_treat, by.x = "ENTREZID", by.y = "GeneID")

DEGs_pos_treat <- DEstatsID_treat[DEstatsID_treat$ENTREZID %in% row.names(dt)[dt[, 1] > 0], 
  c(1, 2, 3, 8:13) ]
DEGs_treat <- DEstatsID_treat[DEstatsID_treat$ENTREZID %in% row.names(dt)[!dt[, 1] == 0], 
  c(1, 2, 3, 8:13) ]

## remove genes without Symbols:
DEGs_pos_treat <- DEGs_pos_treat[complete.cases(DEGs_pos_treat$SYMBOL),]


DEGs_pos_eBayes$DEG_type <- "eBayes"
DEGs_pos_eBayes$DEG_type[DEGs_pos_eBayes$SYMBOL  %in%  DEGs_pos_treat$SYMBOL] <- "Treat"

DEGs_eBayes$DEG_type <- "eBayes_lax"
DEGs_eBayes$DEG_type[DEGs_eBayes$SYMBOL  %in%  DEGs_treat$SYMBOL] <- "Treat_Str"

## remove those with NA as Chr:
DEGs_pos_eBayes <- DEGs_pos_eBayes[complete.cases(DEGs_pos_eBayes$TXCHROM), ]

write.table(DEGs_eBayes,
  "./output/DEGs_bldCells_NK_nonNK_eBayes_Treat_UpDown.txt", row.names = F, sep = "\t")
write.table(DEGs_pos_eBayes,
  "./output/DEGs_bldCells_NK_nonNK_eBayes_Treat.txt", row.names = F, sep = "\t")
write.table(DEstatsID_eBayes,
  "./output/DEStat_bldCells_NK_nonNK_eBayes.txt", row.names = F, sep = "\t")




##================================== NK vs every cell types


design <- model.matrix(~ 0 + dge$samples$celltype_s)
colnames(design) <- c("B", "CD4", "CD8", "Mono", "Neut", "NK")

contr.matrix <- makeContrasts(
  NK_CD4 = NK - CD4,
  NK_CD8 = NK - CD8,
  NK_B = NK - B,
  NK_Mono = NK - Mono,
  NK_Neut = NK - Neut,
  levels = colnames(design))

v <- voom(dgeFiltNorm, design, plot=TRUE)

## Two genes from the combined NK signature gene list do not present in the expression data:
# nkIDs$GeneID[ ! nkIDs$GeneID %in% row.names(dgeFiltNorm$counts)]
## "3803"  "28639"

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)

## First we use eBayes but it gives many DEGs
efit <- eBayes(vfit)
plotSA(efit)
summary(decideTests(efit))

##---- Use TREAT to lower the number of DEGs
tfit <- treat(vfit, lfc=2)
dt <- decideTests(tfit)
summary(dt)


DEStat_NK_CD4   <- topTreat(tfit, coef = "NK_CD4", n = Inf)
DEstatsID_NK_CD4 <- merge(genes, DEStat_NK_CD4, by.x = "ENTREZID", by.y = "GeneID")
DEGs_NK_CD4 <- DEstatsID_NK_CD4[DEstatsID_NK_CD4$ENTREZID %in% row.names(dt)[ ! dt[, "NK_CD4"] == 0], 
  c(1, 2, 3, 8:13) ]
## remove genes without Symbols:
DEGs_NK_CD4 <- DEGs_NK_CD4[complete.cases(DEGs_NK_CD4$SYMBOL) & complete.cases(DEGs_NK_CD4$TXCHROM),]
DEGs_NK_CD4$Comparison <- "NK_CD4"


DEStat_NK_CD8   <- topTreat(tfit, coef = "NK_CD8", n = Inf)
DEstatsID_NK_CD8 <- merge(genes, DEStat_NK_CD8, by.x = "ENTREZID", by.y = "GeneID")
DEGs_NK_CD8 <- DEstatsID_NK_CD8[DEstatsID_NK_CD8$ENTREZID %in% row.names(dt)[ ! dt[, "NK_CD8"] == 0], 
  c(1, 2, 3, 8:13) ]
## remove genes without Symbols:
DEGs_NK_CD8 <- DEGs_NK_CD8[complete.cases(DEGs_NK_CD8$SYMBOL) & complete.cases(DEGs_NK_CD8$TXCHROM),]
DEGs_NK_CD8$Comparison <- "NK_CD8"


DEStat_NK_B   <- topTreat(tfit, coef = "NK_B", n = Inf)
DEstatsID_NK_B <- merge(genes, DEStat_NK_B, by.x = "ENTREZID", by.y = "GeneID")
DEGs_NK_B <- DEstatsID_NK_B[DEstatsID_NK_B$ENTREZID %in% row.names(dt)[ ! dt[, "NK_B"] == 0], 
  c(1, 2, 3, 8:13) ]
## remove genes without Symbols:
DEGs_NK_B <- DEGs_NK_B[complete.cases(DEGs_NK_B$SYMBOL) & complete.cases(DEGs_NK_B$TXCHROM),]
DEGs_NK_B$Comparison <- "NK_B"


DEStat_NK_Mono   <- topTreat(tfit, coef = "NK_Mono", n = Inf)
DEstatsID_NK_Mono <- merge(genes, DEStat_NK_Mono, by.x = "ENTREZID", by.y = "GeneID")
DEGs_NK_Mono <- DEstatsID_NK_Mono[DEstatsID_NK_Mono$ENTREZID %in% row.names(dt)[ ! dt[, "NK_Mono"] == 0], 
  c(1, 2, 3, 8:13) ]
## remove genes without Symbols:
DEGs_NK_Mono <- DEGs_NK_Mono[complete.cases(DEGs_NK_Mono$SYMBOL) & complete.cases(DEGs_NK_Mono$TXCHROM),]
DEGs_NK_Mono$Comparison <- "NK_Mono"


DEStat_NK_Neut   <- topTreat(tfit, coef = "NK_Neut", n = Inf)
DEstatsID_NK_Neut <- merge(genes, DEStat_NK_Neut, by.x = "ENTREZID", by.y = "GeneID")
DEGs_NK_Neut <- DEstatsID_NK_Neut[DEstatsID_NK_Neut$ENTREZID %in% row.names(dt)[ ! dt[, "NK_Neut"] == 0], 
  c(1, 2, 3, 8:13) ]
## remove genes without Symbols:
DEGs_NK_Neut <- DEGs_NK_Neut[complete.cases(DEGs_NK_Neut$SYMBOL) & complete.cases(DEGs_NK_Neut$TXCHROM),]
DEGs_NK_Neut$Comparison <- "NK_Neut"




##==================== define up and down sets for NK genes
nkGenes <- rbind(
  DEGs_NK_CD4[DEGs_NK_CD4$logFC > 0, ],
  DEGs_NK_CD8[DEGs_NK_CD8$logFC > 0, ],
  DEGs_NK_B[DEGs_NK_B$logFC > 0, ],
  DEGs_NK_Mono[DEGs_NK_Mono$logFC > 0, ],
  DEGs_NK_Neut[DEGs_NK_Neut$logFC > 0, ]
)
nkGenes <- nkGenes %>% group_by(SYMBOL) %>%  add_tally() %>%  data.frame()
uu <- unique(nkGenes$SYMBOL[nkGenes$n == 5])[order(unique(nkGenes$SYMBOL[nkGenes$n == 5]))]
# [1] "ADAMTS1"   "ADGRB2"    "AKR1C3"    "ARVCF"     "B3GNT7"    "BAALC"    
# [7] "BNC2"      "C19orf84"  "CABLES1"   "CLIC3"     "COLGALT2"  "DTHD1"    
# [13] "ERBB2"     "FAT4"      "FKBP10"    "GRIK4"     "GZMB"      "HOXA10"   
# [19] "IL12RB2"   "IL2RB"     "KIFC3"     "KIR2DL1"   "KIR3DX1"   "KLRC1"    
# [25] "KLRF1"     "KRT81"     "KRT86"     "LDB2"      "LIM2"      "LINC00298"
# [31] "LINC00299" "LINGO2"    "LRRC43"    "MEIS1"     "MLC1"      "NCAM1"    
# [37] "NCR1"      "NLRP7"     "NMUR1"     "PCDH1"     "PDGFRB"    "PRSS57"   
# [43] "PTGDR"     "RAB27B"    "RAMP1"     "RHOBTB3"   "RNF165"    "SH2D1B"   
# [49] "SPON2"     "SPTSSB"    "TIE1"      "TKTL1"     "TNFRSF11A" "USP43"    
# [55] "XCL2"      "ZMAT4"   

nkGenes_CDs <- intersect(
  DEGs_NK_CD4$SYMBOL[DEGs_NK_CD4$logFC > 8], 
  DEGs_NK_CD8$SYMBOL[DEGs_NK_CD8$logFC > 6])

addCDs <- nkGenes_CDs[!nkGenes_CDs  %in% uu]

uu2 <- c(uu, addCDs)
## This seems to be very specific to T cells: SIGLEC7 --> Mediates inhibition of natural killer cells cytotoxicity
# [1] "COL13A1" "SIGLEC7" "GUCY1A1" "IGFBP7"  "PPP1R9A" "SYK"     "CCNJL"  
# [8] "ATP8B4"  




##------ genes that are down regulated in NK vs at least 4 other cell types:
nonNK_Genes <- rbind(
  DEGs_NK_CD4[DEGs_NK_CD4$logFC < 0, ],
  DEGs_NK_CD8[DEGs_NK_CD8$logFC < 0, ],
  DEGs_NK_B[DEGs_NK_B$logFC < 0, ],
  DEGs_NK_Mono[DEGs_NK_Mono$logFC < 0, ],
  DEGs_NK_Neut[DEGs_NK_Neut$logFC < 0, ]
)

nonNK_Genes <- nonNK_Genes %>% group_by(SYMBOL) %>%  add_tally() %>%  data.frame()

## Three geens are always there: "ID3"   "KCNQ1" "SNN"
nonNK_Genes5 <- nonNK_Genes[nonNK_Genes$n >= 4 , ]
nonNK_Genes5Unique <- nonNK_Genes5[!duplicated(nonNK_Genes5$SYMBOL), ]
# [1] "ECE1-AS1"   "SPINT2"     "BTBD19"     "BTLA"       "MEGF6"     
# [6] "SUSD3"      "GABBR1"     "ARMH1"      "ID3"        "APP"       
# [11] "IL6R"       "IRS1"       "CD82"       "KCNN4"      "KCNQ1"     
# [16] "KRT18"      "LRP6"       "LTB"        "LIMS2"      "PTPRO"     
# [21] "CDH23"      "SPINT1"     "SNN"        "ACTN1"      "ZSWIM1"    
# [26] "ST6GALNAC3" "YBX3P1"     "YBX3"    


##-------- Genes down regulated in at least one other cell types

NKdownSet <- rbind(
  DEGs_NK_CD4[DEGs_NK_CD4$logFC < -8,],
  DEGs_NK_CD8[DEGs_NK_CD8$logFC < -7,],
  DEGs_NK_B[DEGs_NK_B$logFC < -11,],
  DEGs_NK_Mono[DEGs_NK_Mono$logFC < -11,],
  DEGs_NK_Neut[DEGs_NK_Mono$logFC < -12,]
)

NKdownSet <- NKdownSet %>% group_by(SYMBOL) %>%  add_tally() %>%  data.frame()
NKdownSetGenes <- NKdownSet$SYMBOL
names(NKdownSetGenes) <- NKdownSet$Comparison
## 7 genes are duplicated:"EDAR"     "MAL"   "LRRN3" "CACNA1I"     "CD5"     "NOG"    "CD28"
NKdownSetGenes <- NKdownSetGenes[!duplicated(NKdownSetGenes)]
# NKdownSetUnique <- NKdownSet[!duplicated(NKdownSet$SYMBOL), ]
# [1] "EDAR"       "CCR4"       "TSHZ2"      "FBLN7"      "MDS2"      
# [6] "MAL"        "LRRN3"      "ST6GALNAC1" "C14orf132"  "ADTRP"     
# [11] "CACNA1I"    "CD4"        "CD5"        "NOG"        "CD28"      
# [16] "SIT1"       "KCNN4"      "THEMIS"     "LINC01550"  "NPAS2"     
# [21] "NRCAM"      "CD248"      "CACHD1"     "SFRP5"      "CA6"       
# [26] "MMP28"      "REG4"       "CD8A"       "CD8B"       "CPA5"      
# [31] "FCRL1"      "PLD4"       "COL19A1"    "EBF1"       "FCER2"     
# [36] "BLNK"       "PAX5"       "PLEKHG1"    "KLHL14"     "FCRL2"     
# [41] "TCL1A"      "FCRLA"      "CD22"       "VCAN"       "CYBB"      
# [46] "CYP1B1"     "MPEG1"      "FCN1"       "IL1B"       "LGALS2"    
# [51] "LYZ"        "ASGR2"      "GASK1B"     "CPVL"       "CLEC7A"    
# [56] "SIGLEC1"    "TGFBI"      "CD1D"       "CD163"      "CD86"      
# [61] "CD36"       "MAFB"       "S1PR1"      "FBXO4"      "IL1RN"     
# [66] "KCNJ14"     "ACSS2"      "SLC1A7"     "CLMN"       "NAA15"    

##------- Genes down regulated in either CD4 or CD8 cell types

# nonNK_CDs <- unique(
#   DEGs_NK_CD4$SYMBOL[DEGs_NK_CD4$logFC < -7], 
#   DEGs_NK_CD8$SYMBOL[DEGs_NK_CD8$logFC < -7])
# 
# nn2 <- c(nonNK_CDs[! nonNK_CDs  %in% NKdownSetGenes],
#   NKdownSetGenes)


GeneList_GSE60424 <- list(
  DEGs_NK_nonNK = DEGs_eBayes,
  DEGs_NK_CD4 = DEGs_NK_CD4, 
  DEGs_NK_CD8 = DEGs_NK_CD8, 
  DEGs_NK_B = DEGs_NK_B, 
  DEGs_NK_Mono = DEGs_NK_Mono, 
  DEGs_NK_Neut = DEGs_NK_Neut, 
  nonNK_Genes_4CellTypes = nonNK_Genes5,
  Genes_ExtremeDownNK = NKdownSetGenes 
  )

names(GeneList_GSE60424)[1] <- "DEGs_NK_nonNK"

saveRDS(GeneList_GSE60424, "./output/DEGs_List_NK_GSE60424.RDS")

##============ 

logRPKM <- rpkm(dgeFiltNorm, log = T)
dd0 <- merge(genes, logRPKM, by.x = "ENTREZID", by.y = "row.names")

dd0 <- dd0[!duplicated(dd0$SYMBOL), ]
dd0 <- dd0[complete.cases(dd0$SYMBOL), ]
rownames(dd0) <- dd0$SYMBOL
dd0 <- dd0[, 4:ncol(dd0)]

##----- check comDEGs
comDEGs <- intersect(GeneList_GSE107011$DEGs_NK_nonNK$Symbols[GeneList_GSE107011$DEGs_NK_nonNK$logFC > 0], 
  GeneList_GSE60424$DEGs_NK_nonNK$SYMBOL[GeneList_GSE60424$DEGs_NK_nonNK$DEG_type =="Treat_Str" & GeneList_GSE60424$DEGs_NK_nonNK$logFC > 0])

samCol <- c("Green3",  brewer.pal(9, "Purples")[c(3, 4)], "blue", "lightblue", "black")
names(samCol) <- names(table(dgeFiltNorm$samples$celltype_s))

samAnnot <- 
  ComplexHeatmap::columnAnnotation(df = dgeFiltNorm$samples[, "celltype_s", drop = F], col =  list(celltype_s = samCol))

pdf("./figure/Heatmap_DEGsCom_Up_NK_GSE60424.pdf", height = 30, width = 10)
p <- ComplexHeatmap::Heatmap(dd0[comDEGs, ], top_annotation = samAnnot)
draw(p, annotation_legend_side = "right", heatmap_legend_side = "right")
dev.off()


gg <- NKDownData$NKDown[NKDownData$NKDown  %in% rownames(dd0)]
aa <- NKDownData[gg, ]
d0 <- dd0[gg, ]

gAnnot <-
  ComplexHeatmap::rowAnnotation(df = aa[, "Data", drop = F], col =  list(Data = dataCols))


pdf("./figure/Heatmap_DEGs_Down_NK_GSE60424.pdf", height = 18, width = 10)
p <- ComplexHeatmap::Heatmap(d0, left_annotation = gAnnot, top_annotation = samAnnot)
draw(p, annotation_legend_side = "right", heatmap_legend_side = "right")
dev.off()

















