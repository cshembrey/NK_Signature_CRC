
## Xu et al DE analysis: GSE107011
# Note, data is in the TPM format.

dataPath <- "/Users/foroutanm/Documents/data/blood/GEO/GSE107011_RNAseq_bloodCells/"
tpm <-
  read.table(
    paste0(dataPath, "GSE107011_Processed_data_TPM.txt"),
    header = T,
    sep = "\t",
    stringsAsFactors = F,
    check.names = F
  )
head(tpm)

##------------------- Map IDs and generate logTPM
tpm$ENSG <- sapply(tpm[, 1], function(x) unlist(strsplit(x, "[.]"))[[1]])
### remove 45 duplicated names
tpm <- tpm[!duplicated(tpm$ENSG), ]

library(org.Hs.eg.db)               
tpm$GeneIDs = mapIds(
  org.Hs.eg.db,                     
  keys = tpm$ENSG,          
  keytype = "ENSEMBL",              
  column = "ENTREZID",              
  multiVals = "asNA"                
)                                 

## remove rows that do not map to a gene ID
tpm <- tpm[complete.cases(tpm$GeneIDs), ]
tpm <- tpm[ ! duplicated(tpm$GeneIDs), ]

dd <- as.matrix(tpm[, 2:(ncol(tpm)-2)])
row.names(dd) <- tpm$GeneIDs

## Remove PBMC samples
dd <- dd[ , ! grepl("PBMC", colnames(dd))]

logTPM <- log2(dd+1)

write.table(
  logTPM,
  "logTPM_GSE107011_bloodCells_2019.txt",
  row.names = T,
  sep = "\t"
)

##---------------- filter 
kp <- rowSums(logTPM > 2) >= 4
summary(kp)

logTPMFilt <- logTPM[kp, ]

##---------------- annotation data:
ann <- data.frame(
  Samples = colnames(dd),
  CellType = sapply(colnames(dd), function(x)
    unlist(strsplit(x, "[_]"))[[2]])
)

ann$Group[grepl("NK", ann$CellType)] <- "NK"
ann$Group[! grepl("NK", ann$CellType)] <- "nonNK"


##----------------- Do DE analysis through limma-trend
## eBayes with trend=TRUE 
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
  "DEStat_NK_nonNK_GSE107011.txt", row.names = F, sep = "\t")

DEGs <- DEstats[DEstats$GeneIDs %in% row.names(dt)[! dt[, 1] == 0], ]

write.table(DEGs, 
  "DEGs_NK_nonNK_GSE107011.txt", row.names = F, sep = "\t")

##--- only genes with positve log FC

DEGs_pos <- DEGs[DEGs$logFC > 0, ]

write.table(DEGs_pos, 
  "DEGs_positiveLogFC_NK_nonNK_GSE107011.txt", row.names = F, sep = "\t")

