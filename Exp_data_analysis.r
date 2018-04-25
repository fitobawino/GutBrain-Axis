# Pipeline based on Lun et al 2015, It's DE-licious: a recipe for differential expression analyses of RNA-seq...
# Libraries:
#    - edgeR
#    - org.Hs.eg.db
#
# Files:
#    - "annotFile_path" is a cvs file with the following columns: RNAseq_sample_name and sub_structure (group names). Add a mean line.
#    - "rawdata_path" is a cvs file without column names. The first column are the gene names, then the columns of sequencing counts according to RNAseq_sample_name, and the last column is the mean expression value of each gene.
#    - "Genes_path" is a cvs file with the following columns: gene_symbol and entrez_id
#
#

rm(list = ls())
#1. Set variables
annotFile_path <- "/media/fito/Databases/FW_test/Databases/Allen_Institute/RNA-seq/rnaseq_donor10021/SampleAnnot+mean.csv"
rawdata_path   <- "/media/fito/Databases/FW_test/Databases/Allen_Institute/RNA-seq/rnaseq_donor10021/RNAseqCounts+mean.csv"
genes_path     <- "/media/fito/Databases/FW_test/Databases/Allen_Institute/RNA-seq/rnaseq_donor10021/Genes.csv"
save_path      <- "/media/fito/Databases/FW_test/Results/Expression/"
dir.create(save_path)

#2. Load libraries
library(edgeR)
library(org.Hs.eg.db)

#3. Load the data
annotFile <- read.delim(annotFile_path, header = TRUE, sep = ',')
rawdata <- read.delim(rawdata_path, header = FALSE, sep = ',', col.names = c("Genes", levels(annotFile$RNAseq_sample_name)))
genes <- read.delim(genes_path, header = TRUE, sep = ',')

#4. Create an object from a table of counts
y <- DGEList(counts  = rawdata[,2:length(rawdata)], 
             genes   = rawdata[,1], 
             samples = annotFile$RNAseq_sample_name,
             group   = annotFile$sub_structure,
             remove.zeros = TRUE)
dim(y)
m <- match(y$genes$genes, genes$gene_symbol)
y$genes$EntrezGene <- genes$entrez_id[m]

#Remove 'NAs', genes not annotated in entrez
y <- y[complete.cases(y$genes),]
dim(y)

# 5. Normalization
#Remove duplicates
o <- order(rowSums(y$counts), decreasing=TRUE)
d <- duplicated(y$genes)
y <- y[!d,]
dim(y)

# Filter low expressed and normalize
y$samples$lib.size <- colSums(y$counts)
keep <- rowSums(cpm(y) > 0.5) >= 2
y <- y[keep,]
summary(keep)
dim(y)

y <- calcNormFactors(y)
#plotMD(cpm(y, log=TRUE), column = 20) #Mean difference plot
#abline(h=0, col="red", lty=2, lwd=2)

#6. Exploration
table(y$samples$group)
cols <- as.numeric(y$samples$group)+2
barplot(colSums(y$counts), las = 2, col=cols, cex.names=0.5, cex.axis=0.8)
plotMDS(y, col=cols, labels = y$samples$group,  
        top = 500, 
        xlab = "Dim1", 
        ylab = "Dim2",
        #xlim = c(-1, 1),
        #ylim = c(-0.5, 0.5)
        )
        #main = "Multidimensional scaling plot",


#7. GLM approach of experiment design
design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- levels(y$samples$group)
head(design)
fit <- glmQLFit(y, design, robust=TRUE)
head(fit$coefficients)
plotQLDisp(fit)
summary(fit$df.prior)

#8. Compare all unique(y$samples$group) against Mean and perform ANOVA
con <- makeContrasts(
  AnG_i  = AnG_i - Mean,
  AnG_s  = AnG_s - Mean,
  Caudate= Caudate - Mean,
  CbCx   = CbCx - Mean,
  CgG    = CgG - Mean,
  FuG_i  = FuG_its - Mean,
  GP     = GP - Mean,
  GRe    = GRe - Mean,
  Insula = Insula - Mean,
  ITG    = ITG - Mean,
  MFG    = MFG - Mean,
  MTG    = MTG - Mean,
  OrbGyri= OrbGyri - Mean,
  orIFG  = orIFG - Mean,
  PCLa_i = PCLa_i - Mean,
  PCLa_s = PCLa_s - Mean,
  Pcu    = Pcu - Mean,
  pest_V2= pest_V2 - Mean,
  PHG    = PHG - Mean,
  PoG_cs = PoG_cs - Mean,
  PoG_l  = PoG_l - Mean,
  PrG    = PrG - Mean,
  Putame = Putamen - Mean,
  SFG_l  = SFG_l - Mean,
  SFG_m  = SFG_m - Mean,
  SMG_i  = SMG_i - Mean,
  SPL    = SPL - Mean,
  STG    = STG - Mean,
  str_V1 = str_V1 - Mean,
  levels = design)
anov <- glmQLFTest(fit, contrast = con)
topTags(anov)

#9. Save table to file
out_var <- data.frame(row.names = anov$genes$EntrezGene, anov$table)
out_var$FDR <- p.adjust(out_var$PValue, method = "bonferroni")
out_var3 <- subset(out_var, FDR <= 0.05) # Si el FDR es menor o igual a 0.05

write.table(out_var, file = paste0(save_path,"TableS2.txt"), append = TRUE, sep = '\t')
dim(out_var)

#10. Save differentially expressed genes per sub_structure
out_logic <- data.frame(row.names = rownames(out_var))
for (col in names(out_var)[1:29]){
  #out_logic[col] <- abs(out_var[col]) > apply(out_var[col],2,sd)  #>|STD|
  out_logic[col] <- abs(out_var[col]) > 1.5 # > 1.5 FC
}
for (col in names(out_var)[1:29]){
  write(rownames(out_logic)[out_logic[[col]]], 
        file = paste(save_path, "Genes_", col, ".txt", sep = ''), 
        append = TRUE)
}
