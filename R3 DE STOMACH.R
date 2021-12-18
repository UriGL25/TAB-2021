load(file = "C:/Users/oriol/Documents/Uni/4t Curs/Primer semestre/TAB/Definitivas/P4/stomach/rse_gene_stomach.Rdata")

stage <- rse_gene$gdc_cases.diagnoses.tumor_stage

# id list of early tumours
ids.early <-grep(paste("stage i$", "stage ia$","stage ib$","stage ic$", "stage ii$", "stage iia$","stage iib$", "stage iic$",sep="|"), stage)

# id list of late tumours
ids.late <-grep(paste("stage iii$", "stage iiia$", "stage iiib$","stage iiic$", "stage iv$", "stage iva$","stage ivb$", "stage ivc$",sep="|"), stage)

# create an empty column named GROUP
colData(rse_gene)$GROUP <-rep(NA, ncol(rse_gene))
# add early for those patients with tumours at stages i-ii
colData(rse_gene)$GROUP[ids.early] <- "early"
# add late for those patients with tumours at stages iii-iv
colData(rse_gene)$GROUP[ids.late] <- "late"

dataGroups <- data.frame(Stage = rse_gene$gdc_cases.diagnoses.tumor_stage, Group = rse_gene$GROUP)
head(dataGroups)
ggplot(data = dataGroups, mapping = aes(x = Stage, fill = Group)) +
  labs(x = 'Stage', y = 'Number of patients') +
  theme_bw() +
  geom_bar()

naData <- is.na(rse_gene$GROUP)
table(naData)

dim(rse_gene)
# Remove NA
rse_gene <- rse_gene[, !naData]
dim(rse_gene)

table(rse_gene$GROUP)

counts <- assay(rse_gene, "counts")
counts[1:5, 1:2]

phenotype <- colData(rse_gene)
phenotype[1:5, 1:2]

identical(colnames(counts), rownames(phenotype))

annotation <- rowData(rse_gene)
head(annotation)

##NORMALIZATION##

maPlot(counts[,1],
       counts[,2],
       pch=19,
       cex=.5,
       ylim=c(-8,8),
       allCol="darkgray",
       lowess=TRUE,
       xlab=expression(A == log[2] (sqrt(S1/N %.% S2/N))),
       ylab=expression(M == log[2](S1/N)-log[2](S2/N)))
grid(col="black")
title("Raw data")

geneLength <- annotation$bp_length

# RPKM normalization method
counts.rpkm <- t(t(counts/geneLength*1000)/colSums(counts)*1e6)

counts.rpkm[1:5,1:2]

#RPMK#

maPlot(counts.rpkm[,1],
       counts.rpkm[,2],
       pch=19,
       cex=.5,
       ylim=c(-8,8),
       allCol="darkgray",
       lowess=TRUE,
       xlab=expression( A==log[2] (sqrt(S1/N%.%S2/N))),
       ylab=expression(M==log[2](S1/N)-log[2](S2/N)))
grid(col="black")
title("RPKM")

#TMM#
counts.tmm <- normalizeCounts(counts, method = "TMM")
counts.tmm[1:5,1:2]
maPlot(counts.tmm[,1],
       counts.tmm[,2],
       pch=19, cex=.5,
       ylim=c(-8,8),
       allCol="darkgray",
       lowess=TRUE,
       xlab=expression( A==log[2] (sqrt(S1/N%.%S2/N))),
       ylab=expression(M==log[2](S1/N)-log[2](S2/N)))
grid(col="black")
title("TMM")

##Differential analysis##

# stage of each patient
pheno.stage <- subset(phenotype, select=GROUP)
# recreate the counts in a new matrix
counts.adj <- matrix((as.vector(as.integer(counts))), nrow=nrow(counts), ncol=ncol(counts))
rownames(counts.adj) = rownames(counts)
colnames(counts.adj) = colnames(counts)
# check information
identical(colnames(counts.adj), rownames(pheno.stage))
# transform the group variable to factor
pheno.stage$GROUP <- as.factor(pheno.stage$GROUP)
# create the DESeqDataSet input
DEs <- DESeqDataSetFromMatrix(countData = counts.adj,
                              colData = pheno.stage,
                              design = ~ GROUP)
# differential expression analysis
dds <- DESeq(DEs)
# results extracts a result table from a DESEq analysis 
res <- results(dds, pAdjustMethod = "fdr")
head(res)
plotMA(res, ylim=c(-20,20), main='Differentially expressed genes in early vs late pleura cancer')

#Filters
res.padj <- res[(res$padj < 0.001 & !is.na(res$padj)),]
nrow(res.padj)

res.padj_4fold <- res.padj[abs(res.padj$log2FoldChange) > log2(10),] 
nrow(res.padj_4fold)

plotMA(res.padj_4fold, ylim=c(-10,10), main='Most differentially expressed genes')

##POST ANALYSIS##

# we first create the dataframe with the results of the DE
resDF <- as.data.frame(res)
# new column with the gene names
resDF$gene <- rownames(resDF)
# new column with TRUE/FALSE representing if a gene is deferentially expressed or not
resDF$filter <- abs(res$log2FoldChange) > log2(10) & res$padj < 0.001 
table(resDF$filter)

#Volcano#
volcano <- ggplot(resDF, aes(x = log2FoldChange, y = -log10(padj), color=filter)) +
  geom_point(size= 1)
volcano <- volcano + labs (x= 'Effect size: log2(fold-change)', y="-log10(adjusted p-value)", title="Volcano plot") +
  scale_color_manual(values = c("gray", "black")) +
  theme_bw() +
  theme(legend.position = "none")
options(ggrepel.max.overlaps = Inf)
volcano <- volcano + geom_label_repel(data=resDF[abs(resDF$log2FoldChange) > log2(10) & resDF$padj < 0.001 ,], aes(label = as.factor(gene)), alpha = 0.7, size = 2, force = 1.3)  

volcano

##ENRICHMENT ANALYSIS##

# we save a list of DE genes, and we remove the last part of the id
deGenes <- gsub("\\..*", "",rownames(res.padj_4fold))
head(deGenes)
biomart_uses <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
biomart_uses

Entrez_BM <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                   filter="ensembl_gene_id",
                   values=deGenes,
                   mart=biomart_uses,
                   uniqueRows=TRUE)

dim(Entrez_BM)

Entrez_ids <- as.character(na.omit(Entrez_BM$entrezgene_id))

universe <- mappedkeys(org.Hs.egGO)
count.mappedkeys(org.Hs.egGO)

# Set the threshold for the p value to 5%
GOtest <- new("GOHyperGParams",
              geneIds = Entrez_ids,
              universeGeneIds=universe,
              annotation = "org.Hs.eg.db",
              ontology = "BP",
              pvalueCutoff= 0.05,
              conditional = FALSE,
              testDirection = "over") 
# hypergeometric test
GOtestOver <- hyperGTest(GOtest)
GOtestOver

GOresult <- summary(GOtestOver)
head(GOresult)
dim(GOresult)

##WHICH GGPLOT COULD WE USE HERE?##


head(Entrez_BM)