knitr::opts_chunk$set(echo = TRUE, fig.align = "center", fig.width = 6, fig.height = 4)
library(ggplot2)
library(SNPassoc)
library(snpStats)
library(SNPRelate)
library(dplyr)
library(ggrepel)

colorectal.plink <-  read.plink(bed = "C:/Users/oriol/Documents/Uni/4t Curs/Primer semestre/TAB/Definitivas/P3/Project3/data/colorectal/colorectal.bed",
                                bim = "C:/Users/oriol/Documents/Uni/4t Curs/Primer semestre/TAB/Definitivas/P3/Project3/data/colorectal/colorectal.bim",
                                fam = "C:/Users/oriol/Documents/Uni/4t Curs/Primer semestre/TAB/Definitivas/P3/Project3/data/colorectal/colorectal.fam")
names(colorectal.plink)
colorectal.genotype <- colorectal.plink$genotypes
colorectal.genotype

individuals <- colorectal.plink$fam
head(individuals)

annotation <- colorectal.plink$map
head(annotation)

colorectal.phenotype <- read.delim("C:/Users/oriol/Documents/Uni/4t Curs/Primer semestre/TAB/Definitivas/P3/Project3/data/colorectal/colorectal.txt")
head(colorectal.phenotype)

rownames(colorectal.phenotype) <- colorectal.phenotype$id
head(colorectal.phenotype)

identical(rownames(colorectal.phenotype), rownames(colorectal.genotype))

genotype <- colorectal.genotype
phenotype <- colorectal.phenotype

info.snps <- col.summary(genotype)
head(info.snps)

controls <- phenotype$cascon == 0 & !is.na(phenotype$cascon)
genotype.controls <- genotype[controls, ]
info.controls <- col.summary(genotype.controls)

use <- info.snps$Call.rate > 0.95 &
  info.snps$MAF > 0.05 &
  abs(info.controls$z.HWE < 3.3)
mask.snps <- use & !is.na(use)
genotype.qc.snps <- genotype[, mask.snps]
genotype.qc.snps

annotation <- annotation[mask.snps, ]
genotype

genotype.qc.snps

##REPORT SNPS!##

sum(info.snps$Call.rate < 0.95, na.rm = TRUE)
sum(info.snps$MAF < 0.05, na.rm = TRUE)
sum(abs(info.controls$z.HWE > 3.3), na.rm = TRUE)
sum(!mask.snps) #Total of SNPs removed#

info.indv <- row.summary(genotype.qc.snps)
head(info.indv)

genotype.X <- genotype.qc.snps[,annotation$chromosome=="23" & !is.na(annotation$chromosome)]
info.X <- row.summary(genotype.X)
info.X$sex <- phenotype$sex
info.X$id <- phenotype$id

######PROBLEM~#######
ggplot(info.X, aes(y = Heterozygosity, x = id)) +
  geom_point(aes(color=sex), alpha = 0.7) + 
  labs(y = "Heterozygosity", x = "ID", color = "Gender") +
  theme_minimal() + scale_color_manual(values = c("#FFE882", "#4DC4CC"))

##WARNING POTSER PASSA ALGO##
sex.discrep <- (info.X$sex == "male" &
                info.X$Heterozygosity > 0.2) |
                (info.X$sex=="female" &
                 info.X$Heterozygosity < 0.2)
#####


MAF <- col.summary(genotype.qc.snps)$MAF
callmatrix <- !is.na(genotype.qc.snps)
hetExp <- callmatrix %*% (2*MAF*(1-MAF))
hetObs <- with(info.indv,
               Heterozygosity*(ncol(genotype.qc.snps))*Call.rate)
info.indv$hetF <- 1 - (hetObs/hetExp)
head(info.indv)

ggplot(info.indv, aes(x = 1:nrow(info.indv), y = hetF)) +
  geom_point(aes(color = abs(hetF) > 0.1)) +
  geom_hline(yintercept = 0.1, linetype = "dashed") + 
  geom_hline(yintercept = -0.1, linetype = "dashed") + 
  labs(y = "F-Heterozygosity", x = "ID", color = "F-heterozigosity > ±0.1") +
  theme_minimal() + scale_color_manual(values = c("#4DC4CC", "#582602"))

snpgdsBED2GDS("C:/Users/oriol/Documents/Uni/4t Curs/Primer semestre/TAB/Definitivas/P3/Project3/data/colorectal/colorectal.bed",
              "C:/Users/oriol/Documents/Uni/4t Curs/Primer semestre/TAB/Definitivas/P3/Project3/data/colorectal/colorectal.fam",
              "C:/Users/oriol/Documents/Uni/4t Curs/Primer semestre/TAB/Definitivas/P3/Project3/data/colorectal/colorectal.bim",
              out="colorectalGDS"
)

genofile <- snpgdsOpen("colorectalGDS")
set.seed(12345)
snps.qc <- colnames(genotype.qc.snps)
snp.prune <- snpgdsLDpruning(genofile,ld.threshold = 0.2, snp.id=snps.qc)

snps.ibd <- unlist(snp.prune, use.names=FALSE)
ibd <- snpgdsIBDMoM(genofile, kinship = TRUE,
                    snp.id = snps.ibd,
                    num.thread = 1)

ibd.kin <- snpgdsIBDSelection(ibd)
head(ibd.kin)

ibd.kin.thres <- subset(ibd.kin, kinship > 0.1)
head(ibd.kin.thres)


##INDIVIDUALs REMOVED##
ids.rel <- related(ibd.kin.thres)
ids.rel

use <- info.indv$Call.rate > 0.95 &
  abs(info.indv$hetF) < 0.1 &     # or info.inv$Heterozygosity < 0.32
  !sex.discrep &
  !rownames(info.indv)%in%ids.rel
mask.indiv <- use & !is.na(use)
genotype.qc <- genotype.qc.snps[mask.indiv, ]

phenotype.qc <- colorectal.phenotype[mask.indiv, ]
identical(rownames(phenotype.qc), rownames(genotype.qc))

##TOTAL OF INDIVIDUALS WE HAD VS WE END UP KEEPING##
dim(phenotype)
dim(phenotype.qc)

##HOW MANY ARE REMOVED BC OF BAD CALL RATE and for hetero problems (2)##

sum(info.indv$Call.rate < 0.95)
sum(abs(info.indv$hetF)>0.1)

# Number of individuals removed for sex discrepancies
sum(sex.discrep)
# Number of individuals removed to be related with others
length(ids.rel)
# The total number of individuals that do not pass QC
sum(!mask.indiv)

#######THE GWAS ANALYSIS######


gwas <- single.snp.tests(cascon, data = phenotype.qc,
                         snp.data=genotype.qc)
 
gwasStats <- data.frame(SNP=annotation$snp.name, 
                        CHR=annotation$chromosome,
                        BP=annotation$position,
                        P=p.value(gwas, 1))

gwasStats <- subset(gwasStats, !is.na(CHR) & !is.na(P) & CHR!=24 & CHR!=25)
head(gwasStats)

gwas.adj <- snp.rhs.tests(cascon ~ smoke,  data = phenotype.qc,
                          snp.data=genotype.qc, family = "Gaussian")

gwas.adj.Stats <- data.frame(SNP=annotation$snp.name, 
                             CHR=annotation$chromosome,
                             BP=annotation$position,
                             P=p.value(gwas.adj))

gwas.adj.Stats <- subset(gwas.adj.Stats, !is.na(CHR) & !is.na(P) & CHR!=24 & CHR!=25)
head(gwas.adj.Stats)

head(gwasStats)

##MANHATTAN!!!##

nCHR <- length(unique(gwasStats$CHR))
gwasStats$BPcum <- NA
s <- 0
nbp <- c()
for (i in unique(gwasStats$CHR)){
  nbp[i] <- max(gwasStats[gwasStats$CHR == i,]$BP)
  gwasStats[gwasStats$CHR == i,"BPcum"] <- gwasStats[gwasStats$CHR == i,"BP"] + s
  s <- s + nbp[i]
}

axisdf <- gwasStats %>%
  group_by(CHR) %>%
  summarize(center=(max(BPcum) + min(BPcum))/2)

significance <- 1e-04
genomewideline <- 5e-08

##GRAPHIC!!##

manhattanPlot <- ggplot(gwasStats, aes(x = BPcum, y = -log10(P))) +
  geom_point(aes(color=as.factor(CHR)))
manhattanPlot


manhattanPlot <- manhattanPlot +
  scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center)
manhattanPlot

manhattanPlot <- manhattanPlot + labs(x= 'Chromosome', y='p-value (log10)', title='Obesity GWAS study')
manhattanPlot

manhattanPlot <- manhattanPlot +
  geom_hline(yintercept = -log10(significance))
manhattanPlot

manhattanPlot <- manhattanPlot + theme(legend.position="none")
manhattanPlot

mypalette <- c("#E2709A", "#CB4577", "#BD215B", "#970F42", "#75002B")
manhattanPlot <- manhattanPlot +
  scale_color_manual(values = rep(mypalette, length(unique(gwasStats$CHR))))
manhattanPlot

##SNPS ASSOCIATED TO THE DISEASE:###
manhattanPlot <- manhattanPlot +
  geom_label_repel(data=gwasStats[gwasStats$P<significance,], aes(label=as.factor(SNP), alpha=0.7), size=5, force=1.3)
manhattanPlot

