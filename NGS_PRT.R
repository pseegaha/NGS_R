##Dataset packages install order
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19", version = "3.8")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene", version = "3.8") 

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("gdsfmt")
BiocManager::install("SNPRelate")

## NGS analysis of BSgenome.Hsapiens.UCSC.hg19 
library(BSgenome.Hsapiens.UCSC.hg19)
?BSgenome.Hsapiens.UCSC.hg19

genome1 <- BSgenome.Hsapiens.UCSC.hg19
seqlengths(genome1)
plot(seqlengths(genome1))
seqnames(genome1)
genome1$chr1  # same as genome1[["chr1"]]

for (seqname in seqnames(genome1)[1:2]) {
  cat("Checking sequence", seqname, "... ", length(genome1[[seqname]]))
  #seq <- genome1[[seqname]]
  #checkOnlyNsInGaps(seq)
  cat("\n")
}

length(names(genome1$upstream1000))
names(genome1$upstream1000)[1:4]
getSNPcount()
##plot
plot(getSNPcount(),seqlengths(genome1)[1:25])



## Accessing the data stored for TxDb.Hsapiens.UCSC.hg19.knownGene
## Get the locations and alleles of all SNPs on chromosome 22:
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
#available.SNPs()
ch22snps <- getSNPlocs("ch22")
dim(ch22snps)
colnames(ch22snps)
head(ch22snps)
## Convenience wrappers for loading the SNP data:
getSNPcount()
getSNPlocs(c("ch22", "chMT"), as.GRanges=TRUE)



## Extract SNP information for a set of rs ids:
myrsids <- c("rs2639606", "rs75264089", "rs73396229", "rs55871206",
             "rs10932221", "rs56219727", "rs73709730", "rs55838886",
             "rs3734153", "rs79381275", "rs75350930", "rs1516535")
rsidsToGRanges(myrsids)


## Inject the SNPs in hg19 (injectSNPs() "knows" how to map dbSNP chromosome names to UCSC names):
library(BSgenome.Hsapiens.UCSC.hg19)
Hsapiens
Hs2 <- injectSNPs(Hsapiens, "SNPlocs.Hsapiens.dbSNP.20101109")
Hs2
alphabetFrequency(unmasked(Hs2$chr22))
alphabetFrequency(unmasked(Hsapiens$chr22))
## Get the number of nucleotides that were modified by this injection:
neditAt(unmasked(Hs2$chr22), unmasked(Hsapiens$chr22))

## Quality Control(QC) stuff
## Note that dbSNP can assign distinct ids to SNPs located at the same position: rs ids are all distinct
any(duplicated(ch22snps$RefSNP_id)) 
# but some locations are repeated
any(duplicated(ch22snps$loc)) 
# sort by location
ch22snps <- ch22snps[order(ch22snps$loc), ] 
# 777
which(duplicated(ch22snps$loc))[1]  
# rs75258394 and rs78180314 have same locations and alleles
ch22snps[775:778, ]  
## Also note that not all SNP alleles are consistent with the hg19 genome i.e. the alleles reported for a given SNP are not always compatible with the nucleotide found at the SNP location in hg19.
ch1snps <- getSNPlocs("ch1")
all_alleles <- paste(ch1snps$alleles_as_ambig, collapse="")
# 1849438 SNPs on chr1
nchar(all_alleles)  
neditAt(all_alleles, unmasked(Hsapiens$chr1)[ch1snps$loc], fixed=FALSE)
## 3181 SNPs (0.17%) are inconsistent with hg19 chr1.let's check that no SNP falls in an assembly gap:
assemblygaps <- masks(Hsapiens$chr1)$AGAPS
# the assembly gaps
assemblygaps  
## Looping over the assembly gaps:
length(findOverlaps(ch1snps$loc, assemblygaps))  # 0

