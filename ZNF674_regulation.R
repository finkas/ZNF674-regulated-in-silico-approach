library(biomaRt)
library(ggplot2)
library(ggpubr)
library(GenomicRanges)
library(openxlsx)
library(plyr)
library(dplyr)

###################################################################################
################### ChIP-Seq: Filtering of KRAB-regulated Genes ###################
###################################################################################
#Load the BED files, downloaded from 
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78099
#All other KRAB-domain proteins merged in LINUX using BEDtools: 
All_KRAB <- read.table("merged.bed", col.names = c("chr", "start", "end","peak","Intensity", "Strand", "Gene"), sep="\t")

#Extract the maximal Intensity of each ZFP (optional for karyoploteR) and
# annotate it to the data frame
max_intens <- ddply(All_KRAB, .(Gene), summarize, max = max(Intensity, na.rm = T))
All_KRAB$max_intens <- ifelse(All_KRAB$Gene %in% max_intens$Gene, 
                              max_intens$max[match(All_KRAB$Gene, 	
                              max_intens$Gene)], NA)

#Calculating the relative peak intensities for each ZFP by normalizing on the maximum 
All_KRAB$rel_intens <- ifelse(
  All_KRAB$Gene %in% max_intens$Gene, All_KRAB$Intensity/max_intens$max[match(All_KRAB$Gene, max_intens$Gene)], 
  NA)
ZNF674_KRAB <- subset.data.frame(All_KRAB, Gene == "ZNF674")

#Transform peaks to a GRanges object 
KRAB_ranges <- GRanges(seqnames = ZNF674_KRAB$chr, IRanges(
  start = ZNF674_KRAB$start, end = ZNF674_KRAB$end), 
  Intensity = ZNF674_KRAB$Intensity, 
  Gene = ZNF674_KRAB$Gene, 
  rel_intens = ZNF674_KRAB$rel_intens)

#Import GeneHancer Track set from UCSC: 
#Chromosome size: 
# https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
#GRCh37/hg19, Group: regulation, Track: GeneHncer, 
# Table: GH Interactions (DE) (geneHancerInteractionsDpubleElite)
# 'chrom', 'chromStart', 'chromEnd', 'name', 'score', 
# 'geneHancerChrom', 'geneHancerStart', 'geneHancerEnd', 
# 'geneHancerIdentifier', 'geneChrom', 'geneStart', 'geneEnd', 
# 'geneName'
GeneHancer_DE <- read.table("GeneHancer_UCSC_DE.txt", header = T, sep = "\t")

#Transform GeneHancer_DE to a GRanges object
GeneHanDE_ranges <- GRanges(seqnames = GeneHancer_DE$geneHancerChrom, IRanges(
  start = GeneHancer_DE$geneHancerStart, end = GeneHancer_DE$geneHancerEnd), 
  Gene = GeneHancer_DE$geneName, RegEl = GeneHancer_DE$geneHancerIdentifier)

#Find overlaps with the ZNF674 ranges
overlaps_KRAB <- findOverlaps(KRAB_ranges, GeneHanDE_ranges)
overlap_dataKRAB <- data.frame(
  # Index of overlapping ChIP peaks
  KRAB_index = queryHits(overlaps_KRAB),
  # Index of overlapping GeneHancer regions   
  GeneHanDE = subjectHits(overlaps_KRAB),
  # Corresponding gene   
  Gene = GeneHanDE_ranges$Gene[subjectHits(overlaps_KRAB)], 
  # KRAB Regulator
  KRAB_Gene = KRAB_ranges$Gene[queryHits(overlaps_KRAB)], 
  # Extract ChIP-Seq intensity
  RegEl = GeneHancer_DE$geneHancerIdentifier[subjectHits(overlaps_KRAB)],  
  Intensity = KRAB_ranges$Intensity[queryHits(overlaps_KRAB)], 
  Chr = KRAB_ranges@seqnames[queryHits(overlaps_KRAB)], 
  Start = KRAB_ranges@ranges@start[queryHits(overlaps_KRAB)], 
  End = KRAB_ranges@ranges@start[queryHits(overlaps_KRAB)]+
    KRAB_ranges@ranges@width[queryHits(overlaps_KRAB)]
)

#Annotate the ENSG's to the overlapping data.
# Default is the hg19 ENSEMBL version. For genes not known, 
# the hg38 version or the manual curation was applied
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl",
                      host = "https://grch37.ensembl.org")
genENSEMBL <- getBM(attributes=c("hgnc_symbol", "ensembl_gene_id", 
                                 "external_synonym"), mart=ensembl)
ensmebl_hg38 <- useEnsembl(biomart="genes", dataset = "hsapiens_gene_ensembl")
genENSEMBL38 <- getBM(attributes=c("hgnc_symbol", "ensembl_gene_id", 
                                   "external_synonym"), mart=ensmebl_hg38)

#use hg19
overlap_dataKRAB$ensembl_gene_id19 <- ifelse(
  overlap_dataKRAB$Gene %in% genENSEMBL$hgnc_symbol, 
  genENSEMBL$ensembl_gene_id[match(overlap_dataKRAB$Gene, genENSEMBL$hgnc_symbol)], 
  ifelse(overlap_dataKRAB$Gene %in% genENSEMBL$external_synonym, 
         genENSEMBL$ensembl_gene_id[match(overlap_dataKRAB$Gene, 
         genENSEMBL$external_synonym)],
         ifelse(grepl("^ENSG", overlap_dataKRAB$Gene), overlap_dataKRAB$Gene, NA)))

#use hg38
overlap_dataKRAB$ensembl_gene_id38 <- ifelse(
  overlap_dataKRAB$Gene %in% genENSEMBL38$hgnc_symbol, 
  genENSEMBL38$ensembl_gene_id[match(overlap_dataKRAB$Gene, 
  genENSEMBL38$hgnc_symbol)], 
  ifelse(overlap_dataKRAB$Gene %in% genENSEMBL38$external_synonym, 
         genENSEMBL38$ensembl_gene_id[match(overlap_dataKRAB$Gene, 
         genENSEMBL38$external_synonym)],
         ifelse(grepl("^ENSG", overlap_dataKRAB$Gene), overlap_dataKRAB$Gene, NA)))

#use hg19, if not availabe use hg38
overlap_dataKRAB$ensembl_gene_id <- ifelse(
  !is.na(overlap_dataKRAB$ensembl_gene_id19) &
  !is.na(overlap_dataKRAB$ensembl_gene_id38) & 
    overlap_dataKRAB$ensembl_gene_id19 == overlap_dataKRAB$ensembl_gene_id38,
  # they match
  overlap_dataKRAB$ensembl_gene_id38,  
  # only hg38 present
  ifelse(
    is.na(overlap_dataKRAB$ensembl_gene_id19) & 
    !is.na(overlap_dataKRAB$ensembl_gene_id38),overlap_dataKRAB$ensembl_gene_id38,  
    # catch both different & only hg19 present
    ifelse(
      !is.na(overlap_dataKRAB$ensembl_gene_id19), 
      # fallback to hg19 
      overlap_dataKRAB$ensembl_gene_id19, NA_character_  # both NA
    )))

#manually annotate remaining NA's in ensembl_gene_id
NAs <- data.frame(Gene = 	unique(overlap_dataKRAB$Gene[is.na(overlap_dataKRAB$ensembl_gene_id)]))
NAs <- read.xlsx("NAs2.xlsx")
overlap_dataKRAB$ensembl_gene_id <- ifelse(is.na(overlap_dataKRAB$ensembl_gene_id), 
                 NAs$ensembl_gene_id[match(overlap_dataKRAB$Gene, NAs$Gene)],
                                           overlap_dataKRAB$ensembl_gene_id)

overlap_dataKRAB$ensembl_gene_id <- ifelse(is.na(overlap_dataKRAB$ensembl_gene_id), 
                 overlap_dataKRAB$Gene, overlap_dataKRAB$ensembl_gene_id)



###################################################################################
######### Import and Filter for Overexpression in Sample and Male Controls ########
###################################################################################

#Import the data sets extracted from the OUTRIDER data set
sample <- read.xlsx("BGRNA247844.xlsx")
sample$sample <- sample$sampleID
sample$sampleID <- NULL
load('all_unaffMale_Rd2.RData')
all_unaffMale$enseID <- all_unaffMale$geneID
all_unaffMale$geneID <- NULL

#Limit the data sets to the significantly overexpressed exons
sample_over <- subset.data.frame(sample, zScore > 0 & pVal < 0.05)
all_unaffMaleOver <- subset.data.frame(all_unaffMale, zScore > 0 & pVal < 0.05)

#Use the OUTRIDER Annotation table for ensg Annotation
genEx <- read.table("ENSE_to_Gene.tab", 
                    col.names=c("Chr", "start", "end", "enseID", "geneName"), sep="\t")
ENSGs <- getBM(attributes = c("ensembl_exon_id", "ensembl_gene_id", 
                              "external_gene_name"), 
               filters="ensembl_exon_id",values = unique(genEx$enseID),
               mart=ensembl)
sample_over$ensgID <- ifelse(sample_over$enseID %in% ENSGs$ensembl_exon_id, 
                             ENSGs$ensembl_gene_id[match(sample_over$enseID,
                             ENSGs$ensembl_exon_id)], NA)
all_unaffMaleOver$ensgID <- ifelse(all_unaffMaleOver$enseID %in% 
                             ENSGs$ensembl_exon_id, ENSGs$ensembl_gene_id[match(
                             all_unaffMaleOver$enseID,ENSGs$ensembl_exon_id)], NA)

#Annotate KRAB-ZFP regulation in the sample
sample_over$KRABreg <- ifelse(sample_over$ensgID %in% overlap_dataKRAB$ensembl_gene_id, TRUE, FALSE)
sample_overKRAB <- subset.data.frame(sample_over, KRABreg == TRUE)

###################################################################################
## Filtering overexpressed KRAB-regulated Genes/Exons in Sample and MaleControls ##
###################################################################################
#Filtering and Statistics
sample_overKRAB$geneName <- ifelse(sample_overKRAB$ensgID %in% 
                     ENSGs$ensembl_gene_id, ENSGs$external_gene_name[match(
                     sample_overKRAB$ensgID, ENSGs$ensembl_gene_id)], NA)
sample_overKRABsum <- sample_overKRAB %>%
  group_by(geneName, ensgID) %>%
  tally()

#annotate the aberrant ones in the unaffected male, considering the 
# individual with the highest number of aberrant exons
all_unaffMaleOversum <- unique(all_unaffMaleOver[,c("sample", "enseID", "ensgID")])
all_unaffMaleOversum <- all_unaffMaleOversum %>%
  group_by(sample, ensgID) %>%
  tally()
all_unaffMaleOversum <- ddply(all_unaffMaleOversum, .(ensgID), summarize, 
                              max=max(n, na.rm = TRUE))
sample_overKRABsum$nCtrl <- ifelse(sample_overKRABsum$ensgID %in% 	     
       all_unaffMaleOversum$ensgID, all_unaffMaleOversum$max[match(
       sample_overKRABsum$ensgID, all_unaffMaleOversum$ensgID)], 0)

#annotate the total number of RefSeq Exons
ENSGssum <- ENSGs %>%
  group_by(ensembl_gene_id, external_gene_name) %>%
  tally()
sample_overKRABsum$Total <- ifelse(sample_overKRABsum$ensgID %in%
       ENSGssum$ensembl_gene_id, ENSGssum$n[match(sample_overKRABsum$ensgID,
       ENSGssum$ensembl_gene_id)], NA)

#Filtering
sample_overRes <- subset.data.frame(sample_overKRABsum, n/Total > 2/3 & n>nCtrl)
