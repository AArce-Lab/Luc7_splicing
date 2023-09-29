# This script is intended to show how mRNA, spliced and unspliced reads were counted
# It is not intended to be reused directly.
# Count tables are provided in the /data folder.
library("Rsubread")

bamfiles <- c("C:/SMD/Luc7/UMI_Removed/WT_1_umi_1_trimmed.fq.gz.sam.bam_sorted.bam_PCR_removed.bam",
							"C:/SMD/Luc7/UMI_Removed/WT_2_umi_1_trimmed.fq.gz.sam.bam_sorted.bam_PCR_removed.bam",
							"C:/SMD/Luc7/UMI_Removed/WT_3_umi_1_trimmed.fq.gz.sam.bam_sorted.bam_PCR_removed.bam",
							"C:/SMD/Luc7/UMI_Removed/N31_1_umi_1_trimmed.fq.gz.sam.bam_sorted.bam_PCR_removed.bam",
							"C:/SMD/Luc7/UMI_Removed/N31_2_umi_1_trimmed.fq.gz.sam.bam_sorted.bam_PCR_removed.bam",
							"C:/SMD/Luc7/UMI_Removed/N31_4_umi_1_trimmed.fq.gz.sam.bam_sorted.bam_PCR_removed.bam",
							"C:/SMD/Luc7/UMI_Removed/ZNF2_1_umi_1_trimmed.fq.gz.sam.bam_sorted.bam_PCR_removed.bam",
							"C:/SMD/Luc7/UMI_Removed/ZNF2_2_umi_1_trimmed.fq.gz.sam.bam_sorted.bam_PCR_removed.bam",
							"C:/SMD/Luc7/UMI_Removed/ZNF2_3_umi_1_trimmed.fq.gz.sam.bam_sorted.bam_PCR_removed.bam")

#Read in annotations
Introns <- read.delim("C:/SMD/Luc7/annotationFiles/rmMito_ordered_sc3intronsfeb15.saf")
AnnotationFile <- "./annotationFiles/saccharomyces_cerevisiae.20230412.gtf"
ChrAliases <- "./annotationFiles/chrAllias.csv"


AllCounts <-
	Rsubread::featureCounts(
		files = bamfiles,
		isGTFAnnotationFile = TRUE,
		GTF.attrType = "ID",
		GTF.featureType = "gene",
		GTF.attrType.extra = "gene",
		annot.ext = AnnotationFile,
		isPairedEnd = TRUE,
		chrAliases = ChrAliases,
		nthreads = 20,
		nonSplitOnly = FALSE,
		countReadPairs = FALSE,
		allowMultiOverlap = TRUE
	)

saveRDS(AllCounts, file = "./All_Counts.RDS")

#Get 1 NT upstream and downstream of the splice sites
LeftIntrons <- Introns

LeftIntrons$GeneID <- sapply(LeftIntrons$GeneID, paste, "L", sep = ".")

LeftIntrons$IntronName <- Introns$GeneID

LeftIntrons$Start <- Introns$Start - 1

LeftIntrons$End <- Introns$Start


RightIntrons <- Introns

RightIntrons$GeneID <- sapply(RightIntrons$GeneID, paste, "R", sep = ".")

RightIntrons$IntronName <- Introns$GeneID

RightIntrons$Start <- Introns$End

RightIntrons$End <- Introns$End + 1


FlankedIntrons <- rbind(LeftIntrons, RightIntrons)

#Count Unspliced Reads
UnsplicedCounts <- 	Rsubread::featureCounts(
		files = bamfiles,
		annot.ext = FlankedIntrons,
		isPairedEnd = TRUE,
		nonSplitOnly = TRUE,
		nthreads = 20,
		chrAliases = "C:/SMD/Luc7/annotationFiles/chrAllias.csv",
		countReadPairs = FALSE,
		allowMultiOverlap = TRUE,
		minOverlap = 2
	)

saveRDS(UnsplicedCounts, file = "./UnsplicedCounts")

#Get 5 NT outside of either side of intron
LeftIntrons <- Introns

LeftIntrons$GeneID <- sapply(LeftIntrons$GeneID, paste, "L", sep = ".")

LeftIntrons$Start <- Introns$Start - 5

LeftIntrons$End <- Introns$Start - 1


RightIntrons <- Introns

RightIntrons$GeneID <- sapply(RightIntrons$GeneID, paste, "R", sep = ".")

RightIntrons$Start <- Introns$End + 1

RightIntrons$End <- Introns$End + 5


FlankedIntrons <- rbind(LeftIntrons, RightIntrons)

SplicedCounts <- 	Rsubread::featureCounts(
	files = bamfiles,
	annot.ext = FlankedIntrons,
	isPairedEnd = TRUE,
	splitOnly = TRUE,
	nthreads = 20,
	chrAliases = "C:/SMD/Luc7/annotationFiles/chrAllias.csv",
	countReadPairs = FALSE,
	allowMultiOverlap = TRUE,
	minOverlap = 4
)

saveRDS(SplicedCounts, file = "./SplicedCounts")
