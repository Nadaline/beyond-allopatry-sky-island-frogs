<p align="center">
  <img src="assets/logo.png" width="300">
</p>

<h1 align="center">Beyond Allopatry</h1>

<p align="center">
  <b>Isolation by Resistance on Sky-Island Frogs</b>
</p>

<p align="center">
Reproducible analyses accompanying the manuscript.
</p>

---

## Overview

This repository contains all scripts used to reproduce the analyses presented in the manuscript:

> **Beyond Allopatry: Isolation by Resistance on Sky-Island Frogs**

The complete workflow covers raw sequence alignment, SNP calling, population genomics, historical migration inference, and landscape genetics modeling[cite: 1, 2, 3].

---

## Repository Structure

scripts/        # Shell and R scripts for pipeline execution
exploratory/    # Initial data exploration and testing scripts
data/           # Input files (Fastq paths, VCFs, metadata, matrices)
outputs/        # Generated tables, figures, and raw statistics
docs/           # Supplementary files and documentation

---

## Step-by-Step Pipeline & Reproducibility Guide

### Step 1: SNP Calling & Phasing (Bash Protocol)
> **Note:** This script automates alignment using BWA-MEM, cleaning with Picard, and genotyping via GATK v3.8[cite: 3]. Ensure directory variables (contigs_dir, snps_dir) are adjusted to match your environment before execution[cite: 3].
```bash
#!/bin/bash  
# This script is adapted from Ericson et al. 2020 and uses Picard v1.106 & GATK3

contigs_dir=/media/hd1/nadaline/5-correction/consensus 
snps_dir=/media/hd1/nadaline/7-taxon-set/quiriri/SNPs/SNP_calling
samples=("qv33" "qv38" "QR5" "QVA3" "QVA5" "QVA6" "mc48" "mc49" "mc50" "mc51" "mc52" "mc53" "mc54" "mc66" "mc67" "ga56" "ga57" "ga59" "ga60" "ga61" "ga62" "ga63" "ga64" "ga65" "sp68" "sp69" "sp71" "sp72" "sp73" "sp74" "sp77" "sp78" "sp79" "QR10" "QR2" "QR3" "QR4" "QR7" "QVA2" "dzup373" "dzup378" "dzup379" "dzup380" "dzup381" "dzup385" "dzup386" "dzup389" "vt91")

purple="\033[0;45m"
red="\033[0;41m"

# Indexing the outgroup reference
bwa index -p vertebralis -a is vt91.fasta

for i in "${samples[@]}"; do
	contigs="${contigs_dir}/${i}.contigs.fasta"
	if test ! "$snps_dir/$i.MD.bam"; then
		if test ! "$snps_dir/$i.pair.sam"; then
  			echo -e "${red}NÃO FOI encontrado os arquivos sam de $i\033[0m"
  			eval $(echo "bwa mem -t 30 -B 10 -M -R '@RG\tID:$i\tSM:$i\tPL:Illumina' vertebralis /media/hd1/nadaline/1-clean-fastq/${i}/split-adapter-quality-trimmed/${i}-READ1.fastq.gz /media/hd1/nadaline/1-clean-fastq/${i}/split-adapter-quality-trimmed/${i}-READ2.fastq.gz > $i.pair.sam") 
			eval $(echo "bwa mem -t 30 -B 10 -M -R '@RG\tID:$i\tSM:$i\tPL:Illumina' vertebralis /media/hd1/nadaline/1-clean-fastq/$i/split-adapter-quality-trimmed/$i-READ-singleton.fastq.gz > $i.single.sam") 
		else
    		echo -e "${purple}Arquivo sam de $i encontrado, iniciando Picard CleanSam\033[0m"
    		java -jar $snps_dir/picard-tools-1.106/picard.jar CleanSam --INPUT $snps_dir/$i.pair.sam --OUTPUT $snps_dir/$i.pair_CL.sam --VALIDATION_STRINGENCY SILENT
  			java -jar $snps_dir/picard-tools-1.106/picard.jar CleanSam --INPUT $snps_dir/$i.single.sam --OUTPUT $snps_dir/$i.single_CL.sam --VALIDATION_STRINGENCY SILENT
  			samtools sort -o $snps_dir/$i.pair_CL.bam $snps_dir/$i.pair_CL.sam
			samtools sort -o $snps_dir/$i.single_CL.bam $snps_dir/$i.single_CL.sam
			
			# Cleaning intermediate files safely
			if test -f "$snps_dir/$i.pair_CL.bam"; then rm $snps_dir/$i.pair.sam $snps_dir/$i.pair_CL.sam; fi
			if test -f "$snps_dir/$i.single_CL.bam"; then rm $snps_dir/$i.single.sam $snps_dir/$i.single_CL.sam; fi

			if test "$snps_dir/$i.pair_CL.bam"; then
				echo -e "${purple}Arquivo BAM encontrado, iniciando AddOrReplaceReadGroups & MarkDuplicates\033[0m"
				java -Xmx2g -jar $snps_dir/picard-tools-1.106/picard.jar AddOrReplaceReadGroups --INPUT $snps_dir/$i.pair_CL.bam --OUTPUT $snps_dir/$i.pair_RG.bam --SORT_ORDER coordinate --RGPL illumina --RGPU TestXX --RGLB Lib1 --RGID $i --RGSM $i --VALIDATION_STRINGENCY LENIENT
				java -Xmx2g -jar $snps_dir/picard-tools-1.106/picard.jar AddOrReplaceReadGroups --INPUT $snps_dir/$i.single_CL.bam --OUTPUT $snps_dir/$i.single_RG.bam --SORT_ORDER coordinate --RGPL illumina --RGPU TestXX --RGLB Lib1 --RGID $i --RGSM $i --VALIDATION_STRINGENCY LENIENT
				
				java -Xmx2g -jar $snps_dir/picard-tools-1.106/picard.jar MarkDuplicates --INPUT $snps_dir/$i.pair_RG.bam --INPUT $snps_dir/$i.single_RG.bam --OUTPUT $snps_dir/$i.MD.bam --METRICS_FILE $snps_dir/$i.MD.metrics --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 250 --ASSUME_SORTED true --REMOVE_DUPLICATES false
				
				if test -f "$snps_dir/$i.MD.bam"; then rm $snps_dir/$i.single_RG.bam $snps_dir/$i.pair_RG.bam $snps_dir/$i.single_CL.bam $snps_dir/$i.pair_CL.bam; fi
			fi
		fi
	else
		echo -e "${purple}Arquivos MD bam encontrados, iniciando Merge e Variant Calling via GATK3...\033[0m"
		# Merge all bam files into a single matrix (Inputs omitted for readability in overview)
		java -Xmx2g -jar $snps_dir/picard-tools-1.106/picard.jar MergeSamFiles --INPUT ... --OUTPUT $snps_dir/uce-quiriri-merged.bam --SORT_ORDER coordinate
		samtools index $snps_dir/uce-quiriri-merged.bam
		java -Xmx2g -jar $snps_dir/picard-tools-1.106/picard.jar CreateSequenceDictionary --REFERENCE vt91.fasta --OUTPUT $snps_dir/vt91.dict
		samtools faidx vt91.fasta

		# GATK3 Indel Realignment & UnifiedGenotyper
		java -Xmx2g -jar $snps_dir/GenomeAnalysisTK.jar -T RealignerTargetCreator -R vt91.fasta -I $snps_dir/uce-quiriri-merged.bam --minReadsAtLocus 7 -o $snps_dir/uce-quiriri-merged.intervals
		java -Xmx2g -jar $snps_dir/GenomeAnalysisTK.jar -T IndelRealigner -R vt91.fasta -I $snps_dir/uce-quiriri-merged.bam -targetIntervals $snps_dir/uce-quiriri-merged.intervals -LOD 3.0 -o $snps_dir/uce-quiriri-merged_RI.bam
		java -Xmx2g -jar $snps_dir/GenomeAnalysisTK.jar -T UnifiedGenotyper -R vt91.fasta -I $snps_dir/uce-quiriri-merged_RI.bam -gt_mode DISCOVERY -o $snps_dir/uce-quiriri-merged_raw_SNPs.vcf -ploidy 2 -rf BadCigar
		
		# Filtration and Phasing
		java -Xmx2g -jar $snps_dir/GenomeAnalysisTK.jar -T VariantFiltration -R vt91.fasta -V $snps_dir/uce-quiriri-merged_raw_SNPs.vcf --filterExpression "QUAL < 30.0" --filterName "LowQual" -o $snps_dir/uce-quiriri-merged_SNPs_no_indels.vcf -rf BadCigar
		cat $snps_dir/uce-quiriri-merged_SNPs_no_indels.vcf | grep 'PASS\|^#' > $snps_dir/uce-quiriri-merged_SNPs_pass-only.vcf
		java -Xmx2g -jar $snps_dir/GenomeAnalysisTK.jar -T ReadBackedPhasing -R vt91.fasta -I $snps_dir/uce-quiriri-merged_RI.bam --variant $snps_dir/uce-quiriri-merged_SNPs_pass-only.vcf -L $snps_dir/uce-quiriri-merged_SNPs_pass-only.vcf -o $snps_dir/uce-quiriri-merged_SNPs_phased.vcf --phaseQualityThresh 20.0 -rf BadCigar
		break
	fi
done
```

### Step 2: Genetic Structure & TreeMix Export (R Workflow)Note: 
This section contains the complete pipeline to load VCF data, apply quality filtering, run a Principal Coordinate Analysis (PCoA), estimate spatial ancestry proportions using sNMF, and export an LD-pruned dataset formatted for TreeMix.  

```R
# Libraries
library(dartR.popgen)
require(dartR.spatial)
require(dartR.base)
require(vcfR)
require(tess3r)
require(LEA)

# --- 1. Loading and Compliance Checks ---
setwd("/home/nadaline/Documents/Sapinhos/SNP_calling/")
dnv <- read.vcfR("Quiriri-unlike.drop.vcf")
dnv2 <- vcfR2genlight(dnv)
ploidy(dnv2) <- 2
dnv2 <- gl.compliance.check(dnv2)
dnv2 <- gl.add.indmetrics(dnv2, "../Brachy_Metadata.csv")

# --- 2. Quality Control & Filtering ---
dnv3 <- gl.filter.callrate(dnv2, method="loc", threshold=0.7, verbose=0)
dnv4 <- gl.filter.callrate(dnv3, method="ind", threshold=0.7, verbose=0)
dnv5     <- gl.filter.monomorphs(dnv4)

# Optional visualization of missing data profiles
gl.smearplot(dnv5, plot.colors=c("whitesmoke", "mediumorchid", "goldenrod1", "darkseagreen"))

# --- 3. Principal Coordinate Analysis (PCoA) ---
pcoa.dnv <- gl.pcoa(dnv5) 
my.colors <- pals::tol(n=6)
gl.pcoa.plot(x = dnv5, pcoa.dnv, pt.colors = palette(my.colors), ellipse = TRUE)

# --- 4. Sparse Non-Negative Matrix Factorization (sNMF) ---
dnv6 <- dnv5
class(dnv6) <- "genlight"

# Convert genlight object back to VCF, then write to LEA .geno format
gl2vcf(dnv6, plink.bin.path="../../Phd/RADseq/Ameiva/plink_linux_x86_64_20241022/", outpath=getwd(), chr.format="character")
dnv5.vcf <- read.vcfR(file="gl_vcf.vcf")
vcf2geno("./gl_vcf.vcf", "./dnv5.geno")

# Run sNMF algorithms across K=1 to K=6 clusters
obj.snmf = snmf("./dnv5.geno", K=1:6, project="new", repetitions=10, tolerance=0.00001, entropy=TRUE, ploidy=2)
saveRDS(obj.snmf, "brachy-dnv-snmf.rds")

# Evaluate cross-entropy criteria diagnostic plot
plot(obj.snmf, cex=1.2, col="lightblue", pch=19)

# Plotting the Ancestry Matrix Barplot for the selected K (e.g., K=3)
run <- which.min(cross-entropy(obj.snmf, K=3))
my.colors <- c("mediumorchid4", "violetred4", "olivedrab4", "gold3")
barchart(obj.snmf, K=3, run=run, border=NA, space=0, col=my.colors, xlab="Individuals", ylab="Ancestry proportions", main="Ancestry matrix")

# --- 5. Linkage Disequilibrium Pruning & TreeMix Export ---
# Reloading linked variant datasets that include the outgroup (VT91)
dnv_out <- read.vcfR("/media/nadaline/NADALINE II/Backup Linux Doc II/Sapinhos/SNP_calling/Quiriri-linked-out.vcf")
dnv2_out <- vcfR2genlight(dnv_out)
ploidy(dnv2_out) <- 2
dnv2_out <- gl.compliance.check(dnv2_out)
dnv2_out <- gl.add.indmetrics(dnv2_out, "/media/nadaline/NADALINE II/Backup Linux Doc II/Sapinhos/Brachy_Metadata.csv")

# Filter raw data prior to LD report calculations
dnv3_out <- gl.filter.callrate(dnv2_out, method="loc", threshold=0.7, verbose=0)
dnv4_out <- gl.filter.callrate(dnv3_out, method="ind", threshold=0.7, verbose=0)
dnv5_out <- gl.filter.monomorphs(dnv4_out)

# Compute Linkage Disequilibrium maps and filter heavily linked markers
ld_res  <- gl.report.ld.map(dnv5_out, ind.limit=4)
ld_filt <- gl.filter.ld(dnv5_out, ld_res, threshold=0.2)

# Export the final unlinked dataset into a compressed TreeMix archive
gl2treemix(ld_filt, outfile="geo_out_treemix.gz", outpath=getwd())
```
### Step 3: Generalized Dissimilarity Modeling (GDM)
Note: Because resistance and geographic matrices are implemented as raw predictors, an arbitrary placeholder predictor (P) is appended to formatsitepair to properly structure the matrix arrays before running the model.
```R
require(ade4)
require(vegan)
require(gdm)
require(terra)

setwd("/home/nadaline/Documents/Sapinhos/Matriz/")

# Load Pairwise Distance and Environmental Matrices
gen <- read.csv("LandGenReport-pairwise_Gst.Nei.csv", header=T, row.names=1)
top <- read.table("relevo-IBD_resistances.csv", header=T)
cli <- read.table("ann-IBR_resistances.csv", header=T)
gra <- read.table("barriergrasslands_resistances.csv", header=T)
pop <- read.table("/home/nadaline/Documents/Sapinhos/Pop.csv", header=T, row.names=c("QV","GA","MC","SP","QR","TR"))[,2:3]

# Enforce strict population ordering
site <- c("QV","GA","MC","SP","QR","TR")
gen  <- gen[site, site]

# Wipe row/column text attributes for alignment consistency
gen <- gen |> `rownames<-`(NULL) |> `colnames<-`(NULL)
gra <- gra |> `rownames<-`(NULL) |> `colnames<-`(NULL)
cli <- cli |> `rownames<-`(NULL) |> `colnames<-`(NULL)
top <- top |> `rownames<-`(NULL) |> `colnames<-`(NULL)

envdata <- list(as.matrix(cbind(site, top)), as.matrix(cbind(site, gra)), as.matrix(cbind(site, cli)))

# Constructing GDM input data frames
pred    <- cbind(site, pop, as.numeric(c(1,1,1,1,1,1)))
names(pred) <- c("site", "X", "Y", "P")
biores  <- cbind(site, gen)

pairtable := formatsitepair(bioData=biores, bioFormat=3, siteColumn="site", predData=pred, XColumn="X", YColumn="Y", verbose=T)
gdm.tab   := formatsitepair(pairtable, bioFormat=4, siteColumn="site", predData=pred, verbose=T, distPreds=envdata)

# Run GDM Analysis & Variable Importance
res1 <- gdm(gdm.tab, geo=T)
summary(res1)
plot(res1)

# Run Permutation Tests (n=100) for Predictor Significance
modTest_b <- gdm.varImp(gdm.tab, geo=T, nPerm=100)
uu <- data.frame(row.names(modTest_b$`Predictor Importance`), modTest_b$`Predictor Importance`)
us <- uu[rev(order(uu$All.predictors)), ]
barplot(us$All.predictors, names.arg=us[,1], cex.names=0.75, main="Predictor Importance")
```

This repository is still under active development.
