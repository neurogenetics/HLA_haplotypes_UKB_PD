# HLA haplotypes UKBiobank in PD
Investigation of HLA haplotypes in the UK Biobank vs Parkinson's disease

May 12th 2020

Cornelis, Mary, Andy, Mike

#### TL:DR
One of the Parkinson Disease GWAS hits rs504594 (AKA rs112485576) is tagging the DQA1_301 HLA haplotype and this haplotype is confirmed to be the most PD associated HLA haplotype using an improved imputation method specific to HLA haplotypes in the UK Biobank cohort including only European individuals (mainly British). 

#### Longer summary
```
Summary:
In the most recent PD GWAS (PMID: 31701892) a couple hits are in close proximity of the HLA locus.
UKbiobank (https://www.ukbiobank.ac.uk/) is largest semi-public genetic resource out there and 
has imputation HLA haplotypes available for all genotypes participants (~500K).
Here we check:
1) Association between HLA haplotypes and PD (and PD-proxy status => having a parent with PD) 
2) Correlation between HLA haplotypes and the GWAS hits in the HLA locus region

Results:
1) Meta-analyzing association results identified DQA1_301 passing Bonferroni correction 
P=0.00015, OR=0.90, SE=0.0264, this includes => 1529 PD cases, 13404 Proxy PD cases and 347145 controls.
2) HLA haplotype DQA1_301 is highly correclated with rs504594 (number of A alleles) with >0.95 correlation 

Whats next:
Next key is to understand what is special about the DQA1_301 HLA haplotype and how this is somehow protective for PD.

```

## Data and description download
### Get actual data 
Download new data package from UKB

```
module load ukbb/0.1
ukbunpack ukb41967.enc SECRET_KEY
ukbconv ukb41967.enc_ukb csv
```
### Get header of UKB HLA file:
Resource 2182

Name:HLA data headers

Headers associated with imputed HLA values in release version 2.

```
wget -nd biobank.ndph.ox.ac.uk/showcase/showcase/auxdata/ukb_hla_v2.txt
```

### Resource 1520
Name: Imputed HLA values example

Example of imputed HLA loci information for a single participant.

```
wget -nd biobank.ndph.ox.ac.uk/showcase/showcase/examples/eg_hla_impute.dat
```

### PDF description
Resource 182

Name:Imputation of classical HLA types

https://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=182


## Make and prepare analysis files

### UKbiobank HLA data prep

```
# replace comma with tab
sed 's/,/\t/g' ukb41967.csv > ukb41967_tab.txt
sed -i 's/"//g' ukb41967_tab.txt
# replace header, note added SAMPLEID to UKB header
sed '1d' ukb41967_tab.txt > tmpfile; cat ukb_hla_v2.txt tmpfile > ukb41967_tab.txt
rm tmpfile
# ukb41967_tab.txt is the file for later use...
```

### Get phenotypes and covariates

### PD
```
PD.txt => created from UKB field => 42033
```

### PD parent
```
PD_parent_no_PD.txt => created from UKB field => 20107 and 20110 but excluding PD cases
```

### Covariates => PC, age, townsend

#### PC's
```
# creating new UKB PC's 
module load plink
module load R
module load flashpca
module load GCTA

# UKBB_raw_data_no_cousins created using GCTA with 0.125 cut-off
plink --bfile UKBB_raw_data --maf 0.05 --geno 0.05 --hwe 1E-6 --make-bed --out filter --threads 19 --memory 500000
plink --bfile filter --indep-pairwise 50 5 0.5 --out prune --threads 19 --memory 500000
plink --bfile filter --extract prune.prune.in --make-bed --out prune --threads 19 --memory 500000
gcta64 --bfile prune --make-grm --out GRM_matrix --autosome --maf 0.05 --thread-num 20 
gcta64 --grm-cutoff 0.125 --grm GRM_matrix --out GRM_matrix_0125 --make-grm --thread-num 20 
gcta64 --bfile UKBB_raw_data --keep GRM_matrix_0125.grm.id --make-bed --out UKBB_raw_data_no_cousins --thread-num 20

# then make PC's
plink --bfile UKBB_raw_data_no_cousins --maf 0.01 --geno 0.01 --hwe 5e-6 --make-bed --out filter_better \
--exclude range exclusion_regions_hg19.txt --memory 230000 --threads 19
plink --bfile filter_better --indep-pairwise 1000 10 0.02 --out pruned_better_and_more --memory 230000 --threads 19
plink --bfile filter_better --extract pruned_better_and_more.prune.in --make-bed --out input_pca --memory 230000 --threads 19

# 29953 SNPs and 362788 individuals
flashpca --bfile input_pca --suffix _testing_slurm.txt --numthreads 19
mv pcs_testing_slurm.txt pc.txt
```

#### Covariates
```
covariates.txt bases on UKB fields:
# genetic-sex = 22001
# european = 22006
# PC's = 22009
# batch = 22000
# townsend = 189
# Year of birth = 34
# age of recruitment = 21022
```

## Load in R and perform analysis
### Purpose => loading in data to get it all in the right space
```
module load R
R
require(data.table)
library(dplyr)
HLA <- fread("ukb41967_tab.txt",fill=TRUE,header=T)
# set low numbers to missing per recommendations of PDF description HLA UK Biobank
HLA_v2 <- HLA %>% mutate_all(funs(ifelse(.>0 & .<0.7, NA, .)))
PD_case <- fread("PD.txt",header=T)
PD_case$eid <- NULL
proxy_case <- fread("PD_parent_no_PD.txt",header=F)
proxy_case$V2 <- NULL
# COV =>    FID   IID BIRTH_YEAR TOWNSEND AGE_OF_RECRUIT BATCH GENETIC_SEX
COV <- fread("covariates.txt",header=T)
# PC's
PC <- fread("pc.txt",header=T)
PC$IID <- NULL
# data sizes:
dim(HLA_v2)
# [1] 502505    363 
# 502505 samples and 362 HLA haplogroups 
dim(PD_case)
# [1] 2087 PD cases
dim(proxy_case)
# [1] 17987 Proxy cases
dim(PC)
# [1] 362788     11
# 362788 samples and 10 PC's
dim(COV)
# [1] 502616      8
# 502616 samples and known covariates
```

## Process data
### Purpose => process data to get it all in the right and easy to use format
```
# merge and process (!!)
data <- merge(HLA_v2, PC, by.x = "SAMPLEID", by.y = "FID")
datav2 <- merge(data, COV, by.x = "SAMPLEID", by.y = "FID")

# phenotypes
PD <- merge(datav2, PD_case, by.x = "SAMPLEID", by.y = "eid")
mydata2 = select(PD, -381, -382)
PD <- mydata2
proxy <- merge(datav2, proxy_case, by.x = "SAMPLEID", by.y = "V1")
# dim(PD)
[1] 1529 380
# dim(proxy)
[1] 13404  380

# controls 
CONTROL <- anti_join(datav2,PD,by="SAMPLEID")
CONTROL_v2 <- anti_join(CONTROL,proxy,by="SAMPLEID")
dim(CONTROL_v2)
# [1] 347145  380

# 347145 = all need 1/10 for PD and 9/10 proxy
# =34715 for PD and 312430 for proxy

# random order and then subset
random_order <- CONTROL_v2[sample(1:nrow(CONTROL_v2)), ]
controls_PD <- random_order[1:34715,]
controls_proxy <- anti_join(random_order,controls_PD,by="SAMPLEID")
controls_PD$PHENO <- 0
controls_proxy$PHENO <- 0
PD$PHENO <- 1
proxy$PHENO <- 1

# merge with cases
PD_test <- rbind(PD,controls_PD)
Proxy_test <- rbind(proxy,controls_proxy)
```

## Do initial HLA association analysis
### Purpose => perform case/proxy/control regression of HLA haplotypes 

```
# for all 362 HLA haplotypes from UKB, once for PD and once for PD_proxy:

thisFormula1 <- formula(paste("PHENO ~ A_101 + TOWNSEND + AGE_OF_RECRUIT + GENETIC_SEX + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
model1_PD <- glm(thisFormula1, data = PD_test, family = "binomial")
model1_Proxy <- glm(thisFormula1, data = Proxy_test, family = "binomial")
summary(model1_PD)
summary(model1_Proxy)

```

## Meta-analyze HLA association analysis
### Purpose => meta-analyze case/control and proxy/control results

#### Rescale proxy association results
```
module load python/3.7
# python proxy_gwas_gwaxStyle.py --infile short_format.csv --beta-proxy BETA --se-proxy SE --p-proxy P --outfile remove.csv

python proxy_gwas_gwaxStyle.py --infile proxy_HLA_results.csv --beta-proxy Estimate --se-proxy SE --p-proxy P --outfile proxy_HLA_results_rescaled.csv

```

#### Meta analyze results in metal

```
module load metal

# run like this:
metal metal_HLA.txt

```

#### metal file looks like this:
```
#../generic-metal/metal metalAll.txt
#THIS SCRIPT EXECUTES AN ANALYSIS OF EIGHT STUDIES
#THE RESULTS FOR EACH STUDY ARE STORED IN FILES Inputfile1.txt THROUGH Inputfile8.txt
SCHEME  STDERR
AVERAGEFREQ ON
MINMAXFREQ ON
LABEL TotalSampleSize as N # If input files have a column for the sample size labeled as 'N'
# LOAD THE FIRST SEVEN INPUT FILES

# UNCOMMENT THE NEXT LINE TO ENABLE GenomicControl CORRECTION
# GENOMICCONTROL ON

# === DESCRIBE AND PROCESS THE FIRST INPUT FILE ===
MARKER HAPLO
# ALLELE minorAllele majorAllele
FREQ   FREQ_CONTROL
EFFECT Estimate
STDERR SE
PVALUE P
WEIGHT N 
PROCESS toMeta.UKB_HLA_PD_control.tab

# === DESCRIBE AND PROCESS THE SECOND INPUT FILE ===
MARKER HAPLO
# ALLELE minorAllele majorAllele
FREQ   FREQ_CONTROL
EFFECT b_adjusted
STDERR se_adjusted
PVALUE p_derived
WEIGHT N
PROCESS toMeta.UKB_HLA_Proxy_control.tab

OUTFILE HLA_association_UKB_PD .tbl
ANALYZE HETEROGENEITY

QUIT

```

#### Results and forest plot of most significant result

```
# From metal:
## Smallest p-value is 0.0001419 at marker 'DQA1_301'
# Realistic multiple test correction:
# 103 of HLA haplotypes tested in both PD and Proxy + with MAF of >1%
# Bonferroni => 0.05 / 102 = 0.00049
# meaning that only DQA1_301 passes multiple testing correction
```

```
# Forest plot of results
R
require("rmeta")
library(metafor)
data <- read.table("input_forest_HLA.txt", header = T)
labs <- data$DATASET
yi   <- data$BETA
sei  <- data$SE
resFe  <- rma(yi=yi, sei=sei, method="FE")
resRe  <- rma(yi=yi, sei=sei)
print(summary(resFe))
print(summary(resRe))
pdf(file = "HLA_forest_UKB.pdf", width = 8, height = 6)
forest(resFe, xlim=c(resFe$beta-0.9931472,resFe$beta+0.6931472), main="Meta-analysis of DQA1_301 HLA",atransf=exp, xlab=paste("Odds Ratio (95%CI) for SNP",sep=""), slab=labs, mlab="Fixed Effects", col = "red", border = "red", cex=.9)
dev.off()
png(file = "HLA_forest_UKB.png")
forest(resFe, xlim=c(resFe$beta-0.9931472,resFe$beta+0.6931472), main="Meta-analysis of DQA1_301 HLA",atransf=exp, xlab=paste("Odds Ratio (95%CI) for SNP",sep=""), slab=labs, mlab="Fixed Effects", col = "red", border = "red", cex=.9)
dev.off()

```

#### Forest plot of most significant result

<p align="center">
  <img width="600" height="600" src="HLA_forest_UKB.png">
</p>


## Get HLA frequencies
### Purpose => saving the dataframes to get case/proxy/control frequencies of haplotypes

```
# saving dataframe to get frequencies
# note using rounded calls instead of dosages...
# PD
newdata <- PD_test[,c(2:363)]
newdata2 <- round(newdata)
phenos <- PD_test[,c(381)]
IDS <- PD_test[,c(1)]
PD_SAVE <- cbind(IDS,phenos,newdata2)
write.table(PD_SAVE, "PD_testing_data.txt", quote = F, sep = "\t", row.names = F)
# Proxy
newdata <- Proxy_test[,c(2:363)]
newdata2 <- round(newdata)
phenos <- Proxy_test[,c(381)]
IDS <- Proxy_test[,c(1)]
PROXY_SAVE <- cbind(IDS,phenos,newdata2)
write.table(PROXY_SAVE, "Proxy_testing_data.txt", quote = F, sep = "\t", row.names = F)
# frequency calculations:
# mean/2 per group
```

## Check correlation between HLA haplotypes and GWAS variants
### Purpose => just checking if any of the known GWAS variants are tagging the HLA haplotypes of UK Biobank

```
# check correlation / R2 with GWAS hits:
# GWAS hits from https://pdgenetics.shinyapps.io/GWASBrowser/

# 29	rs4140646	6	27738801	LOC100131289
# 30	rs9261484	6	30108683	TRIM40
# 31	rs112485576	6	32578772	HLA-DRB5 # NOTE IS rs504594

# check presence in imputed UKB data:
grep rs4140646 chr6.UKBB.EU.filtered.pvar # <= yes
grep rs9261484 chr6.UKBB.EU.filtered.pvar # <= yes
grep rs504594 chr6.UKBB.EU.filtered.pvar # <= yes

# extract variants from imputed data
./plink2 --pfile chr6.UKBB.EU.filtered --snps rs4140646,rs9261484,rs504594 \
--make-bed --out /data/CARD/UKBIOBANK/HLA/HLA_PD_region_variants

# recode variants for loading into R
plink --bfile HLA_PD_region_variants --recodeA --out HLA_PD_region_variants

# merge  with genotypes HLA
module load R
R
require(data.table)
HLA <- fread("PD_testing_data.txt",header=T)
geno <- fread("HLA_PD_region_variants.raw",header=T)
MERGE <- merge(geno, HLA, by.x="FID",  by.y="SAMPLEID")
write.table(MERGE, "PD_testing_data_WITH_VARIANTS.txt", quote = F, sep = "\t", row.names = F)

# check correlation with all 362 HLA haplotypes vs 3 variants...
# correlation table => correlation_HLA_GWAS_variants.txt

# also note => very nice tool https://biobankengine.shinyapps.io/hla-map/

# results:
# perfect correlation between:
# data:  MERGE$DQA1_301 and MERGE$rs504594_A
# 0.9506786 (Pearson)
# so number of A alleles correlates with DQA1_301 alleles
# NOTE that DQA1_301 is protective
```

### All done

![myImage](https://media.giphy.com/media/XRB1uf2F9bGOA/giphy.gif)
