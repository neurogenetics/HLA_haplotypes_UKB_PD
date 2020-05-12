# HLA_haplotypes_UKB_PD
Investigation of HLA haplotypes in the UK Biobank vs Parkinson's disease


## Data and description download
### get actual data 
Download new data package from UKB

```
module load ukbb/0.1
ukbunpack ukb41967.enc SECRET_KEY
ukbconv ukb41967.enc_ukb csv
```
### get header of UKB HLA file:
Resource 2182

Name:HLA data headers

Headers associated with imputed HLA values in release version 2.

```
wget -nd biobank.ndph.ox.ac.uk/showcase/showcase/auxdata/ukb_hla_v2.txt
```

### Resource 1520
Name:Imputed HLA values example

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
```

### get phenotypes and covariates

### PD

### PD parent


### Covariates => PC, age, townsend



## Load in R and perform analysis

### Load in data
```
module load R
R
require(data.table)
library(dplyr)
HLA <- fread("ukb41967_tab.txt",fill=TRUE,header=T)
PD_case <- fread("PD.txt",header=T)
PD_case$eid <- NULL
proxy_case <- fread("PD_parent_no_PD.txt",header=F)
proxy_case$V2 <- NULL
# COV =>    FID   IID BIRTH_YEAR TOWNSEND AGE_OF_RECRUIT BATCH GENETIC_SEX
COV <- fread("covariates.txt",header=T)
samples_to_keep <- fread("PD_grs.txt",header=T)
samples_to_keep$IID <- NULL
samples_to_keep$PHENO <- NULL
samples_to_keep$CNT <- NULL
samples_to_keep$CNT2 <- NULL
samples_to_keep$SCORE <- NULL
# PC's
PC <- fread("pc.txt",header=T)
PC$IID <- NULL
```

### process data
### Purpose => process data to get it all in the right and easy to use format
```
# merge and process (!!)
data <- merge(HLA, samples_to_keep, by.x = "SAMPLEID", by.y = "FID")
datav2 <- merge(data, COV, by.x = "SAMPLEID", by.y = "FID")
datav3 <- merge(datav2, PC, by.x = "SAMPLEID", by.y = "FID")

# phenotypes
PD <- merge(datav3, PD_case, by.x = "SAMPLEID", by.y = "eid")
mydata2 = select(PD, -381, -382)
PD <- mydata2
proxy <- merge(datav3, proxy_case, by.x = "SAMPLEID", by.y = "V1")
# dim(PD)
[1] 1529 380
# dim(proxy)
[1] 13404  380

# controls 
CONTROL <- anti_join(datav3,PD,by="SAMPLEID")
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

## Do analysis
### Purpose => perform case/proxy/control regression of HLA haplotypes 

```
# for all 362 HLA haplotypes from UKB, once for PD and once for PD_proxy:

thisFormula1 <- formula(paste("PHENO ~ A_101 + TOWNSEND + AGE_OF_RECRUIT + GENETIC_SEX + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"))
model1_PD <- glm(thisFormula1, data = PD_test, family = "binomial")
model1_Proxy <- glm(thisFormula1, data = Proxy_test, family = "binomial")
summary(model1_PD)
summary(model1_Proxy)

```

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

### all done
