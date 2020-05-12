# HLA_haplotypes_UKB_PD
Investigation of HLA haplotypes in the UK Biobank vs Parkinson's disease


## Data and description download
### get actual data 
Download new data package from UKB
module load ukbb/0.1
ukbunpack ukb41967.enc SECRET_KEY
ukbconv ukb41967.enc_ukb csv

### get header of UKB HLA file:
Resource 2182
Name:HLA data headers
Headers associated with imputed HLA values in release version 2.
wget -nd biobank.ndph.ox.ac.uk/showcase/showcase/auxdata/ukb_hla_v2.txt

### Resource 1520
Name:Imputed HLA values example
Example of imputed HLA loci information for a single participant.
wget -nd biobank.ndph.ox.ac.uk/showcase/showcase/examples/eg_hla_impute.dat

### PDF description
Resource 182
Name:Imputation of classical HLA types

