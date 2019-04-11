# Data generation

## Re-compute PCA for modeling - Step 5

Now that we have performed SNP and sample level quality control on our genotype data, we will calculate principal components to be included as covariates in the GWA models. These serve to adjust for any remaining substructure that may confound SNP level association. As with Ancestry filtering we will calculate PCs using the `snpgdsPCA` function from SNPRelate, after performing LD pruning once again on the filtered genotype data set. In this example, we will include the first 10 principal components in GWA models.

```r
## load R data from steps 1-4
load("save.image("GWAS.Steps1-4.Rdata")
```
```r
#Set LD threshold to 0.2
ld.thresh <- 0.2

set.seed(1000)
geno.sample.ids <- rownames(genotype)
snpSUB <- snpgdsLDpruning(genofile, ld.threshold = ld.thresh,
                          sample.id = geno.sample.ids, # Only analyze the filtered samples
                          snp.id = colnames(genotype)) # Only analyze the filtered SNPs
```
```
## SNP pruning based on LD:
## Excluding 204583 SNPs on non-autosomes
## Excluding 0 SNP (monomorphic: TRUE, < MAF: NaN, or > missing rate: NaN)
## Working space: 1401 samples, 656890 SNPs
##  Using 1 (CPU) core
##  Sliding window: 500000 basepairs, Inf SNPs
##  |LD| threshold: 0.2
## Chromosome 1: 8.23%, 5845/71038
## Chromosome 3: 8.08%, 4893/60565
## Chromosome 6: 8.03%, 4352/54176
## Chromosome 12: 8.56%, 3606/42124
## Chromosome 21: 9.41%, 1173/12463
## Chromosome 2: 7.66%, 5647/73717
## Chromosome 4: 8.20%, 4567/55675
## Chromosome 7: 8.49%, 3939/46391
## Chromosome 11: 7.89%, 3489/44213
## Chromosome 10: 7.96%, 3814/47930
## Chromosome 8: 7.65%, 3694/48299
## Chromosome 5: 8.04%, 4514/56178
## Chromosome 14: 8.77%, 2460/28054
## Chromosome 9: 8.21%, 3374/41110
## Chromosome 17: 11.14%, 2222/19939
## Chromosome 13: 8.30%, 2843/34262
## Chromosome 20: 9.39%, 2137/22753
## Chromosome 15: 9.23%, 2390/25900
## Chromosome 16: 9.27%, 2558/27591
## Chromosome 18: 8.87%, 2327/26231
## Chromosome 19: 12.99%, 1491/11482
## Chromosome 22: 10.92%, 1243/11382
## 72578 SNPs are selected in total.
```
```r
snpset.pca <- unlist(snpSUB, use.names=FALSE)
cat(length(snpset.pca),"\n")  #72578 SNPs will be used in PCA analysis

pca <- snpgdsPCA(genofile, sample.id = geno.sample.ids,  snp.id = snpset.pca, num.thread=1)
```
```
## Principal Component Analysis (PCA) on SNP genotypes:
## Excluding 788895 SNPs on non-autosomes
## Excluding 0 SNP (monomorphic: TRUE, < MAF: NaN, or > missing rate: NaN)
## Working space: 1401 samples, 72578 SNPs
##  Using 1 (CPU) core
## PCA: the sum of all working genotypes (0, 1 and 2) = 32714193
## PCA: Wed Jun 24 17:02:49 2015    0%
## PCA: Wed Jun 24 17:34:40 2015    100%
## PCA: Wed Jun 24 17:34:40 2015    Begin (eigenvalues and eigenvectors)
## PCA: Wed Jun 24 17:34:41 2015    End (eigenvalues and eigenvectors)
```
```r
# Find and record first 10 principal components
# pcs will be a N:10 matrix.  Each column is a principal component.
pcs <- data.frame(FamID = pca$sample.id, pca$eigenvect[,1 : 10],
                  stringsAsFactors = FALSE)
colnames(pcs)[2:11]<-paste("pc", 1:10, sep = "")

print(head(pcs))
```
```
##   FamID          pc1          pc2           pc3           pc4
## 1 10002  0.007764870  0.014480384 -0.0006315881  0.0028664643
## 2 10004 -0.012045108 -0.007231015 -0.0030012896 -0.0107972693
## 3 10005 -0.016702930 -0.005347697  0.0144498361 -0.0006151058
## 4 10007 -0.009537235  0.004556977  0.0026835662  0.0166255657
## 5 10008 -0.015392106 -0.002446933  0.0205087909 -0.0057241772
## 6 10009 -0.015123858 -0.002353917  0.0213604518  0.0069156529
##             pc5          pc6           pc7          pc8          pc9
## 1 -0.0188391406  0.009680646  0.0276468057 -0.006645818 -0.023429747
## 2 -0.0077705400 -0.004645751  0.0018061075 -0.003087891 -0.001833242
## 3  0.0345170160  0.038708551  0.0205790788 -0.012265508  0.003592690
## 4 -0.0002363142  0.005514627  0.0159588869  0.027975455  0.029777180
## 5 -0.0039696226  0.005354244 -0.0007269312  0.027014714  0.010672162
## 6  0.0400677558  0.023222478  0.0152485234  0.013296852  0.022746352
##           pc10
## 1  0.010492314
## 2 -0.004538746
## 3 -0.002287043
## 4 -0.007461255
## 5 -0.003352997
## 6  0.013143889
```
```r
# Close GDS file
closefn.gds(genofile)
```
## Genotype imputation - Step 6
In addition to the genotyped SNPs from our study, it is useful to extend the analysis to other known SNPs, that were not typed or were removed by SNP level filtering. In this example, we impute SNPs on chromosome 16.

Performance of genotype imputation requires reference data, which has typed genotypes at the SNPs of interest from similar homogeneous sample. Sources for this data include HapMap and 1000 Genomes.

For this example, we will use 1000 Genomes data, read in from .ped and.info using the `read.pedfile` in from snpStats. Note, that the .info file is similar to the .map file. To specify the column in the .info file with the SNP IDs, we use the `which` argument.

We derive imputation “rules” for the additional SNPs that were not typed in our study using `snp.imputation` based on the genotypes from the 1000 Genomes data. Each rule represents a predictive model for genotypes of untyped SNPs associated with near-by typed SNPs. Using these rules, we can calculate the expected posterior value of the non-typed SNPs using the `impute` function from SNPRelate.

In the last step we remove un-typed SNPs in which we fail to derive imputation “rules”. We also filter out SNPs that have low estimated minor allele frequency, and low imputation accuracy. The latter is based on the R2 value of the model estimated by the `snp.imputation` function.

```r
# Read in 1000g data for given chromosome 16
thougeno <- read.pedfile(onethou.fn$ped, snps = onethou.fn$info, which=1)

# Obtain genotype data for given chromosome
genoMatrix <- thougeno$genotypes

# Obtain the chromosome position for each SNP
support <- thougeno$map
colnames(support)<-c("SNP", "position", "A1", "A2")
head(support)
```
```
##           SNP position A1 A2
## 1 rs140769322    60180  3  2
## 2 rs188810967    60288  2  1
## 3  rs76368850    60291  2  4
## 4 rs185537431    60778  3  1
## 5 rs542544747    60842  2  1
## 6   rs4021615    61349  1  3
```
```r
# Imputation of non-typed 1000g SNPs
presSnps <- colnames(genotype)

# Subset for SNPs on given chromosome
presSnps <- colnames(genotype)
presDatChr <- genoBim[genoBim$SNP %in% presSnps & genoBim$chr==16, ]
targetSnps <- presDatChr$SNP

# Subset 1000g data for our SNPs
# "missing" and "present" are snpMatrix objects needed for imputation rules
is.present <- colnames(genoMatrix) %in% targetSnps

missing <- genoMatrix[,!is.present]
print(missing)             # Almost 400,000 SNPs
```
```
## A SnpMatrix with  99 rows and  377819 columns
## Row names:  CEU_1 ... CEU_99 
## Col names:  rs140769322 ... rs111706106
```
```r
present <- genoMatrix[,is.present]
print(present)                  # Our typed SNPs
```
```
## A SnpMatrix with  99 rows and  20632 columns
## Row names:  CEU_1 ... CEU_99 
## Col names:  rs41340949 ... rs4785775
```
```r
# Obtain positions of SNPs to be used for imputation rules
pos.pres <- support$position[is.present]
pos.miss <- support$position[!is.present]

# Calculate and store imputation rules using snp.imputation()
rules <- snp.imputation(present, missing, pos.pres, pos.miss)
```
```
## SNPs tagged by a single SNP: 82119
## SNPs tagged by multiple tag haplotypes (saturated model): 115769
```
```r
# Remove failed imputations
rules <- rules[can.impute(rules)]
cat("Imputation rules for", length(rules), "SNPs were estimated\n")  
# Imputation rules for 197888 SNPs were estimated
```
```
## Imputation rules for 197888 SNPs were estimated
```

```r
# Quality control for imputation certainty and MAF
# Set thresholds
r2threshold <- 0.7
minor <- 0.01

# Filter on imputation certainty and MAF
rules <- rules[imputation.r2(rules) >= r2threshold]

cat(length(rules),"imputation rules remain after imputations with low certainty were removed\n")  
# 162565 imputation rules remain after imputations with low certainty were removed
```
```r
rules <- rules[imputation.maf(rules) >= minor]
cat(length(rules),"imputation rules remain after MAF filtering\n")  
# 162565 imputation rules remain after MAF filtering
```
```r
# Obtain posterior expectation of genotypes of imputed snps
target <- genotype[,targetSnps]
imputed <- impute.snps(rules, target, as.numeric=FALSE)
print(imputed)  # 162565 SNPs were imputed
```
```
## A SnpMatrix with  1401 rows and  162565 columns
## Row names:  10002 ... 11596 
## Col names:  rs560777354;rs80001234 ... rs62053708
```
```r
# Free some memory in your R session
rm(genoMatrix)
rm(missing)
rm(present)

# Add new imputed, target and rules data to saved results
save.image("GWAS.Steps1-6.Rdata".Rdata")
```
