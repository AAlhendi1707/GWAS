# Genome-wide association analysis
Now that our data is loaded, filtered, and additional SNP genotypes imputed we are ready to perform genome-wide association analysis. This involves regressing each SNP separately on a given trait, adjusted for sample level clinical, environmental, and demographic factors. Due to the large number of SNPs and the generally uncharacterized relationships to the outcome, a simple single additive model will be employed.

The `GWAA` function requires two arguments. The `genodata` argument should specify the entire genotype data object in `SnpMatrix` format. The phenodata argument should be a data frame with a column of sample IDs, corresponding to the row names of genodata, and a columns for the continuous outcome variable. These columns must be named “id” and “phenotype”, respectively. In order to fit the model, genotype data is converted to numeric format using the as function from snpStats. The genotypes of each SNP are then coded as continuous, thereby taking on the value of 0, 1, and 2. For this example, we wish for the value of the genotype to reflect the number of minor alleles. However, following conversion our values will reflect the opposite. To fix this a flip.matrix procedure is included in our `GWAA` function, which can be turned on or off using the `flip` argument.

Due to the large number of models that require fitting, the GWA analysis can be deployed in parallel across multiple processors or machines to reduce the running time. Here we demonstrate two basic methods for performing parallel processing using the doParallel package. This will be carried out differently depending on whether or not the analysis is run on a UNIX based system, though the arguments are the same. The user can specify the number of processes using the `worker` argument (set to 2 by default). Additional arguments include `select.snps` and `nSplits`. The former allows the user to subset the analysis via a vector of SNP IDs. The latter specifies a number of SNP-wise splits that are made to the genotype data. The function runs the GWA analysis on these smaller subsets of the genotype data one at a time. After each subset has finished running the function will print a progress update onto the R console. By default this is set to 10.

## Association analysis of typed SNPs - Step 7
First we create a data frame of phenotype features that is the concatenation of clinical features and the first ten principal components. The HDL feature is normalized using a rank-based inverse normal transform. We then remove variables that we are not including in the analysis, i.e. HDL(non-normalized), LDL, TG, and CAD. Finally, we remove samples with missing normalized HDL data.

```r
## restore the data generated from steps 1-6
load (""Genotype.SNVsfiltered.imputed.target.rultes.Rdata")

## Require GenABEL and GWAA function
library(GenABEL)
source("https://github.com/AAlhendi1707/GWAS/blob/master/R/GWAA.R?raw=true")

# Merge clincal data and principal components to create phenotype table
phenoSub <- merge(clinical,pcs)      # data.frame => [ FamID CAD sex age hdl pc1 pc2 ... pc10 ]

# We will do a rank-based inverse normal transformation of hdl
phenoSub$phenotype <- rntransform(phenoSub$hdl, family="gaussian")

# Show that the assumptions of normality met after transformation
par(mfrow=c(1,2))
hist(phenoSub$hdl, main="Histogram of HDL", xlab="HDL")
hist(phenoSub$phenotype, main="Histogram of Tranformed HDL", xlab="Transformed HDL")
```
<img src="imgs/HDL-transformation.png">

```r
# Remove unnecessary columns from table
phenoSub$hdl <- NULL
phenoSub$ldl <- NULL
phenoSub$tg <- NULL
phenoSub$CAD <- NULL

# Rename columns to match names necessary for GWAS() function
phenoSub <- rename(phenoSub, replace=c(FamID="id"))

# Include only subjects with hdl data
phenoSub<-phenoSub[!is.na(phenoSub$phenotype),]
# 1309 subjects included with phenotype data

print(head(phenoSub))
```

