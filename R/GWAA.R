# ---- gwaa ----

# Genome-wide Association Analysis
# Parallel implementation of linear model fitting on each SNP

GWAA <- function(genodata=genotypes,  phenodata=phenotypes, family = gaussian, filename=NULL,
                 append=FALSE, workers=getOption("mc.cores",2L), flip=TRUE,
                 select.snps=NULL, hosts=NULL, nSplits=10)
{
    if (!require(doParallel)) { stop("Missing doParallel package") }
    
    #Check that a filename was specified
    if(is.null(filename)) stop("Must specify a filename for output.")
    
    #Check that the genotype data is of class 'SnpMatrix'
    if( class(genodata)!="SnpMatrix") stop("Genotype data must of class 'SnpMatrix'.")
    
    #Check that there is a variable named 'phenotype' in phenodata table
    if( !"phenotype" %in% colnames(phenodata))  stop("Phenotype data must have column named 'phenotype'")
    
    #Check that there is a variable named 'id' in phenodata table
    if( !"id" %in% colnames(phenodata)) stop("Phenotype data must have column named 'id'.")
    
    #If a vector of SNPs is given, subset genotype data for these SNPs
    if(!is.null(select.snps)) genodata<-genodata[,which(colnames(genodata)%in%select.snps)]
    
    #Check that there are still SNPs in 'SnpMatrix' object
    if(ncol(genodata)==0) stop("There are no SNPs in the 'SnpMatrix' object.")
    
    #Print the number of SNPs to be checked
    cat(paste(ncol(genodata), " SNPs included in analysis.\n"))
    
    #If append=FALSE than we will overwrite file with column names
    if(!isTRUE(append)) {
        columns<-c("SNP", "Estimate", "Std.Error", "t-value", "p-value")
        write.table(t(columns), filename, row.names=FALSE, col.names=FALSE, quote=FALSE)
    }

    # Check sample counts
    if (nrow(phenodata) != nrow(genodata)) {
        warning("Number of samples mismatch.  Using subset found in phenodata.")
    }

    # Order genodata rows to be the same as phenodata
    genodata <- genodata[phenodata$id,]

    cat(nrow(genodata), "samples included in analysis.\n")

    # Change which allele is counted (major or minor)
    flip.matrix<-function(x) {
        zero2 <- which(x==0)
        two0 <- which(x==2)
        x[zero2] <- 2
        x[two0] <- 0
        return(x)
    }

    nSNPs <- ncol(genodata)
    genosplit <- ceiling(nSNPs/nSplits) # number of SNPs in each subset

    snp.start <- seq(1, nSNPs, genosplit) # index of first SNP in group
    snp.stop <- pmin(snp.start+genosplit-1, nSNPs) # index of last SNP in group

    if (is.null(hosts)) {
        # On Unix this will use fork and mclapply.  On Windows it
        # will create multiple processes on localhost.
        cl <- makeCluster(workers)
    } else {
        # The listed hosts must be accessible by the current user using
        # password-less ssh with R installed on all hosts, all 
        # packages installed, and "rscript" is in the default PATH.
        # See docs for makeCluster() for more information.
        cl <- makeCluster(hosts, "PSOCK")
    }
    show(cl)                            # report number of workers and type of parallel implementation
    registerDoParallel(cl)

    foreach (part=1:nSplits) %do% {
        # Returns a standar matrix of the alleles encoded as 0, 1 or 2
        genoNum <- as(genodata[,snp.start[part]:snp.stop[part]], "numeric")

        # Flip the numeric values of genotypes to count minor allele
        if (isTRUE(flip)) genoNum <- flip.matrix(genoNum)

        # For each SNP, concatenate the genotype column to the
        # phenodata and fit a generalized linear model
        rsVec <- colnames(genoNum)
        res <- foreach(snp.name=rsVec, .combine='rbind') %dopar% {
            a <- summary(glm(phenotype~ . - id, family=family, data=cbind(phenodata, snp=genoNum[,snp.name])))
            a$coefficients['snp',]
        }

        # write results so far to a file
        write.table(cbind(rsVec,res), filename, append=TRUE, quote=FALSE, col.names=FALSE, row.names=FALSE)

        cat(sprintf("GWAS SNPs %s-%s (%s%% finished)\n", snp.start[part], snp.stop[part], 100*part/nSplits))
    }

    stopCluster(cl)
  
    return(print("Done."))
}

