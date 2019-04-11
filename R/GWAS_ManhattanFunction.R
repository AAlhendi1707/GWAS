# ---- manhattan ----
# Receives a data.frame of SNPs with Neg_logP, chr, position, and type.
# Plots Manhattan plot with significant SNPs highlighted.
GWAS_Manhattan <- function(GWAS, col.snps=c("black","gray"),
                           col.detected=c("black"), col.imputed=c("blue"), col.text="black",
                           title="GWAS Tutorial Manhattan Plot", display.text=TRUE,
                           bonferroni.alpha=0.05, bonferroni.adjustment=1000000,
                           Lstringent.adjustment=10000) {
    
    bonferroni.thresh <- -log10(bonferroni.alpha / bonferroni.adjustment)
    Lstringent.thresh <- -log10(bonferroni.alpha / Lstringent.adjustment)
    xscale <- 1000000

    manhat <- GWAS[!grepl("[A-z]",GWAS$chr),]
    
    #sort the data by chromosome and then location
    manhat.ord <- manhat[order(as.numeric(manhat$chr),manhat$position),]
    manhat.ord <- manhat.ord[!is.na(manhat.ord$position),]
    
    ##Finding the maximum position for each chromosome
    max.pos <- sapply(1:21, function(i) { max(manhat.ord$position[manhat.ord$chr==i],0) })
    max.pos2 <- c(0, cumsum(max.pos))                  
    
    #Add spacing between chromosomes
    max.pos2 <- max.pos2 + c(0:21) * xscale * 10
    
    #defining the positions of each snp in the plot
    manhat.ord$pos <- manhat.ord$position + max.pos2[as.numeric(manhat.ord$chr)]
    
    # alternate coloring of chromosomes
    manhat.ord$col <- col.snps[1 + as.numeric(manhat.ord$chr) %% 2]
    
    # draw the chromosome label roughly in the middle of each chromosome band
    text.pos <- sapply(c(1:22), function(i) { mean(manhat.ord$pos[manhat.ord$chr==i]) })
    
    # Plot the data
    plot(manhat.ord$pos[manhat.ord$type=="typed"]/xscale, manhat.ord$Neg_logP[manhat.ord$type=="typed"],
         pch=20, cex=.3, col= manhat.ord$col[manhat.ord$type=="typed"], xlab=NA,
         ylab="Negative Log P-value", axes=F, ylim=c(0,max(manhat$Neg_logP)+1))
    #Add x-label so that it is close to axis
    mtext(side = 1, "Chromosome", line = 1.25)
    
    points(manhat.ord$pos[manhat.ord$type=="imputed"]/xscale, manhat.ord$Neg_logP[manhat.ord$type=="imputed"],
           pch=20, cex=.4, col = col.imputed)
    
    points(manhat.ord$pos[manhat.ord$type=="typed"]/xscale, manhat.ord$Neg_logP[manhat.ord$type=="typed"],
           pch=20, cex=.3, col = manhat.ord$col[manhat.ord$type=="typed"])
    
    axis(2)
    abline(h=0)
    
    SigNifSNPs <- as.character(GWAS[GWAS$Neg_logP > Lstringent.thresh & GWAS$type=="typed", "SNP"])
    
    #Add legend
    legend("topright",c("Bonferroni corrected threshold (p = 5E-8)", "Candidate threshold (p = 5E-6)"),
           border="black", col=c("gray60", "gray60"), pch=c(0, 0), lwd=c(1,1),
           lty=c(1,2), pt.cex=c(0,0), bty="o", cex=0.6)
    
    #Add chromosome number
    text(text.pos/xscale, -.3, seq(1,22,by=1), xpd=TRUE, cex=.8)
    
    #Add bonferroni line
    abline(h=bonferroni.thresh, untf = FALSE, col = "gray60")
    
    #Add "less stringent" line
    abline(h=Lstringent.thresh, untf = FALSE, col = "gray60", lty = 2 )
    
    #Plotting detected genes
    #Were any genes detected?
    if (length(SigNifSNPs)>0){

        sig.snps <- manhat.ord[,'SNP'] %in% SigNifSNPs
        
        points(manhat.ord$pos[sig.snps]/xscale,
               manhat.ord$Neg_logP[sig.snps],
             pch=20,col=col.detected, bg=col.detected,cex=0.5)
      
      text(manhat.ord$pos[sig.snps]/xscale,
           manhat.ord$Neg_logP[sig.snps],
           as.character(manhat.ord[sig.snps,1]), col=col.text, offset=1, adj=-.1, cex=.5)
    }
}
