# ---- map2gene ----
# Returns the subset of SNPs that are within extend.boundary of gene
# using the coords table of gene locations.
map2gene <- function(gene, coords, SNPs, extend.boundary = 5000) {
  coordsSub <- coords[coords$gene == gene,] #Subset coordinate file for spcified gene
  
  coordsSub$start <- coordsSub$start - extend.boundary # Extend gene boundaries
  coordsSub$stop <- coordsSub$stop + extend.boundary
  
  SNPsub <- SNPs[SNPs$position >= coordsSub$start & SNPs$position <= coordsSub$stop &
                 SNPs$chr == coordsSub$chr,] #Subset for SNPs in gene
  
  return(data.frame(SNPsub, gene = gene, stringsAsFactors = FALSE))
}
