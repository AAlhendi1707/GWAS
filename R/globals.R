# ---- globals ---- 
# Customize as needed for file locations

# Modify data.dir to indicate the location of the GWAStutorial files
# Intermediate data files will also be stored in this same location unless you set out.dir
data.dir <- '/Users/ericreed/Desktop/FoulkesLab/SIMFiles'

out.dir <- data.dir                     # may want to write to a separate dir to avoid clutter

# Download files
urlSupport <- "https://www.mtholyoke.edu/courses/afoulkes/Data/GWAStutorial/GWASTutorial_Files.zip"
zipSupport.fn <- sprintf("%s/GWAStutorial_Files.zip", data.dir) 

# Input files
gwas.fn <- lapply(c(bed='bed',bim='bim',fam='fam',gds='gds'), function(n) sprintf("%s/GWAStutorial.%s", data.dir, n))
clinical.fn <- sprintf("%s/GWAStutorial_clinical.csv", data.dir) 
onethou.fn <- lapply(c(info='info',ped='ped'), function(n) sprintf("%s/chr16_1000g_CEU.%s", data.dir, n))
protein.coding.coords.fname <- sprintf("%s/ProCodgene_coords.csv", data.dir)

# Output files
gwaa.fname <- sprintf("%s/GWAStutorialout.txt", out.dir)
gwaa.unadj.fname <- sprintf("%s/GWAStutorialoutUnadj.txt", out.dir)
impute.out.fname <- sprintf("%s/GWAStutorial_imputationOut.csv", out.dir)
CETP.fname <- sprintf("%s/CETP_GWASout.csv", out.dir)

# Working data saved between each code snippet so each can run independently.
# Use save(data, file=working.data.fname(num))
working.data.fname <- function(num) { sprintf("%s/working.%s.Rdata", out.dir, num) }


