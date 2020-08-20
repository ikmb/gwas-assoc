##############################################################################################
###  Genotype-phenotype overlap (GPO): heterozygosity per individual                       ###
##                                     filter when +- 5 standard deviation (SD) from mean  ###
###  2009, David Ellinghaus, d.ellinghaus@ikmb.uni-kiel.de                                 ###
##############################################################################################

## read arguments from command line
file.het <- commandArgs()[4]

## - read data - #
ds.het <- read.table(file=file.het, as.is=T, header=T)
dim(ds.het)

heterozygosity <- 1 - ds.het[,3]/ds.het[,5] 

length(heterozygosity)
summary(heterozygosity)
het.mean <- mean(heterozygosity)
# heterozygosity +- 5SD way from mean
het.sd_5times <- 5 * sd(heterozygosity, na.rm = FALSE)
print(het.sd_5times)

upper_bound  <- het.mean + het.sd_5times
downer_bound <- het.mean - het.sd_5times

output <- paste(file.het, ".outlier.txt",sep="")
sink (file=output)

for (i in 1:length(heterozygosity)) {
  if (heterozygosity[i] < downer_bound || heterozygosity[i] > upper_bound) {  
    cat(ds.het[i,1], ds.het[i,2], heterozygosity[i], "\n")
  }
} 
sink ()

###################################################################
###  END OF SCRIPT  ###############################################
###################################################################
