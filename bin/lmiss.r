###############################
###  2010, David Ellinghaus ###
###############################

rm(list=ls())

file.lmiss = commandArgs()[4]
fraction_miss.thresh = as.numeric(commandArgs()[5])

lmiss = read.table(file.lmiss, as.is=T, header=T)

png(file=paste(file.lmiss, ".2.png", sep=""), width=960, height=960)

hist(lmiss$F_MISS, xlim=c(0,fraction_miss.thresh+0.01), main="All SNPs", xlab="SNP missingness", ylab="Number of SNPs", col="grey", border="black", axes=F)
axis(1, at=seq(0,fraction_miss.thresh+0.01,0.01), tick=T)
axis(2, at=seq(0,2000000,20000), tick=T)
abline(v=fraction_miss.thresh, col="red", lty=2)

dev.off()
