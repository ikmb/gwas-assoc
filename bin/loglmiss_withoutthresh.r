###############################
###  2010, David Ellinghaus ###
###############################

rm(list=ls())

file.lmiss = commandArgs()[4]

lmiss = read.table(file.lmiss, as.is=T, header=T)
lmiss.log10F_MISS = log10(lmiss$F_MISS)

png(file=paste(file.lmiss, ".logscale.1.png", sep=""), width=960, height=960)

hist(lmiss.log10F_MISS, xlim=c(-4,0), main="All SNPs", xlab="SNP missingness", ylab="Number of SNPs", col="grey", border="black", axes=F)
axis(1, at=c(-4,-3,-2,-1,0), labels=c(0.0001,0.001,0.01,0.1,1))
axis(2, at=seq(0,2000000,20000), tick=T)

dev.off()
