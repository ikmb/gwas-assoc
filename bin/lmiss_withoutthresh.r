###############################
###  2010, David Ellinghaus ###
###############################

rm(list=ls())

file.lmiss = commandArgs()[4]

lmiss = read.table(file.lmiss, as.is=T, header=T)

png(file=paste(file.lmiss, ".1.png", sep=""), width=960, height=960)

hist(lmiss$F_MISS, xlim=c(0,1.0), main="All SNPs", xlab="SNP missingness", ylab="Number of SNPs", col="grey", border="black", axes=F)
axis(1, at=seq(0,1.0,0.1), tick=T)
axis(2, at=seq(0,1000000,20000), tick=T)

dev.off()
