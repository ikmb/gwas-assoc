rm(list=ls())

file.maf = commandArgs()[4]
fraction_maf.thresh = as.numeric(commandArgs()[5])

maf = read.table(file.maf, as.is=T, header=T)
maf.log10_MAF = log10(maf$MAF)

png(file=paste(file.maf, ".logscale.2.png", sep=""), width=960, height=960)

hist(maf.log10_MAF, xlim=c(-4,0), main="QCed SNPs", xlab="MAF", ylab="Number of SNPs", col="grey", border="black", axes=F)
axis(1, at=c(-4,-3,-2,log10(fraction_maf.thresh),-1,0), labels=c(0.0001,0.001,0.01,fraction_maf.thresh,0.1,1))
axis(2, at=seq(0,2000000,20000), tick=T)
abline(v=log10(fraction_maf.thresh), col="red", lty=2)

dev.off()
