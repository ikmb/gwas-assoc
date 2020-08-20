###############################
###  2010, David Ellinghaus ###
###############################

rm(list=ls())

file.het   = commandArgs()[4]
file.imiss = commandArgs()[5]

# heterozygosity +- 5SD way from mean
sd.times = 5

imiss = read.table(file.imiss,header=T)
het   = read.table(file.het,header=T)

het = 1 - het$O.HOM./het$N.NM.

png(file=paste(file.het, ".1.png", sep=""), width=960, height=960)

dcols = densCols(imiss$F_MISS,het)
plot(imiss$F_MISS, het, xlim=c(0,1), ylim=c(0,0.5), xlab="Individual missingness", ylab="Heterozygosity rate", col=dcols, pch=20, axes=F)
axis(1, at=seq(0,1,0.1), tick=T)
axis(2, at=seq(0,0.5,0.05), tick=T)
abline(h=mean(het)-(sd.times*sd(het)), lty=2, col="red")
abline(h=mean(het)+(sd.times*sd(het)), lty=2, col="red")

dev.off()
