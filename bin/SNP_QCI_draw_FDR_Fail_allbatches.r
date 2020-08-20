rm(list=ls())

file_entire_collection  <- commandArgs()[4]
file_perbatch <- commandArgs()[5]

# FDR index indexes a 0-based array at [10e-1, 10e-2, ...]
# => actual threshold = 10e-(index+1)
FDR_index_remove_variants  <- as.integer(commandArgs()[6]) + 1

# ------------------------------------------- #
# -- SNP QCI: HWE across entire collection -- #
# ------------------------------------------- #

fdr.hwe.allbatches  <- read.table(file=paste(file_entire_collection, sep=""), header=T)

fdr  <- -log10(fdr.hwe.allbatches$FDR)
pval <- -log10(fdr.hwe.allbatches$HWE_pval_worstbatchremoved)

png(file=paste(file_entire_collection, ".png", sep=""), width=720, height=720)

plot(fdr, fdr.hwe.allbatches$Fail_worstbatchremoved/1000, xlim=c(1,10), ylim=c(0,40), xlab="-log10(FDR threshold)", ylab="# SNPs failed (x 1,000)", pch=17, col="black", type="b", xaxt="n", cex=1.5)
lines(fdr, fdr.hwe.allbatches$Fail_allbatches/1000, type="b", pch=20, col="black", cex=1.5)
axis(1, at=seq(1,10,by=1), labels=fdr)
axis(3, at=seq(1,10,by=1), labels=round(pval, digits=1), xlab="")
mtext("-log10(p-value threshold)", side=3, line=2.5)
legend("topright", legend=c("All batches","Worst batch removed"), pch=c(20,17), col=c("black","black"), text.col=rep("black", 2), bty="y", inset=.05, cex=1.5)

# FDR threshold
lines(x=c(FDR_index_remove_variants,FDR_index_remove_variants), y=c(0,25), lty=2, col="red")
text(fdr[FDR_index_remove_variants], fdr.hwe.allbatches$Fail_allbatches[FDR_index_remove_variants]/1000, labels=c(fdr.hwe.allbatches$Fail_allbatches[FDR_index_remove_variants]/1000), pos=4)

dev.off()

# --------------------------------- #
# -- SNP QCI: HWE for each batch -- #
# --------------------------------- #

fdr.hwe.perbatch    <- read.table(file=paste(file_perbatch, sep=""), header=T)

fdr  <- -log10(fdr.hwe.perbatch$FDR)

png(file=paste(file_perbatch, ".png", sep=""), width=720, height=720)

plot(fdr, fdr.hwe.perbatch$Fail_2plusbatches/1000, xlim=c(1,10), ylim=c(0,25), xlab="-log10(FDR threshold)", ylab="# SNPs failed (x 1,000)", pch=17, col="black", type="b", xaxt="n", cex=1.5)
axis(1, at=seq(1,10,by=1), labels=fdr)
legend("topright", legend=c("Fail in 2+ batches"), pch=c(17), col=c("black"), text.col=rep("black"), bty="y", inset=.05, cex=1.5)
mtext("corresponding -log(p-values thresholds) different for any batch", side=3, line=1.5)

# FDR threshold
lines(x=c(FDR_index_remove_variants,FDR_index_remove_variants), y=c(0,25), lty=2, col="red")
text(fdr[FDR_index_remove_variants], fdr.hwe.perbatch$Fail_2plusbatches[FDR_index_remove_variants]/1000, labels=c(fdr.hwe.perbatch$Fail_2plusbatches[FDR_index_remove_variants]/1000), pos=4)

dev.off()

