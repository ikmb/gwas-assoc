rm(list=ls())

file.prefix  <- commandArgs()[4]
FDR_index_remove_variants  <- as.integer(commandArgs()[5])

# ------------------------------------------- #
# -- SNP QCII: HF across entire collection -- #
# ------------------------------------------- #

fdr.hf  <- read.table(file=paste(file.prefix, ".1.txt", sep=""), header=T)

fdr  <- -log10(fdr.hf$FDR)
pval <- -log10(fdr.hf$HF_pval_worstbatchremoved_cases)

png(file=paste(file.prefix, ".1.png", sep=""), width=720, height=720)

plot(fdr, fdr.hf$Fail_worstbatchremoved_cases/1000, xlim=c(1,10), ylim=c(0,200), xlab="-log10(FDR threshold)", ylab="# SNPs failed (x 1,000)", pch=17, col="black", type="b", yaxt="n", xaxt="n", cex=1.5)
lines(fdr, fdr.hf$Fail_allbatches_cases/1000, type="b", pch=20, col="black", cex=1.5)
axis(1, at=seq(1,10,by=1), labels=fdr)
axis(2, at=seq(0,200,by=20), labels=seq(0,200,by=20))
axis(3, at=seq(1,10,by=1), labels=round(pval, digits=1), xlab="")
mtext("-log10(p-value threshold)", side=3, line=2.5)
legend("topright", legend=c("All batches","Worst batch removed"), pch=c(20,17), col=c("black","black"), text.col=rep("black", 2), bty="y", inset=.05, cex=1.5)

# FDR threshold
lines(x=c(FDR_index_remove_variants,FDR_index_remove_variants), y=c(0,150), lty=2, col="red")
text(fdr[FDR_index_remove_variants], fdr.hf$Fail_worstbatchremoved_cases[FDR_index_remove_variants]/1000, labels=c(fdr.hf$Fail_worstbatchremoved_cases[FDR_index_remove_variants]/1000), pos=4)

dev.off()

