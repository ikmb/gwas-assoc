rm(list=ls())

file.eigenvect   <- commandArgs()[4]
file.annotation  <- commandArgs()[5]

eigenvect <- read.table(file=file.eigenvect, skip=1, col.names=c("FID", "IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "batch"))
annotation       <- read.table(file=file.annotation, header=T, sep="\t")

#####matched <- match(as.character(eigenvect$indivID), as.character(annotation$individualID))

# ------------------------------ #
# -- draw histograms of PC1-6 -- #
# ------------------------------ #

library('sm')
library('plotrix')
#color_code <- c(1:(length(levels(eigenvect$batch))-1)) 
color_code <- c(1:(length(levels(eigenvect$batch)))) 
#legend_code=c("Control", "CD", "UC", "AS", "PS", "PSC")
legend_code= sort(unique(eigenvect$batch))

for (i in seq(3,12)) {
    png(file=paste(file.eigenvect, ".histPC", i-2, ".png", sep=""), width=900, height=720)
    par(mfrow=c(2,2))

    a = eigenvect[eigenvect$batch == legend_code[1], i]
    b = eigenvect[eigenvect$batch == legend_code[2], i]
    #a = eigenvect[eigenvect$batch == "Control", i]
    #b = eigenvect[eigenvect$batch == "CD", i]
    #c = eigenvect[eigenvect$batch == "UC", i]
    #d = eigenvect[eigenvect$batch == "AS", i]
    #e = eigenvect[eigenvect$batch == "PS", i]
    #f = eigenvect[eigenvect$batch == "PSC", i]
    a = a[!is.na(a)] 
    b = b[!is.na(b)]
    #c = c[!is.na(c)]
    #d = d[!is.na(d)]
    #e = e[!is.na(e)]
    #f = f[!is.na(f)]

    print(c(min(a),max(a)))

    #multhist(list(a,b,c,d,e,f), col=color_code, xlab=paste("PC", i-2, sep=""), ylab="Counts", freq=T)
    multhist(list(a,b), col=color_code, xlab=paste("PC", i-2, sep=""), ylab="Counts", freq=T)
    legend("topright", legend=legend_code, fill=color_code)

    #multhist(list(a,b,c,d,e,f), col=color_code, xlab=paste("PC", i-2, sep=""), ylab="Density", freq=F)
    multhist(list(a,b), col=color_code, xlab=paste("PC", i-2, sep=""), ylab="Density", freq=F)
    legend("topright", legend=legend_code, fill=color_code)

    #disease.f <- factor(eigenvect$batch, levels=c("Control", "CD", "UC", "AS", "PS", "PSC"), labels=c("Control", "CD", "UC", "AS", "PS", "PSC"))
    disease.f <- factor(eigenvect$batch, levels=legend_code, labels=legend_code)
    colfill <- c(1:(1+length(levels(disease.f)))) 
    #sm.density.compare(eigenvect[,i], eigenvect$batch, xlab=paste("PC", i-2, sep=""))
    sm.density.compare(eigenvect[,i], eigenvect$batch, xlab=paste("PC", i-2, sep=""), col=color_code)
    legend("topright", levels(disease.f), fill=colfill)
    
    dev.off()
}
