rm(list=ls())

file.eigenvect   <- commandArgs()[4]
file.annotation  <- commandArgs()[5]
eigenvect        <- read.table(file=file.eigenvect, skip=1, col.names=c("FID", "IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "batch"))
#eigenvect <- read.table(file=file.eigenvect, skip=1, col.names=c("FID", "IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", 
#      "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20",
#      "PC21", "PC22", "PC23", "PC24", "PC25", "PC26", "PC27", "PC28", "PC29", "PC30",
#      "PC31", "PC32", "batch"))
annotation       <- read.table(file=file.annotation, header=T, sep="\t")

# ------------------------------ #
# -- draw histograms of PC1-6 -- #
# ------------------------------ #

library('sm')
library('plotrix')
#color_code <- c(1:(length(levels(annotation$diagnosis))-1)) 
color_code <- c(1:(length(levels(annotation$diagnosis)))) 
#legend_code=c("Control", "CD", "UC", "AS", "PS", "PSC")
legend_code= sort(unique(annotation$diagnosis))

for (i in seq(3,8)) {
    png(file=paste(file.eigenvect, ".histPC", i-2, ".png", sep=""), width=900, height=720)
    par(mfrow=c(2,2))

    a = eigenvect[annotation$diagnosis == legend_code[1], i]
    b = eigenvect[annotation$diagnosis == legend_code[2], i]
    #a = eigenvect[annotation$diagnosis == "Control", i]
    #b = eigenvect[annotation$diagnosis == "CD", i]
    #c = eigenvect[annotation$diagnosis == "UC", i]
    #d = eigenvect[annotation$diagnosis == "AS", i]
    #e = eigenvect[annotation$diagnosis == "PS", i]
    #f = eigenvect[annotation$diagnosis == "PSC", i]
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

    #disease.f <- factor(annotation$diagnosis, levels=c("Control", "CD", "UC", "AS", "PS", "PSC"), labels=c("Control", "CD", "UC", "AS", "PS", "PSC"))
    disease.f <- factor(annotation$diagnosis, levels=legend_code, labels=legend_code)
    colfill <- c(1:(1+length(levels(disease.f)))) 
    #sm.density.compare(eigenvect[,i], annotation$diagnosis, xlab=paste("PC", i-2, sep=""))
    sm.density.compare(eigenvect[,i], annotation$diagnosis, xlab=paste("PC", i-2, sep=""), col=color_code)
    legend("topright", levels(disease.f), fill=colfill)
    
    dev.off()
}
