rm(list=ls())

file.eigenvect   <- commandArgs()[4]
file.annotation  <- commandArgs()[5]
#eigenvect        <- read.table(file=file.eigenvect, skip=1, col.names=c("indivID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "batch"))
eigenvect <- read.table(file=file.eigenvect, skip=1, col.names=c("indivID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", 
      "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20",
      "PC21", "PC22", "PC23", "PC24", "PC25", "PC26", "PC27", "PC28", "PC29", "PC30",
      "PC31", "PC32", "batch"))
annotation       <- read.table(file=file.annotation, header=T, sep="\t")

#eigenvect  <- read.table(file="PS_AS_IBD_PSC_SampleQCI_hapmap2_10PC.withoutHapMap.pca.evec", skip=1, col.names=c("indivID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "batch"))
#annotation <- read.table(file="../SNPQCII/PS_AS_IBD_PSC_SNPQCII_annotation.txt", header=T, sep="\t")
#####matched <- match(as.character(eigenvect$indivID), as.character(annotation$individualID))

# ------------------------------ #
# -- draw histograms of PC1-6 -- #
# ------------------------------ #

library('sm')
library('plotrix')
color_code <- c(1:(length(levels(annotation$diagnosis))-1)) 
legend_code=c("Control", "CD", "UC", "AS", "PS", "PSC")

for (i in seq(2,11)) {
    #png(file=paste("PS_AS_IBD_PSC_SampleQCI_hapmap2_10PC.withoutHapMap.pca.evec", ".histPC", i-1, ".png", sep=""), width=900, height=720)
    png(file=paste(file.eigenvect, ".histPC", i-1, ".png", sep=""), width=900, height=720)
    par(mfrow=c(2,2))

    a = eigenvect[annotation$diagnosis == "Control", i]
    b = eigenvect[annotation$diagnosis == "CD", i]
    c = eigenvect[annotation$diagnosis == "UC", i]
    d = eigenvect[annotation$diagnosis == "AS", i]
    e = eigenvect[annotation$diagnosis == "PS", i]
    f = eigenvect[annotation$diagnosis == "PSC", i]
    a = a[!is.na(a)] 
    b = b[!is.na(b)]
    c = c[!is.na(c)]
    d = d[!is.na(d)]
    e = e[!is.na(e)]
    f = f[!is.na(f)]

    print(c(min(a),max(a)))

    multhist(list(a,b,c,d,e,f), col=color_code, xlab=paste("PC", i-1, sep=""), ylab="Counts", freq=T)
    legend("topright", legend=legend_code, fill=color_code)

    multhist(list(a,b,c,d,e,f), col=color_code, xlab=paste("PC", i-1, sep=""), ylab="Density", freq=F)
    legend("topright", legend=legend_code, fill=color_code)

    disease.f <- factor(annotation$diagnosis, levels=c("Control", "CD", "UC", "AS", "PS", "PSC"), labels=c("Control", "CD", "UC", "AS", "PS", "PSC"))
    colfill <- c(1:(1+length(levels(disease.f)))) 
    sm.density.compare(eigenvect[,i], annotation$diagnosis, xlab=paste("PC", i-1, sep=""))
    legend("topright", levels(disease.f), fill=colfill)
    
    dev.off()
}
