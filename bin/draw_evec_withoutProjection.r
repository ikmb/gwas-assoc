rm(list=ls())

file.prefix  <- commandArgs()[4]
file.prefix.withoutHapMap  <- commandArgs()[5]

eigenval  <- read.table(file=paste(file.prefix, ".eval", sep=""), col.names=c("eigenvalues"))
eigenvect <- read.table(file=paste(file.prefix.withoutHapMap, ".pca.evec", sep=""), skip=1, col.names=c("indivID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "batch"))

# ---------------- #
# -- draw PC1-4 -- #
# ---------------- #

png(file=paste(file.prefix.withoutHapMap, ".4PCpairs.png", sep=""), width=900, height=720)
thelabels <- paste("PC", 1:4, "\n", sep="")
pairs(eigenvect[,2:5], col=eigenvect$batch, labels=thelabels)
dev.off()


# -------------------------------------------- #
# -- draw PC1-2 with color code reg batches -- #
# -------------------------------------------- #

min_x = min(eigenvect$PC1)*1.1
max_x = max(eigenvect$PC1)*1.1
min_y = min(eigenvect$PC2)*1.1
max_y = max(eigenvect$PC2)*1.1

# --------------- #
# -- variant 1 -- #
# --------------- #

batches_code <- unique(eigenvect$batch)
color_code <- rgb(runif(length(batches_code)),runif(length(batches_code)),runif(length(batches_code))) 

# -- plot symbols -- #
# batches get '0'
pch_code <- as.vector(unique(eigenvect$batch))
pch_code[pch_code != "HAPMAPCEU" & pch_code != "HAPMAPJPT" & pch_code != "HAPMAPCHB" & pch_code!="HAPMAPYRI"] <- 1
pch_code <- as.integer(pch_code)

# determine color vector for each entry of eigenvect
color_vect <- color_code[match(eigenvect$batch, batches_code)]
# determine symbol vector for each entry of eigenvect
pch_vect <- pch_code[match(eigenvect$batch, batches_code)]

# draw plot
png(file=paste(file.prefix.withoutHapMap, ".2PC.png", sep=""), width=900, height=720)
plot(eigenvect[2:3], xlim=c(min_x,max_x),ylim=c(min_y,max_y), pch=pch_vect, col=color_vect)
legend("topleft", legend=unique(eigenvect$batch), cex=0.7, pch=pch_code, col=color_code, text.col=rep("black", length(unique(eigenvect$batch))), bty="n")
dev.off()



# --------------- #
# -- variant 2 -- #
# --------------- #

#png(file="PS_AS_IBD_PSC_SNPQCI_pruned_hapmap2.2PC.png", width=900, height=720)
#plot(eigenvect[2:3], xlim=c(min_x,max_x),ylim=c(min_y,max_y), pch=rep(1, length(unique(eigenvect$batch))), col=eigenvect$batch)
#####warnings()
#legend("topleft", legend=unique(eigenvect$batch), cex=0.7, pch=rep(1, length(unique(eigenvect$batch))), col=as.numeric(unique(eigenvect$batch)), text.col=rep("black", length(unique(eigenvect$batch))), bty="n")
##legend("topleft", legend=unique(eigenvect$batch), cex=0.7, pch=rep(1, length(unique(eigenvect$batch))), text.col=as.numeric(unique(eigenvect$batch)), col=as.numeric(unique(eigenvect$batch)), bty="n")
#####warnings()
#dev.off()





##### only 0-25 symbols available, otherwise this creates a warning, look at with warnings()
####pch_vector <- c(as.numeric(eigenvect$batch)%%26)
####png(file="PS_AS_IBD_PSC_SNPQCI_pruned_hapmap2.2PC.png", width=900, height=720)
####plot(eigenvect[2:3], xlim=c(min_x,max_x),ylim=c(min_y,max_y), pch=pch_vector, col=eigenvect$batch)
####warnings()
####dev.off()
