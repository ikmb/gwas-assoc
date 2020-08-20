
args = commandArgs(TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript example_script.R <filename><filename_sampleInfo><number>")
}

## Save arguments in variables
filename= args[1]
preQCIMDS_1kG_sample = args[2]
#pci = as.numeric(args[3]) # factor interquartile range
pci = 5 # factor interquartile range

options(stringsAsFactors=F)


pdf(file=paste(filename, ".2PC.pdf", sep=""))
pcs= read.table(paste(filename,".pca.evec",sep=""),h=T,col.names=c("FID","IID","PC1","PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "group"),sep="\t")
print(head(pcs))

samples = read.table(preQCIMDS_1kG_sample,h=T)
samples = samples[match(pcs$IID,samples$ID),]


AFR= c("ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI")
col.AFR=colors()[369:375]
AMR = c("CLM","MXL","PEL","PUR")
col.AMR =colors()[8:11]
EAS = c("CDX","CHB", "CHS", "JPT", "KHV")
col.EAS = colors()[589:593]
EUR=c("CEU","FIN","GBR","IBS","TSI")
col.EUR =colors()[382:386]
SAS=c("BEB","GIH","ITU","PJL","STU")
col.SAS = colors()[493:497]

all = data.frame(pop  = c(AFR,AMR, EAS,EUR,SAS),
                 group= c(rep("AFR", length(AFR)), rep("AMR", length(AMR)), 
                          rep("EAS", length(EAS)), rep("EUR", length(EUR)),
                          rep("SAS", length(SAS))), 
                 col =  c(col.AFR,col.AMR, col.EAS,col.EUR,col.SAS))



col = as.matrix(all$col[match(samples$POP,all$pop)])
pch=rep(16, length(col))
pch[is.na(col)] =3
col[is.na(col)] = colors()[312] 

data =data.frame(pcs[,3:4],
                 group=samples$GROUP, 
                 pop=samples$POP, 
                 col = as.matrix(col), 
                 pch=as.matrix(pch))



medians = cbind( tapply(data[,1],data$group, function(x){median(x,na.rm=T)}),
                 tapply(data[,2],data$group, function(x){median(x,na.rm=T)}))

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
tmp = data[is.na(data$group),]

rest = data[!is.na(data$group),]

plot(rest[,1:2], col=as.matrix(rest$col), pch=rest$pch)
points(tmp[,1:2], col=as.matrix(tmp$col), pch=tmp$pch )
text(medians,labels=rownames(medians))
legend("topright", inset=c(-0.2,0),legend=all$pop, col = as.matrix(all$col), pch=16)

# plot(data[,1:2], col = col, pch = pch, xlab = "PC1", ylab = "PC2")

data = pcs[is.na(data$group),]
pc_one_median = median(data$PC1, na.rm=TRUE ) # median first pc
pc_one_iqr = IQR(data$PC1, na.rm=TRUE) # inter quartile range first pc
pc_two_median =  median(data$PC2, na.rm=TRUE ) # median second pc
pc_two_iqr = IQR(data$PC2, na.rm=TRUE) # inter quartiel range pc

# Plot 1

y = rbind(range(data$PC2), c(pc_two_median - 5*pc_two_iqr,pc_two_median + 5*pc_two_iqr)*1.1)
x = rbind(range(data$PC1), c(pc_one_median - 5*pc_one_iqr, pc_one_median + 5*pc_one_iqr)*1.1) 

xl = c(min(x[,1]), max(x[,2]))
yl = c(min(y[,1]), max(y[,2]))


# mod DE
rect(pc_one_median - 2*pc_one_iqr, pc_two_median - 2*pc_two_iqr, 
     pc_one_median + 2*pc_one_iqr, pc_two_median +2*pc_two_iqr, 
     border="GREEN", lty=2)

rect(pc_one_median - 3*pc_one_iqr, pc_two_median - 3*pc_two_iqr, 
     pc_one_median + 3*pc_one_iqr, pc_two_median +3*pc_two_iqr, 
     border="BLUE", lty=2)
rect(pc_one_median - 5*pc_one_iqr, pc_two_median - 5*pc_two_iqr, 
     pc_one_median + 5*pc_one_iqr, pc_two_median +5*pc_two_iqr, 
     border="RED", lty=2)

legend("bottomright", c("median ± 2*IQR", "median ± 3*IQR", "median ± 5*IQR"), 
       col=c("green", "blue", "red"), 
       lty=c(1,1,1), 
       cex=0.8, 
       bg = "transparent")

out_x = data$PC1 >=pc_one_median + pci*pc_one_iqr | data$PC1 <=pc_one_median - pci*pc_one_iqr # find data within threshold
out_y = data$PC2 >=pc_two_median + pci*pc_two_iqr | data$PC2 <=pc_two_median - pci*pc_two_iqr



write.table(data[out_x|out_y,1:2],paste(filename, ".fail-pca-1KG-qc.txt", sep=""), quote=F, col.names=F, row.names=F)


# dev.off()
graphics.off()
q()
n
