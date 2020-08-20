# 2018-08-27 f.degenhardt Original version
# 2018-08-27 j.kaessens   Remove library path for use with singularity
#            "            Add 1000 Genomes samples as extra parameter

#######################################
# Calculate outliers
#######################################
# .libPaths("~/R/x86_64-redhat-linux-gnu-library/3.2")

calc.outlier =function(data,pci){
  iqr = apply(data[,-(1:2)],2,IQR)
  tmp = c()
  dist =apply(data[,-(1:2)], 2, function(x){return(x-median(x))})
  for(i in c(1)){

    tmp = cbind(tmp,help.calc(dist[,i], dist[,i+1], iqr[c(i,i+1)],pci))

  }
  out = apply(tmp,1,any)
  return(out)
}

help.calc = function(x,y,iqr,pci){
  dist = sqrt(x^2 +y^2)
  out = dist > pci*iqr[1] &  dist > pci*iqr[2]
  return(out)
}

get.origin = function(train,test){
  train$pop= as.matrix(train$pop)
  print(dim(train))
  
  train.median = cbind(tapply(train[,1], train$pop, median),
                       tapply(train[,2], train$pop, median))
  for(i in 1:nrow(test)){
    sqrt = sqrt((train.median[,1] -test[i,1])^2 +(train.median[,2] -test[i,2])^2)
    test$assigned[i] = names(sort(sqrt))[1]

  }
 return(test$assigned)

}

args = commandArgs(TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript example_script.R <filename> <number> <1kg-samples>")
}

## Save arguments in variables
filename= args[1]
pci = as.numeric(args[2])
onekg_samples = args[3]

options(stringsAsFactors=F)


print("Reading pca.evec")
pcs= read.table(paste(filename, ".pca.evec", sep=""),h=T,col.names=c("FID","IID","PC1","PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "group"),sep="\t")
print("Done!")
print(head(pcs))

print("Reading 1kg tables")
samples = read.table(onekg_samples,h=T)
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

pc_columns = dim(pcs)[2]

data =data.frame(pcs[,3:(pc_columns-1)],
                 group=samples$GROUP,
                 pop=samples$POP,
                 col = as.matrix(col),
                 pch=as.matrix(pch))
data$assigned=get.origin(data[!is.na(data$pop),],data)
# write.table(cbind(pcs[,1:2], data$assigned), paste(filename,".origin", sep=""), quote=F, row.names=F)

medians = apply(data[,1:10],2,function(x){
  x=tapply(x,data$group, function(x){return(median(x,na.rm=T))})
  return(x)})




plot.data = function(data, medians, cols, pci,write=F){
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  tmp = data[is.na(data$group),]

  rest = data[!is.na(data$group),]
  x = data[,cols[1]]
  y=  data[,cols[2]]
  plot(rest[,cols], col=as.matrix(rest$col), pch=rest$pch,
       xlim=c(min(x),max(x)),ylim=c(min(y), max(y)))
  points(tmp[,cols], col=as.matrix(tmp$col), pch=tmp$pch )
  text(medians[,cols],labels=rownames(medians))
  legend("topright", inset=c(-0.2,0),legend=all$pop, col = as.matrix(all$col), pch=16)

  # plot(data[,1:2], col = col, pch = pch, xlab = "PC1", ylab = "PC2")
  data = pcs[is.na(data$group),]


  #out =calc.outlier(data[3:(dim(data)[2]-1)], pci)
  #if(length(cols)==2){
  #points(data[,3:(dim(data)[2]-1)][out,cols], col="red", pch=16)}
  #if(write==T){
  #  write.table(data[out,1:2],"fail-pca-1KG-qc.txt", quote=F, col.names=F, row.names=F)
  #}

}

print("DONE1")
print(medians)
#print(data[1:10],1)
pdf(paste(filename,"_pca.pdf", sep=""))
plot.data(data, medians, 1:10,pci, T)
dev.off()
print("DONE2")
library(R.utils)

pdf(paste(filename,"_pca_PCs.pdf", sep=""))
plot.data(data, medians, 1:2,pci, F)
plot.data(data, medians, 3:4,pci, F)
plot.data(data, medians, 5:6,pci, F)
plot.data(data, medians, 7:8,pci, F)
plot.data(data, medians, 9:10,pci, F)
dev.off()
q()
n
