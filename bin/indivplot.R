## Frauke Degenhardt 05092013 - Qualtiy control of gwas (control) data according to Anderson et al. 2010, Nature Protocols
## Jan Kässens       31082018 - Adapation for IKMB QC Pipeline integration
#  Jan Kässens       11012019 - Disabled sex check and outlier tables by request
## Vizualisation of STEP 1: qualitiy control on sample

# Read file
args = commandArgs(TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript example_script.R <file> <x-ped> <y-ped> <imiss> <het>")
}

print(args)

# .libPaths("~/R/x86_64-redhat-linux-gnu-library/3.2")
## Save arguments in variables
file= args[1]
imissfile= args[4]
hetfile=args[5]

library("R.utils")
library("SNPRelate")
library("stats")

ped_x = args[2]
ped_y = args[3]

print("Vor laden 1")

options(stringsAsFactors=F)
X = read.table(ped_x,h=F, colClass="character")
Y = read.table(ped_y,h=F,colClass="character")
fam= X[,1:6]
colnames(fam)=c("FID","IID","","","SEX","PHENO")
fam = data.frame(fam)


print("nach laden 1")

X =X[,-(1:6)]
X = t(apply(X,1,function(x){tmp = x[seq(1,length(x),2)] == x[seq(2,length(x),2)] & x[seq(1,length(x),2)]!="0";
                                  x[tmp]=0; x[!tmp]=2; return(as.numeric(x))}))

Y =Y[,-(1:6)]
Y = t(apply(Y,1,function(x){tmp = x[seq(1,length(x),2)] == x[seq(2,length(x),2)] & x[seq(1,length(x),2)]!="0";
                            x[tmp]=0; x[!tmp]=2; return(as.numeric(x))}))


sumsX = rowSums(X)
sumsY = rowSums(Y)


borderX = (median(sumsX[fam$SEX==2]) + median(sumsX[fam$SEX==1]))/2
borderY = (median(sumsY[fam$SEX==2]) + median(sumsY[fam$SEX==1]))/2

fam$out[fam$SEX==2] = sumsX[fam$SEX==2] < borderX | sumsY[fam$SEX==2] < borderY
fam$out[fam$SEX==1] =  sumsX[fam$SEX==1] > borderX | sumsY[fam$SEX==1] > borderY



write.table(data.frame(FID=as.matrix(fam$FID), IID=as.matrix(fam$IID))[fam$out,], "fail-sexcheck-qc.txt",
            quote = F, row.names = F, col.names =T)

# Color for plotting
print("nach write.table")

pch_sexcheck = as.numeric(fam$SEX)
col_sexcheck =  adjustcolor("grey", alpha=0.5)



# GENOTPYE

imiss=read.table(imissfile,h=T, fill=TRUE)
print("genotype geladen")

# HETEROZYGOSITY

het=read.table(hetfile,h=T,fill=TRUE)
het$meanHet = (het$N.NM. - het$O.HOM.)/het$N.NM.


# IBD

pdf(paste(file,".indiv.pdf",sep = ""))

subplots(nrow=2, ncol=2)
if(length(list.files(pattern = "RData"))!=0){
load(list.files(pattern = "RData"))
ibd.coeff = try(snpgdsIBDSelection(ibd))
}else{
  ibd.coeff="try-error"
  class(ibd.coeff) = "try-error"
}


# Plot



# Plot 1 - GENOTPYE
# Grenze im Plot

maximiss=ifelse(max(imiss$F_MISS)<0.11,0.11, max(imiss$F_MISS))

plot(imiss$F_MISS,
     xlab="Nr. Person",
     ylab="Proportion of missing genotpye",
     main="",
     ylim=c(0,maximiss))

# LINES

abline(h=0.1, col="RED",lty=2)
abline(h=0.05, col="BLUE",lty=2)
abline(h=0.03, col="GREEN",lty=2)
legend("topright",c("3%", "5%", "10%"),
       col=c("green","blue","red"),
       lty = c(1,1,1),
       cex=0.8, bg="transparent")

# Plot 2 - HETEROZYGOSITY
het_mean=mean(het$meanHet)
het_sd= sd(het$meanHet)

plot(het$meanHet,
    xlab="Nr. Person",
    ylab="Heterozygosity rate",
    main="",
    ylim=c(het_mean-10*het_sd, het_mean+10*het_sd))
abline(h=c(het_mean-3*het_sd,
           het_mean+3*het_sd,
           het_mean-2*het_sd,
           het_mean-2*het_sd,
           het_mean-1*het_sd,
           het_mean+1*het_sd),
       col = c(rep("red",2), rep("blue",2), rep ("green",2)),
       lty = rep(2, 6))
legend("topright",
       c("68%", "95%", "99.7%"),
       col=c("green","blue","red"), lty = c(1, 1, 1),
       cex=0.8,
       bg = "transparent")

# Plot 4 - IBD

#if(class(ibd.coeff)!="try-error"){
#  hist(ibd.coeff$kinship,ylim=c(0,100),breaks=100,xlab="Estimated IBD (MLE)",main="")
#
#  ## Grenze
#  abline(v=0.185, lty=2, col="RED")
#  legend("topright", "0.185", col=c("RED"), lty = 1,  cex=0.8, bg="transparent")
#}

# Plot 5 - SEXCHECK
#plot(sumsX, sumsY,
#     col = col_sexcheck,
#     pch = pch_sexcheck,
#     xlab = "sum (homozygous counts X)",
#     ylab = "sum (homozygous counts Y)",
#     cex= 0.8, main="")
#points(sumsX[fam$out], sumsY[fam$out], col = "red", pch = pch_sexcheck[fam$out], bg= "red")
#
#legend("topright",
#       pch = c(1,2) ,
#       col = "grey" ,
#       c("male", "female"),
#       cex = 0.8,
#       bg = "transparent")

dev.off()
quit()
n
