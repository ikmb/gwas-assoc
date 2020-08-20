###############################
###  2009, David Ellinghaus ###
###############################

rm(list=ls())

file.results   <- commandArgs()[4]
disease_name   <- commandArgs()[5]

gwas <- read.table(file.results,header=T,stringsAsFactor=F)

# manhattan
gwas <- gwas[order(gwas$CHR),]
pos.plot <- cpos.abs<-as.numeric(gwas$BP)
cpos.min <- tapply(as.numeric(gwas$BP),gwas$CHR,min)
cpos.max <- tapply(as.numeric(gwas$BP),gwas$CHR,max)
cpos.max.abs <- cpos.max-cpos.min

pos.plot[gwas$CHR==1] <- pos.plot[gwas$CHR==1]-cpos.min[1]+1
for(i in 2:22)
{
  c.ind <- gwas$CHR==i
  pos.plot[c.ind==T] <- sum(cpos.max[1:(i-1)])+pos.plot[c.ind==T]-cpos.min[i]+1
}
#colplot <- ifelse(gwas$CHR%in%c(2,4,6,8,10,12,14,16,18,20,22),grey(.5),grey(.3))
colplot <- ifelse(gwas$CHR%in%c(2,4,6,8,10,12,14,16,18,20,22),"darkblue","lightblue")

jpeg(file=paste(file.results, ".manhattan.jpg", sep=""), quality=100,width=1600,height=700)
plot(pos.plot,-log10(gwas$P),col=colplot,ylab="-log10(p)", xlab="Chromosome",xaxt="n",pch=20,las=1,cex=1.5, main=paste( disease_name, sep=" "))
# ,ylim=c(0,7)
lab.pos <- tapply(pos.plot,gwas$CHR,mean)
axis(1, at=lab.pos, labels=1:22)
axis(2, at=seq(0,40,5), las=1)

genomewideline=-log10(5e-8)
abline(h=genomewideline, col="red", lty=5)
text(500000000, genomewideline + 0.05, labels=c("genome-wide significance"),
    pos=3, cex=1.0, col="red")

dev.off()
