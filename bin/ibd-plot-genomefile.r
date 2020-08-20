#######################################
###  2010, David Ellinghaus         ###
###  2019, Florian Uellendahl-Werth ###
#######################################

rm(list=ls())

file.genome = commandArgs()[4]
inc<-10000
library("stream")
genome<-DSD_ReadCSV(file.genome,sep="",skip=1,take=c(7,8))
wng<-T
lines_n<-0
png(file=paste(file.genome, ".IBD-plot.png", sep=""), width=960, height=960)
plot(get_points(genome,n=1,outofpoints = "warn"), xlim=c(0,1.0), ylim=c(0,1.0), xlab="ZO", ylab="Z1",  pch=20, axes=F)
while (wng) {
  points_tmp<-get_points(genome,n=inc,outofpoints = "ignore")
  points(points_tmp, xlim=c(0,1.0), ylim=c(0,1.0), pch=20)
  if(dim(points_tmp)[1]==0){wng<-F}
  lines_n<-lines_n+inc
  print(lines_n)
  flush.console()
}

axis(1, at=seq(0,1.0,0.2), tick=T)
axis(2, at=seq(0,1.0,0.2), tick=T)
dev.off()
