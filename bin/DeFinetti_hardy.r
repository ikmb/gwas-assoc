###############################
###  2009, David Ellinghaus ###
###############################

## input: data.frame with genotype counts
##         - columns: geno.11: counts for AA
##                    geno.12: counts for AB
##                    geno.22: counts for BB

makeDefinetti <- function(snpSummary, label)
{
  switchStrand <- rep(c(TRUE, FALSE), length = nrow(snpSummary))
  snpSummary$n <- with(snpSummary, geno.11 + geno.12 + geno.22)
  freqs <- with(snpSummary,
                data.frame(freqsA =  (geno.11 + 0.5 * geno.12) / n,
                           freqsAB = geno.12 / (geno.11 + geno.12 + geno.22)))
  freqs$freqsA[switchStrand] <- 1 - freqs$freqsA[switchStrand]

  #print(freqs)
  plot(1, type="n", col="black", xlim = c(0, 1), ylim = c(0,1),
       main=paste("De Finetti diagram", "-  ", label),
       xlab="p", ylab=expression(p[12]),
       xaxs = "i", yaxs = "i", las = 1, xaxt="n", yaxt="n")

  axis(1, seq(0, 1, by = 0.1), tick = T, labels = seq(0,1,by = 0.1))
  axis(2, seq(0, 1, by = 0.1), tick = T, labels = seq(0,1,by = 0.1))
  
  lines(c(0,0.5,1, 0),c(0,1,0, 0), type="l", col="black", ylab=expression(p[12])) # triangle
  with(freqs, points(freqs$freqsA, freqs$freqsAB,
                     cex = 0.5, pch = 20))
  text(x=0.2, y=0.55, labels=expression(p[22]), cex=1.2)
  text(x=0.8, y=0.55, labels=expression(p[11]), cex=1.2)

  # horizontal lines
  for ( i in seq(0,1,by = 0.1))
      lines(c(0,1), c(i,i), lty=3, lwd=0.5, col="black")
  
  # triangle grid 
  for ( i in seq(0,1,by = 0.1))
      lines(c(i,i+0.5-i/2), c(0,1-i), lty=1, lwd=0.5, col="black")
  for ( i in seq(0,1,by = 0.1))
      lines(c(i,i-i/2), c(0,i), lty=1, lwd=0.5, col="black")
  for ( i in seq(0.05,0.95,by = 0.1))
      lines(c(i,i+0.5-i/2), c(0,1-i), lty=3, lwd=0.5, col="black")
  for ( i in seq(0.05,0.95,by = 0.1))
      lines(c(i,i-i/2), c(0,i), lty=3, lwd=0.5, col="black")
  for ( i in seq(0.05,0.95,by = 0.1))
      lines(c(i/2,1-i/2), c(i,i), lty=3, lwd=0.5, col="black")
  
  # ticks
  y <- seq(0.1,0.9,by=.1)
  x <- y/2
  segments(x-.01, y, x, y, col="black")
  text(x = x-.03, y = y, labels = as.character(rev(seq(0.1,0.9,by=.1))), cex=1.3)

  y2 <- rev(seq(0.1,0.9,by=.1))
  x2 <- y/2 + .5
  segments(x2+.01, y2, x2, y2, col="black")
  text(x = x2+.03, y = y2, labels = as.character(seq(0.1,0.9,by=.1)), cex=1.3)
   
  ## HWE curve
  curve(2 * (x - x^2), from = 0, to = 1, col = "red", lwd = 2, add = TRUE)
}

file.hardy              <- commandArgs()[4] # hardy input file
file.out.controls       <- commandArgs()[5] # output filename controls 
file.out.cases          <- commandArgs()[6] # output filename cases 
file.out.cases.controls <- commandArgs()[7] # output filename cases and controls 

hardy <-read.table(file=file.hardy,header=T,stringsAsFactor=F)

# Controls
genos <-matrix(unlist(strsplit(subset(hardy, hardy$TEST == "UNAFF")$GENO, "/")), ncol=3, byrow=T)
genos <-data.frame(geno.11 = as.numeric(genos[,1]),
                geno.12 = as.numeric(genos[,2]),
                geno.22 = as.numeric(genos[,3]))
jpeg(file=paste(file.out.controls, ".jpg", sep=""), quality=100,width=900,height=700)
makeDefinetti(genos, label = "Controls")
dev.off()

# Cases
genos <-matrix(unlist(strsplit(subset(hardy, hardy$TEST == "AFF")$GENO, "/")), ncol=3, byrow=T)
genos <-data.frame(geno.11 = as.numeric(genos[,1]),
                geno.12 = as.numeric(genos[,2]),
                geno.22 = as.numeric(genos[,3]))
jpeg(file=paste(file.out.cases, ".jpg", sep=""), quality=100,width=900,height=700)
makeDefinetti(genos, label = "Cases")
dev.off()

# Cases and Controls
genos <-matrix(unlist(strsplit(subset(hardy, hardy$TEST == "ALL")$GENO, "/")), ncol=3, byrow=T)
genos <-data.frame(geno.11 = as.numeric(genos[,1]),
                geno.12 = as.numeric(genos[,2]),
                geno.22 = as.numeric(genos[,3]))
jpeg(file=paste(file.out.cases.controls, ".jpg", sep=""), quality=100,width=900,height=700)
makeDefinetti(genos, label = "Cases_Controls")
dev.off()
