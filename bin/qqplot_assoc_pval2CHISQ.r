rm(list=ls())
library(snpStats)

file.results    <- commandArgs()[4]
numof.cases     <- as.numeric(commandArgs()[5])
numof.controls  <- as.numeric(commandArgs()[6])
disease_name    <- commandArgs()[7]

# get chisq values
pval <- read.table(file.results, header=T)$P
# calc vector of quantiles
chivect_observed <- qchisq(pval, df=1, ncp=0, lower.tail = FALSE, log.p = FALSE)

# qqplot
retval = qq.chisq(chivect_observed, overdisp=TRUE, slope.one=FALSE, main=basename(file.results), slope.lambda=TRUE)
lambda <- as.numeric(retval[3])
#print(lambda)
#print(numof.cases)
#print(numof.controls)
lambda_1000 <- 1.0 + (lambda - 1.0) * ( (1/numof.cases + 1/numof.controls) / (1/1000 + 1/1000) )

# plot results
#pdf(file=paste(file.results, ".pdf", sep=""), width=6, height=6)
jpeg(file=paste(file.results, ".qqplot.pval2chisq.jpg", sep=""), quality=100,width=700,height=700)
retval = qq.chisq(chivect_observed, x.max=80, overdisp=TRUE, slope.one=FALSE, main=paste( disease_name, ", lambda_1000=", round(lambda_1000, digits=4), ", slope.lambda=T" ,sep="" ) , slope.lambda=TRUE)
dev.off()
