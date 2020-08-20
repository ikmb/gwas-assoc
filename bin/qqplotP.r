rm(list=ls())

file.results    <- commandArgs()[4]
numof.cases     <- as.numeric(commandArgs()[5])
numof.controls  <- as.numeric(commandArgs()[6])
title           <- commandArgs()[7]
qqman           <- commandArgs()[8]
gwas<-read.table(file=file.results,header=T,stringsAsFactor=F)

source(qqman)

# calc vector of quantiles
# version 1, lambda
q <- qchisq(gwas$P, df=1, ncp=0, lower.tail = FALSE, log.p = FALSE)
lambda = median(q)/0.456
print(lambda)

# version 2, lambda
#library(snpMatrix)
#retval = qq.chisq(q, overdisp=TRUE, slope.one=FALSE, main=basename(file.results), slope.lambda=TRUE)
#lambda <- as.numeric(retval[3])
#print(lambda)
lambda_1000 <- 1.0 + (lambda - 1.0) * ( (1/numof.cases + 1/numof.controls) / (1/1000 + 1/1000) )
print(lambda_1000)

## version 3, lambda
#library('GenABEL')
#res <- estlambda(gwas$P, method="median", plot=TRUE, df=1, main="", sub="")
#print(format(res$estimate, digits=3, scientific=FALSE))

# plot QQ
jpeg(file=paste(file.results, ".qqplot.pval.jpg", sep=""), quality=100,width=700,height=700)
qq(gwas$P, title=title, lambda=lambda, lambda_1000=lambda_1000, numof.cases=numof.cases, numof.controls=numof.controls)
dev.off()
