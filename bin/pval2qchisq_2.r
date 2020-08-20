rm(list=ls())

file.or.pval  <- commandArgs()[4]
file.L95.U95  <- commandArgs()[5]
t <- read.table(file=file.or.pval,header=T,stringsAsFactor=F,colClasses="numeric")
q = qchisq(t$P, df=1, ncp=0, lower.tail=FALSE, log.p=FALSE)
# returns string w/o leading and trailing whitespace
#L95 = format( gsub("^\\s+", "", t$OR**(1-1.96/sqrt(q))), digits=4, scientific=FALSE)
#U95 = format( gsub("^\\s+", "", t$OR**(1+1.96/sqrt(q))), digits=4, scientific=FALSE)
L95 = format( t$OR**(1-1.96/sqrt(q)), digits=4, scientific=FALSE)
U95 = format( t$OR**(1+1.96/sqrt(q)), digits=4, scientific=FALSE)
#L95 = format( gsub("[[:blank:]]", "", t$OR**(1-1.96/sqrt(q)), fixed=TRUE), digits=4, scientific=FALSE)
#U95 = format( gsub("[[:blank:]]", "", t$OR**(1+1.96/sqrt(q)), fixed=TRUE), digits=4, scientific=FALSE)
#L95 = format( gsub("\\s+$", "", L95_tmp), digits=4, scientific=FALSE)
#U95 = format( gsub("\\s+$", "", U95_tmp), digits=4, scientific=FALSE)
####trimws(L95, "l")()
####trimws(U95, "l")()
head(L95)
head(U95)
####L95 = math.exp(math.log(oddsratio) - 1.96*se)
####U95 = math.exp(math.log(oddsratio) + 1.96*se)
results <- cbind(as.matrix(L95), as.matrix(U95))
colnames(results) <- c("L95", "U95")
write.table(results, file=file.L95.U95, quote=F, sep="\t", row.names=F, col.names=T, append=F)
