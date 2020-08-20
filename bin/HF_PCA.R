library(snpStats)

# based on template_HF_PCA.CaseControl.R

# Expected args: <plink-basename> <numof_pc> <annotations> <evec> <diagnoses> <outfile>
# where <diagnoses> is a comma-separated list of diagnoses to check.
# <outfile> will be removed if existing
args <- commandArgs(TRUE)

basename <- args[1]
dataset <- read.plink(bed=paste(basename,"bed",sep="."),
                      bim=paste(basename,"bim",sep="."),
                      fam=paste(basename,"fam",sep="."))
nComp <- strtoi(args[2])
annotation <- read.csv(args[3], sep="\t", header=T)
PCA <- read.table(args[4], header=T)
diagnoses <- unlist(strsplit(args[5], ","))
outfile <- args[6]



matched <- match(as.character(PCA$IID), as.character(annotation$individualID))
# PCs start in column #3
U <- svd(cbind(1, PCA[,c(3:(2+nComp))]))$u

basic_indep_anova_test <- function(g, group, select) {
	#basic test of independence

	#choose to use real or permuted genotype
	#    gg <- sample(g, replace=F)
	gg <- g
	ii <- as.numeric(gg) >= 0
	#remove homozygous SNPs
	if(length(unique(gg[select&ii])) <= 1 ){
		result<- list(pval=1, sample=0)
	}
	else{
		y <- gg[select & ii]
		x <- group[select & ii ]
		m <- anova(lm(y~as.factor(x)))
		result <- list(pval=m$P[1], sample=m$Df[2] )
	}
	result
}

adv_indep_anova_test <- function(g, group, select, matched, U) {
	#basic test of independence

	#choose to use real or permuted genotype
	#    gg <- sample(g, replace=F)
	gg <- g[matched]
	ii <- as.numeric(gg) >= 0
	gg[!ii] <- NA
	#remove homozygous SNPs
	if(length(unique(gg[select[matched]&ii])) <= 1 ){
		result<- 1
	}
	else{
		y <- gg
		y <- y-mean(y, na.rm=T);
		sd_y <- sqrt(mean(y*y, na.rm=T)/2)
		y <- y/sd_y
		N <- length(U[,1])
		pc <- colMeans(U*y, na.rm=T)*N
		y_res <- y - rowSums(sweep(U, 2, pc, `*`))
		x <- group[matched][select[matched] & ii ]
		m <- anova(lm(y_res[select[matched] & ii ] ~ as.factor(x)))
		result <-m$P[1]
	}
	result
}

basic_remove_one_test <- function(g, group, select, batch){

	#choose to use real or permuted genotype
	#    gg <- sample(g, replace=F)
	gg <- g
	ii <- as.numeric(gg) >= 0

	result <- c()
	#remove homozygous SNPs
	if(length(unique(gg[select&ii])) <= 1 ){
		result<- rep(1,length(batch))
	}
	else{
		for (b in batch){
			if(length(unique(gg[select&ii&group!=b])) <= 1 ){
				result <- c(result, 1)
			}
			else if(sum(select&ii&group==b)==0 ) {
				result <- c(result, NA)
			}
			else{
				y <- gg[select & ii&group!=b]
				x <- group[select & ii &group!=b]
				m <- anova(lm(y~as.factor(x)))
				result <- c(result, m$P[1])
			}
		}
	}
	c(result)
}

adv_remove_one_test <- function(g, group, select,batch, matched, U){
	#advanced test of independence with adjustment for PCA

	#choose to use real or permuted genotype
	#    gg <- sample(g, replace=F)
	gg <- g[matched]
	ii <- as.numeric(gg) >= 0
	result <- c()
	gg[!ii] <- NA
	#remove homozygous SNPs
	if(length(unique(gg[select[matched]&ii])) <= 1 ){
		result<- rep(1, length(batch));
	}
	else{
		y <- gg
		y <- y-mean(y, na.rm=T);
		sd_y <- sqrt(mean(y*y, na.rm=T)/2)
		y <- y/sd_y
		N <- length(U[,1])
		pc <- colMeans(U*y, na.rm=T)*N
		y_res <- y - rowSums(sweep(U, 2, pc, `*`))

		for (b in batch){
	      batchesLeft <- group[select[matched]&ii&group[matched]!=b]
	      batchesLeft <- unique(batchesLeft[!is.na(batchesLeft)])
	      batchesLeft <- length(batchesLeft)

			if(batchesLeft <= 1 ){
				result <- c(result, 1)
			}
			else if(batchesLeft == 0) {
				result <- c(result, NA)
			}
			else{
				x <- group[matched][select[matched] & ii& group[matched]!=b  ]
				m <- anova(lm(y_res[select[matched] & ii & group[matched]!=b ] ~ as.factor(x)))
				result <- c(result, m$P[1])
			}
		}
		result
	}
}

hf_wrapper <- function(g, annotation, U) {
	g[is.na(g)] <- -1
	pval <- c()
	p_batch <- c()
	i<-1
	group <- annotation$batch
	batch <- sort(unique(group))
	eth	<- "European"

	for(diagnosis in diagnoses) {
		select <- annotation$diagnosis==diagnosis & annotation$ethnicity_predicted ==eth
		ret <- adv_indep_anova_test(g, group, select, matched, U)
		pval[i] <- ret
		ret_batch <- adv_remove_one_test(g, group, select,batch, matched, U)
		p_batch <- c(p_batch, ret_batch)
		i<-i+1
	   gc()
	}
	return(c(length(pval)+length(p_batch), c(pval, p_batch)))
}

genotypes = as(dataset$genotypes, "numeric")
out <- t(apply(genotypes, 2, hf_wrapper, annotation, U))
out <- out[,-1] # first column is vector length

# write to result file
if(file.exists(outfile)) {
   file.remove(outfile)
}

for(row in 1:nrow(out)) {
   cat(dataset$map$chromosome[row],
       dataset$map$snp.name[row],
       dataset$map$position[row],
       dataset$map$allele.1[row],
       out[row,],
       "\n",
       file=outfile, sep="\t", append=TRUE)
}
