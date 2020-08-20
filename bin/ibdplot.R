

args = commandArgs(trailingOnly=T)

file = args[1]
highldexcludes = args[2]
num.samples = as.numeric(args[2])
library("SNPRelate")

cmd = paste("plink --allow-no-sex --noweb",
            "--bfile",file,
            "--exclude", highldexcludes, "--range --indep-pairwise 50 5 0.2",
	    "--out tmp")
system(try(cmd))


cmd = paste("plink --allow-no-sex --noweb",
            "--bfile",file,
            "--extract tmp.prune.in",
	    "--make-bed --out tmp")
system(try(cmd))

cmd = paste("plink --allow-no-sex --noweb",
            "--bfile","tmp",
            "--genome --min 0.185",
	    "--make-bed --out tmp")
system(try(cmd))

bash  = "cat <(awk '{print $1 \"\t\" $2}' tmp.genome) <(awk '{print $3 \"\t\" $4}' tmp.genome) | sort | uniq -c | sort -n >tmp.genome.percent"

write.table(bash, "bash.sh", quote=F, col.names=F, row.names=F)
cmd= "bash bash.sh"
try(system(cmd))
options(stringsAsFactors=F)
tmp = read.table("tmp.genome.percent", h=F, col.names=c("times","FID","IID"),colClasses=c("numeric","character","character"))
tmp$percent= tmp$times/num.samples
tmp = tmp[tmp$percent >=0.9,]
if(nrow(tmp)!=0){
    write.table(tmp[,c("FID","IID")], "fail-IBD-percent-qc.txt", quote=F, row.names=F, col.names=F)
    cmd = paste("plink --allow-no-sex --noweb",
                "--bfile","tmp",
                "--remove","fail-IBD-percent-qc.txt",
                "--genome --min 0.185",
                "--make-bed --out tmp")
    system(try(cmd))
}

cmd = c("cat *genome | awk '{print $1\"\t\" $2}' > tmp.keep;cat *genome | awk '{print $3\"\t\" $4}' >> tmp.keep   ")
try(system(cmd))




cmd = paste("plink --allow-no-sex --noweb",
            "--bfile", "tmp",
            "--extract tmp.prune.in", 
	          "--keep tmp.keep",
            "--make-bed --out", "ibd_pruned")
system(try(cmd))

bed.fn <- "ibd_pruned.bed"
bim.fn <- "ibd_pruned.bim"
fam.fn <- "ibd_pruned.fam"

snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "ibd_pruned.gds")

genofile <- openfn.gds("ibd_pruned.gds")

ibd <- snpgdsIBDMLE(genofile, maf=0.05, missing.rate=0.05, num.thread = 16)
save(ibd, file = "ibd_pruned.RData")

load("ibd_pruned.RData")

ibd.coeff = snpgdsIBDSelection(ibd)
fam = read.table("ibd_pruned.fam",sep="", col.names=c("FID","IID","","","",""))
ID1 = fam$FID[match(ibd.coeff$ID1,fam$IID)]
ID2 = fam$FID[match(ibd.coeff$ID2,fam$IID)]
ibd.coeff = cbind(ID1,ibd.coeff[,1],ID2,ibd.coeff[,-1])
colnames(ibd.coeff)=c("FID1","IID1","FID2","IID2","k0","k1","niter","kinship")

write.table(ibd.coeff, paste(file,".genome",sep=""), quote = F, row.names = F)
write.table(ibd.coeff[(ibd.coeff$kinship> 0.185),], paste(file,"_related.genome",sep=""), quote = F, row.names = F)
png("relatedness.png")
plot(ibd.coeff$k0, ibd.coeff$k1, xlim=c(0,1), ylim=c(0,1),xlab="P(IBD = 0)", ylab="P(IBD = 1)",main ="")
lines(c(0,1), c(1,0), col="red", lty=2)
dev.off()

q()
n
