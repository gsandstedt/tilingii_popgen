###convert multisample VCF into fasta file using vcf-kit 

#install and load packages
install.packages("ape")
install.packages("phangorn")
install.packages("seqinr")

library(ape)
library(phangorn)
library(seqinr)

#read in fasta file and convert to phyDat
tilingii <- read.dna("final.fa", format="fasta")
tilingii_phyDat <- phyDat(tilingii, type = "AA", levels = NULL)
print(tilingii_phyDat)

###final tree for densitree
#load packages for densitree
library(phytools)

#bootstrap
bs <- bootstrap.phyDat(tilingii_phyDat, FUN =function(x) nj(dist.ml(x)), bs=1000)

#root by Dentilobus and drop root tip
bs1 <-root(bs, "DENT")
ntrees<-lapply(bs1,drop.tip,tip="DENT")

#smooth tree
ntrees1<-lapply(ntrees,chronopl, 1, node = "root", CV = FALSE )

#convert to multiphylo
class(ntrees1)<-"multiPhylo"
print(ntrees1)

###write tree to densitree software

write.tree(ntrees1[1:1000], file="final_nj_tree.tre")
