#load SNPRelate
library("SNPRelate")

vcf.fn<-"final.vcf"

#convert vcf to genofile
snpgdsVCF2GDS(vcf.fn, "test2.gds",  method="biallelic.only")
genofile <- snpgdsOpen("test2.gds")
genofile

#text file of population names
pop_code <- scan("pops2.txt", what=character())

#define principal components
head(pop_code)
pca <- snpgdsPCA(genofile, num.thread=2, autosome.only=FALSE)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)

#plot PCA

plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")
head(tab)
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
print(cbind(sample.id, pop_code))


#final PCA with ggplot

library(ggplot2)

tab <- data.frame(sample.id = pca$sample.id,
                  pop = factor(pop_code)[match(pca$sample.id, sample.id)],
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)

p <- ggplot(tab, aes(x=EV1, y=EV2), color = 'sample.id')

p2<-p+geom_point(aes(color = sample.id), size = 4) + theme_bw() + scale_color_manual(values = c("goldenrod1","goldenrod1","skyblue4", "darkorange1","dodgerblue4", "darkorchid4", "darkorchid4", "orange", "midnightblue", "midnightblue","steelblue3", "darkorange","darkorchid", "darkgoldenrod4", "darkgoldenrod4")) 
p2 + theme(axis.text.x = element_text(size=14),
           axis.text.y = element_text(size=14))


