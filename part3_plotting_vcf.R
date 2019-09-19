setwd('~/Google Drive/Moorea/2brad_moorea/Host/')
library(vcfR)
library(adegenet)
library(vegan) 

gl=vcfR2genlight(read.vcfR("~/Google Drive/Moorea/Host/donresult.vcf.gz")) #output from angsd

##some visualizations
# glPlot(gl, posi="topleft")
##white blocks would indicate missing data

##allele frequency spectrum
# myFreq <- glMean(gl)
# hist(myFreq, proba=TRUE, col="gold", xlab="Allele frequencies",
#      main="Distribution of (second) allele frequencies")
# temp <- density(myFreq)
# lines(temp$x, temp$y*1.8,lwd=3)

pops=read.table("part3input_inds2pops_no8.txt",sep="\t") #bams file with a 2nd column describing variable of interest
#pops=read.table("bamscl_year_no8.txt",sep="\t") #another variable I looked at

pop(gl)=pops$V2 #you can have other columns with other variables in other columns, select which one for analysis here
pca=glPca(gl,nf=3,parallel=F) #make pca

quartz()
plot(pca$scores[,1:2],col=transp(as.numeric(as.factor(pop(gl))),0.3),pch=19)
ordispider(pca$scores[,1:2],pop(gl),col=as.numeric(as.factor(pop(gl))),label=T)
ordiellipse(pca$scores[,1:2],pops$V2,label=T,draw="polygon",col="grey90")

scores <- pca$scores
adonis(scores ~ site,method='eu')
#this example helped me figure out how to do this: https://stackoverflow.com/questions/43736314/error-in-g-that-non-conformable-arrays
