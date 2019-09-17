setwd('~/Google Drive/Moorea/2brad_moorea/Host/')
library(vcfR)
library(adegenet)
library(vegan) 

gl=vcfR2genlight(read.vcfR("donresult.vcf.gz")) #output from angsd

pops=read.table("bamscl_site.txt",sep="\t") #bams file with a 2nd column describing variable of interest

pop(gl)=pops$V2
pca=glPca(gl,nf=3,parallel=F)

quartz()
plot(pca$scores[,1:2],col=transp(as.numeric(as.factor(pop(gl))),0.3),pch=19)
ordispider(pca$scores[,1:2],pop(gl),col=as.numeric(as.factor(pop(gl))),label=T)
ordiellipse(pca$scores[,1:2],pops$V2,label=T,draw="polygon",col="grey90")

# scatter(pca)
# pca
# loadingplot(pca,thres=0.0001)
# these <- loadingplot(pca, thres=0.0001)
# bad <- these$var.idx
# bad2 <- as.data.frame(bad)
# baddy <- bad2$bad
# baddy
# write.csv(baddy,"~/Desktop/strict/bads",row.names=F)

#gl@loc.names
#gl@loc.names[10000]
#bad[,1]
#bad2[37,]

#path = "~/Documents/My Data/BRAZIL/Elections/"
#out.file<-""
#file.names <- dir(path, pattern =".txt")
for(i in bad2){
badz <- gl@loc.names[i]
}
str(badz)
baddf <- as.data.frame(badz)
write.csv(baddf,"~/Desktop/depth/baddf.csv")
badz
#removed these from vcf file

gl=vcfR2genlight(read.vcfR("mrstrict_check.recode.vcf"))

pops=read.table("mrstrict_pops_year.txt",sep="\t")
#pops=read.table("mrstrict_pops_site.txt",sep="\t")

pop(gl)=pops$V2
pca=glPca(gl,nf=3,parallel=F)

quartz()
plot(pca$scores[,1:2],col=transp(as.numeric(as.factor(pop(gl))),0.3),pch=19)
ordispider(pca$scores[,1:2],pop(gl),col=as.numeric(as.factor(pop(gl))),label=T)
ordiellipse(pca$scores[,1:2],pops$V2,label=T,draw="polygon",col="grey90")

