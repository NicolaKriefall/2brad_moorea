####PCA from .vcf file####

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

#nice pca
quartz()
plot(pca$scores[,1:2],col=transp(as.numeric(as.factor(pop(gl))),0.3),pch=19)
ordispider(pca$scores[,1:2],pop(gl),col=as.numeric(as.factor(pop(gl))),label=T)
ordiellipse(pca$scores[,1:2],pops$V2,label=T,draw="polygon",col="grey90")

##I didn't do, but useful:
## manually identifying outliers: click on outlier points in the plot.
## adjust n=3 parameter to identify more than 3 points
#outliers=identify(pca$scores[,1:2],labels=gl@ind.names,n=3,cex=0.8)

## re-making PCA without outliers
#gl2=gl[-outliers]
#pca2=glPca(gl2,nf=3,parallel=F)
#colors=as.numeric(as.numeric(as.factor(levels(pop(gl))))) # or use your own colors for populations
#s.class(pca2$scores[],pop(gl2),col=transp(colors,0.5),cstar=1,cellipse=1,clabel=1,axesell=F,grid=F,cpoint=2)

#stats
site <- pops$V2
scores <- pca$scores
adonis(scores ~ site,method='manhattan')
#this example helped me figure out how to do this: https://stackoverflow.com/questions/43736314/error-in-g-that-non-conformable-arrays

# Call:
#   adonis(formula = scores ~ site, method = "manhattan") 
# 
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# site        5    165.86  33.173  1.2078 0.05295  0.268
# Residuals 108   2966.34  27.466         0.94705       
# Total     113   3132.21                 1.00000       

#ggplot pca
library(ggplot2)

#adding zone variable in addition to site
zones <- read.table("part3input_inds2pops_zone.txt") #2nd column is my reef zones instead of sites
zone <- zones$V2 #just selecting second column

#adding more reef info
reefs <- read.table("part3input_inds2pops_reef.txt")
reef <- reefs$V2

p1 <- ggordiplots::gg_ordiplot(ord=scores,groups=site,scaling=0.8,choices=c(1,2),conf=0.7,spiders=FALSE,pt.size=2)
#names(p1)
spid <- p1$df_spiders
#p1$df_spiders
gg <- p1$df_ord
gg$zone <- zone
gg$reef <- reef
quartz()
ggplot(gg,aes(x=x,y=y,group=Group,color=zone,shape=reef))+
  geom_point(size=2)+
  geom_segment(data=p1$df_spiders, aes(x=x, y=y, xend=cntr.x, yend=cntr.y))+
  stat_ellipse(level=0.8,linetype=reef)+
  xlab('PC1 (4.11%)')+
  ylab('PC2 (4.06%)')+
  theme_classic()+
  scale_color_manual(values=c("#ED7953FF","#8405A7FF"),labels=c("Backreef","Forereef"))+
  scale_shape_manual(values=c(21,24,25),labels=c("MNW","MSE","TNW"))
  

#  theme(text=element_text(family="Gill Sans MT"))+ #to change font
#  scale_fill_manual(values=c("coral1","cyan3"),labels=c("Inshore","Offshore"))+
#  labs(shape="Reef zone",color="Reef zone",fill="Reef zone")

#####MDS plot/CCA from IBS matrix:#####
#[from Misha's angsd_ibs_pca.R script]
setwd("~/Google Drive/Moorea/Host")
bams=read.table("bams_no8")[,1] # list of bam files
goods=c(1:length(bams))

## reading table of pairs of replicates (tab-delimited) - skip if there are no clones
#clonepairs=read.table("clonepairs.tab",sep="\t")
#repsa= clonepairs[,1]
#repsb= clonepairs[,2]
## removing "b" replicates
#goods=which(!(bams %in% repsb))

#--------------------
# loading individual to population correspondences
i2p=read.table("part3input_inds2pops_no8.txt",sep="\t") # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list or vcf file.
row.names(i2p)=i2p[,1]
i2p=i2p[goods,]
site=i2p[,2]

# setting up colors for plotting
palette(rainbow(length(unique(site))))
colors=as.numeric(as.factor(site))
colpops=as.numeric(as.factor(sort(unique(site))))

#-------------
# clustering / PCoA based on identity by state (IBS) based on single read resampling
# (for low and/or uneven coverage)

##skipping because removed clones already
#ma = as.matrix(read.table("donresult.ibsMat"))
#hc=hclust(as.dist(ma),"ave")
#plot(hc,cex=0.5)  # this shows how similar clones are

#ma=ma[goods,goods]
#dimnames(ma)=list(bams[goods],bams[goods])
#hc=hclust(as.dist(ma),"ave")
#plot(hc,cex=0.7) # without clones

# performing PCoA and CAP
conds=data.frame(cbind(site))
pp0=capscale(ma~1)
pp=capscale(ma~site,conds)

# significance of by-site divergence
adonis(ma~site,conds) 
# Call:
#   adonis(formula = ma ~ site, data = conds) 
# 
# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
# site        1    0.0615 0.061453   1.052 0.00931  0.171
# Residuals 112    6.5427 0.058417         0.99069       
# Total     113    6.6042                  1.00000    

## eigenvectors
#plot(pp0$CA$eig) 

axes2plot=c(1,2)  
quartz()
library(adegenet) # for transp()
cmd=pp0 #change to pp for CAP, pp0 for MDS
plot(cmd,choices=axes2plot,display="sites",type="n") # choices - axes to display
points(cmd,choices=axes2plot,pch=19,col=transp(colors,alpha=0.7))
#ordihull(cmd,choices= axes2plot,groups= conds$grp,draw="polygon",col=1+as.numeric(unique(as.factor(conds$grp))),label=T)
ordispider(cmd,choices= axes2plot,groups=conds$site,col="grey80")
ordiellipse(cmd,choices= axes2plot,groups= conds$site,draw="polygon",col=colpops,label=T)

## unscaled, to identify outliers
#plot(cmd$CA$u[,axes2plot],pch=19,col=colors)
#ordispider(cmd$CA$u[,axes2plot],groups=conds$site,col="grey80")
#ordiellipse(cmd$CA$u[,axes2plot],groups= conds$site,draw="polygon",col=colpops,label=T)
#identify(cmd$CA$u[,axes2plot],labels=colnames(ma),n=3,cex=0.7)

