setwd("~/Google Drive/Moorea/Host/")
bams=read.table("bams")[,1] # list of bam files
goods=c(1:length(bams))

#--------------------
# loading individual to population correspondences
i2p=read.table("inds2pops_example",sep="\t") # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list.
# I made mine manually in a text file, but I'm sure there's an easier way
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
ma = as.matrix(read.table("part3_allresult_example.ibsMat"))
ma=ma[goods,goods]
dimnames(ma)=list(bams[goods],bams[goods])
hc=hclust(as.dist(ma),"ave")
quartz()
plot(hc,cex=0.5)  # this shows how similar clones are
