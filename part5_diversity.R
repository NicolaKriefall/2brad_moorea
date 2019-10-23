df <- read.table("~/Google Drive/Moorea/2brad_moorea/part5_het_out.txt") #read in data
df$het <- df$V3/(df$V2+df$V3) #heterozygosity calculations
df$V1 <- sub("TO","TNWO",df$V1) #just renaming some sites
df$V1 <- sub("TI","TNWI",df$V1) #just renaming some sites
df$site <- substr(df$V1, 0, 4)
df$site <- as.factor(df$site)
str(df)
plot(het~site,data=df)


