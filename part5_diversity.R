df <- read.table("df")
df$het <- df$V3/(df$V2+df$V3)
df$site <- substr(df$V1, 0, 4)
df[86:137,5] <- substr(df[86:137,5], 0, 2)
str(df)
nocrap <- df[2:137,]
plot(het~site,data=nocrap)
hetz <- nocrap$het
sitez <- nocrap$site
dfmaybe <- data.frame(sitez,hetz)
plot(hetz~sitez,data=dfmaybe)
a1 <- aov(hetz~sitez,data=dfmaybe)
summary(a1)
