# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Sat Jul 8 00:56:05 EDT 2017

################################################################
#   Differential expression analysis with limma
if (!require(Biobase)) {
    source("http://bioconductor.org/biocLite.R")
    biocLite("Biobase")
    library(Biobase)
}
if (!require(GEOquery)) {
    source("http://bioconductor.org/biocLite.R")
    biocLite("GEOquery")
    library(GEOquery)
}
if (!require(limma)) {
    source("http://bioconductor.org/biocLite.R")
    biocLite("limma")
    library(limma)
}
if (!require(apcluster)) {
    source("http://bioconductor.org/biocLite.R")
    biocLite("apcluster")
    library(apcluster)
}
if (!require(affycoretools)) {
    source("http://bioconductor.org/biocLite.R")
    biocLite("affycoretools")
    library(affycoretools)
}
# load series and platform data from GEO

gset <- getGEO("GSE63742", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL1355", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- "000XXX111XXX222XXX333XXX444XXX555XXX666XXX777XXX888XXX999"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G9-G0, G1-G0, G2-G1, G3-G2, G4-G3, G5-G4, G6-G5, G7-G6, G8-G7, G9-G8, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=10000)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","F","Gene.symbol","Gene.title"))
write.table(tT, file=stdout(), row.names=F, sep="\t")


# select 10000 genes and normalize.
selgeneid <- as.character(tT$ID[1:10000])
ex10000 <- ex[match(selgeneid,rownames(ex)),]

# plot to choose the number of genes
qplot(1:30,ex10000[1100,],geom = c("smooth","point"))

ex1500 <- ex10000[1:1500,]
exnorm <- apply(ex1500,1,scale)  # matrix has been transversed
heatmap(exnorm[,1:1000],Rowv = NA)
cov1500 <- (t(exnorm) %*% exnorm)/13
labels <- rep(c("0h","2h","6h","12h","24h","30","36h","72h","120h","168h"),each=3)

############## hierarchical clustering ################
dist = as.dist(1 - cov1500)
hc = hclust(dist,"ave")
plot(hc)
cutree(hc,2)
cut <- cutree(hc,5)
hc5 <- data.frame(cut,names(cut))

############## k-means ################
km = kmeans(cov1500,5,nstart = 10)
km5 <- km$cluster
km5 <- data.frame(km5,cut=names(km5),symbol=as.character(tT$Gene.symbol[1:1500]))
heatmap(exnorm[,order(km5$km5)],Rowv=NA,Colv=NA,labRow=labels)

apclus <-apcluster(cov1000,details = TRUE)
show(apclus)
heatmap(
    apclus,cov1000,margins = c(10,10),Rowv = FALSE, Colv = FALSE
)

group <- factor(hc5$cut)
plotPCA(exnorm[,1:1000],groups=group,groupnames=levels(group))
group2 <- factor(km5$km5)
plotPCA(exnorm,groups=group2,groupnames=levels(group2))

################################################################

#common genes#
selgeneid <- as.character(tT2$ID[match(commongene[1:1000],tT2$Gene.symbol)])
excommon <- ex[match(selgeneid,rownames(ex)),]
excnorm <- apply(excommon,1,scale)  # matrix has been transversed
heatmap(excnorm[,1:1000],Rowv = NA)
covcommon <- (t(exnorm[,1:1000]) %*% exnorm[,1:1000])/13

#hierarchical clustering#
dist = as.dist(1 - cov1000)
hc = hclust(dist,"ave")
plot(hc)
cutree(hc,2)
cut <- cutree(hc,5)
hc5 <- data.frame(cut,names(cut))

############## k-means ################
km = kmeans(covcommon,5)
km5 <- km$cluster
km5 <- data.frame(km5,cut=names(km5))
heatmap(exnorm[,order(commonframe$kmc5)],Rowv=NA,Colv=NA)

apclus <-apcluster(cov1000,details = TRUE)
show(apclus)
heatmap(
    apclus,corr_mat_sample,margins = c(10,10),Rowv = FALSE, Colv = FALSE
)

group <- factor(hc5$cut)
plotPCA(exnorm[,1:1000],groups=group,groupnames=levels(group))
group2 <- factor(commonframe$kmc5)
plotPCA(excnorm[,1:1000],groups=group2,groupnames=levels(group2))


################################################################
#   Boxplot for selected GEO samples
library(Biobase)
library(GEOquery)

# load series and platform data from GEO

gset <- getGEO("GSE63742", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL1355", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# group names for all samples in a series
gsms <- "000XXX111XXX222XXX333XXX444XXX555XXX666XXX777XXX888XXX999"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }
sml <- paste("G", sml, sep="")  #set group names

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# order samples by group
ex <- exprs(gset)[ , order(sml)]
sml <- sml[order(sml)]
fl <- as.factor(sml)
labels <- c("0h","2h","6h","12h","24h","30","36h","72h","120h","168h")

# set parameters and draw the plot
palette(c("#dfeaf4","#f4dfdf","#f2cb98","#dcdaa5","#dff4e4","#f4dff4","#c7ff9d","#f3f388","#ffc56e","#9dffd7", "#AABBCC"))
dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE63742", '/', annotation(gset), " selected samples", sep ='')
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=fl)
legend("topleft", labels, fill=palette(), bty="n")