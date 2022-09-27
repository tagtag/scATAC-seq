#suppose GSE167050_RAW.tar was downloaded from 
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167050
#and untared in the current directory.

files <- list.files("./",pattern="mtx.gz")
files1 <- list.files("./",pattern="features")
require(Matrix) #impot sparse matrix enviroment

library(GenomicFeatures) #import grange tools

#binning scATAC-seq within 200 bp length blocks
L<-200 
j<-5 #repeat the following with changing from j=5 to j=12
#j=5  : CTX1
#j=6  : MGE1
#j=7  : CGE1
#j=8  : LGE1
#j=9  : CTX2
#j=10 : MGE2
#j=11 : CGE2
#j=12 : LGE2

x <- readMM(files[j])
y <- read.csv(files1[j],header=F,sep="\t")
gr <- GRanges(seqnames=y[,1],ranges=IRanges(start=y[,2], end=y[,3]),score=0)
seqlengths(gr) <- getChromInfoFromUCSC("mm10")[match(names(seqlengths(gr)),getChromInfoFromUCSC("mm10")[,1]),2]
gr <- trim(gr)
gr.windows <- tileGenome(seqinfo(gr), tilewidth=L,cut.last.tile.in.chrom=TRUE)
x_all<-NULL
for (i in c(1:dim(x)[2]))
{
    cat(i, " ")
gr$score<-x[,i]
score <-   coverage(gr, weight="score")
gr0.data.binnedAvg <- binnedAverage(gr.windows, score, "binned_score")
M1 <- sparseMatrix(i=(1:length(gr0.data.binnedAvg$binned_score))[gr0.data.binnedAvg$binned_score!=0],j=rep(1,sum(gr0.data.binnedAvg$binned_score!=0)),x=gr0.data.binnedAvg$binned_score[gr0.data.binnedAvg$binned_score!=0],dims=c(length(gr.windows),1))
if(i==1){
x_all <- M1
} else{
x_all <- cbind(x_all, M1)
}
}

# appply SVD to individual samples of eight samples
require(irlba)
A <- t(t(x_all)/colSums(x_all))
SVD <- irlba(A,10)
a <- which(SVD$u!=0,arr.ind=T);SVDp <- spMatrix(dim(SVD$u)[1],dim(SVD$u)[2],i=a[,1],j=a[,2],x=SVD$u[a])  #for j>5, SVDp should be renamed as SVDpn where n=j-5. e.g, for j=6, SVDp should be renamed SVDp1

#apply SVD to convinded matrix
SVD_all <- cbind(SVDp,SVD1p,SVD2p,SVD3p,SVD4p,SVD5p,SVD6p,SVD7p)
SVD_all_all <- irlba(SVD_all,10)
a <- which(SVD_all_all$u!=0,arr.ind=T)
SVD_all_all_p <- spMatrix(dim(SVD_all_all$u)[1],dim(SVD_all_all$u)[2],i=a[,1],j=a[,2],x=SVD_all_all$u[a])

#P-value computation and regions (i) selection
P <- pchisq(scale(SVD_all_all_p[SVD_all_all_p[,2]!=0,2])^2,1,lower.tail=F)
table(p.adjust(P,"BH")<0.01)

#FALSE    TRUE 
#1558701   16469 

#annotation using annotatr
require(annotatr)
sel_ranges <- gr.windows[SVD_all_all_p[,2]!=0][p.adjust(P,"BH")<0.01]
annots <- builtin_annotations()[grep("mm10",builtin_annotations())]
annotations = build_annotations(genome = 'mm10', annotations = annots)

dm_annotated = annotate_regions(
    regions = sel_ranges,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE)

#list gene symbols in selcted regions
data.frame(unique(dm_annotated@elementMetadata@listData$annot@elementMetadata@listData$symbol))

#generate pie chart of regions annotations
TABLE <- table(dm_annotated@elementMetadata@listData$annot@elementMetadata@listData$type)
names(TABLE) <- gsub("mm10_","",names(TABLE))
TABLE <- sort(TABLE) #this is a table that corresponds to pie chart

pdf(file="pie.pdf")
pie(TABLE,init.angle=0,clockwise=T)
dev.off()

#Umap 
X <- t(x_all) %*% SVD_all_all_p #repeat this line from j=5 to j=12 with renaming X to Xn, where n=j-5
X <- rbind(X,X1,X2,X3,X4,X5,X6,X7) 

files2 <- list.files("./",pattern="barcodes.tsv")
num <- NULL
for (i in c(5:length(files2)))
{
    y1 <- read.csv(files2[i],header=F)   
    num <- c(num,dim(y1)[1])
}

col<-NULL
for (i in c(1:length(num)))
{
 col <- c(col,rep(i,num[i]))    
}

set.seed(0)
umap.defaults$n_neighbors <-30
UMAP <- umap(X,config=umap.defaults)


title=c("CTX1","MGE1","CGE1","LGE1","CTX2","MGE2","CGE2","LGE2")
pdf(file="umap.pdf",width=20,height=10)
par(mfrow=c(2,4))
for (i in c(1:8))
{
    plot(UMAP$layout[col==i,],pch=16,cex=0.5,cex.axis=2,cex.lab=1.5,cex.main=2,main=title[i],xlab="UMAP1",ylab="UMAP2")
}
par(mfrow=c(1,1))
dev.off()

#Compute spearman correlations of NIJ between samples
xd <- seq(range(UMAP$layout[,1])[1],range(UMAP$layout[,1])[2],length=11)
yd <- seq(range(UMAP$layout[,2])[1],range(UMAP$layout[,2])[2],length=11)
#repeat the follwing from j=5 to j=12 with renaming cutdx1 and cutdy1 to cutdxn and cutdyn (n=j-4).
cutdx1 <- cut(UMAP$layout[col==1,1],xd)
cutdy1 <- cut(UMAP$layout[col==1,2],xd)
TABLE1 <- table(cutdx1,cutdy1)

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y,method="spearman")
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor)
}

pdf(file="pair.pdf")
pairs(cbind(as.vector(TABLE1),as.vector(TABLE5),as.vector(TABLE2),as.vector(TABLE6),as.vector(TABLE3),as.vector(TABLE7),as.vector(TABLE4),as.vector(TABLE8)),pch=16,labels=title[c(1,5,2,6,3,7,4,8)],lower.panel=panel.cor)
dev.off()

#hierarchical clisutering 
hc <- hclust(-as.dist(cor(cbind(as.vector(TABLE1),as.vector(TABLE5),as.vector(TABLE2),as.vector(TABLE6),as.vector(TABLE3),as.vector(TABLE7),as.vector(TABLE4),as.vector(TABLE8)),method="spearman")),method="av")
pdf(file="hc.pdf")
plot(hc,labels=title[c(1,5,2,6,3,7,4,8)],xlab="")
dev.off()




