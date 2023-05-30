library(ggplot2)
library(reshape)
library(scales)
library(gplots)
require(RColorBrewer)


args <- commandArgs()

my_palette <- colorRampPalette(c( "darkgray","black" ,"yellow"))(n = 1000)

#### heatmap

data3 <- read.csv(args[3],sep="#",header=T,row.names=NULL,quote = NULL)


#defining read size name 15-35 -15 - -35
mat <- t(as.matrix(data3))

n<-seq(15,35,by=1)
n2<- seq(-15,-35,by=-1)
n3<-c(n,n2)
colnames(mat) <- n3

###Defining clustering methods
d2<- as.matrix(mat)

d2[is.na(d2)] <- 0
distCor <- function(x) as.dist(1-cor(t(x),method="pearson"))
hclustAvg <- function(x) hclust(x, method="average")

a1 <- d2
a2 <- a1[rowSums(a1)>0,]


fit <- hclustAvg(distCor(a2))
clusters <- cutree(fit,h=0.3) #model EVEs and TEs

colnames(a2) <- n3
#plotting heatmap
#dev.off()
#h <- heatmap.2(a2,Rowv=TRUE,RowSideColors=as.character(clusters),Colv =NULL,distfun=distCor,hclustfun=hclustAvg, col=my_palette, symbreak=TRUE,margins = c(3,12),keysize=0.2,trace="none",scale="row",key=1)


size = nrow(d2) * 0.11

### number of row too large to generate heatmap
if(nrow(d2) < 200){
	pdf( file = paste(args[4],".profile_correlarion.pdf",sep=""),width = 5, height = size)
	heatmap.2(a2,Rowv=TRUE,Colv=NULL,RowSideColors=as.character(clusters),distfun=distCor,hclustfun=hclustAvg,margins = c(3,12), col=my_palette, symbreak=TRUE,keysize=0.2,trace="none",scale="row",key=1)
	dev.off()
}
###
###  writing clusters files
###
cl = melt(clusters)
cl$names = rownames(cl)
head(cl)


for (i in unique(clusters) ){
    ne = length(cl[cl$value==i,"names"])
    na = (cl[cl$value==i,"names"])
    print(ne)
    
    write.table(na,paste("cluster_",i,"_",ne,"_elements.ids",sep=""),quote = FALSE,col.names = F,row.names = F)
    
}

#329 X 493
