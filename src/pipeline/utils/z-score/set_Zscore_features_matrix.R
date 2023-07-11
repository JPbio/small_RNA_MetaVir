## jpbio jun-2023

args <- commandArgs()
prefix = args[6]

viral <- read.delim(args[3],sep="#",header=F,row.names=NULL,quote = NULL)
viral = t(as.matrix(viral))

nonviral <- read.delim(args[4],sep="#",header=F,row.names=NULL,quote = NULL)
nonviral = t(as.matrix(nonviral))

nohit  <- read.delim(args[5],sep="#",header=F,row.names=NULL,quote = NULL)
nohit = t(as.matrix(nohit))

viral = as.data.frame(viral[-1,])
nonviral = as.data.frame(nonviral[-1,])
nohit = as.data.frame(nohit[-1,])

viral$similarity_label = "viral"
nonviral$similarity_label = "nonviral"
nohit$similarity_label = "nohit"

matrix = rbind(viral,nonviral,nohit)

matrix[,2:50]<- as.data.frame(lapply(matrix[,2:50], function(x) {
  as.numeric(as.character(x))
}))


matrix[,44:46] = matrix[,44:46] / matrix[,ncol(matrix)-1]

matrix[,ncol(matrix) - 2] = matrix[,ncol(matrix) - 2] / matrix[,ncol(matrix)-1]

#matrix[is.infinite(matrix)]=0

n<-seq(15,35,by=1)
n2<- seq(-15,-35,by=-1)

nd = c("dens15to18","dens20to22","dens25to29","ratiosi_pi","ratio_si","dens18to35","length")
#n3<-c(n,n2)
n3<-c("Contigs_ID",c(n,n2,nd),"Similarity_label")
colnames(matrix) <- n3


matrix[is.na(matrix)] <- 0

matrix <- matrix[rowSums(matrix[2:50])>0,]

matrix[, 44:49][matrix[, 44:49] == 0] <- matrix[, 44:49][matrix[, 44:49] == 0] + 0.00001

matrix[,44:49] = log2(matrix[,44:49])

matrix = matrix[,c(1,ncol(matrix),2:(ncol(matrix)-1))]

write.table(matrix, paste0(prefix, "/Zscore_and_features_matrix.tab"), sep="\t", col.names=TRUE, row.names=F, quote=F)