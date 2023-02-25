library(ggplot2)
library(reshape)

args <- commandArgs()
data <- read.csv(args[3],sep="\t",header=T)

data[is.na(data)] <- 0

head(data)
m <-max( max(rowSums(data)), min(rowSums(data)) *-1) * 1.01
b <- melt(data,id.vars="size")

b[b$value < 0,3] <-b[b$value < 0,3] * -1

pdf(paste(args[4],".pdf",sep = ""),width=7,height=4)
ggplot()+ theme_classic(base_size=20)+ geom_bar(data=b,aes(size,value,fill=variable),stat="identity",width=.8)  + xlab("small RNA length (nt)") +ylab(paste("Number of small RNAs \n(",args[5],")",sep="")) + theme(panel.border = element_blank())+ theme(axis.line = element_line(color = 'black')) +
  scale_fill_manual(values = c("#108041", "#3B54A2","#F9A31B","#EA2627")) +
  scale_x_continuous(breaks=c(seq(15,35,by=1)),labels = c("","16","","18","","20","","22","","24","","26","","28","","30","","32","","34",""))

dev.off()

