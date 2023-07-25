library(ggplot2)
library(reshape2)

args <- commandArgs()
data <- read.csv(args[3],sep="\t",header=T)

data[is.na(data)] <- 0

head(data,30)
m <-max( max(rowSums(data[,2:5])), min(rowSums(data[,2:5])) *-1) * 1.01
b <- melt(data,id.vars="X")
b1 <- subset(b,value >= 0)
b2 <- subset(b,value < 0)

pdf(paste(args[4],".pdf",sep = ""),width=7,height=4)

ggplot()+ theme_classic(base_size=20)+ geom_bar(data=b1,aes(X,value,fill=variable),stat="identity",width=.8) + geom_bar(data=b2,aes(X,value,fill=variable),width=.8,stat="identity") + ylim(-m,m) + xlab("small RNA length (nt)") +ylab("Number of small RNAs") + theme(panel.border = element_blank())+ theme(axis.line = element_line(color = 'black')) +
scale_fill_manual(values = c("#108041", "#3B54A2","#F9A31B","#EA2627")) +
scale_x_continuous(breaks=c(seq(15,35,by=1)),labels = c("","16","","18","","20","","22","","24","","26","","28","","30","","32","","34",""))

dev.off()
