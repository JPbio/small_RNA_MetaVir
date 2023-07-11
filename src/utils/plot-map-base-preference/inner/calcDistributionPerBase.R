library(ggplot2)
library(reshape2)

args <- commandArgs()
data <- read.csv(args[3],sep="\t",header=T)

data[is.na(data)] <- 0

head(data)
m <-max( max(rowSums(data)), min(rowSums(data)) *-1) * 1.01
b <- melt(data,id.vars="X")
b1 <- subset(b,value >= 0)
b2 <- subset(b,value < 0)

pdf(paste(args[4],".pdf",sep = ""),width=7,height=4)
ggplot()+ theme_classic(base_size=24)+ geom_bar(data=b1,aes(X,value,fill=variable),stat="identity",width=.8) + geom_bar(data=b2,aes(X,value,fill=variable),width=.8,stat="identity") + scale_fill_brewer(type = "seq", palette = 1) + ylim(-m,m) + xlab("Read size (nt)") +ylab("Number of small RNA reads") + theme(panel.border = element_blank())+ theme(axis.line = element_line(color = 'black'))  
dev.off()
