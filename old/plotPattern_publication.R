
library(ggplot2)
library(reshape)
library(plyr)

args <- commandArgs()
data <- read.csv(args[3],sep="#",header=T)


###
a<- data[,2:22] + data[,23:43]
pos <- melt(data[,2:22])
pos$value[pos$value < 0] <- 0 
neg <- melt(data[,23:43]) 
tot <- melt(a)


pdf(paste(args[4],".pdf",sep = ""),width=6,height=4)
tot.sem <- ddply(tot, c("variable"), summarise,mean=mean(value), sem=sd(value)/sqrt(length(value)))
tot.sem <- transform(tot.sem, lower=mean-sem, upper=mean+sem)
ggplot(tot.sem,aes(x=variable,y=mean)) + geom_bar(stat="identity",fill="grey") + geom_errorbar(aes(ymax=upper,ymin=lower),position=position_dodge(0.9), data=tot.sem) + theme_classic(base_size=20) + xlab("small RNA length (nt)") + ylab("Z-score") +
  theme(panel.border = element_blank())+ theme(axis.line = element_line(color = 'black')) +
  scale_x_discrete(breaks=c("X0","X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","X11","X12","X13","X14","X15","X16","X17","X18","X19","X20"),labels = c("","16","","18","","20","","22","","24","","26","","28","","30","","32","","34",""))
dev.off()
