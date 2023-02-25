
library(ggplot2)
library(reshape)
library(plyr)

args <- commandArgs()
data <- read.csv(args[3],sep="#",header=TRUE)


###
pdf(paste(args[4],".both_strands.pdf",sep = ""),width=6,height=4)


a<- data[,2:22] + data[,23:43]
pos <- melt(data[,2:22])
pos$value[pos$value < 0] <- 0 
neg2 <- melt(data[,23:43] * -1)

neg2$value[neg2$value > 0] <- 0 

pos.sem <- ddply(pos, c( "variable"), summarise,mean=mean(value), sem=sd(value)/sqrt(length(value)))
pos.sem <- transform(pos.sem, lower=mean-sem, upper=mean+sem)
y = function(x) x * -1
neg.sem <- ddply(neg2, c( "variable"), summarise,mean=mean(value), sem=sd(value)/sqrt(length(value)))
neg.sem <- transform(neg.sem, lower=mean-sem, upper=mean+sem)
pos.sem$variable <- neg.sem$variable

print(pos.sem)
print(neg.sem)

ymin <- min(neg.sem$lower)
ymax <-  max(pos.sem$upper) 
m<-max(ymin *-1,ymax)

ggplot(neg.sem,aes(x=variable,y=mean)) + geom_bar(stat="identity",fill="grey") + geom_errorbar(aes(ymax=upper,ymin=lower),position=position_dodge(0.9), data=neg.sem) +  
  geom_bar(data=pos.sem,aes(x=variable,y=mean),stat="identity",fill="grey") + geom_errorbar(aes(ymax=upper,ymin=lower),position=position_dodge(0.9), data=pos.sem) +
  theme_classic(base_size=20) + xlab("small RNA length (nt)") + ylab("Z-score") + ylim(-m,m)+
  theme(panel.border = element_blank())+ theme(axis.line = element_line(color = 'black')) +
  scale_x_discrete(breaks=c("X21","X22","X23","X24","X25","X26","X27","X28","X29","X30","X31","X32","X33","X34","X35","X36","X37","X38","X39","X40","X41"),labels = c("","16","","18","","20","","22","","24","","26","","28","","30","","32","","34",""))
dev.off()
