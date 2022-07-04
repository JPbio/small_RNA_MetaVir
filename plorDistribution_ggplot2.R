library("ggplot2")

args <- commandArgs()
data <- read.csv(args[3],sep="\t",header=T)
d2 <-  data[,c(1,2,3)]

head(data$SENSE)

head(d2)
m <- max(d2[,3],d2[,2]) + 10
d2 <-d2[order(d2[,1]),]

#m <- round(m,digits=-2)

colnames(d2)<- c("Size","SENSE","ANTISENSE")

pdf(paste(args[4],".pdf",sep = ""),width=7,height=4,paper='special')
d2[,3] <- d2[,3] * -1
sink("teste.log")
cat(d2[,1])
cat("\n")
cat(d2[,2])
cat("\n")
cat(d2[,3])
cat("\n")
sink()

ggplot(d2,aes(Size,color=STRAND))+ scale_y_continuous(limits = c(-m, m)) +theme_classic(base_size=24)+ xlab("Read size (nt)") +ylab("Number of reads") + geom_bar(aes(y=SENSE,col="sense"),colour="blue",fill="blue",width=.8,stat="identity") + geom_bar(aes(y=ANTISENSE,col="antisense"),colour="red1",fill="red1",width=.8,stat="identity") + theme(panel.border = element_blank())+ theme(axis.line = element_line(color = 'black'))
#ggplot(d2,aes(Size,y="Number of reads",color=STRAND))+ylab("Number of reads") +theme_classic(base_size=40)  + geom_bar(aes(y=SENSE,col="sense"),colour="blue",fill="navyblue",width=.8,stat="identity") + geom_bar(aes(y=ANTISENSE,col="antisense"),colour="red2",fill="red4",width=.8,stat="identity") 
dev.off()
