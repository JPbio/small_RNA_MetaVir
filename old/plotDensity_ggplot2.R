library("ggplot2")

args <- commandArgs()
data <- read.csv(args[3],sep="\t",header=T)
d2 <-  data[,c(2,4,5)]

head(data$SENSE)

head(d2)
m <- max(d2[,3],d2[,2]) + 10
d2 <-d2[order(d2[,1]),]

#m <- round(m,digits=-2)


pdf(paste(args[4],".pdf",sep = ""),width=7,height=4,paper='special')
d2[,3] <- d2[,3] * -1

 ggplot(d2,aes(START,y=POSITIVE),size=0.5)+ scale_y_continuous(limits = c(-m, m)) +ylab("Small RNA density") +xlab("Reference sequence") + ylim(-as.integer(m),as.integer(m)) + theme_classic(base_size=20)+ geom_density(aes(y=POSITIVE,col="V2"),adjust=5,fill="#0099CC",colour="#0099CC",size=0.1,stat="identity") + geom_density(aes(y=NEGATIVE,col="V3"),colour="#996631",fill="#996631",size=0.1,stat="identity")  + theme(panel.border = element_blank())+ theme(axis.line = element_line(color = 'black'))  



#ggplot(d2,aes(START,y=POSITIVE,color="yellow"),size=0.5) +ylab("Density") +xlab("Reference genome") + ylim(-as.integer(m),as.integer(m)) + theme_bw(base_size=40 )+ geom_density(aes(y=POSITIVE,col="V2"),adjust=5,fill="navyblue",colour="blue",size=0.1,stat="identity") + geom_density(aes(y=NEGATIVE,col="V3"),colour="red2",fill="red4",size=0.1,stat="identity") 
#ggplot(d2,aes(START,y=POSITIVE,color="yellow"),size=0.5) +ylab("Density") +xlab("Reference genome") + ylim(-as.integer(m),as.integer(m)) + theme_bw(base_size=30)+ geom_density(aes(y=POSITIVE,col="V2"),adjust=5,fill="#0080FF",colour="#0080FF",size=0.1,stat="identity") + geom_density(aes(y=NEGATIVE,col="V3"),colour="red2",fill="#FF3333",size=0.1,stat="identity") 
dev.off()
