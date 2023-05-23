args <- commandArgs()
da<-read.table(args[3],sep="\t",header=F)
org=args[4]
org
size=as.numeric(args[5])
size
prefix=args[6]
prefix
org_name=args[7]

pdf(paste(args[6],".positions.",org,".",org_name,".pdf",sep = ""),width=18,height=15)
plot(x=0,type="l",ylim=c(0,8),xlim=c(0,size),main=paste("contigs X ",org_name,sep=""),xlab="Reference")
legend("topright", horiz=TRUE,c("1e-2<x<1e-5","1e-5<x<1e-10","x<1e-10"), lty=c(1,1,1),lwd=c(1,1,1),pch=c("-","-","-"),col=c("firebrick1","gold2","springgreen1"),cex=0.6,title="E-value:",box.lwd = 0.4,box.col = "black",bg = "white", inset=c(0.05, 0),)

abline(v=seq(0,size,by=250),col="gray99",lwd=0.1)

for (j in 1:length(da[,1])){
        field <- strsplit(as.character(da[j,2]), "_l")[[1]]
        e<-as.numeric(da[j,5])
        d=8-j*0.1
        clip(as.numeric(da[j,3]),as.numeric(da[j,4]),10,1)
        if(e < 1e-5 & e> 1e-10){
            abline(h=d,col="gold2",lwd=2)
        }else if(e <= 1e-10){
            abline(h=d,col="springgreen1",lwd=2)
        }else{
            abline(h=d,col="firebrick1",lwd=2)
        }
        clip(-80,size+100,-10,10) 
        text(as.numeric(da[j,3])+2,d+0.04,label=paste(field[1]," (",da[j,6],"nt)",sep=""),cex=0.7)
}

b <- 1:length(da[,4])
t<-cbind(b,as.character(da$V2),as.character(da$V3),as.character(da$V4),as.character(da$V5))

#pushViewport(viewport(x=0.2,y=0.95,width=0.2,height=0.2,just=c("left","top")))
#grid.table(t,gpar.coltext = gpar(cex = 0.4),gpar.rowtext = gpar(cex = 0.4),gp=gpar(fontsize=7))
write.table(t, file = paste(args[6],".legend.",org,".",org_name,".csv",sep = ""), sep = "\t", col.names = NA)


dev.off()


