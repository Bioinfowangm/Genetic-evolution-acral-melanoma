library(ggplot2)
library(ggtext)
library(dplyr)
library(reshape2)
library(gridExtra)
library(ggpubr)
library(gdata)

calc_distance = function(x){
  mutDis = c()
  for (i in 1:nrow(x)){
    dis = NA
    if(i==1){
      if(x$Chr[2]==x$Chr[1]) dis = c(dis,x$Start[2]-x$Start[1])
    }
    else if(i==nrow(x)){
      if(x$Chr[i]==x$Chr[i-1]) dis = c(dis,x$Start[i]-x$Start[i-1])
    }
    else{
      if ((x$Chr[i-1] == x$Chr[i])){
        dis = c(dis,x$Start[i]-x$Start[i-1])
      }
      if( (x$Chr[i+1] == x$Chr[i])){
        dis = c(dis,x$Start[i+1]-x$Start[i])
      }
    }
    mutDis = c(mutDis,min(dis,na.rm=T))
  }
  mutDis = ifelse(mutDis == Inf,100000000,mutDis)
  x$mutDis = mutDis
  return(x)
}

chrlen = read.table("./extdata/f1a_chrLength_hg19.txt",header=T)
InvA_cnr = read.table("./extdata/f4_Case062_InvA.call.cnr",header=T,sep="\t")
InvA_cnr$aStart = chrlen[match(InvA_cnr$chromosome,chrlen$Chr),]$aStart + InvA_cnr$start -1
InvA_cnr$aEnd = chrlen[match(InvA_cnr$chromosome,chrlen$Chr),]$aStart + InvA_cnr$end -1
InvA_cns = read.table("./extdata/f4_Case062_InvA.call.cns",header=T,sep="\t")
InvA_cns$aStart = chrlen[match(InvA_cns$chromosome,chrlen$Chr),]$aStart + InvA_cns$start -1
InvA_cns$aEnd = chrlen[match(InvA_cns$chromosome,chrlen$Chr),]$aStart + InvA_cns$end -1
InvA_cnr$log2 = ifelse(InvA_cnr$log2< -1,(InvA_cnr$log2-8)/9,InvA_cnr$log2)
InvA_cns$log2 = ifelse(InvA_cns$log2< -1,(InvA_cns$log2-8)/9,InvA_cns$log2)
InvA_cns = filter(InvA_cns,(end-start)>10000)

index = c(sample(1:round(nrow(InvA_cnr)),round(nrow(InvA_cnr)/20)),(1:nrow(InvA_cnr))[grepl("CDKN2A|CDKN2B|YAP1|RAF1",InvA_cnr[1:nrow(InvA_cnr),]$gene)])


P_InvA_cns = ggplot(InvA_cnr[index,])+
  geom_vline(data=chrlen[1:22,],aes(xintercept=aEnd),linetype="dashed",color="darkgrey",size=0.15)+
  geom_point(aes((aStart+aEnd)/2,log2),size=0.01,color="gray")+
  geom_hline(data=chrlen[1:22,],aes(yintercept= -0.008695949),color="gray25",size=0.3)+
  geom_segment(aes(x=aStart,xend=aEnd,y=log2,yend=log2),color="purple2",size=0.55,lineend="round",data=InvA_cns)+
  scale_y_continuous(limits=c(-2,6.5),breaks=c(-2,0,3,6),labels=c(0,2,16,128))+
  scale_x_continuous(limits=c(0,3036290862),expand = c(0.01,0.01),breaks=chrlen$midPoint[c(1:18,20,22,23)],labels=c(1:18,20,22,"X"))+
  ylab("CN")+  #annotate("text",x=52000000,y=6,label="PRI",angle=270,size=3.5)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(size=unit(5,"pt"),margin=unit(c(-1,0,0,0),"pt")),
        #axis.text.x = element_text(size=unit(6,"pt"),margin=unit(c(0,0,0,0),"pt")),
        axis.text.y = element_text(size=unit(5,"pt")),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin=unit(c(2,2,0,2),"pt")
  )
##################

InvB_cnr = read.table("./extdata/f4_Case062_InvB.call.cnr",header=T,sep="\t")
InvB_cnr$aStart = chrlen[match(InvB_cnr$chromosome,chrlen$Chr),]$aStart + InvB_cnr$start -1
InvB_cnr$aEnd = chrlen[match(InvB_cnr$chromosome,chrlen$Chr),]$aStart + InvB_cnr$end -1
InvB_cns = read.table("./extdata/f4_Case062_InvB.call.cns",header=T,sep="\t")
InvB_cns$aStart = chrlen[match(InvB_cns$chromosome,chrlen$Chr),]$aStart + InvB_cns$start -1
InvB_cns$aEnd = chrlen[match(InvB_cns$chromosome,chrlen$Chr),]$aStart + InvB_cns$end -1
InvB_cns$log2 = ifelse(InvB_cns$log2< -3,-9.5,InvB_cns$log2)
InvB_cnr$log2 = ifelse(InvB_cnr$log2< -1,(InvB_cnr$log2-8)/9,InvB_cnr$log2)
InvB_cns$log2 = ifelse(InvB_cns$log2< -1,(InvB_cns$log2-8)/9,InvB_cns$log2)
InvB_cns = filter(InvB_cns,(end-start)>10000)

index = c(sample(1:round(nrow(InvB_cnr)),round(nrow(InvB_cnr)/20)),(1:nrow(InvB_cnr))[grepl("CDKN2A|CDKN2B|YAP1|RAF1",InvB_cnr[1:nrow(InvB_cnr),]$gene)])

P_InvB_cns = ggplot(InvB_cnr[index,])+
  geom_vline(data=chrlen[1:22,],aes(xintercept=aEnd),linetype="dashed",color="darkgrey",size=0.15)+
  geom_point(aes((aStart+aEnd)/2,log2),size=0.01,color="gray")+
  geom_hline(data=chrlen[1:22,],aes(yintercept= 0),color="gray25",size=0.3)+
  geom_segment(aes(x=aStart,xend=aEnd,y=log2,yend=log2),color="purple2",size=0.55,lineend="round",data=InvB_cns)+
  scale_y_continuous(limits=c(-2,6.5),breaks=c(-2,0,3,6),labels=c(0,2,16,128))+
  ylab("CN")+
  scale_x_continuous(limits=c(0,3036290862),expand = c(0.01,0.01),breaks=chrlen$midPoint[c(1:18,20,22,23)],labels=c(1:18,20,22,"X"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(size=unit(5,"pt"),margin=unit(c(-1,0,0,0),"pt")),
        #axis.text.x = element_text(size=unit(6,"pt"),margin=unit(c(0,0,0,0),"pt")),
        axis.text.y = element_text(size=unit(5,"pt")),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin=unit(c(0,2,0,2),"pt")
  )
##################
chrlen$col = c(rep(c("white","black"),12))
myplotdat = filter(chrlen,Chr!="chrY")
myrect3 = ggplot(myplotdat) + 
  geom_rect(aes(xmin=aStart,xmax=aEnd,ymin=1-0.5,ymax=1+0.5,fill=col))+
  scale_fill_manual(values=c("black","gray60"))+
  scale_x_continuous(expand = c(0.01,0.01))+
  annotate("text",(myplotdat$aStart+myplotdat$aEnd)/2,y=c(rep(1,17),1.25,0.75,1.25,0.75,1.25,1),label=c(1:22,"X"),color="white",size=1.6)+
  theme(axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        plot.margin=unit(c(-2,1,-2,1),"pt"))

pdf(file="Figure4c.pdf",width=7.2/2.54,height=3.3/2.54)
ggarrange(P_InvA_cns,P_InvB_cns,myrect3,ncol=1,nrow=3,heights=c(1.06,1,0.17),align="v")
dev.off()