library(ggplot2)
library(dplyr)
library(gridExtra)
library(ggpubr)

##### Figure 2a #####
x = data.frame(groups = c("Branch","Trunk","Trunk","Trunk"),
               fill=c("WT","Promoter","upCN","amp"),
               Fraction = c(0,5,7,14)/37)

x$groups = factor(x$groups,levels=c("Trunk","Branch"))
x$fill=factor(x$fill,levels=c("WT","Promoter","upCN","amp"))

pdf(file="Figure2a.pdf",width=4/2.54,height=4.5/2.54)
ggplot(x,aes(groups,Fraction))+
  geom_bar(position="stack",stat="identity",mapping=aes(groups,Fraction,fill=fill),color="gray35")+
  scale_fill_manual(values=c("white","gray96","gray75","gray30"))+
  ylab("Cases with *TERT* alteration (%)")+
  scale_y_continuous(breaks=c(0,0.5,1),labels=c(0,50,100),limits = c(0,1))+
  scale_x_discrete(labels=c("Trunk","Branch"))+
  theme_classic2()+
  theme(
    axis.title.x = element_blank(),
    plot.title = element_blank(),
    axis.title.y = ggtext::element_markdown(size=7),
    axis.text = element_text(size=5.5),
    axis.text.x = element_text(margin=margin(1,0,0,0)),
    axis.line = element_line(linewidth=0.3),
    legend.position="none"
  )
dev.off()


##### Figure 2b #####
Inv1_cnr = read.table("./extdata/f2b_Case113_Inv1.call.cnr",header=T,sep="\t")
Inv1_cns = read.table("./extdata/f2b_Case113_Inv1.call.cns",header=T,sep="\t")
Inv1_cnr$col = ifelse(Inv1_cnr$chromosome=="chr5" & grepl("TERT",Inv1_cnr$gene),"hotpink","gray")
Inv1_cnr$size = ifelse(Inv1_cnr$col=="hotpink",0.6,0.4)
P_Inv1 = ggplot(filter(Inv1_cnr,chromosome=="chr5" & col=="gray" & (1:nrow(Inv1_cnr))%%1 == 0),aes((start+end)/2/1000000,log2)) +
  geom_point(color="gray",size=0.03)+
  geom_segment(aes(x=start/1000000,xend=end/1000000,y=log2,yend=log2),color="purple2",lineend="round",size=0.65,data=filter(Inv1_cns,chromosome=="chr5"))+
  geom_point(data=filter(Inv1_cnr,chromosome=="chr5" & col=="hotpink"),aes((start+end)/2/1000000,log2),color="hotpink",size=0.3)+
  #scale_size_continuous(range=c(0,1))+
  
  #scale_color_manual(values=c("gray","hotpink"))+
  #geom_segment(aes(x=start,xend=start,y=-3,yend=2),linetype="dashed",color="darkorchid1",size=0.35,data=filter(Inv1_cns,chromosome=="chr5"))+
  #geom_vline(xintercept=21994490,linetype="dashed",color="darkorchid1",size=0.35)+
  coord_cartesian(xlim=c(1000000,1300000+250000)/1000000,expand = F,clip = 'on')+
  #annotate("text",x=52000000,y=-0.5,label="Inv1",angle=270,size=3.5)+
  scale_y_continuous(limits=c(-2,1.6),breaks=c(-1,0,1),labels=c(1,2,4))+
  theme_classic2()+
  theme(#axis.title.y=element_text(size=8,margin=unit(c(0,2,0,0),"pt")),
    axis.title.y=element_blank(),    
    axis.line = element_line(linewidth=0.3),
    axis.text.y=element_text(size=5),
    axis.ticks.x=element_blank(),
    axis.text.x=element_blank(),
    axis.title.x=element_blank(),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin=unit(c(1,4,0,5),"pt")
  )

Inv2_cnr = read.table("./extdata/f2b_Case113_Inv2.call.cnr",header=T,sep="\t")
Inv2_cns = read.table("./extdata/f2b_Case113_Inv2.call.cns",header=T,sep="\t")
Inv2_cnr$col = ifelse(Inv2_cnr$chromosome=="chr5" & grepl("TERT",Inv2_cnr$gene),"hotpink","gray")
Inv2_cnr$size = ifelse(Inv2_cnr$col=="hotpink",0.6,0.4)
P_Inv2 = ggplot(filter(Inv2_cnr,chromosome=="chr5" & col=="gray" & (1:nrow(Inv2_cnr))%%1 == 0),aes((start+end)/2/1000000,log2)) +
  geom_point(color="gray",size=0.03)+
  geom_segment(aes(x=start/1000000,xend=end/1000000,y=log2,yend=log2),color="purple2",lineend="round",size=0.65,data=filter(Inv2_cns,chromosome=="chr5"))+
  geom_point(data=filter(Inv2_cnr,chromosome=="chr5" & col=="hotpink"),aes((start+end)/2/1000000,log2),color="hotpink",size=0.3)+
  #scale_size_continuous(range=c(0,1))+
  
  #scale_color_manual(values=c("gray","hotpink"))+
  #geom_segment(aes(x=start,xend=start,y=-3,yend=2),linetype="dashed",color="darkorchid1",size=0.35,data=filter(Inv2_cns,chromosome=="chr5"))+
  #geom_vline(xintercept=21994490,linetype="dashed",color="darkorchid1",size=0.35)+
  coord_cartesian(xlim=c(1000000,1300000+250000)/1000000,expand = F,clip = 'on')+
  #annotate("text",x=52000000,y=-0.5,label="Inv2",angle=270,size=3.5)+
  scale_y_continuous(limits=c(-2,1.6),breaks=c(-1,0,1),labels=c(1,2,4))+
  theme_classic2()+
  theme(#axis.title.y=element_text(size=8,margin=unit(c(0,2,0,0),"pt")),
    axis.title.y=element_blank(),    
    axis.line = element_line(linewidth=0.3),
    axis.text.y=element_text(size=5),
    axis.ticks.x=element_blank(),
    axis.text.x=element_blank(),
    axis.title.x=element_blank(),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin=unit(c(1,4,0,5),"pt")
  )

Met1_cnr = read.table("./extdata/f2b_Case113_Met1.call.cnr",header=T,sep="\t")
Met1_cns = read.table("./extdata/f2b_Case113_Met1.call.cns",header=T,sep="\t")
Met1_cnr$col = ifelse(Met1_cnr$chromosome=="chr5" & grepl("TERT",Met1_cnr$gene),"hotpink","gray")
Met1_cnr$size = ifelse(Met1_cnr$col=="hotpink",0.6,0.4)
P_Met1 = ggplot(filter(Met1_cnr,chromosome=="chr5" & col=="gray" & (1:nrow(Met1_cnr))%%1 == 0),aes((start+end)/2/1000000,log2)) +
  geom_point(color="gray",size=0.03)+  
  geom_segment(aes(x=start/1000000,xend=end/1000000,y=log2,yend=log2),color="purple2",lineend="round",size=0.65,data=filter(Met1_cns,chromosome=="chr5"))+xlab("5p (Mbps)")+
  geom_point(data=filter(Met1_cnr,chromosome=="chr5" & col=="hotpink"),aes((start+end)/2/1000000,log2),color="hotpink",size=0.3)+
  #scale_size_continuous(range=c(0,1))+
  
  #scale_color_manual(values=c("gray","hotpink"))+
  #geom_segment(aes(x=start,xend=start,y=-3,yend=2),linetype="dashed",color="darkorchid1",size=0.35,data=filter(Met1_cns,chromosome=="chr5"))+
  #geom_vline(xintercept=21994490,linetype="dashed",color="darkorchid1",size=0.35)+
  coord_cartesian(xlim=c(1000000,1300000+250000)/1000000,expand = F,clip = 'on')+
  scale_x_continuous(breaks=c(1,1.25,1.5),labels=c(1,1.25,1.5))+
  #annotate("text",x=52000000,y=-0.5,label="Met1",angle=270,size=3.5)+
  scale_y_continuous(limits=c(-2,1.6),breaks=c(-1,0,1),labels=c(1,2,4))+
  theme_classic2()+
  theme(
    #axis.title.y=element_text(size=8,margin=unit(c(0,2,0,0),"pt")),
    axis.title.y=element_blank(),
    axis.line = element_line(linewidth=0.3),
    axis.text.y=element_text(size=5),
    axis.text.x=element_text(size=5),
    axis.title.x=element_text(size=7),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin=unit(c(1,4,2,5),"pt")
  )

pdf(file="Figure2b.pdf",width=2.8/2.54,height=4.2/2.54)
ggarrange(P_Inv1,P_Inv2,P_Met1,ncol=1,nrow=3,heights=c(1,1,1.52),align="v")
dev.off()


