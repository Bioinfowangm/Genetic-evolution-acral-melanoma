library(ggplot2)
library(dplyr)
library(gridExtra)
library(ggpubr)

cytoband = read.table("./extdata/f1a_cytoband_UCSC_update.txt",header=T)
cytoband = filter(cytoband,gieStain != "stalk") %>% 
  mutate(gieStain = replace(gieStain,gieStain=="gvar","gpos100"))

cytoband$gieStain = factor(cytoband$gieStain,levels=c("acen","gpos100","gpos75","gpos50","gpos25","gneg"))
chr5 = filter(cytoband,chrom=="chr5")
P_cytoband11 = ggplot(chr5,aes(xmin=chromStart/1000000,xmax=chromEnd/1000000,ymin=0,ymax=1,fill=gieStain))+geom_rect()+
  scale_fill_manual(values=c("firebrick2","black","gray25","gray50","gray75","gray95"))+
  geom_rect(aes(xmin=0/1000000,xmax=48500000/1000000,ymin=0,ymax=1),col="black",fill=NA,size=0.5)+
  coord_cartesian(xlim=c(-500000,48500000)/1000000,expand = F,clip = 'on')+
  xlab("5p (Mbps)")+
  theme_bw()+
  theme(
    panel.border = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x=element_text(size=5),
    axis.title.x = element_text(size=7),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin=unit(c(-1,1,1,0),"pt")
  )

MIS_cnr = read.table("./extdata/f1d_Case110_MIS.call.cnr",header=T,sep="\t")
MIS_cns = read.table("./extdata/f1d_Case110_MIS.call.cns",header=T,sep="\t")
MIS_cns$log2 = ifelse(MIS_cns$log2 < -1, -1,MIS_cns$log2)
MIS_cnr$col = ifelse(MIS_cnr$chromosome=="chr5" & grepl("TERT|SKP2|RICTOR",MIS_cnr$gene),"hotpink","gray")
MIS_cnr$size = ifelse(MIS_cnr$col=="hotpink",0.6,0.4)
P_MIS_11 = ggplot(filter(MIS_cnr,chromosome=="chr5" & col=="gray" & (1:nrow(MIS_cnr))%%2 == 0),aes((start+end)/2/1000000,log2)) +
  #geom_hline(yintercept=0:6,linetype="dashed")+
  geom_point(color="gray",size=0.02)+
  geom_segment(aes(x=start/1000000,xend=end/1000000,y=log2,yend=log2),color="purple2",lineend="round",size=0.65,data=filter(MIS_cns,chromosome=="chr5"))+
  geom_point(data=filter(MIS_cnr,chromosome=="chr5" & col=="hotpink"),aes((start+end)/2/1000000,log2),color="hotpink",size=0.3)+
  coord_cartesian(xlim=c(-500000,48500000)/1000000,expand = F,clip = 'on')+
  scale_y_continuous(limits=c(-2.5,5),breaks=c(-0,2,4),labels=c(2,8,32))+
  theme_bw()+
  theme(#axis.title.y=element_text(size=8,margin=unit(c(0,2,0,0),"pt")),
    panel.border = element_blank(),
    axis.line.x = element_line(size=0.5),
    axis.title.y=element_blank(),    
    axis.text.y=element_text(size=5),
    axis.ticks.x=element_blank(),
    axis.text.x=element_blank(),
    axis.title.x=element_blank(),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin=unit(c(1,1,0,0),"pt")
  )

Inv1_cnr = read.table("./extdata/f1d_Case110_Inv1.call.cnr",header=T,sep="\t")
Inv1_cns = read.table("./extdata/f1d_Case110_Inv1.call.cns",header=T,sep="\t")
Inv1_cnr$col = ifelse(Inv1_cnr$chromosome=="chr5" & grepl("TERT|SKP2|RICTOR",Inv1_cnr$gene),"hotpink","gray")
Inv1_cnr$size = ifelse(Inv1_cnr$col=="hotpink",0.6,0.4)
P_Inv1_11 = ggplot(filter(Inv1_cnr,chromosome=="chr5" & col=="gray" & (1:nrow(Inv1_cnr))%%2 == 0),aes((start+end)/2/1000000,log2)) +
  geom_point(color="gray",size=0.02)+
  geom_segment(aes(x=start/1000000,xend=end/1000000,y=log2,yend=log2),color="purple2",lineend="round",size=0.65,data=filter(Inv1_cns,chromosome=="chr5"))+
  geom_point(data=filter(Inv1_cnr,chromosome=="chr5" & col=="hotpink"),aes((start+end)/2/1000000,log2),color="hotpink",size=0.3)+
  #coord_cartesian(xlim=c(52000000,180915260)/1000000,expand = F,clip = 'on')+
  coord_cartesian(xlim=c(-500000,48500000)/1000000,expand = F,clip = 'on')+
  #annotate("text",x=52000000,y=-0.5,label="Inv1",angle=270,size=3.5)+
  scale_y_continuous(limits=c(-2.5,5),breaks=c(-0,2,4),labels=c(2,8,32))+
  theme_bw()+
  theme(#axis.title.y=element_text(size=8,margin=unit(c(0,2,0,0),"pt")),
    panel.border = element_blank(),
    axis.line.x = element_line(size=0.5),
    axis.title.y=element_blank(),    
    axis.text.y=element_text(size=5),
    axis.ticks.x=element_blank(),
    axis.text.x=element_blank(),
    axis.title.x=element_blank(),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin=unit(c(1,1,0,0),"pt")
  )

Inv2_cnr = read.table("./extdata/f1d_Case110_Inv2.call.cnr",header=T,sep="\t")
Inv2_cns = read.table("./extdata/f1d_Case110_Inv2.call.cns",header=T,sep="\t")
Inv2_cnr$col = ifelse(Inv2_cnr$chromosome=="chr5" & grepl("TERT|SKP2|RICTOR",Inv2_cnr$gene),"hotpink","gray")
Inv2_cnr$size = ifelse(Inv2_cnr$col=="hotpink",0.6,0.4)
P_Inv2_11 = ggplot(filter(Inv2_cnr,chromosome=="chr5" & col=="gray" & (1:nrow(Inv2_cnr))%%2 == 0),aes((start+end)/2/1000000,log2)) +
  geom_point(color="gray",size=0.02)+
  geom_segment(aes(x=start/1000000,xend=end/1000000,y=log2,yend=log2),color="purple2",lineend="round",size=0.65,data=filter(Inv2_cns,chromosome=="chr5"))+
  geom_point(data=filter(Inv2_cnr,chromosome=="chr5" & col=="hotpink"),aes((start+end)/2/1000000,log2),color="hotpink",size=0.3)+
  #coord_cartesian(xlim=c(52000000,180915260)/1000000,expand = F,clip = 'on')+
  coord_cartesian(xlim=c(-500000,48500000)/1000000,expand = F,clip = 'on')+
  #annotate("text",x=52000000,y=-0.5,label="Inv2",angle=270,size=3.5)+
  scale_y_continuous(limits=c(-2.5,5),breaks=c(-0,2,4),labels=c(2,8,32))+
  theme_bw()+
  theme(#axis.title.y=element_text(size=8,margin=unit(c(0,2,0,0),"pt")),
    panel.border = element_blank(),
    axis.line.x = element_line(size=0.5),
    axis.title.y=element_blank(),    
    axis.text.y=element_text(size=5),
    axis.ticks.x=element_blank(),
    axis.text.x=element_blank(),
    axis.title.x=element_blank(),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin=unit(c(1,1,0,0),"pt")
  )

pdf(file="Figure1d_5q.pdf",width=2.29/2.54,height=4.4/2.54)
ggarrange(P_MIS_11,P_Inv1_11,P_Inv2_11,P_cytoband11,align="v",
          ncol=1,nrow=4,heights=c(1.1,1,1,0.78))
dev.off()