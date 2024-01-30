library(ggplot2)
library(dplyr)
library(gridExtra)
library(ggpubr)

InvA_cnr = read.table("./extdata/f4_Case062_InvA.call.cnr",header=T,sep="\t")
InvA_cnr$log2 = ifelse(InvA_cnr$log2< -1,(InvA_cnr$log2-8)/9,InvA_cnr$log2)
InvA_cns = read.table("./extdata/f4_Case062_InvA.call.cns",header=T,sep="\t")
InvA_cns$log2 = ifelse(InvA_cns$log2< -1,(InvA_cns$log2-8)/9,InvA_cns$log2)

InvA_cnr$col = ifelse(InvA_cnr$chromosome=="chr5" & grepl("TERT",InvA_cnr$gene),"hotpink","gray")
P_InvA = ggplot(filter(InvA_cnr,chromosome=="chr5" & (1:nrow(InvA_cnr))%%1 == 0 ),aes((start/1000000+end/1000000)/2,log2)) +geom_point(size=0.25,color="gray")+
  geom_segment(aes(x=start/1000000,xend=end/1000000,y=log2,yend=log2),color="purple2",lineend="round",size=0.65,data=filter(InvA_cns,chromosome=="chr5"))+
  geom_point(data=filter(InvA_cnr,chromosome=="chr5" & col=="hotpink"),aes((start/1000000+end/1000000)/2,log2),color="hotpink",size=0.3)+
  #geom_hline(yintercept=-1.704,linetype="dashed",color="gray50",size=0.35)+
  coord_cartesian(xlim=c(1000000,1300000+250000)/1000000,expand = F,clip = 'on')+
  #annotate("text",x=55700000,y=6,label="Melanoma1",size=3)+
  scale_y_continuous(breaks=c(-1,0,1),labels=c(1,2,4),limits =c(-1.6,1.6))+
  theme_bw()+
  theme(axis.title.y=element_blank(),
        axis.line.x = element_line(size=0.4),
        axis.text.y=element_text(size=5),
        axis.ticks.x = element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        panel.grid = element_blank(),
        legend.position = "none",
        #panel.grid = element_blank(),
        plot.margin=unit(c(0,4,0,3),"pt"),
        panel.border = element_blank()
  )

MIS_cnr = read.table("./extdata/f4_Case062_MIS.call.cnr",header=T,sep="\t")
MIS_cnr$log2 = ifelse(MIS_cnr$log2< -1,(MIS_cnr$log2-8)/9,MIS_cnr$log2)
MIS_cns = read.table("./extdata/f4_Case062_MIS.call.cns",header=T,sep="\t")
MIS_cns$log2 = ifelse(MIS_cns$log2< -1,(MIS_cns$log2-8)/9,MIS_cns$log2)

MIS_cnr$col = ifelse(MIS_cnr$chromosome=="chr5" & grepl("TERT",MIS_cnr$gene),"hotpink","gray")
P_MIS = ggplot(filter(MIS_cnr,chromosome=="chr5"& (1:nrow(MIS_cnr))%%1 == 0),aes((start/1000000+end/1000000)/2,log2)) +geom_point(size=0.25,color="gray")+
  geom_segment(aes(x=start/1000000,xend=end/1000000,y=log2,yend=log2),lineend = "butt",color="purple2",size=0.65,data=filter(MIS_cns,chromosome=="chr5"))+
  geom_point(data=filter(MIS_cnr,chromosome=="chr5" & col=="hotpink"),aes((start/1000000+end/1000000)/2,log2),color="hotpink",size=0.3)+
  coord_cartesian(xlim=c(1000000,1300000+250000)/1000000,expand = F,clip = 'on')+
  #geom_hline(yintercept=-0.774,linetype="dashed",color="gray50",size=0.35)+
  scale_y_continuous(breaks=c(-1,0,1),labels=c(1,2,4),limits =c(-1.6,1.6))+
  theme_bw()+
  theme(#plot.title =element_text(size=7,margin=margin(1,0,1,0),hjust=0.5),
    plot.title=element_blank(),
    axis.title.y=element_blank(),
    axis.text.y=element_text(size=5),
    axis.ticks.x = element_blank(),
    axis.text.x=element_blank(),
    panel.grid = element_blank(),
    axis.title.x=element_blank(),
    legend.position = "none",
    #panel.grid = element_blank(),
    plot.margin=unit(c(2,4,0,3),"pt"),
    axis.line.x = element_line(size=0.4),
    panel.border = element_blank()
  )

InvB_cnr = read.table("./extdata/f4_Case062_InvB.call.cnr",header=T,sep="\t")
InvB_cnr$log2 = ifelse(InvB_cnr$log2< -1,(InvB_cnr$log2-8)/9,InvB_cnr$log2)
InvB_cns = read.table("./extdata/f4_Case062_InvB.call.cns",header=T,sep="\t")
InvB_cns$log2 = ifelse(InvB_cns$log2< -1,(InvB_cns$log2-8)/9,InvB_cns$log2)

InvB_cnr$col = ifelse(InvB_cnr$chromosome=="chr5" & grepl("TERT",InvB_cnr$gene),"hotpink","gray")
P_InvB = ggplot(filter(InvB_cnr,chromosome=="chr5"& (1:nrow(InvB_cnr))%%1 == 0),aes((start/1000000+end/1000000)/2,log2)) +geom_point(size=0.25,color="gray")+
  geom_segment(aes(x=start/1000000,xend=end/1000000,y=log2,yend=log2),lineend="butt",color="purple2",size=0.65,data=filter(InvB_cns,chromosome=="chr5"))+
  geom_point(data=filter(InvB_cnr,chromosome=="chr5" & col=="hotpink"),aes((start/1000000+end/1000000)/2,log2),color="hotpink",size=0.3)+
  coord_cartesian(xlim=c(1000000,1300000+250000)/1000000,expand = F,clip = 'on')+
  #geom_hline(yintercept=-2.745,linetype="dashed",color="gray50",size=0.35)+
  scale_y_continuous(breaks=c(-1,0,1),labels=c(1,2,4),limits =c(-1.6,1.6))+
  theme_bw()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_text(size=5),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        legend.position = "none",
        #panel.grid = element_blank(),
        plot.margin=unit(c(0,4,0,3),"pt"),
        axis.line.x = element_line(size=0.4),
        panel.border = element_blank()
  )

Met1A_cnr = read.table("./extdata/f4_Case062_Met1A.call.cnr",header=T,sep="\t")
Met1A_cnr$log2 = ifelse(Met1A_cnr$log2< -1,(Met1A_cnr$log2-8)/9,Met1A_cnr$log2)
Met1A_cns = read.table("./extdata/f4_Case062_Met1A.call.cns",header=T,sep="\t")
Met1A_cns$log2 = ifelse(Met1A_cns$log2< -1,(Met1A_cns$log2-8)/9,Met1A_cns$log2)
Met1A_cnr$col = ifelse(Met1A_cnr$chromosome=="chr5" & grepl("TERT",Met1A_cnr$gene),"hotpink","gray")

P_Met1A = ggplot(filter(Met1A_cnr,chromosome=="chr5"& (1:nrow(Met1A_cnr))%%1 == 0),aes((start/1000000+end/1000000)/2,log2)) +geom_point(size=0.25,color="gray")+
  geom_segment(aes(x=start/1000000,xend=end/1000000,y=log2,yend=log2),lineend="butt",color="purple2",size=0.65,data=filter(Met1A_cns,chromosome=="chr5"))+
  geom_point(data=filter(Met1A_cnr,chromosome=="chr5" & col=="hotpink"),aes((start/1000000+end/1000000)/2,log2),color="hotpink",size=0.3)+
  coord_cartesian(xlim=c(1000000,1300000+250000)/1000000,expand = F,clip = 'on')+
  scale_x_continuous(breaks=c(1,1.25,1.5),labels=c(1,1.25,1.5))+
  scale_y_continuous(breaks=c(-1,0,1),labels=c(1,2,4),limits =c(-1.6,1.6))+
  xlab("5p (Mbps)")+
  theme_bw()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_text(size=5),
        panel.grid = element_blank(),
        axis.text.x=element_text(size=5),
        axis.title.x = element_text(size=7,margin=margin(0,0,0,0)),
        legend.position = "none",
        #panel.grid = element_blank(),
        plot.margin=unit(c(0,4,0,3),"pt"),
        axis.line.x = element_line(size=0.4),
        panel.border = element_blank()
  )

pdf(file="Figure4f.pdf",width=2.5/2.54,height=4/2.54)
ggarrange(P_MIS,P_InvA,P_InvB,P_Met1A,ncol=1,nrow=4,heights=c(1.1,1,1,1.6),align="v")
dev.off()