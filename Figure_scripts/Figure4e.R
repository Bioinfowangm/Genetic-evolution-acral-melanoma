library(ggplot2)
library(dplyr)
library(gridExtra)
library(ggpubr)

breakpoints = c(18659478,18893831,20610441,21471422,21905768,24299907,24326299,29126389,29639353,29704611,29724764,29954797,31032372,31092053,31500294,31529089,36886039,40719840,43308638,43716590,45075577,45098983)

InvA_cnr = read.table("./extdata/f4_Case062_InvA.call.cnr",header=T,sep="\t")
InvA_cns = read.table("./extdata/f4_Case062_InvA.call.cns",header=T,sep="\t")
InvA_cnr$col = ifelse(InvA_cnr$chromosome=="chr22" & grepl("EP300",InvA_cnr$gene),"hotpink","gray")
P_InvA = ggplot(filter(InvA_cnr,chromosome=="chr22" & (1:nrow(InvA_cnr))%%6 == 0 ),aes((start+end)/2,log2)) +geom_point(size=0.25,color="gray")+
  geom_vline(xintercept=breakpoints,linetype="dashed",color="black",size=0.12)+
  geom_segment(aes(x=start,xend=end,y=log2,yend=log2),lineend="butt",color="purple2",size=0.65,data=filter(InvA_cns,chromosome=="chr22"))+
  geom_point(data=filter(InvA_cnr,chromosome=="chr22" & col=="hotpink"),aes((start+end)/2,log2),color="hotpink",size=0.3)+
  
  scale_x_continuous(limits=c(15600000,51500000),breaks=c(20000000,30000000,40000000,50000000),labels=c(20,30,40,50),expand=c(0,0))+
  scale_y_continuous(breaks=c(-1,1,3),labels=c(1,4,16),limits =c(-2,3))+
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

MIS_cnr = read.table("extdata/f4_Case062_MIS.call.cnr",header=T,sep="\t")
MIS_cns = read.table("extdata/f4_Case062_MIS.call.cns",header=T,sep="\t")
MIS_cnr$col = ifelse(MIS_cnr$chromosome=="chr22" & grepl("EP300",MIS_cnr$gene),"hotpink","gray")
P_MIS = ggplot(filter(MIS_cnr,chromosome=="chr22"& (1:nrow(MIS_cnr))%%6 == 0),aes((start+end)/2,log2)) +geom_point(size=0.25,color="gray")+  
  geom_vline(xintercept=breakpoints,linetype="dashed",color="black",size=0.12)+
  
  geom_segment(aes(x=start,xend=end,y=log2,yend=log2),lineend = "butt",color="purple2",size=0.65,data=filter(MIS_cns,chromosome=="chr22"))+
  geom_point(data=filter(MIS_cnr,chromosome=="chr22" & col=="hotpink"),aes((start+end)/2,log2),color="hotpink",size=0.3)+
  scale_x_continuous(limits=c(15600000,51500000),breaks=c(20000000,30000000,40000000,50000000),labels=c(20,30,40,50),expand=c(0,0))+
  ggtitle("Chr 22")+
  scale_y_continuous(breaks=c(-1,1,3),labels=c(1,4,16),limits =c(-2,3))+
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
InvB_cns = read.table("./extdata/f4_Case062_InvB.call.cns",header=T,sep="\t")
InvB_cnr$col = ifelse(InvB_cnr$chromosome=="chr22" & grepl("EP300",InvB_cnr$gene),"hotpink","gray")
P_InvB = ggplot(filter(InvB_cnr,chromosome=="chr22"& (1:nrow(InvB_cnr))%%6 == 0),aes((start+end)/2,log2)) +geom_point(size=0.25,color="gray")+  
  geom_vline(xintercept=breakpoints,linetype="dashed",color="black",size=0.12)+
  
  geom_segment(aes(x=start,xend=end,y=log2,yend=log2),lineend="butt",color="purple2",size=0.65,data=filter(InvB_cns,chromosome=="chr22"))+
  geom_point(data=filter(InvB_cnr,chromosome=="chr22" & col=="hotpink"),aes((start+end)/2,log2),color="hotpink",size=0.3)+
  scale_x_continuous(limits=c(15600000,51500000),breaks=c(20000000,30000000,40000000,50000000),labels=c(20,30,40,50),expand=c(0,0))+
  #annotate("text",x=55700000,y=6,label="Melanoma2",size=3)+
  scale_y_continuous(breaks=c(-1,1,3),labels=c(1,4,16),limits =c(-2,3))+
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
Met1A_cns = read.table("./extdata/f4_Case062_Met1A.call.cns",header=T,sep="\t")
Met1A_cnr$col = ifelse(Met1A_cnr$chromosome=="chr22" & grepl("EP300",Met1A_cnr$gene),"hotpink","gray")

P_Met1A = ggplot(filter(Met1A_cnr,chromosome=="chr22"& (1:nrow(Met1A_cnr))%%6 == 0),aes((start+end)/2,log2)) +geom_point(size=0.25,color="gray")+
  geom_segment(aes(x=start,xend=end,y=log2,yend=log2),lineend="butt",color="purple2",size=0.65,data=filter(Met1A_cns,chromosome=="chr22"))+
  geom_vline(xintercept=breakpoints,linetype="dashed",color="black",size=0.12)+
  geom_point(data=filter(Met1A_cnr,chromosome=="chr22" & col=="hotpink"),aes((start+end)/2,log2),color="hotpink",size=0.3)+
  scale_x_continuous(limits=c(15600000,51500000),breaks=c(20000000,30000000,40000000,50000000),labels=c(20,30,40,50),expand=c(0,0))+
  xlab("22q (Mbps)")+
  scale_y_continuous(breaks=c(-1,1,3),labels=c(1,4,16),limits =c(-2,3))+
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

pdf(file="Figure4e.pdf",width=4.6/2.54,height=4/2.54)
ggarrange(P_MIS,P_InvA,P_InvB,P_Met1A,ncol=1,nrow=4,heights=c(1.1,1,1,1.6),align="v")
dev.off()
