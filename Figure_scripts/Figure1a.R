library(ggplot2)
library(ggtext)
library(dplyr)
library(reshape2)
library(gridExtra)
library(ggpubr)
library(gdata)
library(readxl)

##### Input and function #####
# chromosome lengths
chrlen = read.table("./extdata/f1a_chrLength_hg19.txt",header=T)

# CNV Breakpoints
Inv1_BP = read.table("./extdata/f1a_Case101_Invasive1_BP.txt")
colnames(Inv1_BP) = c("Sample","Chr","Start","Chr2","ArmStart","ArmEnd")
Inv1_BP$aStart = chrlen[match(Inv1_BP$Chr,chrlen$Chr),]$aStart + Inv1_BP$Start -1

# copy number data
Inv1_cnr = read.table("./extdata/f1a_Case101_Invasive1.call.cnr",header=T,sep="\t")
Inv1_cnr$aStart = chrlen[match(Inv1_cnr$chromosome,chrlen$Chr),]$aStart + Inv1_cnr$start -1
Inv1_cnr$aEnd = chrlen[match(Inv1_cnr$chromosome,chrlen$Chr),]$aStart + Inv1_cnr$end -1
Inv1_cns = read.table("./extdata/f1a_Case101_Invasive1.call.cns",header=T,sep="\t")
Inv1_cns$aStart = chrlen[match(Inv1_cns$chromosome,chrlen$Chr),]$aStart + Inv1_cns$start -1
Inv1_cns$aEnd = chrlen[match(Inv1_cns$chromosome,chrlen$Chr),]$aStart + Inv1_cns$end -1
Inv1_cnr$log2 = ifelse(Inv1_cnr$log2< -1,(Inv1_cnr$log2-8)/9,Inv1_cnr$log2)
Inv1_cns$log2 = ifelse(Inv1_cns$log2< -1,(Inv1_cnr$log2-8)/9,Inv1_cns$log2)
Inv1_cnr$col = ifelse(grepl("CDK4|MDM2|TERT",Inv1_cnr$gene),"hotpink","gray")
Inv1_cnr$size = ifelse(Inv1_cnr$col=="hotpink",0.6,0.4)

# cytobands
cytoband = read.table("./extdata/f1a_cytoband_UCSC_update.txt",header=T)
cytoband = filter(cytoband,gieStain != "stalk") %>% 
  mutate(gieStain = replace(gieStain,gieStain=="gvar","gpos100"))
cytoband$gieStain = factor(cytoband$gieStain,levels=c("acen","gpos100","gpos75","gpos50","gpos25","gneg"))

##
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

# somatic mutations
Inv1_mut = read_excel("./extdata/f1a_Case101_Invasive1_Mutation.xlsx",sheet=1)
Inv1_mut = calc_distance(Inv1_mut)
Inv1_mut$aStart = chrlen[match(Inv1_mut$Chr,chrlen$Chr),]$aStart + Inv1_mut$Start -1

##### Part1 #####
P_Inv1_cns = ggplot(Inv1_cnr[sample(1:round(nrow(Inv1_cnr)),round(nrow(Inv1_cnr)/12)),])+
  geom_vline(data=chrlen[1:22,],aes(xintercept=aEnd),linetype="dashed",color="darkgrey",size=0.3)+
  geom_point(aes((aStart+aEnd)/2,log2),size=0.05,color="gray")+
  geom_hline(data=chrlen[1:22,],aes(yintercept=0.000411253),color="gray25",size=0.3)+
  geom_segment(aes(x=aStart,xend=aEnd,y=log2,yend=log2),color="purple2",size=0.55,lineend="round",data=Inv1_cns)+
  scale_y_continuous(limits=c(-1.9,6),breaks=c(0,3,6),labels=c(2,16,128))+
  scale_x_continuous(limits=c(0,3036290862),expand = c(0.01,0.01),breaks=chrlen$midPoint[c(1:18,20,22,23)],labels=c(1:18,20,22,"X"))+
  ylab("Copy<br>number")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=unit(5,"pt"),margin=unit(c(0,1,0,-1),"pt")),
        axis.title.x = element_blank(),
        axis.title.y = element_markdown(size=unit(6,"pt"),margin=unit(c(0,0,0,0),"pt")),
        axis.ticks.x = element_blank(),
        plot.margin=unit(c(1,2,0,2),"pt")
  )

Inv1_BP = calc_distance(Inv1_BP)
P_Inv1_BP = ggplot(Inv1_BP)+ 
  geom_vline(data=chrlen[1:22,],aes(xintercept=aEnd),linetype="dashed",color="darkgrey",size=0.3)+
  geom_point(aes(aStart,log10(mutDis)),color="black",size=0.1) +
  ylab("Inter-CNT<br>dist. (bps)")+
  annotate("text",x=-110000000,y=4.1,label="10^4",parse=T,size=1.8,color="gray25")+
  annotate("text",x=-110000000,y=6.1,label="10^6",parse=T,size=1.8,color="gray25")+
  annotate("text",x=-110000000,y=8.1,label="10^8",parse=T,size=1.8,color="gray25")+
  coord_cartesian(xlim=c(-30362909,3036290862+30362909),ylim=c(4,8.1),expand=F,clip="off")+
  scale_y_continuous(expand=c(0.01,0.01),breaks=c(4,6,8))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border=element_rect(colour="black",fill=NA,linewidth=0.25),
        axis.text.x = element_blank(),
        #axis.text.y = element_markdown(size=unit(6,"pt"),margin=unit(c(0,1,0,-1),"pt")),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_markdown(size=unit(6,"pt"),margin=unit(c(0,0,0,0),"pt")),
        #axis.title.y = element_text(size=unit(7,"pt"),margin=unit(c(0,-2,0,0),"pt")),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        plot.margin=unit(c(0,2,0,2),"pt"),
        legend.text = element_text(size=unit(6,"pt"),margin=margin(t=0,b=0,l=-4,r=4))
  )

chrlen$col = c(rep(c("white","black"),12))
myplotdat = filter(chrlen,Chr!="chrY")
pos = (myplotdat$aStart+myplotdat$aEnd)/2
myrect = ggplot(myplotdat) + 
  geom_rect(aes(xmin=aStart,xmax=aEnd,ymin=1-0.5,ymax=1+0.5,fill=col))+
  scale_fill_manual(values=c("black","gray60"))+
  scale_x_continuous(expand = c(0.01,0.01))+
  annotate("text",pos[c(1:18,20,22:23)],y=c(rep(1,21)),label=c(1:18,20,22,"X"),color="white",size=1.8)+
  theme(panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        plot.margin=unit(c(-2,1,-2,1),"pt"))

pdf(file="Figure1a_part1.pdf",width=11/2.54,height=2.4/2.54)
ggarrange(P_Inv1_cns,P_Inv1_BP,myrect,ncol=1,nrow=3,heights=c(1.3,1,0.21),align="v",hjust=-2)
dev.off()


##### Part 2 #####
# plot cytoband
chr5 = filter(cytoband,chrom=="chr5")
P_cytoband5 = ggplot(chr5,aes(xmin=chromStart/1000000,xmax=chromEnd/1000000,ymin=0,ymax=1,fill=gieStain))+geom_rect()+
  scale_fill_manual(values=c("firebrick2","black","gray25","gray50","gray75","gray95"))+
  geom_rect(aes(xmin=0,xmax=49000000/1000000,ymin=0,ymax=1),col="black",fill=NA,size=0.5)+
  coord_cartesian(xlim=c(0,49000000)/1000000,expand = F,clip = 'on')+
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
    plot.margin=unit(c(0,4,1,5),"pt")
  )

chr12 = filter(cytoband,chrom=="chr12")
P_cytoband12 = ggplot(chr12,aes(xmin=chromStart/1000000,xmax=chromEnd/1000000,ymin=0,ymax=1,fill=gieStain))+geom_rect()+
  scale_fill_manual(values=c("firebrick2","black","gray25","gray50","gray75","gray95"))+
  geom_rect(aes(xmin=35800000/1000000,xmax=133851895/1000000,ymin=0,ymax=1),col="black",fill=NA,size=0.5)+
  coord_cartesian(xlim=c(35800000,133851895)/1000000,expand = F,clip = 'on')+
  xlab("12q (Mbps)")+
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
    plot.margin=unit(c(0,4,1,5),"pt")
  )

# plot somatic mutations
Inv1_mut_5 = filter(Inv1_mut,Chr == "chr5" & Start < 49000000)
P_Inv1_5_mut = ggplot(Inv1_mut_5)+ 
  geom_point(aes(aStart,log10(mutDis)),size=0.3) +
  #ylab("Intermut dist.<br>(bps)")+
  #ylab("Intermut dis\n(log10)")+
  coord_cartesian(ylim=c(0,8.1),xlim=c(881626701,881626701+49000000),expand=F,clip="off")+
  annotate("text",x=881626701-5500000,y=4.1,label="10^4",parse=T,size=1.8,color="gray25")+
  annotate("text",x=881626701-5500000,y=0,label=0,parse=T,size=2,color="gray25")+
  annotate("text",x=881626701-5500000,y=8.1,label="10^8",parse=T,size=1.8,color="gray25")+
  #annotate("text",x=881626701-9000000,y=4.1,label="Inter-mut dist.",parse=F,size=2,color="black",angle="90")+
  scale_y_continuous(breaks=c(0,4,8))+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.5),
        axis.title.y = element_blank(),
        #axis.title.y = element_text(size=unit(7,"pt"),margin=unit(c(0,-2,0,0),"pt")),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        plot.margin=unit(c(1,2,0,2),"pt")
  )

Inv1_mut_12 = filter(Inv1_mut,Chr == "chr12" & Start >35800000)
Inv1_mut_12$SNV_pattern2 = factor(as.character(Inv1_mut_12$SNV_pattern),levels=c("C>A","C>G","C>T","T>A","T>C","CC>TT"))
P_Inv1_12_mut = ggplot(Inv1_mut_12)+ 
  geom_point(aes(Start,log10(mutDis),color=SNV_pattern2),size=0.3) +
  ylab("Intermut dis<br>(log<sub>10</sub>)")+
  scale_y_continuous(limits=c(0,8.1),breaks=c(0,4,8),labels=c(0,4,8))+
  scale_x_continuous(limits=c(35800000,133851895),expand = c(0,0))+
  scale_color_manual(values=c("deepskyblue3","black","orangered","mediumpurple1","yellow3","orange"))+
  #guides(colour = guide_legend(nrow = 1,override.aes = list(size = 2)))+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        #axis.text.y = element_text(size=unit(6,"pt"),margin=unit(c(0,1,0,-1),"pt")),
        axis.text.y = element_blank(),
        axis.ticks=element_blank(),
        axis.title.x = element_blank(),
        axis.line.x = element_line(size=0.5),
        #axis.title.y = element_text(size=unit(7,"pt"),margin=unit(c(0,-2,0,0),"pt")),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        plot.margin=unit(c(1,2,0,2),"pt")
  )

# plot copy number
P_Inv1_5_cns = ggplot(filter(Inv1_cnr,chromosome=="chr5" & col=="gray" & (1:nrow(Inv1_cnr))%%2 == 0),aes((start+end)/2/1000000,log2)) +
  #geom_hline(yintercept=0:6,linetype="dashed")+
  ylab("Copy<br>number")+
  geom_point(color="gray",size=0.02)+
  geom_segment(aes(x=start/1000000,xend=end/1000000,y=log2,yend=log2),color="purple2",lineend="round",size=0.6,data=filter(Inv1_cns,chromosome=="chr5"))+
  geom_point(data=filter(Inv1_cnr,chromosome=="chr5" & col=="hotpink"),aes((start+end)/2/1000000,log2),color="hotpink",size=0.3)+
  coord_cartesian(xlim=c(0,49000000)/1000000,expand = F,clip = 'on')+
  scale_y_continuous(limits=c(-2.5,6.1),breaks=c(0,3,6),labels=c(2,16,128))+
  theme_bw()+
  theme(#axis.title.y=element_text(size=8,margin=unit(c(0,2,0,0),"pt")),
    panel.border = element_blank(),
    axis.line.x = element_line(size=0.5),
    axis.title.y = element_markdown(size=unit(6,"pt"),margin=unit(c(0,3,0,-2),"pt")),
    axis.text.y = element_text(size=unit(5,"pt"),margin=unit(c(0,2,0,-3),"pt")),
    axis.ticks.x=element_blank(),
    axis.text.x=element_blank(),
    axis.title.x=element_blank(),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin=unit(c(3,4,0,5),"pt")
  )

P_Inv1_12_cns = ggplot(filter(Inv1_cnr,chromosome=="chr12" & col=="gray" & (1:nrow(Inv1_cnr))%%8 == 0),aes((start+end)/2/1000000,log2)) +
  geom_point(color="gray",size=0.02)+
  geom_segment(aes(x=start/1000000,xend=end/1000000,y=log2,yend=log2),color="purple2",lineend="round",size=0.6,data=filter(Inv1_cns,chromosome=="chr12"))+
  ylab("CN")+
  geom_point(data=filter(Inv1_cnr,chromosome=="chr12" & col=="hotpink"),aes((start+end)/2/1000000,log2),color="hotpink",size=0.3)+
  coord_cartesian(xlim=c(35800000,133851895)/1000000,expand = F,clip = 'on')+
  scale_y_continuous(limits=c(-2.5,6.1),breaks=c(0,3,6),labels=c(2,16,128))+
  theme_bw()+
  theme(#axis.title.y=element_text(size=8,margin=unit(c(0,2,0,0),"pt")),
    panel.border = element_blank(),
    axis.line.x = element_line(size=0.5),
    axis.title.y = element_blank(),
    #axis.text.y = element_text(size=unit(5,"pt"),margin=unit(c(0,1,0,-1),"pt")),
    axis.text.y = element_blank(),
    axis.ticks=element_blank(),
    axis.text.x=element_blank(),
    axis.title.x=element_blank(),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin=unit(c(3,4,0,5),"pt")
  )

pdf(file="Figure1a_part2_5p.pdf",width=3.5/2.54,height=3/2.54)
ggarrange(P_Inv1_5_cns,P_Inv1_5_mut,P_cytoband5,align="v",
          ncol=1,nrow=3,heights=c(1.5,1,1.11))
dev.off()

pdf(file="Figure1a_part2_12q.pdf",width=5.02/2.54,height=3/2.54)
ggarrange(P_Inv1_12_cns,P_Inv1_12_mut,P_cytoband12,align="v",
          ncol=1,nrow=3,heights=c(1.5,1,1.11))
dev.off()
