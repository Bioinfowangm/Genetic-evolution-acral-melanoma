library(gdata)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(ggpattern)
library(readxl)

freq = read_excel("./extdata/f1b_Hailstorms_distribution.xlsx")
freq$chr = factor(freq$chr,levels=c(freq$chr))
freq$col1 = factor(freq$col1)
mychr = ggplot(freq,aes(index,count/37*100))+geom_bar(stat="identity",fill="gray30")+
  scale_x_continuous(expand = c(0.01,0.01),breaks=freq$index,labels=freq$col2)+
  ylab("Cases with\nhailstorm (%)")+
  theme_classic2()+
  theme(
    axis.title.x = element_blank(),
    plot.title=element_text(size=8,hjust=0.5),
    axis.title.y = element_text(size=7),
    axis.line = element_line(size=0.3),
    axis.text = element_text(size=5),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    plot.margin=unit(c(3,3,1,3),"pt")
  )

a = c(1.5,2.5,2); b= c(1.5,1.5,2)
newdat = data.frame(a=c(a,a+2,a+4,a+6,a+8,a+10,a+12,a+14,a+16,a+18,a+20,a+22),
                    b=rep(b,12),
                    t=rep(1:12,each=3))
q = data.frame(x=c(1:12*2,25:27,14:18*2+1,38:39,41),y=1.8,col="blue")
p = data.frame(x=c(1:12*2-1,14:18*2,40),y=1.8,col="red")
arms = rbind(q,p)
myrect1 = ggplot(freq) + 
  geom_rect(aes(xmin=index-0.5,xmax=index+0.5,ymin=1-0.5,ymax=1+0.5,fill=col1))+
  geom_rect(aes(xmin=x-0.5,xmax=x+0.5,ymin=1+0.5,ymax=1+0.85,fill=col),data=arms)+
  scale_fill_manual(values=c("gray60","black",alpha("red",0.3),alpha("blue",0.6)))+
  scale_x_continuous(expand = c(0.01,0.01))+
  scale_y_continuous(limits=c(0.5,2.2),expand = c(0,0))+
  annotate("text",c(1:12*2-0.5,25:27,14:18*2+0.5,0.5+37:38+0.5,0.5+40),y=c(rep(1,23)),label=c(1:22,"X"),color="white",size=1.8)+
  theme(axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        plot.margin=unit(c(-5,0,0,0),"pt"))

pdf(file="Figure1b.pdf",width=10.6/2.54,height=1.8/2.54)
ggarrange(mychr,myrect1,ncol=1,nrow=2,heights=c(1,0.14),align="v")
dev.off()


##### Fig 1c #####
onTree = data.frame(Group = c("Trunk","Branch"),Number = c(59,5))
pdf(file="Figure1c.pdf",width=3.8/2.54,height=5.2/2.54)
onTree$Group = factor(onTree$Group,levels=c("Trunk","Branch"))
ggplot(onTree,aes(Group,Number/64 * 100))+ geom_bar(fill="gray30",stat="identity")+
  ylab("Hailstorms (%)")+ ylim(0,100)+
  theme_classic2()+
  theme(
    axis.title.x = element_blank(),
    #plot.title=element_text(size=8,margin=margin(0,0,0,0),hjust=0.5),
    plot.title=element_blank(),
    axis.title.y = element_text(size=7,margin=margin(0,-1,0,1)),
    axis.text.y = element_text(size=5),
    axis.text.x = element_text(size=6,margin=margin(1,0,0,0)),
    axis.line = element_line(linewidth=0.3),
    plot.margin=unit(c(3,3,1,3),"pt")
  )
dev.off()

