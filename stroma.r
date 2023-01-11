
load('F5_stromal.RData')

source('/home/meisl/bin/bin/bin/source.R')

set.seed(36)

anoS=as.factor(anoS)
#anoS=ordered(as.factor(anoS),levels=c('Mono1','Mono2','Mono3','Macrophage1','Macrophage2','Macrophage3','mDC'))
anoS.pal <- setNames(sample(rainbow(length(levels(anoS)))),levels(anoS));

anoS.palf <- function(n) return(anoS.pal)
a2=scon$plotGraph(groups=anoS,plot.na=F,palette=anoS.palf,size=0.3,alpha=0.2,font.size = c(5, 5.5))
a2



cname=names(anoS)
cname = intersect(cname,names(ssamp))
ano2=data.frame('Cell'=anoS[cname],'SampleType'=ssamp[cname])

# Annotation vs sample
tmp2 <- acast(ano2, Cell ~ SampleType, fun.aggregate=length)
tmp3 <- (sweep(tmp2, 2, colSums(tmp2), FUN='/'))
tmp4 <- melt(tmp3)
head(tmp4)
names(tmp4) <- c('cell', 'sample','pc.of.sample')

tmp4$Group=NULL
tmp4$Group=sample.groups[as.character(tmp4$sample)]


p <- ggplot(na.omit(tmp4),aes(x=cell,y=pc.of.sample,dodge=Group,fill=Group))+geom_boxplot(notch=FALSE,outlier.shape=NA)  +  geom_point(position = position_jitterdodge(jitter.width=0.1),color=adjustcolor(1,alpha=0.3),pch=19,size=0.5)+theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_text(angle = 90, hjust = 0.5))  +xlab("") +ylab("fraction of total stromal cells")+theme(legend.position="top")+
  scale_fill_manual(values=fraction.palette1)
p





library(ggpubr)

df=tmp4
rsig=NULL
for (i in unique(df[,1])){
  tmp=df[df[,1]==i,]

  sig=compare_means(pc.of.sample ~ Group,  data = tmp) #
  sig$cell=i
  rsig=rbind(rsig,sig)
  #sig[sig$p.signif!='ns',]
}



gs = c('SELE','SELP','CLU','ADIRF','EDN1','IGFBP3','PLLP','FBLN5','MCAM','KDR','ENG','CSPG4')




