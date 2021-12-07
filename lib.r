

library(reshape2)
library(xtable)
library(pagoda2)
library(nbHelpers)
library(reshape2)
library(Matrix)
#library(clusterMatch)
library(igraph)
library(ggplot2)
#library(nbHelpers)                  
library(parallel)
#library(fastSave)
library(pagoda2)
#library(nbHelpers)
library(reshape2)
library(pheatmap)
library(conos)
#library('betareg')
library(DESeq2)
library(Matrix)
library(conos)
library(GOstats)
library(RSQLite)
library(biomaRt)
suppressPackageStartupMessages(library(org.Hs.eg.db))

library(ggrepel)

addZScores <- function(df) {
  df$Z <- -qnorm(df$pval/2)
  df$Z[is.na(df$Z)] <- 0
  df$Za <- -qnorm(df$padj/2)
  df$Za[is.na(df$Za)] <- 0
  df$Z <- df$Z * sign(df$log2FoldChange)
  df$Za <- df$Za * sign(df$log2FoldChange)
  
  return(df)
}

runDE=function(count,group,ref){
  print(table(group))
  group = as.factor(group)
  cname=intersect(colnames(count),names(group))
  count=count[,cname]
  group= group [cname]
  cm=as.matrix(count[,names(group)])
  meta <- data.frame(sample.id = colnames(cm), group =group)
  dds1 <- DESeq2::DESeqDataSetFromMatrix(round(cm,0), meta,
                                         design = ~group)
  meta$group <- relevel(meta$group, ref = ref)
  
  dds1 <- DESeq2::DESeq(dds1)
  res1 <- DESeq2::results(dds1, cooksCutoff = FALSE,
                          independentFiltering = FALSE)
  res1 <- as.data.frame(res1)
  res1 <- res1[order(res1$padj, decreasing = FALSE), ]
  
  res1=addZScores(res1)
  
  return(res1)
}




strpart <- function (x, split, n, fixed = FALSE) {
  sapply(strsplit(as.character(x), split, fixed = fixed), "[",n)
}



rbindDEMatrices <- function(mats, cluster.sep.chr) {
  mats <- lapply(names(mats), function(n) {
    rownames(mats[[n]]) <- paste0(n, cluster.sep.chr, rownames(mats[[n]]));
    return(mats[[n]])
  })
  
  return(t(do.call(rbind, mats)))
}

collapseCellsByType <- function(cm, groups, min.cell.count=10) {
  groups <- as.factor(groups);
  cl <- factor(groups[match(rownames(cm),names(groups))],levels=levels(groups));
  # TODO remove dependency on Conos
  tc <- conos:::colSumByFactor(cm,cl);
  tc <- tc[-1,,drop=FALSE]  # omit NA cells
  tc[table(cl)>=min.cell.count,,drop=FALSE]
}


rawMatricesWithCommonGenes=function(raw.mats)
{
  common.genes <- Reduce(intersect, lapply(raw.mats, colnames))
  print(length(common.genes))
  return(lapply(raw.mats, function(x) {
    x[, common.genes]
  }))
}





DotPlot.DE <- function(de1,num=10,ylab='Treat vs control',gn1=NULL,gn2=NULL,orderpvalue=NULL,fc=2,pval=NULL){
  library(ggplot2)
  library(ggrepel)
  
  tmp=de1
  tmp=tmp[!is.na(tmp$pvalue),]
  tmp$score=-log10(tmp$padj)
  if (!is.null(pval)){
    tmp$score=-log10(tmp$pvalue)
  }
  
  gg=data.frame('log2FoldChange'=tmp$log2FoldChange,'score'=tmp$score,'Z'=tmp$Z,'name'=rownames(tmp),
                'padj'=tmp$padj,'pvalue'=tmp$pvalue,'score2'=tmp$score*sign(tmp$log2FoldChange))
  p=ggplot(gg,aes(y=score,x=log2FoldChange))
  p=p+geom_point(shape=".",alpha=1/1,size = 1,color='grey')+theme_classic()
  p=p+labs(y=ylab)
  gg$score2=sign(gg$log2FoldChange)*gg$score
  rownames(gg)=gg$name 
  
  if (is.null(orderpvalue)){
    t1=gg[order(gg$log2FoldChange),]
    t2=gg[order(gg$log2FoldChange,decreasing=TRUE),]
    
  }else{
    t1=gg[order(gg$score2),]
    t2=gg[order(gg$score2,decreasing=TRUE),]
    
  }
  
  t1=t1[t1$log2FoldChange<(-1),]
  n1=rownames(t1)[1:num]
  
  t2=t2[t2$log2FoldChange>1,]
  n2=rownames(t2)[1:num]
  
  if (!is.null(gn1)){
    n1=c(n1,gn1)
  }
  if (!is.null(gn2)){
    n1=c(n2,gn2)
  } 
  
  
  
  p=p+geom_text_repel(
    data = gg[gg$name %in% n1,],
    aes(label = name),
    size = 3,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    color='blue'
  )
  
  
  
  
  p=p+geom_text_repel(
    data = gg[gg$name %in% n2,],
    aes(label = name),
    size = 3,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    color='red1'
  )
  
  
  p=p+ geom_point(data = subset(gg,(log2FoldChange< (-fc) & pvalue < 0.01)),shape=".",alpha=1/1,size = 2,color='blue')
  p=p+ geom_point(data = subset(gg,(log2FoldChange>fc & pvalue < 0.01)),shape=".",alpha=1/1,size = 2,color='red1')
  
  return(p)          
  
}



mergeDat2=function(raw.mats){
  
  genelists <- lapply(raw.mats, function(x) rownames(x))
  str(genelists)
  commongenes <- Reduce(intersect,genelists)
  
  
  matrices2 <- mapply(function(m, name) {
    #   colnames(m) <- paste(name, colnames(m), sep='_');
    m[commongenes,]
  },
  raw.mats,
  names(raw.mats))
  cellNum=unlist(lapply(matrices2, function(x) ncol(x)))
  print('#commongenes')
  print(length(commongenes))
  
  print('#cell numbers')
  print(cellNum)
  bigM2 <- Reduce(cbind, matrices2)
  
  return(bigM2)
}



GOanalysis=function(markers,n){
  
  ENTREZID=unlist(mget(markers, org.Hs.egSYMBOL2EG, ifnotfound=NA))
  ENTREZID=ENTREZID[!is.na(ENTREZID)]
  
  
  for(function_type in c("BP", "CC", "MF")){
    
    param <- new("GOHyperGParams", geneIds=ENTREZID,
                 #universe=universe,
                 annotation="org.Hs.eg.db", ontology=function_type,pvalueCutoff=0.1,
                 conditional=FALSE, testDirection="over")
    hyp <- hyperGTest(param)
    sumTable <- summary(hyp)
    
    
    
    david=sumTable[1:20,]
    david$Pvalue=-log(david[,2])
    termNumber=nrow(david)
    
    
    library(ggplot2)
    p1 <- ggplot(data=david, aes(x=Pvalue, y=Term, size=Count, colour = factor(david$Count)))
    p1 <- p1 + geom_point()
    p1 <- p1 + guides(color = FALSE)
    p1 <- p1 + theme(panel.grid.major = element_line(colour='blue'),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank())
    p1 <- p1 + xlab(paste("-log10(Pvalue)", sep="")) + ylab("")
    p1 <- p1 + labs(title=paste("DAVID:", function_type, sep=""))
    p1 <- p1 + theme(axis.text.x=element_text(size=10, face="plain", colour ='black'))
    p1 <- p1 + theme(axis.text.y=element_text(size=6, face="plain", colour ='black'))
    #p1=p1+theme(axis.title.y=element_text(size=10,face="plain",colour ='black'))
    #p1=p1+theme(legend.position="bottom")
    p1 <- p1 + xlim(min(david$Pvalue), max(david$Pvalue))
    #p1=p1+ylim(-10,+15)
    print(p1)
    ggsave(file=paste(n,'_' ,function_type, ".png", sep=""), scale=0.8, dpi=600, width = 7, height=1+0.25*termNumber) 
  } 
  
}



GOanalysis.term.mouse=function(markers,n){
  
  ENTREZID=unlist(mget(markers, org.Mm.egSYMBOL2EG, ifnotfound=NA))
  ENTREZID=ENTREZID[!is.na(ENTREZID)]
  
  allr=list()
  
  for(function_type in c("BP", "CC", "MF")){
    
    param <- new("GOHyperGParams", geneIds=ENTREZID,
                 #universe=universe,
                 annotation="org.Mm.eg.db", ontology=function_type,pvalueCutoff=0.1,
                 conditional=FALSE, testDirection="over")
    hyp <- hyperGTest(param)
    sumTable <- summary(hyp)
    
    sumTable$p.adjust=p.adjust(sumTable$Pvalue, "BH")
    
    david=sumTable[1:20,]
    david$Pvalue=-log(david[,2])
    termNumber=nrow(david)
    david$genes=apply(david,1,function(x) { paste(names(ENTREZID[ENTREZID %in% get(as.character(x[1]),org.Mm.egGO2ALLEGS)]),collapse = ',' ) } )
    
    allr[[function_type]]=david
    write.table(david,paste(n,'_' ,function_type, ".xls", sep=""),sep='\t',col.names=T,row.names=F,quote=F)
    
    
    library(ggplot2)
    p1 <- ggplot(data=david, aes(x=Pvalue, y=Term, size=Count, colour = factor(david$Count)))
    p1 <- p1 + geom_point()
    p1 <- p1 + guides(color = FALSE)
    p1 <- p1 + theme(panel.grid.major = element_line(colour='blue'),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank())
    
    
    p1 <- p1 + xlab(paste("-log10(Pvalue)", sep="")) + ylab("")
    p1 <- p1 + labs(title=paste("DAVID:", function_type, sep=""))
    p1 <- p1 + theme(axis.text.x=element_text(size=10, face="plain", colour ='black'))
    p1 <- p1 + theme(axis.text.y=element_text(size=6, face="plain", colour ='black'))
    #p1=p1+theme(axis.title.y=element_text(size=10,face="plain",colour ='black'))
    #p1=p1+theme(legend.position="bottom")
    p1 <- p1 + xlim(min(david$Pvalue), max(david$Pvalue))
    #p1=p1+ylim(-10,+15)
    print(p1)
    ggsave(file=paste(n,'_' ,function_type, ".png", sep=""), scale=0.8, dpi=600, width = 7, height=1+0.25*termNumber)
  }
  
  
  saveRDS(allr,paste(n,'.GOterm.rds',sep=''))
  return(allr)
}


connectedBarplot <- function(dat, color=rainbow(nrow(dat)), space=1, alpha=0.2, ...) {  
  b <- barplot(dat, col=color,las=2, space = space, ...)                     
  
  for (i in seq_len(ncol(dat) - 1)) {     
    lines(c(b[i]+0.5, b[i+1]-0.5), c(0, 0)) ## bottom line       
    
    for (j in seq_len(nrow(dat))) {     
      if (j == 1) {                   
        lines(c(b[i]+0.5, b[i+1]-0.5), c(dat[j,i], dat[j,i+1]))                       
        polygon(c(b[i]+0.5, b[i]+0.5, b[i+1]-0.5, b[i+1]-0.5),                        
                c(0, dat[j,i], dat[j,i+1], 0),               
                col=adjustcolor(color[j], alpha.f=alpha))    
      }      
      if (j == 2) {                   
        lines(c(b[i]+0.5, b[i+1]-0.5), c(colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1]))                      
        polygon(c(b[i]+0.5, b[i]+0.5, b[i+1]-0.5, b[i+1]-0.5),                        
                c(dat[1,i], colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1], dat[1,i+1]),                       
                col=adjustcolor(color[j], alpha.f=alpha))    
      }      
      if (j > 2) {                    
        lines(c(b[i]+0.5, b[i+1]-0.5), c(colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1]))                      
        polygon(c(b[i]+0.5, b[i]+0.5, b[i+1]-0.5, b[i+1]-0.5),                        
                c(colSums(dat[1:(j-1),])[i], colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1], colSums(dat[1:(j-1),])[i+1]),              
                col=adjustcolor(color[j], alpha.f=alpha))    
      }      
    }          
  }              
}


