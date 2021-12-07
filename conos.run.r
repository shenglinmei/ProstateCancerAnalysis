
args <- commandArgs(trailingOnly = T);

fin=args[1]
appname <- args[2]


runDE.cluster = function (con, p2, appname) 
{
    #con = readRDS(fc)
    #p2 = readRDS(fp)
    clu = con$clusters$leiden$groups
    folderN = paste(appname, ".DE", sep = "")
    pwd = getwd()
    pwd2 = paste(pwd, "/", folderN, "/", sep = "")
    system(paste("mkdir ", folderN))
    setwd(pwd2)
    DEcaculate2(p2, appname, clu, removeGene = NULL, cutoff = 2, 
        num = 100, GO = NULL)
}

raw2=readRDS(fin)

raw2=lapply(raw2,function(x) x[,colSums(as.matrix(x))>800])

p2lis2=lapply(raw2,function(x) basicP2proc(x,n.cores = 10,min.transcripts.per.cell =0))

#saveRDS(p2lis2,'p2lis2.rds')

c1=runConos(p2lis2,appname)

datraw=raw2
	genelists <- lapply(datraw, function(x) rownames(x))
	str(genelists)
	commongenes <- Reduce(intersect,genelists)
	length(commongenes)


	matrices_raw <- mapply(function(mm, name) {
	 mm[commongenes,]
	},
	datraw,
	names(datraw))


bigM2 <- Reduce(cbind, matrices_raw)

p2=basicP2proc(bigM2,min.cells.per.gene = 0,n.cores = 10)

#p2app = p2app4conos(conos=c1, file=paste(appname,".bin",sep=''), save=TRUE)



 f1=paste(appname,'_p2combined.rds',sep='')
  saveRDS(p2,f1)



runDE.cluster(c1,p2,appname)
