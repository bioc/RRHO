## Author: Jonathan Rosenblatt and Jason Stein 

## Suggest default step size
defaultStepSize <-function(list1, list2){
  n1<- dim(list1)[1]
  n2<- dim(list2)[1]
  result <- ceiling(min(sqrt(c(n1,n2))))	
  return(result)
}	
### Testing:
## list.length <- 100
## list.names <- paste('Gene',1:list.length, sep='')
## gene.list1<- data.frame(list.names, sample(100))
## gene.list2<- data.frame(list.names, sample(100))
## defaultStepSize(gene.list1, gene.list2)


## Compute the overlaps between two *numeric* lists:
numericListOverlap<- function(sample1, sample2, stepsize, alternative, tol=0.5){
  n<- length(sample1)
  
  overlap<- function(a,b) {
    count<-as.integer(sum(as.numeric(sample1[1:a] %in% sample2[1:b])))
    
  switch(alternative,
           enrichment={
             log.pval<- -phyper(q=count-1, m=a, n=n-a+1, k=b, lower.tail=FALSE, log.p=TRUE)         
             signs<- 1L
           },
           two.sided={
             the.mean<- a*b/n
             signs<- sign(count - the.mean)
             if(signs < 0){
               lower<- count 
               upper<- 2*the.mean - count 
             } else{
               lower<- 2*the.mean - count 
               upper<- count 
             }
             
             log.pval<- -log(
               phyper(q=lower+tol, m=a, n=n-a+1, k=b, lower.tail=TRUE) +
                 phyper(q= upper-tol, m=a, n=n-a+1, k=b, lower.tail=FALSE))                               
           })
    
    return(c(counts=count, 
             log.pval=as.numeric(log.pval),
             signs=as.integer(signs)))    
  }
  
  
  indexes<- expand.grid(i=seq(1,n,by=stepsize), j=seq(1,n,by=stepsize))
  overlaps<- apply(indexes, 1, function(x) overlap(x['i'], x['j']))
  
  nrows<- sqrt(ncol(overlaps))
  matrix.counts<- matrix(overlaps['counts',], ncol=nrows)  
  matrix.log.pvals<- matrix(overlaps['log.pval',], ncol=nrows)  
  matrix.signs<- matrix(overlaps['signs',], ncol=nrows)  
  
  return(list(counts=matrix.counts, 
              log.pval=matrix.log.pvals,
              signs= matrix.signs))  
}
### Testing:
# n<- 112
# sample1<- sample(n)
# sample2<- sample(n)  
# .test<- RRHO:::numericListOverlap(sample1, sample2, 10)
# dim(.test$log.oval)
# library(lattice)
# levelplot(.test$counts)
# levelplot(.test$log.pval)
# table(is.na(.test$log.pval))








## Rank Rank Hypergeometric Overlap 
## based on Plaisier et al., Nucleic Acids Research, 2010
RRHO <- function(list1, list2, 
                 stepsize=defaultStepSize(list1, list2), 
                 labels, 
                 alternative,
                 plots=FALSE, 
                 outputdir=NULL, 
                 BY=FALSE,
                 log10.ind=FALSE) {
  ## list 1 is a data.frame from experiment 1 with two columns, 
  ## column 1 is the Gene Identifier, 
  ## column 2 is the signed ranking value (e.g. signed -log(p-value) 
  ##        or fold change)
  ##
  ## list 2 is a data.frame from experiment 2 with two columns, 
  ## column 1 is the Gene Identifier, 
  ## column 2 is the signed ranking value (e.g. signed -log10(p-value) 
  ##    or fold change)
  ## stepsize indicates how many genes to increase by 
  ##    in each algorithm iteration
  
  if (length(list1[,1]) != length(unique(list1[,1])))
    stop('Non-unique gene identifier found in list1')
  if (length(list2[,1]) != length(unique(list2[,1])))
    stop('Non-unique gene identifier found in list2')
  if(plots && (missing(outputdir) || missing(labels)))
    stop('When plots=TRUE, outputdir and labels are required.')
  if(!(alternative=='two.sided' || alternative=='enrichment'))
    stop('Wrong alternative specified.')
  
  result <-list(hypermat=NA, 
                hypermat.counts=NA, 
                hypermat.signs=NA,
                hypermat.by=NA, 
                n.items=nrow(list1), 
                stepsize=stepsize, 
                log10.ind=log10.ind,
                call=match.call()) 
  
  ## Order lists along list2
  list1  <- list1[order(list1[,2],decreasing=TRUE),]
  list2  <- list2[order(list2[,2],decreasing=TRUE),]
  nlist1 <- length(list1[,1])
  nlist2 <- length(list2[,1])
  
  ## Number of genes on the array
  N  <- max(nlist1,nlist2)
  
  .hypermat<- numericListOverlap(list1[,1], list2[,1], stepsize, alternative)
  hypermat<- .hypermat$log.pval
    
  ## Convert hypermat to a vector and Benjamini Yekutieli FDR correct
  
  if(log10.ind) hypermat<- hypermat *log10(exp(1))  
    
  if(BY){
    hypermatvec  <- matrix(hypermat,
                           nrow=nrow(hypermat)*ncol(hypermat),ncol=1)
    hypermat.byvec  <- p.adjust(exp(-hypermatvec),method="BY")
    hypermat.by <- matrix(-log(hypermat.byvec),
                                             nrow=nrow(hypermat),ncol=ncol(hypermat))     
    
    if(log10.ind) hypermat.by<- hypermat.by *log10(exp(1))
    result$hypermat.by<- hypermat.by
  }
  
    
  
  if(plots) {
    try({
    hypermat.signed<- hypermat * .hypermat$signs 
    
    ## Function to plot color bar
    ## Modified from http://www.colbyimaging.com/wiki/statistics/color-bars
    color.bar <- function(lut, min, max=-min, 
                          nticks=11, 
                          ticks=seq(min, max, len=nticks), 
                          title='') {
      scale  <- (length(lut)-1)/(max-min)
      plot(c(0,10), c(min,max), type='n', bty='n', 
           xaxt='n', xlab='', yaxt='n', ylab='')
      mtext(title,2,2.3, cex=0.8)
      axis(2, round(ticks,0), las=1,cex.lab=0.8)
      for (i in 1:(length(lut)-1)) {
        y  <- (i-1)/scale + min
        rect(0,y,10,y+1/scale, col=lut[i], border=NA)
      }
    }
    
    .filename <-paste("RRHOMap", labels[1], "_VS_", labels[2], ".jpg", sep="") 
    jpeg(filename = paste(outputdir,.filename,sep="/"), 
         width=8, height=8, 
         units="in", quality=100, res=150)
    
    jet.colors  <- colorRampPalette(
      c("#00007F", "blue", "#007FFF", "cyan", 
        "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    layout(matrix(c(rep(1,5),2), 1, 6, byrow = TRUE))
    
    image(hypermat.signed, xlab='', ylab='', col=jet.colors(100), 
          axes=FALSE, main="Rank Rank Hypergeometric Overlap Map")
    
        
    mtext(labels[2],2,0.5)
    mtext(labels[1],1,0.5)
    ##mtext(paste("-log(BY P-value) =",max(hypermat.by)),3,0.5,cex=0.5)
    
    finite.ind<- is.finite(hypermat.signed)
    color.bar(jet.colors(100),
              min=min(hypermat.signed[finite.ind], na.rm=TRUE),
              max=max(hypermat.signed[finite.ind], na.rm=TRUE),
              nticks=6,
              title="-log(P-value)")
    
    dev.off()
    
    ## Make a rank scatter plot
    list2ind  <- match(list1[,1],list2[,1])
    list1ind  <- 1:nlist1
    corval  <- cor(list1ind,list2ind,method="spearman")
    .filename <-paste("RankScatter",labels[1],"_VS_",labels[2],".jpg",sep="") 
    jpeg(paste(outputdir,.filename,sep="/"), width=8, 
         height=8, units="in", quality=100, res=150)
    plot(list1ind,list2ind,xlab=paste(labels[1],"(Rank)"), 
         ylab=paste(labels[2],"(Rank)"), pch=20, 
         main=paste(
           "Rank-Rank Scatter (rho = ",signif(corval,digits=3),")"
           ,sep=""), cex=0.5)
    ## TODO: Replace linear fit with LOESS
    model  <- lm(list2ind~list1ind)
    lines(predict(model),col="red",lwd=3)
    dev.off()
    
    ## Make a Venn Diagram for the most significantly associated points
    ## Upper Right Corner (Downregulated in both)
    maxind.ur  <- which(
      max(hypermat.signed[ceiling(nrow(hypermat.signed)/2):nrow(hypermat.signed),
                   ceiling(ncol(hypermat.signed)/2):ncol(hypermat.signed)],
          na.rm=TRUE) == hypermat.signed, 
      arr.ind=TRUE)
    indlist1.ur  <- seq(1,nlist1,stepsize)[maxind.ur[1]]
    indlist2.ur  <- seq(1,nlist2,stepsize)[maxind.ur[2]]
    genelist.ur  <- intersect(
      list1[indlist1.ur:nlist1,1],
      list2[indlist2.ur:nlist2,1])
    ## Lower Right corner (Upregulated in both)
    maxind.lr  <- which(
      max(hypermat.signed[1:(ceiling(nrow(hypermat.signed)/2)-1), 
                   1:(ceiling(ncol(hypermat.signed)/2)-1)],
          na.rm=TRUE) == hypermat.signed, arr.ind=TRUE)
    indlist1.lr  <- seq(1,nlist1,stepsize)[maxind.lr[1]]
    indlist2.lr  <- seq(1,nlist2,stepsize)[maxind.lr[2]]
    genelist.lr  <- intersect(
      list1[1:indlist1.lr,1], 
      list2[1:indlist2.lr,1])
    
    ## Write out the gene lists of overlapping
    .filename <- paste(
      outputdir,"/RRHO_GO_MostDownregulated",labels[1],"_VS_",labels[2],".csv",
      sep="")
    write.table(genelist.ur,.filename,row.names=F,quote=F,col.names=F)
    .filename <- paste(
      outputdir,"/RRHO_GO_MostUpregulated",labels[1],"_VS_",labels[2],".csv",
      sep="")
    write.table(genelist.lr,.filename,row.names=F,quote=F,col.names=F)
    
    .filename <- paste(
      outputdir,"/RRHO_VennMost",labels[1],"__VS__",labels[2],".jpg", 
      sep="")
    jpeg(.filename,width=8.5,height=5,units="in",quality=100,res=150)
    vp1  <- viewport(x=0.25,y=0.5,width=0.5,height=0.9)
    vp2  <- viewport(x=0.75,y=0.5,width=0.5,height=0.9)
    
    pushViewport(vp1)
    h1  <- draw.pairwise.venn(length(indlist1.ur:nlist1),
                              length(indlist2.ur:nlist2),
                              length(genelist.ur), 
                              category=c(labels[1],labels[2]),
                              scaled=TRUE,
                              lwd=c(0,0),
                              fill=c("cornflowerblue", "darkorchid1"),
                              cex=1,
                              cat.cex=1.2,
                              cat.pos=c(0,0),
                              ext.text=FALSE,
                              ind=FALSE,
                              cat.dist=0.01)
    grid.draw(h1)
    grid.text("Down Regulated",y=1)
    upViewport()
    pushViewport(vp2)
    h2  <-  draw.pairwise.venn(length(1:indlist1.lr),
                               length(1:indlist2.lr),
                               length(genelist.lr),
                               category=c(labels[1],labels[2]),
                               scaled=TRUE,
                               lwd=c(0,0),
                               fill=c("cornflowerblue", "darkorchid1"),
                               cex=1,
                               cat.cex=1.2,
                               cat.pos=c(0,0),
                               ext.text=FALSE,
                               main="Negative",
                               ind=FALSE,
                               cat.dist=0.01)
    grid.draw(h2)
    grid.text("Up Regulated",y=1)
    dev.off()
  })
  if(length(h2)==0L) message('Unable to output JPG plots.')
  }



  result$hypermat <- hypermat
  result$hypermat.counts <- .hypermat$counts
  result$hypermat.signs <- .hypermat$signs
  
  return(result)
}
### Testing:
# list.length <- 100
# list.names <- paste('Gene',1:list.length, sep='')
# gene.list1<- data.frame(list.names, sample(list.length))
# gene.list2<- data.frame(list.names, sample(list.length))
# RRHO.example <-  RRHO(gene.list1, gene.list2, alternative = "enrichment")
# library(lattice)
# levelplot(RRHO.example$hypermat)
# RRHO.example <-  RRHO(gene.list1, gene.list2, alternative = "enrichment", BY=TRUE, log10.ind = FALSE)
# RRHO.example <-  RRHO(gene.list1, gene.list2, alternative = "two.sided", BY=TRUE, log10.ind = TRUE)
# levelplot(RRHO.example$hypermat.by)
# RRHO.example <-  
#   RRHO(gene.list1, 
#        gene.list2, 
#        plots=TRUE, outputdir=tempdir(), labels=c("a","b"))




## TODO: Function for FWER control using permutations
pvalRRHO <- function(RRHO.obj, 
                     replications, 
                     stepsize=RRHO.obj$stepsize, 
                     FUN= max){
  ## RRHO.obj <- RRHO.example
  ## FUN<- max
  ## replications<- 100
  ## stepsize <- RRHO.obj$stepsize
  ## result <- list(FUN=FUN, n.items=n.items, stepsize=stepsize)
  ## Note: min(pvals) maps to max(-log(pvals))
  
  interactive.ind<- interactive()
  if(interactive.ind) {
    message('This might take a while')
    pb <- txtProgressBar(min = 0, max = replications, style = 3)
  }
  
  n.items <- RRHO.obj$n.items
  alternative<- RRHO.obj$call$alternative
  log10.ind<- RRHO.obj$log10.ind
  
  result <- list(FUN=FUN, 
                 n.items=n.items, 
                 stepsize=stepsize , 
                 replications= replications, 
                 alternative=alternative,
                 call=match.call())
  
  list.names <- paste('Gene',1:n.items, sep='')
  FUN.vals<- rep(NA, replications)
  for(i in 1:replications){
    ## i<- 1
    ## Generate rankings and compute overlap
    sample1<- data.frame(list.names, sample(n.items))
    sample2<- data.frame(list.names, sample(n.items))	  
    .RRHO<- RRHO(sample1, sample2, stepsize=stepsize, plots=FALSE, BY=FALSE, alternative=alternative, log10.ind=FALSE)
    .clean.result<- na.omit(.RRHO$hypermat)
    FUN.vals[i]<- FUN(.clean.result)
    
    if(interactive.ind) setTxtProgressBar(pb, i)	  
  }
  
  ## Adding a conservative constant in case there were not enough replications.
  FUN.ecdf<- function(x)  min( ecdf(FUN.vals)(x) + 1/replications, 1)
  result$FUN.ecdf<- FUN.ecdf	  
  
  .clean.data<- na.omit(RRHO.obj$hypermat)
  
  FUN.observed<- FUN(.clean.data )
  if(log10.ind) FUN.observed<- FUN.observed / log10(exp(1))
  
  result$pval<- 1-FUN.ecdf(FUN.observed)
  
  if(interactive.ind) close(pb)
  ## Return pvale
  return(result)
}
### Testing:
# pvalRRHO(RRHO.example, replications = 1000)