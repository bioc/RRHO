#Rank Rank Hypergeometric Overlap based on Plaisier et al., Nucleic Acids Research, 2010
#Compares two RRHO maps to determine significance
RRHOComparison <- function(list1, list2, list3, 
                            stepsize, plots=FALSE, labels, outputdir=NULL,
                           log10.ind=FALSE) {
  ## Significance testing of the difference between two RRHO maps.
  ## RRHO map 1: list1 vs list3
  ## RRHO map 2: list2 vs list3
  ##
  ## list 1 is a data.frame from experiment 1 with two columns, column 1 is the Gene Identifier, column 2 is the signed ranking value (e.g. signed -log10(p-value) or fold change)
  ## list 2 is a data.frame from experiment 2 with two columns, column 1 is the Gene Identifier, column 2 is the signed ranking value (e.g. signed -log10(p-value) or fold change)
  ## list 3 is a data.frame from experiment 3 with two columns, column 1 is the Gene Identifier, column 2 is the signed ranking value (e.g. signed -log10(p-value) or fold change).  
  ## stepsize indicates how many genes to increase by in each algorithm iteration
  
  if (length(list1[,1])!=length(unique(list1[,1]))) 
    stop('Non-unique gene identifier found in list1')
  if (length(list2[,1])!=length(unique(list2[,1]))) 
    stop('Non-unique gene identifier found in list2')
  if (length(list3[,1])!=length(unique(list3[,1])))
    stop('Non-unique gene identifier found in list3')
  
  list1 <- list1[order(list1[,2],decreasing=TRUE),]
  list2 <- list2[order(list2[,2],decreasing=TRUE),]
  list3 <- list3[order(list3[,2],decreasing=TRUE),]
  nlist1 <- length(list1[,1])
  nlist2 <- length(list2[,1])
  nlist3 <- length(list3[,1])
  
  ## Number of genes on the array
  N <- max(c(nlist1, nlist2, nlist3))
  
  hypermat1 <- matrix(data=NA, 
                      nrow=length(seq(1,nlist1,stepsize)),
                      ncol=length(seq(1,nlist3,stepsize)))
  OR1 <- matrix(data=NA,
                nrow=length(seq(1,nlist1,stepsize)),
                ncol=length(seq(1,nlist3,stepsize)))
  SE1 <- matrix(data=NA,
                nrow=length(seq(1,nlist1,stepsize)),
                ncol=length(seq(1,nlist3,stepsize)))
  countx <- county <- 0
  ##Loop over the experiments
  for (i in seq(1, nlist1, stepsize)) {
    countx <- countx + 1
    for (j in seq(1,nlist3,stepsize)) {
      county <- county + 1
      ## Parameters for the hypergeometric test
      k <- length(intersect(list1[1:i,1], list3[1:j,1]))
      s <- length(list1[1:i,1])
      M <- length(list3[1:j,1])
      ## Hypergeometric test converting to log10
      ## log10(x) = log_e(x) * log10(e)
      hypermat1[countx, county] <- -phyper(k-1, M, N-M, s, lower.tail=FALSE, log.p=TRUE) 
      contingencytable <- rbind(c(k, M-k),c(s-k, N-s-M+k))
      OR1[countx,county] <- (contingencytable[1,1]*contingencytable[2,2]) /
        (contingencytable[1,2]*contingencytable[2,1])
      
      SE1[countx,county] <- sqrt(1/contingencytable[1,1] +
                                   1/contingencytable[1,2] +
                                   1/contingencytable[2,1] +
                                   1/contingencytable[2,2])
    }
    county<-0
  }
  
  hypermat2 <- matrix(data=NA, 
                      nrow=length(seq(1,nlist2,stepsize)),
                      ncol=length(seq(1,nlist3,stepsize)))
  OR2 <- matrix(data=NA, 
                nrow=length(seq(1,nlist2,stepsize)),
                ncol=length(seq(1,nlist3,stepsize)))
  SE2 <- matrix(data=NA, 
                nrow=length(seq(1,nlist2,stepsize)),
                ncol=length(seq(1,nlist3,stepsize)))
  countx <- county <- 0
  ##Loop over the experiments
  for (i in seq(1, nlist2, stepsize)) {
    countx <- countx + 1
    for (j in seq(1, nlist3, stepsize)) {
      county <- county + 1
      ## Parameters for the hypergeometric test
      k <- length(intersect(list2[1:i,1], list3[1:j,1]))
      s <- length(list2[1:i,1])
      M <- length(list3[1:j,1])
      ## Hypergeometric test converting to log10
      ## log10(x) = log_e(x) * log10(e)
      hypermat2[countx,county] <- -phyper(k-1, M, N-M, s,
                                          lower.tail=FALSE, log.p=TRUE)
      contingencytable <- rbind(c(k,M-k), c(s-k,N-s-M+k))
      OR2[countx, county] <- (contingencytable[1,1]*contingencytable[2,2])/
        (contingencytable[1,2]*contingencytable[2,1])
      
      SE2[countx, county] <- sqrt(1/contingencytable[1,1] +
                                   1/contingencytable[1,2] +
                                   1/contingencytable[2,1] +
                                   1/contingencytable[2,2])
    }
    county<-0
  }
  
  
  ### Now testing of the difference between 1-3 is the same as 2-3.
  ## This is done by testing for the difference 
  ## in log Odds ratio between the two arrays at each location.
  
  ## Normal approximation to the odd ratios ratio
  ## Z = log(OR1)-log(OR2)/sqrt(SE1^2 + SE2^2)
  ORdiff <- log(OR1/OR2)
  SEdiff <- sqrt(SE1^2 + SE2^2)  
  Zdiff <- ORdiff/SEdiff
  
  ## Using a two tailed P-value, need to use log p-values
  ##  Pdiff = (1-pnorm(abs(Zdiff)))*2
  ## Return matrix of log of pvalue of odds ratio test.
  Pdiff <- (-1*pnorm(abs(Zdiff), log.p=TRUE, lower.tail=FALSE)-log(2))
  
  ## Convert hypermat to a vector and Benjamini Yekutieli FDR correct
  Pdiffvec <- matrix(Pdiff, nrow=nrow(Pdiff)*ncol(Pdiff), ncol=1)
  Zdiffvec <- matrix(Zdiff, nrow=nrow(Zdiff)*ncol(Zdiff), ncol=1)
  nonaind <- which(!is.na(Pdiffvec))
  Pdiff.byvec <- matrix(data=NA, nrow=length(Pdiffvec), ncol=1)
  Pdiff.byvec[nonaind] <- 
    sign(Zdiffvec)[nonaind] *
    -log(p.adjust(exp(-Pdiffvec[nonaind]), method="BY"))
  Pdiff.by <- matrix(Pdiff.byvec, nrow=nrow(Pdiff), ncol=ncol(Pdiff))
  
  
  
  if(log10.ind) {
    fact<- log10(exp(1))
    hypermat1<- hypermat1 * fact
    hypermat2<- hypermat2 * fact
    Pdiff<- Pdiff * fact
    Pdiff.by<- Pdiff.by * fact
  }
  
  
  
  
  
  if (plots) {
    ## Function to plot color bar
    ## Modified from http://www.colbyimaging.com/wiki/statistics/color-bars
      color.bar <- function(lut, min, max=-min, 
                            nticks=11, ticks=seq(min, max, len=nticks), title='') {
        scale <- (length(lut)-1)/(max-min)  
        plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='')
        mtext(title,2,2.3,cex=0.8)
        axis(2, round(ticks,0), las=1,cex.lab=0.8)
        for (i in 1:(length(lut)-1)) {
          y <- (i-1)/scale + min
          rect(0, y, 10, y+1/scale, col=lut[i], border=NA)
        }
      }
      
      
      ## Make first comparison plot.
      pdf(paste(
        outputdir,
        paste("RRHOMaps_Diff", labels[1],"__VS__", labels[3],".pdf", sep=""), sep="/")
      )
      jet.colors <- colorRampPalette(
        c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
      layout(matrix(c(rep(1,5),2), 1, 6, byrow = TRUE))
      
      image(hypermat1, xlab='',ylab='',
            col=jet.colors(100),
            axes=FALSE,
            main="Rank Rank Hypergeometric Overlap Map")
      mtext(labels[3],2,0.5)
      mtext(labels[1],1,0.5)
      color.bar(jet.colors(100),min=min(hypermat1,na.rm=TRUE),
                max=max(hypermat1,na.rm=TRUE),nticks=6,title="-log(Nominal P-value)")
      dev.off()
      
      
      ## Make second comparison plot.
      pdf(paste(outputdir, paste("RRHOMaps_Diff",labels[2],"__VS__",labels[3],".pdf",sep=""),sep="/"))
      layout(matrix(c(rep(1,5),2), 1, 6, byrow = TRUE))
      image(hypermat2,xlab='',ylab='',col=jet.colors(100),
            axes=FALSE,main="Rank Rank Hypergeometric Overlap Map")
      mtext(labels[3],2,0.5)
      mtext(labels[2],1,0.5)
      color.bar(jet.colors(100),min=min(hypermat2,na.rm=TRUE),
                max=max(hypermat2,na.rm=TRUE),nticks=6,title="-log(Nominal P-value)")
      dev.off()
      
      
      ## Make diff 
      #Set max and mins of colorbar for in vitro systems used in paper
      minhypermat <- min(Pdiff.by, na.rm=TRUE)
      maxhypermat <- max(Pdiff.by, na.rm= TRUE)

      pdf(paste(outputdir,
                paste("RRHOMaps_DiffMap",
                                labels[1],"__VS__",labels[2],
                                "__VS__",labels[3],".pdf",sep=""),sep="/"))
      layout(matrix(c(rep(1,5),2), 1, 6, byrow = TRUE))
      ##Set NAs as zero
      Pdiff.by[which(is.na(Pdiff.by))] <- 0
      jet.black.colors = colorRampPalette(
        c("#00007F", "blue", "#007FFF", "cyan","black","yellow", "#FF7F00", "red", "#7F0000"))
      
      image(Pdiff.by, xlab='',ylab='',
            col=jet.black.colors(100),
            axes=FALSE,
            main=paste("Rank Rank Hypergeometric Overlap Difference Map\n(",
                       labels[1],"-",labels[3],")-\n(",labels[2],"-",labels[3],")", sep=""),
            zlim=c(minhypermat, maxhypermat))
      
      color.bar(jet.black.colors(100), min=minhypermat,
                max=maxhypermat, nticks=6, title="-log(BY corrected P-value)")
      dev.off()
  }
  
  
  result<- list(hypermat1=hypermat1,
                hypermat2=hypermat2,
                Pdiff=Pdiff,
                Pdiff.by=Pdiff.by)
    
  return(result)
}