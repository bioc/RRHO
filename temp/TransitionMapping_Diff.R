#Rank Rank Hypergeometric Overlap based on Plaisier et al., Nucleic Acids Research, 2010
#Compares two RRHO maps to determine significance
RRHO.Comparison <- function(list1,list2,list3,stepsize,labels,plots=FALSE,outputdir=NULL) {
  ## Significance testing of the difference between two RRHO maps.
  ## RRHO map 1: list1 vs list3
  ## RRHO map 2: list2 vs list3
  ##
  ## list 1 is a data.frame from experiment 1 with two columns, column 1 is the Gene Identifier, column 2 is the signed ranking value (e.g. signed -log10(p-value) or fold change)
  ## list 2 is a data.frame from experiment 2 with two columns, column 1 is the Gene Identifier, column 2 is the signed ranking value (e.g. signed -log10(p-value) or fold change)
  ## list 3 is a data.frame from experiment 3 with two columns, column 1 is the Gene Identifier, column 2 is the signed ranking value (e.g. signed -log10(p-value) or fold change).  
  ## stepsize indicates how many genes to increase by in each algorithm iteration
  
  if (length(list1[,1])!=length(unique(list1[,1]))) 
    stop('Non-unique gene identifier found in list1');
  if (length(list2[,1])!=length(unique(list2[,1]))) 
    stop('Non-unique gene identifier found in list2');
  if (length(list3[,1])!=length(unique(list3[,1])))
    stop('Non-unique gene identifier found in list3');
  
  list1 = list1[order(list1[,2],decreasing=TRUE),];
  list2 = list2[order(list2[,2],decreasing=TRUE),];
  list3 = list3[order(list3[,2],decreasing=TRUE),];
  nlist1 = length(list1[,1]);
  nlist2 = length(list2[,1]);
  nlist3 = length(list3[,1]);
  
  ## Number of genes on the array
  N = max(c(nlist1,nlist2,nlist3));
  
  hypermat1 = matrix(data=NA,nrow=length(seq(1,nlist1,stepsize)),ncol=length(seq(1,nlist3,stepsize)));
  OR1 = matrix(data=NA,nrow=length(seq(1,nlist1,stepsize)),ncol=length(seq(1,nlist3,stepsize)));
  SE1 = matrix(data=NA,nrow=length(seq(1,nlist1,stepsize)),ncol=length(seq(1,nlist3,stepsize)));
  countx = county = 0;
  ##Loop over the experiments
  for (i in seq(1,nlist1,stepsize)) {
    countx = countx + 1;
    for (j in seq(1,nlist3,stepsize)) {
      county = county + 1;
      ## Parameters for the hypergeometric test
      k = length(intersect(list1[1:i,1],list3[1:j,1]));
      s = length(list1[1:i,1]);
      M = length(list3[1:j,1]);
      ## Hypergeometric test converting to log10
      ## log10(x) = log_e(x) * log10(e)
      hypermat1[countx,county] = -phyper(k-1,M,N-M,s,lower.tail=FALSE,log.p=TRUE) * log10(exp(1));
      contingencytable = rbind(c(k,M-k),c(s-k,N-s-M+k));
      OR1[countx,county] = (contingencytable[1,1]*contingencytable[2,2])/(contingencytable[1,2]*contingencytable[2,1]);
      SE1[countx,county] = sqrt(1/contingencytable[1,1] + 1/contingencytable[1,2] + 1/contingencytable[2,1] + 1/contingencytable[2,2]);
    }
    county=0;
  }

  hypermat2 = matrix(data=NA,nrow=length(seq(1,nlist2,stepsize)),ncol=length(seq(1,nlist3,stepsize)));
  OR2 = matrix(data=NA,nrow=length(seq(1,nlist2,stepsize)),ncol=length(seq(1,nlist3,stepsize)));
  SE2 = matrix(data=NA,nrow=length(seq(1,nlist2,stepsize)),ncol=length(seq(1,nlist3,stepsize)));
  countx = county = 0;
  ##Loop over the experiments
  for (i in seq(1,nlist2,stepsize)) {
    countx = countx + 1;
    for (j in seq(1,nlist3,stepsize)) {
      county = county + 1;
      ## Parameters for the hypergeometric test
      k = length(intersect(list2[1:i,1],list3[1:j,1]));
      s = length(list2[1:i,1]);
      M = length(list3[1:j,1]);
      ## Hypergeometric test converting to log10
      ## log10(x) = log_e(x) * log10(e)
      hypermat2[countx,county] = -phyper(k-1,M,N-M,s,lower.tail=FALSE,log.p=TRUE) * log10(exp(1));
      contingencytable = rbind(c(k,M-k),c(s-k,N-s-M+k));
      OR2[countx,county] = (contingencytable[1,1]*contingencytable[2,2])/(contingencytable[1,2]*contingencytable[2,1]);
      SE2[countx,county] = sqrt(1/contingencytable[1,1] + 1/contingencytable[1,2] + 1/contingencytable[2,1] + 1/contingencytable[2,2]);
    }
    county=0;
  }

  ## Normal approximation to the odds ration
  ## Z = log(OR1)-log(OR2)/sqrt(SE1^2+SE2^2)
  ORdiff = log(OR1/OR2);
  SEdiff = sqrt(SE1^2 + SE2^2);  
  Zdiff = ORdiff/SEdiff;
  ## Using a two tailed P-value, need to use log p-values
  ##  Pdiff = (1-pnorm(abs(Zdiff)))*2;
  Pdiff = (-1*pnorm(abs(Zdiff),log.p=TRUE,lower.tail=FALSE)-log(2))*log10(exp(1));
  
  ## Convert hypermat to a vector and Benjamini Yekutieli FDR correct
  Pdiffvec = matrix(Pdiff,nrow=nrow(Pdiff)*ncol(Pdiff),ncol=1);
  Zdiffvec = matrix(Zdiff,nrow=nrow(Zdiff)*ncol(Zdiff),ncol=1);
  nonaind = which(!is.na(Pdiffvec));
  Pdiff.byvec = matrix(data=NA,nrow=length(Pdiffvec),ncol=1);
  Pdiff.byvec[nonaind] = sign(Zdiffvec)[nonaind] * -log10(p.adjust(10^-Pdiffvec[nonaind],method="BY"));
  Pdiff.by = matrix(Pdiff.byvec,nrow=nrow(Pdiff),ncol=ncol(Pdiff));
  
  if (plots) {
    
    ## Function to plot color bar
    ## Modified from http://www.colbyimaging.com/wiki/statistics/color-bars
    color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
      scale = (length(lut)-1)/(max-min);  
      plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='');
      mtext(title,2,2.3,cex=0.8);
      axis(2, round(ticks,0), las=1,cex.lab=0.8);
      for (i in 1:(length(lut)-1)) {
        y = (i-1)/scale + min;
        rect(0,y,10,y+1/scale, col=lut[i], border=NA);
      }
    }
    
    pdf(paste(outputdir,paste("RRHOMaps_Diff",labels[1],"__VS__",labels[3],".pdf",sep=""),sep="/"));
    #png(paste(outputdir,paste("RRHOMaps_Diff",labels[1],"__VS__",labels[3],".png",sep=""),sep="/"),width=10,height=10,units="in",res=300);
    jet.colors = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"));
    layout(matrix(c(rep(1,5),2), 1, 6, byrow = TRUE));
    image(hypermat1,xlab='',ylab='',col=jet.colors(100),axes=FALSE,main="Rank Rank Hypergeometric Overlap Map");
    mtext(labels[3],2,0.5);
    mtext(labels[1],1,0.5);
    color.bar(jet.colors(100),min=min(hypermat1,na.rm=TRUE),max=max(hypermat1,na.rm=TRUE),nticks=6,title="-log10(Nominal P-value)");
    dev.off();

    pdf(paste(outputdir,paste("RRHOMaps_Diff",labels[2],"__VS__",labels[3],".pdf",sep=""),sep="/"));
    #png(paste(outputdir,paste("RRHOMaps_Diff",labels[2],"__VS__",labels[3],".png",sep=""),sep="/"),width=10,height=10,units="in",res=300);
    layout(matrix(c(rep(1,5),2), 1, 6, byrow = TRUE));
    image(hypermat2,xlab='',ylab='',col=jet.colors(100),axes=FALSE,main="Rank Rank Hypergeometric Overlap Map");
    mtext(labels[3],2,0.5);
    mtext(labels[2],1,0.5);
    color.bar(jet.colors(100),min=min(hypermat2,na.rm=TRUE),max=max(hypermat2,na.rm=TRUE),nticks=6,title="-log10(Nominal P-value)");
    dev.off();

    ##Set max and mins of colorbar for in vitro systems used in paper
    minhypermat = -46;
    maxhypermat = 46;
    pdf(paste(outputdir,paste("RRHOMaps_DiffMap",labels[1],"__VS__",labels[2],"__VS__",labels[3],".pdf",sep=""),sep="/"));
    #png(paste(outputdir,paste("RRHOMaps_DiffMap",labels[1],"__VS__",labels[2],"__VS__",labels[3],".png",sep=""),sep="/"),width=10,height=10,units="in",res=300);
    layout(matrix(c(rep(1,5),2), 1, 6, byrow = TRUE));
    ##Set NAs as zero
    Pdiff.by[which(is.na(Pdiff.by))] = 0;
    jet.black.colors = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","black","yellow", "#FF7F00", "red", "#7F0000"));
    image(Pdiff.by,xlab='',ylab='',col=jet.black.colors(100),axes=FALSE,main=paste("Rank Rank Hypergeometric Overlap Difference Map\n",labels[1],"VS",labels[3],"\n",labels[2],"VS",labels[3]),zlim=c(minhypermat,maxhypermat));
    color.bar(jet.black.colors(100),min=minhypermat,max=maxhypermat,nticks=6,title="-log10(BY corrected P-value)");
    dev.off();
  }
  return(list(hypermat1=hypermat1,hypermat2=hypermat2,Pdiff=Pdiff,Pdiff.by=Pdiff.by));
}


##Get the file names as arguments from the command line
args = commandArgs(trailingOnly=T);
origoutputdir = args[1];

rootdir=".";
## output directory
outputdir = paste(origoutputdir,"TransitionMapping",sep="/");
HNPSestandir = paste(rootdir,"SestanDE",sep="/");

options(stringsAsFactors=FALSE);
library(WGCNA);

cat('Comparing HNPs 1 wk vs 4 wk VS In vivo Stage 1 vs Stage 4 to your data progenitor vs differentiated VS In vivo Stage 1 vs Stage 4\n');
## Read in all files
HNPDE = read.csv(paste(HNPSestandir,paste('HNPs_DE_1vs4.csv',sep=""),sep="/"));
SestanDE = read.csv(paste(HNPSestandir,paste('Sestan_DE_1vs4.csv',sep=""),sep="/"));
MyDE = read.csv(paste(outputdir,'YourDatavsInVivo.csv',sep="/"));

## Format lists for RRHO
HNP = data.frame(GeneIdentifier=HNPDE$ID, RankingVal=sign(as.numeric(HNPDE$Beta))*-log10(as.numeric(HNPDE$Pval)));
Sestan = data.frame(GeneIdentifier=SestanDE$ID, RankingVal=sign(as.numeric(SestanDE$Beta))*-log10(as.numeric(SestanDE$Pval)));
MyDE = data.frame(GeneIdentifier=MyDE$ID, RankingVal=sign(as.numeric(MyDE$Beta))*-log10(as.numeric(MyDE$Pval)));

## Use only genes found in all datasets
allgenes = unique(intersect(intersect(HNP$GeneIdentifier,Sestan$GeneIdentifier),MyDE$GeneIdentifier));
HNPmatchind = match(allgenes,HNP$GeneIdentifier);
HNP = HNP[HNPmatchind,];
Sestanmatchind = match(allgenes,Sestan$GeneIdentifier);
Sestan = Sestan[Sestanmatchind,];
Mymatchind = match(allgenes,My$GeneIdentifier);
My = My[Mymatchind,];

#save(Sestan,My, HNP, file='data/lists.RData')


##################################################
##Run difference map
stepsize = 200;
RRHO.Comparison(HNP,My,Sestan,stepsize,c("phNPCs 1 wk PD vs 4 wk PD","My Progenitor vs Diff","In Vivo Stage 1 vs Stage 4"),TRUE,outputdir);


cat('Comparing HNPs 1 wk vs 8 wk VS In vivo Stage 1 vs Stage 6 to your data progenitor vs differentiated VS In vivo Stage 1 vs Stage 6\n');
## Read in all files
HNPDE = read.csv(paste(HNPSestandir,paste('HNPs_DE_1vs8.csv',sep=""),sep="/"));
SestanDE = read.csv(paste(HNPSestandir,paste('Sestan_DE_1vs6.csv',sep=""),sep="/"));
MyDE = read.csv(paste(outputdir,'YourDatavsInVivo.csv',sep="/"));

dim(HNPDE)
dim(SestanDE)
dim(MyDE)

## Format lists for RRHO
HNP = data.frame(GeneIdentifier=HNPDE$ID, RankingVal=sign(as.numeric(HNPDE$Beta))*-log10(as.numeric(HNPDE$Pval)));
Sestan = data.frame(GeneIdentifier=SestanDE$ID, RankingVal=sign(as.numeric(SestanDE$Beta))*-log10(as.numeric(SestanDE$Pval)));
My = data.frame(GeneIdentifier=MyDE$ID, RankingVal=sign(as.numeric(MyDE$Beta))*-log10(as.numeric(MyDE$Pval)));
#My<- MyDE

dim(HNP) ;dim(Sestan); dim(My)

## Use only genes found in all datasets
allgenes = unique(intersect(intersect(HNP$GeneIdentifier,Sestan$GeneIdentifier),My$GeneIdentifier));
HNPmatchind = match(allgenes,HNP$GeneIdentifier);
HNP = HNP[HNPmatchind,];
Sestanmatchind = match(allgenes,Sestan$GeneIdentifier);
Sestan = Sestan[Sestanmatchind,];
Mymatchind = match(allgenes,My$GeneIdentifier);
My = My[Mymatchind,];

dim(HNP); dim(Sestan); dim(My)


##################################################
##Run difference map
stepsize = 200;
RRHO.Comparison(HNP,My,Sestan,stepsize,c("phNPCs 1 wk PD vs 8 wk PD","My Progenitor vs Diff","In Vivo Stage 1 vs Stage 6"),TRUE,outputdir);








