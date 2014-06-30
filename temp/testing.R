#### Testing two sided hypotheses ####
n<- 112
sample1<- sample(n)
sample2<- sample(n)  
# The rank overlap of two independent samples should be about the mean for both one and two sided tests.
.test<- RRHO:::numericListOverlap(sample1, sample2, stepsize=10, alternative="two.sided")
dim(.test$log.pval)
library(lattice)
levelplot(.test$counts)
levelplot(-.test$log.pval)
levelplot(exp(-.test$log.pval))
.test$log.pval
max(exp(-.test$log.pval))
table(is.na(.test$log.pval))


list.length <- 100
list.names <- paste('Gene',1:list.length, sep='')
gene.list1<- data.frame(list.names, sample(list.length))
gene.list2<- data.frame(list.names, sample(list.length))
RRHO.example <-  RRHO(gene.list1, gene.list2, alternative='two.sided', BY = TRUE)
library(lattice)
library(grid)
levelplot(RRHO.example$hypermat)
levelplot(RRHO.example$hypermat.by)
RRHO.example <-  RRHO(alternative = 'enrichment',
       gene.list1,  gene.list2, 
       plots=TRUE, outputdir='/tmp/', labels=c("a","b"))

RRHO.example$call$alternative

RRHO.example <- RRHO(gene.list1, gene.list2, alternative='enrichment',
     plots=TRUE,outputdir='~/Downloads',labels=c('test1','test2'))


pvalRRHO(RRHO.example,1000) 


### Testing website:

data<- read.table('http://systems.crump.ucla.edu/rankrank/demo/rankrank.94aaa8b7716d7f.ranks_true.txt', skip=1)
head(data)
names(data)<- c('UniGene ID','Gene Name','Rank 1','Rank 2','Metric 1','Metric 2')



#### Testing RRHO comparisons ####

# options(stringsAsFactors=FALSE);
library(WGCNA);

## Generate synthetic data:
size<- 1000
list1<- data.frame(GeneIdentifier=paste('gen',1:size, sep=''), RankingVal=-log(runif(size)))
list2<- data.frame(GeneIdentifier=paste('gen',1:size, sep=''), RankingVal=-log(runif(size)))
list3<- data.frame(GeneIdentifier=paste('gen',1:size, sep=''), RankingVal=-log(runif(size)))

comparison.test<-RRHOComparison(list1,list2,list3,
                                stepsize=10, plots=FALSE,
                                labels=c("phNPCs 1 wk PD vs 4 wk PD",
                                         "My Progenitor vs Diff",
                                         "In Vivo Stage 1 vs Stage 4"))

levelplot(comparison.test$Pdiff)
levelplot(comparison.test$Pdiff.by)


(temp.dir<- tempdir())
RRHOComparison(list1,list2,list3,
               stepsize=10,
               labels=c("phNPCs 1 wk PD vs 4 wk PD",
                        "My Progenitor vs Diff",
                        "In Vivo Stage 1 vs Stage 4"),
               plots=TRUE,
               outputdir=temp.dir);


#### Use real data: #####
load(file='data/lists.RData')

##Run difference map
labelssss<- c("phNPCs 1 wk PD vs 8 wk PD","My Progenitor vs Diff","In Vivo Stage 1 vs Stage 6")
stepsize = 200;
RRHO.comp<- RRHOComparison(HNP, My, Sestan, stepsize = stepsize, plots=FALSE);
library(lattice)
levelplot(RRHO.comp$hypermat1)
levelplot(RRHO.comp$Pdiff)
levelplot(RRHO.comp$Pdiff.by)
image(RRHO.comp$Pdiff.by)
table(is.na(RRHO.comp$Pdiff.by))


(temp.dir<- tempdir())
RRHOComparison(HNP ,My,  Sestan, stepsize, plots=TRUE, labels = labelssss, outputdir = temp.dir);