### R code from vignette source 'RRHO.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: RRHO.Rnw:76-85
###################################################
library(RRHO)
# Create "gene" lists:
  list.length <- 100
	list.names <- paste('Gene',1:list.length, sep='')
	gene.list1<- data.frame(list.names, sample(100))
	gene.list2<- data.frame(list.names, sample(100))
# Compute overlap and significance
	RRHO.example <-  RRHO(gene.list1, gene.list2, 
                        BY=TRUE, alternative='enrichment')


###################################################
### code chunk number 2: RRHO.Rnw:88-92
###################################################
# Examine Nominal (-log) pvalues
  image(RRHO.example$hypermat)
# Note: If lattice is available try: 
# levelplot(RRHO.example$hypermat)


###################################################
### code chunk number 3: RRHO.Rnw:95-103
###################################################
# FWER corrected pvalues using 50 random permutations:
pval.testing <- pvalRRHO(RRHO.example, 50)
pval.testing$pval
# The sampling distribution of the minimum 
# of the (-log) nominal p-values:
xs<- seq(0, 10, length=100)
plot(Vectorize(pval.testing$FUN.ecdf)(xs)~xs, 
     xlab='-log(pvalue)', ylab='ECDF', type='S')


###################################################
### code chunk number 4: RRHO.Rnw:108-112
###################################################
# Examine B-Y corrected pvalues 
# Note: probably nothing will be rejected in this 
# example as the data is generated from the null.
  image(RRHO.example$hypermat.by)


###################################################
### code chunk number 5: RRHO.Rnw:120-126
###################################################
m<- 100 ; n<- 100; k<- 50
data<- rhyper(1000, m, n, k)
pvals<- pmin(phyper(data,m,n,k, lower.tail=TRUE), 
             phyper(data,m,n,k, lower.tail=FALSE))
alpha<- 0.05
prop.table(table(pvals<alpha))


###################################################
### code chunk number 6: RRHO.Rnw:129-143
###################################################
getPval<- function(count,m,n,k){
  the.mean<- k*m/(m+n)
  if(count<the.mean){
    lower<- count
    upper<- 2*the.mean-count
  } else{
    lower<- 2*the.mean-count
    upper<- count
  }
  phyper(q=lower, m=m, n=n, k=k, lower.tail=TRUE) +
    phyper(q= upper, m=m, n=n, k=k, lower.tail=FALSE)                  
}
pvals<- sapply(data, getPval, m,n,k)
prop.table(table(pvals<alpha))


###################################################
### code chunk number 7: RRHO.Rnw:154-171
###################################################
size<- 500
list1<- data.frame(
  GeneIdentifier=paste('gen',1:size, sep=''), 
  RankingVal=-log(runif(size)))
list2<- data.frame(
  GeneIdentifier=paste('gen',1:size, sep=''), 
  RankingVal=-log(runif(size)))
list3<- data.frame(
  GeneIdentifier=paste('gen',1:size, sep=''), 
  RankingVal=-log(runif(size)))
rrho.comparison<- RRHOComparison(list1,list2,list3,
                stepsize=10,
                labels=c("phNPCs 1 wk PD vs 4 wk PD",
                         "My Progenitor vs Diff",
                         "In Vivo Stage 1 vs Stage 4"),
                plots=FALSE,
                outputdir=temp.dir);


###################################################
### code chunk number 8: RRHO.Rnw:174-176
###################################################
## The standard RRHO map between list1 and list 3.
image(rrho.comparison$hypermat1)  


###################################################
### code chunk number 9: RRHO.Rnw:178-181
###################################################
## The p-value of the difference between 
# (list1-list3)-(list2-list3).
image(rrho.comparison$Pdiff)


