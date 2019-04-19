#############power#################

res=readRDS('./IMAGE.rds') ###results
sortres <- res[order(res$pvalue ), ] 
loc=sortres[,1]
power=cal_Power(loc,fdr)
pvalue=sortres[,5]

############FDR control####################
res=readRDS('./IMAGE.rds') ###results
sortres <- res[order(res$pvalue ), ] 
loc=sortres[,1]

counts = matrix(nrow=0, ncol=2)
i=1
if(loc[i]>1000)
{
  counts <- rbind(counts, c( 0, 1))
}else{
  counts <- rbind(counts, c( 1, 0))
}
for(i in 2:10000)
{
  if(loc[i]>1000)
  {
    counts <- rbind(counts, c(counts[i-1, 1], counts[i-1, 2]+1))
  }else{
    counts <- rbind(counts, c(counts[i-1, 1]+1, counts[i-1, 2]))
  }
}

FDR_true=numeric()
for(i in 1:10000)
{
  FDR_true[i]=counts[i,2]/i
}


load('pnull.RData')  #######load permutation pvalues

FDR_estimate=FDR(pvalue,pnull,10)

FDR_true=FDR_true[1:length(which(FDR_true<0.2))]
FDR_estimate=FDR_estimate[1:length(FDR_true)]



###########adjust for p-values through genomic control factor#################
adjust_pvalue<-function(pvalue)
{
chisq <- qchisq(1-pvalue,1)
lambda=median(chisq)/qchisq(0.5,1)
chisq=chisq/lambda
pvalue=pchisq( chisq, 1, lower.tail = F)

return(pvalue)

}

###########simulation power#####################
#fdr   power at an fdr of XXX
cal_Power<-function(loc,fdr)
{
  counts = matrix(nrow=0, ncol=2)
  i=1
  if(loc[i]>1000)
  {
    counts <- rbind(counts, c( 0, 1))
  }else{
    counts <- rbind(counts, c( 1, 0))
  }
  for(i in 2:10000)
  {
    if(loc[i]>1000)
    {
      counts <- rbind(counts, c(counts[i-1, 1], counts[i-1, 2]+1))
    }else{
      counts <- rbind(counts, c(counts[i-1, 1]+1, counts[i-1, 2]))
    }
  }
  # counts <- data.frame(
  #                      TPR=counts[, 1]/1000,
  #                       FDR=counts[, 2]/9000)
  TPR=counts[,1]/1000
  FDR=numeric()
  for(i in 1:10000)
  {
    FDR[i]=counts[i,2]/i
  }
  FPR=counts[,2]/9000
  #  idx=which(FDR>0.1)[1]
  idx=which(FDR==min(FDR[which(FDR>fdr)]))
  power=max(TPR[idx])
  return(power)
}

###########simulation/real data FDR##########################
#pmatrix  a vector of permutation p-values
#pvalue  a vector of p-values from the simulation data/real data
#num permutation times
FDR <- function(pvalue,pmatrix,num)
{
  FDR=numeric(length=length(pvalue))
  n=length(pvalue)
  b=pmatrix
  for(i in 1:n)
  {
    num=length(which(b<pvalue[i]))
    FDR[i]=num/(num*i)
  }
  return(FDR)
}