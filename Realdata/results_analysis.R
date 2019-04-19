#########Power########################

res=readRDS('./IMAGE.rds')
idx=which(res[,8]==FALSE)
res[idx,5]=1         ###########set the p-values for unconverged sites as 1


load('pnull.RData')  #######load permutation pvalues
pvalue=sort(pvalue)
fdr=FDR(pvalue,pnull,10)

power=numeric()
for(i in seq(0.01,0.1,0.01))
{
  power=c(power,length(which(fdr<i)))
}

###########Overlap###########################
res1=readRDS('./IMAGE1.rds') #IMAGE results for dataset1
res2=readRDS('./IMAGE2.rds') #IMAGE results for dataset2

idx=which(res1[,8]==FALSE)
res1[idx,5]=1 
idx=which(res2[,8]==FALSE)
res2[idx,5]=1 

res1<-res1[order(res1$pvalue),]
res2<-res2[order(res2$pvalue),]

overlap=numeric()
N=nrow(res1)
n=floor(N/100)
tab=seq(n,n*10,n)

count=1
for(i in tab)
{
  loc1=res1[1:i,1]
  loc2=res2[1:i,1]
  idxx=match(loc1,loc2)
  overlap[cout]=length(which(!is.na(idxx)))/length(idxx)
  count=count+1
}

#########enrich###################################

res=readRDS('./IMAGE.rds')
load('CpG_annotation.RData') #####annotation for CpG islands and other region


load('baboon.RData')


idx=which(res[,8]==FALSE)
res[idx,5]=1

p=res[,5]
p=sort(p)
idx=which(res[,5]<p[7044])   #power at 5% FDR
idxf=which(res[,5]>=p[7044])

t1=0
t2=0
t3=0
t4=0
for(i in 1:21)
{
  idx1=which(chr==i)
  idx2=intersect(idx1,idx)
  cpg=CpGlist[idx2]
  cpg=unique(cpg)##
  t3=t3+length(cpg)
  idxx=match(cpg,island[[i]]) ##### set CpG island  North shore South shore North shelf South shelf Open sea
  t1=t1+length(cpg)-length(which(is.na(idxx)))
  idx2=intersect(idx1,idxf)
  cpg=setdiff(CpGlist[idx2],CpGlist[intersect(idx1,idx)])
  cpg=unique(cpg)###
  t4=t4+length(cpg)
  idxx=match(cpg,island[[i]])
  t2=t2+length(cpg)-length(which(is.na(idxx)))
  
}

d=matrix(c(t1,t3,t2,t4),ncol=2)
fisher.test(d) ########fisher test for enrichment 

###########SNP disrupt CpG#######################
load('baboon.RData')
##disrupt CpG sites###
disrupt=numeric()
d=SNPlist-CpGlist
d=abs(d)
for(i in 1:length(d))
{
  if(d[i]==1)
  {
    disrupt=c(disrupt,i)
  }
}

res=readRDS('IMAGE.rds')
idx=which(res[,8]==FALSE)
res[idx,5]=1

p=res[,5]
p=sort(p)
idx=which(res[,5]<p[7044])   #power at 5% FDR
idxf=which(res[,5]>=p[7044])

t1=length(unique(intersect(disrupt,asm[idx,1])))
t2=length(unique(intersect(disrupt,asm[idxf,1])))

t3=0
t4=0
for(i in 1:21)
{
  idx1=which(chr==i)
  idx2=intersect(idx1,idx)
  cpg=CpGlist[idx2]
  cpg=unique(cpg)##
  t3=t3+length(cpg)
  idx2=intersect(idx1,idxf)
  cpg=setdiff(CpGlist[idx2],CpGlist[intersect(idx1,idx)])
  cpg=unique(cpg)###
  t4=t4+length(cpg)
}

d=matrix(c(t1,t3,t2,t4),ncol=2) 
fisher.test(d) ######fisher test for CpG sites disrupted by SNPs













