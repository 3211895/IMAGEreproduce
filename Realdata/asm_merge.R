name=read.table('fq.txt')    ############names of all individuals ####################
name=as.matrix(name)

for(j in 1:67)  ########67 baboons
{
  tmp=strsplit(name[j],'_')
  name[j]=(tmp[[1]][1])
}

for(i in 1:21) #####21 chromosomes, we call asm for each chromosome.
{
  chr=numeric()
  SNP=numeric()
  CpG=numeric()
  SNPlist=numeric()
  CpGlist=numeric()
  count=1
  for(j in 1:67)
  {
    
    str=paste(name[j],'_part',as.character(i),'.asm',sep='') 
    if(file.exists(str))
    {
      data=read.table(str)
      data=as.matrix(data)
      N=nrow(data)
    }else{
      next
    }
    for(k in 2:N)
    {
      tmp=strsplit(data[k,1],'_')
      chr[count]=390124480-as.numeric(tmp[[1]][2])+1    
      SNP[count]=as.numeric(data[k,2])
      CpG[count]=as.numeric(data[k,6])
      count=count+1
    }
  }
  sSNP=unique(sort(SNP))
  N=length(sSNP)
  for(k in 1:N)
  {
    idx=which(SNP==sSNP[k])
    CpGtmp=unique(sort(CpG[idx]))
    n1=length(CpGtmp)
    n2=length(SNPlist)
    SNPlist[(n2+1):(n2+n1)]=sSNP[k]
    CpGlist[(n2+1):(n2+n1)]=CpGtmp
  }
  r1=matrix(0,ncol=67,nrow=length(CpGlist))
  r2=r1
  y1=r1
  y2=r1
  genotype=r1
  for(j in 1:67)
  {
    
    cat(j)
    str=paste(name[j],'_part',as.character(i),'.asm',sep='')
    if(file.exists(str))
    {
      data=read.table(str)
      data=as.matrix(data)
      N=nrow(data)
    }else{
      next
    }
    for(k in 2:N)
    {
      snp=as.numeric(data[k,2])
      cpg=as.numeric(data[k,6])
      idx=which(SNPlist==snp)
      idx1=which(CpGlist[idx]==cpg)
      if(data[k,3]==data[k,4]&&data[k,3]==data[k,5])
      {
        tmp=strsplit(data[k,7],'-')
        y1[idx[idx1],j]=as.numeric(tmp[[1]][1])
        y2[idx[idx1],j]=0
        r1[idx[idx1],j]=as.numeric(tmp[[1]][1])+as.numeric(tmp[[1]][2])
        r2[idx[idx1],j]=0
        genotype[idx[idx1],j]=0
      }else if(data[k,3]!=data[k,4]&&data[k,3]!=data[k,5]){
        tmp=strsplit(data[k,7],'-')
        y1[idx[idx1],j]=as.numeric(tmp[[1]][1])
        y2[idx[idx1],j]=0
        r1[idx[idx1],j]=as.numeric(tmp[[1]][1])+as.numeric(tmp[[1]][2])
        r2[idx[idx1],j]=0
        genotype[idx[idx1],j]=2
      }else{
        tmp=strsplit(data[k,7],'-')
        tmp2=strsplit(data[k,8],'-')
        y1[idx[idx1],j]=as.numeric(tmp[[1]][1])
        y2[idx[idx1],j]=as.numeric(tmp2[[1]][1])
        r1[idx[idx1],j]=as.numeric(tmp[[1]][1])+as.numeric(tmp[[1]][2])
        r2[idx[idx1],j]=as.numeric(tmp2[[1]][1])+as.numeric(tmp2[[1]][2])
        genotype[idx[idx1],j]=1
      }
    }
    
    
  }
  chr=numeric(length=length(CpGlist))+i
  str2=paste('data_chr',as.character(i),'.RData',sep='')
  save(SNPlist,CpGlist,y1,y2,r1,r2,genotype,file=str2)
  
  
}



CpGList=numeric()
SNPList=numeric()
Chr=numeric()
r1m=matrix(0,ncol=67,nrow=0)
r2m=r1m
y1m=r1m
y2m=r1m
genotypem=r1m

for(i in 1:21)
{
  str=paste('data_chr',as.character(i),'.RData',sep='')
  load(str)
  Chr=c(Chr,chr)
  CpGList=c(CpGList,CpGlist)
  SNPList=c(SNPList,SNPlist)
  r1m=rbind(r1m,r1)
  r2m=rbind(r2m,r2)
  y1m=rbind(y1m,y1)
  y2m=rbind(y2m,y2)
  genotypem=rbind(genotypem,genotype)
}

chr=Chr
ym=y1m+y2m
rm=r1m+r2m
CpGlist=CpGList
SNPlist=SNPList

N=nrow(genotypem)

############Filtering i#####################
num=numeric()
for(i in 1:N)
{
  num[i]=length(which(rm[i,]>0))
}
idx=which(num<20)
y1m=y1m[-idx,]
y2m=y2m[-idx,]
r1m=r1m[-idx,]
r2m=r2m[-idx,]
genotypem=genotypem[-idx,]
ym=ym[-idx,]
rm=rm[-idx,]
CpGlist=CpGlist[-idx]
SNPlist=SNPlist[-idx]
chr=chr[-idx]
N=nrow(genotypem)
############Filtering ii#####################
num=numeric()
num1=numeric()
num2=numeric()
for(i in 1:N)
{
  idx=which(rm[i,]>0)
  ratio=ym[i,idx]/rm[i,idx]
  num[i]=length(idx)
  num1[i]=length(which(ratio<0.1))
  num2[i]=length(which(ratio>0.9))
}
ratio1=num1/num
ratio2=num2/num
idx1=which(ratio1>0.9)
idx2=which(ratio2>0.9)
idx=union(idx1,idx2)
y1m=y1m[-idx,]
y2m=y2m[-idx,]
r1m=r1m[-idx,]
r2m=r2m[-idx,]
genotypem=genotypem[-idx,]
ym=ym[-idx,]
rm=rm[-idx,]
CpGlist=CpGlist[-idx]
SNPlist=SNPlist[-idx]
chr=chr[-idx]
N=nrow(genotypem)

############Filtering iii#####################

a_rdp=numeric()
for(i in 1:N)
{
  a_rdp[i]=sum(rm[i,])/length(which(rm[i,]>0))
}
idx=which(a_rdp<5)
y1m=y1m[-idx,]
y2m=y2m[-idx,]
r1m=r1m[-idx,]
r2m=r2m[-idx,]
genotypem=genotypem[-idx,]
ym=ym[-idx,]
rm=rm[-idx,]
CpGlist=CpGlist[-idx]
SNPlist=SNPlist[-idx]
chr=chr[-idx]
N=nrow(genotypem)

############Filtering iv#####################
maf=numeric()
for(i in 1:N)
{
  idx=which(rm[i,]>0)
  maf[i]=(length(which(genotypem[i,idx]==1))+length(which(genotypem[i,idx]==2))*2)/length(idx)/2
}
idx=union(which(maf<0.05),which(maf>0.95))
y1m=y1m[-idx,]
y2m=y2m[-idx,]
r1m=r1m[-idx,]
r2m=r2m[-idx,]
genotypem=genotypem[-idx,]
ym=ym[-idx,]
rm=rm[-idx,]
CpGlist=CpGlist[-idx]
SNPlist=SNPlist[-idx]
chr=chr[-idx]
N=nrow(genotypem)

geno<-list()
geno[[1]]<-matrix(0,ncol=N,nrow=67)
geno[[2]]<-matrix(0,ncol=N,nrow=67)
for(i in 1:N)
{
  for(j in 1:67)
  {
    if(is.na(genotypem[i,j]))
    {
      geno[[2]][j,i]=0/0
      geno[[1]][j,i]=0/0
    }else if(genotypem[i,j]==1){
      geno[[2]][j,i]=1
    }else if(genotypem[i,j]==2){
      geno[[2]][j,i]=1
      geno[[1]][j,i]=1
    }
  }
}
names(geno) <- c('hap1', 'hap2')
data<-list()
data[[1]]<-matrix(0,ncol=N,nrow=67)
data[[2]]<-matrix(0,ncol=N,nrow=67)
data[[3]]<-matrix(0,ncol=N,nrow=67)
data[[4]]<-matrix(0,ncol=N,nrow=67)
data[[5]]<-matrix(0,ncol=N,nrow=67)
data[[6]]<-matrix(0,ncol=N,nrow=67)
data[[1]]<-t(rm)
data[[2]]<-t(ym)
data[[3]]<-t(r1m)
data[[4]]<-t(r2m)
data[[5]]<-t(y1m)
data[[6]]<-t(y2m)
names(data) <- c('r', 'y', 'r1', 'r2', 'y1', 'y2')




