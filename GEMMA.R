load(file)  #######load the simulation data


##############write the Kinship matrix file##############
STR=character()
for(i in 1:n)
{
  str=character()
  for(j in 1:n)
  {
    str=paste(str,as.character(K[i,j]),sep='\t')
  }
  STR=paste(STR,str,sep='\n')
}
STR=sub('^..','',STR)
Kfile=paste0('K.txt',sep='')
writeLines(STR,con=Kfile,sep="\n")  

genotype=t(geno$hap1+geno$hap2) #######genotype data
ym=t(data$y)
rm=t(data$r) ######methylation data
n=ncol(rm)  ###sample size
m=nrow(rm)  #####number of SNP-CpG pairs

summary<-data.frame()

for(iVar in 1:m)
{
  cat(iVar, '\n')
  tmp1=ym[iVar,]
  tmp2=rm[iVar,]
  y=log2((tmp1+0.01)/(tmp2-tmp1+0.01)) ###normalized data in the form of M-values
  ########write the phenotype file#############################
  str=as.character(y[1])
  for(i in 2:n)
  {
    str=paste(str,as.character(y[i]),sep='\n')
  }
  yfile=paste0('yfile.txt')
  writeLines(str,con=yfile,sep="\n")
  
  #######write the genotype file##########
  x=genotype[iVar,]
  x=t(x)
  
  rs<-rep('rs',1)
  A<-c("A")
  B<-c("C")
  last=cbind(rs,A,B,x)
  gfile=paste0('gfile.txt')
  write.table(last,gfile,quote=F,sep=',',col.names=F,row.names = F)
  ###############GEMMA###########################
  ofile=paste0('gemma')
  str='./gemma -g ' #######path of gemma
  str=paste(str,gfile,' -p ',yfile,' -k ',Kfile,' -lmm 1 -o ',ofile,sep='')
  system(str)
  p=read.table(paste0('./output/',ofile,'.assoc.txt',sep=''))
  p=as.matrix(p)
  pvalue=as.numeric(p[2,11])
  summary <- rbind(summary, data.frame(Loc = iVar, 
                               pvalue = pvalue,Method='GEMMA'))
  
  str=paste('rm ',yfile,sep='')
  system(str)
  str=paste('rm ',gfile,sep='')
  system(str)
  str=paste('rm ./output/',ofile,'.assoc.txt',sep='')
  system(str)
  str=paste('rm ./output/',ofile,'.log.txt',sep='')
  system(str)
  
}

str=paste('rm ',Kfile,sep='')
system(str)
saveRDS(summary, paste0('GEMMA.rds'))
  
  
  