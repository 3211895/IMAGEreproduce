

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
  ###########write the r file################
  str='loc'
  for(i in 1:n)
  {
    str1=paste('id',as.character(i),sep='')
    str=paste(str,str1,sep='\t')
  }
  STR=str
  str='csite_0'
  for(i in 1:n)
  {
    str=paste(str,as.character(rm[iVar,i]),sep='\t')
  }
  STR=paste(STR,str,sep='\n')
  rfile=paste0('rfile.txt')
  writeLines(STR,con=rfile,sep="\n")
  ###########write the y file################
  str='loc'
  for(i in 1:n)
  {
    str1=paste('id',as.character(i),sep='')
    str=paste(str,str1,sep='\t')
  }
  STR=str
  str='csite_0'
  for(i in 1:n)
  {
    str=paste(str,as.character(ym[iVar,i]),sep='\t')
  }
  STR=paste(STR,str,sep='\n')
  yfile=paste0('yfile.txt')
   writeLines(STR,con=yfile,sep="\n")
   
   ###########write the genotype file################
  tmp=genotype[iVar,]
  tmp[which(is.na(tmp))]=0
  str=as.character(tmp[1])
  for(i in 2:n)
  {
    str=paste(str,as.character(tmp[i]),sep='\n')
  }
  pfile=paste0('pred.txt')
  writeLines(str,con=pfile,sep="\n")
  
  ############MACAU##########################
  
  ofile=paste0('macau')
  
  str='./macau -g ' #######path of MACAU
  str=paste0(str,yfile,' -t ',rfile,' -p ',pfile,' -k ')
  str=paste0(str,Kfile,' -c cov_',n,'.txt -bmm -o ',ofile)
  system(str)
  p=read.table(paste0('./output/',ofile,'.assoc.txt',sep=''))
  p=as.matrix(p)
  pvalue=as.numeric(p[2,6])
  p=readLines(paste0('./output/',ofile,'.log.txt',sep=''))
  p=p[16]
  p=strsplit(p,' ')
  p=p[[1]][6]
  
  
  summary <- rbind(summary, data.frame(Loc = iVar, 
                               pvalue = pvalue,
                               Method='MACAU',time=as.numeric(p)*60))
  str=paste('rm ',rfile,sep='')
  system(str)
  str=paste('rm ',yfile,sep='')
  system(str)
  str=paste('rm ',pfile,sep='')
  system(str)
  str=paste('rm ./output/',ofile,'.assoc.txt',sep='')
  system(str)
  str=paste('rm ./output/',ofile,'.log.txt',sep='')
  system(str)
  
  
}

str=paste('rm ',Kfile,sep='')
system(str)
saveRDS(summary, paste0('MACAU.rds'))