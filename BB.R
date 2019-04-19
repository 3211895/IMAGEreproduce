source('simulation_function.R')

load(file)  #######load the simulation data

genotype=t(geno$hap1+geno$hap2)  #######genotype data
ym=t(data$y)
rm=t(data$r)    ######methylation data
n=ncol(rm)  ###sample size
m=nrow(rm)  #####number of SNP-CpG pairs
summary<-data.frame()

for(iVar in 1:m)
{
  cat(iVar, '\n')
  y=rm[iVar,]
  x=ym[iVar,]
  c=genotype[iVar,]
  t1<-system.time(pvalue<-bbfit(x,y,c))
  
  summary <- rbind(summary, data.frame(Loc = iVar, 
                               pvalue = pvalue,
                               Method='BB',time=t1[3]))
}
saveRDS(summary, paste0('BB.rds'))