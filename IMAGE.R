library('IMAGE')

load(file) ##############load the simulation data


summary<-image(geno,data,K)

saveRDS(summary, paste0('IMAGE.rds'))