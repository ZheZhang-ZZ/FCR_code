library(data.table)
library(Amelia)
set.seed(1234)
input<-fread("impute_input.csv")
input<-input[,c('个体号','测定日龄','当天体重','采食量')]
# input<-input[1:1000,]
# start_time <- Sys.time()
imp<-amelia(input,ts='测定日龄',cs='个体号',
            m=30,parallel = "multicore",
            ncpus=30,polytime=2,intercs=T,p2s=0)
# end_time <- Sys.time()  
#  
# end_time - start_time
save(imp,file="amelia_res.RData")
