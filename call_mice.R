library(data.table)
library(mice)
set.seed(1234)
input<-fread("impute_input.csv")
# input<-input[1:1000,]
# start_time <- Sys.time()
imp <- parlmice(input,m=20,
                cluster.seed = 1234,
                maxit=20,method='pmm',
                n.core=20,n.imp.core=1)
# end_time <- Sys.time()
# 
# end_time - start_time
save(imp,file="impute_res.RData")
