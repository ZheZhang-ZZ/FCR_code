library(BGLR)
library(data.table)

input<-fread("inputForBGLR.csv")
nIter=12000
burnIn=2000
saveAt='CORR'

ETA<-list( list(~平均日增重+
                  当天体重+
                  factor(品种品系)+
                  factor(测定站)+
                  factor(性别)+
                  factor(同期组)+
                  etp2+etp3+etp5+
                  etp6+etp8+etp9+
                  etp10+etp11+etp12+
                  etp13+etp14+etp15+
                  etp16+otd2+otd6+otd8+
                  otd9+otd10+otd11+otd12+
                  otd13+otd14+fid5+fid15+
                  fid16,data=input,model='FIXED'),
           list(~factor(个体号), data=input, model='BRR')
)

fm<-BGLR(y=input$采食量,ETA=ETA,nIter=nIter, burnIn=burnIn,saveAt=saveAt,verbose = F)
save(fm,file='fm.RData')
