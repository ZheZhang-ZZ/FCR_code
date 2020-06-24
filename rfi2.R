library(ggplot2)
library(data.table)
library(dplyr)
library(flexmix)
library(MASS)
library(BGLR)
library(orthopolynom)
#library(mice)
par(family='STSong')

setwd("/Users/zhezhang/Desktop/天邦数据/FI/guigang")
ped<-list()
ped_year<-c(2019,2020)
for(i in 1:length(ped_year)){
  ped[[i]]<-read.csv(paste0("ped_",ped_year[i],".csv"),h=T,fill=T)
}
ped<-rbindlist(ped)
ped<-ped[!duplicated(ped$个体号),]

### ped中的大耳号（场内编号）有重复的情况，比如1234出现两次，或者1234和DD.1234
### 这种情况就没法把FIRE产生的个体一一对应到系谱里
### 所以先把所有能对上的都对上，然后看重复对应的ID里哪个距离119日龄最近，就保留

## 第一步，去除大耳号前面的DD/YY/LL/PP
ear<-strsplit(as.character(ped$场内编号),"[.]")
ear<-sapply(ear,function(x)x[length(x)])
#ear<-unlist(ear)
#ear<-ped$场内编号

## 读入数据
data_dir<-"测定站发送美国CG数据/"
files <- list.files(data_dir)

data<-list()
for (file in 1:length(files)){
  data[[file]]<-read.csv(paste0(data_dir,files[file]),h=T,fill=T,
                         col.names=c("测定站","用户自编号","射频耳标",
                                     "进入时间","退出时间","采食时间",
                                     "采食量","当天体重","进入时料重",
                                     "退出时料重","加料料重"),
                         encoding="UTF-8")
  
}

data1<-rbindlist(data)
data1<-unique(data1)
nrow(data1) #4797478
nrow(data1[data1$射频耳标==0,]) #射频耳标=0，意味着猪的耳号识别不出来，这样的有293199条记录

## 将在系谱里的个体数据拎出来

a<-ear[ear%in%as.character(data1$用户自编号)]
length(a) #4559
anyDuplicated(a) # 1181，意味着系谱能对到数据的个体的大耳号有多达1181个重复，比如111071，在ped2019中出现的非常诡异，能对上3头猪
a<-unlist(a)
data2<-droplevels(data1[as.character(data1$用户自编号)%in%ear,]) #4327722

## 时间整理成标准格式
time1 = as.POSIXct(as.character(data2$进入时间),format = "%Y/%m/%d %H:%M")
time2 = as.POSIXct(as.character(data2$退出时间),format = "%Y/%m/%d %H:%M")
data2$进入时间<-time1
data2$退出时间<-time2

## 将数据按照个体号与进入测定站的时间排序
ord<-order(data2$用户自编号,data2$进入时间,decreasing=c(FALSE,FALSE))
data2<-data2[ord,]

## 将个体第一次进入测定站的时间单独拎出来(有一个个体202249，居然在数据里面出现两次，所以必须依靠测定站信息)
onset_time<-data2[!duplicated(data.frame(data2$用户自编号,data2$测定站)),c('用户自编号','测定站','进入时间')]
nrow(onset_time) #3493 这是能匹配到系谱里的个体数，但是由于系谱里有重复ID，导致出现数据里的个体对多个系谱里的ID的情况

##再把数据分为有重复个体的部分和没有重复个体的部分 注意目前为止，
##也只有202249这一个个体出现了在数据中代表两个个体的情况
dup_id<-onset_time[duplicated(onset_time$用户自编号),]
dup_time<-droplevels(onset_time[onset_time$用户自编号%in%dup_id$用户自编号,])
onset_time<-droplevels(onset_time[!onset_time$用户自编号%in%dup_id$用户自编号,])

## 将无重复数据与耳号数据相互进行匹配
b<-match(ear,as.character(onset_time$用户自编号))
onset_time2<-onset_time[b[!is.na(b)],]
ped2<-ped[!is.na(b),]

## 将有重复数据与耳号数据进行匹配
pedInDup<-which(ear%in%dup_time$用户自编号)
dupInPed<-which(dup_time$用户自编号%in%ear)
grid_table<-expand.grid(pedInDup,dupInPed)
grid_table<-grid_table[which(as.character(ear[grid_table$Var1])==
                               as.character(dup_time$用户自编号[grid_table$Var2])),]
dup_time2<-dup_time[grid_table$Var2,]
ped3<-ped[grid_table$Var1,]

## 再把dup数据与非dup数据合并，ped2与ped3合并
onset_time2<-rbind(onset_time2,dup_time2)
ped2<-rbind(ped2,ped3)

## 把进入测定站的日期单独拎出来
date1<-strsplit(as.character(onset_time2$进入时间),split=" ")
date1<-sapply(date1,function(x)x[1])

## 给这些数据加上一些必要的信息
onset_time2$进入日期 <- as.Date(date1,format="%Y-%m-%d")
onset_time2$出生日期 <- as.Date(ped2$出生日期,format="%Y/%m/%d")
onset_time2$品种品系 <- ped2$品种品系
onset_time2$ID <- ped2$个体号

## 计算开测日龄
#onset_time2<-onset_time2[order(onset_time2$用户自编号),]
onset_time2$开测日龄<-onset_time2$进入日期-onset_time2$出生日期

## 对于此时有重复ID的个体，计算其开测日龄距离119的距离,把大于10天的过滤掉
dup_id<-unique(onset_time2$用户自编号[duplicated(onset_time2$用户自编号)])

#把数据分为不重复（onset_time3）与重复两部分
onset_time3<-droplevels(onset_time2[(!onset_time2$用户自编号%in%dup_id),]) # 2438
png(file="距离119日龄的距离.png",width=8, height=5,units='in',res=300)
par(family='STSong')
barplot(table(abs(onset_time3$开测日龄-119)),xlab = "距离119日龄的距离(天)", ylab = "频数") ##看一下这组数据的开测日龄距离119的距离
                                              ##发现在10-11天左右频数就很低了，不妨就把阈值设为10
dev.off()

png(file="距离119日龄的距离2.png",width=8, height=5,units='in',res=300)
par(family='STSong')
barplot(table(onset_time3$开测日龄-119),xlab = "距离119日龄的距离(天)", ylab = "频数") ##看一下这组数据的开测日龄距离119的距离
dev.off()

onset_time3<-droplevels(onset_time3[(abs(119-onset_time3$开测日龄))<=10,]) #2190

dup_time<-droplevels(onset_time2[(onset_time2$用户自编号%in%dup_id),]) #2119
# dup_temp<-dup_time %>% group_by(用户自编号) %>% 
#   summarise(
#     diff = min(abs(119-开测日龄))
#   )

dup_time<-droplevels(dup_time[(abs(119-dup_time$开测日龄))<=10,])

## 此时重复部分仍然还有重复的个体，只能过滤这些个体
dup_id<-unique(data.frame(dup_time$用户自编号,dup_time$测定站)[duplicated(data.frame(dup_time$用户自编号,dup_time$测定站)),])
onset_time4<-dup_time[!dup_time$用户自编号%in%dup_id[,1],] #1030

## 将过滤掉剩下的两部分数据合并
onset_time_final<-rbind(onset_time3,onset_time4) #3220
anyDuplicated(data.frame(onset_time_final$用户自编号,onset_time_final$测定站)) # 0 说明这时候已经没有重复的了

## 将剩下的数据进行整理
a<-paste(data2$用户自编号,data2$测定站,sep=":")
b<-paste(onset_time_final$用户自编号,onset_time_final$测定站,sep=":")
matchID<-match(a,b)
data2$ID <- onset_time_final[matchID,'ID']
data2$开测日龄 <- onset_time_final[matchID,'开测日龄']
data3<-data2[!is.na(data2$ID)] #4003855

id<-droplevels(ped[match(data3$ID,ped$个体号),'个体号'])
sex<-ped[match(data3$ID,ped$个体号),'性别']
line<-ped[match(data3$ID,ped$个体号),'品种品系']
farm<-ped[match(data3$ID,ped$个体号),'本地猪出生猪场']
birth<-ped[match(data3$ID,ped$个体号),'出生日期']
time1 = as.POSIXct(as.character(data3$进入时间),format = "%Y-%m-%d %H:%M")
time2 = as.POSIXct(as.character(data3$退出时间),format = "%Y-%m-%d %H:%M")


data<-data.frame(个体号=id,耳号=data3$用户自编号,性别=sex,
                    出生日期=birth,出生场=farm,品种品系=line,
                    测定站=data3$测定站,开测日龄=as.numeric(data3$开测日龄),
                    进入时间=time1,退出时间=time2,
                    采食时间=data3$采食时间,进入料重=data3$进入时料重,
                    退出料重=data3$退出时料重,采食量=data3$采食量,
                    当天体重=data3$当天体重,喂料=data3$加料料重)

## 进入时间有几条NA数据 根据退出时间与采食时间退出其进入时间
#先把采食时间换算成秒
seconds<-strsplit(as.character(data$采食时间),split=":")
seconds<-sapply(seconds,function(x){
  as.numeric(x[1])*3600+
    as.numeric(x[2])*60+
    as.numeric(x[3])
})
whichNA<-which(is.na(data$进入时间))
data[whichNA,'进入时间']<-data[whichNA,'退出时间']-seconds[whichNA]

ord<-order(data$个体号,data$进入时间,decreasing=c(FALSE,FALSE))
data<-data[ord,]

#要把没有进入时间且没有退出时间的数据过滤掉
data<-droplevels(data[!(is.na(data$进入时间)&is.na(data$退出时间)),])

# 看一下各个测定日龄的测定次数分布
date_in<-strsplit(as.character(data$进入时间),split=" ")
date_in<-sapply(date_in,function(x)x[1]) #进入测定站的日期
date_in <- as.Date(date_in,format="%Y-%m-%d")
birth <- as.Date(data$出生日期,format="%Y/%m/%d")
days<-date_in-birth
data$测定日龄<-days ##4003852

test_frq<-data %>% 
  group_by(个体号,测定日龄) %>%
  summarise(频数=n())

idx<-by(test_frq,test_frq$个体号,function(x)return(1:nrow(x)))
idx<-unlist(idx,use.names=F)
test_frq$起测后x天<-idx

test_frq1<-test_frq %>% 
  group_by(起测后X天) %>%
  summarise(平均测定次数=mean(频数))

png(file="正数测定次数.png",width=8,height=5,units='in',res=300)
par(family='STSong')
barplot(平均测定次数~起测后X天,data=test_frq1)
dev.off()

idx<-by(test_frq,test_frq$个体号,function(x)return(nrow(x):1))
idx<-unlist(idx,use.names=F)
test_frq$结测前x天<-idx

test_frq2<-test_frq %>% 
  group_by(结测前x天) %>%
  summarise(平均测定次数=mean(频数))

png(file="倒数测定次数.png",width=8,height=5,units='in',res=300)
par(family='STSong')
barplot(平均测定次数~结测前X天,data=test_frq2)
dev.off()
## 通过作图我们发现，应该去除每个个体前两天以及最后一天的测定数据

front<-test_frq %>% group_by(个体号) %>% group_map(~ head(.x, 2L),keep=T)
behind<-test_frq %>% group_by(个体号) %>% group_map(~ tail(.x, 1L),keep=T)
front<-rbindlist(front)
behind<-rbindlist(behind)

front_id<-paste(front$个体号,front$测定日龄,sep=":")
behind_id<-paste(behind$个体号,behind$测定日龄,sep=":")
data_id<-paste(data$个体号,data$测定日龄,sep=":")

data<-droplevels(data[!data_id%in%c(front_id,behind_id),])
data$个体号<-factor(data$个体号,levels=as.character(unique(data$个体号)))

#write.csv(data,file="final_data.csv",row.names = F,quote=F,fileEncoding = "utf-8")

#################
### 数据编辑

#重新计算采食时长
seconds<-strsplit(as.character(data$采食时间),split=":")
seconds<-sapply(seconds,function(x){
  as.numeric(x[1])*3600+
    as.numeric(x[2])*60+
    as.numeric(x[3])
})

##按照Casey et al. 2005的方法，先计算7个统计量
#采食量
fiv<-data$采食量*1000

#采食时间
otv<-seconds

#采食效率
frv<-data$采食量 * 1000 / (otv / 60) # 采食效率 g/min

# 下一次进入测定站的料重与本次离开测定站的料重差
get_lwd<-function(x){
  t1<-x$退出料重
  t2<-x$进入料重
  t<-rep(NA,nrow(x))
  t[1:(length(t)-1)]<-(t2[-1]-head(t1,-1))*1000
  return(t)
}
lwd<-by(data,data$个体号,FUN=get_lwd)
lwd<-unlist(lwd, use.names=FALSE)

# 本次进入测定站的料重与上一次离开测定站的料重差
get_fwd<-function(x){
  t1<-x$退出料重
  t2<-x$进入料重
  t<-rep(NA,nrow(x))
  t[2:length(t)]<-(t2[-1]-head(t1,-1))*1000
  return(t)
}
fwd<-by(data,data$个体号,FUN=get_fwd)
fwd<-unlist(fwd, use.names=FALSE)


# 下一次进入测定站的时间与本次离开测定站的时间差
get_ltd<-function(x){
  t1<-x$退出时间
  t2<-x$进入时间
  t<-rep(NA,nrow(x))
  t[1:(length(t)-1)]<-as.numeric(difftime(t2[-1],head(t1,-1),units='secs'))
  return(t)
}
ltd<-by(data,data$个体号,FUN=get_ltd)
ltd<-unlist(ltd, use.names=FALSE)

# 本次进入测定站的时间与上一次离开测定站的时间差
get_ftd<-function(x){
  time1<-x$退出时间
  time2<-x$进入时间
  t<-rep(NA,nrow(x))
  t[2:length(t)]<-as.numeric(difftime(time2[-1],head(time1,-1),units='secs'))
  return(t)
}
ftd<-by(data,data$个体号,FUN=get_ftd)
ftd<-unlist(ftd, use.names=FALSE)

## 看看这几个指标的分布
#bin function
get_bin<-function(x,width){
  bin<-seq(min(x,na.rm=T),max(x,na.rm=T)+width,by=width)
  return(bin)
}
#箱线图函数
box<-function(x){
  ret<-list()
  stat<-boxplot(x,plot=F,na.rm=T)
  remov<-stat$out
  ret$stat<-stat
  ret$remov<-remov
  return(ret)
}

#fiv
et1<-which(fiv < -20)
fiv_remov_id<-et1
#(upper<-quantile(fiv, 3/4,na.rm=T)+1.5*IQR(fiv,na.rm=T)) #415.5 上界
upper<-2000 # 还是按照Casey的标准
et2<-which(fiv>upper)
fiv_remov_id<-append(fiv_remov_id,et2)
fiv2<-fiv[fiv<=upper]
bin<-get_bin(fiv2,5)
freq_fiv2=hist(fiv2, breaks=bin, include.lowest=TRUE, plot=TRUE,xlab="次采食量",ylab="频数")

fiv0<-fiv[otv==0] # 按照Casey的标准，fiv0不能大于20g
et3<-which(fiv%in%fiv0[fiv0>=20] & otv==0)
fiv_remov_id<-append(fiv_remov_id,et3)

#otv
et4<-which(otv < 0)
(upper<-quantile(otv, 3/4,na.rm=T)+1.5*IQR(otv,na.rm=T)) #423 上界
otv_remov_id<-which(otv>upper)
et5<-otv_remov_id
otv2<-otv[otv<=upper]
bin<-get_bin(otv2,5)
freq_otv2=hist(otv2, breaks=bin, include.lowest=TRUE, plot=TRUE)

#frv
et6<-which(fiv>0 & fiv<50 & frv!=Inf & frv>500) # 当FIV在（0，50）之间时，FRV>500的要过滤
frv_remov_id<-et6
et7<-integer(0) #按照et7的定义，我们的数据中没有符合这种情况的
frv2<-frv[which(fiv>=50 & frv!=Inf)] #注意frv是Inf的点并不一定要过滤，在这一步要保留,这一步算的是fiv>=50时的frv分布
#确定上界
(upper<-quantile(frv2, 3/4)+1.5*IQR(frv2)) #95.22206
et8<-which(fiv>=50 & frv!=Inf & frv > upper)
frv2<-frv2[frv2<=upper]
bin<-get_bin(frv2[frv2!=Inf],5)
freq_frv2=hist(frv2[frv2!=Inf], breaks=bin, include.lowest=TRUE, plot=TRUE)
#frv0, 当FRV=0时，要看一下otv的分布
otv_frv<-otv[frv==0]
(upper<-quantile(otv_frv, 3/4)+1.5*IQR(otv_frv)) #94.5 上界
otv_frv2<-otv_frv[otv_frv<=upper]
bin<-get_bin(otv_frv2,5)
freq_otv_frv2=hist(otv_frv2, breaks=bin, include.lowest=TRUE, plot=TRUE)
et9<-which(frv==0 & otv>upper)
frv_remov_id<-append(frv_remov_id,et9)
#frv_lo 当FRV!=0时，看一下frv的分布,最后还是决定用文献的值小于2g/min,这个吃饭速度太慢了
et10<-which(frv!=0 & frv<=2)
frv_remov_id<-c(frv_remov_id,et7,et8,et9,et10)

#lwd
(lower<-quantile(lwd, 1/4,na.rm=T)-1.5*IQR(lwd,na.rm=T)) #-386.2 下界
(upper<-quantile(lwd, 3/4,na.rm=T)+1.5*IQR(lwd,na.rm=T)) #643.4 上界
et11<-which(lwd<lower)
et12<-which(lwd>upper)
lwd2<-lwd[lwd>lower & lwd<upper]
bin<-get_bin(lwd2,5)
freq_lwd2=hist(lwd2, breaks=bin, include.lowest=TRUE, plot=TRUE)
lwd_remov_id<-c(et11,et12)

#fwd
(lower<-quantile(fwd, 1/4,na.rm=T)-1.5*IQR(fwd,na.rm=T)) #-386.2  下界
(upper<-quantile(fwd, 3/4,na.rm=T)+1.5*IQR(fwd,na.rm=T)) #643.4 上界
et13<-which(fwd<lower)
et14<-which(fwd>upper)
fwd2<-fwd[fwd>lower & fwd<upper]
bin<-get_bin(fwd2,5)
freq_fwd2=hist(fwd2, breaks=bin, include.lowest=TRUE, plot=TRUE)
fwd_remov_id<-c(et13,et14)

#ltd与ftd
ltd_remov_id<-which(ltd<0 & !is.na(ltd))
et15<-ltd_remov_id
ftd_remov_id<-which(ftd<0 & !is.na(ftd))
et16<-ftd_remov_id

#把所有这些错误项的id进行合并
remov_id<-unique(c(et1,et2,et3,et4,
                   et5,et6,et7,et8,
                   et9,et10,et11,et12,
                   et13,et14,et15,et16))
length(remov_id)/nrow(data) #0.06928202

## flexmix + robust regression 基于所有数据来对体重进行校正
set.seed(123)
data$测定日龄<-as.numeric(data$测定日龄)
#如果有全0体重是算不出来的，需要过滤，但是暂时先不考虑
#isZeroWeight<-data %>% group_by(个体号) %>% summarise(mean(当天体重))

cluster_weight<-function(a){
  cat(as.character(a$个体号[1]),"\n")
  a$idx<-1:nrow(a)
  res<-list()
  fit <- NULL
  attempt<-0
  b<-a
  a<-droplevels(a[a$当天体重>0,])
  a2<-droplevels(b[b$当天体重<=0,])
  while( is.null(fit) && attempt <= 30 ) {
    attempt <- attempt + 1
  #while(is.null(fit)) {
    
    try(
      fit <- flexmix(当天体重 ~ 测定日龄, data = a, k = 3)
    )
  } 
  # if(is.null(fit)){
  #   a<-droplevels(a[a$当天体重>0,])
  #   while(is.null(fit)) {
  #     try(
  #       fit <- flexmix(当天体重 ~ 测定日龄, data = a, k = 3)
  #     )
  #   }
  # }
  
  m<-fit
  clst<-clusters(m)
  clst_val<-unique(clst)
  clst_val<-sort(clst_val)
  clst_num<-length(clst_val)
  pred119<-predict(m, data.frame(测定日龄=119))
  pred144<-predict(m, data.frame(测定日龄=144))
  pred168<-predict(m, data.frame(测定日龄=168))
  
  res$clst<-data.frame(个体号=a$个体号,测定日龄=a$测定日龄,
                          当天体重=a$当天体重,类别=clst,idx=a$idx)
  if(nrow(a2)>0){
    zero<-data.frame(个体号=a2$个体号,测定日龄=a2$测定日龄,
                        当天体重=a2$当天体重,类别=0,idx=a2$idx)
    res$clst<-rbind(res$clst,zero)
    
    res$clst<-res$clst[order(res$clst[,'idx']),]
  }
  summary<-data.frame(min_weight1=NA,
                      max_weight1=NA,
                      min_weight2=NA,
                      max_weight2=NA,
                      min_weight3=NA,
                      max_weight3=NA,
                      sd1=NA,
                      sd2=NA,
                      sd3=NA,
                      clst1_num=NA,
                      clst2_num=NA,
                      clst3_num=NA,
                      clst1_pred119=NA,
                      clst2_pred119=NA,
                      clst3_pred119=NA,
                      clst1_pred144=NA,
                      clst2_pred144=NA,
                      clst3_pred144=NA,
                      clst1_pred168=NA,
                      clst2_pred168=NA,
                      clst3_pred168=NA,
                      k=clst_num
  )
  
  if(clst_num==1){
    summary[1,2*clst_val[1]-1]<-min(a$当天体重[which(clst==clst_val[1])])
    summary[1,2*clst_val[1]]<-max(a$当天体重[which(clst==clst_val[1])])
    summary[1,clst_val[1]+6]<-sd(a$当天体重[which(clst==clst_val[1])])
    summary[1,clst_val[1]+9]<-sum(clst==clst_val[1])
    summary[1,clst_val[1]+12]<-pred119[[paste('Comp',clst_val[1],sep='.')]]
    summary[1,clst_val[1]+15]<-pred144[[paste('Comp',clst_val[1],sep='.')]]
    summary[1,clst_val[1]+18]<-pred168[[paste('Comp',clst_val[1],sep='.')]]
  }
  if(clst_num==2){
    summary[1,2*clst_val[1]-1]<-min(a$当天体重[which(clst==clst_val[1])])
    summary[1,2*clst_val[2]-1]<-min(a$当天体重[which(clst==clst_val[2])])
    summary[1,2*clst_val[1]]<-max(a$当天体重[which(clst==clst_val[1])])
    summary[1,2*clst_val[2]]<-max(a$当天体重[which(clst==clst_val[2])])
    summary[1,clst_val[1]+6]<-sd(a$当天体重[which(clst==clst_val[1])])
    summary[1,clst_val[2]+6]<-sd(a$当天体重[which(clst==clst_val[2])])
    summary[1,clst_val[1]+9]<-sum(clst==clst_val[1])
    summary[1,clst_val[2]+9]<-sum(clst==clst_val[2])
    summary[1,clst_val[1]+12]<-pred119[[paste('Comp',clst_val[1],sep='.')]]
    summary[1,clst_val[2]+12]<-pred119[[paste('Comp',clst_val[2],sep='.')]]
    summary[1,clst_val[1]+15]<-pred144[[paste('Comp',clst_val[1],sep='.')]]
    summary[1,clst_val[2]+15]<-pred144[[paste('Comp',clst_val[2],sep='.')]]
    summary[1,clst_val[1]+18]<-pred168[[paste('Comp',clst_val[1],sep='.')]]
    summary[1,clst_val[2]+18]<-pred168[[paste('Comp',clst_val[2],sep='.')]]
  }
  if(clst_num==3){
    summary[1,2*clst_val[1]-1]<-min(a$当天体重[which(clst==clst_val[1])])
    summary[1,2*clst_val[2]-1]<-min(a$当天体重[which(clst==clst_val[2])])
    summary[1,2*clst_val[3]-1]<-min(a$当天体重[which(clst==clst_val[3])])
    summary[1,2*clst_val[1]]<-max(a$当天体重[which(clst==clst_val[1])])
    summary[1,2*clst_val[2]]<-max(a$当天体重[which(clst==clst_val[2])])
    summary[1,2*clst_val[3]]<-max(a$当天体重[which(clst==clst_val[3])])
    summary[1,clst_val[1]+6]<-sd(a$当天体重[which(clst==clst_val[1])])
    summary[1,clst_val[2]+6]<-sd(a$当天体重[which(clst==clst_val[2])])
    summary[1,clst_val[3]+6]<-sd(a$当天体重[which(clst==clst_val[3])])
    summary[1,clst_val[1]+9]<-sum(clst==clst_val[1])
    summary[1,clst_val[2]+9]<-sum(clst==clst_val[2])
    summary[1,clst_val[3]+9]<-sum(clst==clst_val[3])
    summary[1,clst_val[1]+12]<-pred119[[paste('Comp',clst_val[1],sep='.')]]
    summary[1,clst_val[2]+12]<-pred119[[paste('Comp',clst_val[2],sep='.')]]
    summary[1,clst_val[3]+12]<-pred119[[paste('Comp',clst_val[3],sep='.')]]
    summary[1,clst_val[1]+15]<-pred144[[paste('Comp',clst_val[1],sep='.')]]
    summary[1,clst_val[2]+15]<-pred144[[paste('Comp',clst_val[2],sep='.')]]
    summary[1,clst_val[3]+15]<-pred144[[paste('Comp',clst_val[3],sep='.')]]
    summary[1,clst_val[1]+18]<-pred168[[paste('Comp',clst_val[1],sep='.')]]
    summary[1,clst_val[2]+18]<-pred168[[paste('Comp',clst_val[2],sep='.')]]
    summary[1,clst_val[3]+18]<-pred168[[paste('Comp',clst_val[3],sep='.')]]
  }
  
  res[['']]<-summary
  
  return(res)
}


plot_flexmix<-function(a){
  a$类别<-as.factor(a$类别)
  ggplot(a, aes(x=测定日龄, y=当天体重,color=类别)) + geom_point()
  
}


# # 测试
# b<-droplevels(data[data$个体号=='DDTGXGG19502923',])
b<-droplevels(data[data$个体号=='DDTGXGG19206728',])
clst_res_temp<-by(b,b$个体号,cluster_weight)
clst_res1_temp<-sapply(clst_res_temp,function(x)x[1])
clst_res2_temp<-sapply(clst_res_temp,function(x)x[2])
clst_res1_temp<-rbindlist(clst_res1_temp)
clst_res2_temp<-rbindlist(clst_res2_temp, idcol='个体号')
plot_flexmix(clst_res1_temp)

# 实操
clst_res<-by(data,data$个体号,cluster_weight)
clst_res1<-sapply(clst_res,function(x)x[1])
clst_res2<-sapply(clst_res,function(x)x[2])
clst_res1<-rbindlist(clst_res1)
clst_res2<-rbindlist(clst_res2, idcol='个体号')
plot_flexmix(droplevels(clst_res1[clst_res1$个体号=="DDTGXGG19502923"])) # 这个个体30个循环都没收敛
plot_flexmix(droplevels(clst_res1[clst_res1$个体号=="DDTGXGG19206728"]))

# 对于两类及以上的个体来说，我们发现其某一类的最小体重值最大也不超过52kg
# 所以，建议还是把最小体重值最小的那一类过滤掉
clst_res3<-clst_res2[clst_res2$k>1,]
min_weight<-apply(clst_res3,1,function(x){
  min_weight<-as.numeric(x[c(2,4,6)])
  min_val<-min(min_weight,na.rm=T)
  return(min_val)
})
summary(min_weight)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.13    9.06    9.63   10.28   10.07   51.23 
#plot_flexmix(droplevels(clst_res1[clst_res1$个体号=="DDTGXGG19205764"]))

# 过滤函数

clst_flt<-function(a){
  clst<-a$类别
  clst_val<-unique(clst)
  clst_val<-sort(clst_val)
  clst_val<-clst_val[clst_val!=0]
  clst_num<-length(clst_val)
  min_weight<-rep(NA,3)
  if(clst_num==2){
    min_weight[clst_val[1]]<-min(a$当天体重[which(clst==clst_val[1])])
    min_weight[clst_val[2]]<-min(a$当天体重[which(clst==clst_val[2])])
  }
  if(clst_num==3){
    min_weight[clst_val[1]]<-min(a$当天体重[which(clst==clst_val[1])])
    min_weight[clst_val[2]]<-min(a$当天体重[which(clst==clst_val[2])])
    min_weight[clst_val[3]]<-min(a$当天体重[which(clst==clst_val[3])])
  }
  
  if(clst_num>1){
    min_weight2<-min_weight[!is.na(min_weight)]
    min_val<-min(min_weight2,na.rm=T)
    tryCatch(
      a$过滤<-clst==clst_val[which(min_weight2==min_val)],
      error=function(e) e, warning=function(w) cat("个体号",as.character(a$个体号[1]),"\n"))
    a$过滤[a$类别==0]<-TRUE
  }else{
    a$过滤<-FALSE
    a$过滤[a$类别==0]<-TRUE
  }
  return(a)
}

#测试

unreg<-by(clst_res1,clst_res1$个体号,function(a){
  clst<-a$类别
  clst_val<-unique(clst)
  clst_val<-sort(clst_val)
  if(length(clst_val)==2 & all(clst_val==c(1,3))){
    return(TRUE)
  }else{
    return(FALSE)
  }
})
which(unreg) #极端值
a<-clst_res1[clst_res1$个体号=='YYTGXGG19117001',]
a<-clst_flt(a)
a<-droplevels(a[!a$过滤,])
plot_flexmix(a)

#过滤实操
clst_flt_res<-by(clst_res1,clst_res1$个体号,clst_flt)
clst_flt_res<-rbindlist(clst_flt_res)
data_weight_cln<-droplevels(data[!clst_flt_res$过滤,])
clst_flt_res2<-droplevels(clst_flt_res[!clst_flt_res$过滤,])

#经过上一步，当类别为3时，还剩下两类
#关于这两类保留哪种还是都保留值得商榷
#先看一下还剩下的两类在144日龄的体重预测差值
two_class<-by(clst_flt_res2,clst_flt_res2$个体号,function(x){
  return(length(unique(x$类别))==2)
})
two_class<-unlist(two_class)
two_class_true<-names(two_class)[two_class]

clst_res4<-droplevels(clst_res2[match(two_class_true,clst_res2$个体号),])
diff144<-apply(clst_res4,1,function(x){
  weight<-c(as.numeric(x[2]), as.numeric(x[4]), as.numeric(x[6]))
  min_weight<-min(weight)
  whichNotMin<-which(weight!=min_weight)
  #return(whichNotMin)
  return(abs(as.numeric(x[whichNotMin[1]+16])-as.numeric(x[whichNotMin[2]+16])))
})
summary(diff144)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00031  0.71048  1.39158  1.85311  2.18671 42.82098

plot_flexmix2<-function(a){
  a$类别<-as.factor(a$类别)
  ggplot(a, aes(x=测定日龄, y=当天体重,color=类别)) +
    facet_wrap(~ 个体号, nrow=3) + geom_point()+
    theme(text = element_text(size=15))
}

round_weight<-clst_res4[which(round(diff144)==14),]
plot_flexmix2(droplevels(clst_res1[clst_res1$个体号%in%round_weight$个体号,]))

plot(density(diff144))
#可以看出当两类之间的差值大于10时就很少了，所以，保留高的那个类别

flt_id<-clst_res4[diff144>=10,'个体号'] #需要执行过滤的个体号
flt_id<-as.data.frame(flt_id)
clst_res5<-clst_res4[diff144>=10,] 

# 获得过滤的类别
flt_grp<-apply(clst_res5,1,function(x){
  weight<-c(as.numeric(x[2]), as.numeric(x[4]), as.numeric(x[6]))
  min_weight<-min(weight)
  whichNotMin<-which(weight!=min_weight)
  mid_weight1<-as.numeric(x[whichNotMin[1]+16])
  mid_weight2<-as.numeric(x[whichNotMin[2]+16])
  if(mid_weight1<mid_weight2){
    return(whichNotMin[1])
  }else{
    return(whichNotMin[2])
  }
})
flt_grp<-paste(as.character(clst_res5$个体号),flt_grp,sep=":")
clst_flt_res2_temp<-paste(as.character(clst_flt_res2$个体号),clst_flt_res2$类别,sep=":")

# 过滤
data_weight_cln2<-droplevels(data_weight_cln[!clst_flt_res2_temp%in%flt_grp,])
clst_flt_res3<-droplevels(clst_flt_res2[!clst_flt_res2_temp%in%flt_grp,])
plot_flexmix(droplevels(clst_flt_res3[clst_flt_res3$个体号=="DDTGXGG19206456"])) #可以看出来只剩一类了
#plot_flexmix(droplevels(clst_flt_res3[clst_flt_res3$个体号=="YYTGXGG19129726"])) #可以看出来只剩一类了

#针对利用flexmix过滤后的体重数据，再用
#robust regression校正一遍
#针对每个个体的数据记录编写校正函数
corr_weight<-function(a){
  #首先利用robus regression把权重小于0.5的数据点过滤掉
  res<-list()
  rr.bisquare <- rlm(当天体重 ~ 测定日龄, data=a, psi = psi.bisquare, maxit=100) #RR模型
  # tryCatch(
  #   rr.bisquare <- rlm(当天体重 ~ 测定日龄, data=a, psi = psi.bisquare, maxit=100), #RR模型
  #   error=function(e) e, warning=function(w) cat("个体号",as.character(a$个体号[1]),"\n"))
  # 
  coef<-summary(rr.bisquare)$coefficients
  biweights <- data.frame(测定日龄=a$测定日龄,当天体重=a$当天体重,残差=rr.bisquare$resid, 权重=rr.bisquare$w)
  biweights2<-biweights[biweights$权重>0.5,]
  
  #然后将去除异常值之后的每日体重求中位数
  biweights3<-biweights2 %>% 
    group_by(测定日龄) %>% 
    summarise(体重 = median(当天体重))
  
  #这时我们构造一个新的数据框，存储从起测到终测每一天的体重，先初始化为NA
  days<-seq(from=a$开测日龄[1],to=a$结测日龄[1])
  pred<-predict(rr.bisquare,data.frame(测定日龄=days))
  pred[match(biweights3$测定日龄,days)]<-biweights3$体重
  pred2<-predict(rr.bisquare,data.frame(测定日龄=c(119,168)))
  biweights4<-data.frame(测定日龄=seq(from=a$开测日龄[1],to=a$结测日龄[1]),
                             原始体重=NA,预测体重=pred)
  biweights4$原始体重[match(biweights3$测定日龄,days)]<-biweights3$体重
  biweights5<-data.frame(个体号=a$个体号[1:2],测定日龄=c(119,168),预测体重=pred2)
  biweights6<-data.frame(个体号=a$个体号[1:2],测定日龄=c((100-coef[1,1])/coef[2,1],(115-coef[1,1])/coef[2,1]),设定体重=c(100,115))
  

  res[['']]<-biweights4
  res$pred119_168<-biweights5
  res$days100_115<-biweights6
  return(res)
}

#将该函数应用于原始数据（需要注意的是，rlm有可能不收敛，这一步在开发程序的时候必须要调整，现在
#我暂时通过提高迭代步数解决了这个问题
ord<-order(data$个体号,data$测定日龄,decreasing=c(FALSE,TRUE),method='radix')
data2<-data[ord,]
offset_date<-data2[!duplicated(data2$个体号),c('个体号','测定日龄')]
data_weight_cln2$结测日龄<-offset_date[match(data_weight_cln2$个体号,offset_date$个体号),'测定日龄']

ord<-order(data$个体号,data$测定日龄,decreasing=c(FALSE,FALSE),method='radix')
data2<-data[ord,]
onset_date<-data2[!duplicated(data2$个体号),c('个体号','测定日龄')]
data_weight_cln2$开测日龄<-onset_date[match(data_weight_cln2$个体号,onset_date$个体号),'测定日龄']

bw<-by(data_weight_cln2,data_weight_cln2$个体号,corr_weight)
bw1<-sapply(bw,function(x)x[1])
bw2<-rbindlist(bw1, idcol='个体号')
bw119_168<-sapply(bw,function(x)x[2])
bw119_168<-rbindlist(bw119_168)
days100_115<-sapply(bw,function(x)x[3])
days100_115<-rbindlist(days100_115)

# 计算动态日增重
calc_dg<-function(a){
  days<-seq(from=min(a$测定日龄),to=max(a$测定日龄))
  #bw<-a$体重[match(days,a$测定日龄)]
  bw<-a$预测体重
  dg<-(bw[-1]-head(bw,-1))
  res<-data.frame(测定日龄=days[-1],日增重=dg)
  return(res)
}
dg<-by(bw2,bw2$个体号,calc_dg)
dg<-rbindlist(dg, idcol='个体号')

## 采食量数据填补

# 首先确定采食量
dfi_data<-data
dfi_data[remov_id,'采食量']<-NA
dfi_data2<-dfi_data %>% 
  group_by(个体号,测定日龄) %>% 
  summarise(采食量 = sum(采食量,na.rm=T))

# 然后是体重
a<-paste(bw2$个体号,bw2$测定日龄,sep=":")
b<-paste(dfi_data2$个体号,dfi_data2$测定日龄,sep=":")
weight<-droplevels(bw2[match(b,a),])

# 首先依据进入测定站的时间以周为单位作为一个同期组（批次）
onset_time<-data[!duplicated(data$个体号),]
batch<-paste0(strftime(onset_time$进入时间, format = "%y"),
              strftime(onset_time$进入时间, format = "%V"))
batch<-data.frame(个体号=onset_time$个体号,同期组=batch)
batch<-batch[match(weight$个体号,batch$个体号),]

# 根据flexmix与robust regression计算得到的平均日增重
weight2<-weight
weight2$原始体重[is.na(weight2$原始体重)]<-weight2$预测体重[is.na(weight2$原始体重)] #把原始体重缺失值给填上
ord<-order(weight2$个体号,weight2$测定日龄,decreasing=c(FALSE,TRUE),method='radix')
weight3<-weight2[ord,]
onset_weight<-weight2[!duplicated(weight2$个体号),]
offset_weight<-weight3[!duplicated(weight3$个体号),]
adg<-data.frame(个体号=onset_weight$个体号,
                   日增重=(offset_weight$原始体重-onset_weight$原始体重)/(offset_weight$测定日龄-onset_weight$测定日龄))
adg<-adg[match(weight$个体号,adg$个体号),]
#plot_flexmix(droplevels(clst_flt_res3[clst_flt_res3$个体号=="YYTGXGG19117127"]))

# 性别
sex<-data.frame(个体号=onset_time$个体号,性别=onset_time$性别)
sex<-sex[match(weight$个体号,sex$个体号),]

# 品种品系
breed<-data.frame(个体号=onset_time$个体号,品种品系=onset_time$品种品系)
breed<-breed[match(weight$个体号,breed$个体号),]

# 测定站
station<-data.frame(个体号=onset_time$个体号,测定站=onset_time$测定站)
station<-station[match(weight$个体号,station$个体号),]

#错误类型数据
for(i in 1:16){
  #if(length(get(paste0('et',i)))==0){cat(i,' ')}
  cat(length(get(paste0('et',i))),' ')
}
# et1 et4 et7 这些都是没有出现的错误，跳过
size<-nrow(data)
etd<-data.frame(个体号=data$个体号,测定日龄=data$测定日龄,
                   采食时间=seconds,采食量=data$采食量*1000,
                   et2=rep(FALSE,size),et3=rep(FALSE, size),
                   et5=rep(FALSE,size),et6=rep(FALSE,size),
                   et8=rep(FALSE,size),et9=rep(FALSE,size),
                   et10=rep(FALSE,size),et11=rep(FALSE,size),
                   et12=rep(FALSE,size),et13=rep(FALSE,size),
                   et14=rep(FALSE,size),et15=rep(FALSE,size),
                   et16=rep(FALSE,size))

etd$et2[et2]<-TRUE
etd$et3[et3]<-TRUE
etd$et5[et5]<-TRUE
etd$et6[et6]<-TRUE
etd$et8[et8]<-TRUE
etd$et9[et9]<-TRUE
etd$et10[et10]<-TRUE
etd$et11[et11]<-TRUE
etd$et12[et12]<-TRUE
etd$et13[et13]<-TRUE
etd$et14[et14]<-TRUE
etd$et15[et15]<-TRUE
etd$et16[et16]<-TRUE

# 错误次数占每天采食次数的百分比
etp<-etd %>% 
  group_by(个体号,测定日龄) %>% 
  summarise(
    etp2 = sum(et2)*100/length(et2),
    etp3 = sum(et3)*100/length(et3),
    etp5 = sum(et5)*100/length(et5),
    etp6 = sum(et6)*100/length(et6),
    etp8 = sum(et8)*100/length(et8),
    etp9 = sum(et9)*100/length(et9),
    etp10 = sum(et10)*100/length(et10),
    etp11 = sum(et11)*100/length(et11),
    etp12 = sum(et12)*100/length(et12),
    etp13 = sum(et13)*100/length(et13),
    etp14 = sum(et14)*100/length(et14),
    etp15 = sum(et15)*100/length(et15),
    etp16 = sum(et16)*100/length(et16)
  )

etp<-etd %>% 
  group_by(个体号,测定日龄) %>% 
  summarise(
    etp2 = sum(et2)*100/length(et2),
    etp3 = sum(et3)*100/length(et3),
    etp5 = sum(et5)*100/length(et5),
    etp6 = sum(et6)*100/length(et6),
    etp8 = sum(et8)*100/length(et8),
    etp9 = sum(et9)*100/length(et9),
    etp10 = sum(et10)*100/length(et10),
    etp11 = sum(et11)*100/length(et11),
    etp12 = sum(et12)*100/length(et12),
    etp13 = sum(et13)*100/length(et13),
    etp14 = sum(et14)*100/length(et14),
    etp15 = sum(et15)*100/length(et15),
    etp16 = sum(et16)*100/length(et16)
  )


# 含有错误的采食时间
otd<-etd %>% 
  group_by(个体号,测定日龄) %>% 
  summarise(
    otd2 = sum(et2*采食时间),
    otd6 = sum(et6*采食时间),
    otd8 = sum(et8*采食时间),
    otd9 = sum(et9*采食时间),
    otd10 = sum(et10*采食时间),
    otd11 = sum(et11*采食时间),
    otd12 = sum(et12*采食时间),
    otd13 = sum(et13*采食时间),
    otd14 = sum(et14*采食时间),
  )

# 含有错误的采食量
fid<-etd %>% 
  group_by(个体号,测定日龄) %>% 
  summarise(
    fid5 = sum(et5*采食量),
    fid15 = sum(et15*采食量),
    fid16 = sum(et16*采食量),
  )

# 不含有错误的采食量
err_all<-unique(c(et1,et2,et3,et4,
                  et5,et6,et7,et8,
                  et9,et10,et11,et12,
                  et13,et14,et15,et16))
etd2<-droplevels(etd[-err_all,])

dfi<-etd %>% 
  group_by(个体号,测定日龄) %>% 
  summarise(
    dfi = sum(采食量),
  )

etp<-as.data.frame(etp)
otd<-as.data.frame(otd)
fid<-as.data.frame(fid)
dfi<-as.data.frame(dfi)


dfi2<-rep(0,nrow(weight))
a<-paste(dfi$个体号,dfi$测定日龄,sep=":")
b<-paste(weight$个体号,weight$测定日龄,sep=":")
dfi2[match(a,b)]<-dfi$dfi

input<-data.frame(个体号=weight$个体号,测定日龄=weight$测定日龄,
                     平均日增重=adg$日增重,当天体重=weight$预测体重,
                     品种品系=breed$品种品系,测定站=station$测定站,
                     性别=sex$性别,同期组=batch$同期组,etp2=etp$etp2,
                     etp3=etp$etp3,etp5=etp$etp5,etp6=etp$etp6,
                     etp8=etp$etp8,etp9=etp$etp9,etp10=etp$etp10,
                     etp11=etp$etp11,etp12=etp$etp12,etp13=etp$etp13,
                     etp14=etp$etp14,etp15=etp$etp15,etp16=etp$etp16,
                     otd2=otd$otd2,otd6=otd$otd6,otd8=otd$otd8,
                     otd9=otd$otd9,otd10=otd$otd10,otd11=otd$otd11,
                     otd12=otd$otd12,otd13=otd$otd13,otd14=otd$otd14,
                     fid5=fid$fid5,fid15=fid$fid15,fid16=fid$fid16,
                     采食量=dfi2
)
input$个体号<-factor(input$个体号,levels=unique(as.character(input$个体号)))

write.csv(input,file="inputForBGLR.csv",row.names = F,quote=F,fileEncoding = "utf-8")

############################################ call BGLR
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


fm<-BGLR(y=input$采食量,ETA=ETA,nIter=nIter, burnIn=burnIn,saveAt=saveAt,verbose = T)

##########################################

## 读取计算完毕的模型
load("fm.RData")

#诊断，都收敛了
varE<-scan('CORRvarE.dat')
plot(varE,type='o',col=2,cex=.5,ylab=expression(var[e]));
abline(h=fm$varE,col=4,lwd=2);
abline(v=fm$burnIn/fm$thin,col=4)

varb<-scan('CORRETA_2_varB.dat')
plot(varb,type='o',col=2,cex=.5,ylab=expression(Varb));
abline(h=fm$ETA[[2]]$varB,col=4,lwd=2);
abline(v=fm$burnIn/fm$thin,col=4)

b<-fread('CORRETA_1_b.dat')
plot(b$otd10,type='o',col=2,cex=.5,ylab=expression(otd10));

# 我们所关注的错误相关的各个固定效应
coef<-tail(fm$ETA[[1]]$b,25)
input_et<-input[,9:(ncol(input)-1)]
input_et<-as.matrix(input_et)
corr_et<-input_et%*%coef
bin<-get_bin(abs(corr_et),100)
freq_otv2=hist(abs(corr_et), breaks=bin, include.lowest=TRUE, plot=TRUE)
adj_dfi<-input$采食量-corr_et
adj_dfi[adj_dfi<0]<-0 #如果校正后还有小于0的数据的话，校正为0
input$校正采食量<-adj_dfi/1000

# 利用robust regression获得119-168日龄的累积采食量
corr_adfi<-function(a){
  tryCatch(
    rr.bisquare <- lm(校正采食量 ~ 测定日龄, data=a),
    error=function(e) e, warning=function(w) cat("个体号",as.character(a$个体号[1]),"\n"))
  
  coef<-summary(rr.bisquare)$coefficients
  days_dfi<-seq(from=119,to=168)
  pred<-predict(rr.bisquare,data.frame(测定日龄=days_dfi))
  pred_dfi<-data.frame(测定日龄=days_dfi,
                       校正采食量=a[match(days_dfi,a$测定日龄),'校正采食量'])
  pred_dfi$校正采食量[is.na(pred_dfi$校正采食量)]<-pred[is.na(pred_dfi$校正采食量)]
  res<-data.frame(个体号=as.character(a$个体号[1]),总采食量=sum(pred_dfi$校正采食量))
  
  return(res)
}

fi119_168<-by(input,input$个体号,corr_adfi)
fi119_168<-rbindlist(fi119_168)
# 基于校正采食量计算料重比
adfi<-fi119_168
adfi$平均日采食量<-adfi$总采食量/(168-119+1)
adg119_168<-bw119_168 %>% group_by(个体号)%>%summarise(增重=预测体重[2]-预测体重[1])
adfi$总增重<-adg119_168$增重[match(adfi$个体号,adg119_168$个体号)]
adfi$平均日增重<-adfi$总增重/(168-119+1)
fcr<-as.data.frame(adfi)
fcr$料重比<-fcr$平均日采食量/fcr$平均日增重

# 达目标体重日龄
days100<-data.frame(个体号=days100_115[days100_115$设定体重==100,]$个体号,
                       测定日龄=round(days100_115[days100_115$设定体重==100,]$测定日龄))
days115<-data.frame(个体号=days100_115[days100_115$设定体重==115,]$个体号,
                       测定日龄=round(days100_115[days100_115$设定体重==115,]$测定日龄))


# 其他的一些与采食行为相关的性状
data_cln<-droplevels(data[-err_all,]) #计算的时候可以把错误类型数据去除，因为反正算的也是平均值
data_cln$采食秒数<-seconds[-err_all] 

afiv<-data_cln %>% 
  group_by(个体号) %>% 
  summarise(
    平均次采食量 = mean(采食量),
  ) # 平均次采食量
afiv<-as.data.frame(afiv)

afrv<-data_cln[data_cln$采食秒数!=0,] %>% 
  group_by(个体号) %>% 
  summarise(
    平均采食速度=mean(采食量/(采食秒数/60)),
  ) # 把采食时间等于0的剔除，防止产生极端值
afrv<-as.data.frame(afrv)

aoid<-data_cln %>% 
  group_by(个体号) %>% 
  summarise(
    平均次采食时间 = mean(采食秒数),
  ) # 平均次采食时间
aoid<-as.data.frame(aoid)

nvd<-data %>% 
  group_by(个体号) %>% 
  summarise(
    平均采食次数 = length(测定日龄)/length(unique(测定日龄)),
  ) # 日平均采食次数
nvd<-as.data.frame(nvd)

# 平均性状的普通BLUP模型
## 将此数据整理成包含各个效应的可用于后续遗传评估的数据
matc<-match(fcr$个体号,input$个体号)
fcr$同期组<-input$同期组[matc]
fcr$性别<-input$性别[matc]
fcr$测定站<-input$测定站[matc]
fcr$品种品系<-input$品种品系[matc]

matc<-match(fcr$个体号,afiv$个体号)
fcr$次平均采食量<-afiv$平均次采食量[matc]
matc<-match(fcr$个体号,afrv$个体号)
fcr$次平均采食速度<-afrv$平均采食速度[matc]
matc<-match(fcr$个体号,aoid$个体号)
fcr$次平均采食时间<-aoid$平均次采食时间[matc]
matc<-match(fcr$个体号,nvd$个体号)
fcr$日平均采食次数<-nvd$平均采食次数[matc]

bw119<-bw119_168[bw119_168$测定日龄==119,]
matc<-match(fcr$个体号,bw119$个体号)
fcr$起测体重<-bw119$预测体重[matc]

matc<-match(fcr$个体号,days100$个体号)
fcr$达100公斤体重日龄<-days100$测定日龄[matc]
matc<-match(fcr$个体号,days115$个体号)
fcr$达115公斤体重日龄<-days115$测定日龄[matc]

fcr$性别[fcr$性别=='混合']<-'公'
fcr$性别<-factor(fcr$性别,levels=c("公","母"))
duroc<-droplevels(fcr[fcr$品种品系=='美系杜洛克',])
duroc<-data.frame(个体号=duroc$个体号,性别=as.numeric(duroc$性别),
                     测定站=duroc$测定站,同期组=duroc$同期组,
                     达100公斤体重日龄=duroc$达100公斤体重日龄,
                     达115公斤体重日龄=duroc$达115公斤体重日龄,
                     平均日采食量=duroc$平均日采食量,
                     平均日增重=duroc$平均日增重,料重比=duroc$料重比,
                     次平均采食量=duroc$次平均采食量,次平均采食速度=duroc$次平均采食速度,
                     次平均采食时间=duroc$次平均采食时间,日平均采食次数=duroc$日平均采食次数,
                     起测体重=duroc$起测体重)

landrace<-droplevels(fcr[fcr$品种品系=='美系长白',])
landrace<-data.frame(个体号=landrace$个体号,性别=as.numeric(landrace$性别),
                        测定站=landrace$测定站,同期组=landrace$同期组,
                        达100公斤体重日龄=landrace$达100公斤体重日龄,
                        达115公斤体重日龄=landrace$达115公斤体重日龄,
                        平均日采食量=landrace$平均日采食量,
                        平均日增重=landrace$平均日增重,料重比=landrace$料重比,
                        次平均采食量=landrace$次平均采食量,次平均采食速度=landrace$次平均采食速度,
                        次平均采食时间=landrace$次平均采食时间,日平均采食次数=landrace$日平均采食次数,
                        起测体重=landrace$起测体重)

yorkshire<-droplevels(fcr[fcr$品种品系=='美系大白',])
yorkshire<-data.frame(个体号=yorkshire$个体号,性别=as.numeric(yorkshire$性别),
                         测定站=yorkshire$测定站,同期组=yorkshire$同期组,
                         达100公斤体重日龄=yorkshire$达100公斤体重日龄,
                         达115公斤体重日龄=yorkshire$达115公斤体重日龄,
                         平均日采食量=yorkshire$平均日采食量,
                         平均日增重=yorkshire$平均日增重,料重比=yorkshire$料重比,
                         次平均采食量=yorkshire$次平均采食量,次平均采食速度=yorkshire$次平均采食速度,
                         次平均采食时间=yorkshire$次平均采食时间,日平均采食次数=yorkshire$日平均采食次数,
                         起测体重=yorkshire$起测体重)

pietrain<-droplevels(fcr[fcr$品种品系=='美系皮特兰',])
pietrain<-data.frame(个体号=pietrain$个体号,性别=as.numeric(pietrain$性别),
                        测定站=pietrain$测定站,同期组=pietrain$同期组,
                        达100公斤体重日龄=pietrain$达100公斤体重日龄,
                        达115公斤体重日龄=pietrain$达115公斤体重日龄,
                        平均日采食量=pietrain$平均日采食量,
                        平均日增重=pietrain$平均日增重,料重比=pietrain$料重比,
                        次平均采食量=pietrain$次平均采食量,次平均采食速度=pietrain$次平均采食速度,
                        次平均采食时间=pietrain$次平均采食时间,日平均采食次数=pietrain$日平均采食次数,
                        起测体重=pietrain$起测体重)

eb5<-droplevels(fcr[fcr$品种品系=='EB5杜洛克',])
eb5<-data.frame(个体号=eb5$个体号,性别=as.numeric(eb5$性别),
                   测定站=eb5$测定站,同期组=eb5$同期组,
                   达100公斤体重日龄=eb5$达100公斤体重日龄,
                   达115公斤体重日龄=eb5$达115公斤体重日龄,
                   平均日采食量=eb5$平均日采食量,
                   平均日增重=eb5$平均日增重,料重比=eb5$料重比,
                   次平均采食量=eb5$次平均采食量,次平均采食速度=eb5$次平均采食速度,
                   次平均采食时间=eb5$次平均采食时间,日平均采食次数=eb5$日平均采食次数,
                   起测体重=eb5$起测体重)

out_dir<-"LMM"
dir.create(out_dir, showWarnings = FALSE)

## 整理系谱
library(pedigree)

add.Inds_new<-function(ped){
  head <- colnames(ped)
  ndams <- match(ped[, 2], ped[, 1])
  ndams <- as.character(unique(ped[is.na(ndams), 2]))
  ndams <- ndams[!is.na(ndams)]
  nsires <- match(ped[, 3], ped[, 1])
  nsires <- as.character(unique(ped[is.na(nsires), 3]))
  nsires <- nsires[!is.na(nsires)]
  nped <- data.frame(matrix(NA, nrow = length(ndams) + length(nsires), 
                            ncol = ncol(ped)))
  colnames(nped) <- colnames(ped)
  nped[, 1] <- c(ndams, nsires)
  ped <- rbind(nped, ped)
  colnames(ped) <- head
  return(ped)
}

get_dmu_dat<-function(name){
  x<-get(name)
  res<-list()
  a<-match(x$个体号,ped$个体号)
  ped_dat<-ped[a,]
  ped_dat<-as.data.frame(ped_dat)
  cols = c(1,28:ncol(ped))
  ped_dat[,cols] <- apply(ped_dat[,cols], 2, function(x) as.character(x))
  ind<-c(ped_dat$个体号,ped_dat$父亲,ped_dat$母亲,ped_dat$父父,ped_dat$父母,ped_dat$母父,ped_dat$母母)
  sire<-c(ped_dat$父亲,ped_dat$父父,ped_dat$母父,ped_dat$父父父,ped_dat$父母父,ped_dat$母父父,ped_dat$母母父)
  dam<-c(ped_dat$母亲,ped_dat$父母,ped_dat$母母,ped_dat$父父母,ped_dat$父母母,ped_dat$母父母,ped_dat$母母母)
  ped_dat2<-cbind(ind,sire,dam)
  ped_dat2<-add.Inds_new(ped_dat2)
  
  ord <- orderPed(ped_dat2)
  ped_dat2<-ped_dat2[order(ord),]
  ped_dat2[,1]<-as.character(ped_dat2[,1])
  ped_dat2[,2]<-as.character(ped_dat2[,2])
  ped_dat2[,3]<-as.character(ped_dat2[,3])
  ord2<-1:nrow(ped_dat2)
  ped_dat2<-cbind(ped_dat2,ord2)
  
  ped_dat2[is.na(ped_dat2)]<-0
  
  name<-toupper(name)
  dir.create(paste0(out_dir,"/",name), showWarnings = FALSE)
  write.table(ped_dat2,file=paste0(out_dir,"/",name,"/","idc.txt"),row.names = F,quote = F,sep="\t")
  phe_ind<-ped_dat2[match(x$个体号,ped_dat2[,1]),4]
  ped_ind<-ped_dat2$ord2
  sire_ind<-ped_dat2[match(ped_dat2[,2],ped_dat2[,1]),4]
  dam_ind<-ped_dat2[match(ped_dat2[,3],ped_dat2[,1]),4]
  ped_dat3<-data.frame(ped_ind,sire_ind,dam_ind,ord2)
  ped_dat3[is.na(ped_dat3)]<-0
  litter<-ped_dat3$dam_ind[match(phe_ind,ped_dat3$ped_ind)]
  
  res$DMU_PED<-ped_dat3
  write.table(ped_dat3,file=paste0(out_dir,"/",name,"/","DMU_PED"),
              row.names = F,col.names=F,quote = F,sep="\t")
  
  parity<-ped$出生胎次[match(x$个体号,ped$个体号)]
  dat_dmu<-data.frame(ID=phe_ind,Mu=rep(1,length(phe_ind)),Pen=x$测定站,
                      Sex=x$性别,Batch=x$同期组,Parity=parity,Dam=litter,
                      DAY100=x$达100公斤体重日龄,
                      DAY115=x$达115公斤体重日龄,
                      ADFI=x$平均日采食量,
                      ADG=x$平均日增重,FCR=x$料重比,
                      FIV=x$次平均采食量,FRV=x$次平均采食速度,
                      OID=x$次平均采食时间,NVD=x$日平均采食次数,
                      BW1=x$起测体重)

  write.table(dat_dmu,file=paste0(out_dir,"/",name,"/","DMU_LMM_DAT"),row.names = F,col.names=F,quote = F,sep="\t")
  
  res$DMU_DAT<-dat_dmu
  return(res)
}


duroc_dmu<-get_dmu_dat(name='duroc')
landrace_dmu<-get_dmu_dat(name='landrace')
yorkshire_dmu<-get_dmu_dat(name='yorkshire')
pietrain_dmu<-get_dmu_dat(name='pietrain')
eb5_dmu<-get_dmu_dat(name='eb5')

cat(paste0(paste(colnames(eb5_dmu$DMU_DAT),collapse="\t"),"\n"),file=paste0(out_dir,"/","DMU_LMM_DAT_HEADER"))

##随机回归模型的数据准备

# 以天作为一个时间单位来进行分析

##日增重/前一天的采食量=料重比
a<-paste(dg$个体号,dg$测定日龄,sep=":")
b<-paste(input$个体号,input$测定日龄,sep=":")
dg_info<-droplevels(input[match(a,b),c('个体号','测定日龄','品种品系','测定站','性别','同期组')])
dg<-merge(dg,dg_info)
dg$测定前日龄<-dg$测定日龄-1
dg$测定前日采食量<-input[match(paste(dg$个体号,dg$测定前日龄,sep=":"),
                        paste(input$个体号,input$测定日龄,sep=":")),'校正采食量']
dg$料重比<-dg$测定前日采食量/(dg$日增重)
input$料重比<-dg[match(paste(input$个体号,input$测定日龄,sep=":"),
                    paste(dg$个体号,dg$测定日龄,sep=":")),'料重比']
input$日增重<-dg[match(paste(input$个体号,input$测定日龄,sep=":"),
                    paste(dg$个体号,dg$测定日龄,sep=":")),'日增重']
input$料重比<-unlist(input$料重比)
input$料重比[(!is.na(input$料重比) & (input$料重比>4.5 | input$料重比<0))]<-NA #料重比应该有个限度
plot(density(input$料重比[!is.na(input$料重比)]))

# Legendre多项式的Phi矩阵
out_dir<-"RRM_DAY"
dir.create(out_dir, showWarnings = FALSE)

get_dmu_dat_leg<-function(breed,name){
  x<-droplevels(input[input$品种品系==breed,])
  leg4coef <- legendre.polynomials(n=4, normalized=TRUE)
  leg4<- as.data.frame(polynomial.values(polynomials=leg4coef,x=scaleX(x$测定日龄, u=-1, v=1)))
  colnames(leg4)<-paste0('phi',1:5)
  x<-cbind(x,leg4)
  dir.create(paste0(out_dir,"/",name), showWarnings = FALSE)
  
  idc<-read.table(paste0("LMM/",name,"/","idc.txt"),h=T,sep="\t")
  phe_ind<-idc$ord2[match(x$个体号,idc$ind)]
  
  litter<-idc$dam[match(x$个体号,idc$ind)]
  litter<-idc$ord2[match(litter,idc$ind)]
  
  parity<-ped$出生胎次[match(x$个体号,ped$个体号)]
  
  dat_dmu<-data.frame(ID=phe_ind,Pen=x$测定站,
                      Sex=as.numeric(x$性别),
                      Batch=x$同期组,Days=x$测定日龄,
                      Dam=litter,Parity=parity,
                      DFI=x$校正采食量,BW=x$当天体重,
                      FCR=x$料重比,lg0=x$phi1,
                      lg1=x$phi2,lg2=x$phi3,
                      lg3=x$phi4,lg4=x$phi5
  )

  dat_dmu[is.na(dat_dmu)]<--999

  
  write.table(dat_dmu,file=paste0(out_dir,"/",name,"/","DMU_RRM_DAT"),
              row.names = F,col.names=F,quote = F,sep="\t")

  return(dat_dmu)
}

duroc_rrm<-get_dmu_dat_leg('美系杜洛克','DUROC')
landrace_rrm<-get_dmu_dat_leg('美系长白','LANDRACE')
yorkshire_rrm<-get_dmu_dat_leg('美系大白','YORKSHIRE')
pietrain_rrm<-get_dmu_dat_leg('美系皮特兰','PIETRAIN')
eb5_rrm<-get_dmu_dat_leg('EB5杜洛克','EB5')
cat(paste0(paste(colnames(eb5_rrm),collapse="\t"),"\n"),file=paste0(out_dir,"/","DMU_RRM_DAT_HEADER"))

# 以周作为一个时间单位来进行分析
out_dir<-"RRM_WK"

get_dmu_dat_week<-function(breed,name){
  
  a<-droplevels(input[input$品种品系==breed,])
  day_start<-min(a$测定日龄)
  day_end<-max(a$测定日龄)
  a$周次<-cut(a$测定日龄,breaks=seq(from=day_start-1,to=day_end+7,by=7),include.lowest = F)
  
  week_dat<-a %>% 
    dplyr::group_by(个体号,周次,.drop=FALSE) %>% 
    dplyr::summarise(
      周平均日采食量=sum(校正采食量)/(n()*1000),
      周平均日增重=(dplyr::last(当天体重)-dplyr::first(当天体重))/n(),
      周记录数=n()
    )
  week_dat<-droplevels(week_dat[week_dat$周记录数>=4,]) ###周记录数至少要大于3条
  week_dat$料重比<-week_dat$周平均日采食量/week_dat$周平均日增重
  week_dat<-droplevels(week_dat[week_dat$料重比>=0 & week_dat$料重比<=4.5,])
  week_dat<-merge(week_dat,fcr[,c('个体号','测定站','性别','同期组')],by='个体号')
  week_dat$周日龄<-as.numeric(gsub("^\\((\\d+),.*","\\1",week_dat$周次))
  week_dat<-week_dat %>% group_by(周日龄) %>% filter(n() >= 50)
  leg4coef<-legendre.polynomials(n=4, normalized=TRUE)
  leg4<-as.data.frame(polynomial.values(polynomials=leg4coef,x=scaleX(week_dat$周日龄, u=-1, v=1)))
  colnames(leg4)<-paste0('phi',1:5)
  
  idc<-read.table(paste0('LMM/',name,"/","idc.txt"),h=T,sep="\t")
  phe_ind<-idc$ord2[match(week_dat$个体号,idc$ind)]
  
  litter<-idc$dam[match(week_dat$个体号,idc$ind)]
  litter<-idc$ord2[match(litter,idc$ind)]
  
  parity<-ped$出生胎次[match(week_dat$个体号,ped$个体号)]
  dat_dmu<-data.frame(ID=phe_ind,Pen=week_dat$测定站,
                      Sex=as.numeric(week_dat$性别),
                      Batch=week_dat$同期组,Days=week_dat$周日龄,
                      Dam=litter,Parity=parity,
                      WDFI=week_dat$周平均日采食量,
                      WADG=week_dat$周平均日增重,
                      FCR=week_dat$料重比,lg0=leg4$phi1,
                      lg1=leg4$phi2,lg2=leg4$phi3,
                      lg3=leg4$phi4,lg4=leg4$phi5
  )
  dir.create(paste0(out_dir,'/',name),showWarnings = F)
  write.table(dat_dmu,file=paste0(out_dir,'/',name,"/","DMU_RRM_DAT"),
              row.names = F,col.names=F,quote = F,sep="\t")
  
  return(dat_dmu)
}

duroc_wk<-get_dmu_dat_week('美系杜洛克','DUROC')
landrace_wk<-get_dmu_dat_week('美系长白','LANDRACE')
yorkshire_wk<-get_dmu_dat_week('美系大白','YORKSHIRE')
pietrain_wk<-get_dmu_dat_week('美系皮特兰','PIETRAIN')
eb5_wk<-get_dmu_dat_week('EB5杜洛克','EB5')
cat(paste0(paste(colnames(eb5_wk),collapse="\t"),"\n"),file=paste0(out_dir,"/","DMU_RRM_DAT_HEADER"))

save.image(file='rfi.RData') 
