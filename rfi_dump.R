
otd<-apply(TureOrFalse,1,function(x){
  val= as.numeric(x[3])
  multip = as.numeric(x[5:length(x)])
  return(val*multip)
})
otd<-t(otd)
otd<-as.data.frame(otd)
colnames(otd)<-paste0("etp",c(2,3,5,6,8:16))

fid<-apply(TureOrFalse,1,function(x){
  val= as.numeric(x[4])
  multip = as.numeric(x[5:length(x)])
  return(val*multip)
})
fid<-t(fid)
fid<-as.data.frame(fid)
colnames(fid)<-paste0("etp",c(2,3,5,6,8:16))

# ######################################################################################
# # 由于内存不足，将数据写出到服务器里计算 Rscript call_mice.R > call_mice.log 2>&1 &
# #填补输入数据
# impute_input<-data.frame(个体号=weight$个体号,测定日龄=weight$测定日龄,
#                             平均日增重=adg$日增重,当天体重=weight$预测体重,
#                             品种品系=breed$品种品系,测定站=station$测定站,
#                             性别=sex$性别,同期组=batch$同期组,采食量=dfi_data$采食量
# )
# write.csv(impute_input,file="impute_input.csv",row.names = F,quote=F,fileEncoding = "utf-8")
# 
# library(data.table)#
# library(mice)
# impute_input<-fread("impute_input.csv")
# impute<-mice(impute_input, m=5, maxit = 10, method = 'pmm', seed = 500)
# impute <- parlmice(impute_input,m=5,
#                    maxit=20,method='pmm',
#                    cluster.seed=1234,
#                    n.core=100,n.imp.core=1)
# save(impute,file="impute_res.RData")
# ######################################################################################

# ## 开测体重
# ord<-order(weight$个体号,weight$测定日龄,decreasing=c(FALSE,TRUE),method='radix')
# weight4<-weight[ord,]
# onset_weight<-weight[!duplicated(weight$个体号),]
# offset_weight<-weight4[!duplicated(weight4$个体号),]
# 
# onset_weight<-cbind(onset_weight,fcr$品种品系)
# plot(density(onset_weight$原始体重))
# 
# ## 结测体重
# onset_weight<-cbind(onset_weight,fcr$品种品系)
# plot(density(onset_weight$原始体重))


## 接下来利用mice包对前期过滤掉的错误数据的采食量进行填充
## 方法为multiple imputation (Jiao et al. 2016)

# 首先确定采食量
dfi_data<-data
dfi_data[remov_id,'采食量']<-NA
dfi_data2<-dfi_data %>% 
  group_by(个体号,测定日龄) %>% 
  summarise(采食量 = sum(采食量))

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
weight2<-droplevels(weight[!is.na(weight$原始体重),])
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


#填补输入数据
impute_input<-data.frame(个体号=weight$个体号,测定日龄=weight$测定日龄,
                            平均日增重=adg$日增重,当天体重=weight$预测体重,
                            品种品系=breed$品种品系,测定站=station$测定站,
                            性别=sex$性别,同期组=batch$同期组,采食量=dfi_data2$采食量
)
write.csv(impute_input,file="impute_input.csv",row.names = F,quote=F,fileEncoding = "utf-8")

######################################################################################
# 由于内存不足，将数据写出到服务器里计算 Rscript call_mice.R > call_mice.log 2>&1 &
library(data.table)#
library(mice)
impute_input<-fread("impute_input.csv")
impute<-mice(impute_input, m=5, maxit = 10, method = 'pmm', seed = 500)
impute <- parlmice(impute_input,m=5,
                   maxit=20,method='pmm',
                   cluster.seed=1234,
                   n.core=100,n.imp.core=1)
save(impute,file="impute_res.RData")
######################################################################################



#wd
stat_res<-box(wd)
wd_remov_id<-which(wd%in%stat_res$remov)
wd2<-wd[!wd%in%stat_res$remov]
bin<-get_bin(wd2,0.1)
freq_wd2=hist(wd2, breaks=bin, include.lowest=TRUE, plot=TRUE)

#之后利用条形图统计5s作为一个区间的频率
bin<-seq(min(seconds2),max(seconds2),by=5)
freq=hist(seconds2, breaks=bin, include.lowest=TRUE, plot=TRUE)
boxplot(seconds2)


## 每天有几条记录？
#先确定每天进入测定站的日期是第几天
onset_date<-strsplit(as.character(data$进入时间),split=" ")
onset_date<-sapply(onset_date,function(x)x[1])
onset_date<-as.POSIXct(onset_date,format="%Y-%m-%d")
data$开测日期<-onset_date
data$测定天数<-as.numeric(data$开测日期-as.POSIXct(as.character(data$出生日期),format="%Y/%m/%d"))-data$开测日龄

record<-data %>% 
  group_by(个体号,测定天数) %>% 
  summarise(n = n())

# 在这一步就发现有的猪一天进入测定站的次数非常多
# 最多居然的一个居然高达197次
# 好多都是进去待个1-2s
# 所以统计一下采食时间的分布
seconds_stat<-boxplot(seconds)
#先利用箱线图过滤掉一线极端值
remov_seconds<-seconds_stat$out
seconds2<-seconds[!seconds%in%remov_seconds]

#之后利用条形图统计5s作为一个区间的频率
bin<-seq(min(seconds2),max(seconds2),by=5)
freq=hist(seconds2, breaks=bin, include.lowest=TRUE, plot=TRUE)
boxplot(seconds2)




ggplot(data, aes(x = 品种品系, y = n)) +
  geom_boxplot() 

ggplot(record, aes(x = as.factor(record$测定天数), y = n)) +
  geom_boxplot() 

#####生长曲线绘制
duroc_dat<-read.csv("duroc_dat.csv",h=T)






##探究一下测定阶段,DDTCZZC17217503结测体重测定有问题，去除其最后一次测定结果
remov<-as.numeric(tail(rownames(duroc_dat[as.character(duroc_dat$个体号)=="DDTCZZC17217503",]),1))
duroc_dat<-droplevels(duroc_dat[-remov,])
duroc_dat2<-duroc_dat[order(duroc_dat$个体号,duroc_dat$进入时间,decreasing=c(FALSE,TRUE)),]

onset_weight<-duroc_dat[!duplicated(duroc_dat$个体号),c('个体号','当天体重')]
colnames(onset_weight)[2]<-'开测体重'

offset_weight<-duroc_dat2[!duplicated(duroc_dat2$个体号),c('个体号','当天体重')]
colnames(offset_weight)[2]<-'结测体重'



stage<-merge(onset_weight,offset_weight,by='个体号')

theme <-theme_get()
theme$text$family <- "STSong"
theme_set(theme)

png(file="原始测定区间.png",width=12, height=8,units='in',family="GB1",res=300)

ggplot(stage) + geom_dumbbell(aes(y = 个体号,
                                  x = 开测体重, xend = 结测体重),
                              colour = "grey60", size = 1,
                              colour_x = "#F7BC08", colour_xend = "#395B74")+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  scale_x_continuous(name="体重阶段(kg)", breaks = seq(30,150,by=10))+ylab("")
dev.off()

####set onset range [60,80], offset range [100,120]
stage1<-cbind(stage,data.frame(开测体重2=60,结测体重2=80))
start<-with(stage1, (开测体重2 <= 开测体重 & 结测体重2 >= 开测体重)|(开测体重2 >= 开测体重 & 结测体重2 <= 结测体重))
stage2<-cbind(stage,data.frame(开测体重2=100,结测体重2=120))
end<-with(stage2, (开测体重2 <= 结测体重 & 结测体重2 >= 结测体重)|(开测体重2 >= 开测体重 & 结测体重2 <= 结测体重))

stage<-stage[start&end,]

png(file="过滤异常值之后所有个体的测定区间.png",width=12, height=8,units='in',family="GB1",res=300)

ggplot(stage) + geom_dumbbell(aes(y = 个体号,
                                  x = 开测体重, xend = 结测体重),
                              colour = "grey60", size = 1,
                              colour_x = "#F7BC08", colour_xend = "#395B74")+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  scale_x_continuous(name="体重阶段", breaks = seq(30,150,by=10))+ylab("")
dev.off()

duroc_dat<-droplevels(duroc_dat[duroc_dat$个体号%in%stage$个体号,])
duroc_dat<-droplevels(duroc_dat[duroc_dat$当天体重>=60 & duroc_dat$当天体重<=120,])

###DDTCZZC17217779 的最后一条记录也有问题，剔除
remov<-which(row.names(duroc_dat)==tail(rownames(duroc_dat[as.character(duroc_dat$个体号)=="DDTCZZC17217779",]),1))
duroc_dat<-droplevels(duroc_dat[-remov,])

duroc_dat2<-duroc_dat[order(duroc_dat$个体号,duroc_dat$进入时间,decreasing=c(FALSE,TRUE)),]


onset_weight<-duroc_dat[!duplicated(duroc_dat$个体号),c('个体号','当天体重')]
colnames(onset_weight)[2]<-'开测体重'

offset_weight<-duroc_dat2[!duplicated(duroc_dat2$个体号),c('个体号','当天体重')]
colnames(offset_weight)[2]<-'结测体重'
stage<-merge(onset_weight,offset_weight,by='个体号')


png(file="最终所有个体的测定区间.png",width=12, height=8,units='in',family="GB1",res=300)
ggplot(stage) + geom_dumbbell(aes(y = 个体号,
                                  x = 开测体重, xend = 结测体重),
                              colour = "grey60", size = 1,
                              colour_x = "#F7BC08", colour_xend = "#395B74")+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  scale_x_continuous(name="体重阶段", breaks = seq(30,150,by=10))
dev.off()

######################################################################
feed<-aggregate(duroc_dat$采食量,by=list(duroc_dat$个体号),sum)
colnames(feed)<-c('个体号','采食量')

date1<-strsplit(as.character(duroc_dat$进入时间),split=" ")
date1<-sapply(date1,function(x)x[1])

date2<-strsplit(as.character(duroc_dat$退出时间),split=" ")
date2<-sapply(date2,function(x)x[1])


duroc_dat$进入日期 <- as.Date(date1,format="%Y-%m-%d")
duroc_dat$退出日期 <- as.Date(date2,format="%Y-%m-%d")

onset_date<-duroc_dat[!duplicated(duroc_dat$个体号),c('个体号','进入日期')]
colnames(onset_date)[2]<-'开测日期'
duroc_dat2<-duroc_dat[order(duroc_dat$个体号,duroc_dat$进入时间,decreasing=c(FALSE,TRUE)),]
offset_date<-duroc_dat2[!duplicated(duroc_dat2$个体号),c('个体号','退出日期')]
colnames(offset_date)[2]<-'结测日期'

nvd<-table(duroc_dat$个体号)
nvd<-data.frame(nvd)
colnames(nvd)<-c('个体号','访问次数')

tpd<-aggregate(duroc_dat$采食时间,by=list(duroc_dat$个体号),function(x)sum(as.numeric(substr(as.character(x),4,5))))
colnames(tpd)<-c('个体号','访问时长')

station<-duroc_dat[!duplicated(duroc_dat$个体号),c('个体号','测定站')]

dmu<-Reduce(function(x, y) merge(x, y, by='个体号'), list(station,onset_weight,offset_weight,onset_date,offset_date,feed,nvd,tpd))

dmu$测定天数<-as.numeric(difftime(dmu$结测日期,dmu$开测日期,units='days'))
dmu$增重<-dmu$结测体重-dmu$开测体重
dmu$平均日增重<-dmu$增重/dmu$测定天数
dmu$平均日采食量<-dmu$采食量/dmu$测定天数
dmu$料重比<-dmu$平均日采食量/dmu$平均日增重
dmu$平均访问次数<-dmu$访问次数/dmu$测定天数
dmu$平均访问时长<-(dmu$访问时长/dmu$访问次数)*dmu$平均访问次数

dmu<-cbind(dmu,ped[match(dmu$个体号,ped$个体号),c('出生日期','性别','本地猪出生猪场','品种品系')])

getSeason <- function(DATE) {
  if(is.na(DATE)){
    return(NA)
  }else{
    month <- as.numeric(format(DATE,'%m'))
    if (month >= 12 | month <= 2){return(4)}
    if (month >= 3 & month <= 5){return(1)}
    if (month >= 6 & month <= 8){return(2)}
    if (month >= 9 & month <= 11){return(3)}
  }
}

ys<-rep(NA,nrow(dmu))

for(i in 1:nrow(dmu)){
  date<-dmu[i,"开测日期"]
  year<-substr(format(date,'%Y'),4,4)
  season<-getSeason(date)
  ys[i]<-paste0(year,season)
}
head(ys)

dmu$年季<-ys

dmu$Sex<-as.numeric(dmu$性别)
pen<-as.numeric(dmu$测定站)

#####pedigree
library(pedigree)

ped_ori<-read.csv("../sire/Duroc_ped.CSV",h=T,colClasses = "character")
ind<-c(ped_ori$个体号,ped_ori$父亲,ped_ori$母亲,ped_ori$父父,ped_ori$父母,ped_ori$母父,ped_ori$母母)
sire<-c(ped_ori$父亲,ped_ori$父父,ped_ori$母父,ped_ori$父父父,ped_ori$父母父,ped_ori$母父父,ped_ori$母母父)
dam<-c(ped_ori$母亲,ped_ori$父母,ped_ori$母母,ped_ori$父父母,ped_ori$父母母,ped_ori$母父母,ped_ori$母母母)
ped<-cbind(ind,sire,dam)
ped<-unique(ped[,1:3])
ped[ped==""]<-NA

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

ped<-add.Inds_new(ped)

ord <- orderPed(ped)
ped<-ped[order(ord),]
#ped<-ped[-1,]
ord2<-0:(nrow(ped)-1)
ped<-cbind(ped,ord2)
ped[,1]<-as.character(ped[,1])
ped[,2]<-as.character(ped[,2])
ped[,3]<-as.character(ped[,3])

ped[is.na(ped)]<-0
dmu[!dmu$个体号%in%ped[,1],1] # nothing

phe_ind<-ped[match(dmu$个体号,ped[,1]),4]
ped_ind<-ped$ord2
sire_ind<-ped[match(ped[,2],ped[,1]),4]
dam_ind<-ped[match(ped[,3],ped[,1]),4]
ped_dmu<-data.frame(ped_ind[-1],sire_ind[-1],dam_ind[-1],ord2[-1])

dat_dmu<-data.frame(ID=phe_ind,Mu=rep(1,length(phe_ind)),Pen=dmu$测定站,
                    Sex=dmu$Sex,ys=dmu$年季,ADG=dmu$平均日增重,
                    ADFI=dmu$平均日采食量,FCR=dmu$料重比,
                    NVD=dmu$平均访问次数,TPD=dmu$平均访问时长)

write.table(ped_dmu,"DMU_PED_DUROC",col.names=F,row.names=F,quote=F)
write.table(dat_dmu,"DMU_DAT_DUROC",col.names=F,row.names=F,quote=F)

cat(paste0(paste(colnames(dat_dmu),collapse="\t"),"\n"),file="DMU_DAT_DUROC_PRO_HEADER")
