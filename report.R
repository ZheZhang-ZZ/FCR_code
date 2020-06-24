######FIRE数据分析报告
library(ggplot2)
setwd('/Users/zhezhang/Desktop/天邦数据/FI/guigang')
#load("rfi.RData")
theme <-theme_get()
theme$text$family <- "STSong"
theme_set(theme)

## 原始数据条数
nrow(data1) 
time = as.POSIXct(as.character(data1$进入时间),format = "%Y/%m/%d %H:%M")
summary(time)
# Min.               1st Qu.                Median                  Mean               3rd Qu. 
# "2019-10-15 11:25:00" "2019-12-12 09:08:00" "2020-01-09 20:49:00" "2020-01-06 21:28:29" "2020-02-03 07:39:00" 
# Max.                  NA's 
# "2020-03-14 23:59:00"                  "36" 
time = as.POSIXct(as.character(data1$退出时间),format = "%Y/%m/%d %H:%M")
summary(time)
# Min.               1st Qu.                Median                  Mean               3rd Qu. 
# "2019-10-15 11:34:00" "2019-12-12 09:11:00" "2020-01-09 20:51:00" "2020-01-06 21:30:29" "2020-02-03 07:41:00" 
# Max.                  NA's 
# "2020-03-15 00:24:00"                  "22" 

## 后续能用的数据条数
nrow(data)  
table(data$品种品系)

## 个体总数
nrow(fcr)

## 各个品种个体数
table(fcr$品种品系)
table(fcr$品种品系,fcr$性别)

## 采食相关数据错误
sapply(paste0("et",1:16),function(x)length(get(x)))
###22号测定站反常
sapply(paste0("et",1:16),function(x)sum(data[get(x),'测定站']==22))
# et1  et2  et3  et4  et5  et6  et7  et8  et9 et10 et11 et12 et13 et14 et15 et16 
# 0    1    0    0 2150  158    0  580    5  211    2  125    2  125 4157 4157 

data_err<-data[err_all,]
plot_station<-table(data_err$测定站)/table(data$测定站)
plot_station<-data.frame(plot_station)
colnames(plot_station)<-c('测定站','错误比例')
png(file="测定站错误数据分布.png",width=12, height=8,units='in',family="GB1",res=300)
ggplot(plot_station, aes(x=测定站, y=错误比例)) + geom_bar(stat="identity")+
  geom_hline(yintercept = nrow(data_err)/nrow(data),colour="red")+
  theme(text = element_text(size=15))

#barplot(table(data_err$测定站)/table(data$测定站),ylim=c(0,0.08),xlab="测定站",ylab="错误比例")
# abline(h=nrow(data_err)/nrow(data),col=4,lwd=2)

dev.off()

###射频耳标识别失败
data_fail<-droplevels(data1[data1$射频耳标==0,])
plot_station<-table(data_fail$测定站)/table(data1$测定站)
plot_station<-data.frame(plot_station)
colnames(plot_station)<-c('测定站','射频耳标识别失败比例')
png(file="测定站射频耳标识别错误数据分布.png",width=12, height=8,units='in',family="GB1",res=300)
ggplot(plot_station, aes(x=测定站, y=射频耳标识别失败比例)) + geom_bar(stat="identity")+
  geom_hline(yintercept = nrow(data_fail)/nrow(data1),colour="red")+
  theme(text = element_text(size=15))


dev.off()

####开测日龄分布
ord<-order(data$个体号,data$测定日龄,decreasing=c(FALSE,TRUE),method='radix')
data2<-data[ord,]
offset_days<-data2[!duplicated(data2$个体号),c('个体号','测定日龄','性别','品种品系')]
offset_days%>%dplyr::group_by(品种品系,性别)%>%dplyr::summarise(数量=dplyr::n())
as.data.frame(offset_days%>%dplyr::group_by(品种品系)%>%dplyr::summarise(平均测定日龄=mean(测定日龄)))
as.data.frame(offset_days%>%dplyr::group_by(品种品系)%>%dplyr::summarise(sd=2*sd(测定日龄)))

ord<-order(data$个体号,data$测定日龄,decreasing=c(FALSE,FALSE),method='radix')
data2<-data[ord,]
onset_days<-data2[!duplicated(data2$个体号),c('个体号','测定日龄','性别','品种品系')]
as.data.frame(onset_days%>%dplyr::group_by(品种品系)%>%dplyr::summarise(平均测定日龄=mean(测定日龄)))
as.data.frame(onset_days%>%dplyr::group_by(品种品系)%>%dplyr::summarise(sd=2*sd(测定日龄)))


png(file="距离119日龄的距离2.png",width=9, height=5,units='in',res=300)
par(family='STSong')
barplot(table(onset_time2$开测日龄-119),
        xlab = "距离119日龄的距离(天)", ylab = "频数") ##看一下这组数据的开测日龄距离119的距离
dev.off()

## 几种典型的体重错误记录
plot_bw<-function(a){
  ggplot(a, aes(x=测定日龄, y=当天体重)) + geom_point()+
    theme(text = element_text(size=15))

}
plot_bw2<-function(a){
  ggplot(a, aes(x=测定日龄, y=当天体重)) +
    facet_wrap(~ 个体号, nrow=2) + geom_point()+
    theme(text = element_text(size=15))
  
}
#比较正常的
png(file="体重_较正常.png",width=12, height=8,units='in',family="GB1",res=300)

plot_bw(droplevels(clst_res1[clst_res1$个体号=="DDTGXGG19501649"]))

dev.off()

#不正常的
abnorm<-clst_res1[clst_res1$个体号=="PPTGXGG19403548" | 
clst_res1$个体号=="LLTGXGG19310376" | 
  clst_res1$个体号=="PPTGXGG19403137" |
  clst_res1$个体号=="DDTGXGG19206456"]
abnorm$个体号<-factor(abnorm$个体号,levels=c("PPTGXGG19403548",
                                       "LLTGXGG19310376",
                                       "PPTGXGG19403137",
                                       "DDTGXGG19206456"))

png(file="体重(校正之前).png",width=12, height=8,units='in',family="GB1",res=300)

plot_bw2(abnorm)

dev.off()

png(file="体重(校正之后).png",width=12, height=8,units='in',family="GB1",res=300)

plot_flexmix2(abnorm)

dev.off()


## 开测日龄
summary(onset_time$测定日龄)
barplot(table(onset_time$测定日龄))
mean(onset_time$测定日龄)
sd(onset_time$测定日龄)

by(onset_time,list(onset_time$品种品系),function(x)mean(x$测定日龄))
by(onset_time,list(onset_time$品种品系),function(x)sd(x$测定日龄))

## 结测日龄
summary(offset_date$测定日龄)
barplot(table(offset_date$测定日龄))
mean(offset_date$测定日龄)
sd(offset_date$测定日龄)

by(offset_date,list(onset_time$品种品系),function(x)mean(x$测定日龄))
by(offset_date,list(onset_time$品种品系),function(x)sd(x$测定日龄))

## 开测日龄与结测日龄差异较大，统一校正到119和168日龄

# weight119<-weight[weight$测定日龄==119,]
# weight168<-weight[weight$测定日龄==168,]
bw119<-bw119_168[bw119_168$测定日龄==119,]
bw168<-bw119_168[bw119_168$测定日龄==168,]

# bw119[!is.na(match(bw119$个体号,weight119$个体号)),'预测体重']<-weight119$预测体重
# bw168[!is.na(match(bw168$个体号,weight168$个体号)),'预测体重']<-weight168$预测体重

bw119<-cbind(bw119,fcr$品种品系)
bw168<-cbind(bw168,fcr$品种品系)

colnames(bw119)[ncol(bw119)]<-'品种品系'
colnames(bw168)[ncol(bw168)]<-'品种品系'

plot_dat<-rbind(bw119,bw168)
plot_dat$品种品系<-mapvalues(plot_dat$品种品系, from = c("美系大白", "美系杜洛克", "美系皮特兰", "EB5杜洛克", "美系长白"), 
                       to = c("大白", "杜洛克", "皮特兰", "EB5", "长白"))
plot_dat$品种品系<-factor(plot_dat$品种品系,levels=c("大白", "长白", "杜洛克", "EB5", "皮特兰"))
png(file="119—168体重.png",width=12, height=8,units='in',family="GB1",res=300)

# ggplot(plot_dat, aes(y=预测体重,x=品种品系,fill=品种品系)) + geom_boxplot()+ 
#   theme(legend.position = "none")+ylab("体重")+
#   facet_wrap(~ 测定日龄, nrow=2,scales='free_y')+
#   theme(text = element_text(size=15))

ggplot(plot_dat, aes(x = 品种品系, y = 预测体重))+geom_violin(aes(fill = 品种品系), trim = FALSE) + 
  geom_boxplot(width = 0.2)+facet_wrap( ~ 测定日龄, nrow=2,scales='free_y')+
  scale_fill_manual(values = c("#00AFBB", "#E7B800","#636363","#FC4E07","#3182bd"))+
  theme(legend.position = "none",text = element_text(size=15))

dev.off()


cor(bw119[bw119$品种品系=='美系大白']$预测体重,
    bw168[bw168$品种品系=='美系大白']$预测体重,
    method = 'spearman') #0.7521227
cor(bw119[bw119$品种品系=='美系杜洛克']$预测体重,
    bw168[bw168$品种品系=='美系杜洛克']$预测体重,
    method = 'spearman') #0.7831125
cor(bw119[bw119$品种品系=='美系皮特兰']$预测体重,
    bw168[bw168$品种品系=='美系皮特兰']$预测体重,
    method = 'spearman') #0.8068745
cor(bw119[bw119$品种品系=='美系长白']$预测体重,
    bw168[bw168$品种品系=='美系长白']$预测体重,
    method = 'spearman') #0.7939808
cor(bw119[bw119$品种品系=='EB5杜洛克']$预测体重,
    bw168[bw168$品种品系=='EB5杜洛克']$预测体重,
    method = 'spearman') #0.803428

####出现体重分层的数据
table(clst_res2$k)
multi_lay<-clst_res2[clst_res2$k>1,]
fcr_ml<-fcr[match(multi_lay$个体号,fcr$个体号),]
ml_dat<-data.frame(测定站=c(as.numeric(names(table(fcr$测定站))),
                         as.numeric(names(table(fcr_ml$测定站)))),
                   个体数=c(table(fcr$测定站)-table(fcr_ml$测定站),table(fcr_ml$测定站)),
                   组=c(rep("总数",length(table(fcr$测定站))),rep("分层个体数",length(table(fcr$测定站)))),
                   比例=rep(table(fcr_ml$测定站)/table(fcr$测定站),2))
ml_dat$组<-factor(ml_dat$组,levels=c('总数','分层个体数'))
png(file="体重分层.png",width=12, height=8,units='in',family="GB1",res=300)

ggplot(ml_dat, aes(x=测定站)) +
  geom_bar(aes(y=个体数, fill=组),stat="identity") +
  geom_line(aes(y = 比例 * 60, group = 1, color = '比例')) +
  geom_point(aes(y = 比例*60, color='比例'))+
  scale_fill_manual(values = c(总数 = "#EFC000FF", 分层个体数 = "#0073C2FF"), name="") +
  scale_colour_manual(values = c(比例 = "#bcbddc"), name="") +
  scale_y_continuous(sec.axis = sec_axis(trans = ~ . / 60)) +theme_bw(base_family="STSong", base_size = 15) +
  theme(legend.title = element_blank(),text = element_text(size=15))
dev.off()

### 性状的描述统计

combine<-data.frame(个体号=fcr$个体号,品种品系=fcr$品种品系,
                       性别=fcr$性别,测定站=fcr$测定站,同期组=fcr$同期组,
                       达115公斤体重日龄=fcr$达115公斤体重日龄,
                       平均日采食量=fcr$平均日采食量,
                       平均日增重=fcr$平均日增重,料重比=fcr$料重比,
                       次平均采食量=fcr$次平均采食量,次平均采食速度=fcr$次平均采食速度,
                       次平均采食时间=fcr$次平均采食时间,日平均采食次数=fcr$日平均采食次数)
size<-nrow(combine)
ntrait<-8
plot_dat<-data.frame(品种品系=rep(combine$品种品系,size*ntrait),
                     性状=rep(c('达115公斤体重日龄','平均日增重','平均日采食量',
                              '料重比','次平均采食量','次平均采食速度',
                              '次平均采食时间','日平均采食次数'),each=size),
                     值=c(combine$达115公斤体重日龄,combine$平均日增重,
                         combine$平均日采食量,combine$料重比,
                         combine$次平均采食量,combine$次平均采食速度,
                         combine$次平均采食速度,combine$日平均采食次数)
                     
                     )
plot_dat$性状<-factor(plot_dat$性状,levels=c('达115公斤体重日龄','平均日增重','平均日采食量',
                                         '料重比','次平均采食量','次平均采食速度',
                                         '次平均采食时间','日平均采食次数'))
plot_dat$品种品系<-mapvalues(plot_dat$品种品系, from = c("美系大白", "美系杜洛克",
                                                 "美系皮特兰", "EB5杜洛克", "美系长白"), 
                            to = c("大白", "杜洛克", "皮特兰", "EB5", "长白"))
plot_dat$品种品系<-factor(plot_dat$品种品系,levels=c("大白", "长白", "杜洛克", "EB5", "皮特兰"))
png(file="性状描述性统计.png",width=12, height=8,units='in',family="GB1",res=300)
par(family='STSong')
ggplot(plot_dat, aes(x = 品种品系, y = 值))+geom_violin(aes(fill = 品种品系), trim = FALSE) + 
  geom_boxplot(width = 0.2)+facet_wrap( ~ 性状,nrow=4,scales='free_y')+
  scale_fill_manual(values = c("#00AFBB", "#E7B800","#636363","#FC4E07","#3182bd"))+
  theme_bw(base_family="STSong", base_size = 15) +
  theme(legend.position = "none",text = element_text(size=15))
dev.off()

avtrait_mean<-combine %>% 
  group_by(品种品系) %>% 
  summarise(
    平均日采食量 = mean(平均日采食量),
    平均日增重 = mean(平均日增重),
    料重比 = mean(料重比),
    次平均采食量 = mean(次平均采食量),
    次平均采食速度 = mean(次平均采食速度),
    次平均采食时间 = mean(次平均采食时间),
    日平均采食次数 = mean(日平均采食次数),
  )

avtrait_se<-combine %>% 
  group_by(品种品系) %>% 
  summarise(
    平均日采食量 = sd(平均日采食量),
    平均日增重= sd(平均日增重),
    料重比 = sd(料重比),
    次平均采食量 = sd(次平均采食量),
    次平均采食速度 = sd(次平均采食速度),
    次平均采食时间 = sd(次平均采食时间),
    日平均采食次数 = sd(日平均采食次数)
  )

avtrait_mean<-as.data.frame(avtrait_mean)
avtrait_se<-as.data.frame(avtrait_se)
row.names(avtrait_mean)<-avtrait_mean[,1]
row.names(avtrait_se)<-avtrait_se[,1]
breed_name<-avtrait_se[,1]
trait_name<-colnames(avtrait_se)[-1]
avtrait_mean<-round(as.numeric(t(avtrait_mean)[-1,]),3)
avtrait_se<-round(as.numeric(t(avtrait_se)[-1,]),3)

avtrait<-paste(avtrait_mean,avtrait_se,sep="±")
avtrait<-matrix(avtrait,nrow=7)
rownames(avtrait)<-trait_name
colnames(avtrait)<-breed_name
write.csv(avtrait,file="平均性状描述统计.csv",quote=F)

adfi_summary<-input %>% 
  group_by(测定日龄,品种品系) %>% 
  summarise(
    日采食量 = mean(校正采食量/1000),
    标准差 = sd(校正采食量/1000)
  )

png(file="采食量描述性统计.png",width=12, height=8,units='in',family="GB1",res=300)

adfi_summary<-droplevels(adfi_summary[complete.cases(adfi_summary),])
ggplot(adfi_summary, aes(x=测定日龄,  y=日采食量, ymin=日采食量-标准差, ymax=日采食量+标准差,fill=品种品系)) +
  geom_bar(stat = "identity") +facet_wrap( ~ 品种品系,nrow=5)+
  geom_errorbar(width=.03, position=position_dodge(0.045))+
  theme(legend.position = 'none',text = element_text(size=15))

dev.off()

#采食高峰
enter_time<-data.frame(品种品系=data$品种品系,
                       进入时间=as.numeric(format(strptime(data$进入时间,"%Y-%m-%d %H:%M:%S"),'%H')))

enter_time_freq<-enter_time %>% 
  group_by(品种品系,进入时间)%>%
  summarise(n=n())%>%
  mutate(频率 = n / sum(n))

png(file="采食量描述性统计.png",width=12, height=8,units='in',family="GB1",res=300)
ggplot(enter_time_freq, aes(x=进入时间,y=频率,fill=品种品系)) +
  geom_bar(stat = "identity",position='dodge') +
  scale_fill_manual(values = c("#00AFBB", "#E7B800","#636363","#FC4E07","#3182bd"))+
  theme_pubclean(base_family="STSong", base_size = 15)+
  scale_x_discrete(limits=0:23,breaks=as.character(0:23),labels=paste(0:23,1:24,sep = "-"))+
  theme(legend.title=element_blank())
dev.off()

