suppressPackageStartupMessages(library(openair))
w<-'work/pal_channel_conversion';a<-'outputs/pal_channel_conversion/analysis'
d<-read.csv(gzfile(file.path(w,'PAL_RB.csv.gz')));d$date<-as.POSIXct(d$date,tz='UTC')
z<-as.data.frame(timeAverage(d,avg.time='month',data.thresh=10,progress=FALSE));z<-z[is.finite(z$value),c('date','value')]
z$date<-format(z$date,'%Y-%m-%d');z$Site<-'PAL_RB';z$Type<-'RB';z$Metric<-'LDSA';z$Phase<-'MetaFull';z$group_id<-'PAL_RB|LDSA|MetaFull'
e<-read.csv(file.path(a,'site_eligibility.csv'));e$valid_months[e$Site=='PAL_RB']<-nrow(z);e$meta_eligible[e$Site=='PAL_RB']<-nrow(z)>=60;write.csv(e,file.path(a,'site_eligibility.csv'),row.names=FALSE)
base<-read.csv(file.path(w,'trend_input_base.csv'));write.csv(rbind(base,z[names(base)]),file.path(w,'trend_input.csv'),row.names=FALSE)
old<-read.csv(file.path(a,'meta_monthly_series.csv'));write.csv(rbind(old[old$Site!='PAL_RB',],z[names(old)]),file.path(a,'meta_monthly_series.csv'),row.names=FALSE)
cat('PAL Meta months:',nrow(z),'\n')
