suppressPackageStartupMessages(library(openair))
w<-'work/network_trends';a<-'outputs/network_trends/analysis'
e<-read.csv(file.path(a,'site_eligibility.csv'));parts<-list()
for (p in list.files(file.path(w,'raw'),full.names=TRUE)) {
 site<-sub('.csv.gz','',basename(p),fixed=TRUE);d<-read.csv(gzfile(p));d$date<-as.POSIXct(d$date,tz='UTC')
 z<-as.data.frame(timeAverage(d,avg.time='month',data.thresh=10,progress=FALSE));z<-z[is.finite(z$value),c('date','value')]
 z$date<-format(z$date,'%Y-%m-%d');z$Site<-site;z$Type<-e$Type[e$Site==site];z$Metric<-'LDSA';z$Phase<-'MetaFull';z$group_id<-paste(site,'LDSA','MetaFull',sep='|');parts[[site]]<-z
 e$valid_months[e$Site==site]<-nrow(z);e$meta_eligible[e$Site==site]<-nrow(z)>=60
 cat(site,'Meta months',nrow(z),'eligible',nrow(z)>=60,'\n')
}
write.csv(e,file.path(a,'site_eligibility.csv'),row.names=FALSE)
base<-read.csv(file.path(w,'trend_input_base.csv'));extra<-do.call(rbind,parts);extra<-extra[names(base)]
write.csv(rbind(base,extra),file.path(w,'trend_input.csv'),row.names=FALSE)
old<-read.csv('outputs/workbook_processing/analysis/meta_monthly_series.csv')
write.csv(rbind(old,extra[names(old)]),file.path(a,'meta_monthly_series.csv'),row.names=FALSE)
