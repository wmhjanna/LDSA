suppressPackageStartupMessages(library(openair));w<-'work/workbook_processing';a<-'outputs/workbook_processing/analysis';all<-read.csv(file.path(a,'full_monthly_series.csv'));rows<-list()
for(site in c('BCN_UB','LND_TR','ATH_SUB')){
 d<-read.csv(gzfile(file.path(w,'meta_raw',paste0(site,'.csv.gz'))));d<-data.frame(date=as.POSIXct(d$Time,tz='UTC'),value=d$LDSA);d<-d[is.finite(d$value)&!is.na(d$date),];d<-d[order(d$date),];d<-d[!duplicated(d$date),]
 x<-as.data.frame(timeAverage(d,avg.time='month',data.thresh=0,progress=FALSE));x$date<-format(x$date,'%Y-%m-%d');y<-all[all$Site==site&all$Metric=='LDSA',];z<-merge(x,y,by='date');err<-max(abs(z$value.x-z$value.y),na.rm=TRUE);stopifnot(err<1e-8,sum(is.finite(x$value))==nrow(y));rows[[site]]<-data.frame(Site=site,months=nrow(y),max_absolute_difference=err,passed=TRUE)
}
write.csv(do.call(rbind,rows),file.path(a,'reference_monthly_aggregation_verification.csv'),row.names=FALSE);print(do.call(rbind,rows))
