suppressPackageStartupMessages(library(openair))
d<-read.csv('@PROJECT_ROOT@/work/network_analysis/GIO_raw.csv');d$date<-as.POSIXct(d$date,tz='UTC')
z<-as.data.frame(timeAverage(d,avg.time='month',data.thresh=10,progress=FALSE));z<-z[is.finite(z$value),c('date','value')];z$date<-format(z$date,'%Y-%m-%d');z$Site<-'GIO_RB';z$Type<-'RB';z$Metric<-'LDSA';z$Phase<-'MetaFull';z$group_id<-'GIO_RB|LDSA|MetaFull';write.csv(z,'@PROJECT_ROOT@/work/network_analysis/gio_meta.csv',row.names=FALSE)
