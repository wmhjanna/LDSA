suppressPackageStartupMessages(library(openair))
w<-'work/workbook_processing';a<-'outputs/workbook_processing/analysis';e<-read.csv('outputs/workbook_processing/analysis/site_eligibility.csv');counts<-list();meta<-list();sensitivity<-list();rawmeta<-list()
for(site in e$Site){
 repeat{q<-tryCatch(read.csv('outputs/workbook_processing/PNSD/QC_summary.csv'),error=function(e)NULL);if(!is.null(q)&&site %in% q$Site)break;Sys.sleep(1)}
 d<-read.csv(gzfile(file.path(w,'meta_raw',paste0(site,'.csv.gz'))));d$date<-as.POSIXct(d$Time,tz='UTC');ty<-e$Type[e$Site==site]
 for(key in c('LDSA','Total_LDSA_17_500nm','Total_LDSA_31.65_500nm','LDSA_strict_primary')){
  x<-data.frame(date=d$date,value=d[[key]]);x<-x[!is.na(x$date)&is.finite(x$value),];if(nrow(x)<3)next
  z<-as.data.frame(timeAverage(x,avg.time='month',data.thresh=10,progress=FALSE));z<-z[is.finite(z$value),c('date','value')];z$date<-format(z$date,"%Y-%m-%d");z$Site<-site;z$Type<-ty
  if(key=='LDSA'){
   counts[[site]]<-nrow(z)
   z$Metric<-'LDSA';z$Phase<-'MetaFull';z$group_id<-paste(site,'LDSA','MetaFull',sep='|');meta[[site]]<-z;rawmeta[[site]]<-z
   cat(site,'legacy 10 percent valid months',nrow(z),'\n')
  }else{
   z$Scenario<-switch(key,Total_LDSA_17_500nm='Uniform17',Total_LDSA_31.65_500nm='Uniform31.65',LDSA_strict_primary='StrictChannels');sensitivity[[paste(site,key)]]<-z
  }
 }
}
repeat{if(file.exists(file.path(w,'trend_input_base.csv')))break;Sys.sleep(1)}
e<-read.csv(file.path(a,'site_eligibility.csv'));e$valid_months<-sapply(e$Site,function(s)counts[[s]]);e$meta_eligible<-e$valid_months>=60;write.csv(e,file.path(a,'site_eligibility.csv'),row.names=FALSE)
base<-read.csv(file.path(w,'trend_input_base.csv'))
z<-do.call(rbind,meta);z<-z[names(base)];write.csv(rbind(base,z),file.path(w,'trend_input.csv'),row.names=FALSE)
write.csv(do.call(rbind,rawmeta),file.path(a,'meta_monthly_series.csv'),row.names=FALSE)
write.csv(do.call(rbind,sensitivity),file.path(w,'sensitivity_monthly_input.csv'),row.names=FALSE)
