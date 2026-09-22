args <- commandArgs(trailingOnly=TRUE)
suppressPackageStartupMessages({library(openair);library(meta)})
d <- read.csv(args[1]); d$date <- as.POSIXct(d$date,tz='UTC')
groups <- split(d,d$group_id); out <- list(); errors <- list()
prior<-file.path(args[2],'trends.csv')
if(file.exists(prior)){
 old<-read.csv(prior);old$percent_CI_valid<-NULL;keys<-paste(old$Site,old$Metric,old$Phase,sep='|');out<-split(old,keys)
}
out<-out[intersect(names(out),names(groups))]
for (id in names(groups)) {
  if(id %in% names(out))next
  set.seed(20260909+sum(utf8ToInt(id)))
  g <- groups[[id]]
  fit <- tryCatch(TheilSen(data.frame(date=g$date,value=g$value),pollutant='value',avg.time='month',deseason=FALSE,slope.percent=TRUE,plot=FALSE,silent=TRUE),error=function(e)e)
  if(inherits(fit,'error')) {errors[[id]]<-data.frame(group_id=id,error=conditionMessage(fit));next}
  r <- fit$data$res2[1,]
  value<-function(name)if(name %in% names(r) && length(r[[name]]))as.numeric(r[[name]][1]) else NA_real_
  out[[id]]<-data.frame(Site=g$Site[1],Type=g$Type[1],Metric=g$Metric[1],Phase=g$Phase[1],n_months=nrow(g),Slope_Unit=value('slope'),CI_Lower_Unit=value('lower'),CI_Upper_Unit=value('upper'),Slope_Percent=value('slope.percent'),CI_Lower_Percent=value('lower.percent'),CI_Upper_Percent=value('upper.percent'),P_value=value('p'),Avg_Conc=value('conc'))
  if(length(out)%%50==0){cat('Fitted',length(out),'of',length(groups),'\n');write.csv(do.call(rbind,out),prior,row.names=FALSE)}
}
result<-do.call(rbind,out)
bad<-is.finite(result$CI_Lower_Percent)&is.finite(result$CI_Upper_Percent)&result$CI_Lower_Percent>result$CI_Upper_Percent
rawbad<-result[bad,]
rawpath<-file.path(args[2],'invalid_percent_CI_raw.csv')
if(file.exists(rawpath)){
 previousbad<-read.csv(rawpath);previousbad$percent_CI_valid<-NULL
 keys<-paste(result$Site,result$Metric,result$Phase,sep='|');oldkeys<-paste(previousbad$Site,previousbad$Metric,previousbad$Phase,sep='|')
 keep<-oldkeys %in% keys[!is.finite(result$CI_Lower_Percent)&!is.finite(result$CI_Upper_Percent)]
 rawbad<-rbind(rawbad,previousbad[keep,names(rawbad)])
}
write.csv(rawbad,rawpath,row.names=FALSE)
result$percent_CI_valid<-is.finite(result$CI_Lower_Percent)&is.finite(result$CI_Upper_Percent)&!bad
result[bad,c('CI_Lower_Percent','CI_Upper_Percent')]<-NA
write.csv(result,file.path(args[2],'trends.csv'),row.names=FALSE)
if(length(errors))write.csv(do.call(rbind,errors),file.path(args[2],'trend_errors.csv'),row.names=FALSE)
e<-read.csv(file.path(args[2],'site_eligibility.csv'))
e$meta_eligible<-tolower(as.character(e$meta_eligible)) %in% c('true','1')
x<-result[result$Metric=='LDSA' & result$Phase=='MetaFull',];x<-merge(x,e[c('Site','meta_eligible')],by='Site')
x$SE<-(x$CI_Upper_Percent-x$CI_Lower_Percent)/(2*qnorm(.975))
x<-x[x$meta_eligible & is.finite(x$SE)&x$SE>0 & is.finite(x$Slope_Percent),]
if(nrow(x)>=2){
 m<-metagen(TE=x$Slope_Percent,seTE=x$SE,studlab=x$Site,subgroup=x$Type,common=FALSE,random=TRUE,method.tau='REML',method.random.ci='classic')
 x$weight_percent<-m$w.random/sum(m$w.random)*100
 write.csv(x,file.path(args[2],'meta_sites.csv'),row.names=FALSE)
 write.csv(data.frame(Type=c(m$subgroup.levels,'Overall'),n=c(m$k.w,m$k),effect=c(m$TE.random.w,m$TE.random),lower=c(m$lower.random.w,m$lower.random),upper=c(m$upper.random.w,m$upper.random),tau2=c(m$tau2.w,m$tau2)),file.path(args[2],'meta_summary.csv'),row.names=FALSE)
 write.csv(data.frame(Q=m$Q.b.random,df=m$df.Q.b.random,p=m$pval.Q.b.random),file.path(args[2],'meta_subgroup_test.csv'),row.names=FALSE)
}
z<-read.csv(file.path(args[2],'meta_summary.csv'));z$estimate_kind<-ifelse(z$n==1,'single site','random effects')
for(t in z$Type[z$n==1]){r<-x[x$Type==t,];z[z$Type==t,c('effect','lower','upper','tau2')]<-list(r$Slope_Percent,r$CI_Lower_Percent,r$CI_Upper_Percent,NA)}
write.csv(z,file.path(args[2],'meta_summary.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(args[2],'R_session_info.txt'))
