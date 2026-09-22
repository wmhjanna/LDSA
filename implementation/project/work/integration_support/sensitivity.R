suppressPackageStartupMessages({library(openair);library(meta)})
w<-'work/integration_support';a<-'outputs/integration_support/analysis'
d<-read.csv(file.path(w,'sensitivity_monthly_input.csv'));d$date<-as.POSIXct(d$date,tz='UTC')
x<-read.csv(file.path(a,'meta_sites.csv'));rows<-list()
for(sc in unique(d$Scenario))for(site in x$Site){
 g<-d[d$Scenario==sc & d$Site==site,];if(nrow(g)<60)next
 set.seed(20260909+sum(utf8ToInt(paste(site,'LDSA','MetaFull',sep='|'))))
 fit<-TheilSen(data.frame(date=g$date,value=g$value),pollutant='value',avg.time='month',deseason=FALSE,slope.percent=TRUE,plot=FALSE,silent=TRUE)
 r<-fit$data$res2[1,];rows[[length(rows)+1]]<-data.frame(Site=site,Type=g$Type[1],Scenario=sc,n=nrow(g),effect=r$slope.percent,low=r$lower.percent,high=r$upper.percent)
}
ss<-do.call(rbind,rows);ss$SE<-(ss$high-ss$low)/(2*qnorm(.975));ss<-ss[is.finite(ss$SE)&ss$SE>0,];write.csv(ss,file.path(a,'sensitivity_site_trends.csv'),row.names=FALSE)
results<-list()
run<-function(g,label,method='classic'){
 fit<-metagen(TE=g$effect,seTE=g$SE,studlab=g$Site,common=FALSE,random=TRUE,method.tau='REML',method.random.ci=method,adhoc.hakn.ci='se')
 data.frame(Scenario=label,n=nrow(g),effect=fit$TE.random,lower=fit$lower.random,upper=fit$upper.random,I2=fit$I2,tau2=fit$tau2)
}
x$effect<-x$Slope_Percent;x$low<-x$CI_Lower_Percent;x$high<-x$CI_Upper_Percent
results[[1]]<-run(x,'Primary');results[[2]]<-run(x,'Primary modified HK','HK');results[[3]]<-run(x[x$Site!='LND_TR',],'Exclude LND_TR')
for(sc in unique(ss$Scenario)){
 g<-ss[ss$Scenario==sc,];results[[length(results)+1]]<-run(g,paste(sc,'same eligible primary sites'))
 if(nrow(g)<nrow(x))results[[length(results)+1]]<-run(x[x$Site %in% g$Site,],paste('Primary matched to',sc))
}
write.csv(do.call(rbind,results),file.path(a,'meta_sensitivity.csv'),row.names=FALSE)
loo<-list();for(site in x$Site)loo[[site]]<-run(x[x$Site!=site,],paste('Exclude',site))
write.csv(do.call(rbind,loo),file.path(a,'meta_leave_one_out.csv'),row.names=FALSE)
print(do.call(rbind,results))
