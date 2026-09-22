suppressPackageStartupMessages(library(meta))
a<-'outputs/record_exclusions/analysis'
x<-read.csv(file.path(a,'meta_sites.csv'));x$effect<-x$Slope_Percent
coverage<-jsonlite::fromJSON('work/record_exclusions/spatial/coverage_selection.json')
run<-function(g,label,method='classic'){
 fit<-metagen(TE=g$effect,seTE=g$SE,studlab=g$Site,common=FALSE,random=TRUE,method.tau='REML',method.random.ci=method,adhoc.hakn.ci='se')
 data.frame(Scenario=label,n=nrow(g),effect=fit$TE.random,lower=fit$lower.random,upper=fit$upper.random,I2=fit$I2,tau2=fit$tau2)
}
rows<-list(run(x,'Primary 65-site network'),run(x,'Primary modified HK','HK'),run(x[x$Site!='LND_TR',],'Exclude LND_TR'),run(x[x$Site %in% coverage$keep,],'Strict 49-site spatial-coverage subset, eligible Meta sites'),run(x[!x$Site %in% coverage$upper_excluded,],'Upper-complete 53-site spatial-coverage subset, eligible Meta sites'))
write.csv(do.call(rbind,rows),file.path(a,'meta_sensitivity.csv'),row.names=FALSE)
loo<-list();for(site in x$Site)loo[[site]]<-run(x[x$Site!=site,],paste('Exclude',site))
write.csv(do.call(rbind,loo),file.path(a,'meta_leave_one_out.csv'),row.names=FALSE)
print(do.call(rbind,rows))
