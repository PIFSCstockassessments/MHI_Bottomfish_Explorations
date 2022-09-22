
#_____________________________________________________________________________________________________________________________
# load packages
	library(data.table)
	library(magrittr)
	library(ggplot2)
	library(ggthemes)


#_____________________________________________________________________________________________________________________________
# set working directory
	proj.dir = "D:/HOME/SAP/2024_Deep7/"


#_____________________________________________________________________________________________________________________________
# bring in indices
	mv_dt = fread(paste0(proj.dir,"VAST/model_runs/dgamma_mv7_vanilla/Index.csv")) %>%
			.[,Category:=factor(Category,levels=paste0("Category_",1:7),labels=c("prfi","etca","etco","prsi","przo","hyqu","apru"))] %>%
			.[Stratum=="Stratum_1"] %>%
			.[,Time:=(2016:2021)[as.numeric(factor(Time))]] %>%
			.[,Model:="mv"] %>%
			setnames(.,"Std. Error for ln(Estimate)","CV") %>%
			.[,.(Model,Category,Time,Estimate,CV)]

	mvpl_dt = fread(paste0(proj.dir,"VAST/model_runs/pl-dgamma_mv7_vanilla/Index.csv")) %>%
			.[,Category:=factor(Category,levels=paste0("Category_",1:7),labels=c("prfi","etca","etco","prsi","przo","hyqu","apru"))] %>%
			.[Stratum=="Stratum_1"] %>%
			.[,Time:=(2016:2021)[as.numeric(factor(Time))]] %>%
			.[,Model:="mvpl"] %>%
			setnames(.,"Std. Error for ln(Estimate)","CV") %>%
			.[,.(Model,Category,Time,Estimate,CV)]

	prfi_dt = fread(paste0(proj.dir,"VAST/model_runs/dgamma_single_prfi/Index.csv")) %>%
			.[,Category:=factor(Category,levels="Category_1",labels="prfi")] %>%
			.[Stratum=="Stratum_1"] %>%
			.[,Time:=(2016:2021)[as.numeric(factor(Time))]] %>%
			.[,Model:="single"] %>%
			setnames(.,"Std. Error for ln(Estimate)","CV") %>%
			.[,.(Model,Category,Time,Estimate,CV)]
	etca_dt = fread(paste0(proj.dir,"VAST/model_runs/dgamma_single_etca/Index.csv")) %>%
			.[,Category:=factor(Category,levels="Category_1",labels="etca")] %>%
			.[Stratum=="Stratum_1"] %>%
			.[,Time:=(2016:2021)[as.numeric(factor(Time))]] %>%
			.[,Model:="single"] %>%
			setnames(.,"Std. Error for ln(Estimate)","CV") %>%
			.[,.(Model,Category,Time,Estimate,CV)]
	etco_dt = fread(paste0(proj.dir,"VAST/model_runs/dgamma_single_etco/simple/Index.csv")) %>%
			.[,Category:=factor(Category,levels="Category_1",labels="etco")] %>%
			.[Stratum=="Stratum_1"] %>%
			.[,Time:=(2016:2021)[as.numeric(factor(Time))]] %>%
			.[,Model:="single"] %>%
			setnames(.,"Std. Error for ln(Estimate)","CV") %>%
			.[,.(Model,Category,Time,Estimate,CV)]
	prsi_dt = fread(paste0(proj.dir,"VAST/model_runs/dgamma_single_prsi/Index.csv")) %>%
			.[,Category:=factor(Category,levels="Category_1",labels="prsi")] %>%
			.[Stratum=="Stratum_1"] %>%
			.[,Time:=(2016:2021)[as.numeric(factor(Time))]] %>%
			.[,Model:="single"] %>%
			setnames(.,"Std. Error for ln(Estimate)","CV") %>%
			.[,.(Model,Category,Time,Estimate,CV)]
	przo_dt = fread(paste0(proj.dir,"VAST/model_runs/dgamma_single_przo/Index.csv")) %>%
			.[,Category:=factor(Category,levels="Category_1",labels="przo")] %>%
			.[Stratum=="Stratum_1"] %>%
			.[,Time:=(2016:2021)[as.numeric(factor(Time))]] %>%
			.[,Model:="single"] %>%
			setnames(.,"Std. Error for ln(Estimate)","CV") %>%
			.[,.(Model,Category,Time,Estimate,CV)]	
	hyqu_dt = fread(paste0(proj.dir,"VAST/model_runs/dgamma_single_hyqu/Index.csv")) %>%
			.[,Category:=factor(Category,levels="Category_1",labels="hyqu")] %>%
			.[Stratum=="Stratum_1"] %>%
			.[,Time:=(2016:2021)[as.numeric(factor(Time))]] %>%
			.[,Model:="single"] %>%
			setnames(.,"Std. Error for ln(Estimate)","CV") %>%
			.[,.(Model,Category,Time,Estimate,CV)]
	apru_dt = fread(paste0(proj.dir,"VAST/model_runs/dgamma_single_apru/Index.csv")) %>%
			.[,Category:=factor(Category,levels="Category_1",labels="apru")] %>%
			.[Stratum=="Stratum_1"] %>%
			.[,Time:=(2016:2021)[as.numeric(factor(Time))]] %>%
			.[,Model:="single"] %>%
			setnames(.,"Std. Error for ln(Estimate)","CV") %>%
			.[,.(Model,Category,Time,Estimate,CV)]
	apru_nof_dt = fread(paste0(proj.dir,"VAST/model_runs/dgamma_single_apru_nofilter/Index.csv")) %>%
			.[,Category:=factor(Category,levels="Category_1",labels="apru")] %>%
			.[Stratum=="Stratum_1"] %>%
			.[,Time:=(2016:2021)[as.numeric(factor(Time))]] %>%
			.[,Model:="single - no filter"] %>%
			setnames(.,"Std. Error for ln(Estimate)","CV") %>%
			.[,.(Model,Category,Time,Estimate,CV)]

	mv_agg_dt = fread(paste0(proj.dir,"VAST/model_runs/dgamma_mv7_vanilla/agg/Index.csv")) %>%
			.[Category=="Category_7"] %>%
			.[,Category:="all"] %>%
			.[Stratum=="Stratum_1"] %>%
			.[,Time:=(2016:2021)[as.numeric(factor(Time))]] %>%
			.[,Model:="mv"] %>%
			setnames(.,"Std. Error for ln(Estimate)","CV") %>%
			.[,.(Model,Category,Time,Estimate,CV)]

	mvpl_agg_dt = fread(paste0(proj.dir,"VAST/model_runs/pl-dgamma_mv7_vanilla/agg/Index.csv")) %>%
			.[Category=="Category_7"] %>%
			.[,Category:="all"] %>%
			.[Stratum=="Stratum_1"] %>%
			.[,Time:=(2016:2021)[as.numeric(factor(Time))]] %>%
			.[,Model:="mvpl"] %>%
			setnames(.,"Std. Error for ln(Estimate)","CV") %>%
			.[,.(Model,Category,Time,Estimate,CV)]

index_dt = rbind(mv_dt,prfi_dt,etca_dt,etco_dt,prsi_dt,przo_dt,hyqu_dt,apru_dt,apru_nof_dt)

p1 = index_dt %>%
	ggplot() +
	ylim(0,NA) +
	facet_wrap(~Category) +
	geom_hline(yintercept=0) +
	geom_path(aes(x=Time,y=Estimate,group=Model,color=Model)) +
	theme_few()
p1

p2 = index_dt %>%
	ggplot() +
	ylim(0,NA) +
	facet_wrap(~Category,scales="free_y") +
	geom_hline(yintercept=0) +
	geom_path(aes(x=Time,y=Estimate,group=Model,color=Model)) +
	theme_few()
p2

p3 =  index_dt %>%
	 .[Model!="single - no filter",.(Estimate=sum(Estimate)),by=.(Model,Time)] %>%
	 ggplot() +
	 ylim(0,NA) +
	 geom_hline(yintercept=0) +
	 geom_path(aes(x=Time,y=Estimate,group=Model,color=Model)) +
	 theme_few()
p3


p4 =  rbind(mv_agg_dt,mvpl_agg_dt) %>% 
	  .[,l95:=exp(log(Estimate)-2*CV)] %>%
	  .[,u95:=exp(log(Estimate)+2*CV)] %>%
	 ggplot() +
	 ylim(0,NA) +
	 geom_hline(yintercept=0) +
	 geom_ribbon(aes(x=Time,ymin=l95,ymax=u95,group=Model,fill=Model),alpha=0.25) +
	 geom_path(aes(x=Time,y=Estimate,group=Model,color=Model)) +
	 theme_few()
p4


p5 = rbind(mv_dt,mvpl_dt) %>% 
	  .[,l95:=exp(log(Estimate)-2*CV)] %>%
	  .[,u95:=exp(log(Estimate)+2*CV)] %>%
	ggplot() +
	ylim(0,NA) +
	facet_wrap(~Category,scales="free_y") +
	geom_hline(yintercept=0) +
	geom_ribbon(aes(x=Time,ymin=l95,ymax=u95,group=Model,fill=Model),alpha=0.25) +
	geom_path(aes(x=Time,y=Estimate,group=Model,color=Model)) +
	theme_few()
p5

