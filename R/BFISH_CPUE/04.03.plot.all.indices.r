

# Nicholas Ducharme-Barth
# 2023/06/04
# Comparison plots

# Copyright (c) 2023 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#_____________________________________________________________________________________________________________________________
# load packages
	library(data.table)
	library(magrittr)
	library(ggplot2)
	library(ggthemes)


#_____________________________________________________________________________________________________________________________
# set working directory
	proj_dir = "D:/HOME/SAP/2024_Deep7/"
	dir_plot = paste0(proj_dir,"VAST/model_runs/comparison_plots/")
	dir.create(dir_plot,recursive = TRUE)

    deep7_code_vec = tolower(c("ETCA","APRU","PRSI","HYQU","PRFI","PRZO","ETCO"))
	deep7_name_vec = c("Ehu (ETCA)", "Lehi (APRU)", "Kalekale (PRSI)", "Hapu'upu'u (HYQU)", "'Opakapaka (PRFI)", "Gindai (PRZO)", "Onaga (ETCO)")

	dir_vec = dir(paste0(proj_dir,"VAST/model_runs/"),recursive=TRUE)
	ests_vec = dir_vec[grep("/index_dt.csv",dir_vec,fixed=TRUE)]
  	ests_vec = ests_vec[grep("/2022_",ests_vec,fixed=TRUE)]	

#_____________________________________________________________________________________________________________________________
# bring in indices
    # design based bfish index data is in kg so convert to lbs (millions)
	design_dt = fread(paste0(proj_dir,"Data/2022_design_based_estimates.csv")) %>%
				.[,.(Model,Category,Time,Estimate,CV)] %>%
                .[,Estimate:=Estimate*2.20462262185] %>%
				.[,Estimate:=Estimate/1000000] %>%
				.[,l95:=exp(log(Estimate)-2*sqrt(log(CV^2+1)))] %>%
				.[,u95:=exp(log(Estimate)+2*sqrt(log(CV^2+1)))] %>%
                .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                .[is.na(Category),Category:="Total"] %>%
                .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]

#_____________________________________________________________________________________________________________________________
# extract other indices
  # order by last modified
  ests_vec = ests_vec[order(sapply(paste0(proj_dir,"VAST/model_runs/",ests_vec),file.mtime))]

  index_dt.list = as.list(rep(NA,length(ests_vec)))
  
  for(i in seq_along(index_dt.list))
  {
    tmp_model = strsplit(strsplit(ests_vec[i],"/")[[1]][2],"_")[[1]]
    tmp_model_name = paste0(tmp_model[1],".",tmp_model[4],".",
                            tmp_model[5],".",
                            paste0(sapply(strsplit(tmp_model[6],"[.]"),function(x)paste0(substr(x,1,1),substr(x,nchar(x),nchar(x)))),collapse=""))
 	index_dt.list[[i]] = fread(paste0(proj_dir,"VAST/model_runs/",ests_vec[i])) %>%
                     .[,Model:=paste0(ifelse(i<10,"0",""),i," ",tmp_model_name)] %>%
                     .[,Category:=deep7_name_vec[match(Category,deep7_code_vec)]] %>%
                     .[is.na(Category),Category:="Total"] %>%
                     .[,Category:=factor(Category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]
 
    rm(list=c("tmp_model","tmp_model_name"))
  }

  index_dt = rbind(design_dt,rbindlist(index_dt.list)) %>%
            .[Model=="design",Model:="00 design"] %>%
			.[,model_number:=as.numeric(sapply(Model,function(x)strsplit(x,"\\s+")[[1]][1]))]

 		p = copy(index_dt) %>%
			ggplot() +
			ylim(0,NA) +
			ylab("Predicted biomass (millions lbs)") +
			xlab("Year") +
			facet_wrap(~Category,scales="free") +
			geom_hline(yintercept=0) +
			geom_ribbon(aes(x=Time,ymin=l95,ymax=u95,group=Model,fill=Model),alpha=0.25) +
			geom_path(aes(x=Time,y=Estimate,group=Model,color=Model),linewidth=1.5) +
			viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			viridis::scale_fill_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			theme_few(base_size=20)
        ggsave(filename=paste0("all_2022_index.png"), plot = p, device = "png", path = dir_plot,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)                        

    	p = copy(index_dt) %>%
			.[,Estimate:=Estimate/mean(Estimate),by=.(Model,Category)] %>%
			.[,l95:=exp(log(Estimate)-2*sqrt(log(CV^2+1)))] %>%
			.[,u95:=exp(log(Estimate)+2*sqrt(log(CV^2+1)))] %>%
            ggplot() +
			ylim(0,NA) +
			ylab("Relative biomass") +
			xlab("Year") +
			facet_wrap(~Category,scales="free") +
			geom_hline(yintercept=0) +
			geom_hline(yintercept=1,linetype="dashed") +
			geom_ribbon(aes(x=Time,ymin=l95,ymax=u95,group=Model,fill=Model),alpha=0.25) +
			geom_path(aes(x=Time,y=Estimate,group=Model,color=Model),linewidth=1.5) +
			viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			viridis::scale_fill_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			theme_few(base_size=20)
        ggsave(filename=paste0("all_2022_index_relative.png"), plot = p, device = "png", path = dir_plot,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)   
   
   
    	p = copy(index_dt) %>%
            ggplot() +
			ylim(0,NA) +
			ylab("Predicted biomass (millions lbs)") +
			xlab("Year") +
			facet_wrap(~Category,scales="free") +
			geom_hline(yintercept=0) +
			geom_path(aes(x=Time,y=Estimate,group=Model,color=Model),linewidth=1.5) +
			viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			viridis::scale_fill_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			theme_few(base_size=20)
        ggsave(filename=paste0("all_2022_index_noRibbon.png"), plot = p, device = "png", path = dir_plot,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)                        

    	p = copy(index_dt) %>%
			# .[model_number %in% c(0,11,15,16)] %>%	
			.[,Estimate:=Estimate/mean(Estimate),by=.(Model,Category)] %>%
			.[,l95:=exp(log(Estimate)-2*sqrt(log(CV^2+1)))] %>%
			.[,u95:=exp(log(Estimate)+2*sqrt(log(CV^2+1)))] %>%
            ggplot() +
			ylim(0,NA) +
			ylab("Relative biomass") +
			xlab("Year") +
			facet_wrap(~Category,scales="free") +
			geom_hline(yintercept=0) +
			geom_hline(yintercept=1,linetype="dashed") +
			geom_path(aes(x=Time,y=Estimate,group=Model,color=Model),linewidth=1.5) +
			viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			viridis::scale_fill_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			theme_few(base_size=20)
        ggsave(filename=paste0("all_2022_index_relative_noRibbon.png"), plot = p, device = "png", path = dir_plot,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)  
