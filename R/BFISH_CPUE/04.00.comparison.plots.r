

# Nicholas Ducharme-Barth
# 12/27/2022
# Comparison plots

# Copyright (c) 2022 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#_____________________________________________________________________________________________________________________________
# load packages
	library(data.table)
	library(magrittr)
	library(ggplot2)
	library(ggthemes)


#_____________________________________________________________________________________________________________________________
# set working directory
	proj.dir = "D:/HOME/SAP/2024_Deep7/"
	dir_plot = paste0(proj.dir,"VAST/model_runs/comparison_plots/")
	dir.create(dir_plot,recursive = TRUE)

#_____________________________________________________________________________________________________________________________
# bring in indices
    # design based bfish index data is in kg so convert to lbs (millions)
	design_dt = fread(paste0(proj.dir,"Data/2021_design_based_estimates.csv")) %>%
				.[,.(Model,Category,Time,Estimate,CV)] %>%
                .[,Estimate:=Estimate*2.20462262185] %>%
				.[,Estimate:=Estimate/1000000] %>%
				.[,l95:=exp(log(Estimate)-2*sqrt(log(CV^2+1)))] %>%
				.[,u95:=exp(log(Estimate)+2*sqrt(log(CV^2+1)))]

#_____________________________________________________________________________________________________________________________
# plot 1) apru single species, poisson-link delta-gamma model
# data filter: 05
# with and without exclusion of the large lehi sample

    # bring in the index_dt from each model fit
    # Note: index_dt is in lbs (millions)

    apru_filter_dt = fread(paste0(proj.dir,"VAST/model_runs/2022-12-27/2021_pldg_apru_05_v_v_TRUE_7.5_FALSE_TRUE_pit_noxval/index_dt.csv")) %>%
                     .[,Model:="model - filter"]
    apru_nofilter_dt = fread(paste0(proj.dir,"VAST/model_runs/2022-12-27/2021_pldg_apru_05_v_v_FALSE_7.5_FALSE_TRUE_pit_noxval/index_dt.csv")) %>%
                     .[,Model:="model - no filter"]

    p1 = rbind(design_dt,apru_filter_dt,apru_nofilter_dt) %>%
            .[Category == "apru"] %>%
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
        ggsave(filename=paste0("p1_apru_filter_comp.png"), plot = p1, device = "png", path = dir_plot,
	  			scale = 1.25, width = 16, height = 9, units = c("in"),
	  			dpi = 300, limitsize = TRUE)


