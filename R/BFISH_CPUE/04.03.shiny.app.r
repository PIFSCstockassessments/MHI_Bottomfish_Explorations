

# Nicholas Ducharme-Barth
# 2023/06/04
# plot indices shiny app
# 1) index
# 2) effects (abundance & catchability)
# 3) influ
# 4) rmse
# 5) spatial predictions
# 6) residuals

# Copyright (c) 2023 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#_____________________________________________________________________________________________________________________________
# load packages
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(data.table)
library(markdown)
library(magrittr)
library(ggplot2)
library(viridis)
library(ggthemes)
# library(DT)


#_____________________________________________________________________________________________________________________________
# set directory
	proj_dir = "D:/HOME/SAP/2024_Deep7/"
  code_dir = "D:/HOME/SAP/Code/MHI_Bottomfish_2023/"

#_____________________________________________________________________________________________________________________________
# set defaults
  shiny.title = "MHI Deep7 BFISH model based indices"
  copyright.year = format(Sys.time(),"%Y")
#_____________________________________________________________________________________________________________________________
# update data
  date_modified = file.mtime(paste0(proj_dir,"VAST/model_runs/comparison_plots/index_dt.csv"))
  date_modified = as.POSIXct(strsplit(as.character(date_modified),"\\s+")[[1]][1])
  today = Sys.time()
  today = as.POSIXct(strsplit(as.character(today),"\\s+")[[1]][1])

  if(date_modified!=today)
  {
    source(paste0(code_dir,"R/BFISH_CPUE/04.04.gather.indices.shiny.r"))
    rm(list=setdiff(ls(),c("code_dir","proj_dir")))
  }

#_____________________________________________________________________________________________________________________________
# load data
	index_dt = fread(file=paste0(proj_dir,"VAST/model_runs/comparison_plots/index_dt.csv")) %>%
            .[,category:=factor(category,levels=c("Total","'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]
	index_summary_dt = fread(file=paste0(proj_dir,"VAST/model_runs/comparison_plots/index_summary_dt.csv")) 
	abundance_effect_dt = fread(file=paste0(proj_dir,"VAST/model_runs/comparison_plots/abundance_effect_dt.csv")) %>%
            .[,species_cd:=factor(species_cd,levels=c("'Opakapaka (PRFI)","Ehu (ETCA)","Onaga (ETCO)","Kalekale (PRSI)","Gindai (PRZO)","Hapu'upu'u (HYQU)","Lehi (APRU)"))]
	plot_empirical_dt = fread(file=paste0(proj_dir,"VAST/model_runs/comparison_plots/plot_empirical_ab_dt.csv"))
	plot_empirical_discrete_dt = fread(file=paste0(proj_dir,"VAST/model_runs/comparison_plots/plot_empirical_discrete_ab_dt.csv"))


#_____________________________________________________________________________________________________________________________
# app options
  start_collapsed = FALSE

#_____________________________________________________________________________________________________________________________
# The app

css <- HTML(
    "#summarytable > .dataTables_wrapper.no-footer > .dataTables_scroll > .dataTables_scrollBody {
        transform:rotateX(180deg);
    }
    #summarytable > .dataTables_wrapper.no-footer > .dataTables_scroll > .dataTables_scrollBody table{
        transform:rotateX(180deg);
    }"
)

ui = dashboardPage(
  header = dashboardHeader(title=shiny.title),
  sidebar = dashboardSidebar(
    br(),
    br(),
    sidebarMenu(id="sidebarmenu",
      menuItem("Introduction", tabName="introduction"),
      menuItem("Summary table", tabName="table"),
      menuItem("Index plots", tabName="index_plots"),
      menuItem("Effect plots", tabName="effect_plots")
    ),

    # Only show these on the plotting tabs - not Introduction and Summary table tabs
    conditionalPanel(condition="input.sidebarmenu == 'index_plots'",
      # category
      awesomeCheckboxGroup(
      inputId = "category",
      label = "Species", 
        choices = unique(as.character(index_dt$category)),
      selected = unique(as.character(index_dt$category))[1],
      inline = FALSE, 
        status = "danger"),
      # scaled
      switchInput(
      inputId = "scale",  
      label = "Rescale index",
      value=TRUE,
      onLabel = "TRUE",
      offLabel = "FALSE",
      onStatus = "success", 
      offStatus = "danger"),
      # scaled
      switchInput(
      inputId = "se",  
      label = "Show uncertainty",
      value=TRUE,
      onLabel = "TRUE",
      offLabel = "FALSE",
      onStatus = "success", 
      offStatus = "danger"),
      # design gray scale
      sliderTextInput(
      inputId = "index_grayscale",
      label = "Choose\ndesign-based index\ngray-scale:", 
        choices = seq(from=0,to=100,by=10),
        selected = c(30),
        force_edges = TRUE,
      grid = TRUE)
    ),
    conditionalPanel(condition="input.sidebarmenu == 'effect_plots'",
      # species
      awesomeCheckboxGroup(
      inputId = "effect_species",
      label = "Species", 
        choices = levels(abundance_effect_dt[!is.na(species_cd)]$species_cd),
      selected = levels(abundance_effect_dt[!is.na(species_cd)]$species_cd),
      inline = FALSE, 
        status = "danger"),
      awesomeRadio(
        inputId = "effect_trans",
        label = "Effect transformation", 
          choices = c("None", "exp()", "Encounter probability"),
        selected = "Encounter probability",
        inline = FALSE, 
          status = "success"),
      # scaled
      switchInput(
      inputId = "effect_emp_scale",  
      label = "Scale empirical\nby gear",
      value=FALSE,
      onLabel = "TRUE",
      offLabel = "FALSE",
      onStatus = "success", 
      offStatus = "danger")
    ),

    br(),
    br(),
    tags$footer(
      div(style="text-align:center",
        tags$p("version 0.0.1"),
        tags$p(paste("Copyright", copyright.year, "NOAA Fisheries, PIFSC SAP"))
      )
    )
  ), # End of sidebar

  body = dashboardBody(
    tags$head(tags$style(HTML('.wrapper {height: auto !important; position:relative; overflow-x:hidden; overflow-y:hidden}') )),
    tags$head(tags$style(css)),
    # Start of main tab stuff
    tabItems(
      # **** Introduction ****
      tabItem(tabName="introduction", h2("Introduction"),
        fluidRow(column(12, includeMarkdown(paste0(code_dir,"R/BFISH_CPUE/shiny/introduction_index.md"))))
      ), # End of introduction tab

      # **** Summary table ****
      tabItem(tabName="table", h2("Summary table"),
        fluidRow(box(title="Model metrics", collapsed=start_collapsed, solidHeader=TRUE, collapsible=TRUE, status="primary", width=12,
         DT::dataTableOutput("summarytable")))
      ), # End of table tab

      # **** Index plots ****
      tabItem(tabName="index_plots", h2("Index plots"),
        fluidRow(
          box(title="Standardized indices", solidHeader=TRUE, collapsible=TRUE, collapsed=start_collapsed, status="primary", width=12,
            plotOutput("index_plots", height="auto"))
        )
      ), # End of index_plots tab

      # **** Effect plots ****
      tabItem(tabName="effect_plots", h2("Effect plots"),
        fluidRow(
          box(title="Estimated abundance effects: continuous covariates", solidHeader=TRUE, collapsible=TRUE, collapsed=start_collapsed, status="primary", width=12,
            p("Select at least one model with defined continuous abundance covariates."),
            plotOutput("ab_effect_plots_c", height="auto"))
        ),
        fluidRow(
          box(title="Estimated abundance effects: discrete covariates", solidHeader=TRUE, collapsible=TRUE, collapsed=start_collapsed, status="primary", width=12,
            p("Select at least one model with defined discrete abundance covariates."),
            plotOutput("ab_effect_plots_d", height="auto"))
        )
      ) # End of effect_plots tab
    ) # End of tabItems
  ) # End of dashboardBody
)

#---------------------------------------------------------------------------
# Server function

server = function(input, output){
  # pixel height for each panel. i.e row height when plotting by species
  height_per_panel = 350

  # # Update model selection choices
  # observeEvent(c(input$data_year,input$error_structure,input$data_treatment,input$q_config,input$ab_config,input$lehi_filter,input$knot_dist),
  # {
  #   tmp_dt = index_dt %>%
  #             .[as.character(data_year) %in% input$data_year| model_number==0] %>%
  #             .[as.character(error_structure) %in% input$error_structure| model_number==0] %>%
  #             .[as.character(data_treatment) %in% input$data_treatment| model_number==0] %>%
  #             .[as.character(q_config) %in% input$q_config| model_number==0] %>%
  #             .[as.character(ab_config) %in% input$ab_config| model_number==0] %>%
  #             .[as.character(lehi_filter) %in% input$lehi_filter| model_number==0] %>%
  #             .[as.character(knot_dist) %in% input$knot_dist| model_number==0]
    

  #   # Method 2
  #   disabled_choices = index_dt[!(model_name_short %in%tmp_dt$model_name_short)]$model_name_short
  #   updatePickerInput(
  #     session = session, inputId = "model_name_short",
  #     choices = index_dt$model_name_short,
  #     choicesOpt = list(
  #       disabled = disabled_choices,
  #       style = ifelse(disabled_choices,
  #                      yes = "color: rgba(119, 119, 119, 0.5);",
  #                      no = "")
  #     )
  #   )
  # }, ignoreInit = TRUE)

  ref_table_reduced = index_summary_dt %>%
                as.data.frame(.)

  output$summarytable = DT::renderDataTable({
    summary_df = index_summary_dt %>%
                 as.data.frame(.,stringsAsFactors=FALSE)
    summary_DT = DT::datatable(summary_df, filter = 'top',rownames=FALSE,
    options = list(scrollX = TRUE, search = list(regex = TRUE, caseInsensitive = FALSE),pageLength = 25))
    return(summary_DT)
  })
  outputOptions(output, "summarytable", suspendWhenHidden = FALSE)

  filtered_table = reactive({
    req(input$summarytable_rows_selected)
    keep_models = c(ref_table_reduced[input$summarytable_rows_selected, ]$model_number)
    return(as.data.frame(index_dt[model_number%in%keep_models],stringsAsFactors=FALSE))  
  })

  # define plots
  output$index_plots = renderPlot({
    input_models = unique(filtered_table()$model_name_short)
    input_species = input$category
    if(length(input_models) < 1 | length(input_species) < 1){
      return()
    }
    plot_dt = index_dt %>%
              .[(model_name_short %in% unique(c("00 design",input_models)))] %>%
              .[category %in% input_species]
    if(nrow(plot_dt) == 0){
      return()
    }
    if(length(input_species)>2)
    {
      n_col = 3
    } else {
      n_col = length(input_species)
    }
    p = plot_dt %>%
      ggplot() +
			ylim(0,NA) +
			xlab("Year") +
			facet_wrap(~category,scales="free",ncol=n_col)
    if(input$scale)
    {
      p = p + ylab("Relative biomass") +
              geom_hline(yintercept=0) + 
              geom_hline(yintercept=1,linetype="dotted")
      if(input$se)
      {
        p = p + geom_ribbon(aes(x=time,ymin=l95_sc,ymax=u95_sc,group=model_name_short,fill=model_name_short),alpha=0.25)
      }
			p = p + geom_path(aes(x=time,y=estimate_sc,group=model_name_short,color=model_name_short),linewidth=1.5)
    } else {
      p = p + ylab("Predicted biomass (millions lbs)") +
              geom_hline(yintercept=0)
      if(input$se)
      {
        p = p + geom_ribbon(aes(x=time,ymin=l95,ymax=u95,group=model_name_short,fill=model_name_short),alpha=0.25)
      }
			p = p + geom_path(aes(x=time,y=estimate,group=model_name_short,color=model_name_short),linewidth=1.5)
    }
    	p = p + scale_colour_manual("Model",values=c(paste0("gray",input$index_grayscale),viridis::viridis_pal(alpha = 1, begin = 0.1, end = 0.8, direction = 1, option = "H")(uniqueN(plot_dt$model_name_short)-1))) +
			scale_fill_manual("Model",values=c(paste0("gray",input$index_grayscale),viridis::viridis_pal(alpha = 1, begin = 0.1, end = 0.8, direction = 1, option = "H")(uniqueN(plot_dt$model_name_short)-1))) +
			theme_few(base_size=20)
			# p = p + viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			# viridis::scale_fill_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			# theme_few(base_size=20)
    return(p)
  },
  height=function(){
    if(length(input$category)>2)
    {
      n_col = 3
    } else {
      n_col = length(input$category)
    }
    return(max(height_per_panel*1.5, (height_per_panel * ceiling(length(input$category) / n_col))))
  })
  # continuous abundance effects
  output$ab_effect_plots_c = renderPlot({
    input_models = unique(filtered_table()$model_name_short)
    input_species = input$effect_species
    input_trans = input$effect_trans
    input_scale = input$effect_emp_scale
    if(length(input_models) < 1 | length(input_species) < 1){
      return()
    }
    plot_dt = abundance_effect_dt %>%
              .[(model_name_short %in% input_models)] %>%
              .[species_cd %in% input_species] %>%
              setnames(.,c("variable_c","value"),c("variable","value_plot")) %>%
              .[,.(model_number,model_name_short,component,species_cd,category,variable,value_plot)] %>%
              na.omit(.) 
    if(input_trans=="exp()")
    {
      plot_dt = plot_dt %>%
                .[,value_link:=exp(value_plot)] %>%
                .[,value_plot:=NULL] %>%
							  .[,value_link_rescale:=value_link/max(value_link),by=.(model_number,species_cd,category)] %>%
                setnames(.,c("value_link_rescale"),c("value_plot")) %>%
                .[,.(model_number,model_name_short,component,species_cd,category,variable,value_plot)]
    } else if(input_trans=="Encounter probability"){
      plot_dt = merge(plot_dt,index_summary_dt[,.(model_number,error_structure)],by="model_number") %>%
                .[,error_1st:=substr(error_structure,1,2)] %>%
                .[,value_link:=as.numeric(NA)] %>%
                .[error_1st=="pl",value_link:=1-exp(-exp(value_plot))] %>%
                .[error_1st=="lg",value_link:=boot::inv.logit(value_plot)] %>%
                .[,value_plot:=NULL] %>%
							  .[,value_link_rescale:=value_link/max(value_link),by=.(model_number,species_cd,category)] %>%
                setnames(.,c("value_link_rescale"),c("value_plot")) %>%
                .[,.(model_number,model_name_short,component,species_cd,category,variable,value_plot)] %>%
                na.omit(.)
    } else {
      plot_dt = plot_dt %>%
                .[,value_link:=value_plot] %>%
                .[,value_plot:=NULL] %>%
							  .[,value_link_rescale:=scales::rescale(value_link),by=.(model_number,species_cd,category)] %>%
                setnames(.,c("value_link_rescale"),c("value_plot")) %>%
                .[,.(model_number,model_name_short,component,species_cd,category,variable,value_plot)]
    }
    if(nrow(plot_dt) == 0){
      return()
    }

    if(input_scale){
      plot_empirical_dt_internal = copy(plot_empirical_dt) %>%
                          setnames(.,"empirical_sc_a","empirical_sc")
    } else {
      plot_empirical_dt_internal = copy(plot_empirical_dt) %>%
                          setnames(.,"empirical_sc_b","empirical_sc")  
    }

    plot_dt$model_name_short = factor(plot_dt$model_name_short,levels=unique(plot_dt$model_name_short))
    p = plot_empirical_dt_internal %>%
    .[gear_type=="camera"] %>%
    .[category %in% unique(plot_dt$category)&species_cd %in%input_species] %>%
    ggplot() +
    ylim(0,NA) +
    ylab("Relative encounter rate") +
	  xlab("Habitat covariate") +
	  facet_wrap(category~species_cd,scales="free",ncol=length(input_species)) +
	  geom_rect(aes(xmin=breaks-x_width,xmax=breaks,ymin=0,ymax=empirical_sc,fill=gear_type),color="white") +
 	  geom_rect(data=plot_empirical_dt_internal[gear_type=="research_fishing" & category %in% unique(plot_dt$category)&species_cd %in%input_species],aes(xmin=breaks,xmax=breaks+x_width,ymin=0,ymax=empirical_sc,fill=gear_type),color="white") +
    geom_line(data=plot_dt,aes(x=variable,y=value_plot,color=model_name_short,group=model_name_short)) +
    geom_hline(yintercept=0) +
    theme_few(base_size=20) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),axis.ticks.y=element_blank(),axis.text.y=element_blank(),strip.text.x = element_text(size = 10)) +
    scale_fill_manual("Gear type",values=c("gray60","gray80")) +
	  viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE)
   
    return(p)
  },
  height=function(){
      n_col = length(input$effect_species)
    return(max(height_per_panel*2.25, (height_per_panel * ceiling(length(input$effect_species) / n_col))))
  })
  # discrete abundance effects
  output$ab_effect_plots_d = renderPlot({
    input_models = unique(filtered_table()$model_name_short)
    input_species = input$effect_species
    input_trans = input$effect_trans
    input_scale = input$effect_emp_scale
    if(length(input_models) < 1 | length(input_species) < 1){
      return()
    }
    plot_dt = abundance_effect_dt %>%
              .[(model_name_short %in% input_models)] %>%
              .[species_cd %in% input_species] %>%
              setnames(.,c("variable_d","value"),c("variable","value_plot")) %>%
              .[,.(model_number,model_name_short,component,species_cd,category,variable,value_plot)] %>%
              na.omit(.) 
    if(input_trans=="exp()")
    {
      plot_dt = plot_dt %>%
                .[,value_link:=exp(value_plot)] %>%
                .[,value_plot:=NULL] %>%
							  .[,value_link_rescale:=value_link/max(value_link),by=.(model_number,species_cd,category)] %>%
                setnames(.,c("value_link_rescale"),c("value_plot")) %>%
                .[,.(model_number,model_name_short,component,species_cd,category,variable,value_plot)]
    } else if(input_trans=="Encounter probability"){
      plot_dt = merge(plot_dt,index_summary_dt[,.(model_number,error_structure)],by="model_number") %>%
                .[,error_1st:=substr(error_structure,1,2)] %>%
                .[,value_link:=as.numeric(NA)] %>%
                .[error_1st=="pl",value_link:=1-exp(-exp(value_plot))] %>%
                .[error_1st=="lg",value_link:=boot::inv.logit(value_plot)] %>%
                .[,value_plot:=NULL] %>%
							  .[,value_link_rescale:=value_link/max(value_link),by=.(model_number,species_cd,category)] %>%
                setnames(.,c("value_link_rescale"),c("value_plot")) %>%
                .[,.(model_number,model_name_short,component,species_cd,category,variable,value_plot)] %>%
                na.omit(.)
    } else {
      plot_dt = plot_dt %>%
                .[,value_link:=value_plot] %>%
                .[,value_plot:=NULL] %>%
							  .[,value_link_rescale:=scales::rescale(value_link),by=.(model_number,species_cd,category)] %>%
                setnames(.,c("value_link_rescale"),c("value_plot")) %>%
                .[,.(model_number,model_name_short,component,species_cd,category,variable,value_plot)]
    }
    if(nrow(plot_dt) == 0){
      return()
    }

    if(input_scale){
      plot_empirical_dt_internal = copy(plot_empirical_discrete_dt) %>%
                          setnames(.,"empirical_sc_a","empirical_sc")
    } else {
      plot_empirical_dt_internal = copy(plot_empirical_discrete_dt) %>%
                          setnames(.,"empirical_sc_b","empirical_sc")  
    }

    plot_dt$model_name_short = factor(plot_dt$model_name_short,levels=unique(plot_dt$model_name_short))
    p = plot_empirical_dt_internal %>%
    ggplot() +
    ylab("Relative encounter rate") +
	  xlab("Island group") +
	  facet_wrap(category~species_cd,scales="free",ncol=7) +
	  geom_bar(aes(x=breaks,y=empirical_sc,fill=gear_type),stat="identity",color="white",position=position_dodge()) +
    geom_point(data=plot_dt,aes(x=variable,y=value_plot,color=model_name_short,group=model_name_short),position="jitter",size=2) +
    geom_hline(yintercept=0) +
    theme_few(base_size=20) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),axis.ticks.y=element_blank(),axis.text.y=element_blank(),strip.text.x = element_text(size = 10)) +
    scale_fill_manual("Gear type",values=c("gray60","gray80")) +
	  viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE)
    
    return(p)
  },
  height=function(){
      n_col = length(input$effect_species)
    return(max(height_per_panel*2.25, (height_per_panel * ceiling(length(input$effect_species) / n_col))))
  })

} # End of server

shinyApp(ui, server)
