

# Nicholas Ducharme-Barth
# 2023/07/28
# Update BFISH data
# Copyright (c) 2023 Nicholas Ducharme-Barth
# You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

#_____________________________________________________________________________________________________________________________
# set working directory
    code_dir = "D:/HOME/SAP/Code/MHI_Bottomfish_2023/"
	proj.dir = proj_dir = "D:/HOME/SAP/2024_Deep7/"

#_____________________________________________________________________________________________________________________________
# define data_flag
    data_flag = "2022_" # includes data through 2022

#_____________________________________________________________________________________________________________________________
# update camera data
    for(data_treatment in 1:5)
    {
        for(all_lengths in c(TRUE,FALSE))
        {
            ls_original = ls()
            if(all_lengths)
            {
                source(paste0(code_dir,"R/BFISH_CPUE/01.01.0",data_treatment,".format.bfish_camera_data.all_lengths.r"))
            } else {
                source(paste0(code_dir,"R/BFISH_CPUE/01.01.0",data_treatment,".format.bfish_camera_data.r"))
            }
            rm(list=setdiff(ls(),ls_original))
        }
    }

#_____________________________________________________________________________________________________________________________
# update research fishing
        for(all_lengths in c(TRUE,FALSE))
        {
            ls_original = ls()
            if(all_lengths)
            {
                source(paste0(code_dir,"R/BFISH_CPUE/01.02.format.bfish_research_fishing_data.all_lengths.r"))
            } else {
                source(paste0(code_dir,"R/BFISH_CPUE/01.02.format.bfish_research_fishing_data.r"))
            }
            rm(list=setdiff(ls(),ls_original))
        }

#_____________________________________________________________________________________________________________________________
# combine
        for(all_lengths in c(TRUE,FALSE))
        {
            ls_original = ls()
            if(all_lengths)
            {
                source(paste0(code_dir,"R/BFISH_CPUE/02.00.00.combine.bfish_research_fishing_and_camera_data.all_lengths.r"))
            } else {
                source(paste0(code_dir,"R/BFISH_CPUE/02.00.00.combine.bfish_research_fishing_and_camera_data.r"))
            }
            rm(list=setdiff(ls(),ls_original))
        }
