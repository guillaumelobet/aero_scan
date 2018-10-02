# Guillaume Lobet - University of Liege
# Root systems image Generator


options(scipen=999) # Disable scientific notation

#-----------------------------------------------------------------------------------------
#--------------------------- GENERAL OPTIONS ---------------------------------------------
#-----------------------------------------------------------------------------------------

# Main directory, where everthing is stored
dir.base <- "/Users/g.lobet/OneDrive - UCL/03_research/aeroscan/phenotype_2_model/"

# Where is ArchiSimple folder
setwd(dir.base) 

#-----------------------------------------------------------------------------------------
#--------------------------- READ THE EXISTING SIMULATION FILES --------------------------
#-----------------------------------------------------------------------------------------
ls <- list.files("simulation_parameters/")

ls <- gsub(".rparam", "", ls)
ls <- gsub(".pparam", "", ls)
ls <- unique(ls)



#-----------------------------------------------------------------------------------------
#--------------------------- GENERATE THE RSMLS FOR EACH SIMULATION ---------------------------------------------
#-----------------------------------------------------------------------------------------

count <- 0
evol <- 0
for(l in ls){
  
  file.copy(paste0(dir.base, "simulation_parameters/", l, ".pparam"), paste0(dir.base, "www/param.pparam"))
  file.copy(paste0(dir.base, "simulation_parameters/", l, ".rparam"), paste0(dir.base, "www/param.rparam"))
  system("www/a-rsml.out")  
  
  ls2 <- list.files("www/")
  ls2 <- ls2[grepl("rootsystem.rsml", ls2)]
  for(l2 in ls2){
    file.rename(paste0(dir.base, "www/", l2), paste0(dir.base, "simulation_rsmls/",gsub("-param", "", l),"-",l2))
  }
  # Counter for tracking advencement
  count <- count+1
  if(count >= length(ls)/200){
    evol <- evol+1
    count <- 0
    message(paste0(evol*(length(ls)/200), " done."))
  }
}  

