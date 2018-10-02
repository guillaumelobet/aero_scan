# Guillaume Lobet - University of Liege
# Root systems Random Generator

# The aim of this script is to run CRootBOX in a batch mode in order to create
# any number of root systems, with any combinaison of parameters.


library(tidyverse)
library(data.table)
library(stringi)
library(plyr)
library(Hmisc)


options(scipen=999) # Disable scientific notation

#-----------------------------------------------------------------------------------------
#--------------------------- GENERAL OPTIONS ---------------------------------------------
#-----------------------------------------------------------------------------------------

# Main directory, where everthing is stored
dir.base <- "/Users/g.lobet/OneDrive - UCL/03_research/aeroscan/phenotype_2_model/"

# Where is ArchiSimple folder
setwd(dir.base) 

# Load custom functions
source("io_function.R")

# User defined parameters
nsimulations <- 10000        # Total number of simulation to run
deviation <- 3            # Deviation to induce from the base dataset
days_to_analyse <- c(7:11, 14: 18)  # These are the days to be used for the matching

delete_old        <- T            # Delete old simulation data
generate_roots    <- T            # Generate root systems using CRootBox
compile_all_data  <- T            # Compile all the data in a single dataframe
plot_data         <- F            # Plot the data at the end of the simulations
match_data        <- F            # Test the root systme matching algorithms


#-----------------------------------------------------------------------------------------
#--------------------------- READ THE INITIAL PARAMETER FILES ---------------------------------------------
#-----------------------------------------------------------------------------------------
# This is done to always start from the some parameters when we induce variations

dataset_init <- read_rparam('original/param.rparam')     
plant_init <- read_pparam('original/param.pparam')       

# Delete the old simulation data
if(delete_old){
  ls <- list.files("simulation_data")
  for(l in ls){
    unlink(paste0("simulation_data/",l))
  }
  
  ls <- list.files("simulation_parameters/")
  for(l in ls){
    unlink(paste0("simulation_parameters/",l))
  }
  
  ls <- list.files("simulation_rsmls/")
  for(l in ls){
    unlink(paste0("simulation_rsmls/",l))
  }
}

#-----------------------------------------------------------------------------------------
#--------------------------- WRITE NEW PARAMETER FILES ---------------------------------------------
#-----------------------------------------------------------------------------------------

if(generate_roots){
  for(k in c(1:nsimulations)){
    
    unid <- stri_rand_strings(1, 10)  # Get a unique ID for the simulation
    
    plant <- plant_init
    dataset <- dataset_init
    
    print("-----------------------")
    print(paste0("Simulation ",k," started"))
    
    #--------------------------- CHANGE PARAMETERS IN THE DATASETS ---------------------------------------------
  
    # # LMAX for TYPE 1
    # tochange <- dataset %>%
    #   filter(type==1 & param == "lmax") %>%
    #   select(val1) %>%
    #   as.numeric()
    # dataset$val1[dataset$type==1 & dataset$param == "lmax"]  <- round(runif(1, tochange/max(deviation, 1e-9), tochange*deviation), 2)
    # 
    # # LMAX for TYPE 2
    # tochange <- dataset %>%
    #   filter(type==2 & param == "lmax") %>%
    #   select(val1) %>%
    #   as.numeric()
    # dataset$val1[dataset$type==2 & dataset$param == "lmax"]  <- round(runif(1, tochange/max(deviation, 1e-9), tochange*deviation), 2)
    # 
    # LN for TYPE 1
    # tochange <- dataset %>%
    #   filter(type==1 & param == "ln") %>%
    #   select(val1) %>%
    #   as.numeric()
    dataset$val1[dataset$type==1 & dataset$param == "ln"]  <- round(runif(1, 0.1, 2.5), 2)
    
    # LA for TYPE 1
    # tochange <- dataset %>%
    #   filter(type==1 & param == "la") %>%
    #   select(val1) %>%
    #   as.numeric()
    dataset$val1[dataset$type==1 & dataset$param == "la"]  <- round(runif(1, 2, 20), 2)
    
    
    # Growth rate for TYPE 1
    # tochange <- dataset %>%
    #   filter(type==1 & param == "r") %>%
    #   select(val1) %>%
    #   as.numeric()
    dataset$val1[dataset$type==1 & dataset$param == "r"]  <- round(runif(1, 0.2, 3.5), 2)
    
    # Growth rate for TYPE 2
    # tochange <- dataset %>%
    #   filter(type==2 & param == "r") %>%
    #   select(val1) %>%
    #   as.numeric()
    dataset$val1[dataset$type==2 & dataset$param == "r"]  <- round(runif(1, 0.2, 5), 2)
    
    # Number of seminal roots
    plant$val1[plant$param == "maxB"]  <- sample(c(1:10), 1)
    
    # Wrtie the data in the parameter file that will be sued nby CRootBox
    write_rparam(dataset, c("www/param.rparam", paste0("simulation_parameters/",unid,"-param.rparam")))
    write_pparam(plant, c("www/param.pparam", paste0("simulation_parameters/",unid,"-param.pparam")))
  
  
    #--------------------------- RUN CODE AND STORE DATA ---------------------------------------------
  
    system("www/a.out")  
    
    # Get the apex data
    ls <- list.files("www")
    ls <- ls[grepl("rootsystem.txt", ls)]
  
    # compile the simulation data
    rs_temp <- NULL
    all_temp <- NULL
    for(l in ls){
      ti <- as.numeric(strsplit(l, "-")[[1]][1])
      if(ti %in% days_to_analyse){
        temp <- fread(paste0("www/",l), header=T) %>%
          arrange(desc(branchID)) %>%
          mutate(node_diff = c(-1, diff(branchID))) %>%
          filter(node_diff < 0) %>%
          select(c(x2,y2,z2)) %>%
          mutate(x=x2, y=y2, z=z2)
    
        temp$time <- ti
        temp$plant <- k
        
        rs_temp <- rbind(rs_temp, temp)
      }
    }
    
    write.csv(rs_temp, file = paste0("simulation_data/plant-",unid,".csv"))
    
    ls2 <- list.files("www/")
    ls2 <- ls2[grepl("rootsystem.rsml", ls2)]
    for(l2 in ls2){
      file.rename(paste0(dir.base, "www/", l2), paste0(dir.base, "simulation_rsmls/",unid,"-",l2))
    }
  }
}



#-----------------------------------------------------------------------------------------
#--------------------------- COMPILE THE DATA FROM ALL SIMULATIONS ---------------------------------------------
#-----------------------------------------------------------------------------------------


if(compile_all_data){

  # COMPILE ALL THE DATA
  print("-------------")
  print("COMPILING ALL DATA")
  ls <- list.files("simulation_data")
  rootsystem <- NULL
  count <- 0
  evol <- 0
  for(l in ls){
    
    # Load each data and compute agregated metrics
    temp <- fread(paste0("simulation_data/",l), header=T) %>%
      rename(c("V1"="id")) %>%
      mutate("id" = gsub(".csv", "", gsub("plant-", "", l))) %>%
      select(id, x, z, time) 
    
    # Divide the root system by slides (along z axis)
    temp$slide <- 1
    temp$slide[temp$z < -2 & temp$z >= -6.5] <- 2
    temp$slide[temp$z < -8.5 & temp$z >= -14.5] <- 3
    temp$slide[temp$z < -14.5 & temp$z >= -20.5] <- 4
    temp$slide[temp$z < -20.5] <- 5
    
    # Get the total number of tips for each time
    temp1 <- ddply(temp, .(id, time), summarise, tips=length(x)) %>%
      mutate(time = paste0("time_", time)) %>%
      dcast(id ~ time, value.var = "tips")
    
    # Get the total number of tips for each time and each slice
    temp <- ddply(temp, .(id, time, slide), summarise, tips=length(x)) %>%
      mutate(time = paste0("time_", time,"_",slide)) %>%
      dcast(id ~ time, value.var = "tips") %>%
      cbind(temp1[-1])
    
    # Check if all the colums existis (for each time * slice combinaison)
    cols <- colnames(temp)
    for(i in c(1:5)){
      for(j in days_to_analyse){
        if(!paste0("time_",j,"_",i) %in% cols) temp[[paste0("time_",j,"_",i)]] <- 0
      }
    }
    
    # Combine data with previous ones
    rootsystem <- rbind(rootsystem, temp)
    
    # Counter for tracking advencement
    count <- count+1
    if(count >= length(ls)/200){
      evol <- evol+1
      count <- 0
      message(paste0(evol*(length(ls)/200), " done. Data size = ",nrow(rootsystem)))
    }
  }
  write_csv(rootsystem, "simulated_data.csv")
  
  # COMPILE ALL THE PARAMETERS
  print("-------------")
  print("------ COMPILING ALL  PARAMETERS")
  ls <- list.files("simulation_parameters/")
  ls <- gsub("-param.rparam", "", ls)
  ls <- gsub("-param.pparam", "", ls)
  ls <- unique(ls)
  params <- NULL
  count <- 0
  evol <- 0
  for(l in ls){
    da <- read_rparam(paste0('simulation_parameters/',l,'-param.rparam'))
    pl <- read_pparam(paste0('simulation_parameters/',l,'-param.pparam'))  
    params <- rbind(params, data.frame(
      id = l,
      # lmax1 = da$val1[da$type==1 & da$param == "lmax"],
      # lmax2 = da$val1[da$type==2 & da$param == "lmax"],
      ln1 = da$val1[da$type==1 & da$param == "ln"],
      la1 = da$val1[da$type==1 & da$param == "la"],
      r1 = da$val1[da$type==1 & da$param == "r"],
      r2 = da$val1[da$type==2 & da$param == "r"],
      maxB = pl$val1[pl$param == "maxB"]
    )) 
    count <- count+1
    if(count >= length(ls)/200){
      evol <- evol+1
      count <- 0
      message(paste0(evol*(length(ls)/200), " done. Data size = ",nrow(params)))
    }
  }
  
  write_csv(params, "simulated_parameters.csv")
}


#-----------------------------------------------------------------------------------------
#--------------------------- PLOT THE DATA ---------------------------------------------
#-----------------------------------------------------------------------------------------


if(plot_data){
  ddply(rootsystem, .(plant, time), summarise, count = length(x)) %>%
    filter(time==max(time)) %>%
    arrange(count) %>%
    mutate(id = c(1:length(count))) %>%
    ggplot(aes(id, count)) +
      geom_col(alpha=0.5, width=0.1) +
      geom_point(size=1) +
      theme_classic()
  
  
  time %>%
    mutate(sumtime = c(cumsum(time))) %>%
    ggplot(aes(plant, time)) + 
      geom_point() + 
      geom_smooth() + 
      theme_classic()
  
  ddply(rootsystem, .(plant, time), summarise, count = length(x)) %>%
    ggplot(aes(time, count, group=plant)) +
      geom_line(alpha=0.2) +
      #geom_point() +
      theme_classic()
}



#-----------------------------------------------------------------------------------------
#--------------------------- MATCH THE DATA ---------------------------------------------
#-----------------------------------------------------------------------------------------

if(match_data){
    
  rootsystem <- fread("simulated_data.csv", header=T)
  
  train_id <- sample(c(1:nrow(rootsystem)), 500)
  train <- rootsystem[train_id,]
  test <- rootsystem[-train_id,]
  
  
  results <- NULL
  x <- rbind(train[,-1], test[,-1])  # make a matrix with the training data and test
  d1 <- as.matrix(dist(x))
  
  for(i in c(1:nrow(train))){
    ind <- which(d1[i,-c(1:nrow(train))] == min(d1[i,-c(1:nrow(train))]))[1]
    results <- rbind(results, data.frame(id = train$id[i], match = test$id[ind]))
  }

  # Test the results values
  test1 <- merge(results, test, by.x="match", by.y="id")%>%
    melt(id.vars=c("id", "match")) %>%
    arrange(id)
  
  train1 <- merge(results, train, by.x="id", by.y="id")%>%
    melt(id.vars=c("id", "match")) %>%
    arrange(id)
  
  dat <- cbind(test1, value2 = train1$value)
  
  regr <- NULL
  for(i in unique(dat$id)){
    temp <- filter(dat, id == i)
    r2 <- round(summary(lm(temp$value ~ temp$value2))$r.squared, 4)
    pearson <- round(rcorr(temp$value, temp$value2, type = "pearson")[[1]][1,2], 4)
    rmse <- sum(sqrt(((temp$value - temp$value2)/temp$value2)^2))
    regr <- rbind(regr, data.frame(id = i, 
                                   r2 = r2,
                                   pearson = pearson,
                                   rmse = rmse))
  }
  
  dat <- merge(dat, regr, by="id") %>%
    filter(r2 > 0.2)
  
  ggplot(dat, aes(value, value2, colour=variable)) + 
    geom_point() + 
    theme_classic() + 
    geom_abline(intercept = 0, slope=1, lty=2) + 
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1)) 
  
  ggplot(dat, aes(r2)) + 
    geom_density(fill="grey") + 
    theme_classic()
    
  
  # Test the params values
  
  params <- fread("simulated_parameters.csv", header=T)
  
  test1 <- merge(results, params, by.x="match", by.y="id")%>%
    melt(id.vars=c("id", "match")) %>%
    arrange(id)
  
  train1 <- merge(results, params, by.x="id", by.y="id")%>%
    melt(id.vars=c("id", "match")) %>%
    arrange(id)
  
  dat <- cbind(test1, value2 = train1$value)
  
  regr <- NULL
  for(i in unique(dat$variable)){
    temp <- filter(dat, variable == i)
    r2 <- round(summary(lm(temp$value ~ temp$value2))$r.squared, 4)
    pearson <- round(rcorr(temp$value, temp$value2, type = "pearson")[[1]][1,2], 4)
    regr <- rbind(regr, data.frame(variable = i, 
                                   r2 = paste0(i,"  ||  r-squared =  ",r2),
                                   pearson = pearson))
  }
  
  merge(dat, regr, by="variable")%>%
    ggplot(aes(value, value2)) + 
      geom_point() + 
      theme_classic() + 
      facet_wrap(~r2, scales="free") +
      stat_smooth(method='lm', se=F) +
      geom_abline(intercept = 0, slope=1, lty=2)  
}





