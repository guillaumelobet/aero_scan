
# Guillaume Lobet - University of Liege
# Root systems Random Generator

# The aim of this script is to match the aeroscan data to the simulated ones


library(tidyverse)
library(data.table)
library(stringi)
library(plyr)
library(Hmisc)
library(reshape2)

options(scipen=999) # Disable scientific notation

#-----------------------------------------------------------------------------------------
#--------------------------- GENERAL OPTIONS ---------------------------------------------
#-----------------------------------------------------------------------------------------

# Main directory, where everthing is stored
dir.base <- "/Users/g.lobet/Dropbox/science/projects/0_aeroscan/phenotype_2_model/"

mod <- 1
mod2 <- 50
# Where is ArchiSimple folder
setwd(dir.base) 


simulated <- read_csv("simulated_data.csv")

params <- read_csv("simulated_parameters.csv") %>%
  select(id, r1)

simulated <- merge(simulated, params, by.x="id", by.y="id")


days_to_analyse <- c(7:11, 14: 18)  # These are the days to be used for the matching
layers <- c(1:5)  # These are the days to be used for the matching

# Get the coluln name patterns for the experimental data
base  <- paste0("time_", days_to_analyse)
base1 <- paste0(rep(base, each=length(layers)), "_", rep(layers, length(days_to_analyse)))
names <- c("id", base, base1)

# Load expeirmenta data
experimental <- read_delim("experimental_data/tips_data.txt", "\t", col_names = names)

# error correction
for(i in c(2:ncol(experimental))){
  experimental[,i] <- experimental[,i] * mod
}


# Add the growth rates
gr_data <- read_delim("experimental_data/growthrate.txt", "\t", col_names = c("id", "growth_rate_primary"))
gr_data <- filter(gr_data, id %in% experimental$id)
experimental <- filter(experimental, id %in% gr_data$id)
experimental <- experimental[order(experimental$id),]
gr_data <- gr_data[order(gr_data$id),]
experimental$r1 <- gr_data$growth_rate_primary

# Synch the colnames
ids <- match(colnames(experimental), colnames(simulated))
simulated <- simulated[,ids]





# train <- experimental[,c(1, 5:11, 62)]
# test <- simulated[,c(1, 5:11, 62)]

train <- experimental[,c(1, 2:11)]
test <- simulated[,c(1, 2:11)]

# give more weight to the growth rate

# train <- train %>% mutate(r1 = r1 * mod2)
train[,c(2:11)] <- train[,c(2:11)] * (c(1:10))
# test <- test %>% mutate(r1 = r1 * mod2)
test[,c(2:11)] <- test[,c(2:11)] * (c(1:10))
# 
# train %>%
#   melt(id.vars = "id") %>%
#   mutate(variable = as.numeric(gsub("time_", "", variable))) %>%
#   ggplot(aes(variable, value, group=id)) + 
#     geom_line(alpha=0.1) + 
#     theme_classic()
# 
# test[c(1:1000),] %>%
#   melt(id.vars = "id") %>%
#   mutate(variable = as.numeric(gsub("time_", "", variable))) %>%
#   ggplot(aes(variable, value, group=id)) + 
#   geom_line(alpha=0.1) + 
#   theme_classic()




results <- NULL
x <- rbind(train[,-1], test[,-1])  # make a matrix with the training data and test
d1 <- as.matrix(dist(x))

for(i in c(1:nrow(train))){
  ind <- which(d1[i,-c(1:nrow(train))] == min(d1[i,-c(1:nrow(train))]))[1]
  results <- rbind(results, data.frame(id = train$id[i], match = test$id[ind], dist = min(d1[i,-c(1:nrow(train))])))
}

ggplot(results, aes(dist)) + 
  geom_density(fill="grey") + 
  theme_classic()

# LOOK AT THE FIT
test1 <- merge(results[,c(1,2)], test, by.x="match", by.y="id")%>%
  melt(id.vars=c("id", "match")) %>%
  arrange(id)

train1 <- merge(results[,c(1,2)], train, by.x="id", by.y="id")%>%
  melt(id.vars=c("id", "match")) %>%
  arrange(id)

dat <- cbind(test1, value2 = train1$value)

ggplot(dat, aes(value, value2))+#, colour=id)) + 
  geom_point(colour="#fbd04a") + 
  stat_smooth(method = "lm", se=F, size=1) + 
  theme_classic() + 
  ylab("Experimental data") + 
  xlab("Model data") + 
  geom_abline(intercept = 0, slope=1, lty=2) + 
  theme(axis.line=element_line(colour="white"),
        axis.text.y=element_text(colour="white"),
        axis.ticks=element_line(colour="white"),
        axis.title.x=element_text(colour="white"),
        axis.title.y=element_text(colour="white"),
        text = element_text(size=20),
        axis.text.x = element_text(angle = 45, hjust = 1, colour="white"),
        panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
        plot.background = element_rect(fill = "transparent",colour = NA))  + 
  ggsave("~/Desktop/reg.png", width=6, height=6, bg="transparent")


# LOOK AT THE ACCURACY
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

dat %>%
  merge(regr, by="id") %>%
  filter(r2 > 0.2) %>%
  select(-c(variable, value, value2)) %>%
  melt(id.vars=c("id","match")) %>%
  ggplot(aes(value)) + 
    geom_density(fill="grey") +
    facet_wrap(~variable, scales="free") + 
    theme_classic()


# Test the params values

params <- read_csv("simulated_parameters.csv")

estimated_params <- merge(results, params, by.x="match", by.y="id") %>%
  melt(id.vars=c("id", "match")) %>% 
  select(-match) %>%
  mutate(type = "exp") %>%
  arrange(id)

all_params <- params %>%
  melt(id.vars=c("id")) %>%
  mutate(type = "sim") %>%
  arrange(id)

data <- rbind(estimated_params, all_params)

ggplot(data, aes(value)) + 
  geom_density(aes(colour=type)) +
  facet_wrap(~variable, scales="free") + 
  theme_classic()

train1 <- merge(results, params, by.x="match", by.y="id")


estimated_params <- merge(results, params, by.x="match", by.y="id")
estimated_params <- merge(estimated_params, train, by.x="id", by.y="id") %>%
  mutate(r1.y = r1.y / mod2)

summary(lm(estimated_params$r1.x ~ estimated_params$r1.y))

ggplot(estimated_params, aes(r1.x, r1.y)) + 
  geom_point() +
  geom_abline(intercept = 0, slope = 1, lty=2) + 
  geom_smooth(method="lm") + 
  theme_classic() + 
  ylab("measured growth rate [cm/day]") + 
  xlab('estimated growth rate [cm/day]')


write_csv(train1, "experimental_data/estimated_parameters.csv")


simulated %>%
  mutate(id1 = c(1:nrow(simulated))) %>%
  filter(id1 < 500) %>%
  melt(id.vars=c("id","id1")) %>%
  mutate(variable = as.character(variable)) %>%
  filter(nchar(variable) < 8) %>%
  mutate(variable = as.numeric(gsub("time_", "", variable))) %>%
  ggplot(aes(variable, value, group=id)) + 
    geom_line(alpha=0.1) +
    theme_classic()
