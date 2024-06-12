library(decisionSupport)
library(dplyr)
library(ggplot2)
library(vctrs)

folder <- "C:/Users/mjimenez/Documents/ABCDR/a Research/CDR_dynamics_AF/SRC_6years"
input_file <- paste0(folder, "/", "CDR_ACSRC_full_cover_yield_values.csv")
input_data <- read.csv(input_file)
runs <- 1e4 # number of runs of the Monte Carlo simulation
milestones <- c(60,90,120,150,180)
milestones_change <- c(59,89,119,149,179)
scenario_milestones <- c(2030, 2035, 2040, 2045, 2050, 2060, 2100)
current_year <- 2024
target_time <- 20
reporting_years <- sort(c(scenario_milestones, current_year+60))-current_year
lenght_rotations <- c(6, 6, 6)
annual_upscale_rate <- 280*100

#Read the values of the model variables. This must be contained in a .csv file
make_variables <- function(est,n=1)
{x <- random(rho=est,n=n)
for(i in colnames(x))assign(i, as.numeric(x[1,i]),envir=.GlobalEnv)}

make_variables(estimate_read_csv(input_file))

#Define each variable as vectors with as many values as the number of years covered by the simulation
time <- 1:n_years # a vector containing the number of each of the years over which the simulation will be run
AGB <- rep(0,n_years)# the stock of aboveground biomass
AGBcarbon <- rep(0,n_years) #the carbon stock in the growing perennial aboveground biomass
biomass_harvest <- rep(0,n_years) #the mass of harvested wood
c_harvest <- rep(0,n_years) #the carbon contained in the harvested wood
c_hwp <- rep(0,n_years) #the carbon contained in the harvested wood accounting for the decay rate
BGB <- # the stock of belowground biomass
  BGBcarbon_capture <- rep(0,n_years) #the absolute annual change of carbon stock in the roots of the growing woody component of the system
dead_BGBcarbon <- rep(0,n_years) #the cumulative amount of BGB carbon that has been grown in the system and died over the whole length of the simulation
annual_soil_carbon_change <- rep(0,n_years) # the average annual change of carbon contained in the soil (excluding roots of the woody component) in tones per hectare
soil_C_stock <- rep(0,n_years)
years_rootstock_rotation <- years_rotation*number_of_rotations

#Simulate the growth dynamics of the coppice stand
#Annual change of carbon captured in Aboveground Woody Biomass
AGBgrowth <- rep(c(biomass_first_rotation, 
                   biomass_first_rotation*biomass_second_rotation, 
                   biomass_first_rotation*biomass_third_rotation), 
                 each = unique(lenght_rotations),
                 length.out = years_rootstock_rotation)
AGBgrowth <- vec_rep(AGBgrowth, times = n_years/years_rootstock_rotation)
AGBgrowth <- vv(AGBgrowth,CV_biomass_prod,n_years)
AGBproduction <- AGBgrowth * wooded_area
AGBcarbon_capture <- AGBproduction * carbon_density

#Carbon captured in Aboveground Woody Biomass
first_rotation_year <- c(1, seq(from = unique(lenght_rotations)+1, to = n_years, by = unique(lenght_rotations)))
harvest_year <- seq(from = unique(lenght_rotations), to = n_years, by = unique(lenght_rotations))

for (i in 1:length(first_rotation_year)){
  AGB[first_rotation_year[i]:harvest_year[i]] <- cumsum(AGBproduction[first_rotation_year[i]:harvest_year[i]])
}
AGBcarbon <- AGB * carbon_density

#Carbon in harvested wood biomass
biomass_harvest[harvest_year] <- cumsum(AGB[harvest_year])
for (i in 1:length(harvest_year)){
  biomass_harvest[c(first_rotation_year[i]:(harvest_year[i]-1))] <- biomass_harvest[harvest_year[1:(length(harvest_year)-1)][i]]
}
c_harvest <- biomass_harvest * carbon_density

# Carbon in harvested wood products, accounting for decay
c_hwp[harvest_year] <- AGBcarbon[harvest_year]
for (i in 1:(n_years/years_rotation-1)){
  c_hwp[(years_rotation*i+1):n_years] <- c_hwp[(years_rotation*i+1):n_years] + AGBcarbon[years_rotation*i] * exp(-log(2)/hl_wood_pulp * 1:(n_years-years_rotation*i))
}

#Carbon captured in Belowground Tree Biomass
BGB <- AGB[1:years_rootstock_rotation] * root_to_shoot_ratio
for (i in 1:(length(BGB)-1)){
  BGB[i+1] <- ifelse(BGB[i+1] < BGB[i], BGB[i], BGB[i+1])
}
BGB <- vec_rep(BGB, times = n_years/years_rootstock_rotation)

BGBcarbon <- AGBcarbon[1:years_rootstock_rotation] * root_to_shoot_ratio
for (i in 1:(length(BGBcarbon)-1)){
  BGBcarbon[i+1] <- ifelse(BGBcarbon[i+1] < BGBcarbon[i], BGBcarbon[i], BGBcarbon[i+1])
}
BGBcarbon <- vec_rep(BGBcarbon, times = n_years/years_rootstock_rotation)

rootstock_turnover_years <- seq(from = years_rootstock_rotation, to = n_years, by = years_rootstock_rotation)
for (i in 2:length(rootstock_turnover_years)){
  current_turnover_year <- rootstock_turnover_years[i]
  dead_BGBcarbon[(current_turnover_year-years_rootstock_rotation+1):current_turnover_year] <- dead_BGBcarbon[current_turnover_year-years_rootstock_rotation] + BGBcarbon[current_turnover_year-years_rootstock_rotation]
  rm(current_turnover_year)
}
#Annual change of carbon captured in Belowground Tree Biomass
BGBcarbon_capture[1] <- AGBcarbon_capture[1] * root_to_shoot_ratio
BGBcarbon_capture[2:n_years] <- diff(BGBcarbon)
BGBcarbon_capture[BGBcarbon_capture < 0] <- 0

#Carbon sequestered in the soil (excluding tree roots)
soil_C_stock[1] <- initial_soil_C_stock
total_biomass_C <- AGBcarbon + BGBcarbon + dead_BGBcarbon
#max_attainable_AGBc <- subset(input_data, input_data$variable == "biomass_third_rotation")$upper * 0.48 * unique(lenght_rotations)
max_attainable_AGBc <- subset(input_data, input_data$variable == "max_att_agb")$upper * 0.48
max_attainable_BGBc <- max_attainable_AGBc * subset(input_data, input_data$variable == "root_to_shoot_ratio")$upper
relative_biomass_C <- total_biomass_C/(max_attainable_AGBc + max_attainable_BGBc + max(dead_BGBcarbon))

for (i in 2:n_years){
  if (max_soil_C_stock > soil_C_stock[i-1]){
    soil_C_stock[i] <- soil_C_stock[i-1] + min(max_soc_seq_rate, (max_soil_C_stock - soil_C_stock[i-1]) * relative_biomass_C[i] * i/n_years)
  } else {
    soil_C_stock[i] <- soil_C_stock[i-1]
  }
}

annual_soil_carbon_change[2:n_years] <- diff(soil_C_stock)
#soil_carbon_sequestration <- cumsum(annual_soil_carbon_change)
soc_perthousand_increase <- annual_soil_carbon_change / soil_C_stock * 1000