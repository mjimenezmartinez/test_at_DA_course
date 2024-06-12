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

#Define the model
AF_CDR <- function(x,varnames)
{
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
  AGBgrowth <- rep(c(biomass_first_rotation, biomass_second_rotation, biomass_third_rotation), 
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
  
  c_harvest[harvest_year] <- biomass_harvest * carbon_density
  
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

  return(list(AGwoodyB_carbon = AGBcarbon, 
              BGwoodyB_carbon = BGBcarbon,
              total_woody_biomass_carbon = AGBcarbon + BGBcarbon,
              soil_carbon_storage = cumsum(annual_soil_carbon_change),
              pew_biomass = biomass_harvest,
              bccs = c_harvest,
              c_products = c_hwp,
              total_woody_biomass_carbon_hwp = AGBcarbon + BGBcarbon + c_hwp,
              total_woody_biomass_carbon_beccs = AGBcarbon + BGBcarbon + c_harvest,
              AF_carbon = cumsum(annual_soil_carbon_change) + AGBcarbon + BGBcarbon,
              lca_carbon = cumsum(annual_soil_carbon_change) + AGBcarbon + BGBcarbon + c_hwp,
              bccs_lca_carbon = cumsum(annual_soil_carbon_change) + AGBcarbon + BGBcarbon + c_harvest,
              annual_AGB_C_capture = AGBcarbon_capture,
              annual_BGB_C_capture = BGBcarbon_capture,
              annual_soil_C_capture = annual_soil_carbon_change,
              annual_AF_C_capture = AGBcarbon_capture + BGBcarbon_capture + annual_soil_carbon_change,
              sequestered_carbon_wood = sum(AGBcarbon_capture),
              sequestered_roots = sum(BGBcarbon_capture),
              sequestered_carbon_mineral_soil = sum(annual_soil_carbon_change),
              sequestered_carbon_soil = sum(annual_soil_carbon_change) + sum(BGBcarbon_capture),
              sequestered_carbon_total = sum(annual_soil_carbon_change) + sum(BGBcarbon_capture) + sum(AGBcarbon_capture),
              four_perthousand_initiative = soc_perthousand_increase)
         )
}

#Run again the Monte Carlo analysis of the model
mcSimulation_results <- decisionSupport::mcSimulation(
  estimate = decisionSupport::estimate_read_csv(input_file),
  model_function = AF_CDR,
  numberOfModelRuns = runs,
  functionSyntax = "plainNames")

#Run again the Monte Carlo analysis, this time returning output files
decisionSupport(inputFilePath = input_file,
                outputPath = paste0(folder, "/", "MCResults", sep=""),
                write_table = TRUE,
                welfareFunction = AF_CDR,
                numberOfModelRuns = runs, #run 10,000 times
                functionSyntax = "plainNames",
                plsrVipAnalysis = FALSE)

MCall <- read.table(file = paste0(folder, "/", "MCResults", "/", "mcSimulationResults.csv"),
                    header = TRUE, sep=",")
colnames(MCall)
################################################################
n_years <- subset(input_data, input_data$variable == "n_years")$upper
landscape_scale_rotation_years <- subset(input_data, input_data$variable == "n_years")$upper / subset(input_data, input_data$variable == "number_of_rotations")$upper
# 1. Carbon capture in one agroforestry plot since establishment on treeless agricultural land
AGCwoody <- MCall[, paste0("AGwoodyB_carbon", 1:n_years)]
#Mean annual carbon balance of the AGB carbon pool
MAC_AGCwoody <- AGCwoody[,1:n_years]/(1:n_years)
CAC_AGCwoody <- AGCwoody
CAC_AGCwoody[,2:n_years] <- AGCwoody[,2:n_years] - AGCwoody[,1:(n_years-1)]
write.csv(AGCwoody, paste0(folder, "/", "AGCwoody.csv"))
quant_stats_unscaledAGCwoody <- rbind(
  apply(AGCwoody[, reporting_years],2,min), 
  apply(AGCwoody[, reporting_years], 2, 
        quantile, probs = c(0.05, 0.25, 0.50, 0.75, 0.95), na.rm = TRUE ),
  apply(AGCwoody[, reporting_years],2,max))
write.csv(quant_stats_unscaledAGCwoody, paste0(folder, "/", "quant_stats_unscaledAGCwoody.csv"))

BGCwoody <- MCall[, paste0("BGwoodyB_carbon", 1:n_years)]
MAC_BGCwoody <- BGCwoody[,1:n_years]/(1:n_years)
CAC_BGCwoody <- BGCwoody
CAC_BGCwoody[,2:n_years] <- BGCwoody[,2:n_years] - BGCwoody[,1:(n_years-1)]
write.csv(BGCwoody, paste0(folder, "/", "BGCwoody.csv"))
quant_stats_unscaledBGCwoody <- rbind(
  apply(BGCwoody[, reporting_years],2,min), 
  apply(BGCwoody[, reporting_years], 2, 
        quantile, probs = c(0.05, 0.25, 0.50, 0.75, 0.95), na.rm = TRUE ),
  apply(BGCwoody[, reporting_years],2,max))
write.csv(quant_stats_unscaledBGCwoody, paste0(folder, "/", "quant_stats_unscaledBGCwoody.csv"))

growing_biom_c <- MCall[, paste0("total_woody_biomass_carbon", 1:n_years)]
MAC_growing_biom_c <- growing_biom_c[,1:n_years]/(1:n_years)
CAC_growing_biom_c <- growing_biom_c
CAC_growing_biom_c[,2:n_years] <- growing_biom_c [,2:n_years] - growing_biom_c [,1:(n_years-1)]
write.csv(growing_biom_c, paste0(folder, "/", "growing_biom_c.csv"))
quant_stats_unscaledgrowing_biom_c <- rbind(
  apply(growing_biom_c[, reporting_years],2,min), 
  apply(growing_biom_c[, reporting_years], 2, 
        quantile, probs = c(0.05, 0.25, 0.50, 0.75, 0.95), na.rm = TRUE ),
  apply(growing_biom_c[, reporting_years],2,max))
write.csv(quant_stats_unscaledgrowing_biom_c, paste0(folder, "/", "quant_stats_growing_biom_c.csv"))

SOC <- MCall[, paste0("soil_carbon_storage", 1:n_years)]
MAC_SOC <- SOC[,1:n_years]/(1:n_years)
CAC_SOC <- SOC
CAC_SOC[,2:n_years] <- SOC[,2:n_years] - SOC[,1:(n_years-1)]
write.csv(SOC, paste0(folder, "/", "SOC.csv"))
quant_stats_unscaledSOC <- rbind(
  apply(SOC[, reporting_years],2,min), 
  apply(SOC[, reporting_years], 2, 
        quantile, probs = c(0.05, 0.25, 0.50, 0.75, 0.95), na.rm = TRUE ),
  apply(SOC[, reporting_years],2,max))
write.csv(quant_stats_unscaledSOC, paste0(folder, "/", "quant_stats_unscaledSOC.csv"))

pew_potential <- MCall[, paste0("pew_biomass", 1:n_years)]
MAC_pew_potential <- pew_potential[,1:n_years]/(1:n_years)
CAC_pew_potential <- pew_potential
CAC_pew_potential[,2:n_years] <- pew_potential[,2:n_years] - pew_potential[,1:(n_years-1)]
write.csv(pew_potential, paste0(folder, "/", "pew_potential.csv"))
quant_stats_unscaled_pew_potential <- rbind(
  apply(pew_potential[, reporting_years],2,min), 
  apply(pew_potential[, reporting_years], 2, 
        quantile, probs = c(0.05, 0.25, 0.50, 0.75, 0.95), na.rm = TRUE ),
  apply(pew_potential[, reporting_years],2,max))
write.csv(quant_stats_unscaled_pew_potential, paste0(folder, "/", "quant_stats_unscaled_pew_potential.csv"))

beccs_potential <- MCall[, paste0("bccs", 1:n_years)]
MAC_beccs_potential <- beccs_potential[,1:n_years]/(1:n_years)
CAC_beccs_potential <- beccs_potential
CAC_beccs_potential[,2:n_years] <- beccs_potential[,2:n_years] - beccs_potential[,1:(n_years-1)]
write.csv(beccs_potential, paste0(folder, "/", "beccs_potential.csv"))
quant_stats_unscaled_beccs_potential <- rbind(
  apply(beccs_potential[, reporting_years],2,min), 
  apply(beccs_potential[, reporting_years], 2, 
        quantile, probs = c(0.05, 0.25, 0.50, 0.75, 0.95), na.rm = TRUE ),
  apply(beccs_potential[, reporting_years],2,max))
write.csv(quant_stats_unscaled_beccs_potential, paste0(folder, "/", "quant_stats_unscaled_beccs_potential.csv"))

c_durable_products <- MCall[, paste0("c_products", 1:n_years)]
MAC_c_durable_products <- c_durable_products[,1:n_years]/(1:n_years)
CAC_c_durable_products <- c_durable_products
CAC_c_durable_products[,2:n_years] <- c_durable_products[,2:n_years] - c_durable_products[,1:(n_years-1)]
write.csv(c_durable_products, paste0(folder, "/", "c_durable_products.csv"))
quant_stats_unscaled_c_durable_products <- rbind(
  apply(c_durable_products[, reporting_years],2,min), 
  apply(c_durable_products[, reporting_years], 2, 
        quantile, probs = c(0.05, 0.25, 0.50, 0.75, 0.95), na.rm = TRUE ),
  apply(c_durable_products[, reporting_years],2,max))
write.csv(quant_stats_unscaled_c_durable_products, paste0(folder, "/", "quant_stats_unscaled_c_durable_products.csv"))

c_hwp <- MCall[, paste0("total_woody_biomass_carbon_hwp", 1:n_years)]
MAC_c_hwp <- c_hwp[,1:n_years]/(1:n_years)
CAC_c_hwp <- c_hwp
CAC_c_hwp[,2:n_years] <- c_hwp[,2:n_years] - c_hwp[,1:(n_years-1)]
write.csv(c_hwp, paste0(folder, "/", "c_hwp.csv"))
quant_stats_unscaled_c_hwp <- rbind(
  apply(c_hwp[, reporting_years],2,min), 
  apply(c_hwp[, reporting_years], 2, 
        quantile, probs = c(0.05, 0.25, 0.50, 0.75, 0.95), na.rm = TRUE ),
  apply(c_hwp[, reporting_years],2,max))
write.csv(quant_stats_unscaled_c_hwp, paste0(folder, "/", "quant_stats_unscaled_c_hwp.csv"))

c_beccs <- MCall[, paste0("total_woody_biomass_carbon_beccs", 1:n_years)]
MAC_c_beccs <- c_beccs[,1:n_years]/(1:n_years)
CAC_c_beccs <- c_beccs
CAC_c_beccs[,2:n_years] <- c_beccs[,2:n_years] - c_beccs[,1:(n_years-1)]
write.csv(c_beccs, paste0(folder, "/", "c_beccs.csv"))
quant_stats_unscaled_c_beccs <- rbind(
  apply(c_beccs[, reporting_years],2,min), 
  apply(c_beccs[, reporting_years], 2, 
        quantile, probs = c(0.05, 0.25, 0.50, 0.75, 0.95), na.rm = TRUE ),
  apply(c_beccs[, reporting_years],2,max))
write.csv(quant_stats_unscaled_c_beccs, paste0(folder, "/", "quant_stats_unscaled_c_beccs.csv"))

AFC <- MCall[, paste0("AF_carbon", 1:n_years)]
MAC_AFC <- AFC[,1:n_years]/(1:n_years)
CAC_AFC <- AFC
CAC_AFC[,2:n_years] <- AFC[,2:n_years] - AFC[,1:(n_years-1)]
write.csv(AFC, paste0(folder, "/", "AFC.csv"))
quant_stats_unscaled_AFC <- rbind(
  apply(AFC[, reporting_years],2,min), 
  apply(AFC[, reporting_years], 2, 
        quantile, probs = c(0.05, 0.25, 0.50, 0.75, 0.95), na.rm = TRUE ),
  apply(AFC[, reporting_years],2,max))
write.csv(quant_stats_unscaled_AFC, paste0(folder, "/", "quant_stats_unscaled_AFC.csv"))

AFClca <- MCall[, paste0("lca_carbon", 1:n_years)]
MAC_AFClca <- AFClca[,1:n_years]/(1:n_years)
CAC_AFClca <- AFClca
CAC_AFClca[,2:n_years] <- AFClca[,2:n_years] - AFClca[,1:(n_years-1)]
write.csv(AFClca, paste0(folder, "/", "AFClca.csv"))
quant_stats_unscaled_AFClca <- rbind(
  apply(AFClca[, reporting_years],2,min), 
  apply(AFClca[, reporting_years], 2, 
        quantile, probs = c(0.05, 0.25, 0.50, 0.75, 0.95), na.rm = TRUE ),
  apply(AFClca[, reporting_years],2,max))
write.csv(quant_stats_unscaled_AFClca, paste0(folder, "/", "quant_stats_unscaled_AFClca.csv"))

AFCbeccs <- MCall[, paste0("bccs_lca_carbon", 1:n_years)]
MAC_AFCbeccs <- AFCbeccs[,1:n_years]/(1:n_years)
CAC_AFCbeccs <- AFCbeccs
CAC_AFCbeccs[,2:n_years] <- AFCbeccs[,2:n_years] - AFCbeccs[,1:(n_years-1)]
write.csv(AFCbeccs, paste0(folder, "/", "AFCbeccs.csv"))
quant_stats_unscaled_AFCbeccs <- rbind(
  apply(AFCbeccs[, reporting_years],2,min), 
  apply(AFCbeccs[, reporting_years], 2, 
        quantile, probs = c(0.05, 0.25, 0.50, 0.75, 0.95), na.rm = TRUE ),
  apply(AFCbeccs[, reporting_years],2,max))
write.csv(quant_stats_unscaled_AFCbeccs, paste0(folder, "/", "quant_stats_unscaled_AFCbeccs.csv"))

MAC_mcSimulation_results <- mcSimulation_results
MAC_mcSimulation_results$y[,paste0("AGwoodyB_carbon", 1:n_years)] <- MAC_AGCwoody
MAC_mcSimulation_results$y[,paste0("BGwoodyB_carbon", 1:n_years)] <- MAC_BGCwoody
MAC_mcSimulation_results$y[,paste0("total_woody_biomass_carbon", 1:n_years)] <- MAC_growing_biom_c
MAC_mcSimulation_results$y[,paste0("soil_carbon_storage", 1:n_years)] <- MAC_SOC
MAC_mcSimulation_results$y[,paste0("pew_biomass", 1:n_years)] <- MAC_pew_potential
MAC_mcSimulation_results$y[,paste0("bccs", 1:n_years)] <- MAC_beccs_potential
MAC_mcSimulation_results$y[,paste0("c_products", 1:n_years)] <- MAC_c_hwp
MAC_mcSimulation_results$y[,paste0("total_woody_biomass_carbon_beccs", 1:n_years)] <- MAC_c_beccs
MAC_mcSimulation_results$y[,paste0("total_woody_biomass_carbon_hwp", 1:n_years)] <- MAC_c_hwp
MAC_mcSimulation_results$y[,paste0("AF_carbon", 1:n_years)] <- MAC_AFC
MAC_mcSimulation_results$y[,paste0("lca_carbon", 1:n_years)] <- MAC_AFClca
MAC_mcSimulation_results$y[,paste0("bccs_lca_carbon", 1:n_years)] <- MAC_AFCbeccs
MAC_mcSimulation_results$y[,paste0("AGwoodyB_carbon", 1:n_years)]
tail(MAC_mcSimulation_results$y[,paste0("AGwoodyB_carbon", 1:n_years)])

CAC_mcSimulation_results <- mcSimulation_results
CAC_mcSimulation_results$y[,paste0("AGwoodyB_carbon", 1:n_years)] <- CAC_AGCwoody
CAC_mcSimulation_results$y[,paste0("BGwoodyB_carbon", 1:n_years)] <- CAC_BGCwoody
CAC_mcSimulation_results$y[,paste0("total_woody_biomass_carbon", 1:n_years)] <- CAC_growing_biom_c
CAC_mcSimulation_results$y[,paste0("soil_carbon_storage", 1:n_years)] <- CAC_SOC
CAC_mcSimulation_results$y[,paste0("pew_biomass", 1:n_years)] <- CAC_pew_potential
CAC_mcSimulation_results$y[,paste0("bccs", 1:n_years)] <- CAC_beccs_potential
CAC_mcSimulation_results$y[,paste0("c_products", 1:n_years)] <- CAC_c_hwp
CAC_mcSimulation_results$y[,paste0("total_woody_biomass_carbon_beccs", 1:n_years)] <- CAC_c_beccs
CAC_mcSimulation_results$y[,paste0("total_woody_biomass_carbon_hwp", 1:n_years)] <- CAC_c_hwp
CAC_mcSimulation_results$y[,paste0("AF_carbon", 1:n_years)] <- CAC_AFC
CAC_mcSimulation_results$y[,paste0("lca_carbon", 1:n_years)] <- CAC_AFClca
CAC_mcSimulation_results$y[,paste0("bccs_lca_carbon", 1:n_years)] <- CAC_AFCbeccs
CAC_mcSimulation_results$y[,paste0("AGwoodyB_carbon", 1:n_years)]
tail(CAC_mcSimulation_results$y[,paste0("AGwoodyB_carbon", 1:n_years)])
###
###### 2. SIMULATE THE TEMPORAL DEVELOPMENT OF AN AGROFORESTRY SYSTEM
######i.e.: every year trees are planted in one additional plot (while trees continue to grow in the previously planted plot/s)
#2.1. Create as many replicates of each simulation run as the length (in years) of the rotation,
#also add a "farm" and a "plot_id" column in order to identify each row
AGCwoody <- AGCwoody %>% slice(rep(1:n(), each = landscape_scale_rotation_years))
n_farms <- floor(nrow(AGCwoody)/landscape_scale_rotation_years)
AGCwoody$farm <- rep(paste0("farm", 1:n_farms), each = landscape_scale_rotation_years)
AGCwoody$plot_id  <- rep(1:landscape_scale_rotation_years, n_farms)
upscale_AGCwoody <- AGCwoody

BGCwoody <- BGCwoody %>% slice(rep(1:n(), each = landscape_scale_rotation_years))
n_farms <- floor(nrow(BGCwoody)/landscape_scale_rotation_years)
BGCwoody$farm <- rep(paste0("farm", 1:n_farms), each = landscape_scale_rotation_years)
BGCwoody$plot_id  <- rep(1:landscape_scale_rotation_years, n_farms)
upscale_BGCwoody <- BGCwoody

growing_biom_c <- growing_biom_c %>% slice(rep(1:n(), each = landscape_scale_rotation_years))
n_farms <- floor(nrow(growing_biom_c)/landscape_scale_rotation_years)
growing_biom_c$farm <- rep(paste0("farm", 1:n_farms), each = landscape_scale_rotation_years)
growing_biom_c$plot_id  <- rep(1:landscape_scale_rotation_years, n_farms)
upscale_growing_biom_c <- growing_biom_c

SOC <- SOC %>% slice(rep(1:n(), each = landscape_scale_rotation_years))
n_farms <- floor(nrow(SOC)/landscape_scale_rotation_years)
SOC$farm <- rep(paste0("farm", 1:n_farms), each = landscape_scale_rotation_years)
SOC$plot_id  <- rep(1:landscape_scale_rotation_years, n_farms)
upscale_SOC <- SOC

pew_potential <- pew_potential %>% slice(rep(1:n(), each = landscape_scale_rotation_years))
n_farms <- floor(nrow(pew_potential)/landscape_scale_rotation_years)
pew_potential$farm <- rep(paste0("farm", 1:n_farms), each = landscape_scale_rotation_years)
pew_potential$plot_id  <- rep(1:landscape_scale_rotation_years, n_farms)
upscale_pew_potential <- pew_potential

beccs_potential <- beccs_potential %>% slice(rep(1:n(), each = landscape_scale_rotation_years))
n_farms <- floor(nrow(beccs_potential)/landscape_scale_rotation_years)
beccs_potential$farm <- rep(paste0("farm", 1:n_farms), each = landscape_scale_rotation_years)
beccs_potential$plot_id  <- rep(1:landscape_scale_rotation_years, n_farms)
upscale_beccs_potential <- beccs_potential

c_durable_products <- c_durable_products %>% slice(rep(1:n(), each = landscape_scale_rotation_years))
n_farms <- floor(nrow(c_durable_products)/landscape_scale_rotation_years)
c_durable_products$farm <- rep(paste0("farm", 1:n_farms), each = landscape_scale_rotation_years)
c_durable_products$plot_id  <- rep(1:landscape_scale_rotation_years, n_farms)
upscale_c_durable_products <- c_durable_products

c_hwp <- c_hwp %>% slice(rep(1:n(), each = landscape_scale_rotation_years))
n_farms <- floor(nrow(c_hwp)/landscape_scale_rotation_years)
c_hwp$farm <- rep(paste0("farm", 1:n_farms), each = landscape_scale_rotation_years)
c_hwp$plot_id  <- rep(1:landscape_scale_rotation_years, n_farms)
upscale_c_hwp <- c_hwp

c_beccs <- c_beccs %>% slice(rep(1:n(), each = landscape_scale_rotation_years))
n_farms <- floor(nrow(c_beccs)/landscape_scale_rotation_years)
c_beccs$farm <- rep(paste0("farm", 1:n_farms), each = landscape_scale_rotation_years)
c_beccs$plot_id  <- rep(1:landscape_scale_rotation_years, n_farms)
upscale_c_beccs <- c_beccs

AFC <- AFC %>% slice(rep(1:n(), each = landscape_scale_rotation_years))
n_farms <- floor(nrow(AFC)/landscape_scale_rotation_years)
AFC$farm <- rep(paste0("farm", 1:n_farms), each = landscape_scale_rotation_years)
AFC$plot_id  <- rep(1:landscape_scale_rotation_years, n_farms)
upscale_AFC <- AFC

AFClca <- AFClca %>% slice(rep(1:n(), each = landscape_scale_rotation_years))
n_farms <- floor(nrow(AFClca)/landscape_scale_rotation_years)
AFClca$farm <- rep(paste0("farm", 1:n_farms), each = landscape_scale_rotation_years)
AFClca$plot_id  <- rep(1:landscape_scale_rotation_years, n_farms)
upscale_AFClca <- AFClca

AFCbeccs <- AFCbeccs %>% slice(rep(1:n(), each = landscape_scale_rotation_years))
n_farms <- floor(nrow(AFCbeccs)/landscape_scale_rotation_years)
AFCbeccs$farm <- rep(paste0("farm", 1:n_farms), each = landscape_scale_rotation_years)
AFCbeccs$plot_id  <- rep(1:landscape_scale_rotation_years, n_farms)
upscale_AFCbeccs <- AFCbeccs
#2.2. Shift each plot (row) one column (year) to the right for as many years it will take to plant trees in that plot
#the objective is to develop an agroforestry system in which clear cuts of wood can be conducted every year
for (i in 1:landscape_scale_rotation_years){
  for (j in 1:landscape_scale_rotation_years){
    if (i != 1){
      upscale_AGCwoody[upscale_AGCwoody$plot_id == i, j+i-1] <- AGCwoody[AGCwoody$plot_id == i, j]
      upscale_BGCwoody[upscale_BGCwoody$plot_id == i, j+i-1] <- BGCwoody[BGCwoody$plot_id == i, j]
      upscale_growing_biom_c[upscale_growing_biom_c$plot_id == i, j+i-1] <- growing_biom_c[growing_biom_c$plot_id == i, j]
      upscale_SOC[upscale_SOC$plot_id == i, j+i-1] <- SOC[SOC$plot_id == i, j]
      upscale_c_durable_products[c_durable_products$plot_id == i, j+i-1] <- c_durable_products[AGCwoody$plot_id == i, j]
      upscale_pew_potential[pew_potential$plot_id == i, j+i-1] <- pew_potential[AGCwoody$plot_id == i, j]
      upscale_beccs_potential[beccs_potential$plot_id == i, j+i-1] <- beccs_potential[AGCwoody$plot_id == i, j]
      upscale_c_hwp[c_hwp$plot_id == i, j+i-1] <- c_hwp[AGCwoody$plot_id == i, j]
      upscale_c_beccs[c_beccs$plot_id == i, j+i-1] <- c_beccs[AGCwoody$plot_id == i, j]
      upscale_AFC[AFC$plot_id == i, j+i-1] <- AFC[AGCwoody$plot_id == i, j]
      upscale_AFClca[AFClca$plot_id == i, j+i-1] <- AFClca[AGCwoody$plot_id == i, j]
      upscale_AFCbeccs[AFCbeccs$plot_id == i, j+i-1] <- AFCbeccs[AGCwoody$plot_id == i, j]
    }
  }
}
for (i in (landscape_scale_rotation_years+1):(landscape_scale_rotation_years*2)){
  for (j in (landscape_scale_rotation_years+1):(landscape_scale_rotation_years*2)){
    if (i != landscape_scale_rotation_years+1){
      upscale_AGCwoody[upscale_AGCwoody$plot_id == i-landscape_scale_rotation_years, j+i-landscape_scale_rotation_years-1] <- AGCwoody[AGCwoody$plot_id == i-landscape_scale_rotation_years, j]
      upscale_BGCwoody[upscale_BGCwoody$plot_id == i-landscape_scale_rotation_years, j+i-landscape_scale_rotation_years-1] <- BGCwoody[BGCwoody$plot_id == i-landscape_scale_rotation_years, j]
      upscale_growing_biom_c[upscale_growing_biom_c$plot_id == i-landscape_scale_rotation_years, j+i-landscape_scale_rotation_years-1] <- growing_biom_c[growing_biom_c$plot_id == i-landscape_scale_rotation_years, j]
      upscale_SOC[upscale_SOC$plot_id == i-landscape_scale_rotation_years, j+i-landscape_scale_rotation_years-1] <- SOC[SOC$plot_id == i-landscape_scale_rotation_years, j]
      upscale_c_durable_products[upscale_c_durable_products$plot_id == i-landscape_scale_rotation_years, j+i-landscape_scale_rotation_years-2] <- c_durable_products[c_durable_products$plot_id == i-landscape_scale_rotation_years, j-1]
      upscale_pew_potential[upscale_pew_potential$plot_id == i-landscape_scale_rotation_years, j+i-landscape_scale_rotation_years-2] <- pew_potential[pew_potential$plot_id == i-landscape_scale_rotation_years, j-1]
      upscale_beccs_potential[upscale_beccs_potential$plot_id == i-landscape_scale_rotation_years, j+i-landscape_scale_rotation_years-2] <- beccs_potential[beccs_potential$plot_id == i-landscape_scale_rotation_years, j-1]
      upscale_c_hwp[upscale_c_hwp$plot_id == i-landscape_scale_rotation_years, j+i-landscape_scale_rotation_years-2] <- c_hwp[c_hwp$plot_id == i-landscape_scale_rotation_years, j-1]
      upscale_c_beccs[upscale_c_beccs$plot_id == i-landscape_scale_rotation_years, j+i-landscape_scale_rotation_years-2] <- c_beccs[c_beccs$plot_id == i-landscape_scale_rotation_years, j-1]
      upscale_AFCbeccs[upscale_AFCbeccs$plot_id == i-landscape_scale_rotation_years, j+i-landscape_scale_rotation_years-2] <- AFCbeccs[AFCbeccs$plot_id == i-landscape_scale_rotation_years, j-1]
      upscale_AFC[upscale_AFC$plot_id == i-landscape_scale_rotation_years, j+i-landscape_scale_rotation_years-2] <- AFC[AFC$plot_id == i-landscape_scale_rotation_years, j-1]
      upscale_AFClca[upscale_AFClca$plot_id == i-landscape_scale_rotation_years, j+i-landscape_scale_rotation_years-2] <- AFClca[AFClca$plot_id == i-landscape_scale_rotation_years, j-1]
    }
  }
}
for (i in (landscape_scale_rotation_years*2+1):n_years){
  for (j in (landscape_scale_rotation_years*2+1):n_years){
    if (i != landscape_scale_rotation_years*2+1){
      upscale_AGCwoody[upscale_AGCwoody$plot_id == i-landscape_scale_rotation_years*2, j+i-landscape_scale_rotation_years*2-1] <- AGCwoody[AGCwoody$plot_id == i-landscape_scale_rotation_years*2, j]
      upscale_BGCwoody[upscale_BGCwoody$plot_id == i-landscape_scale_rotation_years*2, j+i-landscape_scale_rotation_years*2-1] <- BGCwoody[BGCwoody$plot_id == i-landscape_scale_rotation_years*2, j]
      upscale_growing_biom_c[upscale_growing_biom_c$plot_id == i-landscape_scale_rotation_years*2, j+i-landscape_scale_rotation_years*2-1] <- growing_biom_c[growing_biom_c$plot_id == i-landscape_scale_rotation_years*2, j]
      upscale_SOC[upscale_SOC$plot_id == i-landscape_scale_rotation_years*2, j+i-landscape_scale_rotation_years*2-1] <- SOC[SOC$plot_id == i-landscape_scale_rotation_years*2, j]
      upscale_c_durable_products[upscale_c_durable_products$plot_id == i-landscape_scale_rotation_years*2, j+i-landscape_scale_rotation_years*2-2] <- c_durable_products[c_durable_products$plot_id == i-landscape_scale_rotation_years*2, j-1]
      upscale_pew_potential[upscale_pew_potential$plot_id == i-landscape_scale_rotation_years, j+i-landscape_scale_rotation_years-2] <- pew_potential[pew_potential$plot_id == i-landscape_scale_rotation_years, j-1]
      upscale_beccs_potential[upscale_beccs_potential$plot_id == i-landscape_scale_rotation_years, j+i-landscape_scale_rotation_years-2] <- beccs_potential[beccs_potential$plot_id == i-landscape_scale_rotation_years, j-1]
      upscale_c_hwp[upscale_c_durable_products$plot_id == i-landscape_scale_rotation_years*2, j+i-landscape_scale_rotation_years*2-2] <- c_hwp[c_hwp$plot_id == i-landscape_scale_rotation_years*2, j-1]
      upscale_c_beccs[upscale_c_durable_products$plot_id == i-landscape_scale_rotation_years*2, j+i-landscape_scale_rotation_years*2-2] <- c_beccs[c_beccs$plot_id == i-landscape_scale_rotation_years*2, j-1]
      upscale_AFCbeccs[upscale_c_durable_products$plot_id == i-landscape_scale_rotation_years*2, j+i-landscape_scale_rotation_years*2-2] <- AFCbeccs[AFCbeccs$plot_id == i-landscape_scale_rotation_years*2, j-1]
      upscale_AFC[upscale_c_durable_products$plot_id == i-landscape_scale_rotation_years*2, j+i-landscape_scale_rotation_years*2-2] <- AFC[AFC$plot_id == i-landscape_scale_rotation_years*2, j-1]
      upscale_AFClca[upscale_c_durable_products$plot_id == i-landscape_scale_rotation_years*2, j+i-landscape_scale_rotation_years*2-2] <- AFClca[AFClca$plot_id == i-landscape_scale_rotation_years*2, j-1]
    }
  }
}

upscale_AGCwoody <- upscale_AGCwoody[,1:n_years]
upscale_AGCwoody$plot_id  <- rep(1:landscape_scale_rotation_years, n_farms)
upscale_AGCwoody

upscale_BGCwoody <- upscale_BGCwoody[,1:n_years]
upscale_BGCwoody$plot_id  <- rep(1:landscape_scale_rotation_years, n_farms)
upscale_BGCwoody

upscale_growing_biom_c <- upscale_growing_biom_c[,1:n_years]
upscale_growing_biom_c$plot_id  <- rep(1:landscape_scale_rotation_years, n_farms)
upscale_growing_biom_c

upscale_SOC <- upscale_SOC[,1:n_years]
upscale_SOC$plot_id  <- rep(1:landscape_scale_rotation_years, n_farms)
upscale_SOC

upscale_pew_potential <- upscale_pew_potential[,1:n_years]
upscale_pew_potential$plot_id  <- rep(1:landscape_scale_rotation_years, n_farms)
upscale_pew_potential

upscale_beccs_potential <- upscale_beccs_potential[,1:n_years]
upscale_beccs_potential$plot_id  <- rep(1:landscape_scale_rotation_years, n_farms)
upscale_beccs_potential

upscale_c_durable_products <- upscale_c_durable_products[,1:n_years]
upscale_c_durable_products$plot_id  <- rep(1:landscape_scale_rotation_years, n_farms)
upscale_c_durable_products

upscale_c_hwp <- upscale_c_hwp[,1:n_years]
upscale_c_hwp$plot_id  <- rep(1:landscape_scale_rotation_years, n_farms)
upscale_c_hwp

upscale_c_beccs <- upscale_c_beccs[,1:n_years]
upscale_c_beccs$plot_id  <- rep(1:landscape_scale_rotation_years, n_farms)
upscale_c_beccs

upscale_AFC <- upscale_AFC[,1:n_years]
upscale_AFC$plot_id  <- rep(1:landscape_scale_rotation_years, n_farms)
upscale_AFC

upscale_AFClca <- upscale_AFClca[,1:n_years]
upscale_AFClca$plot_id  <- rep(1:landscape_scale_rotation_years, n_farms)
upscale_AFClca

upscale_AFCbeccs <- upscale_AFCbeccs[,1:n_years]
upscale_AFCbeccs$plot_id  <- rep(1:landscape_scale_rotation_years, n_farms)
upscale_AFCbeccs
#2.3 Turn the value to 0's in those cells which represent plots where trees have not yet been planted
boolean_mask <- col(upscale_AGCwoody) < upscale_AGCwoody$plot_id
boolean_mask <- ifelse(boolean_mask == TRUE, 0, 1)
colnames(boolean_mask) <- colnames(upscale_AGCwoody)
final_upscale_AGCwoody <- upscale_AGCwoody[,1:n_years] * boolean_mask[,1:n_years]
final_upscale_AGCwoody$plot_id  <- rep(1:landscape_scale_rotation_years, n_farms)
final_upscale_AGCwoody$farm <- rep(paste0("farm", 1:n_farms), each = landscape_scale_rotation_years)
write.csv(final_upscale_AGCwoody, paste0(folder, "/", "final_upscale_AGCwoody.csv"))

colnames(boolean_mask) <- colnames(upscale_BGCwoody)
final_upscale_BGCwoody <- upscale_BGCwoody[,1:n_years] * boolean_mask[,1:n_years]
final_upscale_BGCwoody$plot_id  <- rep(1:landscape_scale_rotation_years, n_farms)
final_upscale_BGCwoody$farm <- rep(paste0("farm", 1:n_farms), each = landscape_scale_rotation_years)
write.csv(final_upscale_BGCwoody, paste0(folder, "/", "final_upscale_BGCwoody.csv"))

colnames(boolean_mask) <- colnames(upscale_growing_biom_c)
final_upscale_growing_biom_c <- upscale_growing_biom_c[,1:n_years] * boolean_mask[,1:n_years]
final_upscale_growing_biom_c$plot_id  <- rep(1:landscape_scale_rotation_years, n_farms)
final_upscale_growing_biom_c$farm <- rep(paste0("farm", 1:n_farms), each = landscape_scale_rotation_years)
write.csv(final_upscale_growing_biom_c, paste0(folder, "/", "final_upscale_growing_biom_c.csv"))

colnames(boolean_mask) <- colnames(upscale_SOC)
final_upscale_SOC <- upscale_SOC[,1:n_years] * boolean_mask[,1:n_years]
final_upscale_SOC$plot_id  <- rep(1:landscape_scale_rotation_years, n_farms)
final_upscale_SOC$farm <- rep(paste0("farm", 1:n_farms), each = landscape_scale_rotation_years)
write.csv(final_upscale_SOC, paste0(folder, "/", "final_upscale_SOC.csv"))

boolean_mask_harvest <- col(upscale_pew_potential) < upscale_pew_potential$plot_id
boolean_mask_harvest <- ifelse(boolean_mask_harvest == TRUE, 0, 1)
colnames(boolean_mask) <- colnames(upscale_pew_potential)
final_upscale_pew_potential <- upscale_pew_potential[,1:n_years] * boolean_mask_harvest[,1:n_years]
final_upscale_pew_potential$plot_id  <- rep(1:landscape_scale_rotation_years, n_farms)
final_upscale_pew_potential$farm <- rep(paste0("farm", 1:n_farms), each = landscape_scale_rotation_years)
write.csv(final_upscale_pew_potential, paste0(folder, "/", "final_upscale_pew_potential.csv"))

boolean_mask_harvest <- col(upscale_beccs_potential) < upscale_beccs_potential$plot_id
boolean_mask_harvest <- ifelse(boolean_mask_harvest == TRUE, 0, 1)
colnames(boolean_mask) <- colnames(upscale_beccs_potential)
final_upscale_beccs_potential <- upscale_beccs_potential[,1:n_years] * boolean_mask_harvest[,1:n_years]
final_upscale_beccs_potential$plot_id  <- rep(1:landscape_scale_rotation_years, n_farms)
final_upscale_beccs_potential$farm <- rep(paste0("farm", 1:n_farms), each = landscape_scale_rotation_years)
write.csv(final_upscale_beccs_potential, paste0(folder, "/", "final_upscale_beccs_potential.csv"))

boolean_mask_harvest <- col(upscale_c_durable_products) < upscale_c_durable_products$plot_id
boolean_mask_harvest <- ifelse(boolean_mask_harvest == TRUE, 0, 1)
colnames(boolean_mask) <- colnames(upscale_c_durable_products)
final_upscale_c_durable_products <- upscale_c_durable_products[,1:n_years] * boolean_mask_harvest[,1:n_years]
final_upscale_c_durable_products$plot_id  <- rep(1:landscape_scale_rotation_years, n_farms)
final_upscale_c_durable_products$farm <- rep(paste0("farm", 1:n_farms), each = landscape_scale_rotation_years)
write.csv(final_upscale_c_durable_products, paste0(folder, "/", "final_upscale_c_durable_products.csv"))

boolean_mask_harvest <- col(upscale_c_hwp) < upscale_c_hwp$plot_id
boolean_mask_harvest <- ifelse(boolean_mask_harvest == TRUE, 0, 1)
colnames(boolean_mask) <- colnames(upscale_c_hwp)
final_upscale_c_hwp <- upscale_c_hwp[,1:n_years] * boolean_mask_harvest[,1:n_years]
final_upscale_c_hwp$plot_id  <- rep(1:landscape_scale_rotation_years, n_farms)
final_upscale_c_hwp$farm <- rep(paste0("farm", 1:n_farms), each = landscape_scale_rotation_years)
write.csv(final_upscale_c_hwp, paste0(folder, "/", "final_upscale_c_hwp.csv"))

boolean_mask_harvest <- col(upscale_c_beccs) < upscale_c_beccs$plot_id
boolean_mask_harvest <- ifelse(boolean_mask_harvest == TRUE, 0, 1)
colnames(boolean_mask) <- colnames(upscale_c_beccs)
final_upscale_c_beccs <- upscale_c_beccs[,1:n_years] * boolean_mask_harvest[,1:n_years]
final_upscale_c_beccs$plot_id  <- rep(1:landscape_scale_rotation_years, n_farms)
final_upscale_c_beccs$farm <- rep(paste0("farm", 1:n_farms), each = landscape_scale_rotation_years)
write.csv(final_upscale_c_beccs, paste0(folder, "/", "final_upscale_c_beccs.csv"))

boolean_mask_beccs <- col(upscale_beccs_potential) < upscale_beccs_potential$plot_id
boolean_mask_beccs <- ifelse(boolean_mask_beccs == TRUE, 0, 1)
colnames(boolean_mask) <- colnames(upscale_beccs_potential)
final_upscale_beccs <- upscale_beccs_potential[,1:n_years] * boolean_mask_beccs[,1:n_years]
final_upscale_beccs$plot_id  <- rep(1:landscape_scale_rotation_years, n_farms)
final_upscale_beccs$farm <- rep(paste0("farm", 1:n_farms), each = landscape_scale_rotation_years)

boolean_mask_AFC <- col(upscale_AFC) < upscale_AFC$plot_id
boolean_mask_AFC <- ifelse(boolean_mask_AFC == TRUE, 0, 1)
colnames(boolean_mask) <- colnames(upscale_AFC)
final_upscale_AFC <- upscale_AFC[,1:n_years] * boolean_mask_AFC[,1:n_years]
final_upscale_AFC$plot_id  <- rep(1:landscape_scale_rotation_years, n_farms)
final_upscale_AFC$farm <- rep(paste0("farm", 1:n_farms), each = landscape_scale_rotation_years)

boolean_mask_AFClca <- col(upscale_AFClca) < upscale_AFClca$plot_id
boolean_mask_AFClca <- ifelse(boolean_mask_AFClca == TRUE, 0, 1)
colnames(boolean_mask) <- colnames(upscale_AFClca)
final_upscale_AFClca <- upscale_AFClca[,1:n_years] * boolean_mask_AFClca[,1:n_years]
final_upscale_AFClca$plot_id  <- rep(1:landscape_scale_rotation_years, n_farms)
final_upscale_AFClca$farm <- rep(paste0("farm", 1:n_farms), each = landscape_scale_rotation_years)

boolean_mask_AFCbeccs <- col(upscale_AFCbeccs) < upscale_AFCbeccs$plot_id
boolean_mask_AFCbeccs <- ifelse(boolean_mask_AFCbeccs == TRUE, 0, 1)
colnames(boolean_mask) <- colnames(upscale_AFCbeccs)
final_upscale_AFCbeccs <- upscale_AFCbeccs[,1:n_years] * boolean_mask_AFCbeccs[,1:n_years]
final_upscale_AFCbeccs$plot_id  <- rep(1:landscape_scale_rotation_years, n_farms)
final_upscale_AFCbeccs$farm <- rep(paste0("farm", 1:n_farms), each = landscape_scale_rotation_years)
# 2.4 Calculate the carbon captured per hectare of "farmland" (the values expressed in the previous objects only represent carbon captured per hectare of agroforestry plot)
# remember that a farm includes both agroforestry and non-agroforestry plots (until the year equivalent to the length of the rotation. From then on, all the plots of the farm are agroforestry plots)
upscale_average_plot_results <- function (mcResults_farm, group, length_transition, simulation_time){
  outcome_per_hectare <- rowsum(mcResults_farm[,1:simulation_time], mcResults_farm[,group], reorder = FALSE)
  outcome_per_hectare[,1:length_transition] <- sweep(outcome_per_hectare[,1:length_transition], MARGIN=2, 1:length_transition, `/`)
  outcome_per_hectare[,(length_transition+1):simulation_time] <- outcome_per_hectare[,(length_transition+1):simulation_time] / length_transition
  return(outcome_per_hectare)
} 
upscale_plot_results <- function (regional_average, length_transition, simulation_time){
  outcome <- regional_average
  outcome[,1:length_transition] <- sweep(regional_average[,1:length_transition], MARGIN=2, 1:length_transition * annual_upscale_rate, `*`)
  outcome[,(length_transition+1):simulation_time] <- regional_average[,(length_transition+1):simulation_time] * length_transition * annual_upscale_rate
  return(outcome)
}
regional_average_CAI <- function (regional_stock, length_transition, annual_expansion, simulation_time){
  cumulative_outcome <- regional_stock
  cumulative_outcome[,2:n_years] <- regional_stock[,2:simulation_time] - regional_stock[,1:(simulation_time-1)]
  outcome_per_hectare <- cumulative_outcome
  outcome_per_hectare[,1:length_transition] <- sweep(cumulative_outcome[,1:length_transition], MARGIN=2, 1:length_transition * annual_expansion, `/`)
  outcome_per_hectare[,(length_transition+1):simulation_time] <- cumulative_outcome[,(length_transition+1):simulation_time]  / (annual_expansion*length_transition)
  return(outcome_per_hectare)
} 
regional_average_MAI <- function (regional_stock, length_transition, annual_expansion, simulation_time){
  cumulative_outcome <- sweep(regional_stock[,1:n_years], MARGIN=2, 1:n_years, `/`)
  outcome_per_hectare <- cumulative_outcome
  outcome_per_hectare[,1:length_transition] <- sweep(cumulative_outcome[,1:length_transition], 
                                                     MARGIN=2, 1:length_transition * annual_expansion, `/`)
  outcome_per_hectare[,(length_transition+1):simulation_time] <- cumulative_outcome[,(length_transition+1):simulation_time]  / (annual_expansion*length_transition)
  return(outcome_per_hectare)
} 

AGCwoody_stock_transition_average <- upscale_average_plot_results(final_upscale_AGCwoody, "farm", landscape_scale_rotation_years, n_years)
AGCwoody_stock_transition <- upscale_plot_results(AGCwoody_stock_transition_average, landscape_scale_rotation_years, n_years)
AGCwoody_CAC_transition <- AGCwoody_stock_transition
AGCwoody_CAC_transition[,2:n_years] <- AGCwoody_stock_transition[,2:n_years] - AGCwoody_stock_transition[,1:(n_years-1)]
AGCwoody_CAC_transition_average <- regional_average_CAI(AGCwoody_stock_transition, landscape_scale_rotation_years, annual_upscale_rate, n_years)
AGCwoody_MAC_transition <- sweep(AGCwoody_stock_transition[,1:n_years], MARGIN=2, 1:n_years, `/`)
AGCwoody_MAC_transition_average <- regional_average_MAI(AGCwoody_stock_transition, landscape_scale_rotation_years, annual_upscale_rate, n_years)

BGCwoody_stock_transition_average <- upscale_average_plot_results(final_upscale_BGCwoody, "farm", landscape_scale_rotation_years, n_years)
BGCwoody_stock_transition <- upscale_plot_results(BGCwoody_stock_transition_average, landscape_scale_rotation_years, n_years)
BGCwoody_CAC_transition <- BGCwoody_stock_transition
BGCwoody_CAC_transition[,2:n_years] <- BGCwoody_stock_transition[,2:n_years] - BGCwoody_stock_transition[,1:(n_years-1)]
BGCwoody_CAC_transition_average <- regional_average_CAI(BGCwoody_stock_transition, landscape_scale_rotation_years, annual_upscale_rate, n_years)
BGCwoody_MAC_transition <- sweep(BGCwoody_stock_transition[,1:n_years], MARGIN=2, 1:n_years, `/`)
BGCwoody_MAC_transition_average <- regional_average_MAI(BGCwoody_stock_transition, landscape_scale_rotation_years, annual_upscale_rate, n_years)

growing_biom_c_stock_transition_average <- upscale_average_plot_results(final_upscale_growing_biom_c, "farm", landscape_scale_rotation_years, n_years)
growing_biom_c_stock_transition <- upscale_plot_results(growing_biom_c_stock_transition_average, landscape_scale_rotation_years, n_years)
growing_biom_c_CAC_transition <- growing_biom_c_stock_transition
growing_biom_c_CAC_transition[,2:n_years] <- growing_biom_c_stock_transition[,2:n_years] - growing_biom_c_stock_transition[,1:(n_years-1)]
growing_biom_c_CAC_transition_average <- regional_average_CAI(growing_biom_c_stock_transition, landscape_scale_rotation_years, annual_upscale_rate, n_years)
growing_biom_c_MAC_transition <- sweep(growing_biom_c_stock_transition[,1:n_years], MARGIN=2, 1:n_years, `/`)
growing_biom_c_MAC_transition_average <- regional_average_MAI(growing_biom_c_stock_transition, landscape_scale_rotation_years, annual_upscale_rate, n_years)

SOC_stock_transition_average <- upscale_average_plot_results(final_upscale_SOC, "farm", landscape_scale_rotation_years, n_years)
SOC_stock_transition <- upscale_plot_results(SOC_stock_transition_average, landscape_scale_rotation_years, n_years)
SOC_CAC_transition <- SOC_stock_transition
SOC_CAC_transition[,2:n_years] <- SOC_stock_transition[,2:n_years] - SOC_stock_transition[,1:(n_years-1)]
SOC_CAC_transition_average <- regional_average_CAI(SOC_stock_transition, landscape_scale_rotation_years, annual_upscale_rate, n_years)
SOC_MAC_transition <- sweep(SOC_stock_transition[,1:n_years], MARGIN=2, 1:n_years, `/`)
SOC_MAC_transition_average <- regional_average_MAI(SOC_stock_transition, landscape_scale_rotation_years, annual_upscale_rate, n_years)

beccs_stock_transition_average <- upscale_average_plot_results(final_upscale_beccs_potential, "farm", landscape_scale_rotation_years, n_years)
beccs_stock_transition <- upscale_plot_results(beccs_stock_transition_average, landscape_scale_rotation_years, n_years)
beccs_CAC_transition <- beccs_stock_transition
beccs_CAC_transition[,2:n_years] <- beccs_stock_transition[,2:n_years] - beccs_stock_transition[,1:(n_years-1)]
beccs_CAC_transition_average <- regional_average_CAI(beccs_stock_transition, landscape_scale_rotation_years, annual_upscale_rate, n_years)
beccs_MAC_transition <- sweep(beccs_stock_transition[,1:n_years], MARGIN=2, 1:n_years, `/`)
beccs_MAC_transition_average <- regional_average_MAI(beccs_stock_transition, landscape_scale_rotation_years, annual_upscale_rate, n_years)

c_durable_products_stock_transition_average <- upscale_average_plot_results(final_upscale_c_durable_products, "farm", landscape_scale_rotation_years, n_years)
c_durable_products_stock_transition <- upscale_plot_results(c_durable_products_stock_transition_average, landscape_scale_rotation_years, n_years)
c_durable_products_CAC_transition <- c_durable_products_stock_transition
c_durable_products_CAC_transition[,2:n_years] <- c_durable_products_stock_transition[,2:n_years] - c_durable_products_stock_transition[,1:(n_years-1)]
c_durable_products_CAC_transition_average <- regional_average_CAI(c_durable_products_stock_transition, landscape_scale_rotation_years, annual_upscale_rate, n_years)
c_durable_products_MAC_transition <- sweep(c_durable_products_stock_transition[,1:n_years], MARGIN=2, 1:n_years, `/`)
c_durable_products_MAC_transition_average <- regional_average_MAI(c_durable_products_stock_transition, landscape_scale_rotation_years, annual_upscale_rate, n_years)

c_hwp_stock_transition_average <- upscale_average_plot_results(final_upscale_c_hwp, "farm", landscape_scale_rotation_years, n_years)
c_hwp_stock_transition <- upscale_plot_results(c_hwp_stock_transition_average, landscape_scale_rotation_years, n_years)
c_hwp_CAC_transition <- c_hwp_stock_transition
c_hwp_CAC_transition[,2:n_years] <- c_hwp_stock_transition[,2:n_years] - c_hwp_stock_transition[,1:(n_years-1)]
c_hwp_CAC_transition_average <- regional_average_CAI(c_hwp_stock_transition, landscape_scale_rotation_years, annual_upscale_rate, n_years)
c_hwp_MAC_transition <- sweep(c_hwp_stock_transition[,1:n_years], MARGIN=2, 1:n_years, `/`)
c_hwp_MAC_transition_average <- regional_average_MAI(c_hwp_stock_transition, landscape_scale_rotation_years, annual_upscale_rate, n_years)

AFC_stock_transition_average <- upscale_average_plot_results(final_upscale_AFC, "farm", landscape_scale_rotation_years, n_years)
AFC_stock_transition <- upscale_plot_results(AFC_stock_transition_average, landscape_scale_rotation_years, n_years)
AFC_CAC_transition <- AFC_stock_transition
AFC_CAC_transition[,2:n_years] <- AFC_stock_transition[,2:n_years] - AFC_stock_transition[,1:(n_years-1)]
AFC_CAC_transition_average <- regional_average_CAI(AFC_stock_transition, landscape_scale_rotation_years, annual_upscale_rate, n_years)
AFC_MAC_transition <- sweep(AFC_stock_transition[,1:n_years], MARGIN=2, 1:n_years, `/`)
AFC_MAC_transition_average <- regional_average_MAI(AFC_stock_transition, landscape_scale_rotation_years, annual_upscale_rate, n_years)

AFClca_stock_transition_average <- upscale_average_plot_results(final_upscale_AFClca, "farm", landscape_scale_rotation_years, n_years)
AFClca_stock_transition <- upscale_plot_results(AFClca_stock_transition_average, landscape_scale_rotation_years, n_years)
AFClca_CAC_transition <- AFClca_stock_transition
AFClca_CAC_transition[,2:n_years] <- AFClca_stock_transition[,2:n_years] - AFClca_stock_transition[,1:(n_years-1)]
AFClca_CAC_transition_average <- regional_average_CAI(AFClca_stock_transition, landscape_scale_rotation_years, annual_upscale_rate, n_years)
AFClca_MAC_transition <- sweep(AFClca_stock_transition[,1:n_years], MARGIN=2, 1:n_years, `/`)
AFClca_MAC_transition_average <- regional_average_MAI(AFClca_stock_transition, landscape_scale_rotation_years, annual_upscale_rate, n_years)

df_list <- list(AGC_stock = AGCwoody_stock_transition, 
                AGC_CAC = AGCwoody_CAC_transition, 
                AGC_MAC = AGCwoody_MAC_transition,
                AGC_stock_average = AGCwoody_stock_transition_average, 
                AGC_CAC_average = AGCwoody_CAC_transition_average, 
                AGC_MAC_average = AGCwoody_MAC_transition_average,
                BGC_stock = BGCwoody_stock_transition, 
                BGC_CAC = BGCwoody_CAC_transition, 
                BGC_MAC = BGCwoody_MAC_transition,
                BGC_stock_average = BGCwoody_stock_transition_average, 
                BGC_CAC_average = BGCwoody_CAC_transition_average, 
                BGC_MAC_average = BGCwoody_MAC_transition_average,
                growing_biom_c_stock = growing_biom_c_stock_transition, 
                growing_biom_c_CAC = growing_biom_c_CAC_transition, 
                growing_biom_c_MAC = growing_biom_c_MAC_transition,
                growing_biom_c_stock_average = growing_biom_c_stock_transition_average, 
                growing_biom_c_CAC_average = growing_biom_c_CAC_transition_average, 
                growing_biom_c_MAC_average = growing_biom_c_MAC_transition_average,
                SOC_stock = SOC_stock_transition, 
                SOC_CAC = SOC_CAC_transition, 
                SOC_MAC = SOC_MAC_transition,
                SOC_stock_average = SOC_stock_transition_average, 
                SOC_CAC_average = SOC_CAC_transition_average, 
                SOC_MAC_average = SOC_MAC_transition_average,
                beccs_stock = beccs_stock_transition, 
                beccs_CAC = beccs_CAC_transition, 
                beccs_MAC = beccs_MAC_transition,
                beccs_stock_average = beccs_stock_transition_average, 
                beccs_CAC_average = beccs_CAC_transition_average, 
                beccs_MAC_average = beccs_MAC_transition_average,
                c_durable_products_stock = c_durable_products_stock_transition, 
                c_durable_products_CAC = c_durable_products_CAC_transition, 
                c_durable_products_MAC = c_durable_products_MAC_transition,
                c_durable_products_stock_average = c_durable_products_stock_transition_average, 
                c_durable_products_CAC_average = c_durable_products_CAC_transition_average, 
                c_durable_products_MAC_average = c_durable_products_MAC_transition_average,
                c_hwp_stock = c_hwp_stock_transition, 
                c_hwp_CAC = c_hwp_CAC_transition, 
                c_hwp_MAC = c_hwp_MAC_transition,
                c_hwp_stock_average = c_hwp_stock_transition_average, 
                c_hwp_CAC_average = c_hwp_CAC_transition_average, 
                c_hwp_MAC_average = c_hwp_MAC_transition_average,
                AFC_stock = AFC_stock_transition, 
                AFC_CAC = AFC_CAC_transition, 
                AFC_MAC = AFC_MAC_transition,
                AFC_stock_average = AFC_stock_transition_average, 
                AFC_CAC_average = AFC_CAC_transition_average, 
                AFC_MAC_average = AFC_MAC_transition_average,
                AFClca_stock = AFClca_stock_transition, 
                AFClca_CAC = AFClca_CAC_transition, 
                AFClca_MAC = AFClca_MAC_transition,
                AFClca_stock_average = AFClca_stock_transition_average, 
                AFClca_CAC_average = AFClca_CAC_transition_average, 
                AFClca_MAC_average = AFClca_MAC_transition_average)
quant_stats <- lapply(df_list, function(x) x <- rbind(
  apply(x[, reporting_years],2,min), 
  apply(x[, reporting_years], 2, 
        quantile, probs = c(0.05, 0.25, 0.50, 0.75, 0.95), na.rm = TRUE ),
  apply(x[, reporting_years],2,max)))
for(i in 1:length(quant_stats)){
  write.csv(quant_stats[i], paste0(folder, "/", "quant_stats_", names(quant_stats)[i], ".csv"), row.names = c("Min", "5%", "25%", "50%", "75%", "95%", "Max"))
}

upscaled_mcSimulation_results <- mcSimulation_results
upscaled_mcSimulation_results$y[,paste0("AGwoodyB_carbon", 1:n_years)] <- AGCwoody_stock_transition
upscaled_mcSimulation_results$y[,paste0("BGwoodyB_carbon", 1:n_years)] <- BGCwoody_stock_transition
upscaled_mcSimulation_results$y[,paste0("total_woody_biomass_carbon", 1:n_years)] <- growing_biom_c_stock_transition
upscaled_mcSimulation_results$y[,paste0("soil_carbon_storage", 1:n_years)] <- SOC_stock_transition
upscaled_mcSimulation_results$y[,paste0("bccs", 1:n_years)] <- beccs_stock_transition
upscaled_mcSimulation_results$y[,paste0("c_products", 1:n_years)] <- c_durable_products_stock_transition
#upscaled_mcSimulation_results$y[,paste0("total_woody_biomass_carbon_beccs", 1:n_years)] <- c_beccs_stock_transition
upscaled_mcSimulation_results$y[,paste0("field_and_hwp", 1:n_years)] <- c_hwp_stock_transition
upscaled_mcSimulation_results$y[,paste0("AF_carbon", 1:n_years)] <- AFC_stock_transition
upscaled_mcSimulation_results$y[,paste0("lca_carbon", 1:n_years)] <- AFClca_stock_transition
#upscaled_mcSimulation_results$y[,paste0("bccs_lca_carbon", 1:n_years)] <- AFCbeccs_stock_transition
upscaled_mcSimulation_results$y[,paste0("AGwoodyB_carbon", 1:n_years)]
tail(upscaled_mcSimulation_results$y[,paste0("AGwoodyB_carbon", 1:n_years)])

upscaled_mcSimulation_results_per_hectare <- mcSimulation_results
upscaled_mcSimulation_results_per_hectare$y[,paste0("AGwoodyB_carbon", 1:n_years)] <- AGCwoody_stock_transition_average
upscaled_mcSimulation_results_per_hectare$y[,paste0("BGwoodyB_carbon", 1:n_years)] <- BGCwoody_stock_transition_average
upscaled_mcSimulation_results_per_hectare$y[,paste0("total_woody_biomass_carbon", 1:n_years)] <- growing_biom_c_stock_transition_average
upscaled_mcSimulation_results_per_hectare$y[,paste0("soil_carbon_storage", 1:n_years)] <- SOC_stock_transition_average
upscaled_mcSimulation_results_per_hectare$y[,paste0("bccs", 1:n_years)] <- beccs_stock_transition_average
upscaled_mcSimulation_results_per_hectare$y[,paste0("c_products", 1:n_years)] <- c_durable_products_stock_transition_average
#upscaled_mcSimulation_results_per_hectare$y[,paste0("total_woody_biomass_carbon_beccs", 1:n_years)] <- c_beccs_stock_transition
upscaled_mcSimulation_results_per_hectare$y[,paste0("field_and_hwp", 1:n_years)] <- c_hwp_stock_transition_average
upscaled_mcSimulation_results_per_hectare$y[,paste0("AF_carbon", 1:n_years)] <- AFC_stock_transition_average
upscaled_mcSimulation_results_per_hectare$y[,paste0("lca_carbon", 1:n_years)] <- AFClca_stock_transition_average
#upscaled_mcSimulation_results_per_hectare$y[,paste0("bccs_lca_carbon", 1:n_years)] <- AFCbeccs_stock_transition
upscaled_mcSimulation_results_per_hectare$y[,paste0("AGwoodyB_carbon", 1:n_years)]
tail(upscaled_mcSimulation_results_per_hectare$y[,paste0("AGwoodyB_carbon", 1:n_years)])

upscaled_mcSimulation_CAC <- mcSimulation_results
upscaled_mcSimulation_CAC$y[,paste0("AGwoodyB_carbon", 1:n_years)] <- AGCwoody_CAC_transition
upscaled_mcSimulation_CAC$y[,paste0("BGwoodyB_carbon", 1:n_years)] <- BGCwoody_CAC_transition
upscaled_mcSimulation_CAC$y[,paste0("total_woody_biomass_carbon", 1:n_years)] <- growing_biom_c_CAC_transition
upscaled_mcSimulation_CAC$y[,paste0("soil_carbon_storage", 1:n_years)] <- SOC_CAC_transition
upscaled_mcSimulation_CAC$y[,paste0("bccs", 1:n_years)] <- beccs_CAC_transition
upscaled_mcSimulation_CAC$y[,paste0("c_products", 1:n_years)] <- c_durable_products_CAC_transition
#upscaled_mcSimulation_CAC$y[,paste0("total_woody_biomass_carbon_beccs", 1:n_years)] <- c_beccs_MAC_transition
upscaled_mcSimulation_CAC$y[,paste0("field_and_hwp", 1:n_years)] <- c_hwp_CAC_transition
upscaled_mcSimulation_CAC$y[,paste0("AF_carbon", 1:n_years)] <- AFC_CAC_transition
upscaled_mcSimulation_CAC$y[,paste0("lca_carbon", 1:n_years)] <- AFClca_CAC_transition
#upscaled_mcSimulation_CAC$y[,paste0("bccs_lca_carbon", 1:n_years)] <- AFCbeccs_CAC_transition
upscaled_mcSimulation_CAC$y[,paste0("AGwoodyB_carbon", 1:n_years)]
tail(upscaled_mcSimulation_CAC$y[,paste0("AGwoodyB_carbon", 1:n_years)])

upscaled_mcSimulation_CAC_per_hectare <- mcSimulation_results
upscaled_mcSimulation_CAC_per_hectare$y[,paste0("AGwoodyB_carbon", 1:n_years)] <- AGCwoody_CAC_transition_average
upscaled_mcSimulation_CAC_per_hectare$y[,paste0("BGwoodyB_carbon", 1:n_years)] <- BGCwoody_CAC_transition_average
upscaled_mcSimulation_CAC_per_hectare$y[,paste0("total_woody_biomass_carbon", 1:n_years)] <- growing_biom_c_MAC_transition
upscaled_mcSimulation_CAC_per_hectare$y[,paste0("soil_carbon_storage", 1:n_years)] <- SOC_CAC_transition_average
upscaled_mcSimulation_CAC_per_hectare$y[,paste0("bccs", 1:n_years)] <- beccs_CAC_transition_average
upscaled_mcSimulation_CAC_per_hectare$y[,paste0("c_products", 1:n_years)] <- c_durable_products_CAC_transition_average
#upscaled_mcSimulation_CAC_per_hectare$y[,paste0("total_woody_biomass_carbon_beccs", 1:n_years)] <- c_beccs_MAC_transition
upscaled_mcSimulation_CAC_per_hectare$y[,paste0("field_and_hwp", 1:n_years)] <- c_hwp_CAC_transition_average
upscaled_mcSimulation_CAC_per_hectare$y[,paste0("AF_carbon", 1:n_years)] <- AFC_CAC_transition_average
upscaled_mcSimulation_CAC_per_hectare$y[,paste0("lca_carbon", 1:n_years)] <- AFClca_CAC_transition_average
#upscaled_mcSimulation_CAC_per_hectare$y[,paste0("bccs_lca_carbon", 1:n_years)] <- AFCbeccs_CAC_transition
upscaled_mcSimulation_CAC_per_hectare$y[,paste0("AGwoodyB_carbon", 1:n_years)]
tail(upscaled_mcSimulation_CAC_per_hectare$y[,paste0("AGwoodyB_carbon", 1:n_years)])

upscaled_mcSimulation_MAC <- mcSimulation_results
upscaled_mcSimulation_MAC$y[,paste0("AGwoodyB_carbon", 1:n_years)] <- AGCwoody_MAC_transition
upscaled_mcSimulation_MAC$y[,paste0("BGwoodyB_carbon", 1:n_years)] <- BGCwoody_MAC_transition
upscaled_mcSimulation_MAC$y[,paste0("total_woody_biomass_carbon", 1:n_years)] <- growing_biom_c_MAC_transition
upscaled_mcSimulation_MAC$y[,paste0("soil_carbon_storage", 1:n_years)] <- SOC_MAC_transition
upscaled_mcSimulation_MAC$y[,paste0("bccs", 1:n_years)] <- beccs_MAC_transition
upscaled_mcSimulation_MAC$y[,paste0("c_products", 1:n_years)] <- c_durable_products_MAC_transition
#upscaled_mcSimulation_MAC$y[,paste0("total_woody_biomass_carbon_beccs", 1:n_years)] <- c_beccs_CAC_transition
upscaled_mcSimulation_MAC$y[,paste0("field_and_hwp", 1:n_years)] <- c_hwp_MAC_transition
upscaled_mcSimulation_MAC$y[,paste0("AF_carbon", 1:n_years)] <- AFC_MAC_transition
upscaled_mcSimulation_MAC$y[,paste0("lca_carbon", 1:n_years)] <- AFClca_MAC_transition
#upscaled_mcSimulation_MAC$y[,paste0("bccs_lca_carbon", 1:n_years)] <- AFCbeccs_MAC_transition
upscaled_mcSimulation_MAC$y[,paste0("AGwoodyB_carbon", 1:n_years)]
tail(upscaled_mcSimulation_MAC$y[,paste0("AGwoodyB_carbon", 1:n_years)])

upscaled_mcSimulation_MAC_per_hectare <- mcSimulation_results
upscaled_mcSimulation_MAC_per_hectare$y[,paste0("AGwoodyB_carbon", 1:n_years)] <- AGCwoody_MAC_transition_average
upscaled_mcSimulation_MAC_per_hectare$y[,paste0("BGwoodyB_carbon", 1:n_years)] <- BGCwoody_MAC_transition_average
upscaled_mcSimulation_MAC_per_hectare$y[,paste0("total_woody_biomass_carbon", 1:n_years)] <- growing_biom_c_MAC_transition_average
upscaled_mcSimulation_MAC_per_hectare$y[,paste0("soil_carbon_storage", 1:n_years)] <- SOC_MAC_transition_average
upscaled_mcSimulation_MAC_per_hectare$y[,paste0("bccs", 1:n_years)] <- beccs_MAC_transition_average
upscaled_mcSimulation_MAC_per_hectare$y[,paste0("c_products", 1:n_years)] <- c_durable_products_MAC_transition_average
#upscaled_mcSimulation_MAC_per_hectare$y[,paste0("total_woody_biomass_carbon_beccs", 1:n_years)] <- c_beccs_CAC_transition
upscaled_mcSimulation_MAC_per_hectare$y[,paste0("field_and_hwp", 1:n_years)] <- c_hwp_MAC_transition_average
upscaled_mcSimulation_MAC_per_hectare$y[,paste0("AF_carbon", 1:n_years)] <- AFC_MAC_transition_average
upscaled_mcSimulation_MAC_per_hectare$y[,paste0("lca_carbon", 1:n_years)] <- AFClca_MAC_transition_average
#upscaled_mcSimulation_MAC_per_hectare$y[,paste0("bccs_lca_carbon", 1:n_years)] <- AFCbeccs_MAC_transition
upscaled_mcSimulation_MAC_per_hectare$y[,paste0("AGwoodyB_carbon", 1:n_years)]
tail(upscaled_mcSimulation_MAC_per_hectare$y[,paste0("AGwoodyB_carbon", 1:n_years)])