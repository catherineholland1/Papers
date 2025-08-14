#####################################.
## SIMULATE SPATIAL DATA FOR PAPER ##
#####################################.

## PACKAGES ####

library(dplyr)
library(mgcv)
library(compositions)
library(writexl)

##################################.

set.seed(123)

setwd("Data/")

##################################.

## CREATE SQUARE SPATITAL GRID ####

n_side <- 32  
# 32 x 32 = 1024 points

grid <- expand.grid(
  x = seq(0, 1, length.out = n_side),
  y = seq(0, 1, length.out = n_side)
)

n_locations <- nrow(grid)
n_components <- 10  

##################################.

## SIMULATE SPATIAL SURFACES FOR EACH COMPONENT ####

fitted_surfaces <- list()
sim_data_matrix <- matrix(NA, nrow = n_locations, ncol = n_components)

for (i in 1:n_components) {
  # Generate smooth random surface (GAM with thin-plate spline)
  z <- rnorm(n_locations)  # add some noise
  gam_model <- gam(z ~ s(x, y, bs = "tp"), data = cbind(grid, z = z), family = gaussian())
  
  # Simulate one realization from the GAM
  sim_vals <- simulate(gam_model, nsim = 1)[, 1]
  
  # Store simulated raw values
  sim_data_matrix[, i] <- sim_vals
}

##################################.

## APPLY INV CLR TO CONVERT TO PROPORTIONS ####

proportions <- t(apply(sim_data_matrix, 1, function(row) {
  inv_clr <- exp(row) / sum(exp(row))
  return(inv_clr)
}))


##################################.

## CONSTRUCT DATA ####

simulated_data <- data.frame(
  grid,
  proportions
)

# Rename columns to match trees data
names(simulated_data) <- c('X_coord', 'Y_coord', 
                           'Ash_pct', 'Beech_pct', 'Larch_pct', 'Oak_pct', 'ScPine_pct', 
                           'Shadow_pct', 'SBirch_pct', 'SitSpr_pct', 'SwChes_pct', 'Syca_pct')

# Order columns to match trees data
simulated_data <- simulated_data %>%
  select(ends_with('_pct'),
         ends_with('_coord'))

head(simulated_data)
rowSums(simulated_data[,-c(1:2)])

##################################.

## SAVE ####

simulated_data %>%
  write_xlsx(., "simulated_trees.xlsx")

##################################.