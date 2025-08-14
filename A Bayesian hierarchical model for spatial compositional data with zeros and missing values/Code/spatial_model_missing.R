############################################################.
## SPATIAL MULTIVARIATE GAM NIMBLE MODEL - MISSING VALUES ##
############################################################.


## PARALLEL CHAINS ##


## PACKAGES ####

library(dplyr)
library(tidyr)
library(nimble)
library(nimbleHMC)
library(coda)
library(mcmcplots)
library(ggplot2)
library(mgcv)
library(parallel)
library(magrittr)


#######################.

## SET WORKING DIRECTORY ##

setwd("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD Papers : Journals/Spatial/Statistical Modelling Journal/Spatial/")

source("Code/beta-binomial_nimblefunction.R")


#######################.

## DATA ####

tree_data <- readRDS("Data/tree_counts_1000.rds") 


#### SPATIAL LOCATIONS ####

spatial_locations <- tree_data %>%
  select(X_coord, Y_coord)


### TREE TYPES ####

tree_types <- c("larch",
                "oak",
                "sitka_spruce",
                "sycamore")


tree_types_data <- tree_data %>% 
  select(X_coord,
         Y_coord,
         all_of(tree_types),
         total) 


###################.

## RANDOMLY SELECT ROWS ####

N <- nrow(tree_types_data)

N_fit <- 1000

set.seed(451810)
fit_rows <- sample(1:nrow(tree_types_data),N_fit,replace = FALSE)


###################.

## SETUP MISSING DATA ####

selected_data <- tree_types_data[fit_rows,]

missing_tree_data <- selected_data

n_missing <- rep(0,N_fit)


#### SELECT DATA - 1 MISSING ####

missing_1_row <- sample(1:length(fit_rows),200,
                        replace = FALSE)
n_missing[missing_1_row] <- 1

#### SELECT DATA - 2 MISSING ####

missing_2_rows <- sample(setdiff(1:length(fit_rows),
                                 missing_1_row),
                         100,
                         replace = FALSE)
n_missing[missing_2_rows] <- 2

#### SELECT DATA - 3 MISSING ####

missing_3_rows <- sample(setdiff(1:length(fit_rows),
                                 c(missing_1_row,missing_2_rows)),
                         100,
                         replace = FALSE)
n_missing[missing_3_rows] <- 3


#### CREATE NA DATA ####

for(i in 1:N_fit){
  if(n_missing[i]>0){
    missing_tree_data[i,sample(3:6,n_missing[i],replace=FALSE)] <- NA
  }
}

missing_tree_data <- missing_tree_data %>%
  mutate(n_missing = n_missing) %>%
  mutate(across(c(larch:sycamore),
                ~ if_else(is.na(.), TRUE, FALSE),
                .names = "{.col}_missing"))


#### PLOT ####

missing_tree_data %>%
  pivot_longer(larch:sycamore, 
               names_to = "type",
               values_to = "count") %>% 
  mutate(type = factor(type, levels = tree_types)) %>%
  ggplot(., aes(x = X_coord, y = Y_coord, fill = count)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  facet_wrap(~type) +
  labs(x = "x coordinate", y = "y coordinate") +
  theme(legend.position = "bottom") +
  ggtitle("Heatmap for Tree Types - Count")


#### TOTAL COUNT ####

y <- rep(1000,N_fit)


#### TREE COUNTS ####

x_missing = missing_tree_data %>%
  select(larch:sycamore)

x_actual = selected_data %>%
  select(larch:sycamore)


###################.

## CONSTANTS ####

# NUMBER OF TREE TYPES
N_types = length(tree_types)

# NUMBER OF BASIS
N_basis = 400

#######################.

## GAM MODEL ####

dummy_data <- mutate(missing_tree_data,dummy=rbinom(N_fit,1000,0.5)/1000)


larch_gam_bin <- gam(dummy ~ s(X_coord, Y_coord, k=N_basis,bs="tp"), # CHANGED
                     family = "binomial",
                     weights = rep(1000,N_fit),
                     data = dummy_data)

Z <- predict(larch_gam_bin,
             missing_tree_data,
             type="lpmatrix")

Z_full <- predict(larch_gam_bin,
                  tree_types_data,
                  type="lpmatrix")

larch_jagam_bin <- jagam(dummy ~ s(X_coord, Y_coord, k=N_basis,bs="tp"), # CHANGED
                         family = "binomial",
                         weights = rep(1000,N_fit),
                         data = dummy_data,
                         file = "larch.jags")

S1 <- larch_jagam_bin$jags.data$S1

larch_jagam_bin$jags.data$X-Z

#######################.

## NIMBLE MODEL ####

# NUMBER OF CORES
n_cores <- min(detectCores(), 100, 20)

# NUMBER OF CHAINS
nchains = 4

# NUMBER OF CLUSTERS
this_cluster <- makeCluster(nchains)

### SEED ####

initial_seed <- round(runif(nchains, 1234, 8750874))
initial_seed

chain_seed <- round(runif(nchains, 5468, 5632751))
chain_seed


### MODEL FUNCTION ####

run_MCMC <- function(X,
                     initial_seed, chain_seed,
                     x_missing, y,
                     missing_data,
                     Z, S1,
                     N, N_types, N_basis,
                     niter, nburn, nthin) {
  
  library(nimble)
  
  source("beta-binomial_nimblefunction.R")
  
  #### MODEL CODE ####
  
  spatial_model_code <- nimbleCode({
    
    ## EACH TREE TYPE ##
    
    for (k in 1:N_types) {
      
      ## PRIORS ##
      
      for (i in 1:2) {
        sigma[i,k] ~ T(dnorm(0,sd=100),0,)
        lambda[i,k] <- 1/sigma[i,k]^2
        rho[i,k] <- log(lambda[i,k])
        
      }
      
      K1[1:(N_basis - 1), 1:(N_basis-1), k] <- S1[1:(N_basis - 1), 1:(N_basis - 1)] * lambda[1,k] + S1[1:(N_basis - 1), N_basis:((N_basis - 1)*2)] * lambda[2,k]
      
      
      b[1,k] ~ dnorm(0, sd = 10)
      b[2:N_basis,k] ~ dmnorm(zero[1:(N_basis-1)], K1[1:(N_basis - 1), 1:(N_basis-1), k])
      
      
      ## LINEAR PREDICTOR ##
      eta[1:N, k] <- Z[1:N, 1:N_basis] %*% b[1:N_basis, k]
      
      
      ## EXPECTED RESPONSE ##
      mu[1:N, k] <- expit(eta[1:N, k])
      
      ## MODEL COUNTS X ##
      mu_mean[k] <- mean(mu[1:N, k])
      
      for (i in 1:N) {
        
        log(phi[i, k]) <- psi[1, k] +
        psi[2, k]*(mu[i, k]-mu_mean[k]) +
        psi[3, k]*(mu[i, k]-mu_mean[k])^2 +
        psi[4, k]*(mu[i, k]-mu_mean[k])^3
        
      }
      
      psi[1, k] ~ dnorm(2,sd=2)
      
      for(j in 2:4){
        psi[j, k] ~ dnorm(0,sd=1)
      }
      
    }
    
    for (i in 1:N) {
      
      # FIRST TYPE
      x[i, 1] ~ dbetabinomial(mu[i, 1], phi[i, 1], y[i])
      
      for (k in 2:N_types) {
        
        x[i, k] ~ dbetabinomial(mu[i, k], phi[i, k], y[i] - sum(x[i, 1:(k-1)]))
        
      }
    }
    
  })
  
  
  #### DATA ####
  data <- list(x = x_missing)
  
  
  #### CONSTANTS ####
  constants <- list(y = y,
                    Z = Z,
                    S1 = S1,
                    N = N,
                    N_basis = N_basis,
                    N_types = N_types,
                    zero = rep(0, (N_basis - 1)))
  
  
  #### INITIAL VALUES ####
  init_fun <- function(){
    
    x_inits <- matrix(NA, nrow = N, ncol = N_types) 
    for(i in 1:N){
      if(is.na(sum(x_missing[i,]))){
        n_missing <- sum(is.na(x_missing[i,]))
        x_inits[i,is.na(x_missing[i,])] <- rmulti(1,y-sum(x_missing[i,which(!is.na(x_missing[i,]))]),
                                                  prob=rep(1/(n_missing+1),n_missing+1))[1:n_missing]
      }
    }
    
    inits <- list(b = matrix(rnorm(N_basis, 0, 0.05), nrow = N_basis, ncol = N_types),
                  sigma = matrix(runif(2*N_types,0,1), nrow = 2, ncol = N_types),
                  psi = matrix(rnorm(4*N_types,rep(c(2,0,0,0),N_types),0.25), nrow = 4, ncol = N_types),
                  x = x_inits)
    return(inits)
  }
  
  set.seed(initial_seed[X])
  initial_values <- init_fun()
  
  
  #### BUILD MODEL ####
  model <- nimbleModel(code = spatial_model_code,
                       data = data,
                       constants = constants,
                       inits = init_fun())
  
  
  #### COMPILE MODEL ####
  compile_model <- compileNimble(model)
  
  compile_model$calculate()
  
  
  #### CONFIGURE MCMC ####
  conf_model <- configureMCMC(compile_model, 
                              monitors = c("x",
                                           "b",
                                           "lambda",
                                           "sigma",
                                           "rho",
                                           "eta",
                                           "mu",
                                           "phi",
                                           "psi"
                              ),
                              print = TRUE, 
                              useConjugacy = FALSE)
  
  
  #### SAMPLERS ####
  conf_model$removeSamplers(c('sigma', "b","psi"))

  for(j in 1:N_types){
     conf_model$addSampler(target = paste0("b[1:",N_basis,",",j,"]"), type = "AF_slice")
  }
  for(j in 1:N_types){
    conf_model$addSampler(target = paste0("psi[1:",4,",",j,"]"), type = "AF_slice")
  }
  for (node in 1:length(model$expandNodeNames("sigma"))) {
    conf_model$addSampler(target = model$expandNodeNames("sigma")[node], type = "slice")
  }
  
  #### BUILD MCMC ####
  spatial_model <- buildMCMC(conf_model)
  
  
  #### COMPILE MCMC ####
  compile_spatial_model <- compileNimble(spatial_model,
                                         project = model,
                                         resetFunctions = TRUE)
  
  #### RUN MODEL ####
  samples <- runMCMC(compile_spatial_model, 
                     niter =   niter, 
                     nburnin = nburn,
                     thin =   nthin, 
                     samplesAsCodaMCMC = TRUE,
                     setSeed = chain_seed[X])
  
  #return(samples)
  return(list(initial_values = initial_values, samples = samples))
  
}

# system.time({
#   test <- run_MCMC(X=1,
#                    initial_seed = initial_seed,
#                    chain_seed = chain_seed,
#                    x_missing = x_missing,
#                    y = y,
#                    Z = Z,
#                    S1 = S1,
#                    N = N_fit,
#                    N_types = N_types,
#                    N_basis = N_basis,
#                    niter = 200,
#                    nburn = 100,
#                    nthin =   1)
# })




### RUN MODEL ####
run_time <- system.time({samples <- parLapply(cl = this_cluster, 
                                              X = 1:nchains,
                                              fun = run_MCMC, 
                                              
                                              initial_seed = initial_seed,
                                              chain_seed = chain_seed,
                                              
                                              x_missing = x_missing, 
                                              y = y,
                                              
                                              Z = Z,
                                              S1 = S1,
                                              
                                              N = N_fit, 
                                              N_types = N_types,
                                              N_basis = N_basis,
                                              
                                              niter = 20000, 
                                              nburn = 10000, 
                                              nthin =   2)})




# CLOSE CLUSTER
stopCluster(this_cluster)

run_time


#######################.

## OUTPUT ####

save(samples, file = paste0('output/spatial_model_missing_output_', Sys.Date(), '_polynomial_phi.RData'))


samples_list <- lapply(samples, function(chain) chain$samples)


# COMBINE ALL CHAINS
output <- as_tibble(do.call('rbind', samples_list))

N_samples <- nrow(output)

samples_mcmc <- as.mcmc.list(samples_list)

### GELMAN ####

psrf_all <- sapply(1:ncol(output),function(x)gelman.diag(samples_mcmc[,x],
            transform = TRUE,
            autoburnin = FALSE,
            multivariate = FALSE)$psrf)
psrf_all%<>%t()

mean(psrf_all[,1]<=1.05,na.rm=T)


subset_b <- function(chain) {
  cols <- grep("^b", colnames(chain), value = TRUE)
  return(chain[, cols])
}

subset_lambda <- function(chain) {
  cols <- grep("^lambda", colnames(chain), value = TRUE)
  return(chain[, cols])
}

subset_phi <- function(chain) {
  cols <- grep("^phi", colnames(chain), value = TRUE)
  return(chain[, cols])
}

subset_sigma <- function(chain) {
  cols <- grep("^sigma", colnames(chain), value = TRUE)
  return(chain[, cols])
}

subset_psi <- function(chain) {
  cols <- grep("^psi", colnames(chain), value = TRUE)
  return(chain[, cols])
}

subset_mu <- function(chain) {
  cols <- grep("^mu", colnames(chain), value = TRUE)
  return(chain[, cols])
}

#######################.

## PREDICTION ####

x_samples <- array(NA, dim = c(N_fit, N_types, N_samples))
for (i in 1:N_samples) {
  x_row <- output[i, ] %>% select(starts_with("x")) 
  x_samples[, , i] <- array(unlist(x_row), dim = c(N_fit, N_types))
}
x_samples


#######################.

## SAVE PREDICTIONS ####

saveRDS(x_samples, "output/gdm_predictions_polynomial_phi.rds")


#######################.

b_samples <- subset_b(output)%>%unlist()%>%array(dim=c(N_samples,N_basis,N_types))

phi_samples <- subset_phi(output)%>%unlist()%>%array(dim=c(N_samples,N_fit,N_types))

sigma_samples <- subset_sigma(output)%>%unlist()%>%array(dim=c(N_samples,2,N_types))

mu_samples_full <- array(NA,dim=c(N_samples,N,N_types))
for(d in 1:N_types){
  mu_samples_full[,,d] <- expit(b_samples[,,d]%*%t(Z_full))
}

saveRDS(mu_samples_full, "output/mu_samples_full_polynomial.rds")

x_mean_samples_full <- array(NA,dim=c(N_samples,N,N_types))
x_mean_samples_full[,,1] <- mu_samples_full[,,1]
x_mean_samples_full[,,2] <- mu_samples_full[,,2]*(1-x_mean_samples_full[,,1])
for(d in 3:N_types){
  x_mean_samples_full[,,d] <- mu_samples_full[,,d]*(1-apply(x_mean_samples_full[,,1:(d-1)],c(1,2),sum))
}

saveRDS(x_mean_samples_full, "output/x_mean_samples_full_polynomial.rds")


#######################.

## COMPARISON MODELS ####

larch_gam <- gam(larch / total ~ s(X_coord, Y_coord, k = N_basis),
                 family = quasibinomial(link="logit"),
                 weights = rep(1000,N_fit),
                 data = missing_tree_data,
                 optimizer = "efs",
                 method="GCV.Cp")

oak_gam <- gam(oak / total ~ s(X_coord, Y_coord, k = N_basis),
               family = quasibinomial(link="logit"),
               weights = rep(1000,N_fit),
               data = missing_tree_data,
               optimizer = "efs",
               method="GCV.Cp")

sitka_spruce_gam <- gam(sitka_spruce / total ~ s(X_coord, Y_coord, k = N_basis),
                        family = quasibinomial(link="logit"),
                        weights = rep(1000,N_fit),
                        data = missing_tree_data,
                        optimizer = "efs",
                        method="GCV.Cp")


sycamore_gam <- gam(sycamore / total ~ s(X_coord, Y_coord, k = N_basis),
                    family = quasibinomial(link="logit"),
                    weights = rep(1000,N_fit),
                    data = missing_tree_data,
                    optimizer = "efs",
                    method="GCV.Cp")


### PREDICT COMPARISON ####

larch_predict <- predict(larch_gam,
                         newdata = tree_types_data,
                         type = "response")

oak_predict <- predict(oak_gam,
                       newdata = tree_types_data,
                       type = "response")

sitka_spruce_predict <- predict(sitka_spruce_gam,
                                newdata = tree_types_data,
                                type = "response")

sycamore_predict <- predict(sycamore_gam,
                            newdata = tree_types_data,
                            type = "response")


predicted_counts <- cbind(larch = larch_predict,
                          oak = oak_predict,
                          sitka_spruce = sitka_spruce_predict,
                          sycamore = sycamore_predict)


saveRDS(predicted_counts,file="output/predictions_GAM.rds")


#######################.