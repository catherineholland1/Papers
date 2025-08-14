##############################.
## COMPARISON  - GDM VS GAM ##
##############################.

## PACKAGES ####

library(dplyr)
library(tidyr)
library(ggplot2)


#######################.

## SET WORKING DIRECTORY ####

setwd("~/PhD/Spatial - Paper/")

#######################.

## MISSING TREE DATA ####

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


N_types <- length(tree_types)


## RANDOMLY SELECT ROWS ####

N <- nrow(tree_types_data)

N_fit <- 1000

set.seed(451810)
fit_rows <- sample(1:nrow(tree_types_data),N_fit,replace = FALSE)


## SETUP MISSING DATA ####

selected_data <- tree_types_data[fit_rows,]

selected_data %>%
  select(larch:sycamore) %>%
  pivot_longer(., cols = larch:sycamore, values_to = "count", names_to = "type") %>%
  group_by(type) %>%
  summarise(
    mean = mean(count, na.rm = T)
  )

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


#### DATA SUMMARY ####

missing_tree_data %>%
  select(larch:sycamore) %>%
  pivot_longer(., cols = larch:sycamore, values_to = "count", names_to = "type") %>%
  group_by(type) %>%
  summarise(
    n_zero = sum(count == 0, na.rm = TRUE),
    perc_zero = mean(count == 0, na.rm = TRUE) * 100,
    n_na = sum(is.na(count)),
    perc_na = mean(is.na(count)) * 100,
    mean = mean(count, na.rm = T)
  )


missing_tree_data %>%
  group_by(n_missing) %>%
  summarise(n = n()) 


missing_tree_plot <- missing_tree_data %>%
  pivot_longer(larch:sycamore, 
               names_to = "type",
               values_to = "count") %>% 
  mutate(type = case_when(type == 'larch' ~ 'Larch',
                          type == 'oak' ~ 'Oak',
                          type == 'sitka_spruce' ~ 'Sitka spruce',
                          type == 'sycamore' ~ 'Sycamore',
                          TRUE ~ type)) %>%
  mutate(type = factor(type, levels = c('Larch', 'Oak', 'Sitka spruce', 'Sycamore'))) 



missing_tree_plot %>%
  ggplot(aes(x = X_coord, y = Y_coord, fill = count)) +
  geom_tile(data = missing_tree_plot %>% filter(!is.na(count))) +
  scale_fill_viridis_c() +
  
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  facet_wrap(~type) +
  
  geom_tile(data = missing_tree_plot %>% filter(is.na(count)), aes(colour = ""),
            linetype = 0, fill = "#FC8D62") +
  
  labs(x = "x coordinate", y = "y coordinate",
       fill = "count", colour = "missing count") +
  theme(legend.position = "bottom") +
  ggtitle("Heatmap for Modelled Tree Species with Missing Counts") 


# pdf("output/heatmap_missing.pdf", height = 4.75)
# missing_tree_plot %>%
#   ggplot(aes(x = X_coord, y = Y_coord, fill = count)) +
#   geom_tile(data = missing_tree_plot %>% filter(!is.na(count))) +
#   scale_fill_viridis_c() +
#   
#   theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
#   facet_wrap(~type) +
#   
#   geom_tile(data = missing_tree_plot %>% filter(is.na(count)), aes(colour = ""),
#             linetype = 0, fill = "#FC8D62") +
#   
#   labs(x = "", y = "",
#        fill = "Count", colour = "Missing Count") +
#   
#   theme(legend.position = "bottom",
#         legend.margin = margin(t=-15)) 
# dev.off()


#######################.

tree_proportions <- readRDS("Data/tree_proportions.rds")

## LOAD PREDICTIONS ####

### GDM ####

gdm_predictions <- readRDS("output/gdm_predictions_polynomial_phi.rds") %>%
  
  apply(., c(1,2), mean) %>%
  
  as.data.frame() %>%
  rename(larch = "V1",
         oak = "V2",
         sitka_spruce = "V3",
         sycamore = "V4") %>%
  cbind(missing_tree_data %>% select(c(X_coord, Y_coord,
                                       contains("missing")))) %>%
  pivot_longer(., cols = larch:sycamore, 
               names_to = "type", values_to = "predicted") %>%
  left_join(select(tree_proportions,X_coord,Y_coord,larch,oak,sitka_spruce,sycamore) %>% pivot_longer(., cols = larch:sycamore, 
                                            names_to = "type", values_to = "actual"), by = c("X_coord", "Y_coord", "type")) %>%
  
  mutate(predicted = predicted/1000) %>%
  
  mutate(model = "GDM")

x_mean_samples_full <- readRDS("output/x_mean_samples_full_polynomial.rds")
mu_samples_full <- readRDS("output/mu_samples_full_polynomial.rds")

excluded_x_pred_df <- apply(x_mean_samples_full,c(2,3),mean)%>%as.data.frame()%>%
  rename(larch = "V1",
         oak = "V2",
         sitka_spruce = "V3",
         sycamore = "V4") %>%
  cbind(tree_types_data %>% select(c(X_coord, Y_coord,total))) %>%
  anti_join(select(missing_tree_data,X_coord,Y_coord))%>%
  mutate(n_missing=4,larch_missing=TRUE,oak_missing=TRUE,sitka_spruce_missing=TRUE,sycamore_missing=TRUE)%>%
  pivot_longer(., cols = larch:sycamore, 
               names_to = "type", values_to = "predicted") %>%
  left_join(select(tree_proportions,X_coord,Y_coord,larch,oak,sitka_spruce,sycamore) %>% pivot_longer(., cols = larch:sycamore, 
                                                                                                      names_to = "type", values_to = "actual"), by = c("X_coord", "Y_coord", "type")) %>%
  
  mutate(model="GDM")

gdm_predictions <- full_join(excluded_x_pred_df,gdm_predictions)

### GAM ####

gam_predictions_full <- readRDS("output/predictions_GAM.rds")

gam_predictions_400 <- gam_predictions_full[fit_rows,] %>%
  
  as.data.frame() %>%
  # rename(larch = 'V1',
  #        oak = 'V2',
  #        sitka_spruce = 'V3',
  #        sycamore = 'V4') %>%
  cbind(missing_tree_data %>% select(contains("coord"),
                                     contains("missing")))%>%
  pivot_longer(., cols = larch:sycamore, 
               names_to = "type", values_to = "predicted") %>%
  left_join(select(tree_proportions,X_coord,Y_coord,larch,oak,sitka_spruce,sycamore) %>% pivot_longer(., cols = larch:sycamore, 
                                                                                                      names_to = "type", values_to = "actual"), by = c("X_coord", "Y_coord", "type")) %>%
  
  mutate(predicted=predicted,model = "GAM_400")



excluded_gam_preds <- gam_predictions_full%>%as.data.frame()%>%
  cbind(tree_types_data %>% select(c(X_coord, Y_coord,total))) %>%
  anti_join(select(missing_tree_data,X_coord,Y_coord))%>%
  mutate(n_missing=4,larch_missing=TRUE,oak_missing=TRUE,sitka_spruce_missing=TRUE,sycamore_missing=TRUE)%>%
  pivot_longer(., cols = larch:sycamore, 
               names_to = "type", values_to = "predicted") %>%
  left_join(select(tree_proportions,X_coord,Y_coord,larch,oak,sitka_spruce,sycamore) %>% 
              pivot_longer(., cols = larch:sycamore,names_to = "type", values_to = "actual"),
            by = c("X_coord", "Y_coord", "type")) %>%
  mutate(model="GAM_400")

gam_predictions_400 <- full_join(excluded_gam_preds,gam_predictions_400)


#######################.

## COMPARISON DATA ####

comparison_data <- gdm_predictions %>%
  rbind(gam_predictions_400) %>%
  arrange(X_coord, Y_coord)

comparison_data%<>%mutate(naive=actual)%>%group_by(model)%>%
  group_modify(function(x,y){
    x$naive[x$type=="larch"&x$larch_missing] <- mean(x$actual[x$type=="larch"&!x$larch_missing])
    x$naive[x$type=="oak"&x$oak_missing] <- mean(x$actual[x$type=="oak"&!x$oak_missing])
    x$naive[x$type=="sitka_spruce"&x$sitka_spruce_missing] <- mean(x$actual[x$type=="sitka_spruce"&!x$sitka_spruce_missing])
    x$naive[x$type=="sycamore"&x$sycamore_missing] <- mean(x$actual[x$type=="sycamore"&!x$sycamore_missing])
    return(x)
  })

#######################.

## COMPARE PREDICITIONS ####

### PLOT ####

#### PREDICTIONS VS ACTUAL ####

comparison_data %>%
  
  filter((type == "larch" & larch_missing == TRUE) |
           (type == "oak" & oak_missing == TRUE) |
           (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
           (type == "sycamore" & sycamore_missing == TRUE)) %>%
  
  ggplot(., aes(x = actual, y = predicted)) +
  geom_point() +
  
  geom_abline(col = "red") +
  
  facet_grid(type*n_missing ~ model) +
  
  ggtitle("Predicted vs Actual")


comparison_data %>%
  
  filter(model == "GDM" | model == "GAM_400") %>%
  mutate(model = if_else(model == "GAM_400", "GAM", model)) %>%
  
  filter((type == "larch" & larch_missing == TRUE) |
           (type == "oak" & oak_missing == TRUE) |
           (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
           (type == "sycamore" & sycamore_missing == TRUE)) %>%
  
  ggplot(., aes(x = actual, y = predicted)) +
  geom_point() +
  
  geom_abline(col = "red") +
  
  facet_grid(type*n_missing ~ model) +
  
  ggtitle("Predicted vs Actual")


pdf("output/scatterplot_predictions_actual.pdf", height = 4.5)
comparison_data %>%
  filter(model == "GDM" | model == "GAM_400",n_missing<=3) %>%
  mutate(model = if_else(model == "GAM_400", "QB GAMs", model),
         n_missing=case_match(n_missing,
                              1~"1 species missing (n = 200)",
                              2~"2 species missing (n = 400)",
                              3~"3 species missing (n = 600)",
                              4~"All species missing (n = 9224)")) %>%
  
  filter((type == "larch" & larch_missing == TRUE) |
           (type == "oak" & oak_missing == TRUE) |
           (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
           (type == "sycamore" & sycamore_missing == TRUE)) %>%
  mutate(type = case_when(type == 'larch' ~ 'Larch',
                          type == 'oak' ~ 'Oak',
                          type == 'sitka_spruce' ~ 'Sitka spruce',
                          type == 'sycamore' ~ 'Sycamore',
                          TRUE ~ type)) %>%
  
  ggplot(., aes(x = predicted, y = actual,shape=type)) +
  scale_x_continuous(limits=c(0,1))+scale_y_continuous(limits=c(0,1))+
  geom_abline() +
  geom_point() +
  labs(x="Predicted proportion",y="Observed proportion",shape=NULL)+
  facet_grid(model ~ n_missing)+
  #coord_fixed()+
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(vjust=0.8))+
  scale_colour_viridis_d(begin=0,end=0.9)+
  scale_fill_viridis_d(begin=0,end=0.9)+
  scale_shape_manual(values=c(0,1,2,5))
dev.off()


for (i in 1:N_types){
  
  pdf(paste0("output/scatterplot_predictions_actual_", tree_types[i], ".pdf"), height = 5)
  
  plot(comparison_data %>%
         
         filter(model == "GDM" | model == "GAM_400") %>%
         mutate(model = if_else(model == "GAM_400", "GAM", model)) %>%
         
         filter((type == "larch" & larch_missing == TRUE) |
                  (type == "oak" & oak_missing == TRUE) |
                  (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
                  (type == "sycamore" & sycamore_missing == TRUE)) %>%
         
         filter(type == tree_types[i]) %>%
         
         ggplot(., aes(x = actual, y = predicted)) +
         
         geom_abline(col = "darkgrey") +
         
         geom_point() +
         
         facet_grid(~n_missing ~ model) +
         
         labs(x = 'Observed Count', y = 'Predicted Count')
  )
  dev.off()
}


for (i in 1:N_types){
  
  pdf(paste0("output/scatterplot_predictions_actual_", tree_types[i], "-wide.pdf"), height = 3)
  
  plot(comparison_data %>%
         
         filter(model == "GDM" | model == "GAM_400") %>%
         mutate(model = if_else(model == "GAM_400", "GAM", model)) %>%
         mutate(model = factor(model, levels=c('GDM', 'GAM'))) %>%
         
         filter((type == "larch" & larch_missing == TRUE) |
                  (type == "oak" & oak_missing == TRUE) |
                  (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
                  (type == "sycamore" & sycamore_missing == TRUE)) %>%
         
         filter(type == tree_types[i]) %>%
         
         ggplot(., aes(x = actual, y = predicted)) +
         
         geom_abline(col = "darkgrey") +
         
         geom_point() +
         
         facet_grid(model~n_missing) +
         
         labs(x = 'Observed Count', y = 'Predicted Count')
  )
  dev.off()
}



#### RESIDUALS VS PREDICTIONS ####

comparison_data %>%
  
  filter((type == "larch" & larch_missing == TRUE) |
           (type == "oak" & oak_missing == TRUE) |
           (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
           (type == "sycamore" & sycamore_missing == TRUE)) %>%
  
  mutate(residuals = actual - predicted) %>%
  
  ggplot(., aes(x = predicted, y = residuals)) +
  geom_point() +
  facet_grid(type*n_missing ~ model) +
  ggtitle("Residuals vs Predicted")


#### RESID VS ACTUAL ####

comparison_data %>%
  
  filter((type == "larch" & larch_missing == TRUE) |
           (type == "oak" & oak_missing == TRUE) |
           (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
           (type == "sycamore" & sycamore_missing == TRUE)) %>%
  
  mutate(residuals = actual - predicted) %>%
  
  ggplot(., aes(x = actual, y = residuals)) +
  geom_point() +
  facet_grid(type*n_missing ~ model) +
  ggtitle("Residuals vs Actual")


#### HEATMAP ####

comparison_data %>%
  
  filter((type == "larch" & larch_missing == TRUE) |
           (type == "oak" & oak_missing == TRUE) |
           (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
           (type == "sycamore" & sycamore_missing == TRUE)) %>%
  
  pivot_longer(., cols = actual:predicted, values_to = "count", names_to = "value_type") %>%
  mutate(type = case_when(type == 'larch' ~ 'Larch',
                          type == 'oak' ~ 'Oak',
                          type == 'sitka_spruce' ~ 'Sitka spruce',
                          type == 'sycamore' ~ 'Sycamore',
                          TRUE ~ type)) %>%
  mutate(type = factor(type, levels = c('Larch', 'Oak', 'Sitka spruce', 'Sycamore'))) %>%
  
  mutate(model = if_else(value_type == "actual", "actual", model)) %>%
  
  ggplot(., aes(x = X_coord, y = Y_coord, fill = count)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  
  facet_grid(type ~ model) +
  
  labs(x = "x coordinate", y = "y coordinate") +
  theme(legend.position = "bottom") +
  ggtitle("Heatmap for Tree Types - Actual vs Predicted")


comparison_data %>%
  
  filter(model == "GDM" | model == "GAM_400") %>%
  mutate(model = if_else(model == "GAM_400", "GAM", model)) %>%
  
  filter((type == "larch" & larch_missing == TRUE) |
           (type == "oak" & oak_missing == TRUE) |
           (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
           (type == "sycamore" & sycamore_missing == TRUE)) %>%
  
  pivot_longer(., cols = actual:predicted, values_to = "count", names_to = "value_type") %>%
  mutate(type = factor(type, levels = tree_types)) %>%
  
  mutate(model = if_else(value_type == "actual", "actual", model)) %>%
  
  ggplot(., aes(x = X_coord, y = Y_coord, fill = count)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  
  facet_grid(type ~ model) +
  
  labs(x = "x coordinate", y = "y coordinate") +
  theme(legend.position = "bottom") +
  ggtitle("Heatmap for Tree Types - Actual vs Predicted")



pdf("output/heatmap_predictions_actual.pdf", height = 5.75)
comparison_data %>%
  
  filter(model == "GDM" | model == "GAM_400") %>%
  mutate(model = if_else(model == "GAM_400", "GAM", model)) %>%
  
  filter((type == "larch" & larch_missing == TRUE) |
           (type == "oak" & oak_missing == TRUE) |
           (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
           (type == "sycamore" & sycamore_missing == TRUE)) %>%
  
  pivot_longer(., cols = actual:predicted, values_to = "count", names_to = "value_type") %>%
  mutate(type = case_when(type == 'larch' ~ 'Larch',
                          type == 'oak' ~ 'Oak',
                          type == 'sitka_spruce' ~ 'Sitka spruce',
                          type == 'sycamore' ~ 'Sycamore',
                          TRUE ~ type)) %>%
  mutate(type = factor(type, levels = c('Larch', 'Oak', 'Sitka spruce', 'Sycamore'))) %>%
  
  mutate(model = if_else(value_type == "actual", "Observed Counts", model)) %>%
  
  mutate(model = factor(model, levels = c('Observed Counts', 'GDM', 'GAM'))) %>%
  mutate(type = case_when(type == 'sitka_spruce' ~ 'sitka spruce',
                          TRUE ~ type)) %>%
  
  ggplot(., aes(x = X_coord, y = Y_coord, fill = count)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme(axis.text.x=element_text(angle=60, hjust=1)) +
  
  facet_grid(type ~ model) +
  
  labs(x = "", y = "", fill = 'Count') +
  theme(legend.position = "bottom",
        legend.margin = margin(t=-15)) 
dev.off()


for (i in 1:N_types){
  
  pdf(paste0("output/heatmap_predictions_actual-", tree_types[i], ".pdf"), height = 3)
  plot(comparison_data %>%
         
         filter(model == "GDM" | model == "GAM_400") %>%
         mutate(model = if_else(model == "GAM_400", "GAM", model)) %>%
         
         filter(type == tree_types[i]) %>%
         
         # filter((type == "larch" & larch_missing == TRUE) |
         #          (type == "oak" & oak_missing == TRUE) |
         #          (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
         #          (type == "sycamore" & sycamore_missing == TRUE)) %>%
         
         pivot_longer(., cols = actual:predicted, values_to = "count", names_to = "value_type") %>%
         mutate(type = case_when(type == 'larch' ~ 'Larch',
                                 type == 'oak' ~ 'Oak',
                                 type == 'sitka_spruce' ~ 'Sitka spruce',
                                 type == 'sycamore' ~ 'Sycamore',
                                 TRUE ~ type)) %>%
         #mutate(type = factor(type, levels = c('Larch', 'Oak', 'Sitka spruce', 'Sycamore'))) %>%
         
         mutate(model = if_else(value_type == "actual", "Observed Counts", model)) %>%
         
         mutate(model = factor(model, levels = c('Observed Counts', 'GDM', 'GAM'))) %>%
         mutate(type = case_when(type == 'sitka_spruce' ~ 'sitka spruce',
                                 TRUE ~ type)) %>%
         
         ggplot(., aes(x = X_coord, y = Y_coord, fill = count)) +
         geom_tile() +
         scale_fill_viridis_c() +
         theme(axis.text.x=element_text(angle=60, hjust=1)) +
         
         facet_grid(type ~ model) +
         
         labs(x = "", y = "", fill = 'Count') +
         theme( axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
               
               legend.position = "bottom") 
  )
  dev.off()
}


#### DENSITY ####

comparison_data %>%
  
  filter((type == "larch" & larch_missing == TRUE) |
           (type == "oak" & oak_missing == TRUE) |
           (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
           (type == "sycamore" & sycamore_missing == TRUE)) %>%
  
  group_by(type, n_missing, model) %>%  
  mutate(mean_actual = mean(actual, na.rm = TRUE)) %>%
  ungroup() %>% 
  
  ggplot(aes(x = log(predicted), fill = type)) +
  
  geom_density(aes(color = type), alpha = 0.8) +
  
  facet_grid(type*n_missing ~ model, scales = "free") +
  
  geom_vline(aes(xintercept = log(mean_actual), 
                 color = "Original Data"),  
             linetype = "dashed", 
             linewidth = 1) +
  
  scale_color_manual(name = "", values = c("Original Data" = "#FC8D62")) +
  scale_fill_manual(name = "Tree Type", values = c("larch" = "#66C2A5",
                                                   "oak" = "#8DA0CB",
                                                   "sitka_spruce" = "#E78AC3",
                                                   "sycamore" = "#FFD92F")) +
  
  labs(title = "Density Plot of Predicted and Actual Counts - GAM Model",
       x = "Log of Predicted Counts", y = "Density") +
  theme(legend.position = "bottom")


comparison_data %>%
  
  filter(model == "GDM" | model == "GAM_400") %>%
  mutate(model = if_else(model == "GAM_400", "GAM", model)) %>%
  
  filter((type == "larch" & larch_missing == TRUE) |
           (type == "oak" & oak_missing == TRUE) |
           (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
           (type == "sycamore" & sycamore_missing == TRUE)) %>%
  
  group_by(type, n_missing, model) %>%  
  mutate(mean_actual = mean(actual, na.rm = TRUE)) %>%
  ungroup() %>% 
  
  ggplot(aes(x = log(predicted), fill = type)) +
  
  geom_density(aes(color = type), alpha = 0.8) +
  
  facet_grid(type*n_missing ~ model, scales = "free") +
  
  geom_vline(aes(xintercept = log(mean_actual), 
                 color = "Original Data"),  
             linetype = "dashed", 
             linewidth = 1) +
  
  scale_color_manual(name = "", values = c("Original Data" = "#FC8D62")) +
  scale_fill_manual(name = "Tree Type", values = c("larch" = "#66C2A5",
                                                   "oak" = "#8DA0CB",
                                                   "sitka_spruce" = "#E78AC3",
                                                   "sycamore" = "#FFD92F")) +
  
  labs(title = "Density Plot of Predicted and Actual Counts - GAM Model",
       x = "Log of Predicted Counts", y = "Density") +
  theme(legend.position = "bottom")


pdf("output/density_predictions_actual.pdf", height = 10, width = 8)
comparison_data %>%
  
  filter(model == "GDM" | model == "GAM_400") %>%
  mutate(model = if_else(model == "GAM_400", "GAM", model)) %>%
  
  filter((type == "larch" & larch_missing == TRUE) |
           (type == "oak" & oak_missing == TRUE) |
           (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
           (type == "sycamore" & sycamore_missing == TRUE)) %>%
  
  group_by(type, n_missing, model) %>%  
  mutate(mean_actual = mean(actual, na.rm = TRUE)) %>%
  ungroup() %>% 
  
  ggplot(aes(x = log(predicted), fill = type)) +
  
  geom_density(aes(color = type), alpha = 0.8) +
  
  facet_grid(type*n_missing ~ model, scales = "free") +
  
  geom_vline(aes(xintercept = log(mean_actual), 
                 color = "Original Data"),  
             linetype = "dashed", 
             linewidth = 1) +
  
  scale_color_manual(name = "", values = c("Original Data" = "#FC8D62")) +
  scale_fill_manual(name = "Tree Type", values = c("larch" = "#66C2A5",
                                                   "oak" = "#8DA0CB",
                                                   "sitka_spruce" = "#E78AC3",
                                                   "sycamore" = "#FFD92F")) +
  
  labs(title = "Density Plot of Predicted and Actual Counts - GAM Model",
       x = "Log of Predicted Counts", y = "Density") +
  theme(legend.position = "bottom")
dev.off()


### QUANTIFY ####

#### MAE ####

comparison_data %>%
  
  filter((type == "larch" & larch_missing == TRUE) |
           (type == "oak" & oak_missing == TRUE) |
           (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
           (type == "sycamore" & sycamore_missing == TRUE)) %>%
  
  group_by(type, n_missing, model) %>%
  
  summarise(mae = mean(100*abs(predicted - actual))) %>%
  
  pivot_wider(names_from = model, values_from = mae) 


comparison_data %>%
  
  filter((type == "larch" & larch_missing == TRUE) |
           (type == "oak" & oak_missing == TRUE) |
           (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
           (type == "sycamore" & sycamore_missing == TRUE)) %>%
  
  group_by(n_missing, model) %>%
  
  summarise(mae = mean(100*abs(predicted - actual))) %>%
  
  pivot_wider(names_from = model, values_from = mae) 


#### RMSE #####

comparison_data %>%
  
  filter((type == "larch" & larch_missing == TRUE) |
           (type == "oak" & oak_missing == TRUE) |
           (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
           (type == "sycamore" & sycamore_missing == TRUE)) %>%
  
  group_by(type, n_missing, model) %>%
  
  summarise(rmse = 100*sqrt(mean((predicted - actual)^2))) %>%
  
  pivot_wider(names_from = model, values_from = rmse) 


comparison_data %>%
  
  filter((type == "larch" & larch_missing == TRUE) |
           (type == "oak" & oak_missing == TRUE) |
           (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
           (type == "sycamore" & sycamore_missing == TRUE)) %>%
  
  group_by(n_missing, model) %>%
  
  summarise(rmse = 100*sqrt(mean((predicted - actual)^2))) %>%
  
  pivot_wider(names_from = model, values_from = rmse) 




#### MSE #####

comparison_data %>%
  
  filter(model == 'GDM' | model == 'GAM_400') %>%
  mutate(model = if_else(model == 'GAM_400', 'GAM', model)) %>%
  
  filter((type == "larch" & larch_missing == TRUE) |
           (type == "oak" & oak_missing == TRUE) |
           (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
           (type == "sycamore" & sycamore_missing == TRUE)) %>%
  
  group_by(type, n_missing, model) %>%
  
  summarise(mse = mean((predicted - actual)^2),mse0=mean((naive-actual)^2),
            xi=1-mse/mse0) %>% select(-mse,-mse0)%>%
  
  pivot_wider(names_from = model, values_from = xi) 


comparison_data %>%
  
  filter(model == 'GDM' | model == 'GAM_400') %>%
  mutate(model = if_else(model == 'GAM_400', 'GAM', model)) %>%
  
  filter((type == "larch" & larch_missing == TRUE) |
           (type == "oak" & oak_missing == TRUE) |
           (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
           (type == "sycamore" & sycamore_missing == TRUE)) %>%
  
  group_by(n_missing, model) %>%
  
  summarise(mse = mean((predicted - actual)^2),mse0=mean((naive-actual)^2),
            xi=1-mse/mse0) %>% select(-mse,-mse0)%>%
  
  pivot_wider(names_from = model, values_from = xi) 



comparison_data %>%
  
  filter(model == 'GDM' | model == 'GAM_400') %>%
  mutate(model = if_else(model == 'GAM_400', 'GAM', model)) %>%
  
  filter((type == "larch" & larch_missing == TRUE) |
           (type == "oak" & oak_missing == TRUE) |
           (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
           (type == "sycamore" & sycamore_missing == TRUE)) %>%
  
  group_by(n_missing, model) %>%
  
  summarise(mse = mean((predicted - actual)^2),mse0=mean((naive-actual)^2),
            xi=1-mse/mse0) %>% select(-mse,-mse0)%>%
  
  pivot_wider(names_from = model, values_from = xi)

comparison_data %>%
  
  filter(model == 'GDM' | model == 'GAM_400',n_missing<2,n_missing>0) %>%
  mutate(model = if_else(model == 'GAM_400', 'GAM', model)) %>%
  
  filter((type == "larch" & larch_missing == TRUE) |
           (type == "oak" & oak_missing == TRUE) |
           (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
           (type == "sycamore" & sycamore_missing == TRUE)) %>%
  
  group_by(type, model) %>%
  
  summarise(mse = mean((predicted - actual)^2),mse0=mean((naive-actual)^2),
            xi=1-mse/mse0) %>% select(-mse,-mse0)%>%
  
  pivot_wider(names_from = model, values_from = xi)

#######################.

gdm_predictions_full_df <- apply(x_mean_samples_full,c(2,3),mean)%>%as.data.frame()%>%
  rename(larch = "V1",
         oak = "V2",
         sitka_spruce = "V3",
         sycamore = "V4") %>%
  cbind(tree_types_data %>% select(c(X_coord, Y_coord,total)))%>%
  pivot_longer(., cols = larch:sycamore, 
               names_to = "type", values_to = "predicted") %>%
  left_join(select(tree_proportions,X_coord,Y_coord,larch,oak,sitka_spruce,sycamore) %>% pivot_longer(., cols = larch:sycamore, 
                                                                                                      names_to = "type", values_to = "actual"), by = c("X_coord", "Y_coord", "type")) %>%
  
  mutate(model="GDM")

gam_predictions_full_df <- gam_predictions_full%>%as.data.frame()%>%
  cbind(tree_types_data %>% select(c(X_coord, Y_coord,total)))%>%
  pivot_longer(., cols = larch:sycamore, 
               names_to = "type", values_to = "predicted") %>%
  left_join(select(tree_proportions,X_coord,Y_coord,larch,oak,sitka_spruce,sycamore) %>% pivot_longer(., cols = larch:sycamore, 
                                                                                                      names_to = "type", values_to = "actual"), by = c("X_coord", "Y_coord", "type")) %>%
  
  mutate(model="QB")



comparison_data_means <- gdm_predictions_full_df %>%
  rbind(gam_predictions_full_df) %>%
  arrange(X_coord, Y_coord)


for (i in 1:N_types){
  
  pdf(paste0("output/heatmap_predictions_means-", tree_types[i], ".pdf"), height = 3)
  plot(comparison_data_means %>%
         
         filter(model == "GDM" | model == "QB") %>%
         mutate(model = if_else(model == "QB", "QB", model)) %>%
         
         filter(type == tree_types[i]) %>%
         
         # filter((type == "larch" & larch_missing == TRUE) |
         #          (type == "oak" & oak_missing == TRUE) |
         #          (type == "sitka_spruce" & sitka_spruce_missing == TRUE) |
         #          (type == "sycamore" & sycamore_missing == TRUE)) %>%
         
         pivot_longer(., cols = actual:predicted, values_to = "count", names_to = "value_type") %>%
         mutate(type = case_when(type == 'larch' ~ 'Larch',
                                 type == 'oak' ~ 'Oak',
                                 type == 'sitka_spruce' ~ 'Sitka spruce',
                                 type == 'sycamore' ~ 'Sycamore',
                                 TRUE ~ type)) %>%
         #mutate(type = factor(type, levels = c('Larch', 'Oak', 'Sitka spruce', 'Sycamore'))) %>%
         
         mutate(model = if_else(value_type == "actual", "Observed counts", model)) %>%
         
         mutate(model = factor(model, levels = c('Observed counts', 'GDM', 'QB'))) %>%
         mutate(type = case_when(type == 'sitka_spruce' ~ 'sitka spruce',
                                 TRUE ~ type)) %>%
         
         ggplot(., aes(x = X_coord, y = Y_coord, fill = count)) +
         geom_tile() +
         scale_fill_viridis_c() +
         theme(axis.text.x=element_text(angle=60, hjust=1)) +
         
         facet_grid(type ~ model) +
         
         labs(x = "", y = "", fill = 'Count') +
         theme( axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                
                legend.position = "bottom") 
  )
  dev.off()
}

ggplot(gdm_predictions_full_df%>%inner_join(select(filter(missing_tree_data,n_missing==0),X_coord,Y_coord)),
       aes(x=predicted,y=actual))+
  geom_abline()+geom_point()+geom_smooth()

ggplot(gam_predictions_full_df%>%inner_join(select(filter(missing_tree_data,n_missing==0),X_coord,Y_coord)),
       aes(x=predicted,y=actual))+
  geom_abline()+geom_point()+geom_smooth()


pdf("output/heatmap_predictions_means-larch-oak.pdf", height = 4.5)
  plot(comparison_data_means %>%
         
         filter(model == "GDM" | model == "QB") %>%
         mutate(model = if_else(model == "QB", "QB", model)) %>%
         
         filter(type%in%c("larch","oak")) %>%
         
         pivot_longer(., cols = actual:predicted, values_to = "count", names_to = "value_type") %>%
         mutate(type = case_when(type == 'larch' ~ 'Larch',
                                 type == 'oak' ~ 'Oak',
                                 type == 'sitka_spruce' ~ 'Sitka spruce',
                                 type == 'sycamore' ~ 'Sycamore',
                                 TRUE ~ type)) %>%
         #mutate(type = factor(type, levels = c('Larch', 'Oak', 'Sitka spruce', 'Sycamore'))) %>%
         
         mutate(model = if_else(value_type == "actual", "Original data", model)) %>%
         
         mutate(model = factor(model, levels = c('Original data', 'GDM', 'QB'))) %>%
         mutate(type = case_when(type == 'sitka_spruce' ~ 'sitka spruce',
                                 TRUE ~ type)) %>%
         
         ggplot(., aes(x = X_coord, y = Y_coord, fill = count/100,colour=count/100)) +
         geom_tile() +
         scale_fill_viridis_c(name="Proportion",breaks=seq(0,1,by=0.2)) +scale_colour_viridis_c(name="Proportion",breaks=seq(0,1,by=0.2)) +
         theme_bw()+
         theme(axis.text.x=element_text(angle=60, hjust=1)) +
         coord_fixed()+
         facet_grid(type ~ model) +
         labs(x = "", y = "", fill = 'Count') +
         theme( axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                panel.grid = element_blank(),
                legend.position = "bottom",
                legend.key.width = unit(0.5, "in"),
                legend.title = element_text(vjust=0.8)) 
  )
  dev.off()

  
  pdf("output/heatmap_predictions_means-all.pdf", height = 8)
  plot(comparison_data_means %>%
         
         filter(model == "GDM" | model == "QB") %>%
         mutate(model = if_else(model == "QB", "QB", model)) %>%
         
         filter(type%in%c("larch","oak","sitka_spruce","sycamore")) %>%
         
         pivot_longer(., cols = actual:predicted, values_to = "count", names_to = "value_type") %>%
         mutate(type = case_when(type == 'larch' ~ 'Larch',
                                 type == 'oak' ~ 'Oak',
                                 type == 'sitka_spruce' ~ 'Sitka spruce',
                                 type == 'sycamore' ~ 'Sycamore',
                                 TRUE ~ type)) %>%
         #mutate(type = factor(type, levels = c('Larch', 'Oak', 'Sitka spruce', 'Sycamore'))) %>%
         
         mutate(model = if_else(value_type == "actual", "Original data", model)) %>%
         
         mutate(model = factor(model, levels = c('Original data', 'GDM', 'QB'))) %>%
         mutate(type = case_when(type == 'sitka_spruce' ~ 'sitka spruce',
                                 TRUE ~ type)) %>%
         
         ggplot(., aes(x = X_coord, y = Y_coord, fill = count/100,colour=count/100)) +
         geom_tile() +
         scale_fill_viridis_c(name="Proportion",breaks=seq(0,1,by=0.2)) +scale_colour_viridis_c(name="Proportion",breaks=seq(0,1,by=0.2)) +
         #theme_bw()+
         theme(axis.text.x=element_text(angle=60, hjust=1)) +
         coord_fixed()+
         facet_grid(type ~ model) +
         labs(x = "", y = "", fill = 'Count') +
         theme( axis.text.x = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                panel.grid = element_blank(),
                legend.position = "bottom",
                legend.key.width = unit(0.5, "in"),
                legend.title = element_text(vjust=0.8)) 
  )
  dev.off()
  
  

  
comparison_data%>%filter(model=="GDM",type=="larch")%>%group_by(n_missing)%>%
  summarise(larch_missing=sum(larch_missing),
            oak_missing=sum(oak_missing),
            sitka_spruce_missing=sum(sitka_spruce_missing),
            sycamore_missing=sum(sycamore_missing))
