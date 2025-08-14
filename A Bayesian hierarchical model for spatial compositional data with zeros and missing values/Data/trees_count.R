########################.
## TREE DATA - COUNTS ##
########################.

## PACKAGES ####

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

#######################.

setwd("Data/")

#######################.

## DATA ####

data <- read_xlsx("simulated_trees.xlsx")

N = 1000

data_counts <- data %>%
  mutate(total = rowSums(select(., ends_with("_pct")))) %>%
  mutate(across(ends_with("_pct"), ~ floor(N / total * .))) %>%
  mutate(total = rowSums(select(., ends_with("_pct")))) %>%
  rename(ash = "Ash_pct",
         beech = "Beech_pct",
         larch = "Larch_pct",
         oak = "Oak_pct",
         scots_pine = "ScPine_pct",
         silver_birch = "SBirch_pct",
         sitka_spruce = "SitSpr_pct",
         sweet_chestnut = "SwChes_pct",
         sycamore = "Syca_pct",
         shadow = "Shadow_pct")

data_proportions <- data %>%
  mutate(total = rowSums(select(., ends_with("_pct")))) %>%
  mutate(total = rowSums(select(., ends_with("_pct")))) %>%
  rename(ash = "Ash_pct",
         beech = "Beech_pct",
         larch = "Larch_pct",
         oak = "Oak_pct",
         scots_pine = "ScPine_pct",
         silver_birch = "SBirch_pct",
         sitka_spruce = "SitSpr_pct",
         sweet_chestnut = "SwChes_pct",
         sycamore = "Syca_pct",
         shadow = "Shadow_pct")


#######################.

## HEATMAP ####

pdf("output/heatmap_originals.pdf", height = 4)
plot(ggplot(data_proportions%>%
              pivot_longer(cols =larch:sycamore,names_to="type",values_to = "value")%>%
              mutate(type=case_match(type,"larch"~"Larch",
                                     "oak"~"Oak","scots_pine"~"Scots pine",
                                     "shadow"~"Shadow","silver_birch"~"Silver birch",
                                     "sitka_spruce"~"Sitka spruce","sweet_chestnut"~"Sweet chestnut",
                                     "sycamore"~"Sycamore")),
            aes(x=X_coord,y=Y_coord,fill=value/total,colour=value/total))+
       geom_tile()+
       facet_wrap(~type,nrow=2)+
       scale_fill_viridis_c(name="Proportion",breaks=seq(0,1,by=0.2)) +
       scale_colour_viridis_c(name="Proportion",breaks=seq(0,1,by=0.2)) +
       theme(axis.text.x=element_text(angle=60, hjust=1)) +
       coord_fixed()+
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
              legend.title = element_text(vjust=0.8)))
dev.off()


#######################.

## SAVE DATA ####

saveRDS(data_counts, "tree_counts_1000.rds")

saveRDS(data_proportions, "tree_proportions.rds")

#######################.