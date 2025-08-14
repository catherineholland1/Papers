##########################.
## CONFIGURATIN PLOT ##
##########################.

library(dplyr)
library(tidyr)
library(ggplot2)
library(zCompositions)
library(ggh4x)

##################################.

setwd("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD R Code/Forensic Glass/Configurations")

##################################.

## DATA ####

glass <- read.table("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/Documents/PhD R Code/Forensic Glass/Data/database_190310.txt", header = T)

glass$Type <- as.factor(sapply(glass$Name, function(x) substr(x, 1, 1)))


##################################.

## PLOT CONFIGURATIONS ####

pdf('configuration_plot.pdf', height = 3.5)
zPatterns(glass[,4:11], label = 0,
          axis.labels = c("", ""),
          bar.ordered = c(FALSE, FALSE),
          bar.colors = c("#66C2A5", "#FC8D62"),
          cell.colors = c("#8DA0CB", "grey94"),
          cell.labels = c("Absent", "Present"),
          grid.color = "black", grid.lty = "solid",
          legend = TRUE,
          bar.labels = TRUE, 
          show.means = FALSE,
          type.means = "cgm",
          cex.axis = 0.9) 
dev.off()

##################################.

## CONFIGURATIONS ####

config_data <- glass %>%
  group_by(Item) %>%
  mutate(config = case_when(all(Fe == 0 & K != 0) ~ 2,
                            all(Fe != 0 & K == 0) ~ 3,
                            all(Fe == 0 & K == 0) ~ 4,
                            TRUE ~ 1)) %>%
  ungroup() %>%
  mutate(across(Na:Fe, ~ sqrt(.x/O))) %>% 
  rename(Piece = piece)


### PLOT ####

pdf('boxplot_configuration_plot.pdf', height = 6)
config_data %>%
  pivot_longer(cols = c(Na:Fe), names_to = 'element', values_to = 'value') %>%
  mutate(element = factor(element, levels = c('Na', 'Mg', 'Alu', 'Si', 'K', 'Ca', 'Fe')),
         Type = case_when(Type == 'b' ~ 'bulb',
                          Type == 'c' ~ 'car window',
                          Type == 'h' ~ 'headlamp',
                          Type == 'p' ~ 'container',
                          Type == 'w' ~ 'building window'),
         Type = factor(Type, levels = c('bulb', 'car window', 'headlamp', 'container', 'building window')),
         config = case_when(config == '1' ~ 'configuration 1',
                            config == '2' ~ 'configuration 2',
                            config == '3' ~ 'configuration 3',
                            config == '4' ~ 'configuration 4'),
         config = factor(config, levels = c('configuration 1', 'configuration 2', 'configuration 3', 'configuration 4'))) %>%
  ggplot(., aes(x = factor(Type), y = value, fill = factor(Type))) +
  geom_boxplot(outliers = F) +
  facet_grid2(config~element, scales = 'free', independent = 'y') +
  labs(x = '', y = '') +
  scale_x_discrete(labels = c('','','','','')) +
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.1)  # Round to 1 decimal place
  ) +
  scale_fill_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854"), 
                    name = 'Glass Type', 
                    labels = c('bulb', 'car window', 'headlamp', 'container', 'building window')) +
  theme(#axis.text.x=element_text(angle=60, hjust=1),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 8),
        legend.position = 'bottom',
        legend.margin = margin(t = -15))
dev.off()


pdf('boxplot_configuration_2_plot.pdf', height = 3)
config_data %>%
  pivot_longer(cols = c(Na:Fe), names_to = 'element', values_to = 'value') %>%
  mutate(element = factor(element, levels = c('Na', 'Mg', 'Alu', 'Si', 'K', 'Ca', 'Fe')),
         Type = case_when(Type == 'b' ~ 'bulb',
                          Type == 'c' ~ 'car window',
                          Type == 'h' ~ 'headlamp',
                          Type == 'p' ~ 'container',
                          Type == 'w' ~ 'building window'),
         Type = factor(Type, levels = c('bulb', 'car window', 'headlamp', 'container', 'building window')),
         config = case_when(config == '1' ~ 'configuration 1',
                            config == '2' ~ 'configuration 2',
                            config == '3' ~ 'configuration 3',
                            config == '4' ~ 'configuration 4'),
         config = factor(config, levels = c('configuration 1', 'configuration 2', 'configuration 3', 'configuration 4'))) %>%
  filter(config == 'configuration 2') %>%
  
  ggplot(., aes(x = factor(Type), y = value, fill = factor(Type))) +
  geom_boxplot(outliers = F) +
  facet_grid2(~element, scales = 'free', independent = 'y') +
  labs(x = '', y = '') +
  scale_x_discrete(labels = c('','','','','')) +
  scale_y_continuous(
    labels = scales::number_format(accuracy = 0.1)  # Round to 1 decimal place
  ) +
  scale_fill_manual(values = c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854"), 
                    name = 'Glass Type', 
                    labels = c('bulb', 'car window', 'headlamp', 'container', 'building window')) +
  theme(#axis.text.x=element_text(angle=60, hjust=1),
    axis.ticks.x = element_blank(),
    axis.text = element_text(size = 8),
    legend.position = 'bottom',
    legend.margin = margin(t = -15))
dev.off()

##################################.