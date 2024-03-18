
# Data analysis radial growth Colletotrichum species
# Created by Luca Rossini on 6 February 2024
# Last update 7 February 2024
# E-mail: luca.rossini@unitus.it


# Acquisition of the data - File 'DatasetRaggi.csv'

filename <- file.choose()

data_radius = read.csv2(filename, header=T, sep=";", dec=".", 
                        na.string="NA")

head(data_radius)

species <- as.factor(data_radius$Species)
isolate <- as.factor(data_radius$Isolate)
substrate <- as.factor(data_radius$Substrate)
dish <- as.factor(data_radius$Dish)
direction <- as.factor(data_radius$Radius)
weeklyRadius <- as.numeric(data_radius$MeasureRadius_7D)
dailyRadius <- as.numeric(data_radius$MeasureRadius_1D)

  # Check the levels

levels(species)
levels(isolate)
levels(substrate)
levels(dish)
levels(direction)

  # Create sub-datasets for each species for a deeper analysis

dataset_karstii <- data_radius[data_radius$Species == "C_karstii", ]

    # Assign factors and numeric values

isolate_karstii <- as.factor(dataset_karstii$Isolate)
substrate_karstii <- as.factor(dataset_karstii$Substrate)
dish_karstii <- as.factor(dataset_karstii$Dish)
direction_karstii <- as.factor(dataset_karstii$Radius)
weeklyRadius_karstii <- as.numeric(dataset_karstii$MeasureRadius_7D)


dataset_cattleyicola <- data_radius[data_radius$Species == "C_cattleyicola", ]

    # Assign factors and numeric values

isolate_cattleyicola <- as.factor(dataset_cattleyicola$Isolate)
substrate_cattleyicola <- as.factor(dataset_cattleyicola$Substrate)
dish_cattleyicola <- as.factor(dataset_cattleyicola$Dish)
direction_cattleyicola <- as.factor(dataset_cattleyicola$Radius)
weeklyRadius_cattleyicola <- as.numeric(dataset_cattleyicola$MeasureRadius_7D)


dataset_boninense <- data_radius[data_radius$Species == "C_boninense", ]

    # Assign factors and numeric values

isolate_boninense <- as.factor(dataset_boninense$Isolate)
substrate_boninense <- as.factor(dataset_boninense$Substrate)
dish_boninense <- as.factor(dataset_boninense$Dish)
direction_boninense <- as.factor(dataset_boninense$Radius)
weeklyRadius_boninense <- as.numeric(dataset_boninense$MeasureRadius_7D)




# LM - Whole dataset: 'weeklyRadius' variable

    # Check if the dataset needs transformation

library(lattice)
qqmath(weeklyRadius)

    # The answer is yes, so data should be transformed!

library(bestNormalize)
weeklyRadius_Trans <- bestNormalize(weeklyRadius) 
weeklyRadius_Trans <- weeklyRadius_Trans$x.t

library(lme4)

GenLin_weeklyRadius <- lmer(weeklyRadius_Trans ~ species + substrate + 
                              (1 | dish) , data=data_radius)

summary(GenLin_weeklyRadius)

    # Check the dispersion and the model reliability

library(DHARMa)

testDispersion(GenLin_weeklyRadius)
Output_weeklyRadius <- simulateResiduals(fittedModel = GenLin_weeklyRadius, 
                                         plot = T)

library(lattice)
qqmath(GenLin_weeklyRadius)

    # Pairwise comparison - Isolate + Substrate:

library(multcompView)
library(emmeans)

marginal_weeklyRadius = emmeans(GenLin_weeklyRadius, ~ substrate + species)
pairs(marginal_weeklyRadius, adjust="bonferroni")

    # Letters of significance:

library(multcomp)

lettere_weeklyRadius <- cld(marginal_weeklyRadius, alpha=0.05, 
                            Letters=letters, adjust="bonferroni")
lettere_weeklyRadius

    # Pairwise comparison - Species:

marginal_weeklyRadius_species = emmeans(GenLin_weeklyRadius, ~ species)
pairs(marginal_weeklyRadius_species, adjust="bonferroni")

    # Letters of significance:

lettere_weeklyRadius_species <- cld(marginal_weeklyRadius_species, alpha=0.05, 
                                    Letters=letters, adjust="bonferroni")
lettere_weeklyRadius_species

    # Pairwise comparison - Substrate:

marginal_weeklyRadius_substrate = emmeans(GenLin_weeklyRadius, ~ substrate)
pairs(marginal_weeklyRadius_substrate, adjust="bonferroni")

    # Letters of significance:

lettere_weeklyRadius_substrate <- cld(marginal_weeklyRadius_substrate, alpha=0.05, 
                                      Letters=letters, adjust="bonferroni")
lettere_weeklyRadius_substrate




# LM - Only C. boninense to check substrate

    # Check if the dataset needs transformation

qqmath(weeklyRadius_boninense)

    # The answer is yes, so data should be transformed!

weeklyRadius_Trans_boninense <- bestNormalize(weeklyRadius_boninense) 
weeklyRadius_Trans_boninense <- weeklyRadius_Trans_boninense$x.t

GenLin_boninense <- lmer(weeklyRadius_Trans_boninense ~ substrate_boninense +
                                                      (1 | isolate_boninense) +
                                                      (1 | dish_boninense), 
                                                      data=dataset_boninense)

summary(GenLin_boninense)

    # Check the dispersion and the model reliability

testDispersion(GenLin_boninense)
Output_weeklyRadius_boninense <- simulateResiduals(fittedModel = 
                                                    GenLin_boninense, plot = T)

qqmath(GenLin_boninense)

    # Pairwise comparison - Substrate:

marginal_boninense = emmeans(GenLin_boninense, ~ substrate_boninense)
pairs(marginal_boninense, adjust="bonferroni")

    # Letters of significance:

lettere_boninense <- cld(marginal_boninense, alpha=0.05, 
                                  Letters=letters, adjust="bonferroni")
lettere_boninense




# LM - Only karstii to check substrate

    # Check if the dataset needs transformation

qqmath(weeklyRadius_karstii)

    # The answer is yes, so data should be transformed!

weeklyRadius_Trans_karstii <- bestNormalize(weeklyRadius_karstii) 
weeklyRadius_Trans_karstii <- weeklyRadius_Trans_karstii$x.t

    # Here we are going to use the lm with no random effects because their
    # effect is too small and useless

GenLin_karstii <- lm(weeklyRadius_Trans_karstii ~ substrate_karstii,
                                  data=dataset_karstii)

summary(GenLin_karstii)

    # Check the dispersion and the model reliability

testDispersion(GenLin_karstii)
Output_weeklyRadius_karstii <- simulateResiduals(fittedModel = 
                                                    GenLin_karstii, plot = T)

    # Pairwise comparison - Substrate:

marginal_karstii = emmeans(GenLin_karstii, ~ substrate_karstii)
pairs(marginal_karstii, adjust="bonferroni")

    # Letters of significance:

lettere_karstii <- cld(marginal_karstii, alpha=0.05, 
                               Letters=letters, adjust="bonferroni")
lettere_karstii




# LM - Only C. cattleyicola to check substrate

    # Check if the dataset needs transformation

qqmath(weeklyRadius_cattleyicola)

    # The answer is yes, so data should be transformed!

weeklyRadius_Trans_cattleyicola <- bestNormalize(weeklyRadius_cattleyicola) 
weeklyRadius_Trans_cattleyicola <- weeklyRadius_Trans_cattleyicola$x.t

    # Here we are going to use the lm with no random effects because their effect
    # is too small and useless

GenLin_cattleyicola <- lmer(weeklyRadius_Trans_cattleyicola ~ substrate_cattleyicola +
                                    (1 | isolate_cattleyicola) +
                                    (1 | dish_cattleyicola), 
                                    data = dataset_cattleyicola)

summary(GenLin_cattleyicola)

    # Check the dispersion and the model reliability

testDispersion(GenLin_cattleyicola)
Output_weeklyRadius_cattleyicola <- simulateResiduals(fittedModel = 
                                                    GenLin_cattleyicola,plot = T)

    # Pairwise comparison - Substrate:

marginal_cattleyicola = emmeans(GenLin_cattleyicola, ~ substrate_cattleyicola)
pairs(marginal_cattleyicola, adjust="bonferroni")

    # Letters of significance:

lettere_cattleyicola <- cld(marginal_cattleyicola, alpha=0.05, 
                                    Letters=letters, adjust="bonferroni")
lettere_cattleyicola






# Overall plot of the results - Only the general dataset

library(ggplot2)

boxPlot_Complete <- ggplot(data_radius, aes(x=species, y=weeklyRadius,
                                           fill=substrate)) + 
                          geom_boxplot(width=0.5) + 
                          xlab("Isolate") + 
                          ylab("Radius lenght (mm)") + 
                          ggtitle("Radius lenght over the isolates") +
                          theme(plot.title = element_text(hjust=0.5), 
                                text = element_text(size=21)) + 
                          scale_fill_manual(drop = FALSE, name = "Substrate", 
                                labels = c("OA","PDA", "SNA"), 
                                values= alpha(c("green", "blue", 
                                                "orange"), 0.5))  + 
                          scale_x_discrete(guide = guide_axis(angle = 90), labels = 
                                c(expression(italic("Colletotrichum \n boninense")),
                                  expression(italic("Colletotrichum \n cattleyicola")), 
                                  expression(italic("Colletotrichum \n karstii"))))

boxPlot_Complete


boxPlot_species <- ggplot(data_radius, aes(x=species, y=weeklyRadius,
                                           fill=species)) + 
                          geom_boxplot(width=0.3) + 
                          xlab("Species") + 
                          ylab("Radius lenght (mm)") + 
                          ggtitle("Radius lenght over the Species") +
                          theme(plot.title = element_text(hjust=0.5), 
                                text = element_text(size=21)) + 
                          theme(legend.position = "none") +
                          scale_x_discrete(guide = guide_axis(angle = 90),labels = 
                                c(expression(italic("Colletotrichum \n boninense")),
                                  expression(italic("Colletotrichum \n cattleyicola")), 
                                  expression(italic("Colletotrichum \n karstii"))))

boxPlot_species


boxPlot_substrate <- ggplot(data_radius, aes(x=substrate, y=weeklyRadius,
                                           fill=substrate)) + 
                          geom_boxplot(width=0.3) + 
                          xlab("Substrate") + 
                          ylab("Radius lenght (mm)") + 
                          ggtitle("Radius lenght over the substrates") +
                          theme(plot.title = element_text(hjust=0.5), 
                                text = element_text(size=21)) + 
                          theme(legend.position = "none") +
                          scale_x_discrete(labels = c("OA","PDA", "SNA"))

boxPlot_substrate




# Detailed plots - Isolate growth per substrate


    # C. cattleyicola

boxPlotSub_cattleyicola <- ggplot(dataset_cattleyicola, aes(x=Substrate, 
                                                          y=MeasureRadius_7D, 
                                                          fill=Substrate)) + 
                          geom_boxplot(width=0.3) +
                          ylab("Radius lenght (mm)") + 
                          ggtitle("C. cattleyicola") +
                          theme(plot.title = element_text(hjust=0.5, 
                                                          face = "italic"),
                                      text = element_text(size=21)) + 
                          theme(legend.position = "none") +
                          scale_x_discrete(labels = c("OA","PDA", "SNA"))

boxPlotSub_cattleyicola


    # C. karstii

boxPlotSub_karstii <- ggplot(dataset_karstii, aes(x=Substrate, 
                                                  y=MeasureRadius_7D, 
                                                  fill=Substrate)) + 
                          geom_boxplot(width=0.3) +
                          ylab("Radius lenght (mm)") + 
                          ggtitle("C. karstii") +
                          theme(plot.title = element_text(hjust=0.5, 
                                                          face = "italic"), 
                                      text = element_text(size=21)) +
                          theme(legend.position = "none") +
                          scale_x_discrete(labels = c("OA","PDA", "SNA"))

boxPlotSub_karstii


    # C. boninense

boxPlotSub_boninense <- ggplot(dataset_boninense, aes(x=Substrate, 
                                                      y=MeasureRadius_7D, 
                                                      fill=Substrate)) + 
                          geom_boxplot(width=0.3) + 
                          ylab("Radius lenght (mm)") + 
                          ggtitle("C. boninense") +
                          theme(plot.title = element_text(hjust=0.5, 
                                                          face = "italic"), 
                                      text = element_text(size=21)) + 
                          theme(legend.position = "none") +
                          scale_x_discrete(labels = c("OA","PDA", "SNA"))

boxPlotSub_boninense



# Making a grid resuming the plots

library(ggpubr)

gridPlot <- ggarrange(boxPlotSub_boninense, boxPlotSub_karstii, 
                            boxPlotSub_cattleyicola, ncol = 1, nrow = 3)

gridPlot
