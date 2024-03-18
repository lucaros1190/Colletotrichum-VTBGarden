
# Data analysis radial growth Colletotrichum Isolates
# Created by Luca Rossini on 14 November 2023
# Last update 30 November 2023
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

levels(isolate)
levels(substrate)
levels(dish)
levels(direction)

  # Create sub-datasets for each isolate for a deeper analysis

dataset_C1 <- data_radius[data_radius$Isolate == "C1", ][-1]
dataset_C10 <- data_radius[data_radius$Isolate == "C10", ][-1]
dataset_C19 <- data_radius[data_radius$Isolate == "C19", ][-1]
dataset_C2 <- data_radius[data_radius$Isolate == "C2", ][-1]
dataset_C3 <- data_radius[data_radius$Isolate == "C3", ][-1]
dataset_C4 <- data_radius[data_radius$Isolate == "C4", ][-1]
dataset_C5 <- data_radius[data_radius$Isolate == "C5", ][-1]
dataset_C7 <- data_radius[data_radius$Isolate == "C7", ][-1]




# LM - Whole dataset: 'weeklyRadius' variable

    # Check if the dataset needs transformation

library(lattice)
qqmath(weeklyRadius)

    # The answer is yes, so data should be transformed!

library(bestNormalize)
weeklyRadius_Trans <- bestNormalize(weeklyRadius) 
weeklyRadius_Trans <- weeklyRadius_Trans$x.t

library(lme4)

GenLin_weeklyRadius <- lmer(weeklyRadius_Trans ~ isolate + substrate + 
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

marginal_weeklyRadius = emmeans(GenLin_weeklyRadius, ~ substrate + isolate)
pairs(marginal_weeklyRadius, adjust="bonferroni")

    # Letters of significance:

library(multcomp)

lettere_weeklyRadius <- cld(marginal_weeklyRadius, alpha=0.05, 
                            Letters=letters, adjust="bonferroni")
lettere_weeklyRadius

    # Pairwise comparison - Isolate:

marginal_weeklyRadius_isolate = emmeans(GenLin_weeklyRadius, ~ isolate)
pairs(marginal_weeklyRadius_isolate, adjust="bonferroni")

    # Letters of significance:

lettere_weeklyRadius_isolate <- cld(marginal_weeklyRadius_isolate, alpha=0.05, 
                                    Letters=letters, adjust="bonferroni")
lettere_weeklyRadius_isolate

    # Pairwise comparison - Substrate:

marginal_weeklyRadius_substrate = emmeans(GenLin_weeklyRadius, ~ substrate)
pairs(marginal_weeklyRadius_substrate, adjust="bonferroni")

    # Letters of significance:

lettere_weeklyRadius_substrate <- cld(marginal_weeklyRadius_substrate, alpha=0.05, 
                                      Letters=letters, adjust="bonferroni")
lettere_weeklyRadius_substrate




# LM - Only C1 to check substrate

    # Check if the dataset needs transformation

qqmath(dataset_C1$MeasureRadius_7D)

    # The answer is yes, so data should be transformed!

weeklyRadius_Trans_C1 <- bestNormalize(dataset_C1$MeasureRadius_7D) 
weeklyRadius_Trans_C1 <- weeklyRadius_Trans_C1$x.t

GenLin_C1 <- lmer(weeklyRadius_Trans_C1 ~ dataset_C1$Substrate + 
                    (1 | dataset_C1$Dish), data=dataset_C1)

summary(GenLin_C1)

    # Check the dispersion and the model reliability

testDispersion(GenLin_C1)
Output_weeklyRadius_C1 <- simulateResiduals(fittedModel = GenLin_C1, plot = T)

qqmath(GenLin_C1)

    # Pairwise comparison - Isolate + Substrate:

marginal_C1 = emmeans(GenLin_C1, ~ "Substrate")
pairs(marginal_C1, adjust="bonferroni")

    # Letters of significance:

lettere_C1 <- cld(marginal_C1, alpha=0.05, 
                            Letters=letters, adjust="bonferroni")
lettere_C1




# LM - Only C10 to check substrate

    # Check if the dataset needs transformation

qqmath(dataset_C10$MeasureRadius_7D)

    # The answer is yes, so data should be transformed!

weeklyRadius_Trans_C10 <- bestNormalize(dataset_C10$MeasureRadius_7D) 
weeklyRadius_Trans_C10 <- weeklyRadius_Trans_C10$x.t

    # Here we are going to use the lm with no random effects because their
    # effect is too small and useless

GenLin_C10 <- lm(weeklyRadius_Trans_C10 ~ dataset_C10$Substrate, 
                 data=dataset_C10)

summary(GenLin_C10)

    # Check the dispersion and the model reliability

testDispersion(GenLin_C10)
Output_weeklyRadius_C10 <- simulateResiduals(fittedModel = GenLin_C10, plot = T)

    # Pairwise comparison - Substrate:

marginal_C10 = emmeans(GenLin_C10, ~ "Substrate")
pairs(marginal_C10, adjust="bonferroni")

    # Letters of significance:

lettere_C10 <- cld(marginal_C10, alpha=0.05, 
                   Letters=letters, adjust="bonferroni")
lettere_C10


# LM - Only C19 to check substrate

    # Check if the dataset needs transformation

qqmath(dataset_C19$MeasureRadius_7D)

    # The answer is yes, so data should be transformed!

weeklyRadius_Trans_C19 <- bestNormalize(dataset_C19$MeasureRadius_7D) 
weeklyRadius_Trans_C19 <- weeklyRadius_Trans_C19$x.t

    # Here we are going to use the lm with no random effects because their effect
    # is too small and useless

GenLin_C19 <- lm(weeklyRadius_Trans_C19 ~ dataset_C19$Substrate, data=dataset_C19)

summary(GenLin_C19)

    # Check the dispersion and the model reliability

testDispersion(GenLin_C19)
Output_weeklyRadius_C19 <- simulateResiduals(fittedModel = GenLin_C19, plot = T)

    # Pairwise comparison - Substrate:

marginal_C19 = emmeans(GenLin_C19, ~ "Substrate")
pairs(marginal_C19, adjust="bonferroni")

    # Letters of significance:

lettere_C19 <- cld(marginal_C19, alpha=0.05, 
                  Letters=letters, adjust="bonferroni")
lettere_C19




# LM - Only C2 to check substrate

    # Check if the dataset needs transformation

qqmath(dataset_C2$MeasureRadius_7D)

    # The answer is yes, so data should be transformed!

weeklyRadius_Trans_C2 <- bestNormalize(dataset_C2$MeasureRadius_7D) 
weeklyRadius_Trans_C2 <- weeklyRadius_Trans_C2$x.t

GenLin_C2 <- lmer(weeklyRadius_Trans_C2 ~ dataset_C2$Substrate + 
                    (1 | dataset_C2$Dish), data=dataset_C2)

summary(GenLin_C2)

    # Check the dispersion and the model reliability

testDispersion(GenLin_C2)
Output_weeklyRadius_C2 <- simulateResiduals(fittedModel = GenLin_C2, plot = T)

qqmath(GenLin_C2)

    # Pairwise comparison - Substrate:

marginal_C2 = emmeans(GenLin_C2, ~ "Substrate")
pairs(marginal_C2, adjust="bonferroni")

    # Letters of significance:

lettere_C2 <- cld(marginal_C2, alpha=0.05, 
                   Letters=letters, adjust="bonferroni")
lettere_C2




# LM - Only C3 to check substrate

    # Check if the dataset needs transformation

qqmath(dataset_C3$MeasureRadius_7D)

    # The answer is yes, so data should be transformed!

weeklyRadius_Trans_C3 <- bestNormalize(dataset_C3$MeasureRadius_7D) 
weeklyRadius_Trans_C3 <- weeklyRadius_Trans_C3$x.t

GenLin_C3 <- lmer(weeklyRadius_Trans_C3 ~ dataset_C3$Substrate + 
                    (1 | dataset_C3$Dish), data=dataset_C3)

summary(GenLin_C3)

    # Check the dispersion and the model reliability

testDispersion(GenLin_C3)
Output_weeklyRadius_C3 <- simulateResiduals(fittedModel = GenLin_C3, plot = T)

qqmath(GenLin_C3)

    # Pairwise comparison - Substrate:

marginal_C3 = emmeans(GenLin_C3, ~ "Substrate")
pairs(marginal_C3, adjust="bonferroni")

    # Letters of significance:

lettere_C3 <- cld(marginal_C3, alpha=0.05, 
                  Letters=letters, adjust="bonferroni")
lettere_C3




# LM - Only C4 to check substrate

    # Check if the dataset needs transformation

qqmath(dataset_C4$MeasureRadius_7D)

    # The answer is yes, so data should be transformed!

weeklyRadius_Trans_C4 <- bestNormalize(dataset_C4$MeasureRadius_7D) 
weeklyRadius_Trans_C4 <- weeklyRadius_Trans_C4$x.t

GenLin_C4 <- lmer(weeklyRadius_Trans_C4 ~ dataset_C4$Substrate + 
                    (1 | dataset_C4$Dish), data=dataset_C4)

summary(GenLin_C4)

    # Check the dispersion and the model reliability

testDispersion(GenLin_C4)
Output_weeklyRadius_C4 <- simulateResiduals(fittedModel = GenLin_C4, plot = T)

qqmath(GenLin_C4)

    # Pairwise comparison - Isolate + Substrate:

marginal_C4 = emmeans(GenLin_C4, ~ "Substrate")
pairs(marginal_C4, adjust="bonferroni")

    # Letters of significance:

lettere_C4 <- cld(marginal_C4, alpha=0.05, 
                  Letters=letters, adjust="bonferroni")
lettere_C4




# LM - Only C5 to check substrate

    # Check if the dataset needs transformation

qqmath(dataset_C5$MeasureRadius_7D)

    # The answer is yes, so data should be transformed!

weeklyRadius_Trans_C5 <- bestNormalize(dataset_C5$MeasureRadius_7D) 
weeklyRadius_Trans_C5 <- weeklyRadius_Trans_C5$x.t

GenLin_C5 <- lm(weeklyRadius_Trans_C5 ~ dataset_C5$Substrate, data=dataset_C5)

summary(GenLin_C5)

    # Pairwise comparison - Substrate:

marginal_C5 = emmeans(GenLin_C5, ~ "Substrate")
pairs(marginal_C5, adjust="bonferroni")

    # Letters of significance:

lettere_C5 <- cld(marginal_C5, alpha=0.05, 
                  Letters=letters, adjust="bonferroni")
lettere_C5




# LM - Only C7 to check substrate

    # Check if the dataset needs transformation

qqmath(dataset_C7$MeasureRadius_7D)

    # The answer is yes, so data should be transformed!

weeklyRadius_Trans_C7 <- bestNormalize(dataset_C7$MeasureRadius_7D) 
weeklyRadius_Trans_C7 <- weeklyRadius_Trans_C7$x.t

GenLin_C7 <- lmer(weeklyRadius_Trans_C7 ~ dataset_C7$Substrate + 
                    (1 | dataset_C7$Dish), data=dataset_C7)

summary(GenLin_C7)

    # Check the dispersion and the model reliability

testDispersion(GenLin_C7)
Output_weeklyRadius_C7 <- simulateResiduals(fittedModel = GenLin_C7, plot = T)

qqmath(GenLin_C7)

    # Pairwise comparison - Substrate:

marginal_C7 = emmeans(GenLin_C7, ~ "Substrate")
pairs(marginal_C7, adjust="bonferroni")

    # Letters of significance:

lettere_C7 <- cld(marginal_C7, alpha=0.05, 
                  Letters=letters, adjust="bonferroni")
lettere_C7




# Overall plot of the results - Only the general dataset

library(ggplot2)

boxPlot_Complete <- ggplot(data_radius, aes(x=isolate, y=weeklyRadius,
                                           fill=substrate)) + 
                          geom_boxplot(width=0.7) + 
                          xlab("Isolate") + 
                          ylab("Radius lenght (mm)") + 
                          ggtitle("Radius lenght over the isolates") +
                          theme(plot.title = element_text(hjust=0.5), 
                                text = element_text(size=21)) + 
                          scale_fill_manual(drop = FALSE, name = "Substrate", 
                                labels = c("OA","PDA", "SNA"), 
                                values= alpha(c("green", "blue", 
                                                "orange", "purple"), 0.5)) + 
                          scale_x_discrete(labels = c("C1","C10", "C19", "C2", 
                                                      "C3", "C4", "C5", "C7"))

boxPlot_Complete


boxPlot_isolates <- ggplot(data_radius, aes(x=isolate, y=weeklyRadius,
                                           fill=isolate)) + 
                          geom_boxplot(width=0.7) + 
                          xlab("Isolate") + 
                          ylab("Radius lenght (mm)") + 
                          ggtitle("Radius lenght over the isolates") +
                          theme(plot.title = element_text(hjust=0.5), 
                                text = element_text(size=21)) + 
                          theme(legend.position = "none") +
                          scale_x_discrete(labels = c("C1","C10", "C19", "C2", 
                                                      "C3", "C4", "C5", "C7"))

boxPlot_isolates


boxPlot_substrates <- ggplot(data_radius, aes(x=substrate, y=weeklyRadius,
                                            fill=substrate)) + 
                          geom_boxplot(width=0.7) + 
                          xlab("Substrate") + 
                          ylab("Radius lenght (mm)") + 
                          ggtitle("Radius lenght over the substrates") +
                          theme(plot.title = element_text(hjust=0.5), 
                                text = element_text(size=21)) + 
                          theme(legend.position = "none") +
                          scale_x_discrete(labels = c("OA","PDA", "SNA"))

boxPlot_substrates




# Detailed plots - Isolate growth per substrate


    # Isolate C1

boxPlotSub_C1 <- ggplot(dataset_C1, aes(x=Substrate, 
                                        y=MeasureRadius_7D, 
                                        fill=Substrate)) + 
                          geom_boxplot(width=0.7) + 
                          xlab("Substrate") + 
                          ylab("Radius lenght (mm)") + 
                          ggtitle("Isolate C1") +
                          theme(plot.title = element_text(hjust=0.5), 
                                      text = element_text(size=21)) + 
                          theme(legend.position = "none") +
                          scale_x_discrete(labels = c("OA","PDA", "SNA"))

boxPlotSub_C1


    # Isolate C10

boxPlotSub_C10 <- ggplot(dataset_C10, aes(x=Substrate, 
                                        y=MeasureRadius_7D, 
                                        fill=Substrate)) + 
                          geom_boxplot(width=0.7) + 
                          xlab("Substrate") + 
                          ylab("Radius lenght (mm)") + 
                          ggtitle("Isolate C10") +
                          theme(plot.title = element_text(hjust=0.5), 
                                      text = element_text(size=21)) +
                          theme(legend.position = "none") +
                          scale_x_discrete(labels = c("OA","PDA", "SNA"))

boxPlotSub_C10


# Isolate C19

boxPlotSub_C19 <- ggplot(dataset_C19, aes(x=Substrate, 
                                          y=MeasureRadius_7D, 
                                          fill=Substrate)) + 
                          geom_boxplot(width=0.7) + 
                          xlab("Substrate") + 
                          ylab("Radius lenght (mm)") + 
                          ggtitle("Isolate C19") +
                          theme(plot.title = element_text(hjust=0.5), 
                                      text = element_text(size=21)) + 
                          theme(legend.position = "none") +
                          scale_x_discrete(labels = c("OA","PDA", "SNA"))

boxPlotSub_C19


# Isolate C2

boxPlotSub_C2 <- ggplot(dataset_C2, aes(x=Substrate, 
                                        y=MeasureRadius_7D, 
                                        fill=Substrate)) + 
                          geom_boxplot(width=0.7) + 
                          xlab("Substrate") + 
                          ylab("Radius lenght (mm)") + 
                          ggtitle("Isolate C2") +
                          theme(plot.title = element_text(hjust=0.5), 
                                      text = element_text(size=21)) +
                          theme(legend.position = "none") +
                          scale_x_discrete(labels = c("OA","PDA", "SNA"))

boxPlotSub_C2


# Isolate C3

boxPlotSub_C3 <- ggplot(dataset_C3, aes(x=Substrate, 
                                          y=MeasureRadius_7D, 
                                          fill=Substrate)) + 
                          geom_boxplot(width=0.7) + 
                          xlab("Substrate") + 
                          ylab("Radius lenght (mm)") + 
                          ggtitle("Isolate C3") +
                          theme(plot.title = element_text(hjust=0.5), 
                                      text = element_text(size=21)) + 
                          theme(legend.position = "none") +
                          scale_x_discrete(labels = c("OA","PDA", "SNA"))

boxPlotSub_C3


# Isolate C4

boxPlotSub_C4 <- ggplot(dataset_C4, aes(x=Substrate, 
                                        y=MeasureRadius_7D, 
                                        fill=Substrate)) + 
                          geom_boxplot(width=0.7) + 
                          xlab("Substrate") + 
                          ylab("Radius lenght (mm)") + 
                          ggtitle("Isolate C4") +
                          theme(plot.title = element_text(hjust=0.5), 
                                      text = element_text(size=21)) +
                          theme(legend.position = "none") +
                          scale_x_discrete(labels = c("OA","PDA", "SNA"))

boxPlotSub_C4


# Isolate C5

boxPlotSub_C5 <- ggplot(dataset_C5, aes(x=Substrate, 
                                        y=MeasureRadius_7D, 
                                        fill=Substrate)) + 
                          geom_boxplot(width=0.7) + 
                          xlab("Substrate") + 
                          ylab("Radius lenght (mm)") + 
                          ggtitle("Isolate C5") +
                          theme(plot.title = element_text(hjust=0.5), 
                                      text = element_text(size=21)) + 
                          theme(legend.position = "none") +
                          scale_x_discrete(labels = c("OA","PDA", "SNA"))

boxPlotSub_C5


# Isolate C7

boxPlotSub_C7 <- ggplot(dataset_C7, aes(x=Substrate, 
                                          y=MeasureRadius_7D, 
                                          fill=Substrate)) + 
                          geom_boxplot(width=0.7) + 
                          xlab("Substrate") + 
                          ylab("Radius lenght (mm)") + 
                          ggtitle("Isolate C7") +
                          theme(plot.title = element_text(hjust=0.5), 
                                      text = element_text(size=21)) +
                          theme(legend.position = "none") +
                          scale_x_discrete(labels = c("OA","PDA", "SNA"))

boxPlotSub_C7


# Making a grid resuming the plots

library(ggpubr)

gridPlotFirst <- ggarrange(boxPlot_isolates, boxPlot_substrates,
                           ncol = 1, nrow = 2)

gridPlotFirst
  
  
gridPlotSecond <- ggarrange(boxPlotSub_C1, boxPlotSub_C10, boxPlotSub_C19, 
                            boxPlotSub_C2, boxPlotSub_C3, boxPlotSub_C4, 
                            boxPlotSub_C5, boxPlotSub_C7, ncol = 2, nrow = 4)

grid_plotSecond
