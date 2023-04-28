##Scute TNW Analyses##

library(ggplot2)
library(tidyr)
library(dplyr)
library(data.table)
library(reshape2)
library(lubridate)
library(readxl)
library(grid)
library(gridExtra)
library(ggrepel)
library(lme4)
library(SIBER)
library(viridis)
library(car)



skin<-read.csv("LK_SIA_Data.csv") #load skin data
colnames(skin)[1] <- gsub('^...','',colnames(skin)[1])

summary(skin) #check skin import

skin_factors<-c("TurtleID",
                "LifeStage") #select appropriate vars to convert into factors for skin data
skin[skin_factors]<-lapply(skin[skin_factors], as.factor) #convert selected vars to factors


skin_nums<-c("d15N",
             "SCL",
             "d13C")
skin[skin_nums]<-lapply(skin[skin_nums], as.numeric)

summary(skin) #check conversion


#Overall
model_all_C <- aov(d13C ~ TurtleID, data=skin) #same slope but different intercepts for all turtles
model_all_N <- aov(d15N ~ TurtleID, data=skin)


# Use the mean sum of squares to calculate BIV
BIV_all_C <- anova(model_all_C)["TurtleID", "Mean Sq"] # between individual component
BIV_all_C


BIV_all_N <- anova(model_all_N)["TurtleID", "Mean Sq"] # between individual component
BIV_all_N

##################

#Large turtles
skin_Large <- skin %>%
  filter(LifeStage=="Large")

model_Large_C <- aov(d13C ~ TurtleID, data=skin_Large)
model_Large_N <- aov(d15N ~ TurtleID, data=skin_Large)


# Use the mean sum of squares to calculate WIC and BIV
BIV_Large_C <- anova(model_Large_C)["TurtleID", "Mean Sq"] # between individual component

BIV_Large_C


BIV_Large_N <- anova(model_Large_N)["TurtleID", "Mean Sq"] # between individual component

BIV_Large_N


#Small turtles
skin_Small <- skin %>%
  filter(LifeStage=="Small")

model_Small_C <- aov(d13C ~ TurtleID, data=skin_Small)
model_Small_N <- aov(d15N ~ TurtleID, data=skin_Small)


# Use the mean sum of squares to calculate WIC and BIV
BIV_Small_C <- anova(model_Small_C)["TurtleID", "Mean Sq"] # between individual component

BIV_Small_C


BIV_Small_N <- anova(model_Small_N)["TurtleID", "Mean Sq"] # between individual component

BIV_Small_N



