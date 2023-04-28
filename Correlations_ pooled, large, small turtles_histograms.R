

### Analysis associated with Weber publication: "Foraging ecology of Kemp’s ridley (Lepidochelys kempii) turtles in the northeastern Gulf of Mexico: insights from stable isotope analysis"

#load required packages


library("ggpubr")
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
library(forcats)


####Read in data #####
skin<-read.csv("LK_SIA_Data.csv") #load skin data
colnames(skin)[1] <- gsub('^...','',colnames(skin)[1])


summary(skin) #check skin import

skin_factors<-c("TurtleID",
                "LifeStage") #select appropriate variables to convert into factors for skin data
skin[skin_factors]<-lapply(skin[skin_factors], as.factor) #convert selected vars to factors


skin_nums<-c( "SCL",
              "d15N",
              "Weight",
              "BCI",
              "d13C")

skin[skin_nums]<-lapply(skin[skin_nums], as.numeric)
summary(skin)
head(skin)


###testing assumptions of normality to know whether to use pearson or spearman correlations in scatterplots 

# Shapiro-Wilk normality test for SCL
shapiro.test(skin$SCL) # => p = 0.0047    
#not greater than the significance level of 0.05, implying that the distribution ARE significantly different from normal distribution
#Not normally distributed

# Shapiro-Wilk normality test for d13C
shapiro.test(skin$d13C) # => p = 0.24
#NORMALLY DISTRIBUTED

# Shapiro-Wilk normality test for d15N
shapiro.test(skin$d15N) # => p = 0.00000159 
##not greater than the significance level of 0.05, implying that the distribution ARE significantly different from normal distribution
#NOT Normally distributed

shapiro.test(skin$Weight) # => p=0.3791
#Greater than the significance level of 0.05
#NORMALLY DISTRIBUTED

shapiro.test(skin$BCI)
#<0.05, so NOT normally distributed

###### Scatter plots testing for relationships in POOLED data ####

SCL_C <- ggscatter(skin, x = "SCL", y = "d13C", 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "spearman",
                   xlab = "SCL (cm)", ylab = expression(delta^13* C * "  \u2030"))

SCL_N <- ggscatter(skin, x = "SCL", y = "d15N", 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "spearman",
                   xlab = "SCL (cm)", ylab = expression(delta^15* N * "  \u2030"))

Weight_C <- ggscatter(skin, x = "Weight", y = "d13C", 
                      add = "reg.line", conf.int = TRUE, 
                      cor.coef = TRUE, cor.method = "pearson",
                      xlab = "Weight (kg)", ylab = expression(delta^13* C * "  \u2030"))

Weight_N <- ggscatter(skin, x = "Weight", y = "d15N", 
                      add = "reg.line", conf.int = TRUE, 
                      cor.coef = TRUE, cor.method = "spearman",
                      xlab = "Weight (kg)", ylab = expression(delta^15* N * "  \u2030"))

cBCI <- ggscatter(skin, x = "BCI", y = "d13C", 
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = "spearman",
                  xlab = "BCI", ylab = expression(delta^13* C * "  \u2030"))

nBCI <- ggscatter(skin, x = "BCI", y = "d15N", 
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = "spearman",
                  xlab = "BCI", ylab = expression(delta^15* N * "  \u2030"))

#arrange all above plots to be shown together
library(ggpubr)
ggarrange(SCL_C, Weight_C, cBCI, SCL_N, Weight_N, nBCI, 
          labels=c("A", "B", "C", "D", "E", "F"),
          ncol=3, nrow=2)


##### Scatter plots testing for relationships in SMALL TURTLE data ######
############################################################################


small <- skin [which(skin$LifeStage=='Small'),]
summary(small)
head(small)

#check for normality, to know whether to do pearson or spearman correlation for scatterplots

# Shapiro-Wilk normality test for SCL
shapiro.test(small$SCL) # => p = 0.09    
#greater than the significance level of 0.05
#Normally distributed

# Shapiro-Wilk normality test for d13C
shapiro.test(small$d13C) # => p = 0.7001
#Normally distributed

# Shapiro-Wilk normality test for d15N
shapiro.test(small$d15N) # => p = 0.7653
# Normally distributed

shapiro.test(small$Weight) # => p=0.9277
#Greater than the significance level of 0.05
#Normally distributed

shapiro.test(small$BCI) # p = 0.0157
#NOT normally distributed



SCL_C <- ggscatter(small, x = "SCL", y = "d13C", 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "pearson",
                   xlab = "SCL (cm)", ylab = expression(delta^13* C * "  \u2030"))

SCL_N <- ggscatter(small, x = "SCL", y = "d15N", 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "pearson",
                   xlab = "SCL (cm)", ylab = expression(delta^15* N * "  \u2030"))

Weight_C <- ggscatter(small, x = "Weight", y = "d13C", 
                      add = "reg.line", conf.int = TRUE, 
                      cor.coef = TRUE, cor.method = "pearson",
                      xlab = "Weight (kg)", ylab = expression(delta^13* C * "  \u2030"))

Weight_N <- ggscatter(small, x = "Weight", y = "d15N", 
                      add = "reg.line", conf.int = TRUE, 
                      cor.coef = TRUE, cor.method = "pearson",
                      xlab = "Weight (kg)", ylab = expression(delta^15* N * "  \u2030"))

cBCI <- ggscatter(small, x = "BCI", y = "d13C", 
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = "spearman",
                  xlab = "BCI", ylab = expression(delta^13* C * "  \u2030"))

nBCI <- ggscatter(small, x = "BCI", y = "d15N", 
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = "spearman",
                  xlab = "BCI", ylab = expression(delta^15* N * "  \u2030"))

ggarrange(SCL_C, Weight_C, cBCI, SCL_N, Weight_N, nBCI, 
          labels=c("A", "B", "C", "D", "E", "F"),
          ncol=3, nrow=2)

#### Scatter plots testing for relationships in LARGE TURTLE data ######
############################################################################


large <- skin [which(skin$LifeStage=='Large'),]
summary(large)
head(large)


# Shapiro-Wilk normality test for SCL
shapiro.test(large$SCL) # => p = 0.7079   
#greater than the significance level of 0.05
#Normally distributed

# Shapiro-Wilk normality test for d13C
shapiro.test(large$d13C) # => p = 0.1673
#Normally distributed

# Shapiro-Wilk normality test for d15N
shapiro.test(large$d15N) # => p = <0.05
# NOT Normally distributed

shapiro.test(large$Weight) # => p=0.2449
#Greater than the significance level of 0.05
#Normally distributed

shapiro.test(large$BCI) # p = 0.651
#Normally distributed


SCL_C <- ggscatter(large, x = "SCL", y = "d13C", 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "pearson",
                   xlab = "SCL (cm)", ylab = expression(delta^13* C * "  \u2030"))

SCL_N <- ggscatter(large, x = "SCL", y = "d15N", 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "spearman",
                   xlab = "SCL (cm)", ylab = expression(delta^15* N * "  \u2030"))

Weight_C <- ggscatter(large, x = "Weight", y = "d13C", 
                      add = "reg.line", conf.int = TRUE, 
                      cor.coef = TRUE, cor.method = "pearson",
                      xlab = "Weight (kg)", ylab = expression(delta^13* C * "  \u2030"))

Weight_N <- ggscatter(large, x = "Weight", y = "d15N", 
                      add = "reg.line", conf.int = TRUE, 
                      cor.coef = TRUE, cor.method = "spearman",
                      xlab = "Weight (kg)", ylab = expression(delta^15* N * "  \u2030"))

cBCI <- ggscatter(large, x = "BCI", y = "d13C", 
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = "pearson",
                  xlab = "BCI", ylab = expression(delta^13* C * "  \u2030"))

nBCI <- ggscatter(large, x = "BCI", y = "d15N", 
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = "spearman",
                  xlab = "BCI", ylab = expression(delta^15* N * "  \u2030"))


ggarrange(SCL_C, Weight_C, cBCI, SCL_N, Weight_N, nBCI, 
          labels=c("A", "B", "C", "D", "E", "F"),
          ncol=3, nrow=2)


#############   HISTOGRAMS     #############

hist(skin$SCL,
     xlab="Straight Carapace Length (cm)",
     xlim=c(20, 60),
     ylim=c(0, 20),
     breaks=c(20,25,30,35,40,45,50,55,60),
     las=1)


###### Carbon ##########
library(ggplot2)
library(ggplotify)
library(cowplot)


LC<-as.ggplot(~hist(large$d13C,
                    main="Larger Size Class",
                    xlab=expression(delta^13* C * " \u2030"),
                    ylab="Number of individuals",
                    xlim=c(-17,-11),
                    ylim=c(0,26),
                    cex.lab = 1.2,
                    las=1))


SC<-as.ggplot(~hist(small$d13C,
                    main="Smaller Size Class",
                    breaks=c(-17,-16,-15,-14,-13,-12),
                    xlab=expression(delta^13* C * " \u2030"),
                    ylab="Number of individuals",
                    xlim=c(-17,-11),
                    ylim=c(0,7),
                    cex.lab = 1.2,
                    las=1))
########## Nitrogen ########
LN<-as.ggplot(~hist(large$d15N,
                    breaks=c(5,6,7,8,9,10),
                    main=NULL,
                    xlab=expression(delta^15* N * " \u2030"),
                    ylab="Number of individuals",
                    xlim=c(4, 14),
                    cex.lab = 1.2,
                    las=1))

SN<-as.ggplot(~hist(small$d15N,
                    breaks=c(5,6,7,8,9,10),
                    main=NULL,
                    xlab=expression(delta^15* N * " \u2030"),
                    ylab="Number of individuals",
                    xlim=c(4, 14),
                    cex.lab = 1.2,
                    las=1))


FigS1 <- plot_grid(SC, LC, SN, LN, labels = c("A", "B", "C", "D"), ncol=2, nrow=2)
FigS1

ggsave(filename = "FigS1.eps", FigS1, width = 174, height = 220, units = "mm")
#dev.off()

ggsave(filename = "FigS1.jpg", FigS1, width = 174 , height = 220, units = "mm")





