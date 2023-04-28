
### Analysis associated with Weber publication: "Foraging ecology of Kemp’s ridley (Lepidochelys kempii) turtles in the northeastern Gulf of Mexico: insights from stable isotope analysis"

### Boxplots, range, mean, SD for 2 lifestages of C&N data ###

#load required packages


library(tidyverse)
library("DescTools")
library("car")
library("gmodels")
library("ggplot2")
library("qqplotr")
library("dplyr")
library("coin")
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(gridExtra)

initialCR_CN <- read.csv("LK_SIA_Data.csv")
colnames(initialCR_CN)[1] <- gsub('^...','',colnames(initialCR_CN)[1])

weight_BCI_data <- initialCR_CN %>%
  drop_na(Weight) #drop NAs for Weight & BCI (body condition index)

summary(weight_BCI_data)

#means, sd, range of isotopes (C & N), BCI, and weight


#####  pooled data  #######
initialCR_CN %>%
  summarize(count = n(),
            mC = mean(d13C),
            sdC = sd(d13C), 
            rC = range(d13C),
            mN = mean(d15N), 
            sdN = sd(d15N),
            rN = range(d15N))

weight_BCI_data %>% 
  summarize(count = n(),
            mW = mean(Weight), 
            sdW = sd(Weight),
            rW = range(Weight))

weight_BCI_data %>% 
  summarize(count = n(),
            mBCI = mean(BCI), 
            sdBCI = sd(BCI),
            rBCI = range(BCI))

#### metrics by size class (coded here as LifeStage) ####
initialCR_CN %>% group_by(LifeStage) %>%
  summarize(count = n(),
            mC = mean(d13C),
            sdC = sd(d13C), 
            rC = range(d13C),
            mN = mean(d15N), 
            sdN = sd(d15N),
            rN = range(d15N))

weight_BCI_data %>% group_by(LifeStage) %>%
  summarize(count = n(),
            mW = mean(Weight), 
            sdW = sd(Weight),
            rW = range(Weight))

weight_BCI_data %>% group_by(LifeStage) %>%
  summarize(count = n(),
            mBCI = mean(BCI), 
            sdBCI = sd(BCI),
            rBCI = range(BCI))


#visualizing correlations among variables with 'pairs' function, to display histograms (diagonol), pair-wise scatterplots (upper panel), and correlations (lower panel)
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  # text(0.5, 0.5, txt, cex = cex.cor * r)
  text(0.5, 0.5, txt, cex = 3)
}

pairs(select(initialCR_CN, d13C, d15N), 
      diag.panel = panel.hist,
      upper.panel = panel.smooth,
      lower.panel = panel.cor)

pairs(select(initialCR_CN, d13C, d15N))

############ Boxplots code ###########

CR_CN<-read.csv("LK_SIA_Data.csv")

#boxplots by size class
# first Carbon
par(mfrow = c(1,2))
boxplot(d13C ~ LifeStage, data = CR_CN, xlab = "Size Group", ylab = expression(delta^13* C * "  \u2030"))

# second Nitrogen
boxplot(d15N ~ LifeStage, data = CR_CN, xlab = "Size Group", ylab = expression(delta^15* N * "  \u2030"))

#making them prettier:

CR_CN$LifeStage <- factor(CR_CN$LifeStage, c("Small", "Large"))


bpC <- ggplot(CR_CN, aes(x = LifeStage, y = d13C, fill = LifeStage)) +
  stat_boxplot(geom ="errorbar", width = 0.5) +
  geom_boxplot(fill = "grey", outlier.size = 2) +
  ###line below adds the mean, shown as a circle with a cross through it
  stat_summary(fun=mean, geom="point", shape=10, size=5, color="black") +
  theme_bw() + theme(legend.position="none", axis.line = element_line(color='black'),
                     plot.background = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.grid.major = element_blank(),
                     axis.text.x = element_text(color="black", size=12),
                     axis.text.y = element_text(size=12),
                     axis.title.x = element_text(size=14),
                     axis.title.y = element_text(size=14)) +
  xlab (NULL) + ylab(expression(delta^13* C * "  \u2030")) +
  annotate("text", x = -Inf, y = Inf, hjust = 0, vjust = 1,
           label = "MDiff = -0.86, CI: (-1.41, -0.29)", size = 5)

bpN <- ggplot(CR_CN, aes(x = LifeStage, y = d15N, fill = LifeStage)) +
  stat_boxplot(geom ="errorbar", width = 0.5) +
  geom_boxplot(fill = "grey", outlier.size = 2) + 
  stat_summary(fun=mean, geom="point", shape=10, size=5, color="black") + 
  theme_bw() + theme(legend.position="none", axis.line = element_line(color='black'),
                          plot.background = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.grid.major = element_blank(),
                     axis.text.x = element_text(color="black", size=12),
                     axis.text.y = element_text(size=12),
                     axis.title.x = element_text(size=14),
                     axis.title.y = element_text(size=14)) +
  xlab (NULL) + ylab(expression(delta^15* N * "  \u2030")) +
  annotate("text", x = -Inf, y = Inf, hjust = 0, vjust = 1,
           label = "MDiff = 0.73, CI: (0.19, 1.30)", size = 5)
################################################################################
bpC <- ggplot(CR_CN, aes(x = LifeStage, y = d13C, fill = LifeStage)) +
  stat_boxplot(geom ="errorbar", width = 0.5) +
  geom_boxplot(fill = "grey", outlier.size = 2) +
  ###line below adds the mean, shown as a circle with a cross through it
  stat_summary(fun=mean, geom="point", shape=10, size=5, color="black") +
  theme_bw() + theme(legend.position="none", axis.line = element_line(color='black'),
                     plot.background = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.grid.major = element_blank(),
                     axis.text.x = element_text(color="black", size=12),
                     axis.text.y = element_text(size=12),
                     axis.title.x = element_text(size=14),
                     axis.title.y = element_text(size=14)) +
  xlab (NULL) + ylab(expression(delta^13* C * "  \u2030"))

Clab <- bpC + annotate("text", x = 0.5, y = -11, label = expression(italic("Mdiff")~"="~-0.86~","~italic("CI")~":"~"("~-1.41~","~-0.29~")"),
                      size = 3.5, hjust = 0, vjust = 1, color = "black")


bpN <- ggplot(CR_CN, aes(x = LifeStage, y = d15N, fill = LifeStage)) +
  stat_boxplot(geom ="errorbar", width = 0.5) +
  geom_boxplot(fill = "grey", outlier.size = 2) + 
  stat_summary(fun=mean, geom="point", shape=10, size=5, color="black") + 
  theme_bw() + theme(legend.position="none", axis.line = element_line(color='black'),
                     plot.background = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.grid.major = element_blank(),
                     axis.text.x = element_text(color="black", size=12),
                     axis.text.y = element_text(size=12),
                     axis.title.x = element_text(size=14),
                     axis.title.y = element_text(size=14)) +
  xlab (NULL) + ylab(expression(delta^15* N * "  \u2030")) 

Nlab <- bpN + annotate("text", x = 0.5, y = 14, label = expression(italic("Mdiff")~"="~0.73~","~italic("CI")~":"~"("~0.19~","~1.30~")"),
                      size = 3.5, hjust = 0, vjust = 1, color = "black")



library(grid)

bottom<-textGrob("Size Class", gp = gpar(fontsize=16))

fig4 <- grid.arrange(Clab, Nlab, nrow=1, ncol=2, bottom = bottom)
print(fig4)

ggsave(filename = "Fig4.eps", fig4, width = 174, height = 150, units = "mm")
dev.off()

ggsave(filename = "Fig4.jpg", fig4, width = 174, height = 150, units = "mm")

