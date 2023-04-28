
### Analysis associated with Weber publication: "Foraging ecology of Kemp’s ridley (Lepidochelys kempii) turtles in the northeastern Gulf of Mexico: insights from stable isotope analysis"


### Codes for SIBER (Code modified from https://cran.r-project.org/web/packages/SIBER/vignettes/Introduction-to-SIBER.html and SIBER vignettes from Andrew L Jackson) ###

#load required packages
library(ggplot2)
library(MASS)
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
library(rjags)


skin<-read.csv("LK_SIA_Data.csv") #load data
summary(skin) #check skin import

colnames(skin)[1] <- gsub('^...','',colnames(skin)[1])

skin_factors<-c("TurtleID", 
                "LifeStage") #select appropriate vars to convert into factors for skin data
skin[skin_factors]<-lapply(skin[skin_factors], as.factor) #convert selected vars to factors


skin_nums<-c("SCL",
             "d15N",
             "d13C")
skin[skin_nums]<-lapply(skin[skin_nums], as.numeric)

summary(skin) #check conversion

set.seed(1) #set the seed so that the results are comparable, same results each time



#organize & set up dataset for SIBER analysis

skin_siber_dfCN <- skin %>% 
  dplyr::select(d13C, d15N, LifeStage) %>%
  mutate(LifeStage = recode(LifeStage, "Large" = 2, "Small" = 1)) %>%
  drop_na(LifeStage)
skin_siber_dfCN$communityCN <- 1 #doesn't mean anything, siber just needs a community column
names(skin_siber_dfCN)<-c("iso1","iso2","group", "community")
summary(skin_siber_dfCN)
skin_siber_dfCN <- as.data.frame(skin_siber_dfCN)
skin_siber_SmallCN <- skin_siber_dfCN %>%
  filter(group==1)
skin_siber_LargeCN <- skin_siber_dfCN %>%
  filter(group==2)

#group column is Large turtles=2,  Small turtles=1

#create siber objects
skin_siber_objCN <- createSiberObject(skin_siber_dfCN) #pooled data, large & small turtles
skin_siber_Small_objCN <- createSiberObject(skin_siber_SmallCN)
skin_siber_Large_objCN <- createSiberObject(skin_siber_LargeCN)

# Create lists of plotting arguments to be passed onwards to each 
# of the three plotting functions.
community.hulls.args <- list(col = 1, lty = 1, lwd = 3)
group.ellipses.args  <- list(n = 100, p.interval = 0.40, lty = 1, lwd = 3)
group.hulls.args     <- list(lty = 2, col = "grey20")

#plot data, can change T/F to add/subtract ellipses/hulls
par(mfrow=c(1,1))
plotSiberObject(skin_siber_objCN,
                ax.pad = 0.1, #spreads out axes to make sure ellipses don't go outside of axes
                hulls = F, community.hulls.args = community.hulls.args, 
                ellipses = T, group.ellipses.args = group.ellipses.args,
                group.hulls = F, group.hulls.args = group.hulls.args,
                bty = "L",
                iso.order = c(1,2),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030'))
legend <- c("Smaller", "Larger")
legend("topright", inset=0.02, legend, 
       pch = c(1,1,1,2,2,2), col = c(1:4, 1:4), lty = 1)


#calculate summary stats for each group (lifestage), TA = Total Area, SEA = Standard ellipse area, SEAc = Corrected standard ellipse area
group.ML <- groupMetricsML(skin_siber_objCN)
print(group.ML)


#calculate layman metrics for communities - not bayesian
community.ML <- communityMetricsML(skin_siber_objCN) 
print(community.ML)

#fit bayesian model to data
# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many (aka the model will choose every 10th value)
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posterior <- siberMVN(skin_siber_objCN, parms, priors)

# The posterior estimates of the ellipses for each group can be used to calculate the SEA.B for each group.
SEA.B <- siberEllipses(ellipses.posterior)

siberDensityPlot(SEA.B, xticklabels = c("Smaller Size Class", "Larger Size Class"),
                 probs = c(95, 75, 50),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 scl = 1,
                 ylim = c(0, 8),
                 bty = "L",
                 las = 1,
)

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.B), group.ML[3,], col="red", pch = "x", lwd = 2)

# Calculate some credible intervals 
cr.p <- c(0.95, 0.99) # vector of quantiles

# call to hdrcde:hdr using lapply()
SEA.B.credibles <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)

print(SEA.B.credibles) #SEAb (Bayesian standard ellipse areas) and credible intervals

# do similar to get the modes, taking care to pick up multimodal posterior
# distributions if present
SEA.B.modes <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = cr.p, all.modes=T)

print(SEA.B.modes) #SEAb modes

#extract SEA for life stages
SEA_group_modes <- as.data.frame(SEA.B.modes, row.names = "SEA", col.names = c("Small", "Large"))
SEA_group_modes <- as.data.frame(t(SEA_group_modes))
SEA_group_modes

#Prob that Adult ellipse is greater than subadults
Pg1.lt.g2 <- sum( SEA.B[,1] < SEA.B[,2] ) / nrow(SEA.B)
print(Pg1.lt.g2)###comparing 1 with 2

#Overlap in ppm sq
overlap.G1.G2 <- maxLikOverlap("1.1", "1.2", skin_siber_objCN, p = 0.40, n =100)
overlap.G1.G2

#Percent overlap
prop.of.first <- as.numeric(overlap.G1.G2["overlap"] / overlap.G1.G2["area.1"])
print(prop.of.first) ###XX% of ellipse 1 is shared with ellipse 2

prop.of.second <- as.numeric(overlap.G1.G2["overlap"] / overlap.G1.G2["area.2"])
print(prop.of.second)###XX% of eliipse 2 is overlap with ellipse 1

prop.of.both <- as.numeric(overlap.G1.G2["overlap"] / (overlap.G1.G2["area.1"] + overlap.G1.G2["area.2"]))
print(prop.of.both) ####XX% of the total area of both ellipse areas are shared

bayes.overlap.G1.G2 <- bayesianOverlap("1.1", "1.2", ellipses.posterior, 
                                       draws = 10, p.interval = 0.40,
                                       n = 360)
print(bayes.overlap.G1.G2)

##Export Metrics
skin_siber_metrics<-as.data.frame(t(group.ML))
rownames(skin_siber_metrics)<-c("Small", "Large")
skin_siber_metrics$SEAb<-SEA_group_modes[]
names(skin_siber_metrics[,4])<-"SEAb"
write.csv(skin_siber_metrics, file="skin_siber_metrics.csv")


##### Create ellipse plot in ggplot ########

first.plot <- ggplot(data = skin,
                     aes(x = d13C, 
                         y = d15N)) + 
  geom_point(aes(color = LifeStage), size = 5) +
  ylab(expression(paste(delta^{15}, "N (\u2030)"))) +
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  theme(text = element_text(size=16)) + 
  scale_color_brewer(palette="Set1") 

# And print our plot to screen
print(first.plot)

# use our ellipse function to generate the ellipses for plotting

# decide how big an ellipse you want to draw
p.ell <- 0.40

#the following is taken from Jackson 2023 -> https://cran.r-project.org/web/packages/SIBER/vignettes/Introduction-to-SIBER.html
# create our plot based on first.plot above adding the stat_ellipse() geometry.
# We specify thee ellipse to be plotted using the polygon geom, with fill and
# edge colour defined by our column "group", using the normal distribution and
# with a quite high level of transparency on the fill so we can see the points
# underneath. In order to get different ellipses plotted by both columns "group"
# and "community" we have to use the interaction() function to ensure both are
# considered in the aes(group = XYZ) specification. Note also the need to
# specify the scale_fill_viridis_d() function as the mapping of colors for
# points and lines is separate to filled objects and we want them to match.
ellipse.plot <- first.plot + 
  stat_ellipse(aes(group = interaction(LifeStage), 
                   fill = LifeStage, 
                   color = LifeStage), 
               alpha = 0.25, 
               level = p.ell,
               type = "norm",
               geom = "polygon") + 
  scale_color_brewer(palette="Set1")

print(ellipse.plot)

classic.ellipse.plot <- ellipse.plot + theme_classic() + theme(text = element_text(size=18))

# and print to screen
print(classic.ellipse.plot)


