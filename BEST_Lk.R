
library(BEST)
library(dplyr)

### Notes taken from Meredith and Kruschke 2021, https://cran.r-project.org/web/packages/BEST/vignettes/BEST.pdf

### Needs a vector of data values, copy & pasted these from my csv data file ###

# Nitrogen first#

smallN <- c(6.528674887,
8.496798689,
6.323200008,
9.1025748,
7.285840164,
7.748962605,
7.100403502,
6.718640705,
7.844626525,
8.491641849,
6.985143997,
6.449491885,
7.468912274,
5.805489804,
7.21644874,
7.184267211,
6.916571366)

largeN <- c(5.502360256,
            6.917763919,
            5.371208598,
            7.889918106,
            7.504071161,
            5.767274237,
            6.056473187,
            5.11683813,
            5.408238164,
            6.892473972,
            4.737006463,
            7.9204197,
            6.553710366,
            7.729341868,
            5.649014206,
            5.86972964,
            5.721013288,
            5.086761249,
            5.687826086,
            11.1066497,
            6.163061581,
            6.596954296,
            6.201724886,
            6.326619686,
            5.8869493,
            5.829151175,
            6.911729882,
            8.3086094,
            5.978407521,
            8.3086094,
            6.87150297,
            6.543709635,
            6.987735837,
            6.370591438,
            7.228358991,
            5.692247078,
            6.758867617,
            6.887327201,
            7.899182293,
            6.28103763,
            6.125246899,
            11.75968306,
            7.532086228,
            4.263717579,
            6.261026709,
            7.200211876,
            13.29882754)

# Based on previous experience with these sort of trials, we expect reaction times to be approximately 6 secs, but they vary a lot, so we’ll set muM = 6 and muSD = 2. 
# We’ll use the default priors for the other parameters where y = c(y1, y2) 
# sigmaMode = sd(y)
# sigmaSD = sd(y)*5
# nuMean = 30
# nuSD = 30) 
#priors <- list(muM = 6, muSD = 2)

### But we don't the priors...#

# run the model 
BESToutN <- BESTmcmc(smallN, largeN, parallel=FALSE)

# Infer results
plot(BESToutN) # histogram of the posterior distribution of the difference in mean
# Also shown is the mean of the posterior probability, which is an appropriate point estimate of the true difference in means, the 95% Highest Density Interval (HDI), and the posterior probability that the difference is greater than zero. The 95% HDI does not include zero (which means there IS a credible difference) and the probability that the true value is greater than zero is shown as 99.5%.


plot(BESToutN, compVal=0.5, ROPE=c(-0.1,0.1)) # The probability that the difference in reaction times is precisely zero is zero. More interesting is the probability that the difference may be too small to matter. We can define a region of practical equivalence (ROPE) around zero, and obtain the probability that the true value lies therein. For the reaction time example, a difference of ± 0.1 may be too small to matter. The plot shows a high probability that the reaction time increase is > 1. In this case it’s clear that the effect is large, but if most of the probability mass (say, 95%) lay within the ROPE, we would accept the null value for practical purposes.

#I made comp val 0.5

plot(BESToutN, which="sd") # BEST deals appropriately with differences in standard deviations between the samples and departures from normality due to outliers. We can check the difference in standard deviations or the normality parameter with plot 


summary(BESToutN) # Here we have summaries of posterior distributions for the derived parameters: difference in means (muDiff), difference in standard deviations (sigmaDiff) and effect size (effSz).

summary(BESToutN, credMass=0.8, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
        compValeff=1) #As with the plot command, we can set values for compVal and ROPE for each of the parameters of interest

class(BESToutN)
print(BESToutN) #The print function displays the mean, standard deviation and median of the posterior distributions of the parameters in the model, together with a 95% Highest Density Interval

plotPostPred(BESToutN) # Each panel of Figure 5 corresponds to one of the samples, and shows curves produced by selecting 30 random steps in the MCMC chain and plotting the t distribution with the values of µ, σ and ν for that step. Also shown is a histogram of the actual data. We can visually assess whether the model is a reasonably good fit to the sample data (though this is easier for large samples then when n = 6 as here).

plotAll(BESToutN) # puts histograms of all the posterior distributions and the posterior predictive plots onto a single page

names(BESToutN)
meanDiff <- (BESToutN$mu1 - BESToutN$mu2)
meanDiffGTzero <- mean(meanDiff > 0)
meanDiffGTzero


#Carbon

smallC <- c(-14.92919778,
            -15.36107241,
            -14.61685145,
            -16.31469173,
            -14.90466657,
            -15.68898257,
            -13.85466415,
            -14.31819414,
            -15.36660456,
            -15.09062992,
            -15.42320042,
            -14.28737645,
            -15.47037535,
            -16.59904759,
            -12.76703438,
            -14.90642565,
            -15.86119892
)
largeC <- c(-13.19805836,
            -14.17699594,
            -14.59256993,
            -15.64678569,
            -13.84323825,
            -14.60775873,
            -13.74988876,
            -13.82041979,
            -14.58788648,
            -14.99614284,
            -14.8332227,
            -14.55736518,
            -13.28191851,
            -14.30267146,
            -14.46898723,
            -14.74451713,
            -12.50011642,
            -14.60253699,
            -14.38956284,
            -15.79102709,
            -14.26790995,
            -13.65613167,
            -11.43103635,
            -14.24149768,
            -13.0113063,
            -14.47728889,
            -13.21670588,
            -12.55486301,
            -14.74813065,
            -14.46801125,
            -13.05488383,
            -14.21886177,
            -14.73372943,
            -14.93076724,
            -16.79512637,
            -14.64332215,
            -12.63638967,
            -14.45653473,
            -15.50046888,
            -15.01482158,
            -13.85724446,
            -14.15067069,
            -15.17962645,
            -11.86816264,
            -13.12289922,
            -14.32370208,
            -15.35282206
)

BESToutC <- BESTmcmc(smallC, largeC, parallel=FALSE)
plot(BESToutC)

plot(BESToutC, compVal=0.5, ROPE=c(-0.1,0.1)) # The probability that the difference in reaction times is precisely zero is zero. More interesting is the probability that the difference may be too small to matter. We can define a region of practical equivalence (ROPE) around zero, and obtain the probability that the true value lies therein. For the reaction time example, a difference of ± 0.1 may be too small to matter. The plot shows a high probability that the reaction time increase is > 1. In this case it’s clear that the effect is large, but if most of the probability mass (say, 95%) lay within the ROPE, we would accept the null value for practical purposes.

#I made comp val 0.5

plot(BESToutC, which="sd") # BEST deals appropriately with differences in standard deviations between the samples and departures from normality due to outliers. We can check the difference in standard deviations or the normality parameter with plot 


summary(BESToutC) # Here we have summaries of posterior distributions for the derived parameters: difference in means (muDiff), difference in standard deviations (sigmaDiff) and effect size (effSz).

summary(BESToutC, credMass=0.8, ROPEm=c(-0.1,0.1), ROPEsd=c(-0.15,0.15),
        compValeff=1) #As with the plot command, we can set values for compVal and ROPE for each of the parameters of interest

class(BESToutC)
print(BESToutC) #The print function displays the mean, standard deviation and median of the posterior distributions of the parameters in the model, together with a 95% Highest Density Interval

plotPostPred(BESToutC) # Each panel of Figure 5 corresponds to one of the samples, and shows curves produced by selecting 30 random steps in the MCMC chain and plotting the t distribution with the values of µ, σ and ν for that step. Also shown is a histogram of the actual data. We can visually assess whether the model is a reasonably good fit to the sample data (though this is easier for large samples then when n = 6 as here).

plotAll(BESToutC) # puts histograms of all the posterior distributions and the posterior predictive plots onto a single page

names(BESToutC)
meanDiff <- (BESToutC$mu1 - BESToutC$mu2)
meanDiffGTzero <- mean(meanDiff > 0)
meanDiffGTzero


