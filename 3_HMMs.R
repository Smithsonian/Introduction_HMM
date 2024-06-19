# ******************************************************************
# ******************************************************************

# Project: Intro to Animal Mvmts
# Description: Fit Hidden Markov Models with moveHMM.  Code based on exercises developed by Joe Kolowski with sections modified from the moveHMM guide.  Many thanks to Joe and Theo Michelot for the helpful information provided.
# Author: Jared Stabach
# Date: 16 June 2024

# ******************************************************************
# ******************************************************************

### Objectives:   
# 1. Learn the analysis workflow for fitting Hidden Markov Models with moveHMM
# 2. Practice investigation of impact of continuous covariates on state transition probabilities 
# 3. Learn to interpret output from HMMs
# 4. Practice methods for visualizing HMM output both spatially and graphically
# 5. Use model output to subset tracking datasets

### DataSets used in Script:
# White-bearded Wildebeest data from the Athi-Kaputiei Plains in Kenya equipped with Lotek WildCell Satellite GPS collars, resampled to 3hr intervals in "2_CreateTrajectory.R".

# Setup ------------------------------------------------------------

# Remove everything in memory
rm(list = ls())

library(moveHMM)
library(ggmap)
library(tidyverse)


# Data import and formatting ---------------------------------------

# Here I'm bringing in cleaned wildebeest data that has been resampled to 3 hr stepse. This is resampled data we saved from the earlier script, since are investigating the influence of time of day here, and we needed regular sampling throughout the 24 period (even though animals were sampled at hourly intervals during the day). Remember this file now includes only animals from one of the three study areas.

WB.data <- read_rds("Data/wildebeest_3hr.rds")
#WB.data <- read_rds("Data/wildebeest_3hr_adehabitat.rds") # Same thing, except with adehabitat.  Columns are slightly different

# In terms of data prep, a few processing/formatting steps are required: 

# First it is assumed that location error is relatively small. If you are working with a high error data set, you should likely process it before using this package. For example, you may try a state space approach (e.g. crawl, aniMotum). 

# Second, note that the package does not use date/time information, and so assumes the steps are regular. Therefore, data must be organized chronologically.

# Third, variables included in analyses must be numeric.  For example, seasons should be change to 1, 2, 3, and 4 (etc). 

# Fourth, the animal id column must be named "ID".

# Lastly, it seems very helpful to convert the coordinates (currently in meters) to km by dividing by 1000.  All moveHMM vignettes do this conversion of units.  I could not get the functions to work without doing so.

# Here I will assign the data to an object and perform all the necessary reformatting. I'm converting to a dataframe here because this is actually a track object from amt since we used that package for resampling. I'm also creating a column for hour of the day, which we will use below for modeling.

WB.data <- read_rds("Data/wildebeest_3hr.rds") %>%
  as.data.frame() %>% 
  arrange(id, t_) %>%
  select(x = x_,
         y = y_,
         t = t_,
         ID = animal_id,
         temp,
         sex) %>% 
  mutate(hour = lubridate::hour(t),
         x = x/1000,
         y = y/1000,
         sex = sex)


# Create move object and inspect ------------------------------------------

# The function to prepare data to be analyzed in this package is called "prepData()". Here we could have used the raw lat/long values here as well, and then we would indicate the type = "LL".

WB.move <- WB.data  %>% 
  prepData(type = "UTM",
           coordNames = c("x","y"))

class(WB.move)
summary(WB.move)

# We now have movement data for all 12 of our animals and two covariates (temperature and hour of day) to use in analyses. We can also see that the format is a data.frame. This makes it relatively easy to modify our object if we need to. 

names(WB.move)

# Notice that we now have two new columns, "step" and "angle", representing step lengths and turning angles. This is what will be used in all future modeling with this package.

# We can start with an overall plot of our data set. Multiple animals are plotted in separate windows. If ask = TRUE each plot will come separately and you'll need to hit enter after each one, which becomes tedious. Using the argument "compact = T" plots all the tracks at once instead of each animal's track separately. Either way, you will see each animal's plot of step lengths and turning angles (in radians) over time, as well as histograms of step length and turning angle.

plot(WB.move,
     ask = FALSE)

plot(WB.move,
     compact = TRUE,
     ask = FALSE)

# We can plot them separately like this.
plot(WB.move[WB.move$ID == "Kikaya",],
     ask = FALSE)

# We can also investigate the new step and turning angle with our own code.
hist(WB.move$step) 
summary(WB.move$step)

# This indicates one animal moved 23 km in one step and this may be resulting from a gap in the data unfortunately. But this package does not keep track of "bursts" and so there really isn't a way to address this unfortunately.

quantile(WB.move$step, 
         probs = 0.90, 
         na.rm = T) 

# 90% of steps below 1 km. Remember this is 3 hour data.

hist(WB.move$angle)


# Fit Hidden Markov Models - set starting values --------------------------

# To begin our modeling, we need to set starting values for all the relevant parameters that will be estimated so that the optimization algorithm has a starting point. 

# In this case, the initial parameters should be specified in two vectors, stepPar0 (for the step distribution) and anglePar0 (for the angle distribution). We'll start here with some basic default values that often work well for animal movement with two behavioral states. 

# To determine what parameters need to be set, we need to decide which statistical distribution will be used to describe our step lengths and our turning angles. Different distributions have different numbers of parameters. 

# Here we will use a gamma distribution for step lengths, and a von Mises distribution for turning angles. 

# For the gamma distribution, we need to specify a mean and SD. For the von Mises, we need a mean and concentration parameter.

# Zero-inflation must be included in the step length distribution if some steps are of length exactly zero. To do so, another parameter is added to the step distribution: its mass on zero. Let's see if any of these animals have zero step lengths:

# Here I'm looking at the 10 lowest values of step. If "with_ties = T" the function will not count ties for step value and so you'll end up with more than 10 rows in the output.
slice_min(WB.move,
          order_by = step,
          n = 10,
          with_ties = F)

# We do have a single value of zero here, so we'll need to include this extra parameter.

# Here, the initial values are chosen such that they correspond to the commonly observed pattern in 2-state HMMs for animal movement data, with state 1 involving relatively short steps and many turns (hence the choice of a small initial value for the mean of the gamma step length distribution and an initial value of pi radians (180 degrees) for the mean turning angle) and state 2 involving longer steps and fewer turns (hence the choice of a larger initial value for the mean of the gamma step length distribution and an initial value of 0 for the mean turning angle).

# Note that some authors, to ensure starting values did not overly impact results, will repeat model fitting with a range of starting values. This allows one to test the sensitivity of results to starting value decisions.

# Step distributions

# step mean (two parameters: one for each state)
mu0 <- c(0.1, 1) 

# Step standard deviations
sigma0 <- c(0.1, 1)

# Zero-distribution term. These would be higher if there were frequent zeros. 
zeromass0 <- c(0.1, 0.05) 

# Assigning the step distribution starting values
stepPar0 <- c(mu0,
              sigma0,
              zeromass0)

# Turning angle starting values

# Mean turning angle. In radians, pi, or 3.14 represents 180 degrees.
angleMean0 <- c(pi, 0)

# Angle concentration. This reflects the variance in the turning angle distribution.
kappa0 <- c(1, 1) 

# Setting selected parameters
anglePar0 <- c(angleMean0,
               kappa0)


# Let's set up parameters for a 3 behavior model as well, one value for each behavior state.

mu0_3 <- c(0.1, 0.5, 3)
sigma0_3 <- c(0.05, 0.5, 1)
zeromass0_3 <- c(0.05, 0.0001, 0.0001)
stepPar0_3 <- c(mu0_3, sigma0_3, zeromass0_3)
angleMean0_3 <- c( pi, pi, 0)
kappa0_3 <- c(1, 1, 1)
anglePar0_3 <- c(angleMean0_3, kappa0_3)

# For numerical stability, we could decide to standardize the covariate values (in this case temp and hour of the day) before fitting the model. But here we won't be combining them in any models.

# Fit Hidden Markov Models - Model Fitting --------------------------------

# When fitting models we can provide a formula with a covariate to inform transition from one state the next, Here we will investigate whether time of day or temperature influences the probability that animals transition from one behavioral state to the next.

# We also have to specify how many behavioral states we want the model to identify. Here we will indicate "nbStates = 2" to start,  which specifies that we want to fit a two-state model to the data.

# Next we want to define the distributions that we want to use to characterize both the step lengths and turning angles. As noted above, we are going to use a gamma distribution for the former (stepDist = "gamma") and a von-Mises distribution for the latter (angleDist = "vm"). These are the default values but I'm specifying them below. 

# Distribution options for step length are: gamma (“gamma”), Weibull (“weibull”), exponential (“exp”) and log-normal (“lnorm”).

# For turning angle we can choose: von Mises (“vm”) and wrapped-Cauchy (“wrpcauchy”). It is also possible to specify angleDist = "none", if the angles are not modeled.

# Here we will first fit a kind of NULL model, where there are two behavioral states, but no covariates influencing transitions from one state to the next. Note that my formula has an intercept only (~ 1). This follows standard regression formula specification.

WB.null <- fitHMM(data = WB.move, 
                      nbStates = 2, 
                      stepPar0 = stepPar0, 
                      anglePar0 = anglePar0, 
                      stepDist = "gamma",
                      angleDist = "vm",
                      formula = ~ 1)

WB.null

# OK so we have a number of values here summarizing this initial model. We have mean step lengths from each predicted state, along with their SD and zero-mass. State one has a much smaller step length mean around 285m (may vary with each run), with mean for state 2 around 700m. You can see the zero-mass values are very small, which makes sense since we had only one zero step length. 

# Then we have turning angle means and concentrations, with the concentration representing a kind of variance around the mean. State 2 has a mean close to zero. Remember that a relative turn angle of zero means continued movement in the same direction.

# We did not test a regression equation, so we only have an intercept regression coefficient.

# Then finally we have a transition matrix. This shows the average probability of transitioning from state 1 to state 2 is fairly low, at 2.6%. Transitioning from state 2 to state 1 is only slightly higher at 4.6% indicating animals tend to continue in their state rather than frequently switching. This makes sense from a behavioral standpoint and is what we would typically expect in these models.

# Let's take a look at the plots. Here remember as you scroll through that state 1 is orange, state 2 is blue. We can see that path segments are colored according to predicted state in all the full trajectory plots of each animal. And we can see from these that longer distance straight movements are colored as state 2, with more resident/locally focused movements as state 1 (orange). It is likely that here the orange state 1 represents more of a local foraging state, while state 2 (blue) is directed travel. You'll need to hit the back button in the plotting window to see all the plots.

plot(WB.null,
     ask = F)

# We can also see plots of the distributions of step lengths and turning angles for our two states.


# State Assignments and Probabilities -------------------------------------

# To globally decode the state process, the Viterbi algorithm is implemented in the function viterbi(). This function outputs the most likely sequence of states to have generated the observations, under the fitted model. 

WB.states.2 <- viterbi(WB.null)

# Below are the most probable states for the first 25 observations of the first individual. And we can see they are all state 1. 
WB.states.2[1:25]

# This is a very useful output of these models. We now have a predicted state for each location in the dataset. This allows us to look at the proportion of time spent in each state for example.

prop.table(table(WB.states.2))

# So overall, across all animals these wildebeest spend about 64% of their time in this resident, likely intensive foraging state, We could also look at this across individuals, which might indicate range residency vs. more nomadic or even migratory behavior.

# To look at this metric across individuals we can add the state information back to our data frame, as our states vector is the same length as our data. Here I'm making a new object called stateProps which adds in the new state information, then summarizes the proportion of time spent in each state for each animal. Note the use of .by in the mutate function. This is incredibly handy way of creating a new column of data, that is itself based on a grouping of other data.

stateProps <- WB.data %>% 
  
  # add predicted state information to initial data frame
  mutate(state = WB.states.2) %>%
  
  # add new column that is the total locations for each animal
  mutate(locs = n(),
         .by = ID) %>% 
  
  # summarize for each animal, and each state, the proportion of locations.  Using reframe, as summarize has been deprecated in latest version of dplyr
  reframe(stateProp = n()/locs,
          sex = unique(sex),
          .by = c(ID, state)) %>% 
  
  # This reduces our data frame from the same initial size, to one with just the unique rows of information.
  distinct()

# Visually from the maps, PAITA looked to be one of the most resident animals, so I would predict this animal would spend the most time in state 1. And that is indeed the case, as she spends nearly 83% of her time in this state.

# Let's plot this information across individuals. Note that here I'm coloring by sex. The information provided here could allow you to ask questions about sex differences, or other individual differences in movement strategy and movement behavior. We can't include a variable like sex in our models here, because the covariate needs to vary time/location. But we can still make comparisons using the behavior state model results.

# I'm ordered the ID column here by the value of stateprop. This is a convenient way to arrange values in a plot to be ordered according to some other value. 

stateProps %>%
  filter(state == 1) %>% 
  mutate(ID = fct_reorder(ID, 
                          stateProp)) %>% 
  
  ggplot(aes(y = stateProp,
             x = ID, 
             fill = sex)) +
  geom_col(col = "black",
           position = position_dodge()) +
  coord_flip() +
  labs(y = "Prop. time in State 1 (foraging/encamped)") +
  theme_bw()

test <- stateProps %>%
  filter(state == 1) %>% 
  mutate(ID = fct_reorder(ID, 
                          stateProp))

mycols <- c("light blue", "dark blue")

test.or <- test[order(test$stateProp,decreasing = FALSE),]
barplot(test.or$stateProp,
        names.arg = test.or$ID,
        horiz = TRUE,
        xlim=c(0,1),
        xlab="Proportion time in State 1 (Foraging/Encamped)",
        ylab = "Animal ID",
        col = mycols[test.or$sex])
legend(0.6, 3, legend = c("female", "male"), fill = mycols, horiz=FALSE, bty = "n", cex = 1.5)

barplot(test$stateProp,
        names.arg = test$ID,
        horiz = TRUE,
        xlim=c(0,1),
        xlab="Proportion Time in State 1 (Foraging/Encamped)",
        ylab = "Animal ID",
        col = mycols[test$sex])
legend(0.6, 3, legend = c("female", "male"), fill = mycols, horiz=FALSE, bty = "n", cex = 1.5)

        
# So it looks like Kiranto spends the most time in this local foraging state. Let's remind ourselves what his trajectory looks like.

plot(WB.move[WB.move$ID == "Kiranto",],
     ask = FALSE)

# So this is interesting. Just because he spent the most time in state 1 does not mean he spent all his time in the same place. The animal clearly made some long distance movements during the tracking period. While Paita, despite having a somewhat lower proportion of time spent in the resident foraging state, has a much smaller and more defined home range area.

plot(WB.move[WB.move$ID == "Paita",],
     ask = FALSE)

# Differntiating further between these two animals in terms of overall movement strategy, would be better done with something like a net-squared displacement model, or even mean squared displacement. With these HMM models we are more interested in categorizing segments as particular behavioral states.

# It also looks like the 3 wildebeest that spent the least time in foraging/encamped state were males. Ans the wildebeest that spent the MOST time in this state were also males. This kind of analysis could also be split in different times of year or different breeding seasons, and in fact breeding state could be used in the models themselves, as it would change over time. 

# We can also look at the actual probabilities assigned to each location for each behavioral state.

WB.probs2 <- stateProbs(WB.null)

# This is output as a matrix, but we still have one row for each location. Values across a row should add to 1.0.
head(WB.probs2)
nrow(WB.probs2)

# The package guide states that the state with highest probability according to stateProbs() might not be the same as the state in the most probable sequence returned by the Viterbi algorithm. This is because the Viterbi algorithm performs “global decoding”, whereas the state probabilities are “local decoding”. For more details, see Zucchini et al. (2016).

# The function plotStates() can be used to visualize the results of viterbi() and stateProbs(). This plot shows the plots of the most likely state sequence decoded by the Viterbi algorithm, as well as both columns of the matrix of state probabilities, for one individual, “Kiranto”.

plotStates(WB.null,
           animals = "Kiranto",
           ask = FALSE)

# Let's compare this to an animal that spent much more time in state 2, like Sotua

plotStates(WB.null,
           animals = "Sotua",
           ask = FALSE)

# Note that this kind of information could be associated with all kinds of variables like temperature, NDVI, breeding state and so on, to predict what factors influence change in behavioral state.

# Although the package creates automatic plots of trajectories of each animal, with each point colored by predicted state, we have more flexibility if we generate these plots ourselves. Here are some ideas for making custom plots of this information:

WB.data %>% 
  mutate(state = as.factor(WB.states.2)) %>%
  filter(ID == "Sotua") %>% 
  ggplot(aes(x = x,
             y = y,
             fill = state,
             col = state)) +
  geom_path(alpha = 0.5) +
  geom_point(shape = 21,
             alpha = 0.8,
             col = "black") +
  scale_fill_manual(values = c("orange",
                               "cornflowerblue")) +
  scale_color_manual(values = c("orange",
                                "cornflowerblue")) +
  theme_minimal() +
  labs(x = "x",
       y = "y",
       title = "Sotua")





# Model Comparisons and Covariate Testing -----------------------------------------------
# Let's start by looking at the influence of temperature, which was recorded by the GPS collar, on behavioral state transitions. Here I'm indicating the same parameters as above, but now I'm indicating a formula with temperature as a covariate.

wild_m_temp <- fitHMM(data = wild_move, 
                      nbStates = 2, 
                      stepPar0 = stepPar0, 
                      anglePar0 = anglePar0, 
                      formula = ~ temp)


wild_m_temp

# This is the canned function in the package to plot the impacts of our covariate.
plotStationary(wild_m_temp,
               plotCI = TRUE)

# We can see that the probability of transitioning to state 1 (the encamped/foraging state) increases as the temperature increases. Of course I don't know how accurate this temperature sensor is, but at least this gives you a sense of how to test a model with covariates. And one could certainly reason that the animals would be less likely to move longer distances when the temperature was very warm.

# Let's now test a model with time of day. Hour of the day is a challenging variable to model in a regression equation because the 23rd hour is very close to the zero hour, but if it is modeled as a regular number, this connection will be lost. The conversion below allows us to model the hour of the day in a circular format. You can read a succinct description of this here, even though the coding part of the page is not in R:
# https://ianlondon.github.io/blog/encoding-cyclical-features-24hour-time/

wild_m_tod <- fitHMM(data = wild_move, 
                     nbStates = 2, 
                     stepPar0 = stepPar0, 
                     anglePar0 = anglePar0, 
                     formula = ~ cos(2*pi*hour/24) + 
                       sin(2*pi*hour/24))

wild_m_tod

# Let's see what the impact of time of day is.
plotStationary(wild_m_tod,
               plotCI = T)

# This plot is not quite right actually. Read through the sidebar at the end of the script if you want more detail on why, and how to fix it.

# Let's now see which of these three models is a better fit to the data
AIC(wild_m_null,
    wild_m_temp,
    wild_m_tod)

# So the null model is preferred over the temp model, but there is strong evidence that time of day is important in determining behavioral state.

# So now let's end by testing to see if there is evidence that a three state model is favored. We have evidence that time of day is important, so we'll keep this variable. And here we are feeding in starting parameters for the 3 state model, using the same gamme and Von Mises distributions.

wild_m_tod_3 <- fitHMM(data = wild_move, 
                       nbStates = 3, 
                       stepPar0 = stepPar0_3, 
                       anglePar0 = anglePar0_3, 
                       stepDist = "gamma",
                       angleDist = "vm",
                       formula = ~ cos(2*pi*hour/24) + 
                         sin(2*pi*hour/24))


AIC(wild_m_tod, wild_m_tod_3)

# So there is actually substantial evidence here for 3 states.

wild_m_tod_3

# Interestingly, it seems that the orange state 1 may represent a kind of middle level foraging with relatively low displacement. Step length mean is around 280m and turning angle mean is quite large at 2.88. State 2 blue has an even smaller step length mean of around 170m, and a similarly high turning angle. Perhaps this represents highly focal small scale movements, or perhaps even some sort of resting state. State 3 green is representing large scale directional movements with a mean length of 810m and a turning angle near zero.

# Let's see how much time animals tend to spend in each state.
wild_states_3tod <- viterbi(wild_m_tod_3)

prop.table(table(wild_states_3tod))

# So across all individuals, state 1, the kind of middle distance focal foraging movement, is the most common, followed by, interestingly, the longer distance directional movement. I would imagine this varies quite a bit across individuals though and we could certainly replicate our graph above for these three states to look across individuals and sexes.

# Let's see how this looks on a map
plot(wild_m_tod_3,
     ask = FALSE)

# And let's look at the time of day influence on each state probability
plotStationary(wild_m_tod_3)

# It looks like in general our new state 3, which appears to represent large scale directional movements, is most likely to occur between 5 and 7pm, dusk in eastern Africa. State 2, which is the likely focal foraging state is expected most of the day except around 3pm to 9pm,

# Selection of Locations Based on State -----------------------------------

# We'll end by demonstrating how one might remove locations associated with a certain behavioral state. This could be particularly helpful if one is trying to create home ranges for the resident/foraging behavior of an animal, but long-distance travel between ranges is causing unreasonable resident range sizes.

# We first need to assign state information to each location. We've already shown this above, but now we need to preserve the rest of the data. Then I'm multiplying my coordinates by 1000 to return to our original grid values. Here I'm working with the original dataset since I don't need the step and angle information.

# In the 3-state model, state 3 represented directional, long range movements. We'll remove these from our new dataset

wild_data_foraging <- wild_data %>% 
  
  # add predicted state information to initial data frame
  mutate(state = factor(wild_states_3tod),
         x_ = x_ * 1000,
         y_ = y_ * 1000) %>%
  
  # remove hour column
  select(-hour) %>% 
  
  # filter data to only states 1 and 2
  filter(state %in% c(1,2))

# An animal like Ntishya had a few focal areas of use, with substantial travel in between. Let's see if removing state 3 helps us focus in on areas of concentrated foraging.

ntishya_all <- wild_data %>% 
  mutate(state = as.factor(wild_states_3tod)) %>%
  filter(ID == "Ntishya") %>% 
  ggplot(aes(x = x_*1000,
             y = y_*1000,
             col = state,
             fill = state)) +
  geom_path(alpha = 0.5) +
  geom_point(shape = 21,
             alpha = 0.8,
             col = "black") +
  scale_color_manual(values = c("orange",
                                "cornflowerblue",
                                "darkolivegreen")) +
  scale_fill_manual(values = c("orange",
                               "cornflowerblue",
                               "darkolivegreen")) +
  theme_void() +
  labs(x = "x",
       y = "y",
       title = "Ntishya - all data")

ntishya_all

# Let's see what her map looks like without the "travel" state.

ntishya_sub <- wild_data_foraging %>% 
  filter(ID == "Ntishya") %>% 
  ggplot(aes(x = x_,
             y = y_,
             col = state,
             fill = state)) +
  #geom_path(alpha = 0.5) +
  geom_point(shape = 21,
             alpha = 0.8,
             col = "black") +
  scale_color_manual(values = c("orange",
                                "cornflowerblue")) +
  scale_fill_manual(values = c("orange",
                               "cornflowerblue")) +
  theme_void() +
  labs(x = "x",
       y = "y",
       title = "Ntishya - foraging state only")

ntishya_sub

# We can use the cowplot package to combine these plots next to each other.

two_plots <- cowplot::plot_grid(ntishya_all,
                                ntishya_sub,
                                ncol = 2,
                                rel_widths = c(1, 1)
)

two_plots

# This is an example of how one might save a plot like this for a report or publication:

ggsave(plot = two_plots,
       filename = "output/plots/ntishya_2maps.tiff",
       units = "in",
       width = 6.5,
       height = 4)

# This does seem to have helped, but we still have quite a few points in middle areas. This kind of approach might work better for an animal that forages in more focal, selected habitats, compared to a wide-open herd animal grazing like a wildbeest, but at least you can see the approach here.


# PLOTTING AND PREDICTING SIDEBAR - SKIP IF YOU LIKE ####

# In fact in this case the automatic plotting decisions do not work in our case because we did not have any times in the 22nd and 23rd hours. Based on how this function decides what to put on the x axis, we don't get to see the full picture. Because of the circular nature of the time predictor, predicted values in the 23rd hour should approach values in the zero hour. So here I'm going to show how to extract the values we need directly, using other moveHMM functions.

# There are two predicting functions in the package. Those familiar with regression in R will recognize how to use the predict() world of functions. We need to provide a model, which we have, and a set of data to predict values for, using the model equation. We always need to provide a data frame, and it needs to have columns for each covariate we specified in the model, named exactly the same as our data. So here we just need a column with "hour". To make sure we predict state probabilities at each hour, I'm going to specify a sequence of values from 0 to 23, with a total of 100 values in the range.

predValues_hours <- data.frame(hour = seq(from = 0,
                                          to = 23,
                                          length = 100))


# We could use the predictTPM() function to predict the transition probability matrix, but we want to predict the stationary probabilities, that is, the predicted probabilities for any given "hour". For this we need the predictStationary() function. We need to provide our model, and the data we want to predict for.

pred_mod_hour <- predictStationary(m = wild_m_tod, 
                                   newData = predValues_hours, 
                                   returnCI = TRUE, 
                                   alpha = 0.95)

# The result is a list with the following names. The "mle" is our estimate, and the other two are lower and upper confidence limits.
names(pred_mod_hour)

# Now we want to plot this, so we can get a sense of how hour of the day impacts state probabilities. First I'll make a data frame of our input data, along with our predicted values. Here I have to select the first column, because that is associate with state 1. I could plot both states on the plot as well, but in this case they are essentially mirror images of each other. The CIs don't work here because we don't have data from every hour, but this code would work for CIs if you did.

pred_plot_data <- data.frame(hour = predValues_hours$hour,
                             state1prob = pred_mod_hour$mle[,1],
                             lci = pred_mod_hour$lci[,1],
                             uci = pred_mod_hour$lci[,1])

ggplot(pred_plot_data,
       aes(x = hour,
           y = state1prob)) +
  geom_ribbon(aes(ymin = lci, 
                  ymax = uci), 
              alpha = 0.3) +
  geom_line(col = "#D55E00") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "time of day", 
       y = "stationary probability") +
  scale_x_continuous(breaks = seq(0, 24, by = 4)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5), 
                     labels = c("0", "0.5", "1"))

# Remember that state 1 is the foraging, non-traveling state. This indicates the probability of being in this state declines dramatically during the dusk hours, when animals are apparently much more likely be to moving. This is roughly from 4pm to 9pm. 