# ******************************************************************
# ******************************************************************

# Project: Intro to Animal Mvmts
# Description: Create Animal Movement trajectory for modeling
# Author: Jared Stabach
# Date: 15 June 2024

# ******************************************************************
# ******************************************************************

# Setup ------------------------------------------------------------

# Remove anything in memory
rm(list=ls())

library(sf)
library(amt)
library(adehabitatLT)
library(tmap)
library(cowplot)
library(tidyverse)
library(suncalc)

# Set TimeZone and UTM Zone
Timezone1 <- 'UTC'
Timezone2 <- "Africa/Nairobi"

# UTM Zone
LatLong.proj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"  # EPSG:4326
UTMZone.proj <- "+proj=utm +south +zone=37 +units=m +datum=WGS84 +init=EPSG:32737"


# Data Import ------------------------------------------------------
load("./Data/wildebeest_data.rdata")

# Let's load the sf file.  It is essentially a dataframe that is projected (has geographic coordinates)
WB.sf
str(WB.sf)

# Remember that we have projected these data to UTM Zone 37S.  
# These are important items to remember.
# Each dataset will have it's own geographic location.  Important to set this value correctly.

# Check projection
st_crs(WB.sf)


# AMT (Animal Movement Tools): -------------------------------------

# Prepare and Inspect the Track
# Note that we are loading/working with the sf object (projected data)
# These same steps could be done with the non-projected data.frame (WB)

# As a first step, we will convert the sf object to a tibble and pull the projected coordinates from the object.
# We'll then remove the geometry, essentially turning the sf object back into our original dataframe.

WB.data <- WB.sf %>%
  as_tibble() %>% 
  mutate(X = st_coordinates(WB.sf)[ ,1],
         Y = st_coordinates(WB.sf)[ ,2]) %>% 
  dplyr::select(-geometry)

# Create track for a single animal
WB.track.paita <-  WB.data %>%
  filter(animal_id == "Paita") %>%
  make_track(.x = X,
             .y = Y,
             .t = timestamp,
             crs = UTMZone.proj,
             all_cols = TRUE)

# What's the class of this object?
class(WB.track.paita)

# We can see that we have a new kind of object here. While a data frame, it also has the properties of other object classes, which allow new functionality through the amt package. Nothing has really been changed yet in the data itself.

# We can do a simple plot to remind ourselves of this animal's data
plot(WB.track.paita)

# amt has a function to summarize sampling rate information from the resulting track object.
summarize_sampling_rate(WB.track.paita)

# Here we see the median sampling rate is 1 hour, with a minimum of .02 hrs and a max of 211 hours.
# Visualize the # of locations per day, to get a sense of the consistency of the sampling over time.

WB.track.paita %>%
  group_by(date = as_date(t_)) %>%
  summarize(count = n()) %>%
  ggplot(aes(x = date,
             y = count)) +
  geom_col() +
  ylim(c(0,20))

# Add the speed (straight line travel) and net squared displacement to dataframe
WB.track.paita <- WB.track.paita %>%
  add_nsd() %>% 
  mutate(speed = speed(WB.track.paita))

plot(WB.track.paita$t_, sqrt(WB.track.paita$nsd_)/1000, typ="l", xlab = "Tracking Period", ylab = "NSD (km)")

ggplot(data = WB.track.paita,
       aes(x = t_,
           y = sqrt(nsd_)/1000)) +
  geom_line() +
  ggtitle("Paita - displacement (km) over time") +
  ylab("NSD (km)") +
  xlab("Time")

# This highlights the large movements made from the initial point of capture. We can see 7 or 8 relatively large movements, and few smaller scale excursions as well. But by and large this is looking like a very range resident individual, because none of the displacements continue over any significant time periods. 


# AMT: Resample track and summarize steps, 1 animal ----------------

# Here we will resample the track to 3 hr intervals to be able to have a consistent measure of certain movement metrics, like step length, over the entire 24hr period. We could also resample to 1 hour, but that would sacrifice any data during the night, since no 1hr intervals exist during that period. 

# Here we are using a tolerance of 20 min for considering something to be at the 3 hour interval. There is no hard rule for this decision. We will first resample the track, then create a file based on steps.

WB.track.paita.3h <- track_resample(WB.track.paita,
                                       rate = hours(3),
                                       tolerance = minutes(20)) 

# We can see we've gone from 110108 locations to 5426 locations in the resampled track. This is because we've lost a ton of hourly locations from the daytime sampling.


# Convert to Steps -------------------------------------------------

# The next step is to convert the file to a step file. This will give us step lengths and turning angles. We are specifying that we are only keeping bursts (strings of locations) that have at least 3 consecutive locations, the min required to get a turning angle. The steps_by_burst changes the object from one representing points, to one representing steps. 

# I'm also adding time of day, which is a nice built in function, calling another package called "suncalc". You may need to install this package if this code does not work for you. I am not including the crepuscular period here, because it is a relatively short period of the day and with 3 hour intervals we are not likely to sample it well.

WB.stps.paita <- WB.track.paita.3h %>% 
  filter_min_n_burst(min_n = 3) %>%
  steps_by_burst() %>% 
  time_of_day(include.crepuscule = F)

class(WB.stps.paita)

WB.stps.paita %>%
  ggplot(aes(x = sl_)) +
  geom_density(alpha = 0.4) +
  labs(x = "Step length [m]",
       y = "Density") +
  theme_light()

# And we can view these values in a standard data frame format as well. 
WB.stps.paita
names(WB.stps.paita)

# we can also now look at step lengths vs. time of day
WB.stps.paita %>%
  ggplot(aes(x = sl_,
             col = tod_end_ )) +
  geom_density(alpha = 0.4) +
  labs(x = "Step length [m]",
       y = "Density") +
  theme_light()

# There does not appear to be much difference here for step lengths for Paita between the day and night.


# amt: Resample track and summarize steps, all animals -------------

# We can now repeat this workflow but now for all animals together. Here I'm adding nsd in the same code chain that I make the track, and I'm also specifying that we have an id column, so each animal can be kept separate.

WB.tracks <-  WB.data %>%
  make_track(.x = X,
             .y = Y,
             .t = timestamp,
             id = animal_id,
             crs = UTMZone.proj,
             all_cols = TRUE) %>%
  add_nsd()

# This create a huge dataframe of all the animals together
# In the code below I'm indicating that I would like to nest all the data except the id column. 

WB.tracks2 <- WB.tracks %>%  
  nest(data = -"animal_id")

WB.tracks2

# So here we can see that we have an id column, and a data column, where everything else is stored. You can see that essentially each animal id has its own data frame. This operates a bit like the group_by() function, but is an efficient way of actually storing the data, and easily accessing it. As you see below, we can map across these individual dataframes using the purr::map() function, because the data column is a list.

# Run each section of this piped chain to see what happens at each step.

WB.tracks2.3hr <- WB.tracks2  %>%
  mutate(track_rs = map(data, function(x)
    x  %>%  track_resample(rate = hours(3),
                           tolerance = minutes(20)))) %>% 
  dplyr::select(animal_id, track_rs) %>% 
  unnest(cols = track_rs) %>% 
  write_rds("Data/wildebeest_3hr.rds")

# Now we'll start over, and will this time continue to the step of creating bursts for each animal and creating steps. We also have to specify here we want to keep our attribute columns when we convert to steps. Here I'm specifying I want the attribute associated with the end of the step. This would be relevant for example for the point error information, but not for something like sex of the animal.

WB.stps.all <- WB.tracks2  %>%
  mutate(steps = map(data, function(x)
    x  %>%  track_resample(rate = hours(3),
                           tolerance = minutes(20))  %>%
      filter_min_n_burst(min_n = 3) %>%
      steps_by_burst(keep_cols = 'end')))

# Let's have a look at what we may have lost by selecting this sampling rate. This is showing the number of step lengths that were ultimately calculated for each day of the study period. If your animals had periods of sampling at different intervals for example, this would show what time periods were lost from the resampling decision. Here the animals were tracked for different lengths of time, but we have not created any gaps in sampling and we generally have 8 steps per day, which is what we would expect.

WB.stps.all  %>%
  select(animal_id, steps) %>%
  unnest(cols = steps) %>%
  group_by(date = as_date(t2_),
           animal_id) %>%
  summarize(count = n()) %>%
  ggplot(aes(x = date,
             y = count)) +
  geom_col() +
  facet_wrap(~ animal_id) +
  labs(y = "steps per day")

# Finally, we can select name and steps, unnest the new data_frame and create a plot of the step-length distributions.

WB.stps.all  %>%
  select(animal_id, steps) %>%
  unnest(cols = steps) %>%
  ggplot(aes(x = sl_,
             fill = animal_id)) +
  geom_density(alpha = 0.4) +
  labs(x = "Step length [m]",
       y = "Density") +
  xlim(c(0, 5000)) +
  theme_light()

# Even more interesting would be to compare males and females here 

WB.stps.all  %>%
  select(animal_id, steps) %>%
  unnest(cols = steps) %>%
  ggplot(aes(x = sl_,
             fill = sex)) +
  geom_density(alpha = 0.4) +
  labs(x = "Step length [m]",
       y = "Density") +
  xlim(c(0, 2000)) +
  theme_light()


# ADEHABITAT: ------------------------------------------------------

# An animal track in the adehabitatLT package is referred to as a trajectory. Type II trajectories are those for which each location is associated with a recorded time, so that's what we'll use. We'll need to input into this function a dataframe of xy coordinates, and so we'll use the same object we created for amt above, which now has updated projected coordinates. 

# We use the as.ltraj function to create individual animal trajectories and use the infolocs command to include attributes of the data that we want to maintain. The "slsp" argument indicates how to deal with turning angles when successive locations are in the same place. See the help menu for details on this. One big difference from amt is that step lengths and turning angles are created straightaway in this first step.

XY <- WB.data %>% 
  dplyr::select(X:Y) %>% 
  as.data.frame()

WB.traj.raw <-as.ltraj(WB.data %>% 
                          dplyr::select(X:Y) %>% 
                          as.data.frame(),
                        date = WB.data$timestamp,
                        id = WB.data$animal_id, 
                        typeII = TRUE, 
                        infolocs = WB.data[ ,3:6],
                        slsp = "missing")

WB.traj.raw <-as.ltraj(xy = XY,
                       date = WB.data$timestamp,
                       id = WB.data$animal_id,
                       typeII = TRUE, 
                       infolocs = WB.data[ ,3:6],
                       slsp = "missing")

plot(WB.traj.raw[1])

# Let's get a sense of what we have created here
class(WB.traj.raw)

# So this is a list and a trajectory. Let's see what is in the list
length(WB.traj.raw)

# There are 12 elements in the list. These represent each wildebeest. So unlike amt, where everything is inside one large dataframe, we here have a list for each individual. 

# We can get a nice little summary by just running the name of the trajectory object
WB.traj.raw

# We can pull them out individually by indexing.
WB.traj.raw[[1]]

# And we can use indexing for a plot as well. Note in the help file for plot.ltraj that you can modify the output here somewhat.
plot(WB.traj.raw)

# If you have a lot of animals plotting one a time probably makes more sense. 
plot(WB.traj.raw,
     id = c("Kikaya", "Karbolo", "Peria", "Sotua"))

# And we can check to see what our columns are in this new data object that we have for each animal.
names(WB.traj.raw[[1]])
head(WB.traj.raw[[1]])

# So we can see what we have here now. dt: time between locations in seconds; dist: distance between locations, here this is the distance to the NEXT location; R2n: net squared displacement; abs.angle: absolute turning angle; rel.angle: relative turning angle. dx and dy represent the change in the x and y directions.

# Note that there is an additional plotting option in this package to plot various metrics over time, like step distance. 
plotltr(WB.traj.raw[5],
        which = "dist")

# Of course at this point keep in mind that there may be irregular time intervals and that will influence the metrics. And we can look at those times to see how regular they are.

plotltr(WB.traj.raw[5],
        which = "dt")

plotltr(WB.traj.raw[5],
        which = "DOP")


# ADEHABITAT: resampling the data ------------------------------------

# To this point, the function has no idea what our sampling interval was. This is something we must set. These metrics then that have been calculated are raw metrics, and are not necessarily useful, since the time interval between locations may not be constant. And in fact we know it is not. In our case, the data were collected: Every hour from 6 am to 6 pm and every three hours from 6 pm to 6 am. We can actually pull sampling interval data from the trajectory now, since we have times between locations. There is no set function for this as there is for amt.

# Let's look at the median time between locations for each animal in the dataset. Here I'm using a purr::map() to iterate this calculation over each list item. To get the name, unfortunately this is only stored as an attribute of each list, and so it's a bit more tedious to extract. 

map_dfr(1:length(WB.traj.raw), ~
          WB.traj.raw[[.x]] %>% 
          summarize(name  = attr(.,
                                 which = "id"),
                    medianInt_sec = median(dt,
                                           na.rm = T),
                    medianInt_min = medianInt_sec/60,
                    minInt_min = min(dt/60, na.rm = T),
                    maxInt_days = max(dt/86400, na.rm = T)
          ))

# Here I am selecting the first relocation in the dataset, and generally this should work. We can see that this location is almost exactly on the hour.
refda.WB <- WB.traj.raw[[1]]$date[1]

# I could also create a date from scratch like this:
refda.WB2 <- with_tz(as_datetime("2010-10-19 10:00:00"), 
                      "Africa/Nairobi")

# Convert our trajectory into one date frame and get the min date. The ld() function is from adehabitatLT and converts our large list into a single data frame.
myref <- WB.traj.raw %>% 
  ld() %>% 
  pull(date) %>% 
  min() %>% 
  round("mins")


# In the setNA function, dt sets the time lag between successive relocations. tol sets the tolerance. This is the level of imprecision you are willing to put up with and still count a location. Here we are using 10 minutes.

WB.traj.na <- setNA(WB.traj.raw, 
                    date.ref = refda.WB, 
                    dt = 3, 
                    units = "hour") 

is.regular(WB.traj.na)

# Create a regular trajectory
# Round the times of the relocations to the exact time intervals.  

WB.traj.reg <- sett0(WB.traj.na, 
                      date.ref = refda.WB, 
                      dt = 3, 
                      units = "hour")

is.regular(WB.traj.reg)

# Convert to dataframe and Export
WB.traj.reg <- ld(WB.traj.reg)
write_rds(WB.traj.reg, file = "Data/wildebeest_3hr_adehabitat.rds")


# adehabitatLT: Summarize regular trajectory, get metrics ----------------------------------

# Let's now save the summary of our trajectory
Summary.WB.traj <- summary(WB.traj.reg)

# Let's look at this file. 
Summary.WB.traj 

# Note that nb.reloc is the regular expected number of relocations if the animal was sampled at 1hour intervals the whole 24 hours. So this is going to be in our case a number higher than what was actually collected. Let's modify this table now to add some additional useful information. The percentage complete is useful, especially if you have some animals that had trouble getting GPS positions, or if animals were sampled on different intervals, you'd want to see how much data was essentially lost by regularizing the sampling.

# Here I'm calculating the duration of the tracking period, as well as the actual # of records we have, and the percentage of the full sampling intervals that provided data.

# I'm also adding in the actual number of rows in the initial data frame for each animal. This is the actual # of GPS locations collected, after cleaning. I'm getting these values first, before adding them to the table, and I'm getting them from the file before we removed problematic data, since we want to see what we lost in this process.

og_locs <- WB.data %>% 
  as_tibble() %>% 
  summarize(original_locs = n(),
            .by = animal_id)

# Add details to the summary
Summary.WB.traj2 <- Summary.WB.traj %>% 
  mutate(
    Duration = round(difftime(date.end, 
                              date.begin, 
                              units = "days"), 
                     digits = 2),
    Records = nb.reloc - NAs,
    PctComplete = round((Records/nb.reloc)*100,
                         digits = 2),
  ) %>% 
  
  # bring in number of original locations
  left_join(og_locs,
            by = c("id" = "animal_id")) %>% 
  
  # calculate locations lost after regularization
  mutate(lost = original_locs - Records)

Summary.WB.traj2

# Take a look at this file. Most individuals have about a 60% completed records rate. This is essentially the result of the lost locations during the nighttime, since night was sampled at 3hr intervals. If you compare the total # of locations here, with the locations in our initial summary table, you can see how many locations were lost in this process. Unless your sampling interval was highly regular, with no locations collected at odd times, you will always lose some data here. Of course you'll lose more data if some animals had their sampling interval changed, or if some animals were sampled more frequently than others and you had to select a common, lower resolution interval across animals. Here, collections were hourly and were fairly consistent, so we didn't end up losing very much in this process. If we had regularized this to 3 hour intervals as above, we would have lost 2/3 of our daytime locations.

# Keep in mind also though although we do have points scattered during the night (every three hours), these will not be helpful in calculating statistics like step length and turning angle, so the percentage complete is a bit misleading here. We've essentially lost more data than is shown in the lost column, because we have no hourly data during the night (despite having some locations there) that will give us data on movement metrics. We'll see this when we look at the data below as a data frame.

# Let's now use the trajectory data to calculate some movement metrics. We will use the dedicated ld() function from the adehabitatLT package designed to convert the ltraj format to a data frame. We'll save this as a data frame object because we'll use it a few times below.

# Now that the data is at regular intervals, all these movement metrics are comparable and meaningful. It also lets us now look at filtering out potentially unrealistic movement based on speed. We'll get speed by dividing the distance moved in a step (dist) by the time interval, which here should be regular, at 1 hour, or 3600 seconds

WB.df <- WB.traj.reg %>% 
  ld() %>% 
  mutate(speed_m_s = dist/dt)

# Take a look at this dataframe and see all the NA values that have been inserted. So even though these are counted as data, you can see that no metrics like dist and rel.angle are provided.

# We can now look at the distribution of speeds for our animals, trying to identify outliers. Various approaches exist for this, but we can start by simply plotting the data for each animal with a boxplot, which identifies outliers just based on the data distribution.

WB.df %>% 
  filter(!is.na(speed_m_s)) %>% 
  ggplot(aes(x = id,
             y = speed_m_s)) +
  geom_boxplot()

# Let's look at the max speed in the dataset as well
max(WB.df$speed_m_s, na.rm = T)

# We can see across the animals that although there are outlier points in terms of the distribution, the maximum values are around 2 m/s. The outlier from Paita is not here, probably because this was from a 3 hour sampling interval at night that was lost in resampling to 1 hr. The new max speed range translates into about 4.5 mph or 7.2 kph. This is not an unrealistic speed of travel, especially not over 1 hour, and so we don't have any reason to exclude any points at this point. Packages exist to help to filter by speed, like SDLfilter, but ultimately these either require you to input a maximum allowable speed, or they use something like quantiles of the distribution of speeds, like we are doing here with our boxplots.

# We could also plot something like movement speed over time for an animal. Note that this could be compared to where the animal was at each step, for example comparing to distance to human development, distance to road, or habitat type. We'll dig into this a bit more in Week 3 when we attempt to associate movement metrics with animal behavior.

WB.df %>% 
  filter(id == "Karbolo") %>% 
  ggplot(aes(x = date,
             y = speed_m_s)) +
  geom_line()

# Let's now compile some useful stats into a table including average hourly movement distance, total movement distance, max speed, and maximum displacement distance for each animal. Note that we could have made a similar table with data from amt as well, I just decided to do it here once.

WB.moveStats <- WB.df %>% 
  summarize(avg_move_m = round(mean(dist,
                                    na.rm = TRUE),
                               digits = 2), 
            sum_move_km = round(sum(dist/1000,
                                    na.rm = TRUE),
                                digits = 2), 
            max_speed_m_s = round(max(speed_m_s,
                                      na.rm = TRUE),
                                  digits = 2),
            max_disp_km = round(max(sqrt(R2n)/1000,
                                    na.rm = TRUE),
                                digits = 2),
            .by = id)

WB.moveStats

# ADEHABITAT: Visualize Trajectory and Displacement -----------------------------------

# Here we are designing custom plots to look at the trajectory as well as the movement distances and displacement over time.

# Here, we will "loop" over every individuals (by name), creating a separate figure for each animal. We could also open a PDF device before running the loop, and make sure to print() the final map inside the map function to ensure it goes into the PDF. But for here we will just visualize them.

# Note the use of geom_path here, which operates well with a dataframe of X and Y locations that is plotted (not with a spatial sf file). It will connect points in the order that they appear in the dataset. This is why I am arranging the data by date within this map function.

# This is a relatively complicated chunk of code here and is put here as more of a demonstration than anything. But this shows an efficient way of visualizing all your animal tracks, as well as key metrics like step length and displacement over time. You can easily pull sections of this code to generate these plots individually for a single animal. 

# Note also the plots will look much better in the zoom window.

# Assign an animal to "i" to run and test individual sections of code.
# i = unique(gnu_data_fix$name)[1]

map(unique(WB.df$id),
    function (i) {
      
      # Remove the NAs and subset to a single animal. Remember NAs were added when we regularized the trajectory, Note that here "x" is referring to the x coordinate, not to some function value.
      
      my_data <- WB.df %>% 
        filter(id == i,
               !is.na(x))
      
      # Calculate the total days tracked
      
      duration <- my_data %>% 
        summarize(end = max(date),
                  start = min(date),
                  diff = end - start) %>% 
        pull(diff)
      
      
      # # Plot the trajectory
      
      traj_plot <- my_data %>% 
        arrange(date) %>% 
        ggplot(mapping = aes(x = x,
                             y= y)) +
        geom_point(shape = 16,
                   size = 0.5,
                   col = "blue") +
        geom_path(linewidth = 0.05) +
        geom_point(data = my_data %>% 
                     slice_min(date),
                   shape = 17,
                   col = "green",
                   size = 3) +
        geom_point(data = my_data %>% 
                     slice_max(date),
                   shape = 15,
                   col = "red",
                   size = 3) +
        
        labs(title = str_c("Movement of", 
                           i, 
                           sep = " "),
             subtitle = str_c(as_date(min(my_data$date)), 
                              "to", 
                              as_date(max(my_data$date)), 
                              sep = " ")) +
        theme_light()
      
      # Plot the movements over time (Velocity)
      step_plot <- my_data %>% 
        ggplot(aes(x = date,
                   y = dist/1000)) +
        geom_line() +
        labs(x = "Time",
             y = "Distance moved (km)",
             title = "Step lengths over time") +
        theme_light()
      
      
      # Plotting net displacement over steps
      nsd_plot <- my_data %>% 
        ggplot(aes(x = date,
                   y = sqrt(R2n)/1000)) +
        geom_line() +
        labs(x = "Time",
             y = "Distance moved (km)",
             title = "Net Displacement Over Time",
             subtitle = str_c(duration, "days",
                              sep = " ")) +
        theme_light()
      
      # Arrange the two data plots together in one, one column figure
      pg1 <- cowplot::plot_grid(step_plot,
                                nsd_plot,
                                ncol = 1,
                                rel_widths = c(1, 1)
      )
      
      # Arrange the two figure plot above with the main trajectory plot
      pg2 <- traj_plot
      
      cowplot::plot_grid(pg1, 
                         pg2, 
                         ncol = 2,
                         rel_widths = c(1, 1))
      
    }
)

# End Code