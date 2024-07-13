## ----setup, include=FALSE----------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----Setup, message=FALSE, warning=FALSE, echo=TRUE--------------------------------------------------
# Remove items from memory/clean your workspace
rm(list=ls())

# You may need to install these packages first
#install.packages('sf', 'adehabitatLT', 'gridExtra', 'tidyverse')

# Load libraries
library(sf)
library(adehabitatLT)
library(gridExtra)
library(tidyverse)


## ----Timezone, message=FALSE, warning=FALSE, echo=TRUE-----------------------------------------------
# TimeZone
Timezone1 <- 'UTC'
Timezone2 <- "Africa/Nairobi"
 
# Spatial Reference
LatLong.proj <- "EPSG:4326"
UtmZone.proj <- "EPSG:32737"


## ----Data Import, message=FALSE, warning=FALSE, echo=TRUE--------------------------------------------
# Data load/import
load("./Data/wildebeest_data.rdata")

# View your environment or print objects to screen using ls()
# ls()

# Check the projection.
# st_crs(WB.sf)

# AdeHabitatLT requires a dataframe for analyses.
# Let's use the WB.sf dataset because the data have been projected and we'll need the projected coordinates for analyses.  Convert to a 'flat' dataframe.

# Dataframe conversion
WB.data <- WB.sf %>%
  as_tibble() %>%
  mutate(X = st_coordinates(WB.sf)[ ,1],
         Y = st_coordinates(WB.sf)[ ,2]) %>%
  dplyr::select(-geometry)

# Look at the data
head(WB.data)


## ----DOP, message=FALSE, warning=FALSE, echo=TRUE----------------------------------------------------
# Plot the DOP values for confirmation.  
# For this dataset, we will separate 2D and 3D positions and then use a qualitative filter to remove data above a threshold.
# Here I am being more restrictive on 2D positions (dop < 5.0) than 3D positions (dop < 10.0).

# Let's first get a summary of how many records are in the dataset before removing records.  This way we can track how much this filtering impacts the size of the dataset.
val1 <- nrow(WB.data)

# Plot 2D positions
P1.FT2 <- 
  ggplot(WB.data[WB.data$fixType %in% "2",], aes(x = DOP)) +
  geom_histogram(color = "black", fill = "white", bins = 25) +
  labs(title = "GPS Data (2D)", x = "DOP", y = "Frequency") +
  geom_vline(xintercept = 5, color = "red", linetype = "dotted", linewidth = 1) +
  theme_classic()

# Plot 3D positions (most of data)
P1.FT3 <- 
  ggplot(WB.data[WB.data$fixType %in% "3",], aes(x = DOP)) +
  geom_histogram(color = "black", fill = "white", bins = 25) +
  labs(title = "GPS Data (3D)", x = "DOP", y = "Frequency") +
  geom_vline(xintercept = 10, color = "red", linetype = "dotted", linewidth = 1) +
  theme_classic()

# Filter/subset
WB.data <- 
  WB.data %>% 
  filter(
    fixType == 3 & DOP < 10 | fixType == 2 & DOP < 5) # Only accept position with a 3D fixtype and DOP less than 10 or a 2D fixtype and DOP less than 5

# How many records now after filtering?
val2 <- nrow(WB.data)

# Plot again
P2.FT2 <- 
  ggplot(WB.data[WB.data$fixType %in% "2",], aes(x = DOP)) +
  geom_histogram(color = "black", fill = "white", bins = 25) +
  labs(title = "GPS Data (2D) - DOP Filtered", x = "DOP", y = "Frequency") +
  geom_vline(xintercept = 5, color = "red", linetype = "dotted", linewidth = 1) +
  theme_classic()

P2.FT3 <- 
  ggplot(WB.data[WB.data$fixType %in% "3",], aes(x = DOP)) +
  geom_histogram(color = "black", fill = "white", bins = 25) +
  labs(title = "GPS Data (3D) - DOP Filtered", x = "DOP", y = "Frequency") +
  geom_vline(xintercept = 10, color = "red", linetype = "dotted", linewidth = 1) +
  theme_classic()

# Plot the data together
grid.arrange(P1.FT2, P2.FT2, P1.FT3, P2.FT3, ncol = 2)

# What's the percent of data that have been removed?
round((val1-val2)/val1, digits = 4)

# Only 333 records (<1% of data) were removed based on the criteria we used.


## ----Trajectory, message=FALSE, warning=FALSE, echo=TRUE---------------------------------------------
# Extract the X and Y coordinates to input into as.ltraj
XY <- WB.data %>% 
  dplyr::select(X:Y) %>% 
  as.data.frame()

# Create trajectory
WB.traj <-as.ltraj(xy = XY,
                       date = WB.data$timestamp,
                       id = WB.data$animal_id,
                       typeII = TRUE, 
                       infolocs = WB.data[ ,3:6],
                       slsp = "missing")

# Look at the summary of the created object
WB.traj


## ----Resample, message=FALSE, warning=FALSE, echo=TRUE-----------------------------------------------
# Here, we will use the first location in the dataset as the reference date/time.  
# Set reference date/time
refda <- WB.traj[[1]]$date[1]

# Set all null values to NA.  Based on 3 hour interval.
WB.traj.na <- setNA(WB.traj, 
                    date.ref = refda, 
                    dt = 3, 
                    units = "hour") 

# Summarize the new dataset
# WB.traj.na

# Is the trajectory regular?
is.regular(WB.traj.na)

# Create a regular trajectory by rounding the time to the exact time intervals.
WB.traj.reg <- sett0(WB.traj.na, 
                      date.ref = refda, 
                      dt = 3, 
                      units = "hour")

# Summarize the new dataset
WB.traj.reg

# Is the trajectory regular?
is.regular(WB.traj.reg)


## ----Summarize, message=FALSE, warning=FALSE, echo=TRUE----------------------------------------------
# Summarize the movement trajectory
Summary.traj <- summary(WB.traj.reg)

# Note that nb.reloc is the total number of relocations collected at a 3 Hour sampling interval.  Any missing data have been filled with NA, allowing us to calculate the percent complete. 

# Add details to the summary table
Summary.traj <-
  Summary.traj %>% 
  mutate(
    Duration = round(difftime(date.end, 
                              date.begin,
                              units = "days"),
                     digits = 2),
    Records = nb.reloc - NAs,
    PctComplete = round((Records/nb.reloc)*100,
                        digits = 2))

# View the table
Summary.traj


## ----Plotting, message=FALSE, warning=FALSE, echo=TRUE-----------------------------------------------
# Plot all animals together, helpful since all animals are plotted on the same scale.
plot(WB.traj.reg)

# View each animal separately using vector notation
plot(WB.traj.reg[5])

# We can also specify the name of the animals we want to plot.  Here, I've included 4 unique animals to plot.
# unique(WB.data$animal_id)
# plot(WB.traj,
#      id = c("Kikaya", "Karbolo", "Peria", "Sotua"))

# View the columns included in the data object
#names(WB.traj.reg[[1]])
#head(WB.traj.reg[[1]])

# Contents of the ltraj object:
# dt: time between locations in seconds
# dist: distance between the next location
# R2n: net squared displacement
# abs.angle: absolute turning angle
# rel.angle: relative turning angle
# dx and dy represent the change in the x and y directions.

# A variety of additional plotting options also exist in the package.  See help(plotltr).
#plotltr(WB.traj.reg[5], which = "dist") # Distance moved between discrete points (3 hour interval)
#plotltr(WB.traj.reg[5], which = "dt/60") # This shows that the data are regular.  y-axis corrected (dt/60) to show minutes
#plotltr(WB.traj.reg[5], which = "DOP") # DOP.  As expected, all values are < 10 DOP.


## ----Plotting2, message=FALSE, warning=FALSE, echo=TRUE----------------------------------------------
# Here I will use basic R plotting functions, but this same workflow could adopted using GGPLot (A good homework assignment!)

# Convert trajectory to a flat dataframe (all data appended into a single file)
WB.move <- ld(WB.traj.reg)

# Create a reference id so can change easily between individuals
Id.val <- unique(WB.move$id)
 
# Determine which animal you want to plot
i <- 5

# Account for the NAs in the dataset a subset the data to a single individual
WB.sub <- subset(WB.move[!is.na(WB.move$x),], id == Id.val[i]) 

# Setup plot layout with three panels
layout(matrix(c(1,1,2,3), 2, 2, byrow = FALSE), widths=1, heights=c(1,1))

# Calculate the total days tracked
time.diff <- trunc(difftime(WB.sub$date[nrow(WB.sub)],WB.sub$date[1]),units="days")

# Plot the trajectory
plot(WB.sub$x,WB.sub$y,typ="l",xlab="Easting",ylab="Northing",
     main=paste(WB.sub$id[1]," Movement"),frame=FALSE,axes=FALSE,asp=1)
     mtext(paste(format(WB.sub$date[1],"%Y-%m-%d")," to ",format(WB.sub$date[nrow(WB.sub)],"%Y-%m-%d")),cex=0.75) # Just specifying how I want the dates to be reported
     axis(1, labels=TRUE)
     axis(2, labels=TRUE)
     # Color the points
     points(WB.sub$x,WB.sub$y,pch=16,cex=0.5,col="blue")  # All points Blue
     points(WB.sub$x[1],WB.sub$y[1],pch=17,cex=1,col="green") # Starting point Green
     points(WB.sub$x[nrow(WB.sub)],WB.sub$y[nrow(WB.sub)],pch=15,cex=1,col="red") # End point Red

# Plot the movements over time (Velocity)
plot(WB.sub$date, WB.sub$dist/1000, type='l', ylab="Distance moved (km)", xlab="Time", main="Steplengths", frame=FALSE)
    # Calculate the time from release date
    mtext(paste(abs(time.diff)," days"),cex=0.75)

# Plot the net displacement per step
plot(WB.sub$date, sqrt(WB.sub$R2n)/1000, type='l', ylab="Distance (km)", xlab="Days Since Release", main="Net Displacement",frame=FALSE)
    mtext(paste(abs(time.diff)," days"),cex=0.75)


## ----Plotting Loop, message=FALSE, warning=FALSE, echo=FALSE-----------------------------------------
# Loop over each individual (by id) and save result
# Create ID.val
Id.val <- unique(WB.move$id)
 
for (i in 1:length(Id.val)){
  
# Account for the NAs in the dataset a subset the data to a single individual
WB.sub <- subset(WB.move[!is.na(WB.move$x),], id == Id.val[i]) 

# Create the output name of the plot
png(filename = paste0("Output/Plots/",WB.sub$id[i],"_MvmtPlot.png"))
    
# Setup plot layout with three panels
layout(matrix(c(1,1,2,3), 2, 2, byrow = FALSE), widths=1, heights=c(1,1))

# Calculate the total days tracked
time.diff <- trunc(difftime(WB.sub$date[nrow(WB.sub)],WB.sub$date[1]),units="days")

# Plot the trajectory
plot(WB.sub$x,WB.sub$y,typ="l",xlab="Easting",ylab="Northing",
     main=paste(WB.sub$id[1]," Movement"),frame=FALSE,axes=FALSE,asp=1)
     mtext(paste(format(WB.sub$date[1],"%Y-%m-%d")," to ",format(WB.sub$date[nrow(WB.sub)],"%Y-%m-%d")),cex=0.75) # Just specifying how I want the dates to be reported
     axis(1, labels=TRUE)
     axis(2, labels=TRUE)
     # Color the points
     points(WB.sub$x,WB.sub$y,pch=16,cex=0.5,col="blue")  # All points Blue
     points(WB.sub$x[1],WB.sub$y[1],pch=17,cex=1,col="green") # Starting point Green
     points(WB.sub$x[nrow(WB.sub)],WB.sub$y[nrow(WB.sub)],pch=15,cex=1,col="red") # End point Red

# Plot the movements over time (Velocity)
plot(WB.sub$date, WB.sub$dist/1000, type='l', ylab="Distance moved (km)", xlab="Time", main="Steplengths", frame=FALSE)
    # Calculate the time from release date
    mtext(paste(abs(time.diff)," days"),cex=0.75)

# Plot the net displacement per step
plot(WB.sub$date, sqrt(WB.sub$R2n)/1000, type='l', ylab="Distance (km)", xlab="Days Since Release", main="Net Displacement",frame=FALSE)
    mtext(paste(abs(time.diff)," days"),cex=0.75)

dev.off()
}


## ----Metrics, message=FALSE, warning=FALSE, echo=TRUE------------------------------------------------
# First, let's calculate speed, since distance traveled is provided with change in time (dt).  We should also make id and sex factorsr
WB.move <- WB.move %>% 
  mutate(speed = dist/dt) %>%
  mutate(across(c(id,sex), as.factor))

# Create a boxplot of the speeds traveled for each individual 
# Note that I'm filtering the data so that only data that is not NA is included.
# Note also that the units are meters per second
WB.move %>% 
  filter(!is.na(speed)) %>% 
  ggplot(aes(x = id,
             y = speed)) +
  geom_boxplot()

# Create a graph of distance traveled
WB.move %>%
  filter(!is.na(dist)) %>%
  ggplot(aes(x = dist/1000)) +
  geom_histogram(color="black", fill="white", bins = 50) +
  labs(x = "Distance (km)") +
  theme_classic()

# Or we might want the same information, but showing the individual animals
WB.move %>% 
  filter(!is.na(dist)) %>%
  ggplot(aes(x = dist/1000, fill = id)) +
  geom_density(alpha = 0.3) + # alpha determines transparency 
  labs(x = "Distance (km)") +
  theme_classic()

# What is the maximum speed traveled
max(WB.move$speed, na.rm = TRUE)

# What about distance
max(WB.move$dist/1000, na.rm = TRUE) # So a distance of 13.4 km in a 3 hour period.

# Calculate statistics for each animal
WB.moveStats <- WB.move %>% 
  #group_by(id) %>% # Grouping here would give you the same result
  summarize(AvgMove = round(mean(dist/1000, na.rm = TRUE), digits = 2),
            SumMove = round(sum(dist/1000, na.rm = TRUE), digits = 2),
            MaxSpeed = round(max(speed, na.rm = TRUE), digits = 2),
            MaxDisp = round(max(sqrt(R2n)/1000, na.rm = TRUE), digits = 2),
            .by = id)

WB.moveStats


## ----Save, message=FALSE, warning=FALSE, echo=TRUE---------------------------------------------------
write_rds(WB.move, file = "Data/wildebeest_3hr_adehabitat.rds")

