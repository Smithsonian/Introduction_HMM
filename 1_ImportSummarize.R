# ******************************************************************
# ******************************************************************

# Project: Intro to Animal Mvmts
# Description: Read Wildebeest Data
# Author: Jared Stabach
# Date: 13 June 2024

# ******************************************************************
# ******************************************************************

# Remove anything in memory
rm(list=ls())

# Load Packages
library(tidyverse)
library(lubridate)
library(move2)
library(svDialogs)
library(sf)
library(tmap)
library(gt)

# Settings ---------------------------------------------------------


# Set TimeZone and UTM Zone
Timezone1 <- 'UTC'
Timezone2 <- "Africa/Nairobi"
 
# UTM Zone
LatLong.proj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"  # EPSG:4326
UTMZone.proj <- "+proj=utm +south +zone=37 +units=m +datum=WGS84 +init=EPSG:32737"


# ***************************************
# ***************************************

# CSV Download

# ***************************************
# ***************************************

# Create list (so no need to type filename
lf <- list.files(path="./Data/", pattern = '.csv', all.files=FALSE, full.names=TRUE)

WB <- as.data.frame(read_csv(lf[2]))
WB.ref <- as.data.frame(read_csv(lf[1]))

# Fix a few column name issues that result during import
names(WB) <- gsub("-", "_", names(WB))
names(WB) <- gsub(":", "_", names(WB))
names(WB.ref) <- gsub("-", "_", names(WB.ref))


# ***************************************
# ***************************************

# Movebank Download
# Movebank direct download contains more data fields than the uploaded CSV

# ***************************************
# ***************************************

# # Set Movebank Login Details 
# UN <- dlgInput("Enter Movebank UserName: ", Sys.info()[""])$res
# PW <- dlgInput("Enter Movebank Password: ", Sys.info()[""])$res
# 
# # Details
# login <- movebankLogin(username=UN, password=PW)
# 
# # Pull Data from Movebank and convert to a dataframe
# WB <- as.data.frame(getMovebankData(study = "White-bearded wildebeest (Connochaetes taurinus) movements - Kenya", login = login))
# 
# # Reference Data - No need to import here, as data are already subset in Movebank
# WB.ref <- getMovebankReferenceTable(study = "White-bearded wildebeest (Connochaetes taurinus) movements - Kenya", login = login)

# ***************************************
# ***************************************


# Check 
# What tags are included in the study?  How many?
sort(unique(WB$tag_local_identifier))
length(unique(WB$tag_local_identifier))

# Look at Data
head(WB)

# Check the structure of each dataframe
str(WB)
str(WB.ref)

# Check the timezone attribute
tz(WB$timestamp)

# Data Cleaning ----------------------------------------------------

# Create new columns, remove duplicates, filter, and arrange
# Join the reference dataset to include the study area

WB <- WB %>% 
  
  # Pull in the WB.ref dataset to join the study site column
  left_join(WB.ref, by = join_by("tag_local_identifier" == "tag_id")) %>% 
  
  # Rename columns and create timestamp - UTC to local
  # Select only columns of interest
  transmute(id = tag_local_identifier,
            animal_id = individual_local_identifier,
            latitude = location_lat,
            longitude = location_long,
            sex = animal_sex,
            DOP = gps_dop,
            fixType = gps_fix_type_raw,
            temp = external_temperature,
            timestamp = with_tz(timestamp, tz=Timezone2),
            deploy_on = with_tz(deploy_on_date, tz=Timezone2),
            deploy_off = with_tz(deploy_off_date, tz=Timezone2),
            study_site = study_site) %>%
  
  # Correct the data type for specific fields
  mutate(across(c(id,sex,study_site), as.factor)) %>% 
  
  # Make sure no duplicate id and timestamp exist.
  distinct(animal_id, timestamp, .keep_all = TRUE) %>% 
  
  # Remove any records that don't have a timestamp or a Lat/Long location
  filter(!is.na(timestamp),
         !is.na(latitude),
         !is.na(longitude),
         latitude != 0, 
         longitude != 0,
         
         # Grab only the Athi-Kaputiei Data
         study_site == "Athi-Kaputiei Plains",
         
         # And use the deploy on and off dates to further subset
         # (Won't do anything here, as dates already subset)
         timestamp >= deploy_on & timestamp <= deploy_off) %>%
  
  # Remove fields that are now unnecessary.  Dropping the deployment dates and the study area
  dplyr::select(-c(deploy_on, deploy_off, study_site)) %>% 
  
  # Remove extra levels (important since subsetting to a single study area)
  # The droplevels function reassesses what levels are in the data and drops the rest.
  droplevels() %>% 
  
  # Arrange the dataset by id and timestamp
  arrange(id, timestamp)

# Look at the new, condensed dataset
head(WB)
str(WB)

# Now what tags are included?  How many?
sort(unique(WB$id))
length(unique(WB$id))

# For reference, these data were collected:
unique(WB.ref$duty_cycle)


# Visualize --------------------------------------------------------------

# Create sf object
WB.sf <- WB %>% 
  st_as_sf(coords = c('longitude', 'latitude'),
           crs = LatLong.proj) %>% 
  st_transform(UTMZone.proj)

# Plot the data
plot(WB.sf["animal_id"], 
     main = paste0("Wildebeest: Athi-Kaputiei Plains (n = ", length(unique(WB.sf$id)),")"))

# We can do better than that, but at least this confirms what I'd expect
WB.sf %>% 
  ggplot() +
  geom_sf(aes(fill = animal_id),
          alpha = 0.6,
          shape = 21,
          col = "black") +
  theme_bw()

# or separate each individual into its own plot
WB.sf %>% 
  ggplot() +
  geom_sf(aes(fill = animal_id),
          alpha = 0.6,
          shape = 21,
          col = "black") +
  theme_minimal() +
  facet_wrap(~ animal_id)

# This is perhaps the most helpful view for getting a sense of the data collected for each individual, how much their space use varies, and whether there are errant points that may represent actual GPS error. Here fortunately there are no points that can be identified visually as significant outliers. Be sure to use the "zoom" button above the plot to see it in larger size.

# Let's do the same thing using tmap, overlaying the data on a satellite map to investigate potential problems
# tmap camps can be made in "view" or "interactive" mode, or "plot" mode.  Plot mode is better for making stating maps

tmap_mode("view")

# Here we can select a range of basemaps. And you can find a list of available maps HERE. 
# https://leaflet-extras.github.io/leaflet-providers/preview/

# You can also type "providers$" in the console to see all the options, but this won't allow you to easily preview them.

# Here I'm pulling a world satellite map, and a world street map. 
# Note also that in tmap we need to refer to our variables in quotes.

tm_basemap(c("Esri.WorldImagery",
             "OpenStreetMap.HOT")) +
  tm_shape(WB.sf, 
           name = "Wildebeest Locations (Athi-Kaputiei)") +
  tm_dots(size = 0.05,
          title = "Wildebeest ID",
          col = "animal_id",
          alpha = 0.5) 

# With this map, you can zoom in and out to investigate points from a particular animal. You can also toggle on and off the various basemap tiles you have selected using the tile icon in the top left corner. If you plotted other data as well, you can also toggle these layers on and off. 

# Spend some time zooming into the map. You can see a large area in the central north of the study region with no locations. Once you zoom in you can see that this area is developed and that the wildebeest do not venture into this part of the landscape. There very well may be fencing excluding the herds from these areas as well. In the extreme north of the study area is the Kenyan capital Nairobi. This is more obvious when clicking on the street map basemap.

# Again here though, I do not see anything visually that would signal large positional errors. Let's dig into error a bit more in the next section.


# Summarize --------------------------------------------------------------

# Create summary object
wb.Summary <- WB %>% 
  
  summarize(
    Locations = n(),
    Sex = unique(sex),
    Start = min(timestamp),
    End = max(timestamp),
    Duration = round(End - Start, digits = 1),
    .by = animal_id) %>% 
  
  # Arrange results
  arrange(animal_id, Start, desc(Locations))

wb.Summary


# We can make this a little prettier by creating a gt formatted table. 
gt_gnu <- wb.Summary %>% 
  
  # initialize gt table
  gt() %>%
  
  # Make the table easier to read with alternating grey bars
  opt_row_striping() %>%
  
  # add title and subtitle, pulling date of creation
  tab_header(
    title = "White-bearded Wildebeest in Kenya: Tracking Data Summary",
    subtitle = Sys.Date()) %>%
  
  # easy preset date formatting
  fmt_date(
    columns = c(Start, End),
    date_style = 8) %>%
  
  # changing column labels for the table
  cols_label(animal_id = "Wildebeest ID",
             Sex = "Sex",
             Locations = "Total points",
             Start = "First location",
             End = "Last location",
             Duration = "Tracking period (days)") %>%
  
  # center text inside columns
  cols_align(align = "center") 

# Print result
gt_gnu

# Save as html table to send to the project manager, or a shiny app
gtsave(gt_gnu, filename = "summary_gnu_kenya.html")


# Save File --------------------------------------------------------------

#write_rds(WB, file = "Data/wildebeest_data.rds")
#WB <- read_rds(file = "Data/wildebeest_data.rds")

#write_rds(WB.sf, file = "Data/wildebeest_data_sf.rds")

save(WB, WB.sf, file = "Data/wildebeest_data.rdata")
#load("Data/wildebeest_data.rdata")


# End Code