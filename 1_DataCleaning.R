## ----setup, include=FALSE----------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## ----Clean Libraries, message=FALSE, warning=FALSE---------------------------------------------------
# Remove from memory
rm(list=ls())

# You may need to install these packages first
#install.packages('svDialogs', 'tidyverse', 'move2', 'lubridate', 'tmap', 'sf', 'gt')

# Load required libraries
library(svDialogs)
library(tidyverse)
library(move2)
library(lubridate)
library(tmap)
library(sf)
library(gt)


## ----Clean Timezone, message=FALSE, warning=FALSE----------------------------------------------------
# Set TimeZone and UTM Zone
# Other timezones can be found at: https://en.wikipedia.org/wiki/List_of_tz_database_time_zones
Timezone1 <- 'UTC'
Timezone2 <- "Africa/Nairobi"
 
# UTM Zone
LatLong.proj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"  # EPSG:4326
#UtmZone.proj <- "+proj=utm +zone=37 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs" #This is EPSG:32737"
UtmZone.proj <- "EPSG:32737"


## ----Load, message=FALSE, warning=FALSE--------------------------------------------------------------
# We can upload each file directly and convert to a dataframe 
#WB <- as.data.frame(read_csv("./Data/White-bearded wildebeest (Connochaetes taurinus) movements - Kenya.csv"))

# Or list the files so we don't have to type so much
# Create list
lf <- list.files(path="./Data/", pattern = '.csv', all.files=FALSE, full.names=TRUE)
lf

WB <- as.data.frame(read_csv(lf[2]))
WB.ref <- as.data.frame(read_csv(lf[1]))

# Look at the data
# head(WB)
# head(WB.ref)

# Fix a few column name issues
names(WB) <- gsub("-", "_", names(WB))
names(WB) <- gsub(":", "_", names(WB))
names(WB.ref) <- gsub("-", "_", names(WB.ref))

# Check to identify the changes
# head(WB)
# head(WB.ref)

# We could also pull the file directory from Movebank.  To do so, you will need a Movebank UserName and Password.
# This is particularly useful when data continue to be collected (i.e., are actively streaming) on a study.
# Note: When pulling from Movebank, the data will contain a few more data fields than the uploaded CSV

# Set Movebank Login Details
# UN <- dlgInput("Enter Movebank UserName: ", Sys.info()[""])$res
# PW <- dlgInput("Enter Movebank Password: ", Sys.info()[""])$res

# Details
# login <- movebankLogin(username=UN, password=PW)
 
# Pull Data from Movebank and convert to a dataframe
# WB <- as.data.frame(getMovebankData(study = "White-bearded wildebeest (Connochaetes taurinus) movements - Kenya", login = login))
 
# Reference Data - No need to import here, as data are already subset in Movebank
# WB.ref <- getMovebankReferenceTable(study = "White-bearded wildebeest (Connochaetes taurinus) movements - Kenya", login = login)


## ----Verify, message=FALSE, warning=FALSE, results='hide', echo=FALSE--------------------------------
# View the dataset
head(WB)
tail(WB)
head(WB, 10)
WB[1:10,]

# What is the structure of the dataset?
# What is the data type of the time stamp column?
str(WB)
str(WB$timestamp)
str(WB.ref)

# How many rows and columns are there in the movement dataset?
dim(WB)
nrow(WB)
ncol(WB)

# How many tags?
sort(unique(WB$tag_local_identifier))
length(unique(WB$tag_local_identifier))

# How many study sites?  What are there names?
length(unique(WB.ref$study_site))
unique(WB.ref$study_site)

# What is the timezone of the dataset?
tz(WB$timestamp)


## ----Clean Merge, message=FALSE, warning=FALSE, echo=TRUE--------------------------------------------
# Clean the reference file, selecting only the columns that you want to include
WB <- 
  WB %>%

  # Pull in the WB.ref dataset to join the study site column
  left_join(WB.ref, by = join_by("tag_local_identifier" == "tag_id")) %>%

  # Rename columns and create local timestamp - UTC to local
  # Select only columns of interest
  transmute(id = tag_local_identifier,
            animal_id = individual_local_identifier,
            latitude = location_lat,
            longitude = location_long,
            sex = animal_sex,
            DOP = gps_dop,
            fixType = gps_fix_type_raw,
            temp = external_temperature,
            # Converting timestamp with 'with_tz' command
            timestamp = with_tz(timestamp, tz=Timezone2),
            deploy_on = with_tz(deploy_on_date, tz=Timezone2), # Field included in the reference dataset
            deploy_off = with_tz(deploy_off_date, tz=Timezone2), # Field included in the reference dataset
            study_site = study_site) %>% # Field included in the reference dataset

  # The 'across' function allows me to apply a function (as.factor) across multiple fields 
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
         # This is just an example here, as data have already been subset
         timestamp >= deploy_on & timestamp <= deploy_off) %>%

  # Remove fields that are now unnecessary.  Dropping the deployment dates and the study area
  dplyr::select(-c(deploy_on, deploy_off, study_site)) %>%
  
  # Remove extra levels (important since subsetting to a single study area)
  # The droplevels function re-assess what levels are in the data and drops the rest
  droplevels() %>%

  # Arrange the dataset by id and timestamp
  arrange(id, timestamp)

# Look again (yes again!) at your data
head(WB)


## ----Clean Verify, message=FALSE, warning=FALSE, results='hide', echo=FALSE--------------------------
# How many records
dim(WB)
nrow(WB)

# What is the structure of the dataset?
str(WB)

# How many tags?
sort(unique(WB$animal_id))
length(unique(WB$animal_id))

# What is the timezone of the dataset?
tz(WB$timestamp) 


## ----Summarize, message=FALSE, warning=FALSE, echo=TRUE----------------------------------------------
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

# Print Results
wb.Summary

# Now make prettier, saving the results to your Output folder 
gt_gnu <- wb.Summary %>% 
  
  # initialize gt table
  gt() %>%
  
  # Make the table easier to read with alternating grey bars
  opt_row_striping() %>%
  
  # Add title and subtitle, pulling date of creation
  tab_header(
    title = "White-bearded Wildebeest in Kenya: Tracking Data Summary",
    subtitle = Sys.Date()) %>%
  
  # Easy preset date formatting
  fmt_date(
    columns = c(Start, End),
    date_style = 8) %>%
  
  # Change the column labels for the table
  cols_label(animal_id = "Wildebeest ID",
             Sex = "Sex",
             Locations = "Total points",
             Start = "First location",
             End = "Last location",
             Duration = "Tracking period (days)") %>%
  
  # Center text inside columns
  cols_align(align = "center") 

# Print result
gt_gnu

# Save as html table to send to the project manager, or a shiny app
gtsave(gt_gnu, filename = "Output/summary_gnu.html")


## ----Visualize, message=FALSE, warning=FALSE, echo=TRUE----------------------------------------------
# Create very simple plot (non-spatial)
plot(WB$longitude, WB$latitude,
     col = WB$id,
     pch = 16,
     cex = 0.5,
     ylab = 'Northing',
     xlab = 'Easting',
     asp = 1)


## ----Visualize1, message=FALSE, warning=FALSE, echo=TRUE---------------------------------------------
# Convert
WB.sf <- WB %>% 
  st_as_sf(coords = c('longitude', 'latitude'), 
           crs = LatLong.proj) %>% 
  st_transform(UtmZone.proj)

# You could check the coordinate system by:
#st_crs(WB.sf)

# Look at the data
#head(WB.sf)
#str(WB.sf)
class(WB.sf)


## ----Visualize2, message=FALSE, warning=FALSE, echo=TRUE---------------------------------------------
# Plot using basic R function
# plot(WB.sf["animal_id"],
#      main = paste("Wildebeest: Athi-Kaputiei Plains ( n = ", length(unique(WB.sf$animal_id)),")"))

# GGPlot using the spatial object
WB.sf %>%
  ggplot() +
  geom_sf(aes(fill = animal_id),
          alpha = 0.6,
          shape = 21,
          col = "black") +
  scale_fill_discrete(name = "Animal ID") +
  ggtitle(paste("Wildebeest: Athi-Kaputiei Plains (n =", length(unique(WB.sf$animal_id)),")")) +
  coord_sf(datum = st_crs(UtmZone.proj)) + # Note, this line is necessary unless we want the data plotted in Lat/Long
  theme_minimal()

# or use the facet_wrap command to separate each individual into its own plot
# This is perhaps a bit more useful to look for potential errant points for each individual.
# WB.sf %>%
#   ggplot() +
#   geom_sf(aes(fill = animal_id),
#           alpha = 0.6,
#           shape = 21,
#           col = "black") +
#   scale_fill_discrete(name = "Animal ID") +
#   theme_minimal() +
#   facet_wrap(~ animal_id)


## ----Visualize3, message=FALSE, warning=FALSE, echo=TRUE---------------------------------------------
tmap_mode("view")
 
# Select from a range of basemaps:.
# https://leaflet-extras.github.io/leaflet-providers/preview/
# You can also type "providers$" in the console to see all the options, but this won't allow you to easily preview them.

# We'll pull a world satellite map and a world street map.
# Note that the format of tmap requires us to refer to variables in quotes

Athi.Map <-
  tm_basemap(c("Esri.WorldImagery",
             "OpenStreetMap")) +
  tm_shape(WB.sf,
           #name = "Wildebeest Locations (Athi-Kaputiei)") +
           name = paste("Wildebeest: Athi-Kaputiei Plains (n =", length(unique(WB.sf$animal_id)),")")) +
  tm_dots(size = 0.025,
          title = "Wildebeest ID",
          col = "animal_id",
          alpha = 0.5)

# Create the map
Athi.Map

# Save Output
tmap_save(Athi.Map, filename = "Output/Athi_GnuMap.html")


## ----Save, message=FALSE, warning=FALSE, echo=TRUE---------------------------------------------------
# Save individual file
#write_rds(WB, file = "Data/wildebeest_WBonly.rds")
# or use saveRDS(object, file = "my_data.rds")

# Save both files together
save(WB, WB.sf, file = "Data/wildebeest_data.rdata")

