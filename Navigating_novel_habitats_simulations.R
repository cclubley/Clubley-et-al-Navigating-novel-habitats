# The following script carries out the simulations, analysis, and presentation of results (Figures 6 & 7) from the correlated random walk models described and presented in "Navigating novel habitats: assessing the impacts of artificial complexity on space-use by a key intertidal grazer" by Clubley et al. 

library(dplyr)
library(plotrix)
library(raster)
library(SiMRiv)
library(viridis)
library(ggplot2)
library(scales)

# Calculate a resistance value' for each panel type ----------------------

# SiMRiv uses 'resistance' values to define the ability of an organism to move over a certain area. We need to assign a resistance value to each of the four panel types to define the ability of a limpet to move over that panel. We will use the gross distances, calculated from the lab experiment, to inform the resistance values

# Load the data from the laboratory experiments
limpet <- read.csv("./Master_limpet_movement.csv")

# Remove unneeded variables to clean data, and then re-classify the variables
limpet <- limpet[, c(2, 8, 10, 11, 15)]
limpet$Replicate <- as.factor(limpet$Replicate)
limpet$Limpet_ID <- as.character(limpet$Limpet_ID)
limpet$Panel_type <- as.factor(limpet$Panel_type)

# Calculate the gross distance travelled by each individual limpet
limpet <- limpet %>% group_by(Replicate, Limpet_ID, Limpet_size, Panel_type) %>% 
  summarise(Gross_distance=sum(Distance_moved_cm))
# Average by panel type
limpet_avg <- limpet %>% group_by(Panel_type) %>% 
  summarise(mean=mean(Gross_distance))
limpet_avg
#         mean   
# Flat    31.0  
# Ripple  74.1  
# 2.5 cm  43.0  
# 5 cm    31.4 

# Inverse normalise the means:
# We calculate 1 - the normalisation (inversion), because the metric used by SiMRiv is 'resistance', meaning greater resistance produces less movement. Therefore, the panel with the least amount of movement (Flat), should have the highest resistance metric score 
limpet_avg$metric <- 1-(limpet_avg$mean-min(limpet_avg$mean))/(max(limpet_avg$mean)-min(limpet_avg$mean))

# Now normalise the raw gross distances
limpet$normgross <- 1-(limpet$Gross_distance-min(limpet$Gross_distance))/(max(limpet$Gross_distance)-min(limpet$Gross_distance))
# And calculate the standard error of the normalised values so that we can produce a range of resistance values for each panel type
limpet_se <- limpet %>% group_by(Panel_type) %>% 
  summarise(std.error(normgross))

# Calculate the lower and upper bounds for the mean for each tile type
limpet_avg$lower <- limpet_avg$metric-limpet_se$`std.error(normgross)`
limpet_avg$upper <- limpet_avg$metric+limpet_se$`std.error(normgross)`
limpet_avg

# RESISTANCE VALUES FOR EACH PANEL:

#          lower conf.      metric   upper conf.   
# Flat           0.951           1             1             
# Ripple             0           0        0.0545         
# 2.5 cm         0.687       0.721         0.756         
# 5 cm           0.958       0.992             1    

# Clear the environment
rm(limpet, limpet_avg, limpet_se)

# -------------------------------------------------------------------------
# ------------------------------- Control ---------------------------------
# -------------------------------------------------------------------------
# Create a 8 x 16 grid ----------------------------------------------------

# Run a control replicate that uses only the resistance value for the flat panels

# Create a 8 x 16 rectangular grid with a resolution of 0.25 (each grid cell is 0.25 m x 0.25 m - the same size as the concrete panels)
# Define the total extent to match 0.25m x 0.25m cells over 16 x 8 layout
grid <- raster(xmn=0, xmx=4, ymn=0, ymx=2, ncol=16, nrow=8)
res(grid) # Confirm the resolution = 0.25
# Convert to polygons
grid <- rasterToPolygons(grid)
# Plot it to check
plot(grid)

# Assign each polygon within the grid a unique ID
polygonID <- paste0(1:length(grid))
for (i in 1:length(slot(grid, "polygons"))) {
  slot(slot(grid, "polygons")[[i]], "ID") <- polygonID[i]
}

# Add the IDs to the data slot in the grid
slot(grid, "data") <- cbind("ID" = 1:length(grid), slot(grid, "data"))

# Clean the environment
rm(polygonID, i)

# Add resistance values ---------------------------------------------------

# Give the grid polygons resistance attributes based on the inverse normalised average distances for each panel

# Metric values are randomly sampled from within the range of values and assigned to each grid square
grid@data$resistance <- runif(length(grid@data$ID), 0.958, 1)

# Create a resistance raster ----------------------------------------------

# Save the resistance values in the environment
grid_val <- grid@data$resistance

# Convert to a data frame with no duplicates
mapvalues <- as.data.frame(cbind(grid_val))
mapvalues <- data.frame(x = unlist(mapvalues))
row.names(mapvalues) <- NULL
mapvalues$x2 <- as.character(mapvalues$x)
mapvalues <- mapvalues[!duplicated(mapvalues),]

# Create a resistance raster from the grid created earlier
# Include a margin of 2 m to stop random walks from getting trapped in the corner of the grid
resistance <- resistanceFromShape(grid, res=0.25, field="resistance", background=1,
                                  mapvalues=with(mapvalues, setNames(x, x)), margin=2)    
# plot to check
plot(resistance)

# Clear the environment
rm(grid_val, mapvalues)

# Simulate limpet movement ------------------------------------------------

# Based on limpet densities in the field and the size of the model area, we are modelling 1000 limpets 

# Create the limpets
# Max step length set to 0.44 (11 cm in raster resolution) - based on lab observations
# CRW factor set to 0.98 based on prior optimisation
species.list <- lapply(1:1000, function(i) species(state.CRW(0.98)+0.44))

# Run the simulation
set.seed(123)  # for reproducibility
tracks <- simulate(species.list,
                   time = 288, # 24 hours in 5 minute time steps
                   coords = xyFromCell(resistance, sample(1:ncell(resistance), size=1000, replace=TRUE)), # Sample 1000 random cell indices
                   resist=resistance,
                   start.resistance=1)

# Save path length --------------------------------------------------------

# Create an empty list
traj_list <- vector("list", 1000)

# Loop over limpets
for (i in seq_len(1000)) {
  # Compute column indices for each limpet
  cols <- ((i-1)*3+1):(i*3)
  
  # Extract limpet data as a matrix
  limpet_data <- tracks[, cols]
  
  # Apply function
  traj_list[[i]] <- sampleMovement(limpet_data, resolution=1, resist=resistance)
}

all_steplengths <- list()
all_resistance <- list()
limpet_id <- list()

# Loop through results and collect the data
for (i in seq_along(traj_list)) {
  n <- length(traj_list[[i]][["stats"]][["steplengths"]])
  all_steplengths[[i]] <- traj_list[[i]][["stats"]][["steplengths"]]
  all_resistance[[i]] <- traj_list[[i]][["stats"]][["resistance"]]
  limpet_id[[i]] <- rep(i, n)
}

# Combine into vectors
steplengths_vec <- unlist(all_steplengths)
resistance_vec <- unlist(all_resistance)
limpet_vec <- factor(unlist(limpet_id))

# Create data frame
df <- data.frame(
  steplength = steplengths_vec,
  resistance = resistance_vec,
  limpet = limpet_vec
)

# Group the data by limpet and calculate the path length
steplength <- df %>%
  group_by(limpet) %>%
  summarise(path_length=sum(steplength, na.rm=TRUE))

# # Save the data
# write.csv(steplength, "./Control_stats.csv", row.names=FALSE)

# Track plot --------------------------------------------------------------

# Remove every third column (state) as only one state was used in the simulation
tracks <- tracks[, -seq(3, ncol(tracks), by=3)]

# Create a list of data frames, one per individual
tracks <- lapply(1:1000, function(i) {
  df <- data.frame(
    x = tracks[, 2*i - 1],
    y = tracks[, 2*i]
  )
  return(df)
})

# Plot to check
plot(grid, col="azure2")
cols <- magma(1000)
for (i in 1:1000){
  lines(tracks[[i]], col=cols[i], lty=2, lwd=1)
}

# Plot a heat map ---------------------------------------------------------

# Convert the list of tracks to a data frame
tracks_df <- as.data.frame(do.call(rbind, tracks))

# Plot
ggplot()+
  stat_density2d(tracks_df, mapping=aes(x=x, y=y, fill=..density..), geom='raster', contour=FALSE)+
  scale_fill_viridis(option="magma", direction=-1, limits=c(0, 0.5), na.value="white", expand=FALSE, oob=squish) +
  geom_polygon(grid, mapping=aes(x=long, y=lat, group=group), colour="white", fill=NA)

# Save the R environment --------------------------------------------------

# # Save the environment
# save.image(file='./Control.RData')
# load('./Control.RData')

# -------------------------------------------------------------------------
# -------------------------------- Grouped --------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# ----------------------- ~14% coverage, Ripple panel ---------------------
# -------------------------------------------------------------------------
# Create a 8 x 16 grid ----------------------------------------------------

grid <- raster(xmn=0, xmx=4, ymn=0, ymx=2, ncol=16, nrow=8)
grid <- rasterToPolygons(grid)
plot(grid)

ID <- paste0(1:length(grid))
for (i in 1:length(slot(grid, "polygons"))) {
  slot(slot(grid, "polygons")[[i]], "ID") <- ID[i]
}
slot(grid, "data") <- cbind("ID" = 1:length(grid), slot(grid, "data"))

# Sample the grid to create tile sets -------------------------------------

# The first area of coverage is ~14%, which is equal to 18 tiles (14% * 128). 

# Select a random sample of 18 polygons that are next to each other 
# Select one polygon to begin with
gridsample1 <- grid[sample(length(ID), 1),] 

# Set number to 1 for loop
number <- 1
# NB: Sometimes you have to re-run the loop to get grid cells that are next to each other - I am sure there is some way to optimise this...
while (number < 2) { # whilst the number is still 1...
  
  grid <- grid       # re-set the grid
  sample <- grid[sample(length(ID), 1),] # sample 1 grid polygon 
  
  for (i in 1:length(gridsample1)){    # for the length of the grid sample...
    
    if (sample@data[["ID"]] == gridsample1@data[["ID"]][i]){  # if the polygon ID is already present...
      print("Same ID")  # print "same ID"
      number <- number  # number remains at 1
    } 
    
    else if (sample@data[["ID"]] == gridsample1@data[["ID"]][i]+1 |    # if the polygon ID borders an existing polygon ID...
             sample@data[["ID"]] == gridsample1@data[["ID"]][i]-1 |
             sample@data[["ID"]] == gridsample1@data[["ID"]][i]+16 |
             sample@data[["ID"]] == gridsample1@data[["ID"]][i]-16) {
      
      gridsample1 <- rbind(gridsample1, sample) # add the polygon to the grid sample
      gridsample1 <- gridsample1[!duplicated(gridsample1@data),]  # remove any duplicates (was happening in the loop for some reason)
      
    } 
    if (length(gridsample1) == 18) {  # if there are 18 polygons in the sample...
      number <- 2 # number = 2 and the loop is stopped
    }
  } 
}

# Assign the same projection as the original grid
projection(gridsample1) <- projection(grid)
# Remove the sampled polygons form the original grid
grid <- grid-gridsample1
# Plot to check
plot(grid)
plot(gridsample1, col="lightgreen", add=T)

# The remaining polygons become grid_sample2, which will be the flat tiles
gridsample2 <- grid
projection(gridsample2) <- projection(grid)

# Assign each grid sample to a tile type 
flat <- gridsample2 # Grid sample 2 has to be flat
ripple <- gridsample1

# Plot to check
plot(flat, col="white")
plot(ripple, col='lightgreen', add=TRUE)

# Add resistance values ---------------------------------------------------

flat@data$resistance <- runif(length(flat@data$ID), 0.95096189, 1)
ripple@data$resistance <- runif(length(ripple@data$ID), 0, 0.05446582)

crs(ripple) <- crs(flat)
grid_14_ripp <- rbind(flat, ripple)

# Create a resistance raster ----------------------------------------------

flat_val <- flat@data$resistance
ripple_val <- ripple@data$resistance

mapvalues <- as.data.frame(cbind(flat_val, ripple_val))
mapvalues <- data.frame(x = unlist(mapvalues))
row.names(mapvalues) <- NULL
mapvalues$x2 <- as.character(mapvalues$x)
mapvalues <- mapvalues[!duplicated(mapvalues),]

resistance_14_ripp <- resistanceFromShape(grid_14_ripp,    
                                          res=0.25,     
                                          field="resistance",   
                                          mapvalues=with(mapvalues, setNames(x, x2)),    
                                          background=1, margin=2)       

plot(resistance_14_ripp)

# Simulate limpet movement ------------------------------------------------

species.list <- lapply(1:1000, function(i) species(state.CRW(0.98)+0.44))

set.seed(123)
tracks_14_ripp <- simulate(species.list,
                   time = 288, 
                   coords = xyFromCell(resistance_14_ripp, sample(1:ncell(resistance_14_ripp), size=1000, replace=TRUE)),
                   resist=resistance_14_ripp,
                   start.resistance=1)

# Save path length --------------------------------------------------------

traj_list <- vector("list", 1000)

for (i in seq_len(1000)) {
  cols <- ((i-1)*3+1):(i*3)
  limpet_data <- tracks_14_ripp[, cols]
  traj_list[[i]] <- sampleMovement(limpet_data, resolution=1, resist=resistance_14_ripp)
}

all_steplengths <- list()
all_resistance <- list()
limpet_id <- list()

for (i in seq_along(traj_list)) {
  n <- length(traj_list[[i]][["stats"]][["steplengths"]])
  all_steplengths[[i]] <- traj_list[[i]][["stats"]][["steplengths"]]
  all_resistance[[i]] <- traj_list[[i]][["stats"]][["resistance"]]
  limpet_id[[i]] <- rep(i, n)
}

steplengths_vec <- unlist(all_steplengths)
resistance_vec <- unlist(all_resistance)
limpet_vec <- factor(unlist(limpet_id))

df <- data.frame(
  steplength = steplengths_vec,
  resistance = resistance_vec,
  limpet = limpet_vec
)

steplength <- df %>%
  group_by(limpet) %>%
  summarise(path_length=sum(steplength, na.rm=TRUE))

# Track plot --------------------------------------------------------------

tracks_14_ripp <- tracks_14_ripp[, -seq(3, ncol(tracks_14_ripp), by=3)]

tracks_14_ripp <- lapply(1:1000, function(i) {
  df <- data.frame(
    x = tracks_14_ripp[, 2*i - 1],
    y = tracks_14_ripp[, 2*i]
  )
  return(df)
})

plot(flat, col="azure2")
plot(ripple, col="lightyellow", add=TRUE)
cols <- magma(1000)
for (i in 1:1000){
  lines(tracks_14_ripp[[i]], col=cols[i], lty=2, lwd=1)
}

# Plot a heat map ---------------------------------------------------------

tracks_14_ripp_df <- as.data.frame(do.call(rbind, tracks_14_ripp))

ggplot()+
  stat_density2d(tracks_14_ripp_df, mapping=aes(x=x, y=y, fill=..density..), geom='raster', contour=FALSE)+
  scale_fill_viridis(option="magma", direction=-1, limits=c(0, 0.5), na.value="white", expand=FALSE, oob=squish) +
  geom_polygon(grid, mapping=aes(x=long, y=lat, group=group), colour="white", fill=NA)

# Save the R environment --------------------------------------------------

# save.image(file='./R environments/Ripple_Grouped_14percent.RData')
# load('./R environments/Ripple_Grouped_14percent.RData')

# -------------------------------------------------------------------------
# ----------------------- ~14% coverage, 2.5 cm high ----------------------
# -------------------------------------------------------------------------
# Create a 8 x 16 grid ----------------------------------------------------

grid <- raster(xmn=0, xmx=4, ymn=0, ymx=2, ncol=16, nrow=8)
grid <- rasterToPolygons(grid)
plot(grid)

ID <- paste0(1:length(grid))
for (i in 1:length(slot(grid, "polygons"))) {
  slot(slot(grid, "polygons")[[i]], "ID") <- ID[i]
}
slot(grid, "data") <- cbind("ID" = 1:length(grid), slot(grid, "data"))

# Sample the grid to create tile sets -------------------------------------

gridsample1 <- grid[sample(length(ID), 1),] 
number <- 1
while (number < 2) { 
  grid <- grid      
  sample <- grid[sample(length(ID), 1),] 
  for (i in 1:length(gridsample1)){    
    if (sample@data[["ID"]] == gridsample1@data[["ID"]][i]){  
      print("Same ID")  
      number <- number  
    } 
    else if (sample@data[["ID"]] == gridsample1@data[["ID"]][i]+1 |   
             sample@data[["ID"]] == gridsample1@data[["ID"]][i]-1 |
             sample@data[["ID"]] == gridsample1@data[["ID"]][i]+16 |
             sample@data[["ID"]] == gridsample1@data[["ID"]][i]-16) {
      gridsample1 <- rbind(gridsample1, sample) 
      gridsample1 <- gridsample1[!duplicated(gridsample1@data),]  
    } 
    if (length(gridsample1) == 18) {  
      number <- 2 
    }
  } 
}

projection(gridsample1) <- projection(grid)
grid <- grid-gridsample1
plot(grid)
plot(gridsample1, col="lightblue", add=T)

gridsample2 <- grid
projection(gridsample2) <- projection(grid)

flat <- gridsample2 
cm_25 <- gridsample1

plot(flat, col="white")
plot(cm_25, col='lightblue', add=TRUE)

# Add resistance values ---------------------------------------------------

flat@data$resistance <- runif(length(flat@data$ID), 0.95096189, 1)
cm_25@data$resistance <- runif(length(cm_25@data$ID), 0.68657211, 0.75614646)

crs(cm_25) <- crs(flat)
grid_14_cm_25 <- rbind(flat, cm_25)

# Create a resistance raster ----------------------------------------------

flat_val <- flat@data$resistance
cm_25_val <- cm_25@data$resistance

mapvalues <- as.data.frame(cbind(flat_val, cm_25_val))
mapvalues <- data.frame(x = unlist(mapvalues))
row.names(mapvalues) <- NULL
mapvalues$x2 <- as.character(mapvalues$x)
mapvalues <- mapvalues[!duplicated(mapvalues),]

resistance_14_cm_25 <- resistanceFromShape(grid_14_cm_25,    
                                           res=0.25,      
                                           field="resistance",
                                           mapvalues=with(mapvalues, setNames(x, x2)),     
                                           background=1, margin=2)       
plot(resistance_14_cm_25)

# Simulate limpet movement ------------------------------------------------

species.list <- lapply(1:1000, function(i) species(state.CRW(0.98)+0.44))

set.seed(123)  
tracks_14_cm25 <- simulate(species.list,
                           time = 288,
                           coords = xyFromCell(resistance_14_cm_25, sample(1:ncell(resistance_14_cm_25), size=1000, replace=TRUE)),
                           resist=resistance_14_cm_25,
                           start.resistance=1)

# Save path length --------------------------------------------------------

traj_list <- vector("list", 1000)

for (i in seq_len(1000)) {
  cols <- ((i-1)*3+1):(i*3)
  limpet_data <- tracks_14_cm25[, cols]
  traj_list[[i]] <- sampleMovement(limpet_data, resolution=1, resist=resistance_14_cm_25)
}

all_steplengths <- list()
all_resistance <- list()
limpet_id <- list()

for (i in seq_along(traj_list)) {
  n <- length(traj_list[[i]][["stats"]][["steplengths"]])
  all_steplengths[[i]] <- traj_list[[i]][["stats"]][["steplengths"]]
  all_resistance[[i]] <- traj_list[[i]][["stats"]][["resistance"]]
  limpet_id[[i]] <- rep(i, n)
}

steplengths_vec <- unlist(all_steplengths)
resistance_vec <- unlist(all_resistance)
limpet_vec <- factor(unlist(limpet_id))

df <- data.frame(
  steplength = steplengths_vec,
  resistance = resistance_vec,
  limpet = limpet_vec
)

steplength <- df %>%
  group_by(limpet) %>%
  summarise(path_length=sum(steplength, na.rm=TRUE))

# Track plot --------------------------------------------------------------

tracks_14_cm25 <- tracks_14_cm25[, -seq(3, ncol(tracks_14_cm25), by=3)]

tracks_14_cm25 <- lapply(1:1000, function(i) {
  df <- data.frame(
    x = tracks_14_cm25[, 2*i - 1],
    y = tracks_14_cm25[, 2*i]
  )
  return(df)
})

plot(flat, col="azure2")
plot(cm_25, col="lightyellow", add=TRUE)
cols <- magma(1000)
for (i in 1:1000){
  lines(tracks_14_cm25[[i]], col=cols[i], lty=2, lwd=1)
}

# Plot a heat map ---------------------------------------------------------

tracks_14_cm25_df <- as.data.frame(do.call(rbind, tracks_14_cm25))

ggplot()+
  stat_density2d(tracks_14_cm25_df, mapping=aes(x=x, y=y, fill=..density..), geom='raster', contour=FALSE)+
  scale_fill_viridis(option="magma", direction=-1, limits=c(0, 0.5), na.value="white", expand=FALSE, oob=squish) +
  geom_polygon(grid, mapping=aes(x=long, y=lat, group=group), colour="white", fill=NA)

# Save the R environment --------------------------------------------------

# save.image(file='./R environments/cm25_Grouped_14percent.RData')
# load('./R environments/cm25_Grouped_14percent.RData')

# -------------------------------------------------------------------------
# ------------------------ ~14% coverage, 5 cm high -----------------------
# -------------------------------------------------------------------------
# Create a 8 x 16 grid ----------------------------------------------------

grid <- raster(xmn=0, xmx=4, ymn=0, ymx=2, ncol=16, nrow=8)
grid <- rasterToPolygons(grid)
plot(grid)

ID <- paste0(1:length(grid))
for (i in 1:length(slot(grid, "polygons"))) {
  slot(slot(grid, "polygons")[[i]], "ID") <- ID[i]
}
slot(grid, "data") <- cbind("ID" = 1:length(grid), slot(grid, "data"))

# Sample the grid to create tile sets -------------------------------------

gridsample1 <- grid[sample(length(ID), 1),] 
number <- 1
while (number < 2) { 
  grid <- grid       
  sample <- grid[sample(length(ID), 1),] 
  for (i in 1:length(gridsample1)){    
    if (sample@data[["ID"]] == gridsample1@data[["ID"]][i]){  
      print("Same ID")  
      number <- number  
    } 
    else if (sample@data[["ID"]] == gridsample1@data[["ID"]][i]+1 |    
             sample@data[["ID"]] == gridsample1@data[["ID"]][i]-1 |
             sample@data[["ID"]] == gridsample1@data[["ID"]][i]+16 |
             sample@data[["ID"]] == gridsample1@data[["ID"]][i]-16) {
      gridsample1 <- rbind(gridsample1, sample) 
      gridsample1 <- gridsample1[!duplicated(gridsample1@data),]  
      
    } 
    if (length(gridsample1) == 18) {  
      number <- 2 
    }
  } 
}

projection(gridsample1) <- projection(grid)
grid <- grid-gridsample1
plot(grid)
plot(gridsample1, col="orange", add=T)

gridsample2 <- grid
projection(gridsample2) <- projection(grid)

flat <- gridsample2 
cm_5 <- gridsample1

plot(flat, col="white")
plot(cm_5, col='orange', add=TRUE)

# Add resistance values ---------------------------------------------------

flat@data$resistance <- runif(length(flat@data$ID), 0.95096189, 1)
cm_5@data$resistance <- runif(length(cm_5@data$ID), 0.95783566, 1) 

crs(cm_5) <- crs(flat)
grid_14_cm_5 <- rbind(flat, cm_5)

# Create a resistance raster ----------------------------------------------

flat_val <- flat@data$resistance
cm_5_val <- cm_5@data$resistance

mapvalues <- as.data.frame(cbind(flat_val, cm_5_val))
mapvalues <- data.frame(x = unlist(mapvalues))
row.names(mapvalues) <- NULL
mapvalues$x2 <- as.character(mapvalues$x)
mapvalues <- mapvalues[!duplicated(mapvalues),]

resistance_14_cm_5 <- resistanceFromShape(grid_14_cm_5,    
                                          res=0.25,      
                                          field="resistance",   
                                          mapvalues=with(mapvalues, setNames(x, x2)),     
                                          background=1, margin=2)     
plot(resistance_14_cm_5)

# Simulate limpet movement ------------------------------------------------

species.list <- lapply(1:1000, function(i) species(state.CRW(0.98)+0.44))

set.seed(123) 
tracks_14_cm5 <- simulate(species.list,
                           time = 288, 
                           coords = xyFromCell(resistance_14_cm_5, sample(1:ncell(resistance_14_cm_5), size=1000, replace=TRUE)),
                           resist=resistance_14_cm_5,
                           start.resistance=1)

# Save path length --------------------------------------------------------

traj_list <- vector("list", 1000)

for (i in seq_len(1000)) {
  cols <- ((i-1)*3+1):(i*3)
  limpet_data <- tracks_14_cm5[, cols]
  traj_list[[i]] <- sampleMovement(limpet_data, resolution=1, resist=resistance_14_cm_5)
}

all_steplengths <- list()
all_resistance <- list()
limpet_id <- list()

for (i in seq_along(traj_list)) {
  n <- length(traj_list[[i]][["stats"]][["steplengths"]])
  all_steplengths[[i]] <- traj_list[[i]][["stats"]][["steplengths"]]
  all_resistance[[i]] <- traj_list[[i]][["stats"]][["resistance"]]
  limpet_id[[i]] <- rep(i, n)
}

steplengths_vec <- unlist(all_steplengths)
resistance_vec <- unlist(all_resistance)
limpet_vec <- factor(unlist(limpet_id))

df <- data.frame(
  steplength = steplengths_vec,
  resistance = resistance_vec,
  limpet = limpet_vec
)

steplength <- df %>%
  group_by(limpet) %>%
  summarise(path_length=sum(steplength, na.rm=TRUE))

# Track plot --------------------------------------------------------------

tracks_14_cm5 <- tracks_14_cm5[, -seq(3, ncol(tracks_14_cm5), by=3)]

tracks_14_cm5 <- lapply(1:1000, function(i) {
  df <- data.frame(
    x = tracks_14_cm5[, 2*i - 1],
    y = tracks_14_cm5[, 2*i]
  )
  return(df)
})

plot(flat, col="azure2")
plot(cm_5, col="lightyellow", add=TRUE)
cols <- magma(1000)
for (i in 1:1000){
  lines(tracks_14_cm5[[i]], col=cols[i], lty=2, lwd=1)
}

# Plot a heat map ---------------------------------------------------------

tracks_14_cm5_df <- as.data.frame(do.call(rbind, tracks_14_cm5))

ggplot()+
  stat_density2d(tracks_14_cm5_df, mapping=aes(x=x, y=y, fill=..density..), geom='raster', contour=FALSE)+
  scale_fill_viridis(option="magma", direction=-1, limits=c(0, 0.5), na.value="white", expand=FALSE, oob=squish) +
  geom_polygon(grid, mapping=aes(x=long, y=lat, group=group), colour="white", fill=NA)

# Save the R environment --------------------------------------------------

# save.image(file='./R environments/cm5_Grouped_14percent.RData')
# load('./R environments/cm5_Grouped_14percent.RData')

# -------------------------------------------------------------------------
# ----------------------- ~33% coverage, Ripple panel ---------------------
# -------------------------------------------------------------------------
# Create a 8 x 16 grid ----------------------------------------------------

grid <- raster(xmn=0, xmx=4, ymn=0, ymx=2, ncol=16, nrow=8)
grid <- rasterToPolygons(grid)
plot(grid)

ID <- paste0(1:length(grid))
for (i in 1:length(slot(grid, "polygons"))) {
  slot(slot(grid, "polygons")[[i]], "ID") <- ID[i]
}
slot(grid, "data") <- cbind("ID" = 1:length(grid), slot(grid, "data"))

# Sample the grid to create tile sets -------------------------------------

# The second area of coverage is ~33%, which is equal to 42 tiles (33% * 128). 

# Select a random sample of 42 polygons that are next to each other 
# NB: Sometimes you have to re-run the loop to get grid cells that are next to each other - I am sure there is some way to optimise this...
gridsample1 <- grid[sample(length(ID), 1),] 
number <- 1
while (number < 2) { 
  grid <- grid       
  sample <- grid[sample(length(ID), 1),] 
  for (i in 1:length(gridsample1)){    
    if (sample@data[["ID"]] == gridsample1@data[["ID"]][i]){  
      print("Same ID") 
      number <- number 
    } 
    else if (sample@data[["ID"]] == gridsample1@data[["ID"]][i]+1 |    
             sample@data[["ID"]] == gridsample1@data[["ID"]][i]-1 |
             sample@data[["ID"]] == gridsample1@data[["ID"]][i]+16 |
             sample@data[["ID"]] == gridsample1@data[["ID"]][i]-16) {
      gridsample1 <- rbind(gridsample1, sample) 
      gridsample1 <- gridsample1[!duplicated(gridsample1@data),]  
    } 
    if (length(gridsample1) == 42) {  
      number <- 2 
    }
  } 
}

projection(gridsample1) <- projection(grid)
grid <- grid-gridsample1
plot(grid)
plot(gridsample1, col="lightgreen", add=T)

gridsample2 <- grid
projection(gridsample2) <- projection(grid)

flat <- gridsample2 
ripple <- gridsample1

plot(flat, col="white")
plot(ripple, col='lightgreen', add=TRUE)

# Add resistance values ---------------------------------------------------

flat@data$resistance <- runif(length(flat@data$ID), 0.95096189, 1)
ripple@data$resistance <- runif(length(ripple@data$ID), 0, 0.05446582)

crs(ripple) <- crs(flat) 
grid_33_ripp <- rbind(flat, ripple)

# Create a resistance raster ----------------------------------------------

flat_val <- flat@data$resistance
ripple_val <- ripple@data$resistance

mapvalues <- as.data.frame(cbind(flat_val, ripple_val))
mapvalues <- data.frame(x = unlist(mapvalues))
row.names(mapvalues) <- NULL
mapvalues$x2 <- as.character(mapvalues$x)
mapvalues <- mapvalues[!duplicated(mapvalues),]

resistance_33_ripp <- resistanceFromShape(grid_33_ripp,   
                                          res=0.25,     
                                          field="resistance",   
                                          mapvalues=with(mapvalues, setNames(x, x2)),     
                                          background=1, margin=2)       
plot(resistance_33_ripp)

# Simulate limpet movement ------------------------------------------------

species.list <- lapply(1:1000, function(i) species(state.CRW(0.98)+0.44))

set.seed(123)  
tracks_33_ripp <- simulate(species.list,
                          time = 288, 
                          coords = xyFromCell(resistance_33_ripp, sample(1:ncell(resistance_33_ripp), size=1000, replace=TRUE)),
                          resist=resistance_33_ripp,
                          start.resistance=1)

# Save path length --------------------------------------------------------

traj_list <- vector("list", 1000)

for (i in seq_len(1000)) {
  cols <- ((i-1)*3+1):(i*3)
  limpet_data <- tracks_33_ripp[, cols]
  traj_list[[i]] <- sampleMovement(limpet_data, resolution=1, resist=resistance_33_ripp)
}

all_steplengths <- list()
all_resistance <- list()
limpet_id <- list()

for (i in seq_along(traj_list)) {
  n <- length(traj_list[[i]][["stats"]][["steplengths"]])
  all_steplengths[[i]] <- traj_list[[i]][["stats"]][["steplengths"]]
  all_resistance[[i]] <- traj_list[[i]][["stats"]][["resistance"]]
  limpet_id[[i]] <- rep(i, n)
}

steplengths_vec <- unlist(all_steplengths)
resistance_vec <- unlist(all_resistance)
limpet_vec <- factor(unlist(limpet_id))

df <- data.frame(
  steplength = steplengths_vec,
  resistance = resistance_vec,
  limpet = limpet_vec
)

steplength <- df %>%
  group_by(limpet) %>%
  summarise(path_length=sum(steplength, na.rm=TRUE))

# Track plot --------------------------------------------------------------

tracks_33_ripp <- tracks_33_ripp[, -seq(3, ncol(tracks_33_ripp), by=3)]

tracks_33_ripp <- lapply(1:1000, function(i) {
  df <- data.frame(
    x = tracks_33_ripp[, 2*i - 1],
    y = tracks_33_ripp[, 2*i]
  )
  return(df)
})

plot(flat, col="azure2")
plot(ripple, col="lightyellow", add=TRUE)
cols <- magma(1000)
for (i in 1:1000){
  lines(tracks_33_ripp[[i]], col=cols[i], lty=2, lwd=1)
}

# Plot a heat map ---------------------------------------------------------

tracks_33_ripp_df <- as.data.frame(do.call(rbind, tracks_33_ripp))

ggplot()+
  stat_density2d(tracks_33_ripp_df, mapping=aes(x=x, y=y, fill=..density..), geom='raster', contour=FALSE)+
  scale_fill_viridis(option="magma", direction=-1, limits=c(0, 0.5), na.value="white", expand=FALSE, oob=squish) +
  geom_polygon(grid, mapping=aes(x=long, y=lat, group=group), colour="white", fill=NA)

# Save the R environment --------------------------------------------------

# save.image(file='./R environments/Ripple_Grouped_33percent.RData')
# load('./R environments/Ripple_Grouped_33percent.RData')

# -------------------------------------------------------------------------
# ----------------------- ~33% coverage, 2.5 cm high ----------------------
# -------------------------------------------------------------------------
# Create a 8 x 16 grid ----------------------------------------------------

grid <- raster(xmn=0, xmx=4, ymn=0, ymx=2, ncol=16, nrow=8)
grid <- rasterToPolygons(grid)
plot(grid)

ID <- paste0(1:length(grid))
for (i in 1:length(slot(grid, "polygons"))) {
  slot(slot(grid, "polygons")[[i]], "ID") <- ID[i]
}
slot(grid, "data") <- cbind("ID" = 1:length(grid), slot(grid, "data"))

# Sample the grid to create tile sets -------------------------------------

gridsample1 <- grid[sample(length(ID), 1),] 
number <- 1
while (number < 2) { 
  grid <- grid       
  sample <- grid[sample(length(ID), 1),] 
  for (i in 1:length(gridsample1)){    
    
    if (sample@data[["ID"]] == gridsample1@data[["ID"]][i]){  
      print("Same ID")  
      number <- number  
    } 
    else if (sample@data[["ID"]] == gridsample1@data[["ID"]][i]+1 |    
             sample@data[["ID"]] == gridsample1@data[["ID"]][i]-1 |
             sample@data[["ID"]] == gridsample1@data[["ID"]][i]+16 |
             sample@data[["ID"]] == gridsample1@data[["ID"]][i]-16) {
      gridsample1 <- rbind(gridsample1, sample) 
      gridsample1 <- gridsample1[!duplicated(gridsample1@data),]  
    } 
    if (length(gridsample1) == 42) {  
      number <- 2 
    }
  } 
}

projection(gridsample1) <- projection(grid)
grid <- grid-gridsample1
plot(grid)
plot(gridsample1, col="lightblue", add=T)

gridsample2 <- grid
projection(gridsample2) <- projection(grid)

flat <- gridsample2
cm_25 <- gridsample1

plot(flat, col="white")
plot(cm_25, col='lightblue', add=TRUE)

# Add resistance values ---------------------------------------------------

flat@data$resistance <- runif(length(flat@data$ID), 0.95096189, 1)
cm_25@data$resistance <- runif(length(cm_25@data$ID), 0.68657211, 0.75614646)

crs(cm_25) <- crs(flat) 
grid_33_cm_25 <- rbind(flat, cm_25)

# Create a resistance raster ----------------------------------------------

flat_val <- flat@data$resistance
cm_25_val <- cm_25@data$resistance

mapvalues <- as.data.frame(cbind(flat_val, cm_25_val))
mapvalues <- data.frame(x = unlist(mapvalues))
row.names(mapvalues) <- NULL
mapvalues$x2 <- as.character(mapvalues$x)
mapvalues <- mapvalues[!duplicated(mapvalues),]

resistance_33_cm_25 <- resistanceFromShape(grid_33_cm_25,   
                                           res=0.25,      
                                           field="resistance",  
                                           mapvalues=with(mapvalues, setNames(x, x2)),    
                                           background=1, margin=2)       
plot(resistance_33_cm_25)

# Simulate limpet movement ------------------------------------------------

species.list <- lapply(1:1000, function(i) species(state.CRW(0.98)+0.44))

set.seed(123) 
tracks_33_cm25 <- simulate(species.list,
                           time = 288, 
                           coords = xyFromCell(resistance_33_cm_25, sample(1:ncell(resistance_33_cm_25), size=1000, replace=TRUE)),
                           resist=resistance_33_cm_25,
                           start.resistance=1)

# Save path length --------------------------------------------------------

traj_list <- vector("list", 1000)

for (i in seq_len(1000)) {
  cols <- ((i-1)*3+1):(i*3)
  limpet_data <- tracks_33_cm25[, cols]
  traj_list[[i]] <- sampleMovement(limpet_data, resolution=1, resist=resistance_33_cm_25)
}

all_steplengths <- list()
all_resistance <- list()
limpet_id <- list()

for (i in seq_along(traj_list)) {
  n <- length(traj_list[[i]][["stats"]][["steplengths"]])
  all_steplengths[[i]] <- traj_list[[i]][["stats"]][["steplengths"]]
  all_resistance[[i]] <- traj_list[[i]][["stats"]][["resistance"]]
  limpet_id[[i]] <- rep(i, n)
}

steplengths_vec <- unlist(all_steplengths)
resistance_vec <- unlist(all_resistance)
limpet_vec <- factor(unlist(limpet_id))

df <- data.frame(
  steplength = steplengths_vec,
  resistance = resistance_vec,
  limpet = limpet_vec
)

steplength <- df %>%
  group_by(limpet) %>%
  summarise(path_length=sum(steplength, na.rm=TRUE))

# Track plot --------------------------------------------------------------

tracks_33_cm25 <- tracks_33_cm25[, -seq(3, ncol(tracks_33_cm25), by=3)]

tracks_33_cm25 <- lapply(1:1000, function(i) {
  df <- data.frame(
    x = tracks_33_cm25[, 2*i - 1],
    y = tracks_33_cm25[, 2*i]
  )
  return(df)
})

plot(flat, col="azure2")
plot(cm_25, col="lightyellow", add=TRUE)
cols <- magma(1000)
for (i in 1:1000){
  lines(tracks_33_cm25[[i]], col=cols[i], lty=2, lwd=1)
}

# Plot a heat map ---------------------------------------------------------

tracks_33_cm25_df <- as.data.frame(do.call(rbind, tracks_33_cm25))

ggplot()+
  stat_density2d(tracks_33_cm25_df, mapping=aes(x=x, y=y, fill=..density..), geom='raster', contour=FALSE)+
  scale_fill_viridis(option="magma", direction=-1, limits=c(0, 0.5), na.value="white", expand=FALSE, oob=squish) +
  geom_polygon(grid, mapping=aes(x=long, y=lat, group=group), colour="white", fill=NA)

# Save the R environment --------------------------------------------------

# save.image(file='./R environments/cm25_Grouped_33percent.RData')
# load('./R environments/cm25_Grouped_33percent.RData')

# -------------------------------------------------------------------------
# ------------------------ ~33% coverage, 5 cm high -----------------------
# -------------------------------------------------------------------------
# Create a 8 x 16 grid ----------------------------------------------------

grid <- raster(xmn=0, xmx=4, ymn=0, ymx=2, ncol=16, nrow=8)
grid <- rasterToPolygons(grid)
plot(grid)

ID <- paste0(1:length(grid))
for (i in 1:length(slot(grid, "polygons"))) {
  slot(slot(grid, "polygons")[[i]], "ID") <- ID[i]
}
slot(grid, "data") <- cbind("ID" = 1:length(grid), slot(grid, "data"))

# Sample the grid to create tile sets -------------------------------------

gridsample1 <- grid[sample(length(ID), 1),] 
number <- 1
while (number < 2) { 
  grid <- grid       
  sample <- grid[sample(length(ID), 1),]
  for (i in 1:length(gridsample1)){    
    if (sample@data[["ID"]] == gridsample1@data[["ID"]][i]){  
      print("Same ID")  
      number <- number  
    } 
    else if (sample@data[["ID"]] == gridsample1@data[["ID"]][i]+1 |    
             sample@data[["ID"]] == gridsample1@data[["ID"]][i]-1 |
             sample@data[["ID"]] == gridsample1@data[["ID"]][i]+16 |
             sample@data[["ID"]] == gridsample1@data[["ID"]][i]-16) {
      gridsample1 <- rbind(gridsample1, sample) 
      gridsample1 <- gridsample1[!duplicated(gridsample1@data),]  
    } 
    if (length(gridsample1) == 42) {  
      number <- 2 
    }
  } 
}

projection(gridsample1) <- projection(grid)
grid <- grid-gridsample1
plot(grid)
plot(gridsample1, col="orange", add=T)

gridsample2 <- grid
projection(gridsample2) <- projection(grid)

flat <- gridsample2
cm_5 <- gridsample1

plot(flat, col="white")
plot(cm_5, col='orange', add=TRUE)

# Add resistance values ---------------------------------------------------

flat@data$resistance <- runif(length(flat@data$ID), 0.95096189, 1)
cm_5@data$resistance <- runif(length(cm_5@data$ID), 0.95783566, 1) 

crs(cm_5) <- crs(flat) 
grid_33_cm_5 <- rbind(flat, cm_5)

# Create a resistance raster ----------------------------------------------

flat_val <- flat@data$resistance
cm_5_val <- cm_5@data$resistance

mapvalues <- as.data.frame(cbind(flat_val, cm_5_val))
mapvalues <- data.frame(x = unlist(mapvalues))
row.names(mapvalues) <- NULL
mapvalues$x2 <- as.character(mapvalues$x)
mapvalues <- mapvalues[!duplicated(mapvalues),]

resistance_33_cm5 <- resistanceFromShape(grid_33_cm_5,   
                                          res=0.25,      
                                          field="resistance",  
                                          mapvalues=with(mapvalues, setNames(x, x2)),     
                                          background=1, margin=2)     
plot(resistance_33_cm5)

# Simulate limpet movement ------------------------------------------------

species.list <- lapply(1:1000, function(i) species(state.CRW(0.98)+0.44))

set.seed(123) 
tracks_33_cm5 <- simulate(species.list,
                           time = 288,
                           coords = xyFromCell(resistance_33_cm5, sample(1:ncell(resistance_33_cm5), size=1000, replace=TRUE)),
                           resist=resistance_33_cm5,
                           start.resistance=1)

# Save path length --------------------------------------------------------

traj_list <- vector("list", 1000)

for (i in seq_len(1000)) {
  cols <- ((i-1)*3+1):(i*3)
  limpet_data <- tracks_33_cm5[, cols]
  traj_list[[i]] <- sampleMovement(limpet_data, resolution=1, resist=resistance_33_cm5)
}

all_steplengths <- list()
all_resistance <- list()
limpet_id <- list()

for (i in seq_along(traj_list)) {
  n <- length(traj_list[[i]][["stats"]][["steplengths"]])
  all_steplengths[[i]] <- traj_list[[i]][["stats"]][["steplengths"]]
  all_resistance[[i]] <- traj_list[[i]][["stats"]][["resistance"]]
  limpet_id[[i]] <- rep(i, n)
}

steplengths_vec <- unlist(all_steplengths)
resistance_vec <- unlist(all_resistance)
limpet_vec <- factor(unlist(limpet_id))

df <- data.frame(
  steplength = steplengths_vec,
  resistance = resistance_vec,
  limpet = limpet_vec
)

steplength <- df %>%
  group_by(limpet) %>%
  summarise(path_length=sum(steplength, na.rm=TRUE))

# Track plot --------------------------------------------------------------

tracks_33_cm5 <- tracks_33_cm5[, -seq(3, ncol(tracks_33_cm5), by=3)]

tracks_33_cm5 <- lapply(1:1000, function(i) {
  df <- data.frame(
    x = tracks_33_cm5[, 2*i - 1],
    y = tracks_33_cm5[, 2*i]
  )
  return(df)
})

plot(flat, col="azure2")
plot(cm_5, col="lightyellow", add=TRUE)
cols <- magma(1000)
for (i in 1:1000){
  lines(tracks_33_cm5[[i]], col=cols[i], lty=2, lwd=1)
}

# Plot a heat map ---------------------------------------------------------

tracks_33_cm5_df <- as.data.frame(do.call(rbind, tracks_33_cm5))

ggplot()+
  stat_density2d(tracks_33_cm5_df, mapping=aes(x=x, y=y, fill=..density..), geom='raster', contour=FALSE)+
  scale_fill_viridis(option="magma", direction=-1, limits=c(0, 0.5), na.value="white", expand=FALSE, oob=squish) +
  geom_polygon(grid, mapping=aes(x=long, y=lat, group=group), colour="white", fill=NA)

# Save the R environment --------------------------------------------------

# save.image(file='./R environments/cm5_Grouped_33percent.RData')
# load('./R environments/cm5_Grouped_33percent.RData')

# -------------------------------------------------------------------------
# ----------------------- ~50% coverage, Ripple panel ---------------------
# -------------------------------------------------------------------------
# Create a 8 x 16 grid ----------------------------------------------------

grid <- raster(xmn=0, xmx=4, ymn=0, ymx=2, ncol=16, nrow=8)
grid <- rasterToPolygons(grid)
plot(grid)

ID <- paste0(1:length(grid))
for (i in 1:length(slot(grid, "polygons"))) {
  slot(slot(grid, "polygons")[[i]], "ID") <- ID[i]
}
slot(grid, "data") <- cbind("ID" = 1:length(grid), slot(grid, "data"))

# Sample the grid to create tile sets -------------------------------------

# The final area of coverage is ~50%, which is equal to 63 tiles (50% * 128). 

# Select a random sample of 63 polygons that are next to each other 
# NB: Sometimes you have to re-run the loop to get grid cells that are next to each other - I am sure there is some way to optimise this...
gridsample1 <- grid[sample(length(ID), 1),] 
number <- 1
while (number < 2) { 
  grid <- grid       
  sample <- grid[sample(length(ID), 1),] 
  for (i in 1:length(gridsample1)){    
    if (sample@data[["ID"]] == gridsample1@data[["ID"]][i]){  
      print("Same ID")  
      number <- number 
    } 
    else if (sample@data[["ID"]] == gridsample1@data[["ID"]][i]+1 |    
             sample@data[["ID"]] == gridsample1@data[["ID"]][i]-1 |
             sample@data[["ID"]] == gridsample1@data[["ID"]][i]+16 |
             sample@data[["ID"]] == gridsample1@data[["ID"]][i]-16) {
      gridsample1 <- rbind(gridsample1, sample) 
      gridsample1 <- gridsample1[!duplicated(gridsample1@data),]  
    } 
    if (length(gridsample1) == 63) {  
      number <- 2 
    }
  } 
}

projection(gridsample1) <- projection(grid)
grid <- grid-gridsample1
plot(grid)
plot(gridsample1, col="lightgreen", add=T)

gridsample2 <- grid
projection(gridsample2) <- projection(grid)
 
flat <- gridsample2
ripple <- gridsample1

plot(flat, col="white")
plot(ripple, col='lightgreen', add=TRUE)

# Add resistance values ---------------------------------------------------

flat@data$resistance <- runif(length(flat@data$ID), 0.95096189, 1)
ripple@data$resistance <- runif(length(ripple@data$ID), 0, 0.05446582)

crs(ripple) <- crs(flat) 
grid_50_ripp <- rbind(flat, ripple)

# Create a resistance raster ----------------------------------------------

flat_val <- flat@data$resistance
ripple_val <- ripple@data$resistance

mapvalues <- as.data.frame(cbind(flat_val, ripple_val))
mapvalues <- data.frame(x = unlist(mapvalues))
row.names(mapvalues) <- NULL
mapvalues$x2 <- as.character(mapvalues$x)
mapvalues <- mapvalues[!duplicated(mapvalues),]

resistance_50_ripp <- resistanceFromShape(grid_50_ripp,   
                                          res=0.25,      
                                          field="resistance",   
                                          mapvalues=with(mapvalues, setNames(x, x2)),    
                                          background=1, margin=2)    
plot(resistance_50_ripp)

# Simulate limpet movement ------------------------------------------------

species.list <- lapply(1:1000, function(i) species(state.CRW(0.98)+0.44))

set.seed(123) 
tracks_50_ripp <- simulate(species.list,
                          time = 288, 
                          coords = xyFromCell(resistance_50_ripp, sample(1:ncell(resistance_50_ripp), size=1000, replace=TRUE)),
                          resist=resistance_50_ripp,
                          start.resistance=1)

# Save path length --------------------------------------------------------

traj_list <- vector("list", 1000)

for (i in seq_len(1000)) {
  cols <- ((i-1)*3+1):(i*3)
  limpet_data <- tracks_50_ripp[, cols]
  traj_list[[i]] <- sampleMovement(limpet_data, resolution=1, resist=resistance_50_ripp)
}

all_steplengths <- list()
all_resistance <- list()
limpet_id <- list()

for (i in seq_along(traj_list)) {
  n <- length(traj_list[[i]][["stats"]][["steplengths"]])
  all_steplengths[[i]] <- traj_list[[i]][["stats"]][["steplengths"]]
  all_resistance[[i]] <- traj_list[[i]][["stats"]][["resistance"]]
  limpet_id[[i]] <- rep(i, n)
}

steplengths_vec <- unlist(all_steplengths)
resistance_vec <- unlist(all_resistance)
limpet_vec <- factor(unlist(limpet_id))

df <- data.frame(
  steplength = steplengths_vec,
  resistance = resistance_vec,
  limpet = limpet_vec
)

steplength <- df %>%
  group_by(limpet) %>%
  summarise(path_length=sum(steplength, na.rm=TRUE))

# Track plot --------------------------------------------------------------

tracks_50_ripp <- tracks_50_ripp[, -seq(3, ncol(tracks_50_ripp), by=3)]

tracks_50_ripp <- lapply(1:1000, function(i) {
  df <- data.frame(
    x = tracks_50_ripp[, 2*i - 1],
    y = tracks_50_ripp[, 2*i]
  )
  return(df)
})

plot(flat, col="azure2")
plot(ripple, col="lightyellow", add=TRUE)
cols <- magma(1000)
for (i in 1:1000){
  lines(tracks_50_ripp[[i]], col=cols[i], lty=2, lwd=1)
}

# Plot a heat map ---------------------------------------------------------

tracks_50_ripp_df <- as.data.frame(do.call(rbind, tracks_50_ripp))

ggplot()+
  stat_density2d(tracks_50_ripp_df, mapping=aes(x=x, y=y, fill=..density..), geom='raster', contour=FALSE)+
  scale_fill_viridis(option="magma", direction=-1, limits=c(0, 0.5), na.value="white", expand=FALSE, oob=squish) +
  geom_polygon(grid, mapping=aes(x=long, y=lat, group=group), colour="white", fill=NA)

# Save the R environment --------------------------------------------------

# save.image(file='./R environments/Ripple_Grouped_50percent.RData')
# load('./R environments/Ripple_Grouped_50percent.RData')

# -------------------------------------------------------------------------
# ----------------------- ~50% coverage, 2.5 cm high ----------------------
# -------------------------------------------------------------------------
# Create a 8 x 16 grid ----------------------------------------------------

grid <- raster(xmn=0, xmx=4, ymn=0, ymx=2, ncol=16, nrow=8)
grid <- rasterToPolygons(grid)
plot(grid)

ID <- paste0(1:length(grid))
for (i in 1:length(slot(grid, "polygons"))) {
  slot(slot(grid, "polygons")[[i]], "ID") <- ID[i]
}
slot(grid, "data") <- cbind("ID" = 1:length(grid), slot(grid, "data"))

# Sample the grid to create tile sets -------------------------------------

gridsample1 <- grid[sample(length(ID), 1),] 
number <- 1
while (number < 2) { 
  grid <- grid    
  sample <- grid[sample(length(ID), 1),] 
  for (i in 1:length(gridsample1)){   
    if (sample@data[["ID"]] == gridsample1@data[["ID"]][i]){  
      print("Same ID")  
      number <- number  
    } 
    else if (sample@data[["ID"]] == gridsample1@data[["ID"]][i]+1 |   
             sample@data[["ID"]] == gridsample1@data[["ID"]][i]-1 |
             sample@data[["ID"]] == gridsample1@data[["ID"]][i]+16 |
             sample@data[["ID"]] == gridsample1@data[["ID"]][i]-16) {
      gridsample1 <- rbind(gridsample1, sample) 
      gridsample1 <- gridsample1[!duplicated(gridsample1@data),]  
    } 
    if (length(gridsample1) == 63) {  
      number <- 2 
    }
  } 
}

projection(gridsample1) <- projection(grid)
grid <- grid-gridsample1
plot(grid)
plot(gridsample1, col="lightblue", add=T)

gridsample2 <- grid
projection(gridsample2) <- projection(grid)

flat <- gridsample2
cm_25 <- gridsample1

plot(flat, col="white")
plot(cm_25, col='lightblue', add=TRUE)

# Add resistance values ---------------------------------------------------

flat@data$resistance <- runif(length(flat@data$ID), 0.95096189, 1)
cm_25@data$resistance <- runif(length(cm_25@data$ID), 0.68657211, 0.75614646)

crs(cm_25) <- crs(flat)
grid_50_cm_25 <- rbind(flat, cm_25)

# Create a resistance raster ----------------------------------------------

flat_val <- flat@data$resistance
cm_25_val <- cm_25@data$resistance

mapvalues <- as.data.frame(cbind(flat_val, cm_25_val))
mapvalues <- data.frame(x = unlist(mapvalues))
row.names(mapvalues) <- NULL
mapvalues$x2 <- as.character(mapvalues$x)
mapvalues <- mapvalues[!duplicated(mapvalues),]

resistance_50_cm_25 <- resistanceFromShape(grid_50_cm_25,   
                                           res=0.25,    
                                           field="resistance",  
                                           mapvalues=with(mapvalues, setNames(x, x2)),     
                                           background=1, margin=2)       
plot(resistance_50_cm_25)

# Simulate limpet movement ------------------------------------------------

species.list <- lapply(1:1000, function(i) species(state.CRW(0.98)+0.44))

set.seed(123)  
tracks_50_cm25 <- simulate(species.list,
                           time = 288,
                           coords = xyFromCell(resistance_50_cm_25, sample(1:ncell(resistance_50_cm_25), size=1000, replace=TRUE)),
                           resist=resistance_50_cm_25,
                           start.resistance=1)

# Save path length --------------------------------------------------------

traj_list <- vector("list", 1000)

for (i in seq_len(1000)) {
  cols <- ((i-1)*3+1):(i*3)
  limpet_data <- tracks_50_cm25[, cols]
  traj_list[[i]] <- sampleMovement(limpet_data, resolution=1, resist=resistance_50_cm_25)
}

all_steplengths <- list()
all_resistance <- list()
limpet_id <- list()

for (i in seq_along(traj_list)) {
  n <- length(traj_list[[i]][["stats"]][["steplengths"]])
  all_steplengths[[i]] <- traj_list[[i]][["stats"]][["steplengths"]]
  all_resistance[[i]] <- traj_list[[i]][["stats"]][["resistance"]]
  limpet_id[[i]] <- rep(i, n)
}

steplengths_vec <- unlist(all_steplengths)
resistance_vec <- unlist(all_resistance)
limpet_vec <- factor(unlist(limpet_id))

df <- data.frame(
  steplength = steplengths_vec,
  resistance = resistance_vec,
  limpet = limpet_vec
)

steplength <- df %>%
  group_by(limpet) %>%
  summarise(path_length=sum(steplength, na.rm=TRUE))

# Track plot --------------------------------------------------------------

tracks_50_cm25 <- tracks_50_cm25[, -seq(3, ncol(tracks_50_cm25), by=3)]

tracks_50_cm25 <- lapply(1:1000, function(i) {
  df <- data.frame(
    x = tracks_50_cm25[, 2*i - 1],
    y = tracks_50_cm25[, 2*i]
  )
  return(df)
})

plot(flat, col="azure2")
plot(cm_25, col="lightyellow", add=TRUE)
cols <- magma(1000)
for (i in 1:1000){
  lines(tracks_50_cm25[[i]], col=cols[i], lty=2, lwd=1)
}

# Plot a heat map ---------------------------------------------------------

tracks_50_cm25_df <- as.data.frame(do.call(rbind, tracks_50_cm25))

ggplot()+
  stat_density2d(tracks_50_cm25_df, mapping=aes(x=x, y=y, fill=..density..), geom='raster', contour=FALSE)+
  scale_fill_viridis(option="magma", direction=-1, limits=c(0, 0.5), na.value="white", expand=FALSE, oob=squish) +
  geom_polygon(grid, mapping=aes(x=long, y=lat, group=group), colour="white", fill=NA)

# Save the R environment --------------------------------------------------

# save.image(file='./R environments/cm25_Grouped_50percent.RData')
# load('./R environments/cm25_Grouped_50percent.RData')

# -------------------------------------------------------------------------
# ------------------------ ~50% coverage, 5 cm high -----------------------
# -------------------------------------------------------------------------
# Create a 8 x 16 grid ----------------------------------------------------

grid <- raster(xmn=0, xmx=4, ymn=0, ymx=2, ncol=16, nrow=8)
grid <- rasterToPolygons(grid)
plot(grid)

ID <- paste0(1:length(grid))
for (i in 1:length(slot(grid, "polygons"))) {
  slot(slot(grid, "polygons")[[i]], "ID") <- ID[i]
}
slot(grid, "data") <- cbind("ID" = 1:length(grid), slot(grid, "data"))

# Sample the grid to create tile sets -------------------------------------

gridsample1 <- grid[sample(length(ID), 1),] 
number <- 1
while (number < 2) { 
  grid <- grid       
  sample <- grid[sample(length(ID), 1),] 
  for (i in 1:length(gridsample1)){    
    if (sample@data[["ID"]] == gridsample1@data[["ID"]][i]){  
      print("Same ID")  
      number <- number  
    } 
    else if (sample@data[["ID"]] == gridsample1@data[["ID"]][i]+1 |    
             sample@data[["ID"]] == gridsample1@data[["ID"]][i]-1 |
             sample@data[["ID"]] == gridsample1@data[["ID"]][i]+16 |
             sample@data[["ID"]] == gridsample1@data[["ID"]][i]-16) {
      gridsample1 <- rbind(gridsample1, sample) 
      gridsample1 <- gridsample1[!duplicated(gridsample1@data),]  
    } 
    if (length(gridsample1) == 63) {  
      number <- 2 
    }
  } 
}

projection(gridsample1) <- projection(grid)
grid <- grid-gridsample1
plot(grid)
plot(gridsample1, col="orange", add=T)

gridsample2 <- grid
projection(gridsample2) <- projection(grid)

flat <- gridsample2 
cm_5 <- gridsample1

plot(flat, col="white")
plot(cm_5, col='orange', add=TRUE)

# Add resistance values ---------------------------------------------------

flat@data$resistance <- runif(length(flat@data$ID), 0.95096189, 1)
cm_5@data$resistance <- runif(length(cm_5@data$ID), 0.95783566, 1) 

crs(cm_5) <- crs(flat) 
grid_50_cm_5 <- rbind(flat, cm_5)

# Create a resistance raster ----------------------------------------------

flat_val <- flat@data$resistance
cm_5_val <- cm_5@data$resistance

mapvalues <- as.data.frame(cbind(flat_val, cm_5_val))
mapvalues <- data.frame(x = unlist(mapvalues))
row.names(mapvalues) <- NULL
mapvalues$x2 <- as.character(mapvalues$x)
mapvalues <- mapvalues[!duplicated(mapvalues),]

resistance_50_cm_5 <- resistanceFromShape(grid_50_cm_5,  
                                          res=0.25,      
                                          field="resistance",   
                                          mapvalues=with(mapvalues, setNames(x, x2)),     
                                          background=1, margin=2)       
plot(resistance_50_cm_5)

# Simulate limpet movement ------------------------------------------------

species.list <- lapply(1:1000, function(i) species(state.CRW(0.98)+0.44))

set.seed(123)  
tracks_50_cm5 <- simulate(species.list,
                           time = 288, 
                           coords = xyFromCell(resistance_50_cm_5, sample(1:ncell(resistance_50_cm_5), size=1000, replace=TRUE)),
                           resist=resistance_50_cm_5,
                           start.resistance=1)

# Save path length --------------------------------------------------------

traj_list <- vector("list", 1000)

for (i in seq_len(1000)) {
  cols <- ((i-1)*3+1):(i*3)
  limpet_data <- tracks_50_cm5[, cols]
  traj_list[[i]] <- sampleMovement(limpet_data, resolution=1, resist=resistance_50_cm_5)
}

all_steplengths <- list()
all_resistance <- list()
limpet_id <- list()

for (i in seq_along(traj_list)) {
  n <- length(traj_list[[i]][["stats"]][["steplengths"]])
  all_steplengths[[i]] <- traj_list[[i]][["stats"]][["steplengths"]]
  all_resistance[[i]] <- traj_list[[i]][["stats"]][["resistance"]]
  limpet_id[[i]] <- rep(i, n)
}

steplengths_vec <- unlist(all_steplengths)
resistance_vec <- unlist(all_resistance)
limpet_vec <- factor(unlist(limpet_id))

df <- data.frame(
  steplength = steplengths_vec,
  resistance = resistance_vec,
  limpet = limpet_vec
)

steplength <- df %>%
  group_by(limpet) %>%
  summarise(path_length=sum(steplength, na.rm=TRUE))

# Track plot --------------------------------------------------------------

tracks_50_cm5 <- tracks_50_cm5[, -seq(3, ncol(tracks_50_cm5), by=3)]

tracks_50_cm5 <- lapply(1:1000, function(i) {
  df <- data.frame(
    x = tracks_50_cm5[, 2*i - 1],
    y = tracks_50_cm5[, 2*i]
  )
  return(df)
})

plot(flat, col="azure2")
plot(cm_5, col="lightyellow", add=TRUE)
cols <- magma(1000)
for (i in 1:1000){
  lines(tracks_50_cm5[[i]], col=cols[i], lty=2, lwd=1)
}

# Plot a heat map ---------------------------------------------------------

tracks_50_cm5_df <- as.data.frame(do.call(rbind, tracks_50_cm5))

ggplot()+
  stat_density2d(tracks_50_cm5_df, mapping=aes(x=x, y=y, fill=..density..), geom='raster', contour=FALSE)+
  scale_fill_viridis(option="magma", direction=-1, limits=c(0, 0.5), na.value="white", expand=FALSE, oob=squish) +
  geom_polygon(grid, mapping=aes(x=long, y=lat, group=group), colour="white", fill=NA)

# Save the R environment --------------------------------------------------

# save.image(file='./R environments/cm5_Grouped_50percent.RData')
# load('./R environments/cm5_Grouped_50percent.RData')

# -------------------------------------------------------------------------
# -------------------------- ~14% coverage, 3 panel -----------------------
# -------------------------------------------------------------------------
# Create a 8 x 16 grid ----------------------------------------------------

grid <- raster(xmn=0, xmx=4, ymn=0, ymx=2, ncol=16, nrow=8)
grid <- rasterToPolygons(grid)
plot(grid)

ID <- paste0(1:length(grid))
for (i in 1:length(slot(grid, "polygons"))) {
  slot(slot(grid, "polygons")[[i]], "ID") <- ID[i]
}
slot(grid, "data") <- cbind("ID" = 1:length(grid), slot(grid, "data"))

# Sample the grid to create tile sets -------------------------------------

gridsample1 <- grid[sample(length(ID), 1),] 
number <- 1
while (number < 2) { 
  grid <- grid       
  sample <- grid[sample(length(ID), 1),] 
  for (i in 1:length(gridsample1)){    
    if (sample@data[["ID"]] == gridsample1@data[["ID"]][i]){  
      print("Same ID")  
      number <- number 
    } 
    else if (sample@data[["ID"]] == gridsample1@data[["ID"]][i]+1 |    
             sample@data[["ID"]] == gridsample1@data[["ID"]][i]-1 |
             sample@data[["ID"]] == gridsample1@data[["ID"]][i]+16 |
             sample@data[["ID"]] == gridsample1@data[["ID"]][i]-16) {
      gridsample1 <- rbind(gridsample1, sample) 
      gridsample1 <- gridsample1[!duplicated(gridsample1@data),]  
    } 
    if (length(gridsample1) == 6) {  
      number <- 2 
    }
  } 
}
projection(gridsample1) <- projection(grid)
grid <- grid-gridsample1
plot(grid)
plot(gridsample1, col="lightgreen", add=T)

gridsample2 <- grid[sample(length(ID), 1),] 
number <- 1
while (number < 2) {
  grid <- grid       
  sample <- grid[sample(length(grid@data$ID), 1),] 
  for (i in 1:length(gridsample1)){
    if (sample@data[["ID"]] == gridsample1@data[["ID"]][1]){
      print("Same ID")
      number <- number
    }
  }
  for (i in 1:length(gridsample2)){    
    if (sample@data[["ID"]] == gridsample2@data[["ID"]][i]){  
      print("Same ID")  
      number <- number  
    } 
    else if (sample@data[["ID"]] == gridsample2@data[["ID"]][i]+1 |    
             sample@data[["ID"]] == gridsample2@data[["ID"]][i]-1 |
             sample@data[["ID"]] == gridsample2@data[["ID"]][i]+16 |
             sample@data[["ID"]] == gridsample2@data[["ID"]][i]-16) {
      gridsample2 <- rbind(gridsample2, sample) 
      gridsample2 <- gridsample2[!duplicated(gridsample2@data),]  
    } 
    if (length(gridsample2) == 6) {  
      number <- 2 
    }
  } 
}
projection(gridsample2) <- projection(grid)
grid <- grid-gridsample1
grid <- grid-gridsample2
plot(grid)
plot(gridsample1, col="lightgreen", add=T)
plot(gridsample2, col="lightblue", add=T)

gridsample3 <- grid[sample(length(ID), 1),] 
number <- 1
while (number < 2) { 
  grid <- grid       
  sample <- grid[sample(length(grid@data$ID), 1),] 
  for (i in 1:length(gridsample1)){
    if (sample@data[["ID"]] == gridsample1@data[["ID"]][1]){
      print("Same ID")
      number <- number
    }
  }
  for (i in 1:length(gridsample2)){
    
    if (sample@data[["ID"]] == gridsample2@data[["ID"]][1]){
      print("Same ID")
      number <- number
    }
  }
  for (i in 1:length(gridsample3)){    
    if (sample@data[["ID"]] == gridsample3@data[["ID"]][i]){  
      print("Same ID")  
      number <- number  
    } 
    else if (sample@data[["ID"]] == gridsample3@data[["ID"]][i]+1 |    
             sample@data[["ID"]] == gridsample3@data[["ID"]][i]-1 |
             sample@data[["ID"]] == gridsample3@data[["ID"]][i]+16 |
             sample@data[["ID"]] == gridsample3@data[["ID"]][i]-16) {
      gridsample3 <- rbind(gridsample3, sample) 
      gridsample3 <- gridsample3[!duplicated(gridsample3@data),] 
    } 
    if (length(gridsample3) == 6) {  
      number <- 2
    }
  } 
}
projection(gridsample3) <- projection(grid)
grid <- grid-gridsample3

gridsample4 <- grid
projection(gridsample4) <- projection(grid)

flat <- gridsample4 
ripple <- gridsample1
cm_25 <- gridsample2
cm_5 <- gridsample3

plot(flat, col="white")
plot(ripple, col='lightgreen', add=TRUE)
plot(cm_25, col='lightblue', add=TRUE)
plot(cm_5, col='orange', add=TRUE)

# Add resistance values ---------------------------------------------------

flat@data$resistance <- runif(length(flat@data$ID), 0.95096189, 1)
ripple@data$resistance <- runif(length(ripple@data$ID), 0, 0.05446582)
cm_25@data$resistance <- runif(length(cm_25@data$ID), 0.68657211, 0.75614646)
cm_5@data$resistance <- runif(length(cm_5@data$ID), 0.95783566, 1) 

crs(ripple) <- crs(flat) 
crs(cm_25) <- crs(ripple)
crs(cm_5) <- crs(cm_25)
grid_14_3tile <- rbind(flat, ripple, cm_25, cm_5)

# Create a resistance raster ----------------------------------------------

flat_val <- flat@data$resistance
ripple_val <- ripple@data$resistance
cm_25_val <- cm_25@data$resistance
cm_5_val <- cm_5@data$resistance

mapvalues <- as.data.frame(cbind(flat_val, ripple_val, cm_25_val, cm_5_val))
mapvalues <- data.frame(x = unlist(mapvalues))
row.names(mapvalues) <- NULL
mapvalues$x2 <- as.character(mapvalues$x)
mapvalues <- mapvalues[!duplicated(mapvalues),]

resistance_14_3tile <- resistanceFromShape(grid_14_3tile,  
                                           res=0.25,     
                                           field="resistance",   
                                           mapvalues=with(mapvalues, setNames(x, x2)),     
                                           background=1, margin=2)       
plot(resistance_14_3tile)

# Simulate limpet movement ------------------------------------------------

species.list <- lapply(1:1000, function(i) species(state.CRW(0.98)+0.44))

set.seed(123)  
tracks_14_3tile <- simulate(species.list,
                          time = 288,
                          coords = xyFromCell(resistance_14_3tile, sample(1:ncell(resistance_14_3tile), size=1000, replace=TRUE)),
                          resist=resistance_14_3tile,
                          start.resistance=1)

# Save path length --------------------------------------------------------

traj_list <- vector("list", 1000)

for (i in seq_len(1000)) {
  cols <- ((i-1)*3+1):(i*3)
  limpet_data <- tracks_14_3tile[, cols]
  traj_list[[i]] <- sampleMovement(limpet_data, resolution=1, resist=resistance_14_3tile)
}

all_steplengths <- list()
all_resistance <- list()
limpet_id <- list()

for (i in seq_along(traj_list)) {
  n <- length(traj_list[[i]][["stats"]][["steplengths"]])
  all_steplengths[[i]] <- traj_list[[i]][["stats"]][["steplengths"]]
  all_resistance[[i]] <- traj_list[[i]][["stats"]][["resistance"]]
  limpet_id[[i]] <- rep(i, n)
}

steplengths_vec <- unlist(all_steplengths)
resistance_vec <- unlist(all_resistance)
limpet_vec <- factor(unlist(limpet_id))

df <- data.frame(
  steplength = steplengths_vec,
  resistance = resistance_vec,
  limpet = limpet_vec
)

steplength <- df %>%
  group_by(limpet) %>%
  summarise(path_length=sum(steplength, na.rm=TRUE))

# Track plot --------------------------------------------------------------

tracks_14_3tile <- tracks_14_3tile[, -seq(3, ncol(tracks_14_3tile), by=3)]

tracks_14_3tile <- lapply(1:1000, function(i) {
  df <- data.frame(
    x = tracks_14_3tile[, 2*i - 1],
    y = tracks_14_3tile[, 2*i]
  )
  return(df)
})

plot(flat, col="azure2")
plot(ripple, col="lightyellow", add=TRUE)
plot(cm_25, col='lightblue', add=TRUE)
plot(cm_5, col='lightpink', add=TRUE)
cols <- magma(1000)
for (i in 1:1000){
  lines(tracks_14_3tile[[i]], col=cols[i], lty=2, lwd=1)
}

# Plot a heat map ---------------------------------------------------------

tracks_14_3tile_df <- as.data.frame(do.call(rbind, tracks_14_3tile))

ggplot()+
  stat_density2d(tracks_14_3tile_df, mapping=aes(x=x, y=y, fill=..density..), geom='raster', contour=FALSE)+
  scale_fill_viridis(option="magma", direction=-1, limits=c(0, 0.5), na.value="white", expand=FALSE, oob=squish) +
  geom_polygon(grid, mapping=aes(x=long, y=lat, group=group), colour="white", fill=NA)

# Save the R environment --------------------------------------------------

# save.image(file='./R environments/3tile_Grouped_14percent.RData')
# load('./R environments/3tile_Grouped_14percent.RData')

# -------------------------------------------------------------------------
# -------------------------- ~33% coverage, 3 panel -----------------------
# -------------------------------------------------------------------------
# Create a 8 x 16 grid ----------------------------------------------------

grid <- raster(xmn=0, xmx=4, ymn=0, ymx=2, ncol=16, nrow=8)
grid <- rasterToPolygons(grid)
plot(grid)

ID <- paste0(1:length(grid))
for (i in 1:length(slot(grid, "polygons"))) {
  slot(slot(grid, "polygons")[[i]], "ID") <- ID[i]
}
slot(grid, "data") <- cbind("ID" = 1:length(grid), slot(grid, "data"))

# Sample the grid to create tile sets -------------------------------------

gridsample1 <- grid[sample(length(ID), 1),] 
number <- 1
while (number < 2) { 
  grid <- grid       
  sample <- grid[sample(length(ID), 1),] 
  for (i in 1:length(gridsample1)){    
    if (sample@data[["ID"]] == gridsample1@data[["ID"]][i]){ 
      print("Same ID") 
      number <- number 
    } 
    else if (sample@data[["ID"]] == gridsample1@data[["ID"]][i]+1 |    
             sample@data[["ID"]] == gridsample1@data[["ID"]][i]-1 |
             sample@data[["ID"]] == gridsample1@data[["ID"]][i]+16 |
             sample@data[["ID"]] == gridsample1@data[["ID"]][i]-16) {
      gridsample1 <- rbind(gridsample1, sample) 
      gridsample1 <- gridsample1[!duplicated(gridsample1@data),]  
    } 
    if (length(gridsample1) == 14) {  
      number <- 2 
    }
  } 
}
projection(gridsample1) <- projection(grid)
grid <- grid-gridsample1
plot(grid)
plot(gridsample1, col="lightgreen", add=T)

gridsample2 <- grid[sample(length(ID), 1),] 
number <- 1
while (number < 2) { 
  grid <- grid      
  sample <- grid[sample(length(grid@data$ID), 1),] 
  for (i in 1:length(gridsample1)){
    if (sample@data[["ID"]] == gridsample1@data[["ID"]][1]){
      print("Same ID")
      number <- number
    }
  }
  for (i in 1:length(gridsample2)){    
    if (sample@data[["ID"]] == gridsample2@data[["ID"]][i]){ 
      print("Same ID")  
      number <- number 
    } 
    else if (sample@data[["ID"]] == gridsample2@data[["ID"]][i]+1 |    
             sample@data[["ID"]] == gridsample2@data[["ID"]][i]-1 |
             sample@data[["ID"]] == gridsample2@data[["ID"]][i]+16 |
             sample@data[["ID"]] == gridsample2@data[["ID"]][i]-16) {
      gridsample2 <- rbind(gridsample2, sample) 
      gridsample2 <- gridsample2[!duplicated(gridsample2@data),] 
    } 
    if (length(gridsample2) == 14) {  
      number <- 2
    }
  } 
}
projection(gridsample2) <- projection(grid)
grid <- grid-gridsample1
grid <- grid-gridsample2
plot(grid)
plot(gridsample1, col="lightgreen", add=T)
plot(gridsample2, col="lightblue", add=T)

gridsample3 <- grid[sample(length(ID), 1),] 
number <- 1
while (number < 2) { 
  grid <- grid       
  sample <- grid[sample(length(grid@data$ID), 1),] 
  for (i in 1:length(gridsample1)){
    if (sample@data[["ID"]] == gridsample1@data[["ID"]][1]){
      print("Same ID")
      number <- number
    }
  }
  for (i in 1:length(gridsample2)){
    if (sample@data[["ID"]] == gridsample2@data[["ID"]][1]){
      print("Same ID")
      number <- number
    }
  }
  for (i in 1:length(gridsample3)){   
    if (sample@data[["ID"]] == gridsample3@data[["ID"]][i]){  
      print("Same ID") 
      number <- number  
    } 
    else if (sample@data[["ID"]] == gridsample3@data[["ID"]][i]+1 |    
             sample@data[["ID"]] == gridsample3@data[["ID"]][i]-1 |
             sample@data[["ID"]] == gridsample3@data[["ID"]][i]+16 |
             sample@data[["ID"]] == gridsample3@data[["ID"]][i]-16) {
      gridsample3 <- rbind(gridsample3, sample) 
      gridsample3 <- gridsample3[!duplicated(gridsample3@data),]  
    }
    if (length(gridsample3) == 14) {  
      number <- 2 
    }
  } 
}
projection(gridsample3) <- projection(grid)
grid <- grid-gridsample3
gridsample4 <- grid
projection(gridsample4) <- projection(grid)

flat <- gridsample4 
ripple <- gridsample1
cm_25 <- gridsample2
cm_5 <- gridsample3

plot(flat, col="white")
plot(ripple, col='lightgreen', add=TRUE)
plot(cm_25, col='lightblue', add=TRUE)
plot(cm_5, col='orange', add=TRUE)

# Add resistance values ---------------------------------------------------

flat@data$resistance <- runif(length(flat@data$ID), 0.95096189, 1)
ripple@data$resistance <- runif(length(ripple@data$ID), 0, 0.05446582)
cm_25@data$resistance <- runif(length(cm_25@data$ID), 0.68657211, 0.75614646)
cm_5@data$resistance <- runif(length(cm_5@data$ID), 0.95783566, 1) 

crs(ripple) <- crs(flat) 
crs(cm_25) <- crs(ripple)
crs(cm_5) <- crs(cm_25)
grid_33_3tile <- rbind(flat, ripple, cm_25, cm_5)

# Create a resistance raster ----------------------------------------------

flat_val <- flat@data$resistance
ripple_val <- ripple@data$resistance
cm_25_val <- cm_25@data$resistance
cm_5_val <- cm_5@data$resistance

mapvalues <- as.data.frame(cbind(flat_val, ripple_val, cm_25_val, cm_5_val))
mapvalues <- data.frame(x = unlist(mapvalues))
row.names(mapvalues) <- NULL
mapvalues$x2 <- as.character(mapvalues$x)
mapvalues <- mapvalues[!duplicated(mapvalues),]

resistance_33_3tile <- resistanceFromShape(grid_33_3tile,    
                                           res=0.25,      
                                           field="resistance",   
                                           mapvalues=with(mapvalues, setNames(x, x2)),     
                                           background=1, margin=2)       
plot(resistance_33_3tile)

# Simulate limpet movement ------------------------------------------------

species.list <- lapply(1:1000, function(i) species(state.CRW(0.98)+0.44))

set.seed(123) 
tracks_33_3tile <- simulate(species.list,
                            time = 288, 
                            coords = xyFromCell(resistance_33_3tile, sample(1:ncell(resistance_33_3tile), size=1000, replace=TRUE)),
                            resist=resistance_33_3tile,
                            start.resistance=1)

# Save path length --------------------------------------------------------

traj_list <- vector("list", 1000)

for (i in seq_len(1000)) {
  cols <- ((i-1)*3+1):(i*3)
  limpet_data <- tracks_33_3tile[, cols]
  traj_list[[i]] <- sampleMovement(limpet_data, resolution=1, resist=resistance_33_3tile)
}

all_steplengths <- list()
all_resistance <- list()
limpet_id <- list()

for (i in seq_along(traj_list)) {
  n <- length(traj_list[[i]][["stats"]][["steplengths"]])
  all_steplengths[[i]] <- traj_list[[i]][["stats"]][["steplengths"]]
  all_resistance[[i]] <- traj_list[[i]][["stats"]][["resistance"]]
  limpet_id[[i]] <- rep(i, n)
}

steplengths_vec <- unlist(all_steplengths)
resistance_vec <- unlist(all_resistance)
limpet_vec <- factor(unlist(limpet_id))

df <- data.frame(
  steplength = steplengths_vec,
  resistance = resistance_vec,
  limpet = limpet_vec
)

steplength <- df %>%
  group_by(limpet) %>%
  summarise(path_length=sum(steplength, na.rm=TRUE))

# Track plot --------------------------------------------------------------

tracks_33_3tile <- tracks_33_3tile[, -seq(3, ncol(tracks_33_3tile), by=3)]

tracks_33_3tile <- lapply(1:1000, function(i) {
  df <- data.frame(
    x = tracks_33_3tile[, 2*i - 1],
    y = tracks_33_3tile[, 2*i]
  )
  return(df)
})

plot(flat, col="azure2")
plot(ripple, col="lightyellow", add=TRUE)
plot(cm_25, col='lightblue', add=TRUE)
plot(cm_5, col='lightpink', add=TRUE)
cols <- magma(1000)
for (i in 1:1000){
  lines(tracks_33_3tile[[i]], col=cols[i], lty=2, lwd=1)
}

# Plot a heat map ---------------------------------------------------------

tracks_33_3tile_df <- as.data.frame(do.call(rbind, tracks_33_3tile))

ggplot()+
  stat_density2d(tracks_33_3tile_df, mapping=aes(x=x, y=y, fill=..density..), geom='raster', contour=FALSE)+
  scale_fill_viridis(option="magma", direction=-1, limits=c(0, 0.5), na.value="white", expand=FALSE, oob=squish) +
  geom_polygon(grid, mapping=aes(x=long, y=lat, group=group), colour="white", fill=NA)

# Save the R environment --------------------------------------------------

save.image(file='./R environments/3tile_Grouped_33percent.RData')
load('./R environments/3tile_Grouped_33percent.RData')

# -------------------------------------------------------------------------
# -------------------------- ~50% coverage, 3 panel -----------------------
# -------------------------------------------------------------------------
# Create a 8 x 16 grid -----------------------------------------------------------

grid <- raster(xmn = 0, xmx = 4, ymn = 0, ymx = 2, ncol = 16, nrow = 8)
grid <- rasterToPolygons(grid)
plot(grid)

ID <- paste0(1:length(grid))
for (i in 1:length(slot(grid, "polygons"))) {
  slot(slot(grid, "polygons")[[i]], "ID") <- ID[i]
}
slot(grid, "data") <- cbind("ID" = 1:length(grid), slot(grid, "data"))

# Sample the grid to create tile sets -------------------------------------

gridsample1 <- grid[sample(length(ID), 1),] 
number <- 1
while (number < 2) { 
  grid <- grid    
  sample <- grid[sample(length(ID), 1),] 
  for (i in 1:length(gridsample1)){    
    if (sample@data[["ID"]] == gridsample1@data[["ID"]][i]){  
      print("Same ID")  
      number <- number  
    }
    else if (sample@data[["ID"]] == gridsample1@data[["ID"]][i]+1 |    
             sample@data[["ID"]] == gridsample1@data[["ID"]][i]-1 |
             sample@data[["ID"]] == gridsample1@data[["ID"]][i]+16 |
             sample@data[["ID"]] == gridsample1@data[["ID"]][i]-16) {
      gridsample1 <- rbind(gridsample1, sample) 
      gridsample1 <- gridsample1[!duplicated(gridsample1@data),]  
    } 
    if (length(gridsample1) == 21) {  
      number <- 2 
    }
  } 
}
projection(gridsample1) <- projection(grid)
grid <- grid-gridsample1
plot(grid)
plot(gridsample1, col="lightgreen", add=T)

gridsample2 <- grid[sample(length(ID), 1),] 
number <- 1
while (number < 2) {
  grid <- grid      
  sample <- grid[sample(length(grid@data$ID), 1),]
  for (i in 1:length(gridsample1)){
    if (sample@data[["ID"]] == gridsample1@data[["ID"]][1]){
      print("Same ID")
      number <- number
    }
  }
  for (i in 1:length(gridsample2)){   
    if (sample@data[["ID"]] == gridsample2@data[["ID"]][i]){  
      print("Same ID")  
      number <- number  
    } 
    else if (sample@data[["ID"]] == gridsample2@data[["ID"]][i]+1 |    
             sample@data[["ID"]] == gridsample2@data[["ID"]][i]-1 |
             sample@data[["ID"]] == gridsample2@data[["ID"]][i]+16 |
             sample@data[["ID"]] == gridsample2@data[["ID"]][i]-16) {
      gridsample2 <- rbind(gridsample2, sample) 
      gridsample2 <- gridsample2[!duplicated(gridsample2@data),]  
    } 
    if (length(gridsample2) == 21) {  
      number <- 2 
    }
  } 
}
projection(gridsample2) <- projection(grid)
grid <- grid-gridsample1
grid <- grid-gridsample2
plot(grid)
plot(gridsample1, col="lightgreen", add=T)
plot(gridsample2, col="lightblue", add=T)

gridsample3 <- grid[sample(length(ID), 1),] 
number <- 1
while (number < 2) { 
  grid <- grid     
  sample <- grid[sample(length(grid@data$ID), 1),] 
  for (i in 1:length(gridsample1)){
    if (sample@data[["ID"]] == gridsample1@data[["ID"]][1]){
      print("Same ID")
      number <- number
    }
  }
  for (i in 1:length(gridsample2)){
    if (sample@data[["ID"]] == gridsample2@data[["ID"]][1]){
      print("Same ID")
      number <- number
    }
  }
  for (i in 1:length(gridsample3)){   
    if (sample@data[["ID"]] == gridsample3@data[["ID"]][i]){  
      print("Same ID")  
      number <- number  
    } 
    else if (sample@data[["ID"]] == gridsample3@data[["ID"]][i]+1 |    
             sample@data[["ID"]] == gridsample3@data[["ID"]][i]-1 |
             sample@data[["ID"]] == gridsample3@data[["ID"]][i]+16 |
             sample@data[["ID"]] == gridsample3@data[["ID"]][i]-16) {
      gridsample3 <- rbind(gridsample3, sample) 
      gridsample3 <- gridsample3[!duplicated(gridsample3@data),]  
    } 
    if (length(gridsample3) == 21) {  
      number <- 2 
    }
  } 
}
projection(gridsample3) <- projection(grid)
grid <- grid-gridsample3
gridsample4 <- grid
projection(gridsample4) <- projection(grid)
 
flat <- gridsample4 
ripple <- gridsample1
cm_25 <- gridsample2
cm_5 <- gridsample3

plot(flat, col="white")
plot(ripple, col='lightgreen', add=TRUE)
plot(cm_25, col='lightblue', add=TRUE)
plot(cm_5, col='orange', add=TRUE)

# Add resistance values ---------------------------------------------------

flat@data$resistance <- runif(length(flat@data$ID), 0.95096189, 1)
ripple@data$resistance <- runif(length(ripple@data$ID), 0, 0.05446582)
cm_25@data$resistance <- runif(length(cm_25@data$ID), 0.68657211, 0.75614646)
cm_5@data$resistance <- runif(length(cm_5@data$ID), 0.95783566, 1) 

crs(ripple) <- crs(flat) 
crs(cm_25) <- crs(ripple)
crs(cm_5) <- crs(cm_25)
grid_50_3tile <- rbind(flat, ripple, cm_25, cm_5)

# Create a resistance raster -------------------------------------------------------

flat_val <- flat@data$resistance
ripple_val <- ripple@data$resistance
cm_25_val <- cm_25@data$resistance
cm_5_val <- cm_5@data$resistance

mapvalues <- as.data.frame(cbind(flat_val, ripple_val, cm_25_val, cm_5_val))
mapvalues <- data.frame(x = unlist(mapvalues))
row.names(mapvalues) <- NULL
mapvalues$x2 <- as.character(mapvalues$x)
mapvalues <- mapvalues[!duplicated(mapvalues),]

resistance_50_3tile <- resistanceFromShape(grid_50_3tile,    
                                           res=0.25,     
                                           field="resistance",   
                                           mapvalues=with(mapvalues, setNames(x, x2)),     
                                           background=1, margin=2)       
plot(resistance_50_3tile)

# Simulate limpet movement ---------------------------------------------------------

species.list <- lapply(1:1000, function(i) species(state.CRW(0.98)+0.44))

set.seed(123) 
tracks_50_3tile <- simulate(species.list,
                            time = 288, 
                            coords = xyFromCell(resistance_50_3tile, sample(1:ncell(resistance_50_3tile), size=1000, replace=TRUE)),
                            resist=resistance_50_3tile,
                            start.resistance=1)

# Save path length --------------------------------------------------------

traj_list <- vector("list", 1000)

for (i in seq_len(1000)) {
  cols <- ((i-1)*3+1):(i*3)
  limpet_data <- tracks_50_3tile[, cols]
  traj_list[[i]] <- sampleMovement(limpet_data, resolution=1, resist=resistance_50_3tile)
}

all_steplengths <- list()
all_resistance <- list()
limpet_id <- list()

for (i in seq_along(traj_list)) {
  n <- length(traj_list[[i]][["stats"]][["steplengths"]])
  all_steplengths[[i]] <- traj_list[[i]][["stats"]][["steplengths"]]
  all_resistance[[i]] <- traj_list[[i]][["stats"]][["resistance"]]
  limpet_id[[i]] <- rep(i, n)
}

steplengths_vec <- unlist(all_steplengths)
resistance_vec <- unlist(all_resistance)
limpet_vec <- factor(unlist(limpet_id))

df <- data.frame(
  steplength = steplengths_vec,
  resistance = resistance_vec,
  limpet = limpet_vec
)

steplength <- df %>%
  group_by(limpet) %>%
  summarise(path_length=sum(steplength, na.rm=TRUE))

# Track plot --------------------------------------------------------------

tracks_50_3tile <- tracks_50_3tile[, -seq(3, ncol(tracks_50_3tile), by=3)]

tracks_50_3tile <- lapply(1:1000, function(i) {
  df <- data.frame(
    x = tracks_50_3tile[, 2*i - 1],
    y = tracks_50_3tile[, 2*i]
  )
  return(df)
})

plot(flat, col="azure2")
plot(ripple, col="lightyellow", add=TRUE)
plot(cm_25, col='lightblue', add=TRUE)
plot(cm_5, col='lightpink', add=TRUE)
cols <- magma(1000)
for (i in 1:1000){
  lines(tracks_50_3tile[[i]], col=cols[i], lty=2, lwd=1)
}

# Plot a heat map ---------------------------------------------------------

tracks_50_3tile_df <- as.data.frame(do.call(rbind, tracks_50_3tile))

ggplot()+
  stat_density2d(tracks_50_3tile_df, mapping=aes(x=x, y=y, fill=..density..), geom='raster', contour=FALSE)+
  scale_fill_viridis(option="magma", direction=-1, limits=c(0, 0.5), na.value="white", expand=FALSE,oob=squish) +
  geom_polygon(grid, mapping=aes(x=long, y=lat, group=group), colour="white", fill=NA)

# Save the R environment --------------------------------------------------

# save.image(file='./R environments/3tile_Grouped_50percent.RData')
# load('./R environments/3tile_Grouped_50percent.RData')

# -------------------------------------------------------------------------
# --------------------------------- Random --------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# ----------------------- ~14% coverage, Ripple panel ---------------------
# -------------------------------------------------------------------------
# Create a 8 x 16 grid ----------------------------------------------------

grid <- raster(xmn=0, xmx=4, ymn=0, ymx=2, ncol=16, nrow=8)
grid <- rasterToPolygons(grid)
plot(grid)

ID <- paste0(1:length(grid))
for (i in 1:length(slot(grid, "polygons"))) {
  slot(slot(grid, "polygons")[[i]], "ID") <- ID[i]
}
slot(grid, "data") <- cbind("ID" = 1:length(grid), slot(grid, "data"))

# Sample the grid to create tile sets -------------------------------------

gridsample1 <- grid[sample(length(ID), 18),] 
projection(gridsample1) <- projection(grid)
grid <- grid-gridsample1
plot(grid)
plot(gridsample1, col="lightgreen", add=T)

gridsample2 <- grid
projection(gridsample2) <- projection(grid)

flat <- gridsample2 
ripple <- gridsample1

plot(flat, col="white")
plot(ripple, col='lightgreen', add=TRUE)

# Add resistance values ---------------------------------------------------

flat@data$resistance <- runif(length(flat@data$ID), 0.95096189, 1)
ripple@data$resistance <- runif(length(ripple@data$ID), 0, 0.05446582)

crs(ripple) <- crs(flat) 
grid_14_ripp <- rbind(flat, ripple)

# Create a resistance raster ----------------------------------------------

flat_val <- flat@data$resistance
ripple_val <- ripple@data$resistance

mapvalues <- as.data.frame(cbind(flat_val, ripple_val))
mapvalues <- data.frame(x = unlist(mapvalues))
row.names(mapvalues) <- NULL
mapvalues$x2 <- as.character(mapvalues$x)
mapvalues <- mapvalues[!duplicated(mapvalues),]

resistance_14_ripp <- resistanceFromShape(grid_14_ripp,    
                                          res=0.25,      
                                          field="resistance",   
                                          mapvalues=with(mapvalues, setNames(x, x2)),     
                                          background=1, margin=2)       
plot(resistance_14_ripp)

# Simulate limpet movement ------------------------------------------------

species.list <- lapply(1:1000, function(i) species(state.CRW(0.98)+0.44))

set.seed(123)  
tracks_14_ripp <- simulate(species.list,
                            time = 288, 
                            coords = xyFromCell(resistance_14_ripp, sample(1:ncell(resistance_14_ripp), size=1000, replace=TRUE)),
                            resist=resistance_14_ripp,
                            start.resistance=1)

# Save path length --------------------------------------------------------

traj_list <- vector("list", 1000)

for (i in seq_len(1000)) {
  cols <- ((i-1)*3+1):(i*3)
  limpet_data <- tracks_14_ripp[, cols]
  traj_list[[i]] <- sampleMovement(limpet_data, resolution=1, resist=resistance_14_ripp)
}

all_steplengths <- list()
all_resistance <- list()
limpet_id <- list()

for (i in seq_along(traj_list)) {
  n <- length(traj_list[[i]][["stats"]][["steplengths"]])
  all_steplengths[[i]] <- traj_list[[i]][["stats"]][["steplengths"]]
  all_resistance[[i]] <- traj_list[[i]][["stats"]][["resistance"]]
  limpet_id[[i]] <- rep(i, n)
}

steplengths_vec <- unlist(all_steplengths)
resistance_vec <- unlist(all_resistance)
limpet_vec <- factor(unlist(limpet_id))

df <- data.frame(
  steplength = steplengths_vec,
  resistance = resistance_vec,
  limpet = limpet_vec
)

steplength <- df %>%
  group_by(limpet) %>%
  summarise(path_length=sum(steplength, na.rm=TRUE))

# Track plot --------------------------------------------------------------

tracks_14_ripp <- tracks_14_ripp[, -seq(3, ncol(tracks_14_ripp), by=3)]

tracks_14_ripp <- lapply(1:1000, function(i) {
  df <- data.frame(
    x = tracks_14_ripp[, 2*i - 1],
    y = tracks_14_ripp[, 2*i]
  )
  return(df)
})

plot(flat, col="azure2")
plot(ripple, col="lightyellow", add=TRUE)
cols <- magma(1000)
for (i in 1:1000){
  lines(tracks_14_ripp[[i]], col=cols[i], lty=2, lwd=1)
}

# Plot a heat map ---------------------------------------------------------

tracks_14_ripp_df <- as.data.frame(do.call(rbind, tracks_14_ripp))

ggplot()+
  stat_density2d(tracks_14_ripp_df, mapping=aes(x=x, y=y, fill=..density..), geom='raster', contour=FALSE)+
  scale_fill_viridis(option="magma", direction=-1, limits=c(0, 0.5), na.value="white", expand=FALSE, oob=squish) +
  geom_polygon(grid, mapping=aes(x=long, y=lat, group=group), colour="white", fill=NA)

# Save the R environment --------------------------------------------------

# save.image(file='./R environments/Ripple_Random_14percent.RData')
# load('./R environments/Ripple_Random_14percent.RData')

# -------------------------------------------------------------------------
# ----------------------- ~14% coverage, 2.5 cm high ----------------------
# -------------------------------------------------------------------------
# Create a 8 x 16 grid ----------------------------------------------------

grid <- raster(xmn=0, xmx=4, ymn=0, ymx=2, ncol=16, nrow=8)
grid <- rasterToPolygons(grid)
plot(grid)

ID <- paste0(1:length(grid))
for (i in 1:length(slot(grid, "polygons"))) {
  slot(slot(grid, "polygons")[[i]], "ID") <- ID[i]
}
slot(grid, "data") <- cbind("ID" = 1:length(grid), slot(grid, "data"))

# Sample the grid to create tile sets -------------------------------------

gridsample1 <- grid[sample(length(ID), 18),] 
projection(gridsample1) <- projection(grid)
grid <- grid-gridsample1
plot(grid)
plot(gridsample1, col="lightblue", add=T)

gridsample2 <- grid
projection(gridsample2) <- projection(grid)

flat <- gridsample2
cm_25 <- gridsample1

plot(flat, col="white")
plot(cm_25, col='lightblue', add=TRUE)

# Add resistance values ---------------------------------------------------

flat@data$resistance <- runif(length(flat@data$ID), 0.95096189, 1)
cm_25@data$resistance <- runif(length(cm_25@data$ID), 0.68657211, 0.75614646)

crs(cm_25) <- crs(flat) 
grid_14_cm_25 <- rbind(flat, cm_25)

# Create a resistance raster ----------------------------------------------

flat_val <- flat@data$resistance
cm_25_val <- cm_25@data$resistance

mapvalues <- as.data.frame(cbind(flat_val, cm_25_val))
mapvalues <- data.frame(x = unlist(mapvalues))
row.names(mapvalues) <- NULL
mapvalues$x2 <- as.character(mapvalues$x)
mapvalues <- mapvalues[!duplicated(mapvalues),]

resistance_14_cm_25 <- resistanceFromShape(grid_14_cm_25,   
                                           res=0.25,      
                                           field="resistance",   
                                           mapvalues=with(mapvalues, setNames(x, x2)),    
                                           background=1, margin=2)       
plot(resistance_14_cm_25)

# Simulate limpet movement ------------------------------------------------

species.list <- lapply(1:1000, function(i) species(state.CRW(0.98)+0.44))

set.seed(123)
tracks_14_cm25 <- simulate(species.list,
                           time = 288, 
                           coords = xyFromCell(resistance_14_cm_25, sample(1:ncell(resistance_14_cm_25), size=1000, replace=TRUE)),
                           resist=resistance_14_cm_25,
                           start.resistance=1)

# Save path length --------------------------------------------------------

traj_list <- vector("list", 1000)

for (i in seq_len(1000)) {
  cols <- ((i-1)*3+1):(i*3)
  limpet_data <- tracks_14_cm25[, cols]
  traj_list[[i]] <- sampleMovement(limpet_data, resolution=1, resist=resistance_14_cm_25)
}

all_steplengths <- list()
all_resistance <- list()
limpet_id <- list()

for (i in seq_along(traj_list)) {
  n <- length(traj_list[[i]][["stats"]][["steplengths"]])
  all_steplengths[[i]] <- traj_list[[i]][["stats"]][["steplengths"]]
  all_resistance[[i]] <- traj_list[[i]][["stats"]][["resistance"]]
  limpet_id[[i]] <- rep(i, n)
}

steplengths_vec <- unlist(all_steplengths)
resistance_vec <- unlist(all_resistance)
limpet_vec <- factor(unlist(limpet_id))

df <- data.frame(
  steplength = steplengths_vec,
  resistance = resistance_vec,
  limpet = limpet_vec
)

steplength <- df %>%
  group_by(limpet) %>%
  summarise(path_length=sum(steplength, na.rm=TRUE))

# Track plot --------------------------------------------------------------

tracks_14_cm25 <- tracks_14_cm25[, -seq(3, ncol(tracks_14_cm25), by=3)]

tracks_14_cm25 <- lapply(1:1000, function(i) {
  df <- data.frame(
    x = tracks_14_cm25[, 2*i - 1],
    y = tracks_14_cm25[, 2*i]
  )
  return(df)
})

plot(flat, col="azure2")
plot(cm_25, col="lightyellow", add=TRUE)
cols <- magma(1000)
for (i in 1:1000){
  lines(tracks_14_cm25[[i]], col=cols[i], lty=2, lwd=1)
}

# Plot a heat map ---------------------------------------------------------

tracks_14_cm25_df <- as.data.frame(do.call(rbind, tracks_14_cm25))

ggplot()+
  stat_density2d(tracks_14_cm25_df, mapping=aes(x=x, y=y, fill=..density..), geom='raster', contour=FALSE)+
  scale_fill_viridis(option="magma", direction=-1, limits=c(0, 0.5), na.value="white", expand=FALSE, oob=squish) +
  geom_polygon(grid, mapping=aes(x=long, y=lat, group=group), colour="white", fill=NA)

# Save the R environment --------------------------------------------------

# save.image(file='./R environments/cm25_Random_14percent.RData')
# load('./R environments/cm25_Random_14percent.RData')

# -------------------------------------------------------------------------
# ------------------------ ~14% coverage, 5 cm high -----------------------
# -------------------------------------------------------------------------
# Create a 8 x 16 grid ----------------------------------------------------

grid <- raster(xmn=0, xmx=4, ymn=0, ymx=2, ncol=16, nrow=8)
grid <- rasterToPolygons(grid)
plot(grid)

ID <- paste0(1:length(grid))
for (i in 1:length(slot(grid, "polygons"))) {
  slot(slot(grid, "polygons")[[i]], "ID") <- ID[i]
}
slot(grid, "data") <- cbind("ID" = 1:length(grid), slot(grid, "data"))

# Sample the grid to create tile sets -------------------------------------

gridsample1 <- grid[sample(length(ID), 18),] 
projection(gridsample1) <- projection(grid)
grid <- grid-gridsample1
plot(grid)
plot(gridsample1, col="orange", add=T)

gridsample2 <- grid
projection(gridsample2) <- projection(grid)

flat <- gridsample2
cm_5 <- gridsample1

plot(flat, col="white")
plot(cm_5, col='orange', add=TRUE)

# Add resistance values ---------------------------------------------------

flat@data$resistance <- runif(length(flat@data$ID), 0.95096189, 1)
cm_5@data$resistance <- runif(length(cm_5@data$ID), 0.95783566, 1) 

crs(cm_5) <- crs(flat) 
grid_14_cm_5 <- rbind(flat, cm_5)

# Create a resistance raster ----------------------------------------------

flat_val <- flat@data$resistance
cm_5_val <- cm_5@data$resistance

mapvalues <- as.data.frame(cbind(flat_val, cm_5_val))
mapvalues <- data.frame(x = unlist(mapvalues))
row.names(mapvalues) <- NULL
mapvalues$x2 <- as.character(mapvalues$x)
mapvalues <- mapvalues[!duplicated(mapvalues),]

resistance_14_cm_5 <- resistanceFromShape(grid_14_cm_5,    
                                          res=0.25,      
                                          field="resistance",   
                                          mapvalues=with(mapvalues, setNames(x, x2)),     
                                          background=1, margin=2)       
plot(resistance_14_cm_5)

# Simulate limpet movement ------------------------------------------------

species.list <- lapply(1:1000, function(i) species(state.CRW(0.98)+0.44))

set.seed(123)  
tracks_14_cm5 <- simulate(species.list,
                           time = 288, 
                           coords = xyFromCell(resistance_14_cm_5, sample(1:ncell(resistance_14_cm_5), size=1000, replace=TRUE)),
                           resist=resistance_14_cm_5,
                           start.resistance=1)

# Save path length --------------------------------------------------------

traj_list <- vector("list", 1000)

for (i in seq_len(1000)) {
  cols <- ((i-1)*3+1):(i*3)
  limpet_data <- tracks_14_cm5[, cols]
  traj_list[[i]] <- sampleMovement(limpet_data, resolution=1, resist=resistance_14_cm_5)
}

all_steplengths <- list()
all_resistance <- list()
limpet_id <- list()

for (i in seq_along(traj_list)) {
  n <- length(traj_list[[i]][["stats"]][["steplengths"]])
  all_steplengths[[i]] <- traj_list[[i]][["stats"]][["steplengths"]]
  all_resistance[[i]] <- traj_list[[i]][["stats"]][["resistance"]]
  limpet_id[[i]] <- rep(i, n)
}

steplengths_vec <- unlist(all_steplengths)
resistance_vec <- unlist(all_resistance)
limpet_vec <- factor(unlist(limpet_id))

df <- data.frame(
  steplength = steplengths_vec,
  resistance = resistance_vec,
  limpet = limpet_vec
)

steplength <- df %>%
  group_by(limpet) %>%
  summarise(path_length=sum(steplength, na.rm=TRUE))

# Track plot --------------------------------------------------------------

tracks_14_cm5 <- tracks_14_cm5[, -seq(3, ncol(tracks_14_cm5), by=3)]

tracks_14_cm5 <- lapply(1:1000, function(i) {
  df <- data.frame(
    x = tracks_14_cm5[, 2*i - 1],
    y = tracks_14_cm5[, 2*i]
  )
  return(df)
})

plot(flat, col="azure2")
plot(cm_5, col="lightyellow", add=TRUE)
cols <- magma(1000)
for (i in 1:1000){
  lines(tracks_14_cm5[[i]], col=cols[i], lty=2, lwd=1)
}

# Plot a heat map ---------------------------------------------------------

tracks_14_cm5_df <- as.data.frame(do.call(rbind, tracks_14_cm5))

ggplot()+
  stat_density2d(tracks_14_cm5_df, mapping=aes(x=x, y=y, fill=..density..), geom='raster', contour=FALSE)+
  scale_fill_viridis(option="magma", direction=-1, limits=c(0, 0.5), na.value="white", expand=FALSE, oob=squish) +
  geom_polygon(grid, mapping=aes(x=long, y=lat, group=group), colour="white", fill=NA)

# Save the R environment --------------------------------------------------

# save.image(file='./R environments/cm5_Random_14percent.RData')
# load('./R environments/cm5_Random_14percent.RData')

# -------------------------------------------------------------------------
# ----------------------- ~33% coverage, Ripple panel ---------------------
# -------------------------------------------------------------------------
# Create a 8 x 16 grid ----------------------------------------------------

grid <- raster(xmn=0, xmx=4, ymn=0, ymx=2, ncol=16, nrow=8)
grid <- rasterToPolygons(grid)
plot(grid)

ID <- paste0(1:length(grid))
for (i in 1:length(slot(grid, "polygons"))) {
  slot(slot(grid, "polygons")[[i]], "ID") <- ID[i]
}
slot(grid, "data") <- cbind("ID" = 1:length(grid), slot(grid, "data"))

# Sample the grid to create tile sets -------------------------------------

gridsample1 <- grid[sample(length(ID), 42),] 
projection(gridsample1) <- projection(grid)
grid <- grid-gridsample1
plot(grid)
plot(gridsample1, col="lightgreen", add=T)

gridsample2 <- grid
projection(gridsample2) <- projection(grid)

flat <- gridsample2 
ripple <- gridsample1

plot(flat, col="white")
plot(ripple, col='lightgreen', add=TRUE)

# Add resistance values ---------------------------------------------------

flat@data$resistance <- runif(length(flat@data$ID), 0.95096189, 1)
ripple@data$resistance <- runif(length(ripple@data$ID), 0, 0.05446582)

crs(ripple) <- crs(flat) 
grid_33_ripp <- rbind(flat, ripple)

# Create a resistance raster ----------------------------------------------

flat_val <- flat@data$resistance
ripple_val <- ripple@data$resistance

mapvalues <- as.data.frame(cbind(flat_val, ripple_val))
mapvalues <- data.frame(x = unlist(mapvalues))
row.names(mapvalues) <- NULL
mapvalues$x2 <- as.character(mapvalues$x)
mapvalues <- mapvalues[!duplicated(mapvalues),]

resistance_33_ripp <- resistanceFromShape(grid_33_ripp,  
                                          res=0.25,     
                                          field="resistance",   
                                          mapvalues=with(mapvalues, setNames(x, x2)),    
                                          background=1, margin=2)       
plot(resistance_33_ripp)

# Simulate limpet movement ------------------------------------------------

species.list <- lapply(1:1000, function(i) species(state.CRW(0.98)+0.44))

set.seed(123)  
tracks_33_ripp <- simulate(species.list,
                          time = 288, 
                          coords = xyFromCell(resistance_33_ripp, sample(1:ncell(resistance_33_ripp), size=1000, replace=TRUE)),
                          resist=resistance_33_ripp,
                          start.resistance=1)

# Save path length --------------------------------------------------------

traj_list <- vector("list", 1000)

for (i in seq_len(1000)) {
  cols <- ((i-1)*3+1):(i*3)
  limpet_data <- tracks_33_ripp[, cols]
  traj_list[[i]] <- sampleMovement(limpet_data, resolution=1, resist=resistance_33_ripp)
}

all_steplengths <- list()
all_resistance <- list()
limpet_id <- list()

for (i in seq_along(traj_list)) {
  n <- length(traj_list[[i]][["stats"]][["steplengths"]])
  all_steplengths[[i]] <- traj_list[[i]][["stats"]][["steplengths"]]
  all_resistance[[i]] <- traj_list[[i]][["stats"]][["resistance"]]
  limpet_id[[i]] <- rep(i, n)
}

steplengths_vec <- unlist(all_steplengths)
resistance_vec <- unlist(all_resistance)
limpet_vec <- factor(unlist(limpet_id))

df <- data.frame(
  steplength = steplengths_vec,
  resistance = resistance_vec,
  limpet = limpet_vec
)

steplength <- df %>%
  group_by(limpet) %>%
  summarise(path_length=sum(steplength, na.rm=TRUE))

# Track plot --------------------------------------------------------------

tracks_33_ripp <- tracks_33_ripp[, -seq(3, ncol(tracks_33_ripp), by=3)]

tracks_33_ripp <- lapply(1:1000, function(i) {
  df <- data.frame(
    x = tracks_33_ripp[, 2*i - 1],
    y = tracks_33_ripp[, 2*i]
  )
  return(df)
})

plot(flat, col="azure2")
plot(ripple, col="lightyellow", add=TRUE)
cols <- magma(1000)
for (i in 1:1000){
  lines(tracks_33_ripp[[i]], col=cols[i], lty=2, lwd=1)
}

# Plot a heat map ---------------------------------------------------------

tracks_33_ripp_df <- as.data.frame(do.call(rbind, tracks_33_ripp))

ggplot()+
  stat_density2d(tracks_33_ripp_df, mapping=aes(x=x, y=y, fill=..density..), geom='raster', contour=FALSE)+
  scale_fill_viridis(option="magma", direction=-1, limits=c(0, 0.5), na.value="white", expand=FALSE, oob=squish) +
  geom_polygon(grid, mapping=aes(x=long, y=lat, group=group), colour="white", fill=NA)

# Save the R environment --------------------------------------------------

# save.image(file='./R environments/Ripple_Random_33percent.RData')
# load('./R environments/Ripple_Random_33percent.RData')

# -------------------------------------------------------------------------
# ----------------------- ~33% coverage, 2.5 cm high ----------------------
# -------------------------------------------------------------------------
# Create a 8 x 16 grid ----------------------------------------------------

grid <- raster(xmn=0, xmx=4, ymn=0, ymx=2, ncol=16, nrow=8)
grid <- rasterToPolygons(grid)
plot(grid)

ID <- paste0(1:length(grid))
for (i in 1:length(slot(grid, "polygons"))) {
  slot(slot(grid, "polygons")[[i]], "ID") <- ID[i]
}
slot(grid, "data") <- cbind("ID" = 1:length(grid), slot(grid, "data"))

# Sample the grid to create tile sets -------------------------------------

gridsample1 <- grid[sample(length(ID), 42),] 
projection(gridsample1) <- projection(grid)
grid <- grid-gridsample1
plot(grid)
plot(gridsample1, col="lightblue", add=T)

gridsample2 <- grid
projection(gridsample2) <- projection(grid)

flat <- gridsample2 
cm_25 <- gridsample1

plot(flat, col="white")
plot(cm_25, col='lightblue', add=TRUE)

# Add resistance values ---------------------------------------------------

flat@data$resistance <- runif(length(flat@data$ID), 0.95096189, 1)
cm_25@data$resistance <- runif(length(cm_25@data$ID), 0.68657211, 0.75614646)

crs(cm_25) <- crs(flat) 
grid_33_cm_25 <- rbind(flat, cm_25)

# Create a resistance raster ----------------------------------------------

flat_val <- flat@data$resistance
cm_25_val <- cm_25@data$resistance

mapvalues <- as.data.frame(cbind(flat_val, cm_25_val))
mapvalues <- data.frame(x = unlist(mapvalues))
row.names(mapvalues) <- NULL
mapvalues$x2 <- as.character(mapvalues$x)
mapvalues <- mapvalues[!duplicated(mapvalues),]

resistance_33_cm_25 <- resistanceFromShape(grid_33_cm_25,   
                                           res=0.25,     
                                           field="resistance",   
                                           mapvalues=with(mapvalues, setNames(x, x2)),     
                                           background=1, margin=2)       
plot(resistance_33_cm_25)

# Simulate limpet movement ------------------------------------------------

species.list <- lapply(1:1000, function(i) species(state.CRW(0.98)+0.44))

set.seed(123)  
tracks_33_cm25 <- simulate(species.list,
                           time = 288,
                           coords = xyFromCell(resistance_33_cm_25, sample(1:ncell(resistance_33_cm_25), size=1000, replace=TRUE)),
                           resist=resistance_33_cm_25,
                           start.resistance=1)

# Save path length --------------------------------------------------------

traj_list <- vector("list", 1000)

for (i in seq_len(1000)) {
  cols <- ((i-1)*3+1):(i*3)
  limpet_data <- tracks_33_cm25[, cols]
  traj_list[[i]] <- sampleMovement(limpet_data, resolution=1, resist=resistance_33_cm_25)
}

all_steplengths <- list()
all_resistance <- list()
limpet_id <- list()

for (i in seq_along(traj_list)) {
  n <- length(traj_list[[i]][["stats"]][["steplengths"]])
  all_steplengths[[i]] <- traj_list[[i]][["stats"]][["steplengths"]]
  all_resistance[[i]] <- traj_list[[i]][["stats"]][["resistance"]]
  limpet_id[[i]] <- rep(i, n)
}

steplengths_vec <- unlist(all_steplengths)
resistance_vec <- unlist(all_resistance)
limpet_vec <- factor(unlist(limpet_id))

df <- data.frame(
  steplength = steplengths_vec,
  resistance = resistance_vec,
  limpet = limpet_vec
)

steplength <- df %>%
  group_by(limpet) %>%
  summarise(path_length=sum(steplength, na.rm=TRUE))

# Track plot --------------------------------------------------------------

tracks_33_cm25 <- tracks_33_cm25[, -seq(3, ncol(tracks_33_cm25), by=3)]

tracks_33_cm25 <- lapply(1:1000, function(i) {
  df <- data.frame(
    x = tracks_33_cm25[, 2*i - 1],
    y = tracks_33_cm25[, 2*i]
  )
  return(df)
})

plot(flat, col="azure2")
plot(cm_25, col="lightyellow", add=TRUE)
cols <- magma(1000)
for (i in 1:1000){
  lines(tracks_33_cm25[[i]], col=cols[i], lty=2, lwd=1)
}

# Plot a heat map ---------------------------------------------------------

tracks_33_cm_25_df <- as.data.frame(do.call(rbind, tracks_33_cm25))

ggplot()+
  stat_density2d(tracks_33_cm_25_df, mapping=aes(x=x, y=y, fill=..density..), geom='raster', contour=FALSE)+
  scale_fill_viridis(option="magma", direction=-1, limits=c(0, 0.5), na.value="white", expand=FALSE, oob=squish) +
  geom_polygon(grid, mapping=aes(x=long, y=lat, group=group), colour="white", fill=NA)

# Save the R environment --------------------------------------------------

# save.image(file='./R environments/cm25_Random_33percent.RData')
# load('./R environments/cm25_Random_33percent.RData')

# -------------------------------------------------------------------------
# ------------------------ ~33% coverage, 5 cm high -----------------------
# -------------------------------------------------------------------------
# Create a 8 x 16 grid ----------------------------------------------------

grid <- raster(xmn=0, xmx=4, ymn=0, ymx=2, ncol=16, nrow=8)
grid <- rasterToPolygons(grid)
plot(grid)

ID <- paste0(1:length(grid))
for (i in 1:length(slot(grid, "polygons"))) {
  slot(slot(grid, "polygons")[[i]], "ID") <- ID[i]
}
slot(grid, "data") <- cbind("ID" = 1:length(grid), slot(grid, "data"))

# Sample the grid to create tile sets -------------------------------------

gridsample1 <- grid[sample(length(ID), 42),] 
projection(gridsample1) <- projection(grid)
grid <- grid-gridsample1
plot(grid)
plot(gridsample1, col="orange", add=T)

gridsample2 <- grid
projection(gridsample2) <- projection(grid)

flat <- gridsample2 
cm_5 <- gridsample1

plot(flat, col="white")
plot(cm_5, col='orange', add=TRUE)

# Add resistance values ---------------------------------------------------

flat@data$resistance <- runif(length(flat@data$ID), 0.95096189, 1)
cm_5@data$resistance <- runif(length(cm_5@data$ID), 0.95783566, 1) 

crs(cm_5) <- crs(flat) 
grid_33_cm_5 <- rbind(flat, cm_5)

# Create a resistance raster ----------------------------------------------

flat_val <- flat@data$resistance
cm_5_val <- cm_5@data$resistance

mapvalues <- as.data.frame(cbind(flat_val, cm_5_val))
mapvalues <- data.frame(x = unlist(mapvalues))
row.names(mapvalues) <- NULL
mapvalues$x2 <- as.character(mapvalues$x)
mapvalues <- mapvalues[!duplicated(mapvalues),]

resistance_33_cm_5 <- resistanceFromShape(grid_33_cm_5,    
                                          res=0.25,      
                                          field="resistance",   
                                          mapvalues=with(mapvalues, setNames(x, x2)),     
                                          background=1, margin=2)       
plot(resistance_33_cm_5)

# Simulate limpet movement ------------------------------------------------

species.list <- lapply(1:1000, function(i) species(state.CRW(0.98)+0.44))

set.seed(123)  
tracks_33_cm5 <- simulate(species.list,
                           time = 288,
                           coords = xyFromCell(resistance_33_cm_5, sample(1:ncell(resistance_33_cm_5), size=1000, replace=TRUE)),
                           resist=resistance_33_cm_5,
                           start.resistance=1)

# Save path length --------------------------------------------------------

traj_list <- vector("list", 1000)

for (i in seq_len(1000)) {
  cols <- ((i-1)*3+1):(i*3)
  limpet_data <- tracks_33_cm5[, cols]
  traj_list[[i]] <- sampleMovement(limpet_data, resolution=1, resist=resistance_33_cm_5)
}

all_steplengths <- list()
all_resistance <- list()
limpet_id <- list()

for (i in seq_along(traj_list)) {
  n <- length(traj_list[[i]][["stats"]][["steplengths"]])
  all_steplengths[[i]] <- traj_list[[i]][["stats"]][["steplengths"]]
  all_resistance[[i]] <- traj_list[[i]][["stats"]][["resistance"]]
  limpet_id[[i]] <- rep(i, n)
}

steplengths_vec <- unlist(all_steplengths)
resistance_vec <- unlist(all_resistance)
limpet_vec <- factor(unlist(limpet_id))

df <- data.frame(
  steplength = steplengths_vec,
  resistance = resistance_vec,
  limpet = limpet_vec
)

steplength <- df %>%
  group_by(limpet) %>%
  summarise(path_length=sum(steplength, na.rm=TRUE))

# Track plot --------------------------------------------------------------

tracks_33_cm5 <- tracks_33_cm5[, -seq(3, ncol(tracks_33_cm5), by=3)]

tracks_33_cm5 <- lapply(1:1000, function(i) {
  df <- data.frame(
    x = tracks_33_cm5[, 2*i - 1],
    y = tracks_33_cm5[, 2*i]
  )
  return(df)
})

plot(flat, col="azure2")
plot(cm_5, col="lightyellow", add=TRUE)
cols <- magma(1000)
for (i in 1:1000){
  lines(tracks_33_cm5[[i]], col=cols[i], lty=2, lwd=1)
}

# Plot a heat map ---------------------------------------------------------

tracks_33_cm5_df <- as.data.frame(do.call(rbind, tracks_33_cm5))

ggplot()+
  stat_density2d(tracks_33_cm5_df, mapping=aes(x=x, y=y, fill=..density..), geom='raster', contour=FALSE)+
  scale_fill_viridis(option="magma", direction=-1, limits=c(0, 0.5), na.value="white", expand=FALSE, oob=squish) +
  geom_polygon(grid, mapping=aes(x=long, y=lat, group=group), colour="white", fill=NA)

# Save the R environment --------------------------------------------------

# save.image(file='./R environments/cm5_Random_33percent.RData')
# load('./R environments/cm5_Random_33percent.RData')

# -------------------------------------------------------------------------
# ----------------------- ~50% coverage, Ripple panel ---------------------
# -------------------------------------------------------------------------
# Create a 8 x 16 grid ----------------------------------------------------

grid <- raster(xmn=0, xmx=4, ymn=0, ymx=2, ncol=16, nrow=8)
grid <- rasterToPolygons(grid)
plot(grid)

ID <- paste0(1:length(grid))
for (i in 1:length(slot(grid, "polygons"))) {
  slot(slot(grid, "polygons")[[i]], "ID") <- ID[i]
}
slot(grid, "data") <- cbind("ID" = 1:length(grid), slot(grid, "data"))

# Sample the grid to create tile sets -------------------------------------

gridsample1 <- grid[sample(length(ID), 63),] 
projection(gridsample1) <- projection(grid)
grid <- grid-gridsample1
plot(grid)
plot(gridsample1, col="lightgreen", add=T)

gridsample2 <- grid
projection(gridsample2) <- projection(grid)

flat <- gridsample2
ripple <- gridsample1

plot(flat, col="white")
plot(ripple, col='lightgreen', add=TRUE)

# Add resistance values ---------------------------------------------------

flat@data$resistance <- runif(length(flat@data$ID), 0.95096189, 1)
ripple@data$resistance <- runif(length(ripple@data$ID), 0, 0.05446582)

crs(ripple) <- crs(flat) 
grid_50_ripp <- rbind(flat, ripple)

# Create a resistance raster ----------------------------------------------

flat_val <- flat@data$resistance
ripple_val <- ripple@data$resistance

mapvalues <- as.data.frame(cbind(flat_val, ripple_val))
mapvalues <- data.frame(x = unlist(mapvalues))
row.names(mapvalues) <- NULL
mapvalues$x2 <- as.character(mapvalues$x)
mapvalues <- mapvalues[!duplicated(mapvalues),]

resistance_50_ripp <- resistanceFromShape(grid_50_ripp,    
                                          res=0.25,     
                                          field="resistance",   
                                          mapvalues=with(mapvalues, setNames(x, x2)),    
                                          background=1, margin=2)       
plot(resistance_50_ripp)

# Simulate limpet movement ------------------------------------------------

species.list <- lapply(1:1000, function(i) species(state.CRW(0.98)+0.44))

set.seed(123) 
tracks_50_ripp <- simulate(species.list,
                          time = 288, 
                          coords = xyFromCell(resistance_50_ripp, sample(1:ncell(resistance_50_ripp), size=1000, replace=TRUE)),
                          resist=resistance_50_ripp,
                          start.resistance=1)

# Save path length --------------------------------------------------------

traj_list <- vector("list", 1000)

for (i in seq_len(1000)) {
  cols <- ((i-1)*3+1):(i*3)
  limpet_data <- tracks_50_ripp[, cols]
  traj_list[[i]] <- sampleMovement(limpet_data, resolution=1, resist=resistance_50_ripp)
}

all_steplengths <- list()
all_resistance <- list()
limpet_id <- list()

for (i in seq_along(traj_list)) {
  n <- length(traj_list[[i]][["stats"]][["steplengths"]])
  all_steplengths[[i]] <- traj_list[[i]][["stats"]][["steplengths"]]
  all_resistance[[i]] <- traj_list[[i]][["stats"]][["resistance"]]
  limpet_id[[i]] <- rep(i, n)
}

steplengths_vec <- unlist(all_steplengths)
resistance_vec <- unlist(all_resistance)
limpet_vec <- factor(unlist(limpet_id))

df <- data.frame(
  steplength = steplengths_vec,
  resistance = resistance_vec,
  limpet = limpet_vec
)

steplength <- df %>%
  group_by(limpet) %>%
  summarise(path_length=sum(steplength, na.rm=TRUE))

# Track plot --------------------------------------------------------------

tracks_50_ripp <- tracks_50_ripp[, -seq(3, ncol(tracks_50_ripp), by=3)]

tracks_50_ripp <- lapply(1:1000, function(i) {
  df <- data.frame(
    x = tracks_50_ripp[, 2*i - 1],
    y = tracks_50_ripp[, 2*i]
  )
  return(df)
})

plot(flat, col="azure2")
plot(ripple, col="lightyellow", add=TRUE)
cols <- magma(1000)
for (i in 1:1000){
  lines(tracks_50_ripp[[i]], col=cols[i], lty=2, lwd=1)
}

# Plot a heat map ---------------------------------------------------------

tracks_50_ripp_df <- as.data.frame(do.call(rbind, tracks_50_ripp))

ggplot()+
  stat_density2d(tracks_50_ripp_df, mapping=aes(x=x, y=y, fill=..density..), geom='raster', contour=FALSE)+
  scale_fill_viridis(option="magma", direction=-1, limits=c(0, 0.5), na.value="white", expand=FALSE, oob=squish) +
  geom_polygon(grid, mapping=aes(x=long, y=lat, group=group), colour="white", fill=NA)

# Save the R environment --------------------------------------------------

# save.image(file='./R environments/Ripple_Random_50percent.RData')
# load('./R environments/Ripple_Random_50percent.RData')

# -------------------------------------------------------------------------
# ----------------------- ~50% coverage, 2.5 cm high ----------------------
# -------------------------------------------------------------------------
# Create a 8 x 16 grid ----------------------------------------------------

grid <- raster(xmn=0, xmx=4, ymn=0, ymx=2, ncol=16, nrow=8)
grid <- rasterToPolygons(grid)
plot(grid)

ID <- paste0(1:length(grid))
for (i in 1:length(slot(grid, "polygons"))) {
  slot(slot(grid, "polygons")[[i]], "ID") <- ID[i]
}
slot(grid, "data") <- cbind("ID" = 1:length(grid), slot(grid, "data"))

# Sample the grid to create tile sets -------------------------------------

gridsample1 <- grid[sample(length(ID), 63),] 
projection(gridsample1) <- projection(grid)
grid <- grid-gridsample1
plot(grid)
plot(gridsample1, col="lightblue", add=T)

gridsample2 <- grid
projection(gridsample2) <- projection(grid)

flat <- gridsample2 
cm_25 <- gridsample1

plot(flat, col="white")
plot(cm_25, col='lightblue', add=TRUE)

# Add resistance values ---------------------------------------------------

flat@data$resistance <- runif(length(flat@data$ID), 0.95096189, 1)
cm_25@data$resistance <- runif(length(cm_25@data$ID), 0.68657211, 0.75614646)

crs(cm_25) <- crs(flat)
grid_50_cm_25 <- rbind(flat, cm_25)

# Create a resistance raster ----------------------------------------------

flat_val <- flat@data$resistance
cm_25_val <- cm_25@data$resistance

mapvalues <- as.data.frame(cbind(flat_val, cm_25_val))
mapvalues <- data.frame(x = unlist(mapvalues))
row.names(mapvalues) <- NULL
mapvalues$x2 <- as.character(mapvalues$x)
mapvalues <- mapvalues[!duplicated(mapvalues),]

resistance_50_cm_25 <- resistanceFromShape(grid_50_cm_25,   
                                           res=0.25,     
                                           field="resistance",   
                                           mapvalues=with(mapvalues, setNames(x, x2)),    
                                           background=1, margin=2)       
plot(resistance_50_cm_25)

# Simulate limpet movement ------------------------------------------------

species.list <- lapply(1:1000, function(i) species(state.CRW(0.98)+0.44))

set.seed(123) 
tracks_50_cm25 <- simulate(species.list,
                           time = 288, 
                           coords = xyFromCell(resistance_50_cm_25, sample(1:ncell(resistance_50_cm_25), size=1000, replace=TRUE)),
                           resist=resistance_50_cm_25,
                           start.resistance=1)

# Save path length --------------------------------------------------------

traj_list <- vector("list", 1000)

for (i in seq_len(1000)) {
  cols <- ((i-1)*3+1):(i*3)
  limpet_data <- tracks_50_cm25[, cols]
  traj_list[[i]] <- sampleMovement(limpet_data, resolution=1, resist=resistance_50_cm_25)
}

all_steplengths <- list()
all_resistance <- list()
limpet_id <- list()

for (i in seq_along(traj_list)) {
  n <- length(traj_list[[i]][["stats"]][["steplengths"]])
  all_steplengths[[i]] <- traj_list[[i]][["stats"]][["steplengths"]]
  all_resistance[[i]] <- traj_list[[i]][["stats"]][["resistance"]]
  limpet_id[[i]] <- rep(i, n)
}

steplengths_vec <- unlist(all_steplengths)
resistance_vec <- unlist(all_resistance)
limpet_vec <- factor(unlist(limpet_id))

df <- data.frame(
  steplength = steplengths_vec,
  resistance = resistance_vec,
  limpet = limpet_vec
)

steplength <- df %>%
  group_by(limpet) %>%
  summarise(path_length=sum(steplength, na.rm=TRUE))

# Track plot --------------------------------------------------------------

tracks_50_cm25 <- tracks_50_cm25[, -seq(3, ncol(tracks_50_cm25), by=3)]

tracks_50_cm25 <- lapply(1:1000, function(i) {
  df <- data.frame(
    x = tracks_50_cm25[, 2*i - 1],
    y = tracks_50_cm25[, 2*i]
  )
  return(df)
})

plot(flat, col="azure2")
plot(cm_25, col="lightyellow", add=TRUE)
cols <- magma(1000)
for (i in 1:1000){
  lines(tracks_50_cm25[[i]], col=cols[i], lty=2, lwd=1)
}

# Plot a heat map ---------------------------------------------------------

tracks_50_cm_25_df <- as.data.frame(do.call(rbind, tracks_50_cm25))

ggplot()+
  stat_density2d(tracks_50_cm_25_df, mapping=aes(x=x, y=y, fill=..density..), geom='raster', contour=FALSE)+
  scale_fill_viridis(option="magma", direction=-1, limits=c(0, 0.5), na.value="white", expand=FALSE, oob=squish) +
  geom_polygon(grid, mapping=aes(x=long, y=lat, group=group), colour="white", fill=NA)

# Save the R environment --------------------------------------------------

# save.image(file='./R environments/cm25_Random_50percent.RData')
# load('./R environments/cm25_Random_50percent.RData')

# -------------------------------------------------------------------------
# ------------------------ ~50% coverage, 5 cm high -----------------------
# -------------------------------------------------------------------------
# Create a 8 x 16 grid ----------------------------------------------------

grid <- raster(xmn=0, xmx=4, ymn=0, ymx=2, ncol=16, nrow=8)
grid <- rasterToPolygons(grid)
plot(grid)

ID <- paste0(1:length(grid))
for (i in 1:length(slot(grid, "polygons"))) {
  slot(slot(grid, "polygons")[[i]], "ID") <- ID[i]
}
slot(grid, "data") <- cbind("ID" = 1:length(grid), slot(grid, "data"))

# Sample the grid to create tile sets -------------------------------------

gridsample1 <- grid[sample(length(ID), 63),] 
projection(gridsample1) <- projection(grid)
grid <- grid-gridsample1
plot(grid)
plot(gridsample1, col="orange", add=T)

gridsample2 <- grid
projection(gridsample2) <- projection(grid)

flat <- gridsample2 
cm_5 <- gridsample1

plot(flat, col="white")
plot(cm_5, col='orange', add=TRUE)

# Add resistance values ---------------------------------------------------

flat@data$resistance <- runif(length(flat@data$ID), 0.95096189, 1)
cm_5@data$resistance <- runif(length(cm_5@data$ID), 0.95783566, 1) 

crs(cm_5) <- crs(flat) 
grid_50_cm_5 <- rbind(flat, cm_5)

# Create a resistance raster ----------------------------------------------

flat_val <- flat@data$resistance
cm_5_val <- cm_5@data$resistance

mapvalues <- as.data.frame(cbind(flat_val, cm_5_val))
mapvalues <- data.frame(x = unlist(mapvalues))
row.names(mapvalues) <- NULL
mapvalues$x2 <- as.character(mapvalues$x)
mapvalues <- mapvalues[!duplicated(mapvalues),]

resistance_50_cm_5 <- resistanceFromShape(grid_50_cm_5,    
                                          res=0.25,      
                                          field="resistance",   
                                          mapvalues=with(mapvalues, setNames(x, x2)),     
                                          background=1, margin=2)      
plot(resistance_50_cm_5)

# Simulate limpet movement ------------------------------------------------

species.list <- lapply(1:1000, function(i) species(state.CRW(0.98)+0.44))

set.seed(123)  
tracks_50_cm5 <- simulate(species.list,
                           time = 288, 
                           coords = xyFromCell(resistance_50_cm_5, sample(1:ncell(resistance_50_cm_5), size=1000, replace=TRUE)),
                           resist=resistance_50_cm_5,
                           start.resistance=1)

# Save path length --------------------------------------------------------

traj_list <- vector("list", 1000)

for (i in seq_len(1000)) {
  cols <- ((i-1)*3+1):(i*3)
  limpet_data <- tracks_50_cm5[, cols]
  traj_list[[i]] <- sampleMovement(limpet_data, resolution=1, resist=resistance_50_cm_5)
}

all_steplengths <- list()
all_resistance <- list()
limpet_id <- list()

for (i in seq_along(traj_list)) {
  n <- length(traj_list[[i]][["stats"]][["steplengths"]])
  all_steplengths[[i]] <- traj_list[[i]][["stats"]][["steplengths"]]
  all_resistance[[i]] <- traj_list[[i]][["stats"]][["resistance"]]
  limpet_id[[i]] <- rep(i, n)
}

steplengths_vec <- unlist(all_steplengths)
resistance_vec <- unlist(all_resistance)
limpet_vec <- factor(unlist(limpet_id))

df <- data.frame(
  steplength = steplengths_vec,
  resistance = resistance_vec,
  limpet = limpet_vec
)

steplength <- df %>%
  group_by(limpet) %>%
  summarise(path_length=sum(steplength, na.rm=TRUE))

# Track plot --------------------------------------------------------------

tracks_50_cm5 <- tracks_50_cm5[, -seq(3, ncol(tracks_50_cm5), by=3)]

tracks_50_cm5 <- lapply(1:1000, function(i) {
  df <- data.frame(
    x = tracks_50_cm5[, 2*i - 1],
    y = tracks_50_cm5[, 2*i]
  )
  return(df)
})

plot(flat, col="azure2")
plot(cm_5, col="lightyellow", add=TRUE)
cols <- magma(1000)
for (i in 1:1000){
  lines(tracks_50_cm5[[i]], col=cols[i], lty=2, lwd=1)
}

# Plot a heat map ---------------------------------------------------------

tracks_50_cm_5_df <- as.data.frame(do.call(rbind, tracks_50_cm5))

ggplot()+
  stat_density2d(tracks_50_cm_5_df, mapping=aes(x=x, y=y, fill=..density..), geom='raster', contour=FALSE)+
  scale_fill_viridis(option="magma", direction=-1, limits=c(0, 0.5), na.value="white", expand=FALSE, oob=squish) +
  geom_polygon(grid, mapping=aes(x=long, y=lat, group=group), colour="white", fill=NA)

# Save the R environment --------------------------------------------------

# save.image(file='./R environments/cm5_Random_50percent.RData')
# load('./R environments/cm5_Random_50percent.RData')

# -------------------------------------------------------------------------
# -------------------------- ~14% coverage, 3 panel -----------------------
# -------------------------------------------------------------------------
# Create a 8 x 16 grid ----------------------------------------------------

grid <- raster(xmn=0, xmx=4, ymn=0, ymx=2, ncol=16, nrow=8)
grid <- rasterToPolygons(grid)
plot(grid)

ID <- paste0(1:length(grid))
for (i in 1:length(slot(grid, "polygons"))) {
  slot(slot(grid, "polygons")[[i]], "ID") <- ID[i]
}
slot(grid, "data") <- cbind("ID" = 1:length(grid), slot(grid, "data"))

# Sample the grid to create tile sets -------------------------------------

gridsample1 <- grid[sample(length(ID), 6),] 
projection(gridsample1) <- projection(grid)
grid <- grid-gridsample1
plot(grid)
plot(gridsample1, col="lightgreen", add=T)

gridsample2 <- grid[sample(length(grid), 6),] 
projection(gridsample2) <- projection(grid)
grid <- grid-gridsample2
plot(grid)
plot(gridsample1, col="lightgreen", add=T)
plot(gridsample2, col="lightblue", add=T)

gridsample3 <- grid[sample(length(grid), 6),] 
projection(gridsample3) <- projection(grid)
grid <- grid-gridsample3

gridsample4 <- grid
projection(gridsample4) <- projection(grid)

flat <- gridsample4 
ripple <- gridsample1
cm_25 <- gridsample2
cm_5 <- gridsample3

plot(flat, col="white")
plot(ripple, col='lightgreen', add=TRUE)
plot(cm_25, col='lightblue', add=TRUE)
plot(cm_5, col='orange', add=TRUE)

# Add resistance values ---------------------------------------------------

flat@data$resistance <- runif(length(flat@data$ID), 0.95096189, 1)
ripple@data$resistance <- runif(length(ripple@data$ID), 0, 0.05446582)
cm_25@data$resistance <- runif(length(cm_25@data$ID), 0.68657211, 0.75614646)
cm_5@data$resistance <- runif(length(cm_5@data$ID), 0.95783566, 1) 

crs(ripple) <- crs(flat) 
crs(cm_25) <- crs(ripple)
crs(cm_5) <- crs(cm_25)
grid_14_3tile <- rbind(flat, ripple, cm_25, cm_5)

# Create a resistance raster ----------------------------------------------

flat_val <- flat@data$resistance
ripple_val <- ripple@data$resistance
cm_25_val <- cm_25@data$resistance
cm_5_val <- cm_5@data$resistance

mapvalues <- as.data.frame(cbind(flat_val, ripple_val, cm_25_val, cm_5_val))
mapvalues <- data.frame(x = unlist(mapvalues))
row.names(mapvalues) <- NULL
mapvalues$x2 <- as.character(mapvalues$x)
mapvalues <- mapvalues[!duplicated(mapvalues),]

resistance_14_3tile <- resistanceFromShape(grid_14_3tile,   
                                           res=0.25,      
                                           field="resistance",   
                                           mapvalues=with(mapvalues, setNames(x, x2)),    
                                           background=1, margin=2)       
plot(resistance_14_3tile)

# Simulate limpet movement ------------------------------------------------

species.list <- lapply(1:1000, function(i) species(state.CRW(0.98)+0.44))

set.seed(123)  
tracks_14_3tile <- simulate(species.list,
                          time = 288, 
                          coords = xyFromCell(resistance_14_3tile, sample(1:ncell(resistance_14_3tile), size=1000, replace=TRUE)),
                          resist=resistance_14_3tile,
                          start.resistance=1)

# Save path length --------------------------------------------------------

traj_list <- vector("list", 1000)

for (i in seq_len(1000)) {
  cols <- ((i-1)*3+1):(i*3)
  limpet_data <- tracks_14_3tile[, cols]
  traj_list[[i]] <- sampleMovement(limpet_data, resolution=1, resist=resistance_14_3tile)
}

all_steplengths <- list()
all_resistance <- list()
limpet_id <- list()

for (i in seq_along(traj_list)) {
  n <- length(traj_list[[i]][["stats"]][["steplengths"]])
  all_steplengths[[i]] <- traj_list[[i]][["stats"]][["steplengths"]]
  all_resistance[[i]] <- traj_list[[i]][["stats"]][["resistance"]]
  limpet_id[[i]] <- rep(i, n)
}

steplengths_vec <- unlist(all_steplengths)
resistance_vec <- unlist(all_resistance)
limpet_vec <- factor(unlist(limpet_id))

df <- data.frame(
  steplength = steplengths_vec,
  resistance = resistance_vec,
  limpet = limpet_vec
)

steplength <- df %>%
  group_by(limpet) %>%
  summarise(path_length=sum(steplength, na.rm=TRUE))

# Track plot --------------------------------------------------------------

tracks_14_3tile <- tracks_14_3tile[, -seq(3, ncol(tracks_14_3tile), by=3)]

tracks_14_3tile <- lapply(1:1000, function(i) {
  df <- data.frame(
    x = tracks_14_3tile[, 2*i - 1],
    y = tracks_14_3tile[, 2*i]
  )
  return(df)
})

plot(flat, col="azure2")
plot(ripple, col='lightyellow', add=TRUE)
plot(cm_25, col='lightblue', add=TRUE)
plot(cm_5, col='lightpink', add=TRUE)
cols <- magma(1000)
for (i in 1:1000){
  lines(tracks_14_3tile[[i]], col=cols[i], lty=2, lwd=1)
}

# Plot a heat map ---------------------------------------------------------

tracks_14_3tile_df <- as.data.frame(do.call(rbind, tracks_14_3tile))

ggplot()+
  stat_density2d(tracks_14_3tile_df, mapping=aes(x=x, y=y, fill=..density..), geom='raster', contour=FALSE)+
  scale_fill_viridis(option="magma", direction=-1, limits=c(0, 0.5), na.value="white", expand=FALSE, oob=squish) +
  geom_polygon(grid, mapping=aes(x=long, y=lat, group=group), colour="white", fill=NA)

# Save the R environment --------------------------------------------------

# save.image(file='./R environments/3tile_Random_14percent.RData')
# load('./R environments/3tile_Random_14percent.RData')

# -------------------------------------------------------------------------
# -------------------------- ~33% coverage, 3 panel -----------------------
# -------------------------------------------------------------------------
# Create a 8 x 16 grid ----------------------------------------------------

grid <- raster(xmn=0, xmx=4, ymn=0, ymx=2, ncol=16, nrow=8)
grid <- rasterToPolygons(grid)
plot(grid)

ID <- paste0(1:length(grid))
for (i in 1:length(slot(grid, "polygons"))) {
  slot(slot(grid, "polygons")[[i]], "ID") <- ID[i]
}
slot(grid, "data") <- cbind("ID" = 1:length(grid), slot(grid, "data"))

# Sample the grid to create tile sets -------------------------------------

gridsample1 <- grid[sample(length(ID), 14),] 
projection(gridsample1) <- projection(grid)
grid <- grid-gridsample1
plot(grid)
plot(gridsample1, col="lightgreen", add=T)

gridsample2 <- grid[sample(length(grid), 14),] 
projection(gridsample2) <- projection(grid)
grid <- grid-gridsample2
plot(grid)
plot(gridsample1, col="lightgreen", add=T)
plot(gridsample2, col="lightblue", add=T)

gridsample3 <- grid[sample(length(grid), 14),] 
projection(gridsample3) <- projection(grid)
grid <- grid-gridsample3

gridsample4 <- grid
projection(gridsample4) <- projection(grid)

flat <- gridsample4 
ripple <- gridsample1
cm_25 <- gridsample2
cm_5 <- gridsample3

plot(flat, col="white")
plot(ripple, col='lightgreen', add=TRUE)
plot(cm_25, col='lightblue', add=TRUE)
plot(cm_5, col='orange', add=TRUE)

# Add resistance values ---------------------------------------------------

flat@data$resistance <- runif(length(flat@data$ID), 0.95096189, 1)
ripple@data$resistance <- runif(length(ripple@data$ID), 0, 0.05446582)
cm_25@data$resistance <- runif(length(cm_25@data$ID), 0.68657211, 0.75614646)
cm_5@data$resistance <- runif(length(cm_5@data$ID), 0.95783566, 1) 

crs(ripple) <- crs(flat) 
crs(cm_25) <- crs(ripple)
crs(cm_5) <- crs(cm_25)
grid_33_3tile <- rbind(flat, ripple, cm_25, cm_5)

# Create a resistance raster ----------------------------------------------

flat_val <- flat@data$resistance
ripple_val <- ripple@data$resistance
cm_25_val <- cm_25@data$resistance
cm_5_val <- cm_5@data$resistance

mapvalues <- as.data.frame(cbind(flat_val, ripple_val, cm_25_val, cm_5_val))
mapvalues <- data.frame(x = unlist(mapvalues))
row.names(mapvalues) <- NULL
mapvalues$x2 <- as.character(mapvalues$x)
mapvalues <- mapvalues[!duplicated(mapvalues),]

resistance_33_3tile <- resistanceFromShape(grid_33_3tile,   
                                           res=0.25,      
                                           field="resistance",  
                                           mapvalues=with(mapvalues, setNames(x, x2)),     
                                           background=1, margin=2)      
plot(resistance_33_3tile)

# Simulate limpet movement ------------------------------------------------

species.list <- lapply(1:1000, function(i) species(state.CRW(0.98)+0.44))

set.seed(123)  
tracks_33_3tile <- simulate(species.list,
                            time = 288, 
                            coords = xyFromCell(resistance_33_3tile, sample(1:ncell(resistance_33_3tile), size=1000, replace=TRUE)),
                            resist=resistance_33_3tile,
                            start.resistance=1)

# Save path length --------------------------------------------------------

traj_list <- vector("list", 1000)

for (i in seq_len(1000)) {
  cols <- ((i-1)*3+1):(i*3)
  limpet_data <- tracks_33_3tile[, cols]
  traj_list[[i]] <- sampleMovement(limpet_data, resolution=1, resist=resistance_33_3tile)
}

all_steplengths <- list()
all_resistance <- list()
limpet_id <- list()

for (i in seq_along(traj_list)) {
  n <- length(traj_list[[i]][["stats"]][["steplengths"]])
  all_steplengths[[i]] <- traj_list[[i]][["stats"]][["steplengths"]]
  all_resistance[[i]] <- traj_list[[i]][["stats"]][["resistance"]]
  limpet_id[[i]] <- rep(i, n)
}

steplengths_vec <- unlist(all_steplengths)
resistance_vec <- unlist(all_resistance)
limpet_vec <- factor(unlist(limpet_id))

df <- data.frame(
  steplength = steplengths_vec,
  resistance = resistance_vec,
  limpet = limpet_vec
)

steplength <- df %>%
  group_by(limpet) %>%
  summarise(path_length=sum(steplength, na.rm=TRUE))

# Track plot --------------------------------------------------------------

tracks_33_3tile <- tracks_33_3tile[, -seq(3, ncol(tracks_33_3tile), by=3)]

tracks_33_3tile <- lapply(1:1000, function(i) {
  df <- data.frame(
    x = tracks_33_3tile[, 2*i - 1],
    y = tracks_33_3tile[, 2*i]
  )
  return(df)
})

plot(flat, col="azure2")
plot(ripple, col='lightyellow', add=TRUE)
plot(cm_25, col='lightblue', add=TRUE)
plot(cm_5, col='lightpink', add=TRUE)
cols <- magma(1000)
for (i in 1:1000){
  lines(tracks_33_3tile[[i]], col=cols[i], lty=2, lwd=1)
}

# Plot a heat map ---------------------------------------------------------

tracks_33_3tile_df <- as.data.frame(do.call(rbind, tracks_33_3tile))

ggplot()+
  stat_density2d(tracks_33_3tile_df, mapping=aes(x=x, y=y, fill=..density..), geom='raster', contour=FALSE)+
  scale_fill_viridis(option="magma", direction=-1, limits=c(0, 0.5), na.value="white", expand=FALSE, oob=squish) +
  geom_polygon(grid, mapping=aes(x=long, y=lat, group=group), colour="white", fill=NA)

# Save the R environment --------------------------------------------------

# save.image(file='./R environments/3tile_Random_33percent.RData')
# load('./R environments/3tile_Random_33percent.RData')

# -------------------------------------------------------------------------
# -------------------------- ~50% coverage, 3 panel -----------------------
# -------------------------------------------------------------------------
# Create a 8 x 16 grid ----------------------------------------------------

grid <- raster(xmn=0, xmx=4, ymn=0, ymx=2, ncol=16, nrow=8)
grid <- rasterToPolygons(grid)
plot(grid)

ID <- paste0(1:length(grid))
for (i in 1:length(slot(grid, "polygons"))) {
  slot(slot(grid, "polygons")[[i]], "ID") <- ID[i]
}
slot(grid, "data") <- cbind("ID" = 1:length(grid), slot(grid, "data"))

# Sample the grid to create tile sets -------------------------------------

gridsample1 <- grid[sample(length(ID), 21),] 
projection(gridsample1) <- projection(grid)
grid <- grid-gridsample1
plot(grid)
plot(gridsample1, col="lightgreen", add=T)

gridsample2 <- grid[sample(length(grid), 21),] 
projection(gridsample2) <- projection(grid)
grid <- grid-gridsample2
plot(grid)
plot(gridsample1, col="lightgreen", add=T)
plot(gridsample2, col="lightblue", add=T)

gridsample3 <- grid[sample(length(grid), 21),] 
projection(gridsample3) <- projection(grid)
grid <- grid-gridsample3

gridsample4 <- grid
projection(gridsample4) <- projection(grid)

flat <- gridsample4 
ripple <- gridsample1
cm_25 <- gridsample2
cm_5 <- gridsample3

plot(flat, col="white")
plot(ripple, col='lightgreen', add=TRUE)
plot(cm_25, col='lightblue', add=TRUE)
plot(cm_5, col='orange', add=TRUE)

# Add resistance values ---------------------------------------------------

flat@data$resistance <- runif(length(flat@data$ID), 0.95096189, 1)
ripple@data$resistance <- runif(length(ripple@data$ID), 0, 0.05446582)
cm_25@data$resistance <- runif(length(cm_25@data$ID), 0.68657211, 0.75614646)
cm_5@data$resistance <- runif(length(cm_5@data$ID), 0.95783566, 1) 

crs(ripple) <- crs(flat) 
crs(cm_25) <- crs(ripple)
crs(cm_5) <- crs(cm_25)
grid_50_3tile <- rbind(flat, ripple, cm_25, cm_5)

# Create a resistance raster ----------------------------------------------

flat_val <- flat@data$resistance
ripple_val <- ripple@data$resistance
cm_25_val <- cm_25@data$resistance
cm_5_val <- cm_5@data$resistance

mapvalues <- as.data.frame(cbind(flat_val, ripple_val, cm_25_val, cm_5_val))
mapvalues <- data.frame(x = unlist(mapvalues))
row.names(mapvalues) <- NULL
mapvalues$x2 <- as.character(mapvalues$x)
mapvalues <- mapvalues[!duplicated(mapvalues),]

resistance_50_3tile <- resistanceFromShape(grid_50_3tile,   
                                           res=0.25,     
                                           field="resistance",  
                                           mapvalues=with(mapvalues, setNames(x, x2)),     
                                           background=1, margin=2)       
plot(resistance_50_3tile)

# Simulate limpet movement ------------------------------------------------

species.list <- lapply(1:1000, function(i) species(state.CRW(0.98)+0.44))

set.seed(123)  
tracks_50_3tile <- simulate(species.list,
                            time = 288, 
                            coords = xyFromCell(resistance_50_3tile, sample(1:ncell(resistance_50_3tile), size=1000, replace=TRUE)),
                            resist=resistance_50_3tile,
                            start.resistance=1)

# Save path length --------------------------------------------------------

traj_list <- vector("list", 1000)

for (i in seq_len(1000)) {
  cols <- ((i-1)*3+1):(i*3)
  limpet_data <- tracks_50_3tile[, cols]
  traj_list[[i]] <- sampleMovement(limpet_data, resolution=1, resist=resistance_50_3tile)
}

all_steplengths <- list()
all_resistance <- list()
limpet_id <- list()

for (i in seq_along(traj_list)) {
  n <- length(traj_list[[i]][["stats"]][["steplengths"]])
  all_steplengths[[i]] <- traj_list[[i]][["stats"]][["steplengths"]]
  all_resistance[[i]] <- traj_list[[i]][["stats"]][["resistance"]]
  limpet_id[[i]] <- rep(i, n)
}

steplengths_vec <- unlist(all_steplengths)
resistance_vec <- unlist(all_resistance)
limpet_vec <- factor(unlist(limpet_id))

df <- data.frame(
  steplength = steplengths_vec,
  resistance = resistance_vec,
  limpet = limpet_vec
)

steplength <- df %>%
  group_by(limpet) %>%
  summarise(path_length=sum(steplength, na.rm=TRUE))

# Track plot --------------------------------------------------------------

tracks_50_3tile <- tracks_50_3tile[, -seq(3, ncol(tracks_50_3tile), by=3)]

tracks_50_3tile <- lapply(1:1000, function(i) {
  df <- data.frame(
    x = tracks_50_3tile[, 2*i - 1],
    y = tracks_50_3tile[, 2*i]
  )
  return(df)
})

plot(flat, col="azure2")
plot(ripple, col='lightyellow', add=TRUE)
plot(cm_25, col='lightblue', add=TRUE)
plot(cm_5, col='lightpink', add=TRUE)
cols <- magma(1000)
for (i in 1:1000){
  lines(tracks_50_3tile[[i]], col=cols[i], lty=2, lwd=1)
}

# Plot a heat map ---------------------------------------------------------

tracks_50_3tile_df <- as.data.frame(do.call(rbind, tracks_50_3tile))

ggplot()+
  stat_density2d(tracks_50_3tile_df, mapping=aes(x=x, y=y, fill=..density..), geom='raster', contour=FALSE)+
  scale_fill_viridis(option="magma", direction=-1, limits=c(0, 0.5), na.value="white", expand=FALSE, oob=squish) +
  geom_polygon(grid, mapping=aes(x=long, y=lat, group=group), colour="white", fill=NA)

# Save the R environment --------------------------------------------------

# save.image(file='./R environments/3tile_Random_50percent.RData')
# load('./R environments/3tile_Random_50percent.RData')
