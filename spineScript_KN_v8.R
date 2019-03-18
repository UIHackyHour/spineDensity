
library(dplyr)          # to install these , enter install.packages("package name") into the console in RStudio. replace package name with each package name in quotes, ie install.packages("dplyr")
library(ggplot2)
library(lattice)
library(latticeExtra)
library(cluster)
library(stats)
library(readr)


fileChosen <- file.choose()     # opens file dialog to open 
filePath <- dirname(fileChosen) # gets the directory name that that file is in
setwd(filePath)           # sets the working directory. important for saving files later
fileList <- list.files(filePath)  # a list of files in the directory to go through
fileList <- grep(".csv", fileList, value = TRUE) # finds only lists .csv files in case there are other types
data_all <- data.frame()  # initializes a data.frame that will store all of the data for the experiment

for(i in 1:length(fileList)){  ##add back in 'i' when adding for loop back in
  fileName <- unlist(strsplit(fileList[i], "[.]"))[1]   # removes ".csv" from file name
  fileLoc <- file.path(filePath, fileList[i])           # creates a path to each file (differs from FileChosen as the loop runs)
  df <- read_csv(fileLoc, col_types = cols()) # reads the csv file, suppresses output information about columns
  colnames(df) <- gsub("-", "_", colnames(df)) # R doesn't like dashes in column names, replaces "-" with "_" for all colnames of df
  df <- df[complete.cases(df$SOMA_DISTANCE),] # removes cases where there is no soma distance data
  df <- df[complete.cases(df$RAYBURST_VOLUME),] # removes cases where there is no data
  df <- df[complete.cases(df$MAX_DTS),] # removes cases where there is no  data
  df$file <- fileName # adds a column with repeated information about which file that spine came from
  total_length <- df %>% group_by(SECTION_NUMBER)  %>% summarise(section_length=max(SECTION_LENGTH)) %>%ungroup() %>% summarise(total_length=sum(section_length)) %>% as.double() 
  # ^^this uses the dplyr package to take the data.frame and group it by section. Then it finds the length of each section before adding them all together
  total_spines <- nrow(df) # number of rows = number of spines
  density_overall <- total_spines/total_length  # calculate density of segment
  density_mushroom <- sum(as.numeric(df$TYPE=="mushroom"))/total_length  # counts the number of "mushroom" rows and divides by length for mushroom density
  density_thin <- sum(as.numeric(df$TYPE=="thin"))/total_length  # same as above, with "thin"
  density_stubby <- sum(as.numeric(df$TYPE=="stubby"))/total_length # same as above, with "stubby"
  agn <- agnes(df[c("X", "Y", "Z")], metric = "euclidean", method = "average") # runs the agnes, computes agglomerative hierarchical clustering, "average" = UPGMA
  dist <- as.matrix(dist(df[c("X", "Y", "Z")])) # creates a distance matrix of soma distance
  df$nn_dist <- apply(dist,2, function(x) sort(x)[2])   # finds nearest neighbor for each spine
  df$nn2_dist <- apply(dist,2, function(x) sort(x)[3])  # finds second nearest neighbor
  df$nn3_dist <- apply(dist,2, function(x) sort(x)[4])  # finds third
  df$exp_nn <- 1/density_overall  # calculates the expected nearest neighbor based on overall density
  df$nn_clustered <- apply(dist,2, function(x) sum(as.numeric(x<df$exp_nn))) # how many spines are within the expected distance
  df$is_it_2clustered <- as.numeric(df$nn_dist<df$exp_nn) # if nearest neighbor is closer than expected neighbor (given the density), 1 = TRUE, 0 = FALSE
  df$is_it_3clustered <- as.numeric(df$nn2_dist<df$exp_nn)
  df$is_it_4clustered <- as.numeric(df$nn3_dist<df$exp_nn)
  dendrite_ac <- agn$ac # defines agglomeration coefficent of the dendrite

  ac <- list() #initializes list for ac's from next for loop
  for(j in 1:15000){
    test_data <- data.frame(sample(df$X), sample(df$Y), sample(df$Z)) # randomize the X's, Y's, and Z's to make a "biologically plausible" dataframe.
    test_dist <- as.matrix(dist(test_data)) #creates distance matrix for random sample
    assign("test_cluster", agnes(test_data, metric = "euclidean", method = "average")) #runs UPGMA on sample and labels the agn output "test_cluster"
    ac <- rbind(ac,test_cluster$ac) # add row to ac using the 'test_cluster' ac
  }
  # mean density of the simulation
  m <- mean(unlist(ac))
  # standard deviation of the simulation
  s <- sd(ac)
  # calculated zscore (assuming simulation data are normally distributed which appears true)
  df$zscore <- (dendrite_ac - m) / s
  cScore <- sum(as.numeric(ac<dendrite_ac))/15000 # average (divide by 15000 samples) how many times random ac is smaller than dendrite_ac, the smaller the value, the more "clustering"
  df$c_score <- cScore #add a row to df with the cScore of the dendrite
  df$density_overall <- density_overall # add these data to df
  df$density_mushroom <- density_mushroom
  df$density_thin <- density_thin
  df$density_stubby <- density_stubby
  data_all <- rbind(data_all, df) # adds df as next row in the data_all file
} # end of for loop, all files in folder should be added to data_all

data_all$animal_num <- lapply(data_all$file, function(x) unlist(strsplit(x, "-"))[2]) # this pulls out an animal number from the file number
data_all$retro_label <- substring(data_all$animal_num, nchar(data_all$animal_num), nchar(data_all$animal_num))  #pulls off 'L' or 'N' from ID to indicate whether cell was retro-gradely labeled
#data_all$PDB <- substring(data_all$animal_num, nchar(data_all$animal_num), nchar(data_all$animal_num))
data_all$PDB <- lapply(data_all$file, function(x) unlist(strsplit(x, "-"))[3]) # pulls dendrite location out of name and adds column
data_all$PDB <- replace(data_all$PDB, data_all$PDB=="p", "prox")
data_all$PDB <- replace(data_all$PDB, data_all$PDB=="d", "dist")
data_all$PDB <- replace(data_all$PDB, data_all$PDB=="b", "basal")
data_all$retro_label <- replace(data_all$retro_label, data_all$retro_label=="L", "labeled") 
data_all$retro_label <- replace(data_all$retro_label, data_all$retro_label=="N", "not labeled")
data_all$stack <- lapply(data_all$file, function(x) unlist(strsplit(x, "-"))[4]) #gives letter ID of different imaging days
data_all$stack <- unlist(data_all$stack)
data_all$animal_num <- lapply(data_all$animal_num, function(x) unlist(strsplit(x, "L"))[1]) #removes letter from behind animal ID name
data_all$animal_num <- lapply(data_all$animal_num, function(x) unlist(strsplit(x, "N"))[1])
data_all$RAYBURST_VOLUME <- as.numeric(data_all$RAYBURST_VOLUME)
data_all$MAX_DTS <- as.numeric(data_all$MAX_DTS)

n <- readline(prompt="Enter group name:") # gives prompt on console to enter animal group
data_all$group <- n # adds the group to the master data file

data_all$animal_num <- unlist(data_all$animal_num) #changes animal number to a vector rather than a list, which is important for executing the following task
data_all$PDB <- unlist(data_all$PDB) # same as above

TYPE_ave <- data_all %>% group_by(group, animal_num, retro_label, TYPE) %>% summarise(Cscore = mean(c_score), zscore = mean(zscore), nn = mean(nn_dist), nn2 = mean(nn2_dist), nn3 = mean(nn3_dist), density_overall = mean(density_overall), density_mushroom = mean(density_mushroom), density_thin = mean(density_thin), density_stubby = mean(density_stubby), spine_vol_overall = mean(RAYBURST_VOLUME, na.rm = TRUE), spine_length_overall = mean(MAX_DTS, na.rm = TRUE))
PDB_TYPE_ave <- data_all %>% group_by(group, animal_num, stack, retro_label, PDB, TYPE) %>% summarise(Cscore = mean(c_score), zscore = mean(zscore), nn = mean(nn_dist), nn2 = mean(nn2_dist), nn3 = mean(nn3_dist), density_overall = mean(density_overall), density_mushroom = mean(density_mushroom), density_thin = mean(density_thin), density_stubby = mean(density_stubby), spine_vol_overall = mean(RAYBURST_VOLUME, na.rm = TRUE), spine_length_overall = mean(MAX_DTS, na.rm = TRUE))

PDB_ave <- data_all %>% group_by(group, animal_num, retro_label, PDB) %>% summarise(Cscore = mean(c_score), zscore = mean(zscore), nn = mean(nn_dist), nn2 = mean(nn2_dist), nn3 = mean(nn3_dist), density_overall = mean(density_overall), density_mushroom = mean(density_mushroom), density_thin = mean(density_thin), density_stubby = mean(density_stubby), spine_vol_overall = mean(RAYBURST_VOLUME, na.rm = TRUE), spine_length_overall = mean(MAX_DTS, na.rm = TRUE))
# ^first groups by treatment group, then by animal, then by day of image (i.e. -F), then by dendrite location
# summarise then finds averages for each group for all various pieces of data

animal_ave <- data_all %>% group_by(group, animal_num) %>% summarise(Cscore = mean(c_score), zscore = mean(zscore), nn = mean(nn_dist), nn2 = mean(nn2_dist), nn3 = mean(nn3_dist), density_overall = mean(density_overall), density_mushroom = mean(density_mushroom), density_thin = mean(density_thin), density_stubby = mean(density_stubby), spine_vol_overall = mean(RAYBURST_VOLUME, na.rm = TRUE), spine_length_overall = mean(MAX_DTS, na.rm = TRUE))
# groups by treatment then animal, then finds averages per animal for various data

labeled_ave <- data_all %>% group_by(group, animal_num, retro_label) %>% summarise(Cscore = mean(c_score), zscore = mean(zscore), nn = mean(nn_dist), nn2 = mean(nn2_dist), nn3 = mean(nn3_dist), density_overall = mean(density_overall), density_mushroom = mean(density_mushroom), density_thin = mean(density_thin), density_stubby = mean(density_stubby), spine_vol_overall = mean(RAYBURST_VOLUME, na.rm = TRUE), spine_length_overall = mean(MAX_DTS, na.rm = TRUE))
# groups by treatment then animal, then finds averages per animal for various data



dir.create("analysis")  #creates a directory to create a file on the computer
savePath <- paste(filePath,"/", "analysis", collapse = "/", sep="") # creates the path where files can be saved
setwd(savePath) # sets working directory to the path created above
write.csv(data_all, "data_all.csv")
write.csv(PDB_ave, "PDB averages.csv")
write.csv(animal_ave, "Animal Averages.csv")
write.csv(TYPE_ave, "spine_type_averages.csv")
write.csv(PDB_TYPE_ave, "PDB_spine_type_averages.csv")
write.csv(labeled_ave, "labeled_averages.csv")

