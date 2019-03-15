# modeling Yadav et al. 2012 cluster analysis; tailored towards clustering of attachment points
# initially coded by S. Johnson- updates and 3D analysis of spine heads added by K. Nett


library(dplyr)          # to install these , enter install.packages("package name") into the console in RStudio. replace package name with each package name in quotes, ie install.packages("dplyr")
library(ggplot2)
library(lattice)
library(latticeExtra)
library(cluster)
library(stats)
library(readr)
library(bio3d) #install.packages("bio3d", dependencies = TRUE)
library(DECIPHER)
#how to install DECIPHER - a is enter in response to prompt for "all"
#source("https://bioconductor.org/biocLite.R")
#biocLite("DECIPHER")
#a


fileChosen <- file.choose()     # opens file dialog to open 
filePath <- dirname(fileChosen) # gets the directory name that that file is in
setwd(filePath)           # sets the working directory. important for saving files later
fileList <- list.files(filePath)  # a list of files in the directory to go through
fileList <- grep(".csv", fileList, value = TRUE) # finds only lists .csv files in case there are other types
data_all <- data.frame()  # initializes a data.frame that will store all of the data for the experiment




for(i in 1:length(fileList)){  ##add back in 'i' when adding for loop back in
  fileName <- unlist(strsplit(fileList[1], "[.]"))[1]   # removes ".csv" from file name
  fileLoc <- file.path(filePath, fileList[1])           # creates a path to each file (differs from FileChosen as the loop runs)
  df <- read_csv(fileLoc, col_types = cols()) # reads the csv file, suppresses output information about columns
  colnames(df) <- gsub("-", "_", colnames(df)) # R doesn't like dashes in column names, replaces "-" with "_" for all colnames of df
  df <- df[complete.cases(df$SOMA_DISTANCE),] # removes cases where there is no soma distance data
  df <- df[complete.cases(df$RAYBURST_VOLUME),] # removes cases where there is no data
  df <- df[complete.cases(df$MAX_DTS),] # removes cases where there is no  data
  df$file <- fileName # adds a column with repeated information about which file that spine came from
  total_length <- df %>% group_by(SECTION_NUMBER)  %>% summarise(section_length=max(SECTION_LENGTH)) %>%ungroup() %>% summarise(total_length=sum(section_length)) %>% as.double() 
  # ^^this uses the dplyr package to take the data.frame and group it by section. Then it finds the length of each section before adding them all together
  
  total_spines <- as.numeric(nrow(df)) # number of rows = number of spines
  density_overall <- total_spines/total_length  # calculate density of segment
  density_mushroom <- sum(as.numeric(df$TYPE=="mushroom"))/total_length  # counts the number of "mushroom" rows and divides by length for mushroom density
  density_thin <- sum(as.numeric(df$TYPE=="thin"))/total_length  # same as above, with "thin"
  density_stubby <- sum(as.numeric(df$TYPE=="stubby"))/total_length # same as above, with "stubby"
  df <- df[order(df$SOMA_DISTANCE),]  #order data by column SOMA_DISTANCE so that coordinates pulled are in spine attachment order
  
  data_coords <- data.frame(df$X, df$Y, df$Z) # creates data frame with spine head coordinates

# Shane's original clustering analysis using spine attachment distance from soma  
    
  agn <- agnes(df$SOMA_DISTANCE,metric = "euclidean", method = "average") # runs the agnes, computes agglomerative hierarchical clustering, "average" = UPGMA
  dist_ac <- as.matrix(dist(df$SOMA_DISTANCE))  #creates a distance matrix of soma distances
  df$nn_dist_ac <- apply(dist_ac,2, function(x) sort(x)[2])   # finds nearest neighbor for each spine
  df$nn2_dist_ac <- apply(dist_ac,2, function(x) sort(x)[3])  # finds second nearest neighbor
  df$nn3_dist_ac <- apply(dist_ac,2, function(x) sort(x)[4])  # finds third
  dendrite_ac <- agn$ac # defines agglomeration coefficent of the dendrite
  
  ac_test <- list() #initializes list for ac's from next for loop  
  possible_dist <- seq(1,total_length, by=0.01) # create list of possible soma distances by 0.01 increments to pull from
  
   for(k in 1:5000){
    test_data <- sample(possible_dist, total_spines) # take a sample from all possible locations on dendrite to match total number of spines
    test_dist <- as.matrix(dist(test_data)) #creates distance matrix for random sample
    assign("test_cluster", agnes(test_data,metric = "euclidean", method = "average")) #runs UPGMA on sample and labels the agn output "test_cluster"
    ac_test <- rbind(ac_test,test_cluster$ac) # add row to ac using the 'test_cluster' ac
  } 
  
  cScore_ac_1D <- sum(as.numeric(ac_test<dendrite_ac))/5000 # average (divide by 1000 samples) how many times random ac is smaller than dendrite_ac, the smaller the value, the more "clustering"
  df$c_score_ac_1D <- cScore_ac_1D #add a row to df with the cScore of the dendrite
#end Shane's original analysis
  

  
#1D clustering (from spine attachment points)
  dist_1D <- as.matrix(dist(df$SOMA_DISTANCE)) # creates a distance matrix using distance from soma, same as dist_ac, but keeping separate for different analysis for now
  UPGMA_obsv_1D <- IdClusters(dist_1D, method = "UPGMA", cutoff=0.75, showPlot=TRUE)  #runs cluster analysis with cutoff from Yadav paper on observed spines
  UPGMA_obsv_1D$rn <- rownames(UPGMA_obsv_1D) # adds a column with row names to keep spine ID, not sure if this step is necessary
  cluster_all_1D <- cbind(UPGMA_obsv_1D, df$SOMA_DISTANCE)  # adds cluster number, spine ID(row number) and corresponding distance from soma together
  cluster_all_1D <- cluster_all_1D[order(cluster_all_1D$cluster),] #orders data by cluster number
  cluster_freq_1D <- table(cluster_all_1D$cluster) #creates a table with how many spines are in each cluster
  cluster_freq_1D <- as.data.frame(cluster_freq_1D) #converts above table to a data frame
  cluster_freq_1D$Freq <- as.numeric(cluster_freq_1D$Freq) # turns the "cluster" column to a number
  cluster_freq_1D$is_clustered <- as.numeric(cluster_freq_1D$Freq > 1) #returns 1 if there is >1 spine in a cluster and 0 if not-- therefore all 1s reflect a true "cluster" since it has more than 1 spine in ity 
  num_clusters_1D <- sum(cluster_freq_1D$is_clustered) #count how many clusters are on this segment
  spines_in_cluster_1D <- cluster_freq_1D %>% group_by(is_clustered) %>% summarise(num_clus_spines_1D = sum(Freq)) # count how many spines are in the true clusters (>1 spines) and how many are not
  spines_clustered_1D <- list() #initialize list for number of spines that are clustered
  spines_clustered_1D <- rbind(spines_clustered_1D, spines_in_cluster_1D[1,2]) #add number of spines clustered to list
  spines_clustered_1D <- as.matrix(spines_clustered_1D) #converts # of spines clustered to matrix
  spines_clustered_1D <- as.numeric(spines_clustered_1D) #converts to numerical form
  
  spines_not_1D <- as.numeric(total_spines - spines_clustered_1D) #define and calculate number of spines not clustered (total spines minus number of spines clustered)
  
  test_spines_clustered_all_1D <- data.frame() #initialize dataframe to contain # of spines in a cluster for all random samples
  test_num_clusters_all_1D <- data.frame()
  
# 1D random spines for-loop
  for(j in 1:500){
    random_spines_1D <- sample(possible_dist, total_spines) # take a sample from all possible locations on dendrite to match total number of spines
    test_dist_1D <- as.matrix(dist(random_spines_1D)) #creates distance matrix for random sample
    UPGMA_test_1D <- IdClusters(test_dist_1D, method="UPGMA", cutoff=0.75, showPlot=FALSE) #runs UPGMA with cutoff (same as with obsv but with random sample)
    UPGMA_test_1D$rn <- rownames(UPGMA_test_1D) #add row name (spine ID) to cluster number  
    cluster_all_test_1D <- cbind(UPGMA_test_1D, random_spines_1D) #add cluster number/spine ID to randomly generated soma distances
    cluster_all_test_1D <- cluster_all_test_1D[order(cluster_all_test_1D$cluster),] #order by cluster number
    cluster_freq_test_1D <- table(cluster_all_test_1D$cluster) #generate table with number of spines in each cluster  
    cluster_freq_test_1D <- as.data.frame(cluster_freq_test_1D) #change table to data frame 
    cluster_freq_test_1D$is_clustered <- as.numeric(cluster_freq_test_1D$Freq > 1) #ask whether cluster has >1 spines in it (1 for yes, 0 for no)
    num_clusters_test_1D <- sum(cluster_freq_test_1D$is_clustered) # add how many 1s (or how many clusters have >1 spines) to get true number of clusters
    num_clusters_test_1D[is.na(num_clusters_test_1D)] <- 0 #changes possible NA from 0 clusters to 0
    test_num_clusters_all_1D <- rbind(test_num_clusters_all_1D, num_clusters_test_1D) #adds total number of clusters to running list
    spines_in_cluster_test_1D <- cluster_freq_test_1D %>% group_by(is_clustered) %>% summarise(num_clusters_test_1D = sum(Freq)) #count how many spines that are in a true cluster or not
    spines_clustered_test_1D <- spines_in_cluster_test_1D[2,2] # define how many spines are in a cluster
    spines_not_test_1D <- as.numeric(total_spines - spines_clustered_test_1D) # calculate how many spines are not clustered
    spines_clustered_test_1D[is.na(spines_clustered_test_1D)] <- 0  #returns 0 instead of Na if no spines are clustered in the random sample
    test_spines_clustered_all_1D <- rbind(test_spines_clustered_all_1D, spines_clustered_test_1D) #add number of clustered spines to running list 
    
  } # end of random spines 1D for-loop
  
  test_spines_clustered_all_1D <- as.matrix(test_spines_clustered_all_1D) #changes data.frame to matrix (not sure if this is neccesary but it works)
  test_spines_clustered_all_1D <- as.numeric(test_spines_clustered_all_1D) #changes all vales to numeric (again, not sure if neccessary)
  
  
  std_test_1D <- sd(test_spines_clustered_all_1D) #generates STD of total number of spines clustered in entire random sample
  mean_test_1D <- mean(test_spines_clustered_all_1D) #generates average number of spines in a cluster in random data
  curve_dnorm_1D <- dnorm(test_spines_clustered_all_1D, mean_test_1D, std_test_1D) #gives probability density function, or height of probability distribution at each point(height = frequency)
  std_curve_1D <- sd(curve_dnorm_1D)
  mean_curve_1D <- mean(curve_dnorm_1D)
  
  Cscore_1D <- pnorm(spines_clustered_1D, mean_curve_1D, sd=sqrt(std_curve_1D)) #calculates the probability that given the random sample distribution (curve mean+SD) the total number of spines observed is higher
  #therefore, the closer to 1, higher probability that a given number has more clustered spines than the random normal distribution and vice versa
  
 # add all 1D data to all data 
  cluster_data_1D <- data.frame()  
  cluster_data_1D <- rbind(cluster_data_1D, num_clusters_1D)
  cluster_data_1D <- cbind(cluster_data_1D, spines_clustered_1D)
  cluster_data_1D <- cbind(cluster_data_1D, spines_not_1D)
  cluster_data_1D <- cbind(cluster_data_1D, Cscore_1D)
  cluster_data_1D[is.na(cluster_data_1D)] <- 0
  colnames(cluster_data_1D) <- c("# of clusters - 1D", "spines clustered - 1D", "spines not clustered - 1D", "Cscore - 1D")
  
  
# 3D clustering analysis  
  dist_3D <- as.matrix(dist.xyz(data_coords)) # creates distance matrix of spine head coordinates
  UPGMA_obsv_3D <- IdClusters(dist_3D, method = "UPGMA", cutoff=0.75, showPlot=TRUE) #run cluster analysis with cutoff used in Yadav paper
  #gives cluster number associated with which spine (i.e. spines 25 and 26 are in cluster 1)
  UPGMA_obsv_3D$rn <- rownames(UPGMA_obsv_3D)   #adds a column with rownows to keep spine ID, not sure if this step is neccessary
  cluster_all_3D <- cbind(UPGMA_obsv_3D, data_coords)   #get the coordinates for each spine ID, if kept in row name/number order, will correctly correspond to each spine
  cluster_all_3D <- cluster_all_3D[order(cluster_all_3D$cluster),] #sort data frame by cluster number
  cluster_freq_3D <- table(cluster_all_3D$cluster) #create a table counting how many times each cluster Variable occurs (i.e. how many spines in each cluster)
  cluster_freq_3D <- as.data.frame(cluster_freq_3D)
  cluster_freq_3D$Freq <- as.numeric(cluster_freq_3D$Freq)
  cluster_freq_3D$is_clustered <- as.numeric(cluster_freq_3D$Freq >1) # create column where 1 means there is more than one spine in a cluster or 0 if just 1
  num_clusters_3D <- sum(cluster_freq_3D$is_clustered) #count how many 1s to determine how many clusters (spines > 1) in the segment
  spines_in_cluster_3D <- cluster_freq_3D %>% group_by(is_clustered) %>% summarise(num_clus_spines_3D = sum(Freq))
  spines_clustered_3D <- list()
  spines_clustered_3D <- rbind(spines_clustered_3D, spines_in_cluster_3D[1,2])
  spines_clustered_3D <- as.matrix(spines_clustered_3D) #conerts num of spines clustered to matrix
  spines_clustered_3D <- as.numeric(spines_clustered_3D) #converts to numerical form
  
  spines_not_3D <- as.numeric(total_spines - spines_clustered_3D)
  
  test_spines_clustered_all_3D <- data.frame()
  test_num_clusters_all_3D <- data.frame()

  
# 3D random spines for loop
  for(j in 1:100){
    test_data_X <- data.frame(sample(df$X), df$Y, df$Z) # randomize the X's, Y's, and Z's to make a "biologically plausible" dataframe.
    colnames(test_data_X) <- c( "x", "Y", "Z")
    test_data_Y <- data.frame(df$X, sample(df$Y), df$Z)
    colnames(test_data_Y) <- c( "x", "Y", "Z")
    test_data_Z <- data.frame(df$X, df$Y, sample(df$Z))
    colnames(test_data_Z) <- c( "x", "Y", "Z")
    test_data <- rbind(test_data_X, test_data_Y, test_data_Z)
    test_data_final <-data.frame(sample_n(test_data, total_spines))
    test_dist_3D <- as.matrix(dist(test_data_final)) #creates distance matrix for random sample

    UPGMA_test_3D <- IdClusters(test_dist_3D, method = "UPGMA", cutoff=0.75, showPlot=TRUE) #run cluster analysis with cutoff used in Yadav paper
    #gives cluster number associated with which spine (i.e. spines 25 and 26 are in cluster 1)
    UPGMA_test_3D$rn <- rownames(UPGMA_test_3D)   #adds a column with rownows to keep spine ID, not sure if this step is neccessary
    cluster_all_test_3D <- cbind(UPGMA_test_3D, test_data_final)   #get the coordinates for each spine ID, if kept in row name/number order, will correctly correspond to each spine
    cluster_all_test_3D <- cluster_all_test_3D[order(cluster_all_test_3D$cluster),] #sort data frame by cluster number
    cluster_freq_test_3D <- table(cluster_all_test_3D$cluster) #create a table counting how many times each cluster Variable occurs (i.e. how many spines in each cluster)
    cluster_freq_test_3D <- as.data.frame(cluster_freq_test_3D)
    cluster_freq_test_3D$is_clustered <- as.numeric(cluster_freq_test_3D$Freq >1) # create column where 1 means there is more than one spine in a cluster or 0 if just 1
    num_clusters_test_3D <- sum(cluster_freq_test_3D$is_clustered)
    num_clusters_test_3D[is.na(num_clusters_test_3D)] <- 0
    test_num_clusters_all_3D <- rbind(test_num_clusters_all_3D, num_clusters_test_3D)
    spines_in_cluster_test_3D <- cluster_freq_test_3D %>% group_by(is_clustered) %>% summarise(num_clus_spines_test_3D = sum(Freq))
    spines_clustered_test_3D <- spines_in_cluster_test_3D[2,2]
    
    spines_not_test_3D <- as.numeric(spines_in_cluster_test_3D[1,2])
    spines_clustered_test_3D[is.na(spines_clustered_test_3D)] <- 0
    test_spines_clustered_all_3D <- rbind(test_spines_clustered_all_3D, spines_clustered_test_3D)
    
  } # end of random spines 3D for-loop
  
  test_spines_clustered_all_3D <- as.numeric(as.matrix(test_spines_clustered_all_3D))

  
  std_test_3D <- sd(test_spines_clustered_all_3D)
  mean_test_3D <- mean(test_spines_clustered_all_3D)
  curve_dnorm_3D <- dnorm(test_spines_clustered_all_3D, mean_test_3D, std_test_3D)
  std_curve_3D <- sd(curve_dnorm_3D)
  mean_curve_3D <- mean(curve_dnorm_3D)
  
  Cscore_3D <- pnorm(spines_clustered_3D, mean_curve_3D, std_curve_3D)


  # add 3D data to data frame
  cluster_data_3D <- data.frame()
  cluster_data_3D <- rbind(cluster_data_3D, num_clusters_3D)
  cluster_data_3D <- cbind(cluster_data_3D, spines_clustered_3D)
  cluster_data_3D <- cbind(cluster_data_3D, spines_not_3D)
  cluster_data_3D <- cbind(cluster_data_3D, Cscore_3D)
  cluster_data_3D[is.na(cluster_data_3D)] <- 0
  colnames(cluster_data_3D) <- c("# of clusters-3D", "spines clustered-3D", "spines not clustered-3D", "Cscore-3D")
  
  
  
  df$density_overall <- density_overall # add these data to df
  df$density_mushroom <- density_mushroom
  df$density_thin <- density_thin
  df$density_stubby <- density_stubby
  data_all <- rbind(data_all, df) # adds df as next row in the data_all file
  
  
  data_all <- cbind(data_all, cluster_data_1D)
  data_all <- cbind(data_all, cluster_data_3D)

  
} # end of file for-loop 

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

TYPE_ave <- data_all %>% group_by(group, animal_num, retro_label, TYPE) %>% summarise(Cscore_1D = mean(c_score_1D), nn_1D = mean(nn_dist_1D), nn2_1D = mean(nn2_dist_1D), nn3_1D = mean(nn3_dist_1D), Cscore_3D = mean(c_score_3D), nn_3D = mean(nn_dist_3D), nn2_3D = mean(nn2_dist_3D), nn3_3D = mean(nn3_dist_3D), density_overall = mean(density_overall), density_mushroom = mean(density_mushroom), density_thin = mean(density_thin), density_stubby = mean(density_stubby), spine_vol_overall = mean(RAYBURST_VOLUME, na.rm = TRUE), spine_length_overall = mean(MAX_DTS, na.rm = TRUE))
PDB_TYPE_ave <- data_all %>% group_by(group, animal_num, stack, retro_label, PDB, TYPE) %>% summarise(Cscore_1D = mean(c_score_1D), nn_1D = mean(nn_dist_1D), nn2_1D = mean(nn2_dist_1D), nn3_1D = mean(nn3_dist_1D), Cscore_3D = mean(c_score_3D), nn_3D = mean(nn_dist_3D), nn2_3D = mean(nn2_dist_3D), nn3_3D = mean(nn3_dist_3D), density_overall = mean(density_overall), density_mushroom = mean(density_mushroom), density_thin = mean(density_thin), density_stubby = mean(density_stubby), spine_vol_overall = mean(RAYBURST_VOLUME, na.rm = TRUE), spine_length_overall = mean(MAX_DTS, na.rm = TRUE))

PDB_ave <- data_all %>% group_by(group, animal_num, retro_label, PDB) %>% summarise(Cscore_1D = mean(c_score_1D), nn_1D = mean(nn_dist_1D), nn2_1D = mean(nn2_dist_1D), nn3_1D = mean(nn3_dist_1D), Cscore_3D = mean(c_score_3D), nn_3D = mean(nn_dist_3D), nn2_3D = mean(nn2_dist_3D), nn3_3D = mean(nn3_dist_3D), density_overall = mean(density_overall), density_mushroom = mean(density_mushroom), density_thin = mean(density_thin), density_stubby = mean(density_stubby), spine_vol_overall = mean(RAYBURST_VOLUME, na.rm = TRUE), spine_length_overall = mean(MAX_DTS, na.rm = TRUE))
# ^first groups by treatment group, then by animal, then by day of image (i.e. -F), then by dendrite location
# summarise then finds averages for each group for all various pieces of data

animal_ave <- data_all %>% group_by(group, animal_num) %>% summarise(Cscore_1D = mean(c_score_1D), nn_1D = mean(nn_dist_1D), nn2_1D = mean(nn2_dist_1D), nn3_1D = mean(nn3_dist_1D), Cscore_3D = mean(c_score_3D), nn_3D = mean(nn_dist_3D), nn2_3D = mean(nn2_dist_3D), nn3_3D = mean(nn3_dist_3D), density_overall = mean(density_overall), density_mushroom = mean(density_mushroom), density_thin = mean(density_thin), density_stubby = mean(density_stubby), spine_vol_overall = mean(RAYBURST_VOLUME, na.rm = TRUE), spine_length_overall = mean(MAX_DTS, na.rm = TRUE))
# groups by treatment then animal, then finds averages per animal for various data

labeled_ave <- data_all %>% group_by(group, animal_num, retro_label) %>% summarise(Cscore_1D = mean(c_score_1D), nn_1D = mean(nn_dist_1D), nn2_1D = mean(nn2_dist_1D), nn3_1D = mean(nn3_dist_1D), Cscore_3D = mean(c_score_3D), nn_3D = mean(nn_dist_3D), nn2_3D = mean(nn2_dist_3D), nn3_3D = mean(nn3_dist_3D), density_overall = mean(density_overall), density_mushroom = mean(density_mushroom), density_thin = mean(density_thin), density_stubby = mean(density_stubby), spine_vol_overall = mean(RAYBURST_VOLUME, na.rm = TRUE), spine_length_overall = mean(MAX_DTS, na.rm = TRUE))
# groups by treatment then animal, then finds averages per animal for various data



#dir.create("analysis")  #creates a directory to create a file on the computer
#savePath <- paste(filePath,"/", "analysis", collapse = "/", sep="") # creates the path where files can be saved
#setwd(savePath) # sets working directory to the path created above
#write.csv(data_all, "data_all.csv")
#write.csv(PDB_ave, "PDB averages.csv")
#write.csv(animal_ave, "Animal Averages.csv")
#write.csv(TYPE_ave, "spine_type_averages.csv")
#write.csv(PDB_TYPE_ave, "PDB_spine_type_averages.csv")
#write.csv(labeled_ave, "labeled_averages.csv")

