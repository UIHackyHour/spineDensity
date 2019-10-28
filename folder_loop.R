#loop attempt with raw_data_cleanup.R

library(dplyr)           
library(ggplot2)        
library(lattice)
library(latticeExtra)
library(cluster)
library(stats)
library(readr)
library(tidyr)
library(bio3d)
library(DECIPHER)
library(svDialogs)
library(rowr)
library(tcltk)


parent.folder <- tk_choose.dir(default = "", caption = "select folder")
sub.folders <- list.dirs(parent.folder, recursive = TRUE)[-1]
setwd(parent.folder)
saveFolder <- dirname(parent.folder)

for(i in 1:length(sub.folders)) {
setwd(parent.folder)
  fileChosen <- sub.folders[i]




#for(i in 1:length(sub.folders)) {
#fileChosen <-  sub.folders[i]


# Initialize Code: Select file(s) and define group

#fileChosen <- file.choose()                                             # opens a folder dialog box to select first file in folder containing data from 1 treatment group 
filePath   <- dirname(fileChosen)                                       # gets the directory name that that file is in
setwd(filePath)                                                         # sets the working directory. important for saving files later
fileList   <- list.files(fileChosen)                                      # list the files in the directory to go through
fileList   <- grep(".csv", fileList, value = TRUE)                      # finds only lists .csv files in case there are other types
data_all   <- data.frame()                                              # initializes data frame to store data throughout the for-loop
#rat_ID     <- dlgInput("Enter Rat ID #", Sys.info()["user"])$res       # prompt given in console to enter the animal treatment group (i.e. Sal-D0)



# File list for-loop: Apply code to each file in fileList

for(i in 1:length(fileList))  {                                         # peform code inside {} for each .csv file in the file list
  fileName     <- unlist(strsplit(fileList[i], "[.]"))[1]               # remove '.csv' from the file name
  fileLoc      <- file.path(fileChosen, fileList[i])                      # create path to file (will differ from initial FileChosen as loop progresses)
  df           <- read_csv(fileLoc, col_types = cols())                 # read the csv file, uses column names from NeuronStudio file
  colnames(df) <- gsub("-", "_", colnames(df))                          # replaces dashes ('-') w/ underscores ('_') for column names (R doesn't like dashes)
  df           <- df[complete.cases(df$SOMA_DISTANCE) ,]                # removes cases where there is no soma distance data
  df           <- df[complete.cases(df$RAYBURST_VOLUME), ]              # " " spine volume data
  df           <- df[complete.cases(df$MAX_DTS), ]                      # " " length data
  df$file      <- fileName                                              # add a column identifying the file name for each spine (row)
  ## Current NeuronStudio data should only have 1 section/dendrite; however, it is possible
  ## to have multiple if they are not linked together. The next 6 lines calculate total
  ## length when there is more than 1 section
  
  total_length <- df %>%                                                # dplyr package: define total length of dendritic segment
    group_by(SECTION_NUMBER)  %>%                                       # group data frame by section number
    summarise(section_length = max(SECTION_LENGTH)) %>%                 # find the length of each section
    ungroup() %>%                                                       # ungroup previous grouping by section number
    summarise(total_length = sum(section_length)) %>%                   # sum section length of each section
    as.double()                                                         # allows for 64 bit storage (increase precision with more significant digits)
  

  df$RAYBURST_VOLUME[df$RAYBURST_VOLUME == 0]      <- NA       #if no data, make NA
  df$HEAD_DIAMETER[df$HEAD_DIAMETER == 0]          <- NA
  df$MAX_DTS[df$MAX_DTS == 0]                      <- NA
  
  df$HEAD_DIAMETER[df$HEAD_DIAMETER > 1.5]                        <- NA      #mushroom HD cutoff
  df$MAX_DTS[df$MAX_DTS > 3.00]                                   <- NA      #mushroom/thin length cutoff
  df$MAX_DTS[df$TYPE == "stubby"       & df$MAX_DTS > 0.80]       <- NA      #stubby length cutoff
  df$HEAD_DIAMETER[df$TYPE == "stubby" & df$HEAD_DIAMETER > 0.97] <- NA      #stubby HD cutoff
  df$HEAD_DIAMETER[df$TYPE == "thin"   & df$HEAD_DIAMETER > 1.22] <- NA      #thin HD cutoff

  df$RAYBURST_VOLUME[df$TYPE == "mushroom" & df$RAYBURST_VOLUME > 0.60] <- NA   #cut-offs from 2SD about mean in Harris et al. 1992 paper, pyramidal CA1 neurons
  df$RAYBURST_VOLUME[df$TYPE == "thin"     & df$RAYBURST_VOLUME > 0.10] <- NA   
  df$RAYBURST_VOLUME[df$TYPE == "stubby"   & df$RAYBURST_VOLUME > 0.05] <- NA

  total_spines     <- count(df$ID[!is.na(df$ID)])                       # number of rows = number of spines
  density_overall  <- total_spines/total_length                         # calculate density of dendritic segment
  density_mushroom <- sum(as.numeric(df$TYPE == "mushroom"),            # find total number of mushroom spines and
                          na.rm = TRUE)/total_length                    # divide by total length to find mushroom density
  density_thin     <- sum(as.numeric(df$TYPE == "thin"),                # as above for thin spines
                          na.rm = TRUE)/total_length                    # ""
  density_stubby   <- sum(as.numeric(df$TYPE == "stubby"),              # as above for stubby spines
                          na.rm = TRUE)/total_length                    # ""
  
  
  #add data to running list in data_master file
  df$density_overall  <- density_overall # add these data to df
  df$density_mushroom <- density_mushroom
  df$density_thin     <- density_thin
  df$density_stubby   <- density_stubby
  df$animal_num       <- lapply(df$file, function(x) unlist(strsplit(x, "-"))[2])             # pulls animal ID out of file name, with retro label attached
  df$retro_label      <- substring(df$animal_num, nchar(df$animal_num), nchar(df$animal_num)) # adds column reported retro-label
  df$retro_label      <- replace(df$retro_label, df$retro_label == "L", 1)                    # adds labeled for L and not labeled for N
  df$retro_label      <- replace(df$retro_label, df$retro_label == "N", 0)
  df$retro_label      <- as.numeric(df$retro_label)
  df$location         <- lapply(df$file, function(x) unlist(strsplit(x, "-"))[3])
  df$location         <- as.character(df$location)
  df$animal_num       <- lapply(df$animal_num, function(x) unlist(strsplit(x, "L"))[1])       # removes L or N from animal number
  df$animal_num       <- lapply(df$animal_num, function(x) unlist(strsplit(x, "N"))[1])       # removes L or N from animal number
  df$animal_num       <- as.numeric(df$animal_num)
  
  
  
  df$density_overall[df$density_overall > 4.00]   <- NA
  df$density_mushroom[df$density_mushroom > 1.00] <- NA
  df$density_thin[df$density_thin > 4.00]         <- NA
  df$density_stubby[df$density_stubby > 1.00]     <- NA
  
  

  
  
  #df[!complete.cases(df),]                         <- NA
  
    data_all            <- rbind(data_all, df)                                                  # adds df as next row in the data_master file
  
}                                                                                             # end of file for-loop 


rat_ID <- data_all$animal_num[1]
data_all$MAX_DTS<- as.numeric(data_all$MAX_DTS)
data_all$RAYBURST_VOLUME <- as.numeric(data_all$RAYBURST_VOLUME)
data_all$HEAD_DIAMETER <- as.numeric(data_all$HEAD_DIAMETER)
#A_A

data_mast_A_A <- data_all %>%
  group_by(file) %>%
  summarise(den_ov_A_A   = mean(density_overall, na.rm = TRUE),
            den_mush_A_A = mean(density_mushroom, na.rm = TRUE),
            den_thin_A_A = mean(density_thin, na.rm = TRUE),
            den_stub_A_A = mean(density_stubby, na.rm = TRUE),
            vol_ov_A_A   = mean(RAYBURST_VOLUME, na.rm = TRUE),
            len_ov_A_A   = mean(MAX_DTS, na.rm = TRUE),
            hd_ov_A_A    = mean(HEAD_DIAMETER, na.rm = TRUE)) %>%
  ungroup()

data_mast_A_A <- data_mast_A_A %>%
  summarise(den_ov   = mean(den_ov_A_A, na.rm = TRUE),
            den_mush = mean(den_mush_A_A, na.rm = TRUE),
            den_thin = mean(den_thin_A_A, na.rm = TRUE),
            den_stub = mean(den_stub_A_A, na.rm = TRUE),
            vol_ov   = mean(vol_ov_A_A, na.rm = TRUE),
            len_ov   = mean(len_ov_A_A, na.rm = TRUE),
            hd_ov    = mean(hd_ov_A_A, na.rm = TRUE)) %>%
  ungroup()

data_TYPE_A_A <- data_all %>%
  group_by(file, TYPE) %>%
  summarise(vol_TYPE_A_A = mean(RAYBURST_VOLUME, na.rm = TRUE),
            len_TYPE_A_A = mean(MAX_DTS, na.rm = TRUE),
            hd_TYPE_A_A  = mean(HEAD_DIAMETER, na.rm = TRUE)) %>%
  ungroup()

data_TYPE_A_A <- data_TYPE_A_A %>%
  group_by(TYPE)               %>%
  summarise(vol_TYPE = mean(vol_TYPE_A_A, na.rm = TRUE),
            len_TYPE = mean(len_TYPE_A_A, na.rm = TRUE),
            hd_TYPE  = mean(hd_TYPE_A_A, na.rm = TRUE))  %>%
  ungroup()

data_mush_A_A <- subset(data_TYPE_A_A, data_TYPE_A_A$TYPE == 'mushroom', select = 2:4)
data_thin_A_A <- subset(data_TYPE_A_A, data_TYPE_A_A$TYPE == 'thin',     select = 2:4)
data_stub_A_A <- subset(data_TYPE_A_A, data_TYPE_A_A$TYPE == 'stubby',   select = 2:4)

label_location <- c('A_A')

data_mast_A_A <- cbind.fill(label_location, 
                            data_mast_A_A, 
                            data_mush_A_A, 
                            data_thin_A_A, 
                            data_stub_A_A, 
                            fill = 0)


#N_A

data_mast_NL_A <- data_all %>%
  group_by(file, retro_label)%>%
  summarise(den_ov_NL_A   = mean(density_overall, na.rm = TRUE),
            den_mush_NL_A = mean(density_mushroom, na.rm = TRUE),
            den_thin_NL_A = mean(density_thin, na.rm = TRUE),
            den_stub_NL_A = mean(density_stubby, na.rm = TRUE),
            vol_ov_NL_A   = mean(RAYBURST_VOLUME, na.rm = TRUE),
            len_ov_NL_A   = mean(MAX_DTS, na.rm = TRUE),
            hd_ov_NL_A    = mean(HEAD_DIAMETER, na.rm = TRUE))   %>%
  ungroup()

data_mast_NL_A <- data_mast_NL_A %>%
  group_by(retro_label) %>%
  summarise(den_ov   = mean(den_ov_NL_A, na.rm = TRUE),
            den_mush = mean(den_mush_NL_A, na.rm = TRUE),
            den_thin = mean(den_thin_NL_A, na.rm = TRUE),
            den_stub = mean(den_stub_NL_A, na.rm = TRUE),
            vol_ov   = mean(vol_ov_NL_A, na.rm = TRUE),
            len_ov   = mean(len_ov_NL_A, na.rm = TRUE),
            hd_ov    = mean(hd_ov_NL_A, na.rm = TRUE))    %>%
  ungroup()

data_mast_N_A <- subset(data_mast_NL_A, data_mast_NL_A$retro_label == 0, select = 2:8)

data_TYPE_NL_A <- data_all %>%
  group_by(file, retro_label, TYPE) %>%
  summarise(vol_TYPE_NL_A = mean(RAYBURST_VOLUME, na.rm = TRUE),
            len_TYPE_NL_A = mean(MAX_DTS, na.rm = TRUE),
            hd_TYPE_NL_A  = mean(HEAD_DIAMETER, na.rm = TRUE))   %>%
  ungroup()

data_TYPE_NL_A <- data_TYPE_NL_A %>%
  group_by(retro_label, TYPE) %>%
  summarise(vol_TYPE = mean(vol_TYPE_NL_A, na.rm = TRUE),
            len_TYPE = mean(len_TYPE_NL_A, na.rm = TRUE),
            hd_TYPE = mean(hd_TYPE_NL_A, na.rm = TRUE))  %>%
  ungroup()


data_mush_N_A <- subset(data_TYPE_NL_A, data_TYPE_NL_A$retro_label == 0 & data_TYPE_NL_A$TYPE == 'mushroom', select= 3:5)
data_thin_N_A <- subset(data_TYPE_NL_A, data_TYPE_NL_A$retro_label == 0 & data_TYPE_NL_A$TYPE == 'thin',     select= 3:5)
data_stub_N_A <- subset(data_TYPE_NL_A, data_TYPE_NL_A$retro_label == 0 & data_TYPE_NL_A$TYPE == 'stubby',   select= 3:5)

label_location <- c('N_A')

data_mast_N_A <- cbind.fill(label_location, 
                            data_mast_N_A, 
                            data_mush_N_A, 
                            data_thin_N_A, 
                            data_stub_N_A, 
                            fill = 0)


#L_A

data_mast_L_A <- subset(data_mast_NL_A, data_mast_NL_A$retro_label == 1,  select = 2:8)
data_mush_L_A <- subset(data_TYPE_NL_A, data_TYPE_NL_A$retro_label == 1 & data_TYPE_NL_A$TYPE == 'mushroom', select = 3:5)
data_thin_L_A <- subset(data_TYPE_NL_A, data_TYPE_NL_A$retro_label == 1 & data_TYPE_NL_A$TYPE == 'thin',     select = 3:5)
data_stub_L_A <- subset(data_TYPE_NL_A, data_TYPE_NL_A$retro_label == 1 & data_TYPE_NL_A$TYPE == 'stubby',   select = 3:5)


label_location <- c('L_A')

data_mast_L_A <- cbind.fill(label_location, 
                            data_mast_L_A, 
                            data_mush_L_A, 
                            data_thin_L_A, 
                            data_stub_L_A, 
                            fill = 0)

# for A_b

data_mast_A_bpd <- data_all %>%
  group_by(file, location) %>%
  summarise(den_ov_A_b   = mean(density_overall, na.rm = TRUE),
            den_mush_A_b = mean(density_mushroom, na.rm = TRUE),
            den_thin_A_b = mean(density_thin, na.rm = TRUE),
            den_stub_A_b = mean(density_stubby, na.rm = TRUE),
            vol_ov_A_b   = mean(RAYBURST_VOLUME, na.rm = TRUE),
            len_ov_A_b   = mean(MAX_DTS, na.rm = TRUE),
            hd_ov_A_b    = mean(HEAD_DIAMETER, na.rm = TRUE)) %>%
  ungroup()

data_mast_A_bpd <- data_mast_A_bpd %>%
  group_by(location) %>%
  summarise(den_ov   = mean(den_ov_A_b, na.rm = TRUE),
            den_mush = mean(den_mush_A_b, na.rm = TRUE),
            den_thin = mean(den_thin_A_b, na.rm = TRUE),
            den_stub = mean(den_stub_A_b, na.rm = TRUE),
            vol_ov   = mean(vol_ov_A_b, na.rm = TRUE),
            len_ov   = mean(len_ov_A_b, na.rm = TRUE),
            hd_ov    = mean(hd_ov_A_b, na.rm = TRUE)) %>%
  ungroup()

data_mast_A_b <- subset(data_mast_A_bpd, data_mast_A_bpd$location == 'b', select = 2:8)


data_TYPE_A_bpd <- data_all %>%
  group_by(file, location, TYPE) %>%
  summarise(vol_TYPE_A_b = mean(RAYBURST_VOLUME, na.rm = TRUE),
            len_TYPE_A_b = mean(MAX_DTS, na.rm = TRUE),
            hd_TYPE_A_b  = mean(HEAD_DIAMETER, na.rm = TRUE)) %>%
  ungroup()

data_TYPE_A_bpd <- data_TYPE_A_bpd %>%
  group_by(location, TYPE) %>%
  summarise(vol_TYPE = mean(vol_TYPE_A_b, na.rm = TRUE),
            len_TYPE = mean(len_TYPE_A_b, na.rm = TRUE),
            hd_TYPE  = mean(hd_TYPE_A_b, na.rm = TRUE)) %>%
  ungroup()

data_mush_A_b <- subset(data_TYPE_A_bpd, data_TYPE_A_bpd$location == 'b' & data_TYPE_A_bpd$TYPE == 'mushroom', select = 3:5)
data_thin_A_b <- subset(data_TYPE_A_bpd, data_TYPE_A_bpd$location == 'b' & data_TYPE_A_bpd$TYPE == 'thin',     select = 3:5)
data_stub_A_b <- subset(data_TYPE_A_bpd, data_TYPE_A_bpd$location == 'b' & data_TYPE_A_bpd$TYPE == 'stubby',   select = 3:5)

label_location <- c('A_b')

data_mast_A_b <- cbind.fill(label_location, 
                            data_mast_A_b, 
                            data_mush_A_b, 
                            data_thin_A_b, 
                            data_stub_A_b,
                            fill = 0)

#A_p

data_mast_A_p <- subset(data_mast_A_bpd, data_mast_A_bpd$location == 'p', select = 2:8)
data_mush_A_p <- subset(data_TYPE_A_bpd, data_TYPE_A_bpd$location == 'p' & data_TYPE_A_bpd$TYPE == 'mushroom', select = 3:5)
data_thin_A_p <- subset(data_TYPE_A_bpd, data_TYPE_A_bpd$location == 'p' & data_TYPE_A_bpd$TYPE == 'thin',     select = 3:5)
data_stub_A_p <- subset(data_TYPE_A_bpd, data_TYPE_A_bpd$location == 'p' & data_TYPE_A_bpd$TYPE == 'stubby',   select = 3:5)

label_location <- c('A_p')

data_mast_A_p <- cbind.fill(label_location, 
                            data_mast_A_p, 
                            data_mush_A_p, 
                            data_thin_A_p, 
                            data_stub_A_p, 
                            fill = 0)

#A_d

data_mast_A_d <- subset(data_mast_A_bpd, data_mast_A_bpd$location == 'd',  select = 2:8)
data_mush_A_d <- subset(data_TYPE_A_bpd, data_TYPE_A_bpd$location == 'd' & data_TYPE_A_bpd$TYPE == 'mushroom', select = 3:5)
data_thin_A_d <- subset(data_TYPE_A_bpd, data_TYPE_A_bpd$location == 'd' & data_TYPE_A_bpd$TYPE == 'thin',     select = 3:5)
data_stub_A_d <- subset(data_TYPE_A_bpd, data_TYPE_A_bpd$location == 'd' & data_TYPE_A_bpd$TYPE == 'stubby',   select = 3:5)

label_location <- c('A_d')

data_mast_A_d <- cbind.fill(label_location, 
                            data_mast_A_d, 
                            data_mush_A_d, 
                            data_thin_A_d, 
                            data_stub_A_d, 
                            fill = 0)

#for N_b

data_mast_NL_bpd <- data_all %>%
  group_by(file, retro_label, location) %>%
  summarise(den_ov_NL_bpd   = mean(density_overall, na.rm = TRUE),
            den_mush_NL_bpd = mean(density_mushroom, na.rm = TRUE),
            den_thin_NL_bpd = mean(density_thin, na.rm = TRUE),
            den_stub_NL_bpd = mean(density_stubby, na.rm = TRUE),
            vol_ov_NL_bpd   = mean(RAYBURST_VOLUME, na.rm = TRUE),
            len_ov_NL_bpd   = mean(MAX_DTS, na.rm = TRUE),
            hd_ov_NL_bpd    = mean(HEAD_DIAMETER, na.rm = TRUE)) %>%
  ungroup()

data_count_master <- data_mast_NL_bpd[,1:3]

count_N_b <- data.frame(nrow(subset(data_count_master, retro_label == 0 & location =='b')))
count_N_p <- data.frame(nrow(subset(data_count_master, retro_label == 0 & location =='p')))
count_N_d <- data.frame(nrow(subset(data_count_master, retro_label == 0 & location =='d')))
count_N_A <- sum(count_N_b, count_N_p, count_N_d)

count_L_b <- data.frame(nrow(subset(data_count_master, retro_label == 1 & location =='b')))
count_L_p <- data.frame(nrow(subset(data_count_master, retro_label == 1 & location =='p')))
count_L_d <- data.frame(nrow(subset(data_count_master, retro_label == 1 & location =='d')))
count_L_A <- sum(count_L_b, count_L_p, count_L_d)

count_A_A <- sum(count_N_A, count_L_A)
count_A_b <- sum(count_N_b, count_L_b)
count_A_p <- sum(count_N_p, count_L_p)
count_A_d <- sum(count_N_d, count_L_d)



data_count_master <- cbind.fill(count_A_A,
                                count_A_b,
                                count_A_p,
                                count_A_d,
                                count_N_A,
                                count_N_b, 
                                count_N_p, 
                                count_N_d,
                                count_L_A,
                                count_L_b, 
                                count_L_p, 
                                count_L_d, 
                                fill = 0)

data_count_master            <- t(data_count_master)
row.names(data_count_master) <- NULL


data_mast_NL_bpd <- data_mast_NL_bpd %>%
  group_by(retro_label, location) %>%
  summarise(den_ov   = mean(den_ov_NL_bpd, na.rm = TRUE),
            den_mush = mean(den_mush_NL_bpd, na.rm = TRUE),
            den_thin = mean(den_thin_NL_bpd, na.rm = TRUE),
            den_stub = mean(den_stub_NL_bpd, na.rm = TRUE),
            vol_ov   = mean(vol_ov_NL_bpd, na.rm = TRUE),
            len_ov   = mean(len_ov_NL_bpd, na.rm = TRUE),
            hd_ov   = mean(hd_ov_NL_bpd, na.rm = TRUE)) %>%
  ungroup()

data_mast_N_b <- subset(data_mast_NL_bpd, data_mast_NL_bpd$retro_label == 0 & data_mast_NL_bpd$location == 'b', select = 3:9)

data_TYPE_NL_bpd <- data_all %>%
  group_by(file, retro_label, location, TYPE) %>%
  summarise(vol_TYPE_NL_bpd = mean(RAYBURST_VOLUME, na.rm = TRUE),
            len_TYPE_NL_bpd = mean(MAX_DTS, na.rm = TRUE),
            hd_TYPE_NL_bpd  = mean(HEAD_DIAMETER, na.rm = TRUE)) %>%
  ungroup()

data_TYPE_NL_bpd <- data_TYPE_NL_bpd %>%
  group_by(retro_label, location, TYPE) %>%
  summarise(vol_TYPE = mean(vol_TYPE_NL_bpd, na.rm = TRUE),
            len_TYPE = mean(len_TYPE_NL_bpd, na.rm = TRUE),
            hd_TYPE  = mean(hd_TYPE_NL_bpd, na.rm = TRUE)) %>%
  ungroup()

data_mush_N_b <- subset(data_TYPE_NL_bpd, 
                        data_TYPE_NL_bpd$retro_label == 0 & 
                          data_TYPE_NL_bpd$location    == 'b' & 
                          data_TYPE_NL_bpd$TYPE        == 'mushroom', 
                        select = 4:6)

data_thin_N_b <- subset(data_TYPE_NL_bpd, 
                        data_TYPE_NL_bpd$retro_label == 0 & 
                          data_TYPE_NL_bpd$location    == 'b' & 
                          data_TYPE_NL_bpd$TYPE        == 'thin', 
                        select = 4:6)

data_stub_N_b <- subset(data_TYPE_NL_bpd, 
                        data_TYPE_NL_bpd$retro_label == 0 & 
                          data_TYPE_NL_bpd$location    == 'b' & 
                          data_TYPE_NL_bpd$TYPE        == 'stubby', 
                        select = 4:6)


label_location <- c('N_b')

data_mast_N_b <- cbind.fill(label_location, 
                            data_mast_N_b, 
                            data_mush_N_b,
                            data_thin_N_b,
                            data_stub_N_b, 
                            fill = 0)

#N_p

data_mast_N_p <- subset(data_mast_NL_bpd, data_mast_NL_bpd$retro_label == 0 & data_mast_NL_bpd$location == 'p', select = 3:9)


data_mush_N_p <- subset(data_TYPE_NL_bpd, 
                        data_TYPE_NL_bpd$retro_label == 0 & 
                          data_TYPE_NL_bpd$location    == 'p' & 
                          data_TYPE_NL_bpd$TYPE        == 'mushroom', 
                        select = 4:6)


data_thin_N_p <- subset(data_TYPE_NL_bpd, 
                        data_TYPE_NL_bpd$retro_label == 0 & 
                          data_TYPE_NL_bpd$location    == 'p' & 
                          data_TYPE_NL_bpd$TYPE        == 'thin', 
                        select = 4:6)


data_stub_N_p <- subset(data_TYPE_NL_bpd, 
                        data_TYPE_NL_bpd$retro_label == 0 & 
                          data_TYPE_NL_bpd$location    == 'p' & 
                          data_TYPE_NL_bpd$TYPE        == 'stubby', 
                        select = 4:6)


label_location <- c('N_p')

data_mast_N_p <- cbind.fill(label_location, 
                            data_mast_N_p, 
                            data_mush_N_p, 
                            data_thin_N_p, 
                            data_stub_N_p, 
                            fill = 0)

#N_d

data_mast_N_d <- subset(data_mast_NL_bpd, 
                        data_mast_NL_bpd$retro_label == 0 & 
                          data_mast_NL_bpd$location    == 'd', 
                        select = 3:9)

data_mush_N_d <- subset(data_TYPE_NL_bpd, 
                        data_TYPE_NL_bpd$retro_label == 0 & 
                          data_TYPE_NL_bpd$location    == 'd' & 
                          data_TYPE_NL_bpd$TYPE        == 'mushroom', 
                        select = 4:6)

data_thin_N_d <- subset(data_TYPE_NL_bpd, 
                        data_TYPE_NL_bpd$retro_label == 0 & 
                          data_TYPE_NL_bpd$location    == 'd' & 
                          data_TYPE_NL_bpd$TYPE        == 'thin', 
                        select = 4:6)

data_stub_N_d <- subset(data_TYPE_NL_bpd, 
                        data_TYPE_NL_bpd$retro_label == 0 & 
                          data_TYPE_NL_bpd$location    == 'd' & 
                          data_TYPE_NL_bpd$TYPE        == 'stubby', 
                        select = 4:6)

label_location <- c('N_d')

data_mast_N_d <- cbind.fill(label_location, 
                            data_mast_N_d, 
                            data_mush_N_d, 
                            data_thin_N_d, 
                            data_stub_N_d, 
                            fill = 0)

#L_b

data_mast_L_b <- subset(data_mast_NL_bpd, 
                        data_mast_NL_bpd$retro_label == 1 & 
                          data_mast_NL_bpd$location    == 'b', 
                        select = 3:9)


data_mush_L_b <- subset(data_TYPE_NL_bpd, 
                        data_TYPE_NL_bpd$retro_label == 1 & 
                          data_TYPE_NL_bpd$location    == 'b' & 
                          data_TYPE_NL_bpd$TYPE        == 'mushroom', 
                        select = 4:6)


data_thin_L_b <- subset(data_TYPE_NL_bpd, 
                        data_TYPE_NL_bpd$retro_label == 1 & 
                          data_TYPE_NL_bpd$location    == 'b' & 
                          data_TYPE_NL_bpd$TYPE        == 'thin', 
                        select = 4:6)


data_stub_L_b <- subset(data_TYPE_NL_bpd, 
                        data_TYPE_NL_bpd$retro_label == 1 & 
                          data_TYPE_NL_bpd$location    == 'b' & 
                          data_TYPE_NL_bpd$TYPE        == 'stubby', 
                        select = 4:6)


label_location <- c('L_b')
data_mast_L_b <- cbind.fill(label_location, 
                            data_mast_L_b, 
                            data_mush_L_b, 
                            data_thin_L_b, 
                            data_stub_L_b, 
                            fill = 0)

#L_p

data_mast_L_p <- subset(data_mast_NL_bpd, 
                        data_mast_NL_bpd$retro_label == 1 & 
                          data_mast_NL_bpd$location    == 'p', 
                        select = 3:9)


data_mush_L_p <- subset(data_TYPE_NL_bpd, 
                        data_TYPE_NL_bpd$retro_label == 1 & 
                          data_TYPE_NL_bpd$location    == 'p' & 
                          data_TYPE_NL_bpd$TYPE        == 'mushroom', 
                        select = 4:6)


data_thin_L_p <- subset(data_TYPE_NL_bpd, 
                        data_TYPE_NL_bpd$retro_label == 1 & 
                          data_TYPE_NL_bpd$location    == 'p' & 
                          data_TYPE_NL_bpd$TYPE        == 'thin', 
                        select = 4:6)


data_stub_L_p <- subset(data_TYPE_NL_bpd, 
                        data_TYPE_NL_bpd$retro_label == 1 & 
                          data_TYPE_NL_bpd$location    == 'p' & 
                          data_TYPE_NL_bpd$TYPE        == 'stubby', 
                        select = 4:6)



label_location <- c('L_p')

data_mast_L_p <- cbind.fill(label_location, 
                            data_mast_L_p, 
                            data_mush_L_p, 
                            data_thin_L_p, 
                            data_stub_L_p, 
                            fill = 0)

#L_d

data_mast_L_d <- subset(data_mast_NL_bpd, 
                        data_mast_NL_bpd$retro_label == 1 & 
                          data_mast_NL_bpd$location    == 'd', 
                        select = 3:9)


data_mush_L_d <- subset(data_TYPE_NL_bpd, 
                        data_TYPE_NL_bpd$retro_label == 1 & 
                          data_TYPE_NL_bpd$location    == 'd' & 
                          data_TYPE_NL_bpd$TYPE        == 'mushroom', 
                        select = 4:6)


data_thin_L_d <- subset(data_TYPE_NL_bpd, 
                        data_TYPE_NL_bpd$retro_label == 1 & 
                          data_TYPE_NL_bpd$location    == 'd' & 
                          data_TYPE_NL_bpd$TYPE        == 'thin', 
                        select = 4:6)


data_stub_L_d <- subset(data_TYPE_NL_bpd, 
                        data_TYPE_NL_bpd$retro_label == 1 & 
                          data_TYPE_NL_bpd$location    == 'd' & 
                          data_TYPE_NL_bpd$TYPE        == 'stubby', 
                        select = 4:6)


label_location <- c('L_d')

data_mast_L_d <- cbind.fill(label_location, 
                            data_mast_L_d, 
                            data_mush_L_d, 
                            data_thin_L_d, 
                            data_stub_L_d, 
                            fill = 0)




data_master <- rbind(data_mast_A_A, 
                     data_mast_A_b, 
                     data_mast_A_p, 
                     data_mast_A_d, 
                     data_mast_N_A, 
                     data_mast_N_b, 
                     data_mast_N_p, 
                     data_mast_N_d, 
                     data_mast_L_A, 
                     data_mast_L_b, 
                     data_mast_L_p, 
                     data_mast_L_d)

data_master <- cbind.fill (rat_ID, 
                           data_count_master, 
                           data_master)

colnames(data_master) <- c('rat_ID', 
                           'count',
                           'group',
                           'den_ov', 
                           'den_mush', 
                           'den_thin', 
                           'den_stub', 
                           'vol_ov', 
                           'len_ov', 
                           'hd_ov', 
                           'vol_mush', 
                           'len_mush', 
                           'hd_mush', 
                           'vol_thin', 
                           'len_thin', 
                           'hd_thin',
                           'vol_stub', 
                           'len_stub', 
                           'hd_stub')



savePath_all <- paste0(saveFolder, "/data_all") # creates the path where files can be saved
setwd(savePath_all)
write.csv(data_all, paste0(rat_ID, "_data_all.csv"))
#setwd(fileLoc) # sets working directory to the path created above
#dir.create(savePath, "data_all") #create a directory to create a file inside
setwd(saveFolder)
savePath_master <- paste0(saveFolder, "/data_master")
setwd(savePath_master)
write.csv(data_master, paste0(rat_ID, "_data_master.csv"))

}


winDialog("ok", "Code complete")