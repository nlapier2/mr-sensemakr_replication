# Setup -------------------------------------------------------------------

# Cleans workspace and loads package
#rm(list = ls())

## You should create a folder called "packages" first
#lib <-  "../packages/"

## now change the libPaths to include the default + the folder
#current.lib <- .libPaths()
#.libPaths(c(lib, current.lib))

## checks if it works
#.libPaths()

# File location
#file <- "/u/home/b/blhill/code/mr_ukb/UKBB_features.tsv"
#file <- "/u/project/sriram/ukbiobank/data/mr_ukb_split/data/UKBB_features_expanded.tsv"
file <- "ukbb_data_with_prs.tsv"

# Loads data table for fast loading
library(data.table)

# dplyr for some manipulations
library(dplyr)


data <- fread(file, stringsAsFactors = FALSE)

# creates hdl, ldl, triglycerides, chd
data[,hdl:= HDL_cholesterol]
data[,ldl:= LDL_direct]
data[,triglycerides:= Triglycerides]
data[,chd:=as.numeric(derived_chd)][,summary(chd)]
data[,prs_hdl:= PRS_hdl]
data[,prs_ldl:= PRS_ldl]

# creates t2d data
#data[,t2d:= ifelse(diabetes == "Yes", 1, ifelse(diabetes == "No", 0, NA))]
#data[, summary(t2d)]


# creates chd (coronary heart disease) data
# data[,chd:= ifelse(derived_chd == "True", 1, ifelse(derived_chd == "False", 0, NA))]
# data[, summary(chd)]


data[, summary(age)]
data[, .(mean_age = mean(age, na.rm = TRUE), 
         sd_age = sd(age, na.rm = TRUE))]


data[, summary(bmi)]
data[, .(mean_bmi = mean(bmi, na.rm = T), sd_bmi = sd(bmi, na.rm = T))]

#data[, table(frequency_alcohol)/nrow(data)]

data[, table(alcohol_frequency)/nrow(data)]

data[, table(smoking_status)]
data[, table(smoking_status)/nrow(data)]

# creates alcohol
data[,alcohol:= alcohol_frequency]

# creates smoking
data[,smoking:= smoking_status ]

# reorder
out_data <- data[, 
       list(EID, 
            assessment_centre,
            genotype_batch, 
            genotype_plate, 
            genotype_well,
            ethnicity,
            sex,
            age,
            prs_hdl,
            prs_ldl,
            hdl,
            chd,
            ldl,
            triglycerides,
            bmi,
#            t2d,
            alcohol, 
            smoking,
            PC1,
            PC2,
            PC3,
            PC4,
            PC5,
            PC6,
            PC7,
            PC8,
            PC9,
            PC10,
            PC11,
            PC12,
            PC13,
            PC14,
            PC15,
            PC16,
            PC17,
            PC18,
            PC19,
            PC20
           )]

saveRDS(object = out_data, file = "processed_data/unfiltered_data.rds")

# filter 4
filter_file <- "/u/project/sriram/ukbiobank/data/geno/cal/filter4.fam"
filter <- read.table(filter_file)

filter4_data <- out_data[out_data$EID %in% filter$V1, ] 

saveRDS(object = filter4_data, file = "processed_data/filter4_data.rds")
#
