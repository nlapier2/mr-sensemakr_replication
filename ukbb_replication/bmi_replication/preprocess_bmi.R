# Setup -------------------------------------------------------------------

# Cleans workspace and loads package
rm(list = ls())

## You should create a folder called "packages" first
lib <-  "../packages/"

## now change the libPaths to include the default + the folder
current.lib <- .libPaths()
.libPaths(c(lib, current.lib))

## checks if it works
.libPaths()

# File location
file <- "/u/home/b/blhill/code/mr_ukb/UKBB_features.tsv"

# Loads data table for fast loading
library(data.table)

# dplyr for some manipulations
library(dplyr)

data <- fread(file, stringsAsFactors = FALSE)


# creates t2d data
data[,t2d:= ifelse(diabetes == "Yes", 1, ifelse(diabetes == "No", 0, NA))]
data[, summary(t2d)]


data[, summary(age)]
data[, .(mean_age = mean(age, na.rm = TRUE), 
         sd_age = sd(age, na.rm = TRUE))]


data[, summary(height)]
data[, .(mean_age = mean(height, na.rm = TRUE), 
         sd_age = sd(height, na.rm = TRUE))]
data  %>% group_by(sex)  %>% summarise(mean = mean(height, na.rm = T),
                                      sd = sd(height, na.rm = T))


data[, summary(weight)]
data[, .(mean_age = mean(weight, na.rm = TRUE), 
         sd_age = sd(weight, na.rm = TRUE))]
data  %>% group_by(sex)  %>% summarise(mean = mean(weight, na.rm = T),
                                      sd = sd(weight, na.rm = T))

data[, summary(bmi)]
data[, .(mean_bmi = mean(bmi, na.rm = T), sd_bmi = sd(bmi, na.rm = T))]

data[, table(frequency_alcohol)/nrow(data)]

data[, table(alcohol_frequency)/nrow(data)]

data[, table(smoking_status)]
data[, table(smoking_status)/nrow(data)]

data[,sbp:= systolic_blood_pressure][,summary(sbp)]

data[,.(mean_sbp = mean(sbp, na.rm = TRUE), 
        sd_sbp = sd(sbp, na.rm = TRUE))]

data[,dbp:= diastolic_blood_pressure][,summary(dbp)]

data[,.(mean_dbp = mean(dbp, na.rm = T),
       sd_dbp = sd(dbp, na.rm = T))]

data[,summary(pulse_rate)]

data[ ,.(mean_pr = mean(pulse_rate, na.rm = T),
      sd_pr = sd(pulse_rate, na.rm = T))]

data[,chd:=as.numeric(derived_chd)][,summary(chd)]

data[,angina:=as.numeric(derived_Angina)][,summary(angina)]

data[,stroke:=as.numeric(derived_Stroke)][,summary(stroke)]

data[,heart_attack:=as.numeric(`derived_Heart attack`)][,summary(heart_attack)]

data[,hbp:=as.numeric(`derived_High blood pressure`)][,summary(hbp)]

data[,table(average_income)/nrow(data)]

data %>% group_by(sex, average_income) %>% summarise(n = n())  %>% mutate(freq = n/sum(n))

data[,summary(townsend )]
#hist(data$townsend, breaks= 100, col = "darkblue")

data[,.(mean_tws = mean(townsend, na.rm = T),
        sd_tws = sd(townsend, na.rm = T))]

data  %>% group_by(sex) %>% summarise(mean = mean(townsend, na.rm = TRUE),
                                     sd = sd(townsend, na.rm = TRUE))
data[,prs:= PRS_split1][is.na(prs), prs:= PRS_split2 ]

# creates alcohol
data[,alcohol:= alcohol_frequency]

# creates smoking
data[,smoking:= smoking_status ]

data[,meds:= as.numeric(derived_medication_cholesterol_blood_pressure_diabetes)]

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
            income = average_income,
            prs,
            bmi,
            t2d,
            townsend,
            alcohol, 
            smoking,
            weight,
            body_fat,
            sbp,
            dbp,
            hbp,
            meds,
            pulse_rate,
            chd,
            angina,
            stroke,
            heart_attack,
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

saveRDS(object = out_data, file = "../processed_data/out_data.rds")

# filter 4
filter_file <- "/u/project/sriram/ukbiobank/data/geno/cal/filter4.fam"
filter <- read.table(filter_file)

filter4_data <- out_data[out_data$EID %in% filter$V1, ] 

saveRDS(object = filter4_data, file = "../processed_data/filter4_data.rds")
#
