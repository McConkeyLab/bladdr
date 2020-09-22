
# Goal --------------------------------------------------------------------
    ## Create Datasets that will be used in other scripts; create ones that are
    ## easy to read in later and clean

    ## for now - no data about amplifications. either Y, N, M

# Set Workplace -----------------------------------------------------------

library(tidyverse)  ##includes ggplot2, dplyr
library(readxl)
library(lubridate)


# Find Duplicates Function ------------------------------------------------
## returns all duplicates, not just second
## needs to be in columnn called sample to work

## re-run this function until you only get "FALSE"
remove_first_duplicates <- function(data) {
    find <- data[duplicated(data$sample),]
    others <- data$sample %in% find$sample
    duplicates <- data[others,]
    other_data <- data[!others,]
    keepers <- duplicates[duplicated(duplicates$sample),]
    end_data <- rbind(other_data, keepers)
    print(table(duplicated(end_data$sample)))
    end_data
}


# Allegheny Aim 1 ---------------------------------------------------------
a1 <- read_excel("./b_datasets/a1-lib_prep.xlsx",
                       sheet = 1, skip = 3, col_names = T,
                       na = "")
a1 <- a1[,c(3, 4, 5, 6, 7, 8, 9, 15, 18)] %>%
    `colnames<-`(c("sample", "nano_conc", "protein", "salt", "tape_conc",
                  "conc_200", "dv200", "input_rna", "success"))

a1 <- a1 %>% mutate(., tech = "bridget", project = "a1",
             success = sub("MAYBE", "M", a1$success))

a1 <- a1 %>% mutate(., success = sub("weird", NA, a1$success))

## need to remove duplicated samples - keeping the last one that was done
## duplicated only marks 2nd duplicate as a duplicate, this dataset is in the
## order that the library preps were done


a1 <- remove_first_duplicates(a1) %>% remove_first_duplicates(.)

## all duplicates removed - only keeping LAST library made.



# Allegheny Aim 2 ---------------------------------------------------------
a2_total <- read_excel("./b_datasets/a2-lib_prep.xlsx",
                       sheet = 1,
                       na = "") %>%
    mutate(ffpe_date = ymd(ffpe_date)) %>%
    mutate(rin = ifelse(tape_conc < 25, rin == "NA", rin),
           project = "allegheny aim 2")


a2_sub <- a2_total$success %in% c("Y", "N")
a2 <- a2_total[a2_sub,]

remove(a2_sub)


# ADAPT-BLADDER -----------------------------------------------------------
ab_total <- read_excel("./b_datasets/adapt_bladder-library_prep.xlsx",
                       sheet = 1,
                       col_types = c("text", "date", "text", rep("numeric", 8),
                                     "text"),
                       na = "") %>%
    mutate(ffpe_date = ymd(ffpe_date)) %>%
    mutate(rin = ifelse(tape_conc < 25, rin == "NA", rin),
           project = "ADAPT-BLADDER")

ab_sub <- ab_total$success %in% c("Y", "N")
ab <- ab_total[ab_sub,]

remove(ab_sub)


# GBCI --------------------------------------------------------------------
gbci_total <- read_excel("./b_datasets/gbci-library_prep.xlsx", sheet = 1,
                         col_types = c("text", "date", "text",
                                       rep("numeric", 8), "text"),
                         na = "") %>%
    mutate(ffpe_date = ymd(ffpe_date)) %>%
    mutate(rin = ifelse(tape_conc < 25, rin == "NA", rin),
           project = "MIBC gbci")


gbci_sub <- gbci_total$success %in% c("Y", "N")
gbci <- gbci_total[gbci_sub,]

remove(gbci_sub)


# All Together Now --------------------------------------------------------

libraries_sub <- rbind(a1, a2, ab, gbci) %>%
    mutate(., kit = as.factor(kit),
           success = as.factor(success),
           project = as.factor(project),
           tech = "bridget")

libraries_total <- rbind(a1_total, a2_total, ab_total, gbci_total) %>%
    mutate(., kit = as.factor(kit),
           success = as.factor(success),
           project = as.factor(project),
           tech = "bridget")

libraries_binary <- libraries_sub %>%
    mutate(., success = case_when(success == "Y" ~ 1, success == "N" ~ 0))


# Remove Smaller Datasets (optional) --------------------------------------
remove(a1, a1_total)
remove(a2, a2_total)
remove(ab, ab_total)
remove(gbci, gbci_total)

# Create Train, Test, Validate --------------------------------------------
## want equal amounts of failed + succeeded for train dataset.
## NOT FINAL - JUST MESSING AROUND HERE.


train_n <- libraries_sub[libraries_sub$success == "N",]
set.seed(100)
train_y <- libraries_sub[!libraries_sub$success == "N",] %>%
    slice_sample(., n = nrow(train_n))
train <- rbind(train_n, train_y) %>% slice_sample(., prop = 1.0) %>%
    mutate(., success_binary = case_when(success == "Y" ~ 1,
                                         success == "N" ~ 0))
remove(train_n, train_y)


# Write Datasets for Later Use --------------------------------------------
write.table(libraries_total, "libraries_total")
write.table(libraries_binary, "libraries_binary")
write.table(train, "training_data")

