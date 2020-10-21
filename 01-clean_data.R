
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
                ## removes anything that was duplicated
        keepers <- duplicates[duplicated(duplicates$sample),]
        end_data <- rbind(other_data, keepers)
        print(table(duplicated(end_data$sample)))
        end_data
}

remove_duplicates <- function(data) {
        data[!duplicated(data$sample),]
}


# Allegheny Aim 1 ---------------------------------------------------------
a1 <- read_excel("./b_datasets/a1-lib_prep.xlsx",
                 sheet = 1, skip = 3, col_names = T,
                 na = "")

## want to remove any samples that were speed-vacced
find <- grepl("speed", a1$DF)
a1 <- a1[!find,]

## need to find any samples that were not diluted for lib prep
## ("5.0 uL straight from tube")

find <- grepl(".*tube.*", a1$DF)

a1 <- a1[,c(3, 4, 5, 6, 7, 8, 9, 15, 18)] %>%
        `colnames<-`(c("sample", "nano_conc", "protein", "salt", "tape_conc",
                       "conc_200", "dv200", "input_rna", "success"))

a1 <- mutate(a1, input_rna = case_when(find == T ~ "undiluted",
                                       TRUE ~ as.character(input_rna)),
             tech = "bridget",
             project = "a1",
             kit = "roche",
             dv200 = dv200/100,
             success = sub("MAYBE", "M", success),
             success = sub("weird", NA, success))


a1 <- remove_duplicates(a1)



# Allegheny Aim 2 ---------------------------------------------------------
a2 <- read_excel("./b_datasets/a2-lib_prep.xlsx",
                 sheet = 1, skip = 5, col_names = T, na = "")

find <- grepl(".*tube.*", a2$DF)

a2 <- a2[,c(3, 4, 5, 6, 7, 8, 9, 17, 19)]  %>%
        `colnames<-`(c("sample", "nano_conc", "protein", "salt", "tape_conc",
                       "conc_200", "dv200", "input_rna", "success"))

a2 <- a2 %>% mutate(., tech = "bridget",
               project = "a2",
               kit = "qiagen",
               input_rna = case_when(find == T ~ "undiluted",
                                       TRUE ~ as.character(input_rna)))

a2_1 <- remove_duplicates(a2)
## 266 unique observations


# ADAPT-BLADDER -----------------------------------------------------------
ab <- read_excel("./b_datasets/ab-lib_prep.xlsx",
                 sheet = 1, skip = 4, col_names = T, na = "")

find <- grepl(".*tube.*", ab$DF)


ab <- ab[,c(2, 3, 4, 5, 6, 7, 8, 9, 17, 19)]  %>%
        `colnames<-`(c("sample", "timepoint", "nano_conc", "protein", "salt",
                       "tape_conc",
                       "conc_200", "dv200", "input_rna", "success"))

ab <- ab %>% mutate(., tech = "bridget",
               project = "ab",
               kit = "roche",
               sample = paste(sample, timepoint, sep = "_"),
               dv200 = dv200/100,
               input_rna = case_when(find == T ~ "undiluted",
                                       TRUE ~ as.character(input_rna)))


ab <- select(ab, !timepoint)

## no duplicates to remove



# GBCI --------------------------------------------------------------------
gbci <- read_excel("./b_datasets/gbci-lib_prep.xlsx", sheet = 1,
                   skip = 3, col_names = T, na = "")

find <- grepl(".*tube.*", gbci$DF)


gbci <- gbci[,c(3, 4, 5, 6, 7, 8, 9, 17)] %>%
        `colnames<-`(c("sample", "nano_conc", "protein", "salt",
                       "tape_conc",
                       "conc_200", "dv200", "success"))

gbci <- gbci %>% mutate(kit = "roche",
               tech = "bridget",
               input_rna = 20,
               project = "gbci",
               input_rna = case_when(find == T ~ "undiluted",
                                       TRUE ~ as.character(input_rna)))


gbci <- remove_duplicates(gbci)
## 23 unique observations


# CIS ---------------------------------------------------------------------
cis <- read_excel("./b_datasets/cis-lib_prep.xlsx", sheet = 1, skip = 4,
                  col_names = T, na = "")

find <- grepl(".*tube.*", cis$DF)

cis <- cis[,c(2, 3, 4, 5, 6, 7, 8, 16, 18)] %>%
        `colnames<-`(c("sample", "nano_conc", "protein", "salt",
                       "tape_conc",
                       "conc_200", "dv200", "input_rna", "success"))

cis <- cis %>% mutate(kit = "roche",
               tech = "bridget",
               project = "cis",
               dv200 = dv200/100,
               input_rna = case_when(find == T ~ "undiluted",
                                           TRUE ~ as.character(input_rna)))


## no duplicates to remove


# CS ----------------------------------------------------------------------
cs <- read_excel("./b_datasets/cs-lib_prep.xlsx", sheet = 1, skip = 3,
                 col_names = T, na = "")

cs <- cs[,c(2, 3, 4, 5, 6, 7, 8, 14, 17)] %>%
        `colnames<-`(c("sample", "nano_conc", "protein", "salt",
                       "tape_conc",
                       "conc_200", "dv200", "input_rna", "success"))

cs <- cs %>% mutate(.,kit = "roche",
               tech = "bridget",
               project = "cs",
               success = sub("MAYBE/Y|MAYBE", "M", success),
               dv200 = dv200/100)

cs <- remove_duplicates(cs)
## 49 unique observations - cores that were isolated alone (A,B, etc)
## and then the patients that were combined should count as separate
## samples (since they kind of are)

## no samples speed-vacced or taken directly from tube

# WC samples --------------------------------------------------------------
wc <- read_excel("./b_datasets/wc-lib_prep.xlsx", sheet = 1, skip = 2,
                 col_names = T, na = "")

wc <- wc[,c(3, 4, 5, 6, 7, 8, 9, 17)] %>%
        `colnames<-`(c("sample", "nano_conc", "protein", "salt",
                       "tape_conc",
                       "conc_200", "dv200", "success")) %>%
        mutate(.,kit = "roche",
               tech = "bridget",
               project = "wc",
               input_rna = 20)

## no duplicates, no samples spun down or taken directly from tube

# MSK Roche ---------------------------------------------------------------
msk_roche <- read_excel("./b_datasets/megan_data/MSK Ampliseq Roche Lib Master.xlsx",
                        sheet = 1, col_names = T, skip = 1)

msk_roche <- msk_roche[, c(2, 3, 4, 5, 6, 7, 8, 11)] %>%
        'colnames<-'(c("sample", "nano_conc", "protein", "salt", "tape_conc",
                       "conc_200", "dv200", "success")) %>%
        mutate(., kit = "roche",
               tech = "megan",
               project = "msk",
               input_rna = c(rep(20, 20), rep(5, 8), rep(10, 8), rep(20,3)))

msk_roche <- remove_duplicates(msk_roche)




# MSK Qiagen --------------------------------------------------------------
msk_qiagen <- read_excel("./b_datasets/megan_data/MSK Ampliseq Qiagen Lib Master_updated.xlsx",
                         sheet = 2, col_names = T, skip = 2, na = "")

findh20 <- msk_qiagen$H2O == 0
## finding all instances where the sample was undiluted, needs to be done
## before we delete this column

msk_qiagen <- msk_qiagen[, c(3, 4, 5, 6, 7, 8, 10, 11, 12, 26, 29)] %>%
        'colnames<-'(c("sample", "id", "type", "nano_conc", "protein", "salt",
                       "tape_conc", "conc_200", "dv200", "input_rna",
                       "success"))

msk_qiagen <- msk_qiagen %>% mutate(., kit = "qiagen",
               tech = "megan",
               project = "msk",
               success = sub("^[Mm].*","M", success),
               success = sub(".*<200bp.*|undetected", NA, success),
               success = sub("^Y$|.*barcodes.*", "Y", success))

find <- grepl(".*max.*", msk_qiagen$type)
msk_qiagen <- mutate(msk_qiagen, input_rna = case_when(find == T ~ "undiluted",
                                         findh20 == T ~ "undiluted",
                                         TRUE ~ as.character(input_rna)))

find <- grepl(".*new.*", msk_qiagen$type)
msk_qiagen <- mutate(msk_qiagen,
                     sample = case_when(find == T ~ paste0(sample, "b"),
                                        TRUE ~ sample))
## all 'repeat' samples, but with new RNA

find <- grepl(".*barcode.*", msk_qiagen$type)
msk_qiagen <- msk_qiagen[!find,]
## removing all repeats that were repeated because of the wrong barcode

find <- grepl("[Ss]peed.*", msk_qiagen$type)
msk_qiagen <- msk_qiagen[!find,]
## removing all samples that were speed-vacced


msk_qiagen <- remove_duplicates(msk_qiagen)

msk_qiagen <- select(msk_qiagen, !c(id, type))

# Chemo -------------------------------------------------------------------
chemo <- read_excel("./b_datasets/megan_data/Chemo_QC.xlsx", sheet = 1,
                    skip = 2, col_names = T, na = "")

chemo <- chemo[,c(2, 5, 6, 7, 14, 15, 16, 22, 25)]  %>%
        'colnames<-'(c("sample", "nano_conc", "protein", "salt", "tape_conc",
                       "conc_200", "dv200", "input_rna", "success"))

chemo <- chemo %>% mutate(., kit = "roche",
               tech = "megan",
               project = "wc_chemo",
               success = sub("MAYBE", "M", success))

## no duplicates to delete, no samples speed-vacced or pulled straight from
## tube


# Eva Roche ---------------------------------------------------------------
eva_roche <- read_excel("./b_datasets/megan_data/Eva AmpliSeq Lib Master.xlsx",
                        sheet = 1, skip = 2, col_names = T, na = "")

eva_roche <- eva_roche[,c(2, 3, 4, 5, 6, 7, 8, 14, 17)] %>%
        'colnames<-'(c("sample", "nano_conc", "protein", "salt", "tape_conc",
                       "conc_200", "dv200", "input_rna", "success"))

eva_roche <- eva_roche %>% mutate(., kit = "roche",
               tech = "megan",
               project = "eva_roche",
               success = sub("MAYBE", "M", success))

eva_roche <- remove_duplicates(eva_roche)


# Eva Qiagen --------------------------------------------------------------
eva_qiagen <- read_excel("./b_datasets/megan_data/Eva AmpliSeq Lib Master.xlsx",
                        sheet = 3, skip = 2, col_names = T, na = "")

findh20 <- eva_qiagen$H2O == 0

eva_qiagen <- eva_qiagen[,c(3, 4, 5, 6, 7, 8, 9, 10, 16, 19)] %>%
        'colnames<-'(c("sample", "type", "nano_conc", "protein", "salt",
                       "tape_conc", "conc_200", "dv200", "input_rna",
                       "success"))

eva_qiagen <- eva_qiagen %>% mutate(., kit = "qiagen",
               tech = "megan",
               project = "eva_qiagen",
               success = sub("[Mm].*", "M", success),
               input_rna = case_when(findh20 == T ~ "undiluted",
                                           TRUE ~ as.character(input_rna)))

eva_qiagen <- remove_duplicates(eva_qiagen)

eva_qiagen <- select(eva_qiagen, !type)


# Eva NE ------------------------------------------------------------------
eva_ne <- read_excel("./b_datasets/megan_data/Eva NE Lib Prep Master.xlsx",
                     sheet = 1, skip = 2, col_names = T, na = "")

findh20 <- eva_ne$H2O == 0

eva_ne <- eva_ne[,c(2, 3, 4, 5, 6, 8, 9, 10, 19, 22)] %>%
        'colnames<-'(c("sample", "type", "nano_conc", "protein", "salt",
                       "tape_conc", "conc_200", "dv200", "input_rna",
                       "success"))

eva_ne <- eva_ne %>%  mutate(., kit = NA,
               tech = "megan",
               project = "eva_ne",
               input_rna = case_when(findh20 == T ~ "undiluted",
                                           TRUE ~ as.character(input_rna)),
               success = sub("No", "N", success),
               success = sub("Maybe", "M", success))

eva_ne <- remove_duplicates(eva_ne)

eva_ne <- select(eva_ne, !type)



# SCC ---------------------------------------------------------------------
scc <- read_excel("./b_datasets/megan_data/SCC QC Master.xlsx", sheet = 2,
                  skip = 2, col_names = T, na = "")

scc <- scc[,c(2, 5, 6, 7, 9, 10, 11, 25, 28)] %>%
        'colnames<-'(c("sample", "nano_conc", "protein", "salt",
                       "tape_conc", "conc_200", "dv200", "input_rna",
                       "success"))

scc <- scc %>% mutate(., kit = "qiagen",
               tech = "megan",
               project = "scc",
               success = sub("Maybe", "M", success))

# WC BCG ------------------------------------------------------------------
bcg <- read_excel("./b_datasets/megan_data/WC BCG Pairs Ampliseq Lib Master.xlsx",
                  sheet = 2, skip = 2, col_names = T, na = "")

find <- bcg$RNA == 3.5

bcg <- bcg[,c(2, 7, 8, 9, 10, 11, 17, 20)] %>%
        'colnames<-'(c("sample", "nano_conc", "protein", "salt",
                       "tape_conc", "dv200", "input_rna",
                       "success"))

bcg <- bcg %>% mutate(., kit = "qiagen",
               tech = "megan",
               project = "scc",
               dv200 = dv200/100,
               input_rna = case_when(find == T ~ "undiluted",
                                     TRUE ~ as.character(input_rna)),
               success = sub("[Mm][Aa][Yy][Bb][Ee].*", "M", success),
               success = sub("^[Yy]es$", "Y", success),
               success = sub(".*<200bp.*|.*>200bp.*|MARKER", NA, success),
               success = sub("No", "N", success),
               conc_200 = NA)

bcg <- remove_duplicates(bcg)


# All Together Now --------------------------------------------------------
total <- rbind(a1, a2, ab, bcg, chemo, cis, cs, eva_ne, eva_qiagen,
               eva_roche, gbci, msk_qiagen, msk_roche, scc, wc)

nas <- data.frame(sample = sum(is.na(total$sample)),
                  nano_conc = sum(is.na(total$nano_conc)),
                  protein = sum(is.na(total$protein)),
                  salt = sum(is.na(total$salt)),
                  tape_conc = sum(is.na(total$tape_conc)),
                  conc_200 = sum(is.na(total$conc_200)),
                  dv200 = sum(is.na(total$dv200)),
                  input_rna = sum(is.na(total$input_rna)),
                  success = sum(is.na(total$success)),
                  kit = sum(is.na(total$kit)))

complete_total <- total[complete.cases(total),]

binary_complete <- complete_total[!complete_total$success == "M",]


# Clean Up ----------------------------------------------------------------
remove(a1, a2, ab, bcg, chemo, cis, cs, eva_ne, eva_qiagen, eva_roche, gbci,
       msk_qiagen, msk_roche, scc, wc, find, findh20)

# Create Train, Test, Validate --------------------------------------------
## want equal amounts of failed + succeeded for train dataset.
## NOT FINAL - JUST MESSING AROUND HERE.

set.seed(100)

train <- slice_sample(binary_complete, prop = 0.75) %>%
        mutate(., success_binary = case_when(success == "Y" ~ 1,
                                             success == "N" ~ 0))

train_n <- binary_complete[binary_complete$success == "N",]
set.seed(100)
train_y <- binary_complete[!binary_complete$success == "N",] %>%
        slice_sample(., n = nrow(train_n))
train <- rbind(train_n, train_y) %>% slice_sample(., prop = 1.0) %>%
        mutate(., success_binary = case_when(success == "Y" ~ 1,
                                             success == "N" ~ 0))
remove(train_n, train_y)




# Write Datasets for Later Use --------------------------------------------
write.table(complete_total, "complete_total.txt")
write.table(train, "training_data.txt")
write.table(binary_complete, "binary_complete.txt")

