### this preprocess code file redefines ICU encounter & sepsis encounter
### first define ICU encounter to be the 1st day with CAM
### then define sepsis encounter by applying Rhee's definition
### all patients are included first together and later keep sepsis ICU only 
### (on the level of encounter)

rm(list=ls())

library(dplyr)        # alternatively, this also loads %>%
library(tableone)
library(eeptools)    # age_calc to calculate age
library(readxl)
library(scales)
library(kableExtra)
library(magrittr)  # %<>%
library(runner)    # calculate on rolling/running windows

## haplogroup (analysis group)
h1_tmp = read_excel("../Mito Delirium BioVU Data/Genetics/CAM_Haplogroups.xlsx", sheet="White") %>%
    rename(grid=GRID, haplogroup=`analysis group`) %>%
    dplyr::select(grid, haplogroup) %>%
    mutate(race="White")
h2_tmp = read_excel("../Mito Delirium BioVU Data/Genetics/CAM_Haplogroups.xlsx", sheet="Black") %>%
    rename(grid=GRID, haplogroup=`analysis group`) %>%
    dplyr::select(grid, haplogroup) %>%
    mutate(race="Black")
h3_tmp = read_excel("../Mito Delirium BioVU Data/Genetics/CAM_Haplogroups.xlsx", sheet="HIspanic") %>%
    rename(grid=GRID, haplogroup=`analysis group`) %>%
    dplyr::select(grid, haplogroup) %>%
    mutate(race="Hispanic")
h4_tmp = read_excel("../Mito Delirium BioVU Data/Genetics/CAM_Haplogroups.xlsx", sheet="Drop") %>%
    rename(grid=GRID, haplogroup=`analysis group`) %>%
    dplyr::select(grid, haplogroup) %>%
    mutate(race="Drop")
haplo_group = rbind(h1_tmp, h2_tmp, h3_tmp, h4_tmp) %>%
    mutate(haplogroup = if_else(haplogroup=="L0", "Black Other", haplogroup)) %>%
    mutate(grid=as.character(grid),
           haplogroup=as.factor(haplogroup),
           race=as.factor(race))

stopifnot(length(unique(haplo_group$grid))==dim(haplo_group)[1])
rm(h1_tmp, h2_tmp, h3_tmp, h4_tmp)

## date of death (DOD)
dod_data = read_excel("../Mito Delirium BioVU Data/Demographics/Date of Death data.xlsx", 
                      sheet="FINAL DATE OF DEATH DATA") %>%
    rename(grid=GRID, dod=`DOD or proxy date`, dod_determine=`how DoD determined`,
           deceased=`edited DECEASED`) %>%
    mutate(dod=as.Date(dod), grid=as.character(grid)) %>%
    dplyr::select(grid, dod, dod_determine, deceased)

## all cohort subjects
all_grids = read_excel("../datafile/sepsis_grids_20200106.xlsx", sheet="All GRIDs") %>%
    filter(category == "Have Haplogroup") %>%   # only use grids having haplogroup
    dplyr::select(grid, gender, dob, ethnicity) %>%
    left_join(haplo_group, by="grid") %>%
    left_join(dod_data, by="grid")

rm(dod_data, haplo_group)

## adm & dc dates (include worst sofa score per admission)
adm_dates = read.csv("../datafile/sepsis_compare_20191217.csv") %>%
    mutate(grid = as.character(grid),
           adm_date = as.Date(adm_date),
           dc_date = as.Date(dc_date)) %>%
    rename(sofa_worst = sofa) %>%
    dplyr::select(grid, adm_date, dc_date, sofa_worst) %>%
    filter(grid %in% all_grids$grid)

## daily status
daily_status = read.csv("../datafile/daily_status_20190925.csv") %>%
    mutate(dt = as.Date(dt),
           grid = as.character(grid)) %>%
    dplyr::select(grid, dt, status.today) %>%
    filter(grid %in% all_grids$grid) %>%
    mutate(status.today = ifelse(status.today=="Unknown: conflicting CAM", "Delirious", status.today))

## daily lab (includes sofa, rass, creatinine, platelet ...)
daily_lab = read.csv("../datafile/daily_sofa_score_20191010.csv") %>%
    mutate(grid = as.character(grid)) %>%
    filter(grid %in% all_grids$grid) %>%
    rename(dt = lab_date) %>%
    mutate(dt = as.Date(dt)) %>%
    filter(day>=0)    # !!95 records have bugs of day<0; day=lab_date-adm_date+1

## create all days within an encounter
encounters = all_grids %>%
    left_join(adm_dates, by=c("grid")) %>%
    group_by(grid, adm_date) %>%
    do(data.frame(., dt=seq(.$adm_date, .$dc_date, by='1 day'))) %>% # final rows=223451
    ungroup() %>%
    left_join(daily_status, by=c("grid", "dt")) %>%
    left_join(daily_lab, by=c("grid", "dt")) %>%
    mutate(cam_yn = ifelse(status.today %in% c("Normal", "Delirious", "Comatose"), 1, 0),
           rass_yn = ifelse(is.na(rass), 0,1)) %>%    # 1-known 0-missing
    group_by(grid, adm_date) %>%
    mutate(cam_next = lead(cam_yn, 1),
           flag_end = case_when(
               dt!=dc_date & cam_yn==0 & cam_next==0 ~ 0,    # 2 consecutive unk cams
               dt==dc_date & cam_yn==0 & rass_yn==0 ~ 0,    # dc date have no cam & rass
               TRUE ~ 1    # 0 non-ICU day; 1 ICU day
           )) %>%
    ungroup() %>%
    # remove encounters with no cam
    group_by(grid, adm_date) %>%
    filter(1 %in% cam_yn) %>%
    ungroup()

rm(adm_dates, daily_status, daily_lab)

# define 1st ICU day - 1st day in an encounter with cam known
icu_start = encounters %>%
    group_by(grid, adm_date) %>%
    filter(cam_yn==1) %>%
    slice(c(1)) %>%
    ungroup() %>%
    dplyr::select(grid, adm_date, dt) %>%
    rename(icu_start_date=dt)
stopifnot(sum(is.na(icu_start$icu_start_date))==0)

# define last ICU day - 
# the day before 2 or more consecutive days without CAM, 
# or the day before a discharge day with neither CAM nor RASS
icu_end = encounters %>%
    dplyr::select(grid, adm_date, flag_end, dt) %>%
    left_join(icu_start, by=c("grid", "adm_date")) %>%
    filter(icu_start_date<=dt) %>%    # end day needs to be after start day
    group_by(grid, adm_date) %>%
    do (
        if (0 %in% .$flag_end) slice(., 1:(which.max(.$flag_end==0)-1))    # get first non-NA days group
        else slice(., 1:n())    # all observed, no 0
    ) %>%
    slice(c(n())) %>%    # then last non-NA day
    ungroup() %>%
    rename(icu_last_date = dt) %>%
    dplyr::select(-c(flag_end, icu_start_date))
stopifnot(sum(is.na(icu_end$icu_last_date))==0)
stopifnot(all(icu_start$icu_start_date <= icu_end$icu_last_date))    # all start dates<=end dates

# define ICU encounters
icu_daily = encounters %>%
    left_join(icu_start, by=c("grid", "adm_date")) %>%
    left_join(icu_end, by=c("grid", "adm_date"))

rm(icu_start, icu_end)

# create encounter level data
icu_encounters = icu_daily %>%
    group_by(grid, adm_date, dc_date, icu_start_date) %>%
    slice(c(1)) %>%
    ungroup() %>%
    dplyr::select(grid, adm_date, dc_date, icu_start_date)

#---- Apply Rhee's def to define sepsis encounters
# 1. old def: starts from adm date
source("rhee_infection_old.R")
rhee_infection_old = FUN_rhee_infection_old(icu_encounters)
source("rhee_organ_dysfunction_old.R")
sepsis_rhee_old = FUN_rhee_organ_old(rhee_infection_old)

# 2. new def: starts from 1st cam date
source("rhee_infection_new.R")
rhee_infection_new = FUN_rhee_infection_new(icu_encounters)
source("rhee_organ_dysfunction_new.R")
sepsis_rhee_new = FUN_rhee_organ_new(rhee_infection_new)

## save above partial results
# filename = paste0("preprocess_part1_define_sepsis_encnts.RData")
# save.image(filename)

################################################################################
########################## Saved part1 results #################################
################################################################################

#---- compare old def & new def
## 2-by-2 table for the old & new definitions (encounter level)
## old: icu encounters start from adm date
## new: icu encounters start from 1st cam date

# old - obtained from `rhee_infection_old.R` & `rhee_organ_dysfunction_old.R`
sepsis_rhee_old %<>% 
    dplyr::select(grid, adm_date, rhee)
sepsis_encnts_old = icu_daily %>%
    left_join(sepsis_rhee_old, by=c("grid", "adm_date")) %>%
    mutate(rhee = ifelse(is.na(rhee), 0, rhee))
rhee_old = sepsis_encnts_old %>%
    group_by(grid, adm_date) %>%
    slice(c(1)) %>%
    ungroup() %>%
    dplyr::select(grid, adm_date, rhee)
# new - obtained from `rhee_infection_new.R` & `rhee_organ_dysfunction_new.R`
sepsis_rhee_new %<>% 
    dplyr::select(grid, adm_date, rhee)
sepsis_encnts_new = icu_daily %>%
    left_join(sepsis_rhee_new, by=c("grid", "adm_date")) %>%
    mutate(rhee = ifelse(is.na(rhee), 0, rhee))
rhee_new = sepsis_encnts_new %>%
    group_by(grid, adm_date) %>%
    slice(c(1)) %>%
    ungroup() %>%
    dplyr::select(grid, adm_date, rhee)

rm(rhee_infection_old, sepsis_rhee_old, rhee_infection_new, sepsis_rhee_new)

## 2-by-2 table for old def vs new def
sepsis_compare = rhee_old %>%
    inner_join(rhee_new, by=c("grid", "adm_date")) %>%
    rename(rhee_old = rhee.x, rhee_new = rhee.y) 

# sepsis_compare %>% 
#     with(., addmargins(xtabs( ~ rhee_old + rhee_new, addNA = T))) %>% 
#     kable(caption="Table for sepsis ICU definition old vs. new") %>% 
#     kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = F)

icu_daily %<>%
    left_join(sepsis_compare, by=c("grid", "adm_date"))

#---- covariates ------------
## neurological exclusion
# neuro damage data (remove encounters with bad icd codes)
# ICD9/10 https://www.medicalbillingandcodingonline.com/icd-cm-codes/
bad_icd9 = c(430,431,432,851,852,853,854)  # icd9:xxx.xx 5 numeric digits
bad_icd10 = c("I60","I61","I62","I63","S02.0","S02.1","S02.91","S06","S07","S09.0","S09.8","S09.90")  # icd10:A00.000C 7digits
neuro_icd9 = read_excel("../Mito Delirium BioVU Data/Neuro Exclusions/neuro_icd9_V2.xlsx") %>%
    as.data.frame() %>%
    rename(grid=GRID, dt=CODE_DATE, icd9_code=CODE) %>%
    mutate(dt = as.Date(dt), code_group=floor(icd9_code),
           icd9_neuro=ifelse(code_group %in% bad_icd9, 1,0)) %>%  # 1: has damage; 0: no damage
    dplyr::select(grid, dt, icd9_code, icd9_neuro) %>%
    distinct(grid, dt, .keep_all = TRUE)
inverse_starts_with = function(standard, sample) {
    return(startsWith(sample, standard))
}
neuro_icd10 = read_excel("../Mito Delirium BioVU Data/Neuro Exclusions/neuro_icd10_V2.xlsx") %>%
    as.data.frame() %>%
    rename(grid=GRID, dt=CODE_DATE, icd10_code=CODE) %>%
    mutate(dt = as.Date(dt),
           icd10_neuro=ifelse(TRUE %in% sapply(bad_icd10, FUN=inverse_starts_with, icd10_code), 1,0)) %>% # 1:damage; 0:no damage
    dplyr::select(grid, dt, icd10_code, icd10_neuro) %>%
    distinct(grid, dt, .keep_all = TRUE)

## save above partial results
# filename = paste0("preprocess_part1.1_neuro_exclude.RData")
# save.image(filename)


################################################  comorbidity score ################################################## 
all_icd_files = list.files(path="../Mito Delirium BioVU Data/Elixhauser Comorbidities/", pattern=".xlsx$", full.names=T)
get_comorbid_name = function(file_name, sep1, sep2) {
    file_name_tmp = strsplit(file_name, sep1)[[1]][length(strsplit(file_name, sep1)[[1]])]  # take last item
    file_name_new = strsplit(file_name_tmp, sep2)[[1]][1]
    return(file_name_new)
}
comorbid_names = unlist(lapply(as.character(all_icd_files), FUN=get_comorbid_name, sep1="/", sep2=".xlsx"))
comorbid_data = list()
for (i in 1:length(all_icd_files)) {
    cat(i, comorbid_names[i], "\n")
    tmpdata = read_excel(all_icd_files[i], sheet=2)
    tmpdata[comorbid_names[i]] = 1
    comorbid_data[[i]] = na.omit(tmpdata)
}
merge_comorbid = comorbid_data %>%
    Reduce(function(x, y) full_join(x, y, by=c("GRID", "CODE_DATE")), .) %>%
    rename(grid = GRID,
           dt = CODE_DATE) %>%
    filter(grid %in% all_grids$grid) %>%
    mutate(dt = as.Date(dt)) %>%
    arrange(grid, dt) %>%  # sort
    replace(is.na(.), 0)

comorbid_grids = icu_daily %>%
    dplyr::select(grid, dt) %>%
    full_join(merge_comorbid, by=c("grid", "dt")) %>%
    arrange(grid, dt) %>%
    replace(is.na(.), 0)

# EI based on all previous history
EIscore_prev_mat = comorbid_grids %>%
    group_by(grid) %>%
    mutate_at(.vars = 3:33, .funs=cummax) %>%
    ungroup() %>%
    mutate(EIscore_prev = rowSums(.[3:33])) %>%
    right_join(icu_daily, by=c("grid", "dt")) %>%
    dplyr::select(grid, dt, EIscore_prev, all_of(comorbid_names))

# EI based on 1-year history
EIscore_1yr_mat = comorbid_grids %>%
    group_by(grid) %>%
    mutate_at(.vars = 3:33, .funs=list(~sum_run(x=., k=365, idx=dt))) %>%
    ungroup() %>%
    mutate_if(is.numeric, ~1 * (. > 0)) %>% # replace all >1 = 1
    mutate(EIscore_1yr = rowSums(.[3:33])) %>%
    right_join(icu_daily, by=c("grid", "dt")) %>%
    dplyr::select(grid, dt, EIscore_1yr) #, all_of(comorbid_names)

## preexisting dementia (look at history of a subject <= adm)
daily_dementia = read_excel("../Mito Delirium BioVU Data/Dementia/Dementia.xlsx", sheet="Initial Dementia Code date") %>%
    bind_cols(dementia_raw=rep(1, nrow(.))) %>%
    rename(grid=GRID, dt=CODE_DATE) %>%
    mutate(dt = as.Date(dt)) %>%
    right_join(icu_daily, by=c("grid", "dt")) %>%
    mutate(dementia_raw = ifelse(is.na(dementia_raw), 0,1)) %>%
    dplyr::select(grid, dt, dementia_raw) %>%
    group_by(grid) %>%  # look at previous history of a subject
    mutate(dementia_adm = cumsum(dementia_raw)) %>%
    ungroup()

## drugs
tmp_drugs1 = read.csv("../Mito Delirium BioVU Data/Data/grid_date_med1.csv") %>%
    as.data.frame() %>%
    dplyr::select(GRID, DRUG_DATE, DRUG_CLASS)
tmp_drugs2 = read.csv("../Mito Delirium BioVU Data/Data/grid_date_med2.csv") %>%
    as.data.frame() %>%
    dplyr::select(GRID, DRUG_DATE, DRUG_CLASS)
tmp_drugs3 = read.csv("../Mito Delirium BioVU Data/Data/grid_date_med3.csv") %>%
    as.data.frame() %>%
    dplyr::select(GRID, DRUG_DATE, DRUG_CLASS)
daily_drugs = rbind(tmp_drugs1, tmp_drugs2, tmp_drugs3) %>%
    filter(DRUG_CLASS %in% c("opiate", "benzo", "propofol", "dexmedetomidine")) %>%
    rename(grid=GRID, dt=DRUG_DATE) %>%
    mutate(dt = as.Date(dt),
           grid = as.character(grid)) %>%
    group_by(grid, dt, DRUG_CLASS) %>%
    summarise(count = n())  %>%
    tidyr::spread(DRUG_CLASS, count, fill = 0) %>%  # long to wide format
    ungroup() %>%
    mutate_if(is.numeric, ~1 * (. > 0)) # >0 to be 1
rm(tmp_drugs1, tmp_drugs2, tmp_drugs3)


## merge daily
daily_mat = icu_daily %>%
    left_join(daily_dementia, by=c("grid", "dt")) %>%    # same rows
    left_join(EIscore_prev_mat, by=c("grid", "dt")) %>%  # same rows
    left_join(EIscore_1yr_mat, by=c("grid", "dt")) %>%   # same rows
    left_join(daily_drugs, by=c("grid", "dt")) %>%
    left_join(neuro_icd9, by=c("grid", "dt")) %>%
    left_join(neuro_icd10, by=c("grid", "dt")) %>%
    tidyr::replace_na(list(icd9_neuro=0, icd10_neuro=0)) %>%
    mutate(neuro_yn=ifelse(icd9_neuro==0 & icd10_neuro==0, 0,1),
           neuro_yn=factor(neuro_yn, levels=c(0,1), labels=c("no neuro issues", "had neuro issues")),
           female = ifelse(gender=="F", 1, 0),
           # have.haplogroup = ifelse(have.haplogroup=="yes",1,0),
           dob = as.Date(dob),
           age_adm = as.numeric(adm_date - dob) / 365,
           stay_duration = as.numeric(age_calc(adm_date, dc_date, units="days"))+1,
           status.today = factor(status.today),
           pressor_bin = ifelse(is.na(pressor), 0,1),
           bilirubin_bin = ifelse(is.na(bilirubin), 0,1),
           status.today = forcats::fct_explicit_na(status.today, "No CAM or RASS"),
           diff_2EIscores = EIscore_prev - EIscore_1yr,
           vent = ifelse(is.na(vent), 0, vent),
           opiate = ifelse(is.na(opiate), 0, opiate),
           benzo = ifelse(is.na(benzo), 0, benzo),
           propofol = ifelse(is.na(propofol), 0, propofol),
           dexmedetomidine = ifelse(is.na(dexmedetomidine), 0, dexmedetomidine),
           dementia_adm = ifelse(is.na(dementia_adm), 0, dementia_adm),
           dementia_adm = factor(dementia_adm)) %>%
    # remove encounter with any neuro damage
    group_by(grid, adm_date) %>%
    mutate(neuro_days = sum(neuro_yn=="has neuro issues")) %>%  # get sum of neuro_yn by encounter
    ungroup() %>%
    filter(neuro_days == 0) %>%  
    # add num of days per encounter (grouped by grid+adm_date)
    add_count(grid, adm_date, name="num_days") %>%
    # add num of unknown days per encounter
    group_by(grid, adm_date) %>%
    mutate(unknown_days = sum(status.today %in% c("Unknown: RASS only", "Unknown: conflicting CAM", "No CAM or RASS")),
           unknown_prop = unknown_days / num_days) %>%
    ungroup()

rm(daily_dementia, EIscore_prev_mat, EIscore_1yr_mat, daily_drugs, neuro_icd9, neuro_icd10)

sepsis_dailymat = daily_mat %>%
    filter(rhee_new == 1)  # 54877 days, 3953 encounters
nosepsis_dailymat = daily_mat %>%
    filter(rhee_new == 0)  # 178988 days, 21523 encounters
stopifnot(dim(sepsis_dailymat)[1] + dim(nosepsis_dailymat)[1] == dim(daily_mat)[1])

## save sepsis_dailymat,nosepsis_dailymat
# save(sepsis_dailymat,nosepsis_dailymat,file = paste0("preprocess_dailymat_raw.RData"))

## save above partial results
# filename = paste0("preprocess_part1.2_sepsis_nonsepsis_dailymat_raw.RData")
# save.image(filename)

