### this reads results generated from `preprocess_part1.2_sepsis_nonsepsis_dailymat_raw.RData`

rm(list=ls())

library(dplyr)        # alternatively, this also loads %>%
library(tableone)
library(eeptools)    # age_calc to calculate age
library(readxl)
library(scales)
library(kableExtra)
library(magrittr)  # %<>%
library(runner)    # calculate on rolling/running windows

# sepsis_dailymat,nosepsis_dailymat
load("preprocess_dailymat_raw.RData")

# sepsis_dailymat <- nosepsis_dailymat

#### impute daily status ----
## check pattern by tables
nk_str = c("Comatose", "Delirious", "Normal")
unk_str = c("Unknown: RASS only", "No CAM or RASS")

cam = sepsis_dailymat %>%
    dplyr::select(grid, adm_date, dc_date, icu_start_date, icu_last_date, dt, status.today, rass) %>%
    mutate(status.today = ifelse(status.today %in% unk_str, NA, as.character(status.today))) %>%
    group_by(grid, adm_date) %>%
    mutate(yesterday = lag(status.today, 1),
           tomorrow = lead(status.today, 1)) %>%
    ungroup()

## impute
# 1. a single middle day with t-1 & t+1 same status
s1 = cam %>%
    filter(dt!=dc_date & yesterday==tomorrow & is.na(status.today)) %>%
    mutate(status.today_new = yesterday) %>%
    dplyr::select(grid, dt, status.today_new)

## other scenario: borrow info from rass & t-1 cam & t+1 cam
# 2.1 last day with known rass
FUN_draw_probs2 = function(yesterday_value, rass_value, approach=2) {
    df_tmp = cam %>%
        filter(yesterday==yesterday_value & !is.na(status.today) & dt==dc_date & rass==rass_value)
    if (dim(df_tmp)[1]==0) {  # no records for day t-1 comatose & day t rass=-3
        status.today_new = yesterday_value
        cat("Unseen combinations: yesterday =", yesterday_value, "; rass_today =", rass_value, "\n")
        return (status.today_new)
    }
    coma_prob = sum(df_tmp$status.today=="Comatose") / dim(df_tmp)[1]
    deli_prob = sum(df_tmp$status.today=="Delirious") / dim(df_tmp)[1]
    norm_prob = sum(df_tmp$status.today=="Normal") / dim(df_tmp)[1]
    nk_str = c("Comatose", "Delirious", "Normal")
    if (approach==2) {  # assign mode
        status.today_new = nk_str[which.max(c(coma_prob, deli_prob, norm_prob))]
    }
    else {  # sample from prob dist
        status.today_new = sample(c("Comatose", "Delirious", "Normal"), 1, 
                                  prob=c(coma_prob, deli_prob, norm_prob))
    }
    return (status.today_new)
}

set.seed(1234)
s2 = cam %>%
    filter(dt==dc_date & dc_date==icu_last_date & is.na(status.today) & !is.na(rass)) %>%
    rowwise() %>%  # conduct each row operation
    mutate(status.today_new = FUN_draw_probs2(yesterday, rass)) %>%
    ungroup() %>%
    dplyr::select(grid, dt, status.today_new)


# 3.1 middle day with known rass but conflict before & after
FUN_draw_probs3 = function(yesterday_value, tomorrow_value, rass_value, approach=2) {
    df_tmp = cam %>%
        filter(yesterday==yesterday_value & !is.na(status.today) & 
                   tomorrow==tomorrow_value & rass==rass_value)
    if (dim(df_tmp)[1]==0) {  # no records for day t-1 comatose & day t rass=-3
        status.today_new = "Normal" #yesterday_value  # all impute as normal instead of yesterday
        cat("Unseen combinations: yesterday =", yesterday_value, "; rass_today =", rass_value, 
            "; tomorrow =", tomorrow_value, "\n")
        return (status.today_new)
    }
    coma_prob = sum(df_tmp$status.today=="Comatose") / dim(df_tmp)[1]
    deli_prob = sum(df_tmp$status.today=="Delirious") / dim(df_tmp)[1]
    norm_prob = sum(df_tmp$status.today=="Normal") / dim(df_tmp)[1]
    if (approach==2) {  # assign mode
        status.today_new = nk_str[which.max(c(coma_prob, deli_prob, norm_prob))]
    }
    else {  # sample from prob dist
        status.today_new = sample(c("Comatose", "Delirious", "Normal"), 1, 
                                  prob=c(coma_prob, deli_prob, norm_prob))
    }
    return (status.today_new)
}

set.seed(1234)
s3 = cam %>%
    filter(dt!=dc_date & icu_start_date<=dt & dt<=icu_last_date & is.na(status.today) & !is.na(rass) 
           & yesterday!=tomorrow) %>%
    rowwise() %>%  # conduct each row operation
    mutate(status.today_new = FUN_draw_probs3(yesterday, tomorrow, rass)) %>%
    ungroup() %>%
    dplyr::select(grid, dt, status.today_new)


# 3.2 middle day with unknown rass and conflict before & after
FUN_draw_probs4 = function(yesterday_value, tomorrow_value, approach=2) {
    df_tmp = cam %>%
        filter(yesterday==yesterday_value & !is.na(status.today) & 
                   tomorrow==tomorrow_value)
    if (dim(df_tmp)[1]==0) {  # no records for day t-1 comatose & day t rass=-3
        status.today_new = yesterday_value
        cat("Unseen combinations: yesterday =", yesterday_value, "; tomorrow =", tomorrow_value, "\n")
        return (status.today_new)
    }
    coma_prob = sum(df_tmp$status.today=="Comatose") / dim(df_tmp)[1]
    deli_prob = sum(df_tmp$status.today=="Delirious") / dim(df_tmp)[1]
    norm_prob = sum(df_tmp$status.today=="Normal") / dim(df_tmp)[1]
    if (approach==2) {  # assign mode
        status.today_new = nk_str[which.max(c(coma_prob, deli_prob, norm_prob))]
    }
    else {  # sample from prob dist
        status.today_new = sample(c("Comatose", "Delirious", "Normal"), 1, 
                                  prob=c(coma_prob, deli_prob, norm_prob))
    }
    return (status.today_new)
}

set.seed(1234)
s4 = cam %>%
    filter(dt!=dc_date & icu_start_date<=dt & dt<=icu_last_date & is.na(status.today) & is.na(rass)
           & yesterday!=tomorrow) %>%
    rowwise() %>%  # conduct each row operation
    mutate(status.today_new = FUN_draw_probs4(yesterday, tomorrow)) %>%
    ungroup() %>%
    dplyr::select(grid, dt, status.today_new)

status_new_mat = cam %>%
    left_join(s1, by=c("grid","dt")) %>%
    left_join(s2, by=c("grid","dt")) %>%
    left_join(s3, by=c("grid","dt")) %>%
    left_join(s4, by=c("grid","dt")) %>%
    mutate(status.today_new = coalesce(status.today, status.today_new.x, status.today_new.y,
                                       status.today_new.x.x, status.today_new.y.y)) %>%
    dplyr::select(grid, dt, status.today_new)

# add status.today_new to sepsis_dailymat
sepsis_dailymat %<>%
    left_join(status_new_mat, by=c("grid","dt"))

rm(s1, s2, s3, s4, status_new_mat)

#### impute SOFA ----
# liver sofa: last obs carried forward
# impute 0 if all missing in an encounter
library(zoo)
sepsis_dailymat %<>%
    mutate(sofa_liver_new = sofa_liver) %>%
    group_by(grid, adm_date) %>%
    tidyr::fill(sofa_liver_new, .direction="downup") %>%
    ungroup() %>%
    mutate(sofa_liver_new = ifelse(is.na(sofa_liver_new), 0, sofa_liver_new))

# Cardiovascular sofa, Coagulation sofa, Renal sofa
# impute a single unknown day with t-1 & t+1 same status
FUN_impute_single = function(df, today_var) {
    today_new = df %>%
        group_by(grid, adm_date) %>%
        mutate(yesterday = lag(!!as.name(today_var), 1),
               tomorrow = lead(!!as.name(today_var), 1),
               today_new = if_else(yesterday==tomorrow & is.na(!!as.name(today_var)), 
                                   yesterday, !!as.name(today_var))) %>%
        ungroup() %>%
        dplyr::select(grid, dt, today_new)
    return(today_new)  # type: tbl_df
}

# for now only impute one single unknown day
sofa_cardio_new = FUN_impute_single(sepsis_dailymat, "sofa_cardio") %>%
    rename(sofa_cardio_new = today_new)
sofa_coagul_new = FUN_impute_single(sepsis_dailymat, "sofa_coagulation") %>%
    rename(sofa_coagul_new = today_new)
sofa_renal_new = FUN_impute_single(sepsis_dailymat, "sofa_renal") %>%
    rename(sofa_renal_new = today_new)

# add sofa_new to sepsis_dailymat
sepsis_dailymat %<>%
    left_join(sofa_cardio_new, by=c("grid","dt")) %>%
    left_join(sofa_coagul_new, by=c("grid","dt")) %>%
    left_join(sofa_renal_new, by=c("grid","dt")) %>%
    group_by(grid, adm_date) %>%
    mutate(status.pre = lag(status.today_new, 1),
           benzo.pre = lag(benzo, 1),
           opiate.pre = lag(opiate, 1),
           propofol.pre = lag(propofol, 1), 
           dex.pre = lag(dexmedetomidine, 1)) %>%
    ungroup()

rm(sofa_cardio_new, sofa_coagul_new, sofa_renal_new)

## save above partial results
# filename = paste0("preprocess_part2.1_sepsis_dailymat_imputed.RData")
# save.image(filename)
# filename = paste0("preprocess_part2.2_nosepsis_dailymat_imputed.RData")
# save.image(filename)

## move to preprocess_redefine_14days.R if want the first 14 days
## else, move on

################################################################################
########################## Saved part2 results #################################
################################################################################

## sepsis ICU encounters (only keep from ICU start date to ICU end date)
sepsisICU = sepsis_dailymat %>%
    group_by(grid, adm_date) %>%
    filter(icu_start_date<=dt & dt<=icu_last_date) %>%
    ungroup() %>%
    add_count(grid, adm_date, name="num_icu_days") %>%
    group_by(grid, adm_date) %>%
    do(data.frame(., na_cnt_cardio = sum(is.na(.$sofa_cardio_new)),
                  na_cnt_coagul = sum(is.na(.$sofa_coagul_new)),
                  na_cnt_renal = sum(is.na(.$sofa_renal_new)))) %>%
    ungroup() %>%
    # remove encounters with all missing (cardio, coagul, renal) sofa
    filter(na_cnt_cardio < num_icu_days) %>%
    filter(na_cnt_coagul < num_icu_days) %>%
    filter(na_cnt_renal < num_icu_days)

## impute cardio/coagul/renal sofa with *linear interpolation*
sepsisICU %<>%
    group_by(grid, adm_date) %>%
    mutate(sofa_cardio_imp = na.approx(sofa_cardio_new, rule = 2),
           sofa_coagul_imp = na.approx(sofa_coagul_new, rule = 2),
           sofa_renal_imp = na.approx(sofa_renal_new, rule = 2)) %>%
    ungroup() %>%
    mutate(sofa_cardio_imp = round(sofa_cardio_imp),
           sofa_coagul_imp = round(sofa_coagul_imp),
           sofa_renal_imp = round(sofa_renal_imp)) %>%
    mutate(
        sofa_new = case_when(
            is.na(sofa_respiration) + is.na(sofa_coagul_imp) + is.na(sofa_liver_new) + 
                is.na(sofa_cardio_imp) + is.na(sofa_cns) + is.na(sofa_renal_imp) != 6 ~
                coalesce(sofa_liver_new, 0) + coalesce(sofa_coagul_imp, 0) + 
                coalesce(sofa_cns, 0) + coalesce(sofa_cardio_imp, 0) + 
                coalesce(sofa_respiration, 0) + coalesce(sofa_renal_imp, 0)
        )) %>%
    mutate(sofa_cardio_imp = factor(sofa_cardio_imp),
           sofa_coagul_imp = factor(sofa_coagul_imp),
           sofa_renal_imp = factor(sofa_renal_imp),
           sofa_liver_new = factor(sofa_liver_new),
           vent = factor(vent))

## sofa at first icu day for each encounter
sepsisICU %<>%
    dplyr::group_by(grid, adm_date) %>%
    arrange(dt) %>%
    mutate(sofa_icu_start = sofa_new[1]) %>%  # sepsisICU starts from first ICU day
    ungroup()

## sofa at admission for each encounter
# daily_mat %<>% 
#   group_by(grid, adm_date) %>%
#   slice(c(1)) %>%  # get the first day of each encounter
#   ungroup() %>%
#   dplyr::select(grid, adm_date, sofa) %>%
#   rename(sofa_adm = sofa) %>%
#   right_join(daily_mat, by=c("grid", "adm_date"))

stopifnot(sum(is.na(sepsisICU$sofa_cardio_imp)) == 0)
stopifnot(sum(is.na(sepsisICU$sofa_coagul_imp)) == 0)
stopifnot(sum(is.na(sepsisICU$sofa_renal_imp)) == 0)


## in-icu death: encounter-level
tmp_dod = sepsisICU %>%
    group_by(grid, adm_date) %>%  # add adm_date to make it into encounter-level
    slice(n()) %>%  # get the last icu date of a patient seen in all encounters
    ungroup() %>%
    dplyr::select(grid, adm_date, dc_date, icu_last_date, dod, deceased) %>%
    mutate(in_icu_death = case_when(
        icu_last_date<dc_date ~ 0,
        dod<=icu_last_date ~ 1,
        icu_last_date==dc_date & deceased=="null" ~ 0,
        dod>icu_last_date ~ 0
    )) %>%
    mutate(in_icu_death = factor(in_icu_death)) %>%
    dplyr::select(grid, adm_date, in_icu_death)

death_manual <- read.delim("../datafile/in_ICU_death_manual_review.txt")
for (i in 1:dim(death_manual)[1]) {
    tmp_dod$in_icu_death[which(tmp_dod$grid==death_manual$grid[i])] = ifelse(
        death_manual$manual.review.in.ICU.death[i]=="Yes", 1, 0)
}
# R206551853 - alive since had drugs on 10/10/12
tmp_dod$in_icu_death[which(tmp_dod$grid=="R206551853")] = 0
## check last date in system
# which(is.na(tmp_dod$in_icu_death))
# sepsisICU %>%
#   filter(grid=="R257565449") %>%
#   dplyr::select(grid, adm_date, dc_date, icu_last_date)
tmp_dod$in_icu_death[which(tmp_dod$grid=="R220237780")] = 0
tmp_dod$in_icu_death[which(tmp_dod$grid=="R248606672")] = 0
tmp_dod$in_icu_death[which(tmp_dod$grid=="R257565449")] = 0

stopifnot(sum(is.na(tmp_dod$in_icu_death)) == 0)


sepsisICU %<>%
    left_join(tmp_dod, by=c("grid", "adm_date"))  # add dod to sepsisICU
rm(tmp_dod)

stopifnot(sepsisICU$status.today_new %in% c("Normal", "Delirious", "Comatose"))

## analysis mat: include normal/delirious/comatose daily records
# after impute, all status should be known
analysis_mat = sepsisICU %>%
    filter(status.today_new %in% c("Normal", "Delirious", "Comatose")) %>%  # after impute
    add_count(grid, adm_date, name="num_kn_days") %>%  # num of known days
    mutate(is_delirium = ifelse(status.today_new=="Delirious", 1, 0),
           is_coma = ifelse(status.today_new=="Comatose", 1, 0)) %>%
    group_by(grid, adm_date) %>%
    mutate(day_id = row_number(),
           is_delirium_ever = ifelse("Delirious" %in% status.today_new, 1, 0),
           is_coma_ever = ifelse("Comatose" %in% status.today_new, 1, 0),
           is_delirium_ever = factor(is_delirium_ever),
           is_coma_ever = factor(is_coma_ever)) %>%  # number days within encounters
    do(data.frame(., day_benzo = sum(.$benzo), 
                  day_opiate = sum(.$opiate),
                  day_propofol = sum(.$propofol),
                  day_dex = sum(.$dexmedetomidine),
                  day_delirium = sum(.$is_delirium), 
                  day_coma = sum(.$is_coma))) %>%
    ungroup() %>%
    mutate(day_delirium = ifelse(day_delirium==0, NA, day_delirium),
           day_coma = ifelse(day_coma==0, NA, day_coma))  # only need days when day>=1


# daily-level
analysis_daily = analysis_mat %>%
    # filter(race=="White" | race=="Black") %>% # only use white & black #UPDATED on 2022-06-16: not filter race
    mutate(
        # race = forcats::fct_drop(race),  # remove null levels
        # race = forcats::fct_relevel(race, "White", "Black"),
        # haplogroup = forcats::fct_drop(haplogroup),
        # haplogroup = forcats::fct_relevel(haplogroup, "H", "I", "W", "X", "J", "T", "UK", "White Other", 
        #                                   "L3", "L1", "L2", "Black Other"),
        is_WhiteOther = ifelse(haplogroup=="White Other", 1, 0),
        is_H = ifelse(haplogroup=="H", 1, 0),
        is_I = ifelse(haplogroup=="I", 1, 0),
        is_W = ifelse(haplogroup=="W", 1, 0),
        is_X = ifelse(haplogroup=="X", 1, 0),
        is_J = ifelse(haplogroup=="J", 1, 0),
        is_T = ifelse(haplogroup=="T", 1, 0),
        is_UK = ifelse(haplogroup=="UK", 1, 0),
        is_BlackOther = ifelse(haplogroup=="Black Other", 1, 0),
        is_L1 = ifelse(haplogroup=="L1", 1, 0),
        is_L2 = ifelse(haplogroup=="L2", 1, 0),
        is_L3 = ifelse(haplogroup=="L3", 1, 0))  # re-arrange levels

daily_model_deli = analysis_daily %>%
    filter(!(adm_date==icu_start_date & dt==adm_date)) %>%
    filter(status.today_new!="Comatose") %>%
    mutate(delirium = ifelse(status.today_new=="Delirious", 1, 0))

daily_model_coma = analysis_daily %>%
    filter(!(adm_date==icu_start_date & dt==adm_date)) %>%
    filter(status.today_new!="Delirious") %>%
    mutate(comatose = ifelse(status.today_new=="Comatose", 1, 0))

# encounter-level
analysis_adm = analysis_daily %>%
    group_by(grid, adm_date) %>%
    slice(c(1)) %>%  # get the first ICU day
    ungroup()

analysis_enc <- analysis_daily %>% 
    filter(status.today_new != "Comatose") %>%  # remove coma for now
    add_count(grid, adm_date, name="num_icu_days_wo_coma") %>%
    group_by(grid, adm_date) %>%
    mutate(enc_id = cur_group_id(),  # encounter id (cluster)
           vent_ever = ifelse(1 %in% vent, 1, 0),  # whether vent or not within an encounter
           mean_cardio_sofa = mean(as.numeric(as.character(sofa_cardio_imp))),
           mean_coagul_sofa = mean(as.numeric(as.character(sofa_coagul_imp))),
           mean_liver_sofa = mean(as.numeric(as.character(sofa_liver_new))),
           mean_renal_sofa = mean(as.numeric(as.character(sofa_renal_imp))),
           day_delirium=ifelse(is.na(day_delirium), 0, day_delirium),
           prop_delirium_wocoma = day_delirium / num_icu_days_wo_coma,
           prop_delirium_wocoma_bin = ifelse(prop_delirium_wocoma>=0.5, 1, 0)) %>%
    slice(1) %>% 
    ungroup()

# subject-level  # don't plan to use subject-level covars; only use encounter & daily
analysis_subj = analysis_daily %>%
    group_by(grid) %>%
    slice(c(1)) %>%
    ungroup()

## post-preprocessing: check consecutive unknown patterns
source("consecutive_unk_seq.R")
tmp_cardio_res = FUN_unk_seq(sepsisICU, "sofa_cardio_new")
tmp_coagul_res = FUN_unk_seq(sepsisICU, "sofa_coagul_new")
tmp_renal_res = FUN_unk_seq(sepsisICU, "sofa_renal_new")
# missing patterns (0 as known and non-zero value as consecutive unknown days)
unk_seq_cardio = tmp_cardio_res$consec_unkdays
unk_seq_coagul = tmp_coagul_res$consec_unkdays
unk_seq_renal = tmp_renal_res$consec_unkdays
# daily sofa for encounters with unknown days
unk_daily_cardio = tmp_cardio_res$myvar_daily
unk_daily_coagul = tmp_coagul_res$myvar_daily
unk_daily_renal = tmp_renal_res$myvar_daily

# ## check encounters with consecutive missing
# myvar = "sofa_renal_new"
# a = sepsisICU %>%
#   dplyr::select(grid, adm_date, icu_start_date, icu_last_date, dt, !!as.name(myvar)) %>%
#   group_by(grid, adm_date) %>%
#   filter(NA %in% !!as.name(myvar)) %>%
#   slice(n()) %>%
#   ungroup()
# 
# b = sepsisICU %>%
#   dplyr::select(grid, adm_date, icu_start_date, icu_last_date, dt, !!as.name(myvar)) %>%
#   group_by(grid, adm_date) %>%
#   filter(NA %in% !!as.name(myvar)) %>%
#   ungroup()
# 
# c = b %>%
#   group_by(grid, adm_date) %>%
#   filter(4 %in% !!as.name(myvar)) %>%
#   filter(1 %in% !!as.name(myvar)) %>%
#   ungroup()
# 
# ## encounters for David to look into
# # for CV sofa
# R298483492, 2010-06-03  # 0-10-3, 3-5-0
# 
# # for coa sofa
# R211962586, 2007-06-21  # 0-6-0
# R296672677, 2006-11-04  # 0-1-4
# 
# # for renal sofa
# R211232481, 2011-11-26  # last 11 days missing
# R293693829, 2012-04-03  # 1-1-4

# filename = paste0("preprocess_2020-11-17.RData")
# save.image(filename)

analysis_daily2 <- analysis_daily %>% 
    dplyr::select(-adm_id, -day, -map, -data_type, -flag_end, -sofa_cardio_new, 
                  -sofa_coagul_new, -sofa_renal_new)

################################################################################
########################## Saved final results #################################
################################################################################

# data.table::fwrite(analysis_daily2, file = "sepsis_icu_daily_imputed_2022-06-16.csv")
# data.table::fwrite(analysis_daily2, file = "nosepsis_icu_daily_imputed_2022-06-16.csv")
