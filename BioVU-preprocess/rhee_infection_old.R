### this file is included in `preprocess_redefine.R`
### to apply Rhee's definition to get sepsis ICU

library(tidyverse)
library(readxl)
library(magrittr)
library(lubridate)

FUN_rhee_infection_old = function(cam_visits) {

  #. Data import and clean ----------------------------------------------------
  static_raw <- read_csv("../Mito Delirium BioVU Data/Data/Samuels_Delirium_STATIC_20180718.csv")
  changed_grid <- read_csv("../datafile/changed_grid_dob_20190924.csv")
  # cam_visits <- read_csv("../datafile/icu_visits_20200609.csv")  # xtan: use ICU_start_date not adm_date
  blood_raw <- read_excel("../Mito Delirium BioVU Data/Phenotype data/culture_merge.xlsx",
                          sheet = "Blood Culture Days")
  
  names(static_raw) <- str_to_lower(names(static_raw))
  names(blood_raw) <- str_to_lower(names(blood_raw))
  names(changed_grid) <- str_to_lower(names(changed_grid))
  
  med_raw <- NULL
  for (i in 1:3) {
    med_raw %<>% 
      bind_rows(
        read_csv(
          paste0("../Mito Delirium BioVU Data/Data/grid_date_med", i, ".csv"),
          col_names = c("grid", "drug_date", "drug_name", "drug_class", 
                        "drug_route1", "drug_route2", "drug_route3"),  # could have multiple routes (tho mostly single)
          skip = 1
        ) 
      ) 
  }
  # med_raw %>% count(drug_class)
  abx_raw <- med_raw %>% 
    filter(drug_class == "antibiotic") %>% 
    select(-drug_class) 
  
  #> clean up blood culture data -------------------------------------------------
  ## convert messed-up GRIDs and dates 
  # length(unique(blood_raw$grid))
  # sum(unique(blood_raw$grid) %in% unique(static_raw$grid)) 
  # sum(unique(blood_raw$grid) %in% unique(changed_grid$old_grid))
  # sum(unique(blood_raw$grid) %in% unique(changed_grid$updated_grid))
  
  blood_raw1 <- blood_raw %>% 
    left_join(changed_grid, by = c("grid" = "old_grid")) %>% 
    mutate(
      blood_date = if_else(!is.na(updated_grid),
                           as_date(`blood culture code_date`) - old_dob + updated_dob,
                           as_date(`blood culture code_date`)),
      grid = if_else(!is.na(updated_grid), updated_grid, grid)
    ) %>% 
    distinct(grid, blood_date) %>% 
    arrange(grid, blood_date) # some duplicates
  
  ## blood culture day within 1 day of hospital admissiond date 
  blood_raw2 <- sqldf::sqldf('SELECT * 
                           FROM cam_visits as t1
                           INNER JOIN blood_raw1 as t2
                           ON t1.grid = t2.grid AND blood_date <= dc_date AND blood_date BETWEEN adm_date-1 AND adm_date+1') %>%
    setNames(make.unique(names(.))) %>%  # make unique column names
    as_tibble() %>% 
    select(grid, adm_date, dc_date, icu_start_date, blood_date) %>% 
    mutate(day = as.numeric(blood_date - adm_date) + 1) %>%   # xtan: from blood_date to icu_start_date
    select(grid, adm_date, dc_date, icu_start_date, blood_date, day) 
  
  # blood_raw2 %>%
  #   group_by(grid, blood_date) %>%
  #   summarise(n = n_distinct(adm_date)) %>%
  #   filter(n > 1) # no overlap
  # blood_raw2 %>% 
  #   distinct(grid, adm_id)
  
  #> clean up antibiotic data and find QADs ---------------------------------------
  # length(unique(abx_raw$grid))
  # sum(unique(abx_raw$grid) %in% unique(static_raw$grid))  # check whether any new updated GRID
  # sum(unique(abx_raw$grid) %in% unique(changed_grid$old_grid))
  # sum(unique(abx_raw$grid) %in% unique(changed_grid$updated_grid))
  # abx_raw %>% 
  #   count(drug_name) %>% 
  #   print(n = 100)
  # abx_raw %>% 
  #   count(drug_route1, drug_route2, drug_route3)
  # abx_raw %>% 
  #   filter(drug_name == 'VANCOMYCIN') %>% 
  #   count(drug_route1, drug_route2, drug_route3)
  
  ## convert messed-up GRIDs and simplify drug routes
  abx_raw1 <- abx_raw %>% 
    left_join(changed_grid, by = c("grid" = "old_grid")) %>% 
    mutate(
      drug_date = if_else(!is.na(updated_grid),
                          as.Date(drug_date) - old_dob + updated_dob,
                          as.Date(drug_date)),
      grid = if_else(!is.na(updated_grid), updated_grid, grid),
      ivm = if_else(drug_route1 %in% c("IV", "IM") | drug_route2 %in% c("IV", "IM") | drug_route3 %in% c("IV", "IM"), 1, 0),
      po = if_else(drug_route1 %in% 'PO' | drug_route2 %in% 'PO' | drug_route3 %in% 'PO', 1, 0)
    ) %>% 
    distinct(grid, drug_date, drug_name, drug_route1, drug_route2, drug_route3, ivm, po) %>% 
    arrange(grid, drug_date) # some duplicates
  
  # abx_raw1 %>% 
  #   count(ivm, drug_route1, drug_route2, drug_route3)
  # abx_raw1 %>% 
  #   count(po, drug_route1, drug_route2, drug_route3)
  # abx_raw1 %>% 
  #   count(po, ivm)
  # abx_raw1 %>% 
  #   filter(drug_name == 'VANCOMYCIN') %>% 
  #   count(po, ivm)
  # abx_raw1 %>% 
  #   distinct(grid, drug_date, drug_name, ivm, po)
  
  ## filter by GRID and identify new QAD
  abx_raw2 <- abx_raw1 %>% 
    semi_join(blood_raw2, by = "grid") %>% 
    distinct(grid, drug_date, drug_name, ivm, po)
  
  abx_raw3 <- sqldf::sqldf('SELECT *
        FROM abx_raw2 as t1
        left join abx_raw2 as t2
        ON t1.grid = t2.grid AND t1.drug_name == t2.drug_name AND t2.drug_date BETWEEN t1.drug_date-2 AND t1.drug_date-1') %>%
    setNames(make.unique(names(.))) %>%  # make unique column names
    as_tibble() %>% 
    mutate(new_abx = case_when(
      drug_name != 'VANCOMYCIN' & is.na(grid.1) ~ 1,
      drug_name != 'VANCOMYCIN' & !is.na(grid.1) ~ 0,
      drug_name == 'VANCOMYCIN' & is.na(grid.1) ~ 1,
      drug_name == 'VANCOMYCIN' & !is.na(grid.1) & ivm*ivm.1 == 0 & po*po.1 == 0 ~ 1,
      drug_name == 'VANCOMYCIN' & !is.na(grid.1) ~ 0)) %>%
    group_by(grid, drug_date, drug_name, ivm, po) %>% 
    summarise(new_abx = min(new_abx)) %>% 
    ungroup()
  
  # abx_raw3 %>%
  #   count(new_abx)
  # abx_raw3 %>%
  #   filter(drug_name == 'VANCOMYCIN') %>%
  #   count(new_abx)
  
  
  ## NEW QAD within 2 days of blood culture day 
  first_qad <- sqldf::sqldf('SELECT * 
                           FROM abx_raw3 as t1
                           INNER JOIN blood_raw2 as t2
                           ON t1.grid = t2.grid AND drug_date <= dc_date AND drug_date BETWEEN blood_date-2 AND blood_date+2') %>% 
    setNames(make.unique(names(.))) %>%  # make unique column names
    as_tibble() %>% 
    select(-grid.1) %>% 
    filter(new_abx == 1) %>% 
    distinct(grid, adm_date, dc_date, icu_start_date, blood_date, drug_date) %>% 
    arrange(grid, blood_date, drug_date) %>% 
    rename(first_qad_date = drug_date) #39,191  -> xtan: 13,336
  # first_qad %>%
  #   distinct(grid, first_qad_date) #33,476  -> xtan: 11,151
  
  
  ## check 4 consecutive days after first QAD date 
  qads <- sqldf::sqldf('SELECT * 
             FROM first_qad as t1
             INNER JOIN abx_raw3 as t2
             ON t1.grid = t2.grid AND drug_date BETWEEN first_qad_date AND first_qad_date + 4') %>% 
    setNames(make.unique(names(.))) %>%  # make unique column names
    as_tibble() %>%
    select(-grid.1, -po) %>%
    arrange(grid, adm_date, dc_date, icu_start_date, blood_date, first_qad_date, drug_date, drug_name) %>% 
    group_by(grid, adm_date, dc_date, icu_start_date, blood_date, first_qad_date, drug_name) %>% 
    # new abx initiated in current QAD sequence
    mutate(new_abx_this = if_else(cumsum(new_abx) > 0, 1, 0)) %>% 
    ungroup() 
  # qads %>% 
  #   count(new_abx_this, new_abx)
  # qads %>% 
  #   filter(drug_date > dc_date)
  # qads %>% 
  #   filter(first_qad_date > dc_date)
  
  
  # calculate # of calendar days, # of drug days, max gap for each QAD sequence
  qads_grp <- qads %>%  
    filter(new_abx_this == 1) %>%  # remove abx not initiated in current QAD sequence since they are not counted as QADs.
    distinct(grid, adm_date, dc_date, icu_start_date, blood_date, first_qad_date, drug_date) %>%
    group_by(grid, adm_date, dc_date, icu_start_date, blood_date, first_qad_date) %>%
    summarise(max_day = as.numeric(max(drug_date - first_qad_date) + 1),
              n_day = n_distinct(drug_date),
              max_gap = if_else(max_day == 1, 0, as.numeric(max(drug_date - lag(drug_date), na.rm = T)))
    ) %>% 
    ungroup()
  
  # qads_grp %>% 
  #   count(max_day, n_day, max_gap)
  
  # filter by >= 4 QADs, worst case, 3 same antibiotics in 5 with one day gap inbetween
  qads_ge4d <- qads_grp %>% 
    filter(max_day >= 4, n_day >= 3,  max_gap <= 2) 
  # qads_ge4d %>% 
  #   count(max_day, n_day, max_gap)
  
  # for <= 4 calendar days, consider death
  qads_l4d <- qads_grp %>%  
    filter(max_day <= 3, max_gap <= 1) %>% 
    mutate(last_qad = first_qad_date + max_day - 1)
  # qads_l4d %>% 
  #   count(max_day, n_day, max_gap)
  
  #> check death  data --------------------------------
  # static_raw %>% 
  #   select(grid, dod, death_flag) %>% 
  #   Hmisc::describe()
  
  static_raw1 <- static_raw %>% 
    select(grid, dod, death_flag) %>% 
    left_join(changed_grid, by = c("grid" = "old_grid")) %>% 
    mutate(
      dod = if_else(!is.na(updated_grid),
                    as.Date(dod, format="%m/%d/%Y") - old_dob + updated_dob,
                    as.Date(dod, format="%m/%d/%Y")),
      grid = if_else(!is.na(updated_grid), updated_grid, grid)) %>% 
    distinct(grid, dod, death_flag) %>% 
    arrange(grid, death_flag, dod)
  
  ## check data
  # static_raw1 %>% 
  #   group_by(grid) %>% 
  #   filter(n() > 1)
  # static_raw1 %>% 
  #   filter(!is.na(dod), death_flag != 1)
  # static_raw1 %>% 
  #   filter(is.na(dod), death_flag != 0)
  # static_raw1 %>% 
  #   group_by(grid) %>% 
  #   filter(n_distinct(dod, na.rm = T) > 1)
  
  ## check QAD for these two GRIDs
  # qads %>% 
  #   filter(grid %in% c("R213850517", "R297625095"))
  ## All we want is death date now
  death_date <- static_raw1 %>% 
    filter(!is.na(dod)) %>%  
    group_by(grid) %>% 
    summarise(dod = min(dod)) %>% 
    ungroup()
  
  
  qads_l4d_death <- sqldf::sqldf('SELECT * 
             FROM qads_l4d as t1
             INNER JOIN death_date as t2
             ON t1.grid = t2.grid AND dod BETWEEN last_qad AND last_qad + 1') %>% 
    setNames(make.unique(names(.))) %>%  # make unique column names
    as_tibble() %>% 
    select(-last_qad, -grid.1)
  
  
  # filter by having at least one new IV/IM within 2 days window of blood date.
  qads_ivm <- qads %>%  
    filter(new_abx_this == 1) %>%
    filter(new_abx == 1, ivm == 1, abs(drug_date - blood_date) <= 2) %>%
    distinct(grid, adm_date, dc_date, icu_start_date, blood_date, first_qad_date)
  
  # combining the two criteria above
  first_qad1 <- qads_ge4d %>% 
    bind_rows(qads_l4d_death) %>% 
    inner_join(qads_ivm, by = c("grid", "adm_date", "dc_date", "icu_start_date", "blood_date",  "first_qad_date"))
  
  # first_qad1 %>% 
  #   distinct(grid, first_qad_date)
  # first_qad1 %>% 
  #   distinct(grid, adm_id, blood_date)
  # first_qad1 %>% 
  #   distinct(grid, adm_id)
  
  #. Now only need to check those QAD sequence with gap = 2 days ----------------
  # first_qad1 %>% 
  #   count(max_day, n_day, max_gap)
  # 
  # first_qad1 %>% 
  #   filter(max_gap == 2, n_day == 3, max_day == 4)
  ## not 4 QAD, new drug only after gap
  # qads %>% 
  #   filter(new_abx_this == 1) %>% 
  #   filter(grid == 'R200459372', first_qad_date == ymd(20140403))
  ## 4 QAD, only 1 drug
  # qads %>% 
  #   filter(new_abx_this == 1) %>% 
  #   filter(grid == 'R200617056', first_qad_date == ymd(20091024))
  # 
  # first_qad1 %>% 
  #   filter(max_gap == 2, n_day == 3, max_day == 5)
  ## 4 QAD, only one drug
  # qads %>% 
  #   filter(new_abx_this == 1) %>% 
  #   filter(grid == 'R200388325', first_qad_date == ymd(20081113))
  ## not 4 QAD, new drug only after gap
  # qads %>% 
  #   filter(new_abx_this == 1) %>% 
  #   filter(grid == 'R201758046', first_qad_date == ymd(20100207)) 
  ## 4 QAD, same drug after gap
  # qads %>% 
  #   filter(new_abx_this == 1) %>% 
  #   filter(grid == 'R204357412', first_qad_date == ymd(20140319))
  # first_qad1 %>% 
  #   filter(max_gap == 2, n_day == 4, max_day == 5)
  ## 4 QAD, same drug after gap
  # qads %>% 
  #   filter(new_abx_this == 1) %>% 
  #   filter(grid == 'R200077040', first_qad_date == ymd(20140511))
  ## not 4 QAD, new drug only after gap
  # qads %>% 
  #   filter(new_abx_this == 1) %>% 
  #   filter(grid == 'R200124824', first_qad_date == ymd(20080418))
  ## ABX sequences with 2-day gaps
  qads_gap2 <- qads %>% 
    filter(new_abx_this == 1) %>% 
    semi_join(
      first_qad1 %>% 
        filter(max_gap == 2),
      by = c("grid", "adm_date", "dc_date", "icu_start_date", "blood_date", "first_qad_date"))
  # find the day after gap, check whether it's only new abx on this day
  # for a sequence, if any of the day after gap had only new abx, then it's not a 4 day QAD sequence
  qads_gap_rm <- qads_gap2 %>% 
    select(grid:drug_date) %>% 
    distinct() %>% 
    group_by(grid, adm_date, dc_date, icu_start_date, blood_date, first_qad_date) %>% 
    mutate(gap = drug_date - lag(drug_date)) %>% 
    ungroup() %>% 
    filter(gap == 2) %>% 
    left_join(select(qads_gap2, -ivm, -new_abx_this), 
              by = c("grid", "adm_date", "dc_date", "icu_start_date", "blood_date", "first_qad_date", 'drug_date')) %>% 
    group_by(grid, adm_date, dc_date, icu_start_date, blood_date, first_qad_date, drug_date) %>%
    summarise(all_new = if_else(sum(new_abx) == n(), 1, 0)) %>% 
    ungroup() %>% 
    group_by(grid, adm_date, dc_date, icu_start_date, blood_date, first_qad_date) %>% 
    summarise(any_gap_all_new = sum(all_new)) %>% 
    ungroup() %>% 
    filter(any_gap_all_new == 1) # remove 581 sequences
  
  # remove those unqualified sequences
  rhee_infection <- first_qad1 %>% 
    anti_join(qads_gap_rm,
              by = c("grid", "adm_date", "dc_date", "icu_start_date", "blood_date", "first_qad_date"))  %>% 
    distinct(grid, adm_date, dc_date, icu_start_date, blood_date)
  # rhee_infection %>% distinct(grid, adm_date) #16975  -> xtan: 5673
  
  return(rhee_infection)
}

# write_csv(rhee_infection, "../datafile/rhee_infection_old_20200609.csv")
