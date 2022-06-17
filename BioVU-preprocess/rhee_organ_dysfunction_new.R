library(tidyverse)
library(readxl)
library(magrittr)
library(lubridate)

FUN_rhee_organ_new = function(rhee_infection) {
  
  # rhee_infection <- read_csv("../datafile/rhee_infection_old_20200609.csv")
  changed_grid <- read_csv("../datafile/changed_grid_dob_20190924.csv")
  static_raw <- read_csv("../Mito Delirium BioVU Data/Data/Samuels_Delirium_STATIC_20180718.csv") 
  names(static_raw) <- str_to_lower(names(static_raw))
  
  #. Vasopressor initiation ------------------
  med_raw <- NULL
  for (i in 1:3) {
    med_raw %<>% 
      bind_rows(
        read_csv(
          paste0("../Mito Delirium BioVU Data/Data/grid_date_med", i, ".csv"),
          col_names = c("grid", "drug_date", "drug_name", "drug_class", 
                        "drug_route1", "drug_route2", "drug_route3"),
          skip = 1
        ) 
      ) 
  }
  # med_raw %>% count(drug_class)
  # med_raw %>% filter(drug_class == "pressor") %>% count(drug_name)
  pressor_raw <- med_raw %>% 
    filter(drug_class == "pressor", drug_name != "DOBUTAMINE") %>% # Dobutamine is not listed as vasopressor in Rhee definition
    select(-drug_class) 
  # pressor_raw %>% 
  #   select(drug_date, drug_name, starts_with('drug_route')) %>% 
  #   Hmisc::describe()
  
  
  #> convert messed-up GRIDs and dates ---
  # length(unique(pressor_raw$grid))
  # sum(unique(pressor_raw$grid) %in% unique(static_raw$grid))  # check whether any new updated GRID
  # sum(unique(pressor_raw$grid) %in% unique(changed_grid$old_grid))
  # sum(unique(pressor_raw$grid) %in% unique(changed_grid$updated_grid))
  
  
  pressor_raw1 <- pressor_raw %>% 
    left_join(changed_grid, by = c("grid" = "old_grid")) %>% 
    mutate(
      drug_date = if_else(!is.na(updated_grid),
                          as.Date(drug_date) - old_dob + updated_dob,
                          as.Date(drug_date)),
      grid = if_else(!is.na(updated_grid), updated_grid, grid)
    ) %>% 
    distinct(grid, drug_date, drug_name, drug_route1, drug_route2, drug_route3) %>% 
    arrange(grid, drug_date) # one duplicate
  # sum(pressor_raw1$grid %in% changed_grid$old_grid)
  # pressor_raw1 %>% 
  #   count(drug_route1, drug_route2, drug_route3)
  
  # keep IV vasopressors only
  # In the supp. We only included intravenous administrations of vasopressors, and excluded vasopressors that were clearly single bolus injections rather than continuous infusions.
  pressor_raw1 %<>% 
    filter(drug_route1 == "IV" | drug_route2 == "IV" | drug_route3 == "IV")
  
  
  #> find initiation of new vasopressors (no same vasopressor on the prior calendar day) ---
  pressor_new_date <- sqldf::sqldf('SELECT *
        FROM pressor_raw1 as t1
        left join pressor_raw1 as t2
        ON t1.grid = t2.grid AND t1.drug_name == t2.drug_name AND t2.drug_date == t1.drug_date - 1') %>% 
    setNames(make.unique(names(.))) %>%  # make unique column names
    as_tibble() %>% 
    filter(is.na(grid.1)) %>% 
    distinct(grid, drug_date) 
  
  pressor_new <- sqldf::sqldf('SELECT * 
              FROM rhee_infection as t1
             INNER JOIN pressor_new_date as t2 
             ON t1.grid = t2.grid AND drug_date <= dc_date AND drug_date BETWEEN blood_date-2 AND blood_date+2') %>% 
    setNames(make.unique(names(.))) %>%  # make unique column names
    as_tibble() %>% 
    select(-grid.1) %>% 
    distinct(grid, adm_date, dc_date, icu_start_date, blood_date) %>% 
    mutate(new_pressor = 1)
  
  #. Mechnical ventilation initiation -----------------
  vent_raw <- read_excel("../Mito Delirium BioVU Data/Phenotype data/ventilation_days.xlsx")
  names(vent_raw) <- str_to_lower((names(vent_raw)))
  # vent_raw %>% Hmisc::describe()
  
  #> convert messed-up GRIDs and dates ---
  # length(unique(vent_raw$grid))
  # sum(unique(vent_raw$grid) %in% static_raw$grid)
  # sum(unique(vent_raw$grid) %in% changed_grid$old_grid)
  # sum(unique(vent_raw$grid) %in% changed_grid$updated_grid)
  vent_raw1 <- vent_raw %>% 
    left_join(changed_grid, by = c("grid" = "old_grid")) %>% 
    mutate(
      code_date = if_else(!is.na(updated_grid),
                          as_date(code_date) - old_dob + updated_dob,
                          as_date(code_date)),
      grid = if_else(!is.na(updated_grid), updated_grid, grid)
    ) %>% 
    distinct(grid, code_date) %>% 
    arrange(grid, code_date)
  # sum(vent_raw1$grid %in% changed_grid$old_grid)
  
  #> find initiation of mechanical ventilation (no ventilation on the prior calendar day) ---
  vent_new_date <- sqldf::sqldf('SELECT *
        FROM vent_raw1 as t1
        left join vent_raw1 as t2
        ON t1.grid = t2.grid AND t2.code_date == t1.code_date - 1') %>% 
    setNames(make.unique(names(.))) %>%  # make unique column names
    as_tibble() %>% 
    filter(is.na(grid.1)) %>% 
    distinct(grid, code_date) 
  
  vent_new <- sqldf::sqldf('SELECT * 
              FROM rhee_infection as t1
             INNER JOIN vent_new_date as t2 
             ON t1.grid = t2.grid AND code_date <= dc_date AND code_date BETWEEN blood_date-2 AND blood_date+2') %>% 
    setNames(make.unique(names(.))) %>%  # make unique column names
    as_tibble() %>% 
    select(-grid.1) %>% 
    distinct(grid, adm_date, dc_date, icu_start_date, blood_date) %>%
    mutate(new_vent = 1)
  
  #. lactate ---------------------
  file_names <- list.files("../Mito Delirium BioVU Data/Lab values/Lactate",
                           full.names = T)
  
  #> out of range value ---
  lactate_oor <- NULL
  for (file in file_names) {
    readfile = read_excel(file, sheet = 2)
    if (dim(readfile)[1] != 0) {
      lactate_oor %<>%  
        bind_rows(readfile)
    }
  }
  
  # Hmisc::describe(lactate_oor) #238 rows
  # lactate_oor %>% 
  #   distinct(`LAC (mmol/L)`) %>% 
  #   pull(`LAC (mmol/L)`) %>% 
  #   str_view_all( "[[:digit:]]*\\.*[[:digit:]]+")
  lactate_raw <- lactate_oor %>% 
    rename(oor_value = `LAC (mmol/L)`) %>% 
    mutate(`LAC (mmol/L)` = as.numeric(str_extract(oor_value, "[[:digit:]]*\\.*[[:digit:]]+"))
    ) 
  # lactate_raw %>% 
  #   distinct(oor_value, `LAC (mmol/L)`)
  
  #> normal value ---
  for (file in file_names) {
    lactate_raw <- lactate_raw %>% 
      bind_rows(read_excel(file, sheet = 1))
  }
  # Hmisc::describe(lactate_raw)
  names(lactate_raw) <- tolower(names(lactate_raw))
  
  #> convert messed-up GRIDs and dates ---
  # length(unique(lactate_raw$grid))
  # sum(unique(lactate_raw$grid) %in% unique(static_raw$grid))  # check whether any new updated GRID
  # sum(unique(lactate_raw$grid) %in% unique(changed_grid$old_grid)) # no OLD GRID
  # sum(unique(lactate_raw$grid) %in% unique(changed_grid$updated_grid))
  # lactate_raw %>% 
  #   filter(!grid %in% static_raw$grid, !grid %in% changed_grid$updated_grid)
  
  #> find >= 2 days
  lactate_date <- lactate_raw %>% 
    filter(`lac (mmol/l)` >= 2.0) %>% 
    distinct(grid, lab_date) %>% 
    mutate(lab_date = as_date(lab_date))
  lactate_ge2 <- sqldf::sqldf('SELECT * 
              FROM rhee_infection as t1
                         INNER JOIN lactate_date as t2 
                         ON t1.grid = t2.grid AND lab_date <= dc_date AND lab_date BETWEEN blood_date-2 AND blood_date+2') %>% 
    setNames(make.unique(names(.))) %>%  # make unique column names
    as_tibble() %>% 
    select(-grid.1) %>% 
    distinct(grid, adm_date, dc_date, icu_start_date, blood_date) %>%
    mutate(lactate_ge2 = 1)
  
  #. renal dysfunction ------------------------
  #> creatinine doubling ----
  file_names <- list.files("../Mito Delirium BioVU Data/Lab values/Creatinine",
                           full.names = T)
  
  #> out of range value ---
  creatinine_oor <- NULL
  for (file in file_names) {
    creatinine_oor <- creatinine_oor %>% 
      bind_rows(read_excel(file, sheet = 2))
  }
  # Hmisc::describe(creatinine_oor) #650 rows
  # creatinine_oor %>% 
  #   distinct(`Creat mg/dL`) %>% 
  #   pull(`Creat mg/dL`) %>% 
  #   str_view_all( "[[:digit:]]*\\.*[[:digit:]]+")
  creatinine_raw <- creatinine_oor %>% 
    rename(oor_value = `Creat mg/dL`) %>% 
    mutate(`Creat mg/dL` = case_when(
      str_detect(oor_value, ",") ~ as.numeric(str_replace(oor_value, ",", ".")),
      str_detect(oor_value, "-") ~ str_extract_all(oor_value, "[[:digit:]]*\\.*[[:digit:]]+") %>% sapply(function(x) mean(as.numeric(x))),
      T ~ as.numeric(str_extract(oor_value, "[[:digit:]]*\\.*[[:digit:]]+"))
    )) 
  # creatinine_raw %>% 
  #   distinct(oor_value, `Creat mg/dL`)
  
  #> normal value ---
  for (file in file_names) {
    creatinine_raw <- creatinine_raw %>% 
      bind_rows(read_excel(file, sheet = 1))
  }
  # Hmisc::describe(creatinine_raw)
  names(creatinine_raw) <- tolower(names(creatinine_raw))
  
  #> convert messed-up GRIDs and dates ---
  # length(unique(creatinine_raw$grid))
  # sum(unique(creatinine_raw$grid) %in% unique(static_raw$grid))  # check whether any new updated GRID
  # sum(unique(creatinine_raw$grid) %in% unique(changed_grid$old_grid)) # no OLD GRID
  # sum(unique(creatinine_raw$grid) %in% unique(changed_grid$updated_grid))
  # creatinine_raw %>% 
  #   filter(!grid %in% static_raw$grid, !grid %in% changed_grid$updated_grid)
  creatinine_raw %<>% 
    filter(!is.na(`creat mg/dl`)) %>% 
    mutate(lab_date = as_date(lab_date))
  # ggplot(creatinine_raw) +
  #   geom_histogram(aes(x = `creat mg/dl`)) +
  #   scale_y_log10()
  # ggplot(creatinine_raw) +
  #   geom_boxplot(aes(x ="", y = `creat mg/dl`))
  # summary(creatinine_raw$`creat mg/dl`) 
  # creatinine_raw %>% 
  #   filter(`creat mg/dl` > 50)
  
  
  
  
  #> calculate EGFR---------------
  ## get sex, race from static_raw and convert messed-up GRIDs
  # static_raw %>% 
  #   select(sex, race) %>% 
  #   Hmisc::describe()
  demo_raw <- static_raw %>% 
    select(grid, sex, race, dob) %>%
    left_join(changed_grid, by = c("grid" = "old_grid")) %>% 
    mutate(
      grid = if_else(!is.na(updated_grid), updated_grid, grid),
      dob = if_else(!is.na(updated_dob), updated_dob, mdy(dob))) %>% 
    distinct(grid, sex, race, dob) %>% 
    arrange(grid)
  # sum(demo_raw$grid %in% changed_grid$old_grid)
  # demo_raw %>% 
  #   Hmisc::describe()
  # demo_raw %>% 
  #   group_by(grid) %>% 
  #   filter(n_distinct(race) > 1) %>% 
  #   print(n = 100)
  # demo_raw %>% 
  #   group_by(grid) %>% 
  #   filter(n_distinct(sex) > 1) %>% 
  #   print(n = 100)
  # demo_raw %>% 
  #   group_by(grid) %>% 
  #   filter(n_distinct(dob) > 1) %>% 
  #   print(n = 100)
  demo <- demo_raw %>% 
    group_by(grid) %>% 
    mutate(
      race_b = case_when(
        sum(str_detect(race, "B")) > 0 ~ "Black",
        sum(race == "U") != n() ~ "White or other")) %>% 
    ungroup() %>% 
    distinct(grid, sex, race_b, dob) %>% 
    mutate(
      sex = case_when(
        sex != "U" ~ sex))
  # demo %>% Hmisc::describe()
  
  
  ## Calculate age at a given reference date 
  calc_age <- function(birthDate, refDate = Sys.Date()) {
    period <- as.period(interval(birthDate, refDate),
                        unit = "year")
    period$year
  }
  
  egfr_raw <- creatinine_raw %>% 
    select(-oor_value, -lab_time) %>% 
    left_join(demo, by = "grid") %>% 
    mutate(
      age = calc_age(dob, lab_date),
      egfr = case_when(
        race_b == "Black" & sex == "F" & `creat mg/dl` <= 0.7 ~ 166*(`creat mg/dl`/0.7)^(-0.329)*(0.993)^age,
        race_b == "Black" & sex == "F" & `creat mg/dl` >  0.7 ~ 166*(`creat mg/dl`/0.7)^(-1.209)*(0.993)^age,
        race_b == "Black" & sex == "M" & `creat mg/dl` <= 0.9 ~ 163*(`creat mg/dl`/0.9)^(-0.411)*(0.993)^age,
        race_b == "Black" & sex == "M" & `creat mg/dl` >  0.9 ~ 163*(`creat mg/dl`/0.9)^(-1.209)*(0.993)^age,
        race_b == "White or other" & sex == "F" & `creat mg/dl` <= 0.7 ~ 144*(`creat mg/dl`/0.7)^(-0.329)*(0.993)^age,
        race_b == "White or other" & sex == "F" & `creat mg/dl` >  0.7 ~ 144*(`creat mg/dl`/0.7)^(-1.209)*(0.993)^age,
        race_b == "White or other" & sex == "M" & `creat mg/dl` <= 0.9 ~ 141*(`creat mg/dl`/0.9)^(-0.411)*(0.993)^age,
        race_b == "White or other" & sex == "M" & `creat mg/dl` >  0.9 ~ 141*(`creat mg/dl`/0.9)^(-1.209)*(0.993)^age
      ),
      egfr_60 = if_else(egfr > 60, 60, egfr)
    )
  # egfr_raw %>% 
  #   Hmisc::describe()
  # ggplot(egfr_raw) +
  #   geom_histogram(aes(x = egfr)) 
  # egfr_raw %>% 
  #   filter(egfr > 60) # more than half
  # egfr_raw %>% 
  #   filter(egfr == Inf)
  
  #> exclude end-stage kidney disease --------------------------
  esrd_raw <- read_tsv("../Mito Delirium BioVU Data/Phenotype data/ESRD_ICD.txt")
  names(esrd_raw) <- str_to_lower(names(esrd_raw))
  esrd_raw %<>% 
    mutate(esrd_code_date = mdy(esrd_code_date))
  # esrd_raw %>% 
  #   Hmisc::describe() # 4705 pts
  
  
  # check for messed-up GRID
  # length(unique(creatinine_raw$grid))
  # sum(unique(esrd_raw$grid) %in% unique(static_raw$grid))  # check whether any new updated GRID
  # sum(unique(esrd_raw$grid) %in% unique(changed_grid$old_grid)) # no OLD GRID
  # sum(unique(esrd_raw$grid) %in% unique(changed_grid$updated_grid))
  
  rhee_infection_ex <- sqldf::sqldf('SELECT * 
              FROM rhee_infection as t1
              INNER JOIN esrd_raw as t2 
              ON t1.grid = t2.grid AND esrd_code_date BETWEEN blood_date - 2 AND blood_date + 2') %>% # on 12/4/2019 meeting, we decided to use same time window (within 2 days of BC day) for end-stage renal disease code 
    setNames(make.unique(names(.))) %>%  # make unique column names
    as_tibble() %>% 
    distinct(grid, adm_date, dc_date, icu_start_date, blood_date)
  rhee_infection1 <- rhee_infection  %>% 
    anti_join(rhee_infection_ex)
  
  #> identify creatinine doubling/edfr decline infections --------------------  -> xtan: use icu_start_date instead
  # find baseline for each hospitaliztion, then filter to find creatinine doubling within 2 day window
  renal_dysfct <- sqldf::sqldf('SELECT * 
                               FROM rhee_infection1 as t1
                               INNER JOIN egfr_raw as t2 
                               ON t1.grid = t2.grid AND lab_date BETWEEN icu_start_date-3 AND dc_date') %>% 
    setNames(make.unique(names(.))) %>%  # make unique column names
    as_tibble() %>% 
    select(grid, adm_date, dc_date, blood_date, icu_start_date, lab_date, `creat mg/dl`, sex, race_b, age, egfr, egfr_60) %>% 
    group_by(grid, adm_date) %>% 
    mutate(creat_bl = min(`creat mg/dl`),
           egfr_bl = max(egfr_60)) %>% 
    ungroup() %>% 
    group_by(grid, adm_date, dc_date, icu_start_date, blood_date) %>% 
    filter(lab_date >= blood_date - 2, 
           lab_date <= blood_date + 2,
           `creat mg/dl` >= 2*creat_bl | egfr_60 <= 0.5*egfr_bl) %>% 
    ungroup()
  ## check range of creatinine
  # renal_dysfct %>% 
  #   select(`creat mg/dl`:egfr_bl) %>% 
  #   Hmisc::describe()
  # renal_dysfct %>% 
  #   filter(`creat mg/dl` < 1.2)
  # renal_dysfct %>% 
  #   filter(`creat mg/dl` < 1.2) %>% 
  #   distinct(grid, adm_id, adm_date, dc_date, creat_bl)
  ## check egfr
  # renal_dysfct %>% 
  #   filter(egfr > 60)
  # get distinct BC dates
  renal_dysfct_date <- renal_dysfct %>% 
    distinct(grid, adm_date, dc_date, icu_start_date, blood_date) %>% 
    mutate(renal_dysfct  = 1)
  
  # check End Stage Renal Disease Status
  # sqldf::sqldf('SELECT * 
  #               FROM renal_dysfct_date as t1
  #               INNER JOIN esrd_raw as t2 
  #               ON t1.grid = t2.grid AND adm_date > esrd_code_date') %>% 
  #   setNames(make.unique(names(.))) %>%  # make unique column names
  #   as_tibble() %>% 
  #   distinct(grid) #133 pts
  # esrd_raw %>% 
  #   filter(grid == "R200617056") %>% 
  #   arrange(esrd_code_date)
  
  #. bilirubin doubling --------------------------------------------------------
  file_names <- list.files("../Mito Delirium BioVU Data/Lab values/Bilirubin",
                           full.names = T)
  
  #> out of range value ---
  bilirubin_oor <- NULL
  for (file in file_names) {
    bilirubin_oor <- bilirubin_oor %>% 
      bind_rows(read_excel(file, sheet = "out of range"))
  }
  # Hmisc::describe(bilirubin_oor) #857 rows
  # bilirubin_oor %>% 
  #   distinct(`Tbil (mg/dL)`) %>% 
  #   pull(`Tbil (mg/dL)`) %>% 
  #   str_view_all( "[[:digit:]]*\\.*[[:digit:]]+")
  bilirubin_raw <- bilirubin_oor %>% 
    rename(oor_value = `Tbil (mg/dL)`) %>% 
    mutate(`Tbil (mg/dL)` = if_else(
      str_detect(oor_value, "-"),
      str_extract_all(oor_value, "[[:digit:]]*\\.*[[:digit:]]+") %>% sapply(function(x) mean(as.numeric(x))),
      as.numeric(str_extract(oor_value, "[[:digit:]]*\\.*[[:digit:]]+"))
    )) 
  # bilirubin_raw %>% 
  #   distinct(oor_value, `Tbil (mg/dL)`)
  
  #> normal value ---
  for (file in file_names) {
    bilirubin_raw <- bilirubin_raw %>% 
      bind_rows(read_excel(file, sheet = 1))
  }
  # Hmisc::describe(bilirubin_raw)
  names(bilirubin_raw) <- tolower(names(bilirubin_raw))
  
  #> convert messed-up GRIDs and dates ---
  # length(unique(bilirubin_raw$grid))
  # sum(unique(bilirubin_raw$grid) %in% unique(static_raw$grid))  # check whether any new updated GRID
  # sum(unique(bilirubin_raw$grid) %in% unique(changed_grid$old_grid)) # no OLD GRID
  # sum(unique(bilirubin_raw$grid) %in% unique(changed_grid$updated_grid))
  # bilirubin_raw %>% 
  #   filter(!grid %in% static_raw$grid, !grid %in% changed_grid$updated_grid)
  bilirubin_raw %<>% 
    filter(!is.na(`tbil (mg/dl)`)) %>% 
    mutate(lab_date = as_date(lab_date))
  
  # ggplot(bilirubin_raw) +
  #   geom_histogram(aes(x = `tbil (mg/dl)`)) +
  #   scale_y_log10()
  # summary(bilirubin_raw$`tbil (mg/dl)`) 
  # bilirubin_raw %>% 
  #   filter(`tbil (mg/dl)` > 50)
  # find baseline for each hospitaliztion, then filter to find bilirubin doubling within 2 day window 
  bilirubin_dbl <- sqldf::sqldf('SELECT * 
              FROM rhee_infection as t1
              INNER JOIN bilirubin_raw as t2 
              ON t1.grid = t2.grid AND lab_date BETWEEN icu_start_date-3 AND dc_date') %>% 
    setNames(make.unique(names(.))) %>%  # make unique column names
    as_tibble() %>% 
    select(grid, adm_date, dc_date, icu_start_date, blood_date, lab_date, `tbil (mg/dl)`) %>% 
    group_by(grid, adm_date) %>% 
    mutate(bil_bl = min(`tbil (mg/dl)`)) %>% 
    ungroup() %>% 
    group_by(grid, adm_date, dc_date, icu_start_date, blood_date) %>% 
    filter(lab_date >= blood_date - 2, 
           lab_date <= blood_date + 2,
           `tbil (mg/dl)` >= 2.0,
           `tbil (mg/dl)` >= 2*bil_bl) %>% 
    ungroup() 
  ## check range of bilirubin
  # bilirubin_dbl %>% 
  #   select(bil_bl, `tbil (mg/dl)`) %>% 
  #   Hmisc::describe()
  # bilirubin_dbl %>% 
  #   filter(`tbil (mg/dl)` > 50)
  # bilirubin double dates
  bilirubin_dbl_date <- bilirubin_dbl %>% 
    distinct(grid, adm_date, dc_date, icu_start_date, blood_date) %>% 
    mutate(bil_dbl  = 1)
  
  
  #. platelet decline ---------------------
  #> out of range value ---
  file_names <- list.files("../Mito Delirium BioVU Data/Lab values/Platelet",
                           pattern = "range.xlsx$",
                           full.names = T)
  platelet_oor <- NULL
  for (file in file_names) {
    platelet_oor <- platelet_oor %>% 
      bind_rows(read_excel(file))
  }
  # Hmisc::describe(platelet_oor) # 1462 rows
  # platelet_oor %>% 
  #   distinct(`Plt-Ct (thou/uL)`) %>% 
  #   pull(`Plt-Ct (thou/uL)`) %>% 
  #   str_view_all( "[[:digit:]]+") 
  platelet_raw <- platelet_oor %>% 
    rename(oor_value = `Plt-Ct (thou/uL)`) %>% 
    mutate(`Plt-Ct (thou/uL)` = if_else(
      str_detect(oor_value, "-"),
      str_extract_all(oor_value, "[[:digit:]]+") %>% sapply(function(x) mean(as.numeric(x))),
      as.numeric(str_extract(oor_value, "[[:digit:]]+"))
    )) 
  # platelet_raw %>% 
  #   distinct(oor_value, `Plt-Ct (thou/uL)`)
  
  #> normal value ---
  file_names <- list.files("../Mito Delirium BioVU Data/Lab values/Platelet",
                           pattern = "labs.xlsx$",
                           full.names = T)
  for (file in file_names) {
    platelet_raw <- platelet_raw %>% 
      bind_rows(read_excel(file))
  }
  # Hmisc::describe(platelet_raw)
  names(platelet_raw) <- tolower(names(platelet_raw))
  
  #> convert messed-up GRIDs and dates ---
  # length(unique(platelet_raw$grid))
  # sum(unique(platelet_raw$grid) %in% unique(static_raw$grid))  # check whether any new updated GRID
  # sum(unique(platelet_raw$grid) %in% unique(changed_grid$old_grid)) # No OLD GRIDS
  # sum(unique(platelet_raw$grid) %in% unique(changed_grid$updated_grid))
  # platelet_raw %>%
  #   filter(!grid %in% static_raw$grid, !grid %in% changed_grid$updated_grid)
  platelet_raw %<>%
    filter(!is.na(`plt-ct (thou/ul)`)) %>%
    mutate(lab_date = as_date(lab_date))
  
  # ggplot(platelet_raw) +
  #   geom_histogram(aes(x = `plt-ct (thou/ul)`)) +
  #   scale_y_log10()
  # ggplot(creatinine_raw) +
  #   geom_boxplot(aes(x ="", y = `creat mg/dl`))
  # summary(platelet_raw$`plt-ct (thou/ul)`) 
  # platelet_raw %>% 
  #   filter(`plt-ct (thou/ul)` > 5000)
  
  # find baseline for each hospitaliztion, then filter to find platelet doubling within 2 day window
  platelet_dcl <- sqldf::sqldf('SELECT *
                              FROM rhee_infection as t1
                              INNER JOIN platelet_raw as t2
                              ON t1.grid = t2.grid AND lab_date BETWEEN icu_start_date-3 AND dc_date') %>%
    setNames(make.unique(names(.))) %>%  # make unique column names
    as_tibble() %>%
    select(grid, adm_date, dc_date, icu_start_date, blood_date, lab_date, `plt-ct (thou/ul)`) %>%
    group_by(grid, adm_date) %>%
    mutate(plt_bl = max(`plt-ct (thou/ul)`)) %>%
    ungroup() %>%
    group_by(grid, adm_date, dc_date, icu_start_date, blood_date) %>%
    filter(lab_date >= blood_date - 2,
           lab_date <= blood_date + 2,
           plt_bl > 100, # Kaveh think Rhee paper has a typo of platelets, should be 100x10^3 cells/uL
           `plt-ct (thou/ul)` < 100,
           `plt-ct (thou/ul)` <= 0.5*plt_bl) %>%
    ungroup()
  ## check range of creatinine
  # platelet_dcl %>% 
  #   select(plt_bl, `plt-ct (thou/ul)`) %>% 
  #   Hmisc::describe()
  # platelet_dcl %>% 
  #   filter(plt_bl > 1000)
  # platelet_dcl %>% 
  #   filter(plt_bl > 1000) %>% 
  #   distinct(grid, adm_id, adm_date, dc_date, plt_bl)
  # get distinct BC dates
  platelet_dcl_date <- platelet_dcl %>%
    distinct(grid, adm_date, dc_date, icu_start_date, blood_date) %>%
    mutate(plt_dcl  = 1)
  
  #. Put everything together for Rhee definition ---------------------------------
  ## keep the first blood date with the worst organ dysfunction value for each encounter
  rhee_sepsis <- rhee_infection %>% 
    left_join(pressor_new) %>% 
    left_join(vent_new) %>% 
    left_join(renal_dysfct_date) %>% 
    left_join(bilirubin_dbl_date) %>% 
    left_join(platelet_dcl_date) %>% 
    left_join(lactate_ge2) %>% 
    mutate(rhee = if_else(is.na(new_pressor | new_vent | renal_dysfct | bil_dbl | plt_dcl | lactate_ge2), 0, 1)) %>% 
    group_by(grid, adm_date, dc_date, icu_start_date) %>% 
    slice(which.max(rhee)) %>% 
    ungroup()
  
  # rhee_sepsis %>% 
  #   count(rhee)
  
  return(rhee_sepsis)
}

# write_csv(rhee_sepsis, "../datafile/sepsis_rhee_new_20200609.csv")
