# BioVU Data Analysis

The goal of the project is to study the association between mtDNA haplogroups and delirium in sepsis patients.  


## Statistical modeling (regression and mediation analysis)

### Input data
- Haplogroup: `Mito Delirium BioVU Data/Genetics/CAM_Haplogroups.xlsx`
- Date of death (DOD): `Mito Delirium BioVU Data/Demographics/Date of Death data.xlsx` (sheet `FINAL DATE OF DEATH DATA`)
- All cohort subjects: `datafile/sepsis_grids_20200106.xlsx` (sheet `All GRIDs`)
- Admission & discharge dates (include worst sofa score per admission): `datafile/sepsis_compare_20191217.csv`
- Daily CAM status: `datafile/daily_status_20190925.csv`
- Daily lab (includes sofa, rass, creatinine, platelet ...): `datafile/daily_sofa_score_20191010.csv`
- Neuro damage data (remove encounters with "bad" icd codes): 
  - ICD9: `Mito Delirium BioVU Data/Neuro Exclusions/neuro_icd9_V2.xlsx`
  - ICD10: `Mito Delirium BioVU Data/Neuro Exclusions/neuro_icd10_V2.xlsx`
- Comorbidity score: `Mito Delirium BioVU Data/Elixhauser Comorbidities/*.xlsx`
- Daily dementia: `Mito Delirium BioVU Data/Dementia/Dementia.xlsx` (sheet `Initial Dementia Code date`)
- Medications
  - `Mito Delirium BioVU Data/Data/grid_date_med1.csv`
  - `Mito Delirium BioVU Data/Data/grid_date_med2.csv`
  - `Mito Delirium BioVU Data/Data/grid_date_med3.csv`
- In ICU death manual review: `datafile/in_ICU_death_manual_review.txt`

### Output data
- Data dictionary: `data_arxiv/analysis_daily_dict.xlsx`
- Daily-level data: `data_arxiv/analysis_daily.rds`
- Encounter-level data: `grid` and `adm_date` in daily-level data uniquely determine an encounter
- Subject-level data: `grid` in daily-level data uniquely determine a subject

### Functions
- TODO
- 

