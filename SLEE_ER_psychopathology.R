### Load libraries ####
library(dplyr) #for dataframe management
library(purrr) #for dataframe management
library(readr) #to read csvs
library(FSA) #for function Dunn tests ie dunnTest()
library(lme4) #for linear regression
library(lmerTest) #for linear regression p values
library(rsq) #for variance explained for mixed effect linear regression models
library(nortest) #for function Anderson-Darling tests ie ad.test
library(psych) #for correlation matrices ie corr.test()
library(lavaan) #for SEM
library(misty) #for grand mean centering variables ie center()
library(bruceR) #for conditional process modeling

### Prevent use of scientific notation ####
options(scipen=999)

## PREP DATA ####

### Read in raw data from release 5.1 ####
#### Gender data ####
gish_y_gi <- read_csv("stress_er_psychopathology/gish_y_gi.csv")

#### DERS-P for emotion regulation ####
# Note: warning will flag 18 problems total for row 4938 where the following
# columns are NA: "mh_p_ders__catast_003", "mh_p_ders__catast_004",
# "mh_p_ders__catast_005", "mh_p_ders__catast_006", "mh_p_ders__catast_007"
# "mh_p_ders__catast_008", "mh_p_ders__catast_009", "mh_p_ders__catast_010",
# "mh_p_ders__catast_011", "mh_p_ders__catast_012", "mh_p_ders__distract_002"
# "mh_p_ders__distract_003", "mh_p_ders__distract_004", "mh_p_ders__negscnd_003",
# "mh_p_ders__negscnd_004", "mh_p_ders__negscnd_005","mh_p_ders__negscnd_006",
# "mh_p_ders__negscnd_007" 
mh_p_ders <- read_tsv('/shared/abcd-6-0-release-tabulated-prod/mh_p_ders.tsv')
# problems(mh_p_ders)
# colnames(mh_p_ders[problems(mh_p_ders)$col])

#### BPM for youth-report psychopathology symptoms ####
# Note: warning will flag 6 problems total from rows 29049, 47184, 49676, 50048,
# 54499, and 75969 from 5 total columns, none of which are bpm internalizing or
# externalizing scores
mh_y_bpm <- read_tsv('/shared/abcd-6-0-release-tabulated-prod/mh_y_bpm.tsv')
# problems(mh_y_bpm)
# unique(problems(mh_y_bpm)$row)
# unique(problems(mh_y_bpm)$col)
# colnames(mh_y_bpm[,unique(problems(mh_y_bpm)$col)])
# problems(mh_y_bpm)[which(problems(mh_y_bpm)$col==which(colnames(mh_y_bpm)=="mh_y_bpm__int_sum")),]
# problems(mh_y_bpm)[which(problems(mh_y_bpm)$col==which(colnames(mh_y_bpm)=="mh_y_bpm__ext_sum")),]

#### Longitudinal tracking data ####
# Note: warning will flag 38 problems total from 28 rows all from column 6 ie
# ab_g_dyn__visit_days so okay to ignore this warning
ab_g_dyn <- read_tsv('/shared/abcd-6-0-release-tabulated-prod/ab_g_dyn.tsv')
# problems(ab_g_dyn)
# unique(problems(ab_g_dyn)$row)
# unique(problems(ab_g_dyn)$col)
# colnames(ab_g_dyn)[6]

#### LES (youth-reported) ####
#### LES for exposure to negative life events
# Note: warning will flag 30 problems total from 7 rows and 19 columns, none of
# which are mh_y_ple__exp__bad_count__v01 ie year 3 sum of bad events
mh_y_ple <- read_tsv('/shared/abcd-6-0-release-tabulated-prod/mh_y_ple.tsv')
# problems(mh_y_ple)
# unique(problems(mh_y_ple)$row)
# unique(problems(mh_y_ple)$col)
# colnames(mh_y_ple)[unique(problems(mh_y_ple)$col)]

### Prepare longitudinal tracking data for analysis ####
trackdata <- ab_g_dyn %>%
  # rename age and site columns
  rename(age = ab_g_dyn__visit_age, site = ab_g_dyn__design_site) %>%
  # select only relevant columns
  select(participant_id,session_id,age,site)

# see all unique values (can visually check for NA or errors)
map(trackdata,unique)

### Prepare gender data for analysis ####

#### Identify gender groups ####
genderdata <- gish_y_gi %>%
  # keep only data from year 3 and 4 follow-up visits
  # Before this step, n should be 49083. After this step, n should be 15064.
  filter(eventname=="3_year_follow_up_y_arm_1"|
           eventname=="4_year_follow_up_y_arm_1") %>%
  # convert numeric sex values to human-readable character strings
  mutate(sex_details = case_when(
    kbi_sex_assigned_at_birth==1 ~ "male",
    kbi_sex_assigned_at_birth==2 ~ "female",
    kbi_sex_assigned_at_birth==999 ~ "dont_know",
    kbi_sex_assigned_at_birth==777 ~ "refuse"
  )) %>%
  # set "don't know" or "refuse" to be NA for sex
  mutate(sex = case_when(
    sex_details=="male" ~ "male",
    sex_details=="female" ~ "female",
    sex_details=="dont_know" ~ NA_character_,
    sex_details=="refuse" ~ NA_character_
  )) %>%
  # make "male" reference level for sex
  mutate(sex = relevel(as.factor(sex), ref="male")) %>%
  # convert numeric gender values to human-readable character strings
  mutate(gender = case_when(
    kbi_gender==1 ~ "boy",
    kbi_gender==2 ~ "girl",
    kbi_gender==3 ~ "nb",
    kbi_gender==999 ~ "dont_understand",
    kbi_gender==777 ~ "refuse"
  )) %>%
  # convert numeric values for trans identity to human-readable character strings
  mutate(trans = case_when(
    kbi_y_trans_id==1 ~ "yes",
    kbi_y_trans_id==2 ~ "maybe",
    kbi_y_trans_id==3 ~ "no",
    kbi_y_trans_id==4 ~ "dont_understand",
    kbi_y_trans_id==777 ~ "refuse"
  )) %>%
  # combine gender and trans identity information to make five gender groups
  mutate(gender_details = case_when(
    kbi_gender==777 ~ "refuse", #refuse to provide gender
    kbi_y_trans_id==777 ~ "refuse", #refuse to say whether trans
    kbi_gender==999 ~ "dont_understand", #don't understand gender item
    kbi_y_trans_id==4 ~ "dont_understand", #don't understand trans item
    kbi_gender==1 & kbi_y_trans_id==1 ~ "trans_boy", #boy and yes trans
    kbi_gender==1 & kbi_y_trans_id==2 ~ "trans_boy", #boy and maybe trans
    kbi_gender==1 & kbi_y_trans_id==3 ~ "cis_boy", #boy and not trans
    kbi_gender==2 & kbi_y_trans_id==1 ~ "trans_girl", #girl and yes trans
    kbi_gender==2 & kbi_y_trans_id==2 ~ "trans_girl", #girl and maybe trans
    kbi_gender==2 & kbi_y_trans_id==3 ~ "cis_girl", #girl and not trans
    kbi_gender==3 & kbi_y_trans_id==1 ~ "nb", #nonbinary/other gender and yes trans
    kbi_gender==3 & kbi_y_trans_id==2 ~ "nb", #nonbinary/other gender and maybe trans
    kbi_gender==3 & kbi_y_trans_id==3 ~ "nb" #nonbinary/other gender and not trans
  )) %>%
  # combine all not cis groups due to sample size to make one gender diverse group (gd)
  mutate(genderid = case_when(
    gender_details=="trans_boy" ~ "gd",
    gender_details=="trans_girl" ~ "gd",
    gender_details=="nb" ~ "gd",
    gender_details=="cis_boy" ~ "cis_boy",
    gender_details=="cis_girl" ~ "cis_girl",
    gender_details=="refuse" ~ "refuse",
    gender_details=="dont_understand" ~ "dont_understand"
  )) %>%
  
  ### to get number of gd youth in year 3 vs year 4 and to get number
  ### of subjects who said refuse or don't know to gender question at year 4
  ### for methods, uncomment and run just the line below
  # group_by(eventname,genderid) %>% count() %>% print()
  ### to get number of subjects who said refuse or don't know to sex question at
  ### year 4, uncomment and run just the line below
  # group_by(eventname,sex_details) %>% count() %>% print()
  
  # remove subjects who refused to answer and/or did not understand gender
  # or trans questions. Before this step, n should be 15064. After this step,
# n should be 14495.
filter(genderid!="refuse",
       genderid!="dont_understand") %>%
  # make genderid a factor (automatically uses cisboy as reference)
  mutate(genderid = as.factor(genderid)) %>%
  # make columns for participant id and session so consistent with other data
  mutate(session_id = case_when(eventname=="3_year_follow_up_y_arm_1" ~ "ses-03A",
                                eventname=="4_year_follow_up_y_arm_1" ~ "ses-04A"),
         participant_id = str_replace(src_subject_id, "NDAR_INV", "sub-")) %>%
  # keep only columns relevant to analysis
  select(participant_id,session_id,
         sex,sex_details,
         genderid,gender_details
  ) 

# see all unique values (can visually check for NA or errors)
map(genderdata,unique)

#### Count number of subjects per gender group per data collection year ####
genderdata %>% 
  group_by(session_id) %>% 
  count(genderid)

### Prepare BPM data for analysis ####
bpmdata <- mh_y_bpm %>%
  # keep only year 4 data. Before this step, n = 94624 and after this step, n = 9664
  filter(session_id=="ses-04A") %>%
  # select only columns relevant to analysis
  select(participant_id,session_id,
         mh_y_bpm__int_sum,
         mh_y_bpm__ext_sum
  ) %>%
  # rename subscale columns to be more human-readable and shorter
  rename(bpm_ext = mh_y_bpm__ext_sum,
         bpm_int = mh_y_bpm__int_sum
  ) %>%
  # make internalizing and externalizing scores numeric
  mutate(bpm_ext = as.numeric(na_if(bpm_ext, "n/a")),
         bpm_int = as.numeric(na_if(bpm_int, "n/a"))) %>%
  # remove subjects with NA for bpm internalizing or externalizing. Before this
  # step, n = 9664 and after this step n = 9296
  filter(!is.na(bpm_int),
         !is.na(bpm_ext))

# see all unique values (can visually check for NA or errors)
map(bpmdata,unique)

### Prepare DERS-P data for analysis ####

#### Create cumulative score ####
dersdata <- mh_p_ders %>%
  # keep only year 3 data. Before this step, n = 33218 and after this step, n = 10140
  filter(session_id=="ses-03A") %>%
  # keep only relevant columns
  select(-c("mh_p_ders_dtt","mh_p_ders_age","mh_p_ders_lang",
            "mh_p_ders__attun_mean","mh_p_ders__attun_nm",
            "mh_p_ders__catast_mean","mh_p_ders__catast_nm",
            "mh_p_ders__distract_mean","mh_p_ders__distract_nm",
            "mh_p_ders__negscnd_mean","mh_p_ders__negscnd_nm")) %>%
  # remove subjects who refused to answer one or more items. Before this step,
  # n should be 10140 After this step, n should be 9821
  filter(!if_any(everything(), ~ . == 777)) %>%
  # remove subjects who answered don't know to one or more items. Before this step,
  # n should be 9821 After this step, n should be 9821
  filter(!if_any(everything(), ~ . == 999)) %>%
  # remove subjects with na in any row. Before this step, n should be 9821
  # After this step, n should still be 9821
  drop_na() %>%
  # add column to reverse score "my child is clear about their feelings"
  mutate(rev_mh_p_ders__attun_001 = 
           case_when(mh_p_ders__attun_001 == 1 ~ 5,
                     mh_p_ders__attun_001 == 2 ~ 4,
                     mh_p_ders__attun_001 == 3 ~ 3,
                     mh_p_ders__attun_001 == 4 ~ 2,
                     mh_p_ders__attun_001 == 5 ~ 1)) %>%
  # add column to reverse score "my child pays attention to how they feel"
  mutate(rev_mh_p_ders__attun_002 = 
           case_when(mh_p_ders__attun_002 == 1 ~ 5,
                     mh_p_ders__attun_002 == 2 ~ 4,
                     mh_p_ders__attun_002 == 3 ~ 3,
                     mh_p_ders__attun_002 == 4 ~ 2,
                     mh_p_ders__attun_002 == 5 ~ 1)) %>%
  # add column to reverse score "my child is attentive to their feelings"
  mutate(rev_mh_p_ders__attun_003 = 
           case_when(mh_p_ders__attun_003 == 1 ~ 5,
                     mh_p_ders__attun_003 == 2 ~ 4,
                     mh_p_ders__attun_003 == 3 ~ 3,
                     mh_p_ders__attun_003 == 4 ~ 2,
                     mh_p_ders__attun_003 == 5 ~ 1)) %>%
  # add column to reverse score "my child knows exactly how they are feeling"
  mutate(rev_mh_p_ders__attun_004 = 
           case_when(mh_p_ders__attun_004 == 1 ~ 5,
                     mh_p_ders__attun_004 == 2 ~ 4,
                     mh_p_ders__attun_004 == 3 ~ 3,
                     mh_p_ders__attun_004 == 4 ~ 2,
                     mh_p_ders__attun_004 == 5 ~ 1)) %>%
  # add column to reverse score "my child cares about what they are feeling"
  mutate(rev_mh_p_ders__attun_005 = 
           case_when(mh_p_ders__attun_005 == 1 ~ 5,
                     mh_p_ders__attun_005 == 2 ~ 4,
                     mh_p_ders__attun_005 == 3 ~ 3,
                     mh_p_ders__attun_005 == 4 ~ 2,
                     mh_p_ders__attun_005 == 5 ~ 1)) %>%
  # add column to reverse score "when my child is upset, they acknowledge their emotions"
  mutate(rev_mh_p_ders__attun_006 = 
           case_when(mh_p_ders__attun_006 == 1 ~ 5,
                     mh_p_ders__attun_006 == 2 ~ 4,
                     mh_p_ders__attun_006 == 3 ~ 3,
                     mh_p_ders__attun_006 == 4 ~ 2,
                     mh_p_ders__attun_006 == 5 ~ 1)) %>%
  # add column to reverse score "when my child is upset, they know that
  # they can find a way to eventually feel better"
  mutate(rev_mh_p_ders__catast_006 = 
           case_when(mh_p_ders__catast_006 == 1 ~ 5,
                     mh_p_ders__catast_006 == 2 ~ 4,
                     mh_p_ders__catast_006 == 3 ~ 3,
                     mh_p_ders__catast_006 == 4 ~ 2,
                     mh_p_ders__catast_006 == 5 ~ 1)) %>%
  # add column to reverse score "when my child is upset, they feel like
  # they can remain in control of their behaviors"
  mutate(rev_mh_p_ders__catast_007 = 
           case_when(mh_p_ders__catast_007 == 1 ~ 5,
                     mh_p_ders__catast_007 == 2 ~ 4,
                     mh_p_ders__catast_007 == 3 ~ 3,
                     mh_p_ders__catast_007 == 4 ~ 2,
                     mh_p_ders__catast_007 == 5 ~ 1))

ders_cols_to_sum <- dersdata %>%
  select(-all_of(c("participant_id","session_id","mh_p_ders__attun_001",
                   "mh_p_ders__attun_002","mh_p_ders__attun_003","mh_p_ders__attun_004",
                   "mh_p_ders__attun_005","mh_p_ders__attun_006","mh_p_ders__catast_006",
                   "mh_p_ders__catast_007"))) %>% colnames()
print(ders_cols_to_sum)
# The following 29 columns should be part of the list to be summed:
# "rev_mh_p_ders__attun_001","rev_mh_p_ders__attun_002","rev_mh_p_ders__attun_003",
# "rev_mh_p_ders__attun_004","rev_mh_p_ders__attun_005","rev_mh_p_ders__attun_006",
# "mh_p_ders__catast_001","mh_p_ders__catast_002","mh_p_ders__catast_003",
# "mh_p_ders__catast_004","mh_p_ders__catast_005","rev_mh_p_ders__catast_006,
# "rev_mh_p_ders__catast_007","mh_p_ders__catast_008","mh_p_ders__catast_009",
# "mh_p_ders__catast_010","mh_p_ders__catast_011","mh_p_ders__catast_012",
# "mh_p_ders__distract_001","mh_p_ders__distract_002","mh_p_ders__distract_003",
# "mh_p_ders__distract_004","mh_p_ders__negscnd_001","mh_p_ders__negscnd_002",
# "mh_p_ders__negscnd_003","mh_p_ders__negscnd_004","mh_p_ders__negscnd_005",
# "mh_p_ders__negscnd_006","mh_p_ders__negscnd_007"

dersdata <- dersdata %>%
  # add column to sum across all items (using using reverse-scored versions of
  # eight items above) and make one cumulative score
  mutate(ders_total = rowSums(
    across(all_of(ders_cols_to_sum)))) %>%
  # select only columns relevant to analysis
  select(participant_id,session_id,ders_total)

# see all unique values (can visually check for NA or errors)
map(dersdata,unique)

### Prepare LES data for analysis ####
#### Identify and prepare relevant columns ####
ledata <- mh_y_ple  %>%
  # keep only year 2 data. Before this step, n = 56129 and after this step, n = 10946
  filter(session_id=="ses-02A") %>%
  # code specific items
  mutate(family_died = 
           case_when(is.na(mh_y_ple_001) ~ NA, #NA for experienced or not
                     mh_y_ple_001==777 ~ NA, # decline to answer
                     mh_y_ple_001==999 ~ NA, # don't know
                     mh_y_ple_001==0 ~ 0, # did not experience
                     mh_y_ple_001==1 & mh_y_ple__exp_001==2 ~ 1, # experienced and mostly bad
                     mh_y_ple_001==1 & mh_y_ple__exp_001==1 ~ 0, # experienced and mostly good
                     mh_y_ple_001==1 & mh_y_ple__exp_001==999 ~ 0, # experienced and don't know
                     mh_y_ple_001==1 & mh_y_ple__exp_001==444 ~ NA, # experienced and not applicable
                     mh_y_ple_001==1 & mh_y_ple__exp_001=="n/a" ~ NA # experienced and NA
           )) %>%
  mutate(family_injured = 
           case_when(is.na(mh_y_ple_002) ~ NA, #NA for experienced or not
                     mh_y_ple_002==777 ~ NA, # decline to answer
                     mh_y_ple_002==999 ~ NA, # don't know
                     mh_y_ple_002==0 ~ 0, # did not experience
                     mh_y_ple_002==1 & mh_y_ple__exp_002==2 ~ 1, # experienced and mostly bad
                     mh_y_ple_002==1 & mh_y_ple__exp_002==1 ~ 0, # experienced and mostly good
                     mh_y_ple_002==1 & mh_y_ple__exp_002==999 ~ 0, # experienced and don't know
                     mh_y_ple_002==1 & mh_y_ple__exp_002==444 ~ NA, # experienced and not applicable
                     mh_y_ple_002==1 & mh_y_ple__exp_002=="n/a" ~ NA # experienced and NA
           )) %>%
  mutate(saw_crime_accident = 
           case_when(is.na(mh_y_ple_003) ~ NA, #NA for experienced or not
                     mh_y_ple_003==777 ~ NA, # decline to answer
                     mh_y_ple_003==999 ~ NA, # don't know
                     mh_y_ple_003==0 ~ 0, # did not experience
                     mh_y_ple_003==1 & mh_y_ple__exp_003==2 ~ 1, # experienced and mostly bad
                     mh_y_ple_003==1 & mh_y_ple__exp_003==1 ~ 0, # experienced and mostly good
                     mh_y_ple_003==1 & mh_y_ple__exp_003==999 ~ 0, # experienced and don't know
                     mh_y_ple_003==1 & mh_y_ple__exp_003==444 ~ NA, # experienced and not applicable
                     mh_y_ple_003==1 & mh_y_ple__exp_003=="n/a" ~ NA # experienced and NA
           )) %>%
  mutate(friend_died = 
           case_when(is.na(mh_y_ple_004) ~ NA, #NA for experienced or not
                     mh_y_ple_004==777 ~ NA, # decline to answer
                     mh_y_ple_004==999 ~ NA, # don't know
                     mh_y_ple_004==0 ~ 0, # did not experience
                     mh_y_ple_004==1 & mh_y_ple__exp_004==2 ~ 1, # experienced and mostly bad
                     mh_y_ple_004==1 & mh_y_ple__exp_004==1 ~ 0, # experienced and mostly good
                     mh_y_ple_004==1 & mh_y_ple__exp_004==999 ~ 0, # experienced and don't know
                     mh_y_ple_004==1 & mh_y_ple__exp_004==444 ~ NA, # experienced and not applicable
                     mh_y_ple_004==1 & mh_y_ple__exp_004=="n/a" ~ NA # experienced and NA
           )) %>%
  mutate(lost_friend = 
           case_when(is.na(mh_y_ple_005) ~ NA, #NA for experienced or not
                     mh_y_ple_005==777 ~ NA, # decline to answer
                     mh_y_ple_005==999 ~ NA, # don't know
                     mh_y_ple_005==0 ~ 0, # did not experience
                     mh_y_ple_005==1 & mh_y_ple__exp_005==2 ~ 1, # experienced and mostly bad
                     mh_y_ple_005==1 & mh_y_ple__exp_005==1 ~ 0, # experienced and mostly good
                     mh_y_ple_005==1 & mh_y_ple__exp_005==999 ~ 0, # experienced and don't know
                     mh_y_ple_005==1 & mh_y_ple__exp_005==444 ~ NA, # experienced and not applicable
                     mh_y_ple_005==1 & mh_y_ple__exp_005=="n/a" ~ NA # experienced and NA
           )) %>%
  mutate(friend_sick_injured = 
           case_when(is.na(mh_y_ple_006) ~ NA, #NA for experienced or not
                     mh_y_ple_006==777 ~ NA, # decline to answer
                     mh_y_ple_006==999 ~ NA, # don't know
                     mh_y_ple_006==0 ~ 0, # did not experience
                     mh_y_ple_006==1 & mh_y_ple__exp_006==2 ~ 1, # experienced and mostly bad
                     mh_y_ple_006==1 & mh_y_ple__exp_006==1 ~ 0, # experienced and mostly good
                     mh_y_ple_006==1 & mh_y_ple__exp_006==999 ~ 0, # experienced and don't know
                     mh_y_ple_006==1 & mh_y_ple__exp_006==444 ~ NA, # experienced and not applicable
                     mh_y_ple_006==1 & mh_y_ple__exp_006=="n/a" ~ NA # experienced and NA
           )) %>%
  mutate(parent_financial = 
           case_when(is.na(mh_y_ple_007) ~ NA, #NA for experienced or not
                     mh_y_ple_007==777 ~ NA, # decline to answer
                     mh_y_ple_007==999 ~ NA, # don't know
                     mh_y_ple_007==0 ~ 0, # did not experience
                     mh_y_ple_007==1 & mh_y_ple__exp_007==2 ~ 1, # experienced and mostly bad
                     mh_y_ple_007==1 & mh_y_ple__exp_007==1 ~ 0, # experienced and mostly good
                     mh_y_ple_007==1 & mh_y_ple__exp_007==999 ~ 0, # experienced and don't know
                     mh_y_ple_007==1 & mh_y_ple__exp_007==444 ~ NA, # experienced and not applicable
                     mh_y_ple_007==1 & mh_y_ple__exp_007=="n/a" ~ NA # experienced and NA
           )) %>%
  mutate(family_drug_alcohol = 
           case_when(is.na(mh_y_ple_008) ~ NA, #NA for experienced or not
                     mh_y_ple_008==777 ~ NA, # decline to answer
                     mh_y_ple_008==999 ~ NA, # don't know
                     mh_y_ple_008==0 ~ 0, # did not experience
                     mh_y_ple_008==1 & mh_y_ple__exp_008==2 ~ 1, # experienced and mostly bad
                     mh_y_ple_008==1 & mh_y_ple__exp_008==1 ~ 0, # experienced and mostly good
                     mh_y_ple_008==1 & mh_y_ple__exp_008==999 ~ 0, # experienced and don't know
                     mh_y_ple_008==1 & mh_y_ple__exp_008==444 ~ NA, # experienced and not applicable
                     mh_y_ple_008==1 & mh_y_ple__exp_008=="n/a" ~ NA # experienced and NA
           )) %>%
  mutate(sick = 
           case_when(is.na(mh_y_ple_009) ~ NA, #NA for experienced or not
                     mh_y_ple_009==777 ~ NA, # decline to answer
                     mh_y_ple_009==999 ~ NA, # don't know
                     mh_y_ple_009==0 ~ 0, # did not experience
                     mh_y_ple_009==1 & mh_y_ple__exp_009==2 ~ 1, # experienced and mostly bad
                     mh_y_ple_009==1 & mh_y_ple__exp_009==1 ~ 0, # experienced and mostly good
                     mh_y_ple_009==1 & mh_y_ple__exp_009==999 ~ 0, # experienced and don't know
                     mh_y_ple_009==1 & mh_y_ple__exp_009==444 ~ NA, # experienced and not applicable
                     mh_y_ple_009==1 & mh_y_ple__exp_009=="n/a" ~ NA # experienced and NA
           )) %>%
  mutate(injured = 
           case_when(is.na(mh_y_ple_010) ~ NA, #NA for experienced or not
                     mh_y_ple_010==777 ~ NA, # decline to answer
                     mh_y_ple_010==999 ~ NA, # don't know
                     mh_y_ple_010==0 ~ 0, # did not experience
                     mh_y_ple_010==1 & mh_y_ple__exp_010==2 ~ 1, # experienced and mostly bad
                     mh_y_ple_010==1 & mh_y_ple__exp_010==1 ~ 0, # experienced and mostly good
                     mh_y_ple_010==1 & mh_y_ple__exp_010==999 ~ 0, # experienced and don't know
                     mh_y_ple_010==1 & mh_y_ple__exp_010==444 ~ NA, # experienced and not applicable
                     mh_y_ple_010==1 & mh_y_ple__exp_010=="n/a" ~ NA # experienced and NA
           )) %>%
  mutate(parent_argue = 
           case_when(is.na(mh_y_ple_011) ~ NA, #NA for experienced or not
                     mh_y_ple_011==777 ~ NA, # decline to answer
                     mh_y_ple_011==999 ~ NA, # don't know
                     mh_y_ple_011==0 ~ 0, # did not experience
                     mh_y_ple_011==1 & mh_y_ple__exp_011==2 ~ 1, # experienced and mostly bad
                     mh_y_ple_011==1 & mh_y_ple__exp_011==1 ~ 0, # experienced and mostly good
                     mh_y_ple_011==1 & mh_y_ple__exp_011==999 ~ 0, # experienced and don't know
                     mh_y_ple_011==1 & mh_y_ple__exp_011==444 ~ NA, # experienced and not applicable
                     mh_y_ple_011==1 & mh_y_ple__exp_011=="n/a" ~ NA # experienced and NA
           )) %>%
  mutate(parent_lost_job = 
           case_when(is.na(mh_y_ple_012) ~ NA, #NA for experienced or not
                     mh_y_ple_012==777 ~ NA, # decline to answer
                     mh_y_ple_012==999 ~ NA, # don't know
                     mh_y_ple_012==0 ~ 0, # did not experience
                     mh_y_ple_012==1 & mh_y_ple__exp_012==2 ~ 1, # experienced and mostly bad
                     mh_y_ple_012==1 & mh_y_ple__exp_012==1 ~ 0, # experienced and mostly good
                     mh_y_ple_012==1 & mh_y_ple__exp_012==999 ~ 0, # experienced and don't know
                     mh_y_ple_012==1 & mh_y_ple__exp_012==444 ~ NA, # experienced and not applicable
                     mh_y_ple_012==1 & mh_y_ple__exp_012=="n/a" ~ NA # experienced and NA
           )) %>%
  mutate(parent_away =
           case_when(is.na(mh_y_ple_013) ~ NA, #NA for experienced or not
                     mh_y_ple_013==777 ~ NA, # decline to answer
                     mh_y_ple_013==999 ~ NA, # don't know
                     mh_y_ple_013==0 ~ 0, # did not experience
                     mh_y_ple_013==1 & mh_y_ple__exp_013==2 ~ 1, # experienced and mostly bad
                     mh_y_ple_013==1 & mh_y_ple__exp_013==1 ~ 0, # experienced and mostly good
                     mh_y_ple_013==1 & mh_y_ple__exp_013==999 ~ 0, # experienced and don't know
                     mh_y_ple_013==1 & mh_y_ple__exp_013==444 ~ NA, # experienced and not applicable
                     mh_y_ple_013==1 & mh_y_ple__exp_013=="n/a" ~ NA # experienced and NA
           )) %>%
  mutate(family_arrest =
           case_when(is.na(mh_y_ple_014) ~ NA, #NA for experienced or not
                     mh_y_ple_014==777 ~ NA, # decline to answer
                     mh_y_ple_014==999 ~ NA, # don't know
                     mh_y_ple_014==0 ~ 0, # did not experience
                     mh_y_ple_014==1 & mh_y_ple__exp_014==2 ~ 1, # experienced and mostly bad
                     mh_y_ple_014==1 & mh_y_ple__exp_014==1 ~ 0, # experienced and mostly good
                     mh_y_ple_014==1 & mh_y_ple__exp_014==999 ~ 0, # experienced and don't know
                     mh_y_ple_014==1 & mh_y_ple__exp_014==444 ~ NA, # experienced and not applicable
                     mh_y_ple_014==1 & mh_y_ple__exp_014=="n/a" ~ NA # experienced and NA
           )) %>%
  mutate(family_mental_emotional = 
           case_when(is.na(mh_y_ple_015) ~ NA, #NA for experienced or not
                     mh_y_ple_015==777 ~ NA, # decline to answer
                     mh_y_ple_015==999 ~ NA, # don't know
                     mh_y_ple_015==0 ~ 0, # did not experience
                     mh_y_ple_015==1 & mh_y_ple__exp_015==2 ~ 1, # experienced and mostly bad
                     mh_y_ple_015==1 & mh_y_ple__exp_015==1 ~ 0, # experienced and mostly good
                     mh_y_ple_015==1 & mh_y_ple__exp_015==999 ~ 0, # experienced and don't know
                     mh_y_ple_015==1 & mh_y_ple__exp_015==444 ~ NA, # experienced and not applicable
                     mh_y_ple_015==1 & mh_y_ple__exp_015=="n/a" ~ NA # experienced and NA
           )) %>%
  mutate(sibling_left = 
           case_when(is.na(mh_y_ple_016) ~ NA, #NA for experienced or not
                     mh_y_ple_016==777 ~ NA, # decline to answer
                     mh_y_ple_016==999 ~ NA, # don't know
                     mh_y_ple_016==0 ~ 0, # did not experience
                     mh_y_ple_016==1 & mh_y_ple__exp_016==2 ~ 1, # experienced and mostly bad
                     mh_y_ple_016==1 & mh_y_ple__exp_016==1 ~ 0, # experienced and mostly good
                     mh_y_ple_016==1 & mh_y_ple__exp_016==999 ~ 0, # experienced and don't know
                     mh_y_ple_016==1 & mh_y_ple__exp_016==444 ~ NA, # experienced and not applicable
                     mh_y_ple_016==1 & mh_y_ple__exp_016=="n/a" ~ NA # experienced and NA
           )) %>%
  mutate(victim_crime_violence_assault = 
           case_when(is.na(mh_y_ple_017) ~ NA, #NA for experienced or not
                     mh_y_ple_017==777 ~ NA, # decline to answer
                     mh_y_ple_017==999 ~ NA, # don't know
                     mh_y_ple_017==0 ~ 0, # did not experience
                     mh_y_ple_017==1 & mh_y_ple__exp_017==2 ~ 1, # experienced and mostly bad
                     mh_y_ple_017==1 & mh_y_ple__exp_017==1 ~ 0, # experienced and mostly good
                     mh_y_ple_017==1 & mh_y_ple__exp_017==999 ~ 0, # experienced and don't know
                     mh_y_ple_017==1 & mh_y_ple__exp_017==444 ~ NA, # experienced and not applicable
                     mh_y_ple_017==1 & mh_y_ple__exp_017=="n/a" ~ NA # experienced and NA
           )) %>%
  mutate(parent_separate_divorce = 
           case_when(is.na(mh_y_ple_018) ~ NA, #NA for experienced or not
                     mh_y_ple_018==777 ~ NA, # decline to answer
                     mh_y_ple_018==999 ~ NA, # don't know
                     mh_y_ple_018==0 ~ 0, # did not experience
                     mh_y_ple_018==1 & mh_y_ple__exp_018==2 ~ 1, # experienced and mostly bad
                     mh_y_ple_018==1 & mh_y_ple__exp_018==1 ~ 0, # experienced and mostly good
                     mh_y_ple_018==1 & mh_y_ple__exp_018==999 ~ 0, # experienced and don't know
                     mh_y_ple_018==1 & mh_y_ple__exp_018==444 ~ NA, # experienced and not applicable
                     mh_y_ple_018==1 & mh_y_ple__exp_018=="n/a" ~ NA # experienced and NA
           )) %>%
  mutate(parent_trouble_law = 
           case_when(is.na(mh_y_ple_019) ~ NA, #NA for experienced or not
                     mh_y_ple_019==777 ~ NA, # decline to answer
                     mh_y_ple_019==999 ~ NA, # don't know
                     mh_y_ple_019==0 ~ 0, # did not experience
                     mh_y_ple_019==1 & mh_y_ple__exp_019==2 ~ 1, # experienced and mostly bad
                     mh_y_ple_019==1 & mh_y_ple__exp_019==1 ~ 0, # experienced and mostly good
                     mh_y_ple_019==1 & mh_y_ple__exp_019==999 ~ 0, # experienced and don't know
                     mh_y_ple_019==1 & mh_y_ple__exp_019==444 ~ NA, # experienced and not applicable
                     mh_y_ple_019==1 & mh_y_ple__exp_019=="n/a" ~ NA # experienced and NA
           )) %>%
  mutate(new_school = 
           case_when(is.na(mh_y_ple_020) ~ NA, #NA for experienced or not
                     mh_y_ple_020==777 ~ NA, # decline to answer
                     mh_y_ple_020==999 ~ NA, # don't know
                     mh_y_ple_020==0 ~ 0, # did not experience
                     mh_y_ple_020==1 & mh_y_ple__exp_020==2 ~ 1, # experienced and mostly bad
                     mh_y_ple_020==1 & mh_y_ple__exp_020==1 ~ 0, # experienced and mostly good
                     mh_y_ple_020==1 & mh_y_ple__exp_020==999 ~ 0, # experienced and don't know
                     mh_y_ple_020==1 & mh_y_ple__exp_020==444 ~ NA, # experienced and not applicable
                     mh_y_ple_020==1 & mh_y_ple__exp_020=="n/a" ~ NA # experienced and NA
           )) %>%
  mutate(family_moved = 
           case_when(is.na(mh_y_ple_021) ~ NA, #NA for experienced or not
                     mh_y_ple_021==777 ~ NA, # decline to answer
                     mh_y_ple_021==999 ~ NA, # don't know
                     mh_y_ple_021==0 ~ 0, # did not experience
                     mh_y_ple_021==1 & mh_y_ple__exp_021==2 ~ 1, # experienced and mostly bad
                     mh_y_ple_021==1 & mh_y_ple__exp_021==1 ~ 0, # experienced and mostly good
                     mh_y_ple_021==1 & mh_y_ple__exp_021==999 ~ 0, # experienced and don't know
                     mh_y_ple_021==1 & mh_y_ple__exp_021==444 ~ NA, # experienced and not applicable
                     mh_y_ple_021==1 & mh_y_ple__exp_021=="n/a" ~ NA # experienced and NA
           )) %>%
  mutate(parent_jail = 
           case_when(is.na(mh_y_ple_022) ~ NA, #NA for experienced or not
                     mh_y_ple_022==777 ~ NA, # decline to answer
                     mh_y_ple_022==999 ~ NA, # don't know
                     mh_y_ple_022==0 ~ 0, # did not experience
                     mh_y_ple_022==1 & mh_y_ple__exp_022==2 ~ 1, # experienced and mostly bad
                     mh_y_ple_022==1 & mh_y_ple__exp_022==1 ~ 0, # experienced and mostly good
                     mh_y_ple_022==1 & mh_y_ple__exp_022==999 ~ 0, # experienced and don't know
                     mh_y_ple_022==1 & mh_y_ple__exp_022==444 ~ NA, # experienced and not applicable
                     mh_y_ple_022==1 & mh_y_ple__exp_022=="n/a" ~ NA # experienced and NA
           )) %>%
  mutate(new_step_parent = 
           case_when(is.na(mh_y_ple_023) ~ NA, #NA for experienced or not
                     mh_y_ple_023==777 ~ NA, # decline to answer
                     mh_y_ple_023==999 ~ NA, # don't know
                     mh_y_ple_023==0 ~ 0, # did not experience
                     mh_y_ple_023==1 & mh_y_ple__exp_023==2 ~ 1, # experienced and mostly bad
                     mh_y_ple_023==1 & mh_y_ple__exp_023==1 ~ 0, # experienced and mostly good
                     mh_y_ple_023==1 & mh_y_ple__exp_023==999 ~ 0, # experienced and don't know
                     mh_y_ple_023==1 & mh_y_ple__exp_023==444 ~ NA, # experienced and not applicable
                     mh_y_ple_023==1 & mh_y_ple__exp_023=="n/a" ~ NA # experienced and NA
           )) %>%
  mutate(parent_new_job = 
           case_when(is.na(mh_y_ple_024) ~ NA, #NA for experienced or not
                     mh_y_ple_024==777 ~ NA, # decline to answer
                     mh_y_ple_024==999 ~ NA, # don't know
                     mh_y_ple_024==0 ~ 0, # did not experience
                     mh_y_ple_024==1 & mh_y_ple__exp_024==2 ~ 1, # experienced and mostly bad
                     mh_y_ple_024==1 & mh_y_ple__exp_024==1 ~ 0, # experienced and mostly good
                     mh_y_ple_024==1 & mh_y_ple__exp_024==999 ~ 0, # experienced and don't know
                     mh_y_ple_024==1 & mh_y_ple__exp_024==444 ~ NA, # experienced and not applicable
                     mh_y_ple_024==1 & mh_y_ple__exp_024=="n/a" ~ NA # experienced and NA
           )) %>%
  mutate(new_sibling = 
           case_when(is.na(mh_y_ple_025) ~ NA, #NA for experienced or not
                     mh_y_ple_025==777 ~ NA, # decline to answer
                     mh_y_ple_025==999 ~ NA, # don't know
                     mh_y_ple_025==0 ~ 0, # did not experience
                     mh_y_ple_025==1 & mh_y_ple__exp_025==2 ~ 1, # experienced and mostly bad
                     mh_y_ple_025==1 & mh_y_ple__exp_025==1 ~ 0, # experienced and mostly good
                     mh_y_ple_025==1 & mh_y_ple__exp_025==999 ~ 0, # experienced and don't know
                     mh_y_ple_025==1 & mh_y_ple__exp_025==444 ~ NA, # experienced and not applicable
                     mh_y_ple_025==1 & mh_y_ple__exp_025=="n/a" ~ NA # experienced and NA
           )) %>%
  # year 3 only: mh_y_ple__exp_026, mh_y_ple__exp_027, mh_y_ple__exp_028, mh_y_ple__exp_029, mh_y_ple__exp_030, mh_y_ple__exp_031
  # year 4 only: mh_y_ple__exp_032, mh_y_ple__exp_033
  # no item for good/bad effect for mh_y_ple__exp_034 so not included in cumulative score
  filter(!is.na(family_died)) %>% #before = 10946, after = 10329
  filter(!is.na(family_injured)) %>% #after = 10093
  filter(!is.na(saw_crime_accident)) %>% #after = 9687
  filter(!is.na(friend_died)) %>% #after = 9681
  filter(!is.na(lost_friend)) %>% #after = 9599
  filter(!is.na(friend_sick_injured)) %>% #after = 9546
  filter(!is.na(parent_financial)) %>% #after = 9516
  filter(!is.na(family_drug_alcohol)) %>% #after = 9429
  filter(!is.na(sick)) %>% #after = 9410
  filter(!is.na(injured)) %>% #after = 9378
  filter(!is.na(parent_argue)) %>% #after = 9343
  filter(!is.na(parent_lost_job)) %>% #after = 9245
  filter(!is.na(parent_away)) %>% #after = 9073
  filter(!is.na(family_arrest)) %>% #after = 8999
  filter(!is.na(family_mental_emotional)) %>% #after = 8941
  filter(!is.na(sibling_left)) %>% #after = 8861
  filter(!is.na(victim_crime_violence_assault)) %>% #after = 8859
  filter(!is.na(parent_separate_divorce)) %>% #after = 8760
  filter(!is.na(parent_trouble_law)) %>% #after = 8746
  filter(!is.na(new_school)) %>% #after = 8598
  filter(!is.na(family_moved)) %>% #after = 8492
  filter(!is.na(parent_jail)) %>% #after = 8480
  filter(!is.na(new_step_parent)) %>% #after = 8452
  filter(!is.na(parent_new_job)) %>% #after = 8367
  filter(!is.na(new_sibling)) %>% #after = 8333
  mutate(total_bad_le = rowSums(across(c(family_died,family_injured,saw_crime_accident,
                                         friend_died,lost_friend,friend_sick_injured,
                                         parent_financial,family_drug_alcohol,sick,injured,
                                         parent_argue,parent_lost_job,parent_away,family_arrest,
                                         family_mental_emotional,sibling_left,
                                         victim_crime_violence_assault,parent_separate_divorce,
                                         parent_trouble_law,new_school,family_moved,parent_jail,
                                         new_step_parent,parent_new_job,new_sibling
  )))) %>%
  # remove rows with NA. before this step, n = 8333 and after n = 8333
  filter(!is.na(total_bad_le)) %>%
  # select only columns relevant to analysis ie participant ID, time point,
  # and total number of life events experienced as bad
  select(participant_id,session_id,total_bad_le) %>%
  mutate(total_bad_le = as.numeric(total_bad_le))

# see all unique values (can visually check for NA or errors)
map(ledata,unique)

### Combine all data for analysis into one data frame ####
analysis_data <- genderdata %>%
  # keep only gender data from year 4
  filter(session_id=="ses-04A") %>%
  # remove session_id
  select(-session_id) %>%
  # add longitudinal tracking data based on subject ID and data collection year.
  left_join(trackdata %>% filter(session_id=="ses-04A") %>% select(-session_id),
            by = c("participant_id")) %>%
  # add DERS-P data based on subject ID and data collection year. before this
  left_join(dersdata %>% select(-session_id),by=c("participant_id")) %>%
  # add LES data based on subject ID and data collection year
  left_join(ledata %>% select(-session_id),by=c("participant_id")) %>%
  # add BPM data based on subject ID and data collection year
  left_join(bpmdata %>% select(-session_id),by=c("participant_id")) %>%
  # remove subjects with NA for LES. Before this step, n = 4612 and
  # after this step, n = 3455
  filter(!(is.na(total_bad_le))) %>%
  # remove subjects with NA for year 3 DERS. Before this step, n = 3455 and
  # after this step, n = 3294
  filter(!(is.na(ders_total))) %>%
  # remove subjects with NA for year 4 BPM internalizing. Before this step, 
  # n = 3294 and after this step, n = 3189
  filter(!(is.na(bpm_int))) %>%
  # remove subjects with NA for year 4 BPM externalizing. Before this step, 
  # n = 3189 and after this step, n = 3189
  filter(!(is.na(bpm_ext))) %>%
  # Z score continuous variables
  mutate(across(
    c(age, ders_total,
      total_bad_le,
      bpm_int, bpm_ext,
    ),
    # Z-score and center continuous variables
    ~ as.numeric(scale(.,center=TRUE,scale=TRUE)),
    .names = "Z_{.col}"
  )) %>%
  # make categorical variables factors
  mutate(across(c(genderid, gender_details, sex, sex_details, site), as.factor))

# see all unique values (can visually check for NA or errors)
map(analysis_data,unique)

### Get general overview of all data ####

#### Get age for participants at each time point
ab_g_dyn %>% 
  filter(participant_id %in% analysis_data$participant_id) %>%
  select(participant_id,session_id,ab_g_dyn__visit_age) %>%
  group_by(session_id) %>%
  summarize(mean_age = round(mean(ab_g_dyn__visit_age),3),
            sd_age = round(sd(ab_g_dyn__visit_age),3)) %>%
  as.data.frame()

#### Determine how many subjects in each gender group ####
analysis_data %>%
  group_by(genderid) %>%
  # group_by(gender_details) %>%
  count()


#### Determine how many subjects in each combination of gender and sex group ####
analysis_data %>%
  group_by(genderid, sex) %>%
  count()

## PLOTS ####

##### Graph of LES vs DERS by gender ####
ggplot(analysis_data, 
       aes(x=total_bad_le,y=ders_total, fill=genderid)) +
  geom_smooth(aes(linetype=genderid),method="lm",color="black", se=TRUE) +
  scale_linetype_manual(values = c("cis_boy"="31",
                                   "cis_girl"="11",
                                   "gd"="solid")) +
  scale_shape_manual(values=c(21,22,23)) +
  scale_colour_grey(start=0.9,end=0) +
  scale_fill_manual(values=c("grey40","grey85","black")) +
  scale_x_continuous(expand = c(0,0),
                     breaks = seq(0,20, by = 2),
                     limits=c(0,20)) +
  coord_cartesian(ylim = c(29, 130)) +
  scale_y_continuous(breaks=seq(0,130,10)) +
  guides(
    shape = guide_legend(override.aes = list(size = 3)),
    line = guide_legend(override.aes = list(size = 2))
  ) +
  theme_classic() +
  theme(legend.key.width = unit(0.5, "in"))
# Save graph
# ggsave("les_vs_ders_gender.tiff",width=8.5,height=6,unit="in",path="figures")

##### Graph of LES vs DERS by sex ####
# Note: 24 subjects with sex as NA (but with gender data) will not be plotted
ggplot(analysis_data[-which(is.na(analysis_data$sex)),], 
       aes(x=total_bad_le,y=ders_total, fill=sex)) +
  geom_smooth(aes(linetype=sex),method="lm",color="black", se=TRUE) +
  scale_linetype_manual(values = c("male"="31",
                                   "female"="11")) +
  scale_shape_manual(values=c(21,22)) +
  scale_colour_grey(start=0.9,end=0) +
  scale_fill_manual(values=c("grey40","grey85")) +
  scale_x_continuous(expand = c(0,0),
                     breaks = seq(0,20, by = 2),
                     limits=c(0,20)) +
  coord_cartesian(ylim = c(29, 130)) +
  scale_y_continuous(breaks=seq(0,130,10)) +
  guides(
    shape = guide_legend(override.aes = list(size = 3)),
    line = guide_legend(override.aes = list(size = 2))
  ) +
  theme_classic() +
  theme(legend.key.width = unit(0.5, "in"))
# Save graph
# ggsave("les_vs_ders_sex.tiff",width=8.5,height=6,unit="in",path="figures")

##### Graphs of LES vs BPM by gender ####
outcome_list <- c("bpm_int","bpm_ext")
bpm_les_gender_plot_list <- list()
for (outcome in outcome_list) {
  # Create plot
  bpm_les_plot <-
    ggplot(aes(x=total_bad_le,y=.data[[outcome]],
               fill = genderid
    ),data=analysis_data) +
    geom_smooth(aes(linetype=genderid),method="lm",color="black", se=TRUE) +
    scale_linetype_manual(values = c("cis_boy"="31",
                                     "cis_girl"="11",
                                     "gd"="solid")) +
    scale_shape_manual(values=c(21,22,23)) +
    scale_colour_grey(start=0.9,end=0) +
    scale_fill_manual(values=c("grey40","grey85","black")) +
    scale_x_continuous(expand = c(0,0),
                       breaks = seq(0,20, by = 2),
                       limits=c(0,20)) +
    scale_y_continuous(expand = c(0,0),
                       breaks = seq(0,20, by = 2),
                       limits=c(0,14)) +
    theme_classic()
  bpm_les_gender_plot_list[[outcome]] <- bpm_les_plot
  # Save plot
  # ggsave(paste0("les_vs_",outcome,"_bygender.tiff"),
  # width=8.5,height=6,units = "in",path="figures")
}
bpm_les_gender_plot_list

##### Graphs of LES vs BPM by sex ####
# Note: 24 subjects with sex as NA (but with gender data) will not be plotted
bpm_les_sex_plot_list <- list()
for (outcome in outcome_list) {
  # Create plot
  bpm_les_plot <-
    ggplot(aes(x=total_bad_le,y=.data[[outcome]],
               fill = sex
    ),data=analysis_data[-which(is.na(analysis_data$sex)),]) +
    geom_smooth(aes(linetype=sex),method="lm",color="black", se=TRUE) +
    scale_linetype_manual(values = c("male"="31",
                                     "female"="11")) +
    scale_shape_manual(values=c(21,22)) +
    scale_colour_grey(start=0.9,end=0) +
    scale_fill_manual(values=c("grey40","grey85")) +
    scale_x_continuous(expand = c(0,0),
                       breaks = seq(0,20, by = 2),
                       limits=c(0,20)) +
    scale_y_continuous(expand = c(0,0),
                       breaks = seq(0,20, by = 2),
                       limits=c(0,14)) +
    coord_cartesian(ylim = c(0,14)) +
    theme_classic()
  bpm_les_sex_plot_list[[outcome]] <- bpm_les_plot
  # Save plot
  # ggsave(paste0("les_vs_",outcome,"_bysex.tiff"),
  # width=8.5,height=6,units = "in",path="figures")
}
bpm_les_sex_plot_list

##### Graphs of DERS vs BPM by gender ####
bpm_ders_gender_plot_list <- list()
for (outcome in outcome_list) {
  # Create plot
  bpm_ders_plot <- ggplot(aes(x=ders_total,y=.data[[outcome]],
                              linetype=genderid,
                              shape = genderid,
                              fill = genderid
  ),data=analysis_data) +
    geom_smooth(method="lm",se=TRUE,color="black") +
    scale_linetype_manual(values = c("cis_boy"="31",
                                     "cis_girl"="11",
                                     "gd"="solid")) +
    scale_shape_manual(values=c(21,22,23)) +
    scale_fill_manual(values=c("grey40","grey85","black")) +
    scale_colour_grey(start=0.9,end=0) +
    scale_x_continuous(expand = c(0,0),
                       breaks=seq(20,130,10),
                       limits = c(29,130)) +
    scale_y_continuous(expand = c(0,0),
                       breaks = seq(0,20, by = 2),
                       limits=c(0,14)) +
    coord_cartesian(ylim = c(0,14)) +
    theme_classic() +
    guides(shape = guide_legend(override.aes = list(size = 2)))
  bpm_ders_gender_plot_list[[outcome]] <- bpm_ders_plot
  # Save plot
  # ggsave(paste0("ders_vs_",outcome,"_bygender.tiff"),
  #        width=8.3,height=6,units = "in",path="figures")
}
bpm_ders_gender_plot_list

##### Graphs of DERS vs BPM by sex ####
# Note: 24 subjects with sex as NA (but with gender data) will not be plotted
bpm_ders_sex_plot_list <- list()
for (outcome in outcome_list) {
  # Create plot
  bpm_ders_plot <- ggplot(aes(x=ders_total,y=.data[[outcome]],
                              linetype=sex,
                              shape = sex,
                              fill = sex
  ),data=analysis_data[-which(is.na(analysis_data$sex)),]) +
    geom_smooth(method="lm",se=TRUE,color="black") +
    scale_linetype_manual(values = c("male"="31",
                                     "female"="11")) +
    scale_shape_manual(values=c(21,22,23)) +
    scale_fill_manual(values=c("grey40","grey85")) +
    scale_colour_grey(start=0.9,end=0) +
    scale_x_continuous(expand = c(0,0),
                       breaks=seq(20,130,10),
                       limits = c(29,130)) +
    scale_y_continuous(expand = c(0,0),
                       breaks = seq(0,20, by = 2),
                       limits=c(0,14)) +
    coord_cartesian(ylim = c(0,14)) +
    theme_classic() +
    guides(shape = guide_legend(override.aes = list(size = 2)))
  bpm_ders_sex_plot_list[[outcome]] <- bpm_ders_plot
  # Save plot
  # ggsave(paste0("ders_vs_",outcome,"_bysex.tiff"),
  #        width=8.3,height=6,units = "in",path="figures")
}
bpm_ders_sex_plot_list

## STEP ONE: BASIC STATS ####

### Get summary stats for each variable ####
sumstats <- 
  analysis_data %>% 
  # group_by(genderid) %>% #group by cis boy, cis girl, or nb and get n for each group
  # group_by(gender_details) %>% #group by cis boy, cis girl, nb, trans boy, or trans girl 
  filter(!is.na(sex)) %>% #remove 24 subjects with NA sex
  group_by(sex) %>% #group by male or female
  # group_by(sex_details) %>% #group by sex including subjects with don't know or refuse for sex
  summarise(
    n = n(),
    across(
      c("age","total_bad_le","ders_total","bpm_int","bpm_ext"),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd = ~sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    )
  ) %>%
  # transpose so groups are columns and summary stats are rows
  t() %>% 
  # make data frame so values are not strings
  as.data.frame() %>%
  # make first row ie group names into column names
  set_names(.[1, ]) %>%
  # remove first row which is now column names
  slice(-1) %>%
  # make all values numeric
  mutate(across(everything(), as.numeric)) %>%
  # round all values to two decimal places
  mutate(across(where(is.numeric), ~ round(.x, 2)))
sumstats

### Assess normality of distributions of each variable from all data for methods ####
# p-value < 0.05 suggests data is not normally distributed
# walk() applies a function, in this case shapiro.test, to a list
# cat() makes it print the variable name before each test
walk(c("total_bad_le","ders_total","bpm_int","bpm_ext"), 
     ~ {
       cat("Variable:", .x, "\n")
       print(shapiro.test(analysis_data[[.x]]))
     })

### Create basic histogram for each variable ####
# Store histograms in list in case want to view later
variable_histograms <- 
  map(c("total_bad_le", 
        "ders_total", 
        "bpm_int","bpm_ext"
        # "cbcl_int","cbcl_ext"
  ), 
  # note that !!sym(.x) turns the variables in the list above into arguments
  # that can be passed to ggplot
  ~ ggplot(analysis_data, aes(x = !!sym(.x))) + 
    geom_histogram() +
    ggtitle(.x)
  )
# Print all histograms from stored list
print(variable_histograms)

### Create correlation matrix for all variables from year 4 data for results and table 1 ####
# Make correlation matrix
corrmat <- 
  analysis_data %>%
  # select relevant columns
  select(c(age,
           total_bad_le,
           ders_total,
           bpm_int,bpm_ext
           # cbcl_int,cbcl_ext
  )) %>% 
  # run correlation tests for all pairs of variables, adjust using fdr
  corr.test(adjust="fdr")
# print correlation matrix, correlation coefficients, and p-values
print(corrmat,digits=3)
# See fdr-adjusted p-value for each correlation test
# note that list gives p values above diagonal, going across down columns above diagonal
round(corrmat$p.adj,5)

### Kruskal-Wallis (non-parametric version of one-way ANOVA) and Dunn test ie ####
### pairwise Wilcoxon tests to determine whether variables differ based on gender

#### Age
kruskal.test(age ~ genderid, data = analysis_data) 

#### LES  
kruskal.test(total_bad_le ~ genderid, data = analysis_data)
dunnTest(analysis_data$total_bad_le, analysis_data$genderid, method = "bh")

#### DERS  
kruskal.test(ders_total ~ genderid, data = analysis_data)
dunnTest(analysis_data$ders_total, analysis_data$genderid, method = "bh")

#### BPM internalizing 
kruskal.test(bpm_int ~ genderid, data = analysis_data)
dunnTest(analysis_data$bpm_int, analysis_data$genderid, method = "bh")

#### BPM externalizing
kruskal.test(bpm_ext ~ genderid, data = analysis_data)
dunnTest(analysis_data$bpm_ext, analysis_data$genderid, method = "bh")


### Mann-Whitney U (non-parametric version of two-sample t-test) to determine ####
### whether variables differ based on sex 

#### Age 
wilcox.test(age ~ sex, data = analysis_data)

#### LES
wilcox.test(total_bad_le ~ sex, data = analysis_data)

#### DERS
wilcox.test(ders_total ~ sex, data = analysis_data)

#### BPM internalizing
wilcox.test(bpm_int ~ sex, data = analysis_data)

#### BPM externalizing 
wilcox.test(bpm_ext ~ sex, data = analysis_data)

## STEP TWO: MEDIATING EFFECT OF ER ON BPM ~ LES #### 

### Simple mediation model for bpm internalizing ####

bpm_int_gender_model4 <- PROCESS(
  analysis_data,
  y = "Z_bpm_int",
  x = "Z_total_bad_le",
  meds = c("Z_ders_total"),
  covs = c("Z_age"), # neither gender nor sex as covariate
  clusters = "site",
  cov.path = c("both"),
  nsim = 5000,
  seed = 1234,
  digits = 5)

### Simple mediation model for bpm externalizing ####
bpm_ext_gender_model4 <- PROCESS(
  analysis_data,
  y = "Z_bpm_ext",
  x = "Z_total_bad_le",
  meds = c("Z_ders_total"),
  covs = c("Z_age"), # neither gender nor sex as covariate
  clusters = "site",
  cov.path = c("both"),
  nsim = 5000,
  seed = 1234,
  digits = 5)

## STEP THREE: MODERATING EFFECT OF GENDER OR SEX ON MEDIATION ####

### Moderated mediation model (Hayes model 15) to test whether gender ####
### moderates mediating effect of DERS on relationship between LES and BPM
### internalizing 
bpm_int_gender_model15 <- PROCESS(
  analysis_data,
  y = "Z_bpm_int",
  x = "Z_total_bad_le",
  meds = c("Z_ders_total"),
  mods = c("genderid"),
  covs = c("Z_age"),
  clusters = "site",
  # mod.path = c("x-y","x-m","m-y"), #only sig interaction for gender*les so run line below
  mod.path = c("x-y"),
  cov.path = c("both"),
  nsim = 5000,
  seed = 1234,
  digits = 5)

# are the conditional direct effects [c'] of X on Y significant different for
# different groups? Use Z test to find out. 
# Z = (beta1 - beta2 / (sqrt(SE1^2 + SE2^2)))
# cis boy (beta1) vs cis girl (beta2): 
#     Z = (0.05920      - 0.15855     )/(sqrt((0.02333^2)+(0.02460)^2)) = -2.930374
# cis boy (beta1) vs gd (beta2):
#     Z = (0.05920      - 0.22014     )/(sqrt((0.02333^2)+(0.05864)^2)) = -2.550129
# cis girl (beta1) vs gd (beta2):
#     Z = (0.15855      - 0.22014     )/(sqrt((0.02460^2)+(0.05864)^2)) = -0.968534
# to go from Z score to p-value, find probability of being outside absolute value
# of Z score (because don't know if beta1 is smaller or larger than beta2) and 
# then multiply that by 2 because two-tailed test. Can use default settings of
# mean = 0 and sd = 1 in pnorm function because that is true of Z scores
# cis boy vs cis girl: Z = -2.930374, so pnorm(-abs(-2.930374))*2 = 0.003385543
# cis boy vs gd: Z = -2.550129, so pnorm(-abs(-2.550129))*2 = 0.01076831
# cis girl vs gd: Z = -0.968534, so pnorm(-abs(-0.968534))*2 = 0.3327777
# Finally, we need to fdr correct for multiple tests:
# p.adjust(c(0.003385543,0.01076831,0.3327777),method="fdr")
# So final p-values rounded to three places are:
# cis boy vs cis girl: p = .010
# cis boy vs gd: p = .016
# cis girl vs gd: p = .333

### Moderated mediation model (Hayes model 15) to test whether gender ####
### moderates mediating effect of DERS on relationship between LES and BPM
### externalizing  
bpm_ext_gender_model15 <- PROCESS(
  analysis_data,
  y = "Z_bpm_ext",
  x = "Z_total_bad_le",
  meds = c("Z_ders_total"),
  mods = c("genderid"),
  covs = c("Z_age"),
  clusters = "site",
  mod.path = c("x-y","x-m","m-y"), # no sig interactions
  cov.path = c("both"),
  nsim = 5000,
  seed = 1234,
  digits = 5)

### Moderated mediation model (Hayes model 15) to test whether sex ####
### moderates mediating effect of DERS on relationship between LES and BPM
### internalizing 
bpm_int_sex_model15 <- PROCESS(
  analysis_data,
  y = "Z_bpm_int",
  x = "Z_total_bad_le",
  meds = c("Z_ders_total"),
  mods = c("sex"),
  covs = c("Z_age"),
  clusters = "site",
  # mod.path = c("x-y","x-m","m-y"), #sig interaction only for sex*les on bpm and sex*ders on bpm so run line below
  mod.path = c("x-y","m-y"),
  cov.path = c("both"),
  nsim = 5000,
  seed = 1234,
  digits = 5)

### Moderated mediation model (Hayes model 15) to test whether sex ####
### moderates mediating effect of DERS on relationship between LES and BPM
### externalizing 
bpm_ext_sex_model15 <- PROCESS(
  analysis_data,
  y = "Z_bpm_ext",
  x = "Z_total_bad_le",
  meds = c("Z_ders_total"),
  mods = c("sex"),
  covs = c("Z_age"),
  clusters = "site",
  mod.path = c("x-y","x-m","m-y"), #no sig interactions
  cov.path = c("both"),
  nsim = 5000,
  seed = 1234,
  digits = 5)