### Set working directory ####
setwd("C:/Users/Kate Scheuer/OneDrive - UW/Desktop/Lab/aces_emotionreg_cbcl")

### Load libraries ####
library(tidyverse) # for dplyr and associated functions
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
gish_y_gi <- read_csv("data/gish_y_gi.csv")

#### DERS-P for emotion regulation ####
mh_p_ders <- read_csv("data/mh_p_ders.csv")

#### CBCL for parent-report psychopathology symptoms ####
# Note: warning will flag 226 problems total from rows 9941, 14373, 14660, 
# 16373, 17952, 21439, 23991, 28582, 31278, 36700, 38672, 40583, 42272, and 
# 46007 and from 20 total columns all ending in _m indicating they pertain to 
# missing values, so it is okay to ignore this warning
# problems(mh_p_cbcl)
# unique(problems(mh_p_cbcl)$row)
# colnames(mh_p_cbcl[,unique(problems(mh_p_cbcl)$col)])
mh_p_cbcl <- read_csv("data/mh_p_cbcl.csv")

#### BPM for youth-report psychopathology symptoms ####
mh_y_bpm <- read_csv("data/mh_y_bpm.csv")

#### Longitudinal tracking data ####
# Note: warning will flag row 79003 and column 5 ie rel_birth_id which contains
# NA, so it is okay to ignore this warning
# problems(abcd_y_lt)
# abcd_y_lt[problems(abcd_y_lt)$row,problems(abcd_y_lt)$col]
abcd_y_lt <- read_csv("data/abcd_y_lt.csv")

#### LES (youth-reported) ####
#### LES for exposure to negative life events
# Note: warning will flag 76 problems total from 47 unique rows and columns 109
# and 110 ie ple_homeless_fu_y and ple_homeless_fu2_y, neither of which is 
# relevant for this analysis, so it is okay to ignore this warning
# problems(mh_y_le)
# unique(problems(mh_y_le)$row)
# colnames(mh_y_le[,unique(problems(mh_y_le)$col)])
mh_y_le <- read_csv("data/mh_y_le.csv")

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
  # keep only columns relevant to analysis
  select(src_subject_id,eventname,
         sex,sex_details,
         genderid,gender_details
  ) %>%
  
  ### to get percentage gd youth in year 3 vs year 4 for methods and to get number
  ### of subjects who said refuse or don't know to gender question at year 4
  ### for methods, uncomment and run just the line below
  # group_by(eventname,genderid) %>% count() %>% print()
  ### to get number of subjects who said refuse or don't know to sex question at
  ### year 4 for methods, uncomment and run just the line below
  # group_by(eventname,sex_details) %>% count() %>% print()
  
  # remove subjects who refused to answer and/or did not understand gender
  # or trans questions. Before this step, n should be 15064. After this step,
# n should be 14495.
filter(genderid!="refuse",
       genderid!="dont_understand") %>%
  # make genderid a factor (automatically uses cisboy as reference)
  mutate(genderid_refcisboy = as.factor(genderid)) %>%
  # make another column that uses cis_girl as reference instead
  mutate(genderid_refcisgirl = relevel(genderid_refcisboy,"cis_girl"))

#### Count number of subjects per gender group per data collection year ####
# 3 year follow-up: 5113 cis boys, 4462 cis girls, 308 gd
# 4 year follow-up: 2393 cis boys, 1964 cis girls, 255 gd
genderdata %>% 
  group_by(eventname) %>% 
  count(genderid)

### to get number of subjects with CBCL data in year 4 who refused to answer
### at least one item (okay if missing language)
mh_p_cbcl %>% group_by(eventname) %>% count()
mh_p_cbcl %>% group_by(eventname) %>% 
  filter(if_any(all_of(c("cbcl_scr_syn_internal_r","cbcl_scr_syn_external_r")), is.na)) %>% count()

### Prepare CBCL data for analysis ####
cbcldata <- mh_p_cbcl %>%
  # select only columns relevant to analysis
  select(src_subject_id,eventname,
         cbcl_scr_syn_internal_r,
         cbcl_scr_syn_external_r
  ) %>%
  # rename subscale columns to be more human-readable and shorter
  rename(cbcl_ext = cbcl_scr_syn_external_r,
         cbcl_int = cbcl_scr_syn_internal_r
  )  

### to get number of subjects with CBCL data in year 4 who refused to answer
### at least one item (okay if missing language)
mh_y_bpm %>% group_by(eventname) %>% count()
mh_y_bpm %>% group_by(eventname) %>% 
  filter(if_any(all_of(c("bpm_y_scr_internal_r","bpm_y_scr_external_r")), is.na)) %>% count()


### Prepare BPM data for analysis ####
bpmdata <- mh_y_bpm %>%
  # select only columns relevant to analysis
  select(src_subject_id,eventname,
         bpm_y_scr_internal_r,
         bpm_y_scr_external_r
  ) %>%
  # rename subscale columns to be more human-readable and shorter
  rename(bpm_ext = bpm_y_scr_external_r,
         bpm_int = bpm_y_scr_internal_r
  ) 

### Prepare DERS-P data for analysis ####

### to get number of subjects with DERS data in year 3 who refused to answer
### at least one item (okay if missing language)
mh_p_ders %>% group_by(eventname) %>% count()
mh_p_ders %>% group_by(eventname) %>% 
  filter(if_any(-all_of(c("ders_p_select_language___1")), ~ . == 777)) %>% count()

#### Create cumulative score ####
dersdata <- mh_p_ders %>%
  # remove subjects who refused to answer one or more items. Before this step,
  # n should be 14708. After this step, n should be 14225.
  filter(!if_any(everything(), ~ . == 777)) %>%
  # remove subjects with na in any row. Before this step, n should be 14225.
  # After this step, n should still be 14225.
  drop_na() %>%
  # add column to reverse score "my child pays attention to how he/she feels"
  mutate(rev_ders_attn_awareness_p = 
           case_when(ders_attn_awareness_p == 1 ~ 5,
                     ders_attn_awareness_p == 2 ~ 4,
                     ders_attn_awareness_p == 3 ~ 3,
                     ders_attn_awareness_p == 4 ~ 2,
                     ders_attn_awareness_p == 5 ~ 1)) %>%
  # add column to reverse score "my child is attentive to his/her feelings"
  mutate(rev_ders_feelings_attentive_p = 
           case_when(ders_feelings_attentive_p == 1 ~ 5,
                     ders_feelings_attentive_p == 2 ~ 4,
                     ders_feelings_attentive_p == 3 ~ 3,
                     ders_feelings_attentive_p == 4 ~ 2,
                     ders_feelings_attentive_p == 5 ~ 1)) %>%
  # add column to reverse score "my child cares about what he/she is feeling"
  mutate(rev_ders_feelings_care_p = 
           case_when(ders_feelings_care_p == 1 ~ 5,
                     ders_feelings_care_p == 2 ~ 4,
                     ders_feelings_care_p == 3 ~ 3,
                     ders_feelings_care_p == 4 ~ 2,
                     ders_feelings_care_p == 5 ~ 1)) %>%
  # add column to reverse score "when my child is upset, he/she acknowledges
  # his/her emotions"
  mutate(rev_ders_upset_ack_p = 
           case_when(ders_upset_ack_p == 1 ~ 5,
                     ders_upset_ack_p == 2 ~ 4,
                     ders_upset_ack_p == 3 ~ 3,
                     ders_upset_ack_p == 4 ~ 2,
                     ders_upset_ack_p == 5 ~ 1)) %>%
  # add column to reverse score "my child is clear about his/her feelings"
  mutate(rev_ders_clear_feelings_p = 
           case_when(ders_clear_feelings_p == 1 ~ 5,
                     ders_clear_feelings_p == 2 ~ 4,
                     ders_clear_feelings_p == 3 ~ 3,
                     ders_clear_feelings_p == 4 ~ 2,
                     ders_clear_feelings_p == 5 ~ 1)) %>%
  # add column to reverse score "my child knows exactly how he/she is feeling"
  mutate(rev_ders_feelings_know_p = 
           case_when(ders_feelings_know_p == 1 ~ 5,
                     ders_feelings_know_p == 2 ~ 4,
                     ders_feelings_know_p == 3 ~ 3,
                     ders_feelings_know_p == 4 ~ 2,
                     ders_feelings_know_p == 5 ~ 1)) %>%
  # add column to reverse score "when my child is upset, he/she feels like
  # he/she can remain in control of his/her behaviors"
  mutate(rev_ders_upset_behavior_control_p = 
           case_when(ders_upset_behavior_control_p == 1 ~ 5,
                     ders_upset_behavior_control_p == 2 ~ 4,
                     ders_upset_behavior_control_p == 3 ~ 3,
                     ders_upset_behavior_control_p == 4 ~ 2,
                     ders_upset_behavior_control_p == 5 ~ 1)) %>%
  # add column to reverse score "when my child is upset, he/she knows that
  # he/she can find a way to eventually feel better"
  mutate(rev_ders_upset_better_p = 
           case_when(ders_upset_better_p == 1 ~ 5,
                     ders_upset_better_p == 2 ~ 4,
                     ders_upset_better_p == 3 ~ 3,
                     ders_upset_better_p == 4 ~ 2,
                     ders_upset_better_p == 5 ~ 1)) 

ders_cols_to_sum <- dersdata %>%
  select(-all_of(c(
    "src_subject_id",
    "eventname",
    "ders_p_select_language___1",
    "ders_attn_awareness_p",
    "ders_feelings_attentive_p",
    "ders_feelings_care_p",
    "ders_upset_ack_p",
    "ders_clear_feelings_p",
    "ders_feelings_know_p",
    "ders_upset_behavior_control_p",
    "ders_upset_better_p"
  ))) %>%
  colnames()
print(ders_cols_to_sum)
# The following 29 columns should be part of the list to be summed:
# "ders_emotion_overwhelm_p","ders_upset_angry_p","ders_upset_ashamed_p",
# "ders_upset_behavior_p","ders_upset_concentrate_p","ders_upset_control_p",
# "ders_upset_depressed_p","ders_upset_difficulty_p","ders_upset_embarrassed_p",
# "ders_upset_emotion_overwhelm_p","ders_upset_esteem_p",
# "ders_upset_feel_better_p","ders_upset_fixation_p","ders_upset_focus_p",
# "ders_upset_guilty_p","ders_upset_irritation_p",
# "ders_upset_long_time_better_p","ders_upset_lose_control_p",
# "ders_upset_out_control_p","ders_upset_time_p","ders_upset_weak_p",
# "rev_ders_attn_awareness_p","rev_ders_feelings_attentive_p",
# "rev_ders_feelings_care_p","rev_ders_upset_ack_p","rev_ders_clear_feelings_p",
# "rev_ders_feelings_know_p","rev_ders_upset_behavior_control_p","rev_ders_upset_better_p"

dersdata <- dersdata %>%
  # add column to sum across all items (using using reverse-scored versions of
  # eight items above) and make one cumulative score
  mutate(ders_total = rowSums(
    across(all_of(ders_cols_to_sum)))) %>%
  # select only columns relevant to analysis
  select(src_subject_id,eventname,ders_total)


#### See all unique values for each column in DERS-P data ####
# see all unique values (can visually check for NA or errors)
map(dersdata,unique)
# explicitly check for any NA values
which(is.na(dersdata))

#### Provide summary statistics for DERS-P data by data collection year ####
# year 3 follow-up: mean = 55.2, sd = 17.5, n = 9727
# year 4 follow-up: mean = 53.7, sd = 16.8, n = 4498
dersdata %>%
  group_by(eventname) %>%
  summarise(
    mean_ders = mean(ders_total, na.rm = TRUE),
    sd_ders = sd(ders_total, na.rm = TRUE),
    n = n()
  )

### Prepare LES data for analysis ####

### to get number of subjects with LES data in year 3 with NA for at least one
### item other than homelessness and knowing someone who attempted suicide
### (which were only asked in year 4)
mh_y_le %>% group_by(eventname) %>% count()
mh_y_le %>% group_by(eventname) %>% 
  filter(if_any(all_of(c("ple_died_y","ple_injured_y","ple_crime_y",
                         "ple_friend_y","ple_friend_injur_y",
                         "ple_financial_y","ple_sud_y","ple_ill_y",
                         "ple_injur_y","ple_argue_y","ple_job_y",
                         "ple_away_y","ple_arrest_y","ple_friend_died_y",
                         "ple_mh_y","ple_sib_y","ple_victim_y","ple_separ_y",
                         "ple_law_y","ple_school_y","ple_move_y","ple_jail_y",
                         "ple_step_y","ple_new_job_y","ple_new_sib_y",
                         # year 3 items
                         "ple_hit_y", "ple_hospitalized_y", "ple_deported_y", 
                         "ple_foster_care_y", "ple_shot_y", "ple_lockdown_y"
                         # year 4 items (note subjects not asked whether ple_suicide_y was good or bad)
                         # "ple_homeless_y","ple_suicide_y
  )), is.na)) %>% count()


#### Identify and prepare relevant columns ####
ledata <- mh_y_le %>% 
  # remove rows with NA in any of the main items asking about whether event
  # was or was not experienced. Note: items below were asked at all time points 
  # when LES was administered, ie year one, two, three, and four follow-up visits.
  # Before this step, n should be 49151. After this step, n should be 37203.
  filter(!if_any(c(
    "ple_died_y","ple_injured_y","ple_crime_y",
    "ple_friend_y","ple_friend_injur_y",
    "ple_financial_y","ple_sud_y","ple_ill_y",
    "ple_injur_y","ple_argue_y","ple_job_y",
    "ple_away_y","ple_arrest_y","ple_friend_died_y",
    "ple_mh_y","ple_sib_y","ple_victim_y","ple_separ_y",
    "ple_law_y","ple_school_y","ple_move_y","ple_jail_y",
    "ple_step_y","ple_new_job_y","ple_new_sib_y"), is.na)) %>%
  # remove rows corresponding to year three follow-up visits if any of the items
  # added in year three are NA. Before this step, n should be 37203. After this
  # step, n should still be 37203.
  filter(
    !(eventname == "3_year_follow_up_y_arm_1" & 
        # items added in year 3 are below
        if_any(c("ple_hit_y", "ple_hospitalized_y", "ple_deported_y", 
                 "ple_foster_care_y", "ple_shot_y", "ple_lockdown_y"), is.na))) %>%
  # remove rows corresponding to year four follow-up visits if any of the items
  # added in year three or year four are NA. Before this step, n should be 37203. After this
  # step, n should be 37006.
  filter(
    !(eventname == "4_year_follow_up_y_arm_1" & 
        # items added in year 4 are below
        if_any(c("ple_hit_y", "ple_hospitalized_y", "ple_deported_y", 
                 "ple_foster_care_y", "ple_shot_y", "ple_lockdown_y",
                 "ple_homeless_y"
                 # subjects not asked whether suicide was good/bad
                 # "ple_suicide_y"
        ), is.na))) %>%
  # Note: ple_y_ss_total_bad ie total number of events subject reported
  # experiencing as bad does *not* include items added in years 3 and 4. To get
  # total number of events experienced as bad in year 3, need to manually
  # identify the number of these items experienced as bad
  mutate(ple_hit_bad = case_when(is.na(ple_hit_y) ~ NA, #NA
                                 ple_hit_y==0 ~ 0, #not experienced
                                 ple_hit_y==1 & ple_hit_fu_y!=2 ~ 0, #experienced but not bad
                                 ple_hit_y==1 & ple_hit_fu_y==2 ~ 1 #experienced as bad
  )) %>%
  mutate(ple_hospitalized_bad = case_when(is.na(ple_hospitalized_y) ~ NA, #NA
                                          ple_hospitalized_y==0 ~ 0, #not experienced
                                          ple_hospitalized_y==1 & ple_hospitalized_fu_y!=2 ~ 0, #experienced but not bad
                                          ple_hospitalized_y==1 & ple_hospitalized_fu_y==2 ~ 1 #experienced as bad
  )) %>%
  mutate(ple_deported_bad = case_when(is.na(ple_deported_y) ~ NA, #NA
                                      ple_deported_y==0 ~ 0, #not experienced
                                      ple_deported_y==1 & ple_deported_fu_y!=2 ~ 0, #experienced but not bad
                                      ple_deported_y==1 & ple_deported_fu_y==2 ~ 1 #experienced as bad
  )) %>%
  mutate(ple_foster_care_bad = case_when(is.na(ple_foster_care_y) ~ NA, #NA
                                         ple_foster_care_y==0 ~ 0, #not experienced
                                         ple_foster_care_y==1 & ple_foster_care_fu_y!=2 ~ 0, #experienced but not bad
                                         ple_foster_care_y==1 & ple_foster_care_fu_y==2 ~ 1 #experienced as bad
  )) %>%
  mutate(ple_shot_bad = case_when(is.na(ple_shot_y) ~ NA, #NA
                                  ple_shot_y==0 ~ 0, #not experienced
                                  ple_shot_y==1 & ple_shot_fu_y!=2 ~ 0, #experienced but not bad
                                  ple_shot_y==1 & ple_shot_fu_y==2 ~ 1 #experienced as bad
  )) %>%
  mutate(ple_lockdown_bad = case_when(is.na(ple_lockdown_y) ~ NA, #NA
                                      ple_lockdown_y==0 ~ 0, #not experienced
                                      ple_lockdown_y==1 & ple_lockdown_fu_y!=2 ~ 0, #experienced but not bad
                                      ple_lockdown_y==1 & ple_lockdown_fu_y==2 ~ 1 #experienced as bad
  )) %>%
  mutate(ple_homeless_bad = case_when(is.na(ple_homeless_y) ~ NA, #NA
                                      ple_homeless_y==0 ~ 0, #not experienced
                                      ple_homeless_y==1 & ple_homeless_fu_y!=2 ~ 0, #experienced but not bad
                                      ple_homeless_y==1 & ple_homeless_fu_y==2 ~ 1 #experienced as bad
  )) %>%
  # Note: no item asking whether knowing someone who attempting suicide was good or bad
  # mutate(ple_suicide_bad = case_when(is.na(ple_suicide_y) ~ NA, #NA
  #                                ple_suicide_y==0 ~ 0, #not experienced
  #                                ple_suicide_y==1 ~ 1 #experienced (assumed bad)
  #                                )) %>%
  mutate(total_bad_le = case_when(
    # in year 1, total bad life events = ple_y_ss_total_bad
    eventname=="1_year_follow_up_y_arm_1" ~ ple_y_ss_total_bad,
    # in year 2, total bad life events = ple_y_ss_total_bad
    eventname=="2_year_follow_up_y_arm_1" ~ ple_y_ss_total_bad,
    # in year 3, total bad life events = ple_y_ss_total_bad +
    # six additional items
    eventname=="3_year_follow_up_y_arm_1" ~ ple_y_ss_total_bad+
      ple_hit_bad+ple_hospitalized_bad+ple_deported_bad+
      ple_foster_care_bad+ple_shot_bad+ple_lockdown_bad,
    # in year 4, total bad life events = ple_y_ss_total_bad +
    # six additional items added in year 3 and two additional items
    # added in year 4
    eventname=="4_year_follow_up_y_arm_1" ~ ple_y_ss_total_bad+
      ple_hit_bad+ple_hospitalized_bad+ple_deported_bad+
      ple_foster_care_bad+ple_shot_bad+ple_lockdown_bad+
      ple_homeless_bad
    # +ple_suicide_bad
  )) %>%
  # select only columns relevant to analysis ie participant ID, time point,
  # and total number of life events experienced as bad
  select(src_subject_id,eventname,total_bad_le) 


#### Provide summary statistics for LES data by data collection year ####
# year 3 follow-up: mean = 2.47, sd = 2.39, n = 10308
# year 4 follow-up: mean = 2.34, sd = 2.24, n = 4542
ledata %>% 
  group_by(eventname) %>% 
  summarise(
    mean_les = mean(total_bad_le, na.rm = TRUE),
    sd_les = sd(total_bad_le, na.rm = TRUE),
    n = n()
  )

### Combine all data for analysis into one data frame ####
alldata <- genderdata %>%
  # add longitudinal tracking data based on subject ID and data collection year
  left_join(select(abcd_y_lt,
                   c(src_subject_id,eventname,
                     site_id_l,interview_age)),
            by=c("src_subject_id","eventname")) %>%
  # rename age and site columns
  rename(age = interview_age, site = site_id_l) %>%
  # add DERS-P data based on subject ID and data collection year
  left_join(dersdata,by=c("src_subject_id","eventname")) %>%
  # add LES data based on subject ID and data collection year
  left_join(ledata,by=c("src_subject_id","eventname")) %>%
  # add CBCL data based on subject ID and data collection year
  left_join(cbcldata,by=c("src_subject_id","eventname")) %>%
  # add BPM data based on subject ID and data collection year
  left_join(bpmdata,by=c("src_subject_id","eventname")) %>%
  mutate(across(
    c(age, ders_total, total_bad_le, 
      cbcl_int, cbcl_ext,
      bpm_int, bpm_ext),
    # Z-score and center continuous variables
    ~ as.numeric(scale(.,center=TRUE,scale=TRUE)),
    .names = "Z_{.col}"
  )) %>%
  # make genderid, sex, and site factors rather than characters
  mutate(across(c(genderid, gender_details, sex, site), as.factor))

### Get general overview of all data ####
#### Get raw number and percentage for each gender group and data collection year ####
# year 3 follow-up: cis boy = 5113 (51.7%), cis girl = 4462 (45.1%), gd = 308 (3.12%)
# year 4 follow-up: cis boy = 2393 (51.9%), cis girl = 1964 (42.6%), gd = 255 (5.53%)
# note that this does not include individuals who didn't understand or refused item
alldata %>% 
  group_by(eventname) %>% 
  count(genderid) %>% 
  mutate(percentage = n / sum(n) * 100)

#### Create separate data frame for just data from year 4 follow-up visit ####
# n = 4612
yr4data <- alldata %>% filter(eventname=="4_year_follow_up_y_arm_1")

#### Create separate data frame for just data from year 3 follow-up visit ####
# n = 9883
yr3data <- alldata %>% filter(eventname=="3_year_follow_up_y_arm_1")

#### Determine how many subjects switched gender groups between years 3 and 4 ####
gender_change <-
  # start with year 4 data
  yr4data %>%
  # select only genderid and subject id columns
  select(src_subject_id,genderid) %>%
  # rename genderid column to show corresponds to genderid in year 4
  rename(genderid_yr4 = genderid) %>%
  # add genderid from year 3 for subjects present in year 4 data
  left_join(select(yr3data,
                   c(src_subject_id,genderid)),
            by=c("src_subject_id")) %>%
  # rename genderid column to show correspond to genderid in year 3
  rename(genderid_yr3 = genderid) %>%
  # remove subjects with the same gender id in years 3 and 4. before this step
  # n = 4612 and after this step n = 202 so 4410 subjects had a stable gender
  # identity between years 3 and 4
  filter(genderid_yr4 != genderid_yr3) %>%
  group_by(genderid_yr3,genderid_yr4) %>% 
  count()
gender_change
# cis boy -> cis girl = 2, cis boy -> gd = 25, cis girl -> cis boy = 1, 
# cis girl -> gd = 154, gd -> cis boy = 6, gd -> cis girl = 14

### Combine year 4 and year 3 data to have option to use year 3 as x variable ####
#### Make combined year 3 and year 4 data for analysis ####
analysis_data <- yr4data %>%
  # add year 4 specifier
  rename(
    yr4_age = age,
    yr4_total_bad_le = total_bad_le,
    yr4_ders_total = ders_total,
    yr4_cbcl_int = cbcl_int,
    yr4_cbcl_ext = cbcl_ext,
    yr4_bpm_int = bpm_int,
    yr4_bpm_ext = bpm_ext,
    Z_yr4_age = Z_age,
    Z_yr4_total_bad_le = Z_total_bad_le,
    Z_yr4_ders_total = Z_ders_total,
    Z_yr4_cbcl_int = Z_cbcl_int,
    Z_yr4_cbcl_ext = Z_cbcl_ext,
    Z_yr4_bpm_int = Z_bpm_int,
    Z_yr4_bpm_ext = Z_bpm_ext
  ) %>%
  # add year 3 data
  left_join(select(yr3data, 
                   c(src_subject_id,
                     age,
                     total_bad_le,
                     ders_total,
                     cbcl_int,
                     cbcl_ext,
                     bpm_int,
                     bpm_ext,
                     Z_age,
                     Z_total_bad_le,
                     Z_ders_total,
                     Z_cbcl_int,
                     Z_cbcl_ext,
                     Z_bpm_int,
                     Z_bpm_ext
                   )),
            by="src_subject_id") %>%
  # name year 3 data to have year 3 specifier
  rename(
    yr3_age = age,
    yr3_total_bad_le = total_bad_le,
    yr3_ders_total = ders_total,
    yr3_cbcl_int = cbcl_int,
    yr3_cbcl_ext = cbcl_ext,
    yr3_bpm_int = bpm_int,
    yr3_bpm_ext = bpm_ext,
    Z_yr3_age = Z_age,
    Z_yr3_total_bad_le = Z_total_bad_le,
    Z_yr3_ders_total = Z_ders_total,
    Z_yr3_cbcl_int = Z_cbcl_int,
    Z_yr3_cbcl_ext = Z_cbcl_ext,
    Z_yr3_bpm_int = Z_bpm_int,
    Z_yr3_bpm_ext = Z_bpm_ext
  ) %>%
  # remove subjects without LES or DERS in year 3 or without CBCL or BPM in 
  # year 4. before this step, n = 4612 and after this step, n = 3763 
  filter(!is.na(yr3_ders_total), !is.na(yr3_total_bad_le),
         !is.na(yr4_cbcl_int), !is.na(yr4_cbcl_ext),
         !is.na(yr4_bpm_int), !is.na(yr4_bpm_ext))

#### See type of each column ####
str(analysis_data)

#### See all unique values for each column ####
map(analysis_data,unique)

#### Determine how many subjects in each more detailed gender group ####
# cis boy = 1955, cis girl = 1614, nb = 144, trans boy = 33, trans girl = 17
analysis_data %>%
  group_by(gender_details) %>%
  count()

#### Determine how many subjects in each combination of gender and sex group ####
# cis boy and male = 1943, cis boy and female = 2, cis boy and NA sex = 10,
# cis girl and male = 1, cis girl and female = 1604, cis girl and NA = 9,
# gd and male = 22, gd and female = 169, gd and NA sex = 3
analysis_data %>%
  group_by(genderid, sex) %>%
  count()

## BASIC STATS ####

### Get summary stats for each variable ####
sumstats <- 
  analysis_data %>% 
  # group_by(genderid) %>% #group by cis boy, cis girl, or nb and get n for each group
  # group_by(eventname) %>% #get average age at each time point (abstract)
  # group_by(gender_details) %>% #group by cis boy, cis girl, nb, trans boy, or trans girl for supplemental table 1
  # filter(!is.na(sex)) %>% #remove 22 subjects with NA sex (all don't know or refuse) for table 3 column 1
  # group_by(sex) %>% #group by male or female for table 3 columns 2 and 3
  # group_by(sex_details) %>% #group by sex including subjects with don't know or refuse for sex
  summarise(
    n = n(),
    across(
      c("yr4_age",
        # "yr3_age",
        "yr3_total_bad_le","yr3_ders_total",
        "yr4_bpm_int","yr4_cbcl_int","yr4_bpm_ext","yr4_cbcl_ext"),
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
# walk() applies a function, in this case ad.test, to a list
# cat() makes it print the variable name before each test
# note: Shapiro-Wilk test only allows up to 5000 samples so using Anderson-
# Darling test.
# All variables are non-normally distributed.
walk(c("yr3_total_bad_le", 
       "yr3_ders_total", 
       "yr4_bpm_int",
       "yr4_cbcl_int",
       "yr4_bpm_ext",
       "yr4_cbcl_ext"), 
     ~ {
       cat("Variable:", .x, "\n")
       print(ad.test(analysis_data[[.x]]))
     })

### Create basic histogram for each variable ####
# Store histograms in list in case want to view later
variable_histograms <- 
  map(c("yr3_total_bad_le", 
        "yr3_ders_total",  
        "yr4_bpm_int", 
        "yr4_cbcl_int", 
        "yr4_bpm_ext", 
        "yr4_cbcl_ext"), 
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
# All variables are significantly correlated with each other except for age.
corrmat <- 
  analysis_data %>%
  # select relevant columns
  select(c(yr4_age,
           yr3_total_bad_le,
           yr3_ders_total,
           yr4_bpm_int,
           yr4_cbcl_int,
           yr4_bpm_ext,
           yr4_cbcl_ext)) %>% 
  # run correlation tests for all pairs of variables, adjust using fdr
  corr.test(adjust="fdr")
# print correlation matrix, correlation coefficients, and p-values
# note that values are already rounded to the decimal place which is showing
print(corrmat,digits=3)
# See fdr-adjusted p-value for each correlation test
# note that list gives p values above diagonal, going across rows (ie not down
# columns) 
round(corrmat$p.adj,5)

## PLOTS ####

##### Graph of LES vs DERS by gender ####
ggplot(analysis_data, 
       aes(x=yr3_total_bad_le,y=yr3_ders_total, fill=genderid)) +
  # geom_point(aes(color=genderid, shape = genderid),size=2) +
  geom_smooth(aes(linetype=genderid),method="lm",color="black", se=TRUE) +
  scale_linetype_manual(values = c("cis_boy"="31",
                                   "cis_girl"="11",
                                   "gd"="solid")) +
  scale_shape_manual(values=c(21,22,23)) +
  scale_colour_grey(start=0.9,end=0) +
  scale_fill_manual(values=c("grey40","grey85","black")) +
  scale_x_continuous(expand = c(0,0),
                     breaks = seq(0,19, by = 1),
                     limits=c(-.5,19.5)) +
  # scale_y_continuous(expand = c(0,0),
  #                    breaks=seq(20,130,10),
  #                    limits = c(25,131)) +
  guides(
    shape = guide_legend(override.aes = list(size = 3)),
    line = guide_legend(override.aes = list(size = 2))
  ) +
  theme_classic() +
  theme(legend.key.width = unit(0.5, "in"))
# Save graph
# ggsave("les_vs_ders_scatter_bygender.tiff",width=8.5,height=6,unit="in",path="figures")

##### Graph of LES vs DERS by sex ####
# Note: 22 subjects with sex as NA (but with gender data) will not be plotted
ggplot(analysis_data[-which(is.na(analysis_data$sex)),], 
       aes(x=yr3_total_bad_le,y=yr3_ders_total, fill=sex)) +
  # geom_point(aes(color=sex, shape = sex),size=2) +
  geom_smooth(aes(linetype=sex),method="lm",color="black", se=TRUE) +
  scale_linetype_manual(values = c("male"="31",
                                   "female"="11")) +
  scale_shape_manual(values=c(21,22)) +
  scale_colour_grey(start=0.9,end=0) +
  scale_fill_manual(values=c("grey40","grey85")) +
  scale_x_continuous(expand = c(0,0),
                     breaks = seq(0,19, by = 1),
                     limits=c(-.5,19.5)) +
  # scale_y_continuous(expand = c(0,0),
  #                    breaks=seq(20,130,10),
  #                    limits = c(25,131)) +
  guides(
    shape = guide_legend(override.aes = list(size = 3)),
    line = guide_legend(override.aes = list(size = 2))
  ) +
  theme_classic() +
  theme(legend.key.width = unit(0.5, "in"))
# Save graph
# ggsave("les_vs_ders_scatter_bysex.tiff",width=8.5,height=6,unit="in",path="figures")

##### Graphs of LES vs CBCL by gender ####
outcome_list <- c("yr4_cbcl_int","yr4_cbcl_ext")
cbcl_les_plot_list <- list()
for (outcome in outcome_list) {
  # Create plot
  cbcl_les_plot <-
    ggplot(aes(x=yr3_total_bad_le,y=.data[[outcome]],
               fill = genderid
    ),data=analysis_data) +
    # geom_point(aes(color=genderid, shape = genderid), size=2) +
    geom_smooth(aes(linetype=genderid),method="lm",color="black", se=TRUE) +
    scale_linetype_manual(values = c("cis_boy"="31",
                                     "cis_girl"="11",
                                     "gd"="solid")) +
    scale_shape_manual(values=c(21,22,23)) +
    scale_colour_grey(start=0.9,end=0) +
    scale_fill_manual(values=c("grey40","grey85","black")) +
    scale_x_continuous(expand = c(0,0),
                       breaks = seq(0,19, by = 1),
                       limits=c(-.5,19.5)) +
    # scale_y_continuous(expand = c(0,0),
    #                    breaks=seq(0,60,by=10),
    #                    limits=c(0,52)) +
    theme_classic()
  cbcl_les_plot_list[[outcome]] <- cbcl_les_plot
  # Save plot
  # ggsave(paste0("les_vs_",outcome,"_bygender.tiff"),
  # width=8.5,height=6,units = "in",path="figures")
}
cbcl_les_plot_list

##### Graphs of LES vs CBCL by sex ####
# Note: 22 subjects with sex as NA (but with gender data) will not be plotted
outcome_list <- c("yr4_cbcl_int","yr4_cbcl_ext")
cbcl_les_plot_list <- list()
for (outcome in outcome_list) {
  # Create plot
  cbcl_les_plot <-
    ggplot(aes(x=yr3_total_bad_le,y=.data[[outcome]],
               fill = sex
    ),data=analysis_data[-which(is.na(analysis_data$sex)),]) +
    # geom_point(aes(color=sex, shape = sex), size=2) +
    geom_smooth(aes(linetype=sex),method="lm",color="black", se=TRUE) +
    scale_linetype_manual(values = c("male"="31",
                                     "female"="11")) +
    scale_shape_manual(values=c(21,22)) +
    scale_colour_grey(start=0.9,end=0) +
    scale_fill_manual(values=c("grey40","grey85")) +
    scale_x_continuous(expand = c(0,0),
                       breaks = seq(0,19, by = 1),
                       limits=c(-.5,19.5)) +
    # scale_y_continuous(expand = c(0,0),
    #                    breaks=seq(0,60,by=10),
    #                    limits=c(0,52)) +
    theme_classic()
  cbcl_les_plot_list[[outcome]] <- cbcl_les_plot
  # Save plot
  # ggsave(paste0("les_vs_",outcome,"_bysex.tiff"),
  # width=8.5,height=6,units = "in",path="figures")
}
cbcl_les_plot_list

##### Graphs of LES vs BPM by gender ####
outcome_list <- c("yr4_bpm_int","yr4_bpm_ext")
bpm_les_plot_list <- list()
for (outcome in outcome_list) {
  # Create plot
  bpm_les_plot <-
    ggplot(aes(x=yr3_total_bad_le,y=.data[[outcome]],
               fill = genderid
    ),data=analysis_data) +
    # geom_point(aes(color=genderid, shape = genderid),size=2) +
    geom_smooth(aes(linetype=genderid),method="lm",color="black", se=TRUE) +
    scale_linetype_manual(values = c("cis_boy"="31",
                                     "cis_girl"="11",
                                     "gd"="solid")) +
    scale_shape_manual(values=c(21,22,23)) +
    scale_colour_grey(start=0.9,end=0) +
    scale_fill_manual(values=c("grey40","grey85","black")) +
    scale_x_continuous(expand = c(0,0),
                       breaks = seq(0,19, by = 1),
                       limits=c(-.5,19.5)) +
    # scale_y_continuous(expand = c(0,0),
    #                    breaks=seq(0,15,by=5),
    #                    limits=c(0,15)) +
    theme_classic()
  bpm_les_plot_list[[outcome]] <- bpm_les_plot
  # Save plot
  # ggsave(paste0("les_vs_",outcome,"_bygender.tiff"),
  # width=8.5,height=6,units = "in",path="figures")
}
bpm_les_plot_list

##### Graphs of LES vs BPM by sex ####
# Note: 22 subjects with sex as NA (but with gender data) will not be plotted
outcome_list <- c("yr4_bpm_int","yr4_bpm_ext")
bpm_les_plot_list <- list()
for (outcome in outcome_list) {
  # Create plot
  bpm_les_plot <-
    ggplot(aes(x=yr3_total_bad_le,y=.data[[outcome]],
               fill = sex
    ),data=analysis_data[-which(is.na(analysis_data$sex)),]) +
    geom_point(aes(color=sex, shape = sex),size=2) +
    geom_smooth(aes(linetype=sex),method="lm",color="black", se=TRUE) +
    scale_linetype_manual(values = c("male"="31",
                                     "female"="11")) +
    scale_shape_manual(values=c(21,22)) +
    scale_colour_grey(start=0.9,end=0) +
    scale_fill_manual(values=c("grey40","grey85")) +
    scale_x_continuous(expand = c(0,0),
                       breaks = seq(0,19, by = 1),
                       limits=c(-.5,19.5)) +
    # scale_y_continuous(expand = c(0,0),
    #                    breaks=seq(0,15,by=5),
    #                    limits=c(0,15)) +
    theme_classic()
  bpm_les_plot_list[[outcome]] <- bpm_les_plot
  # Save plot
  # ggsave(paste0("les_vs_",outcome,"_bysex.tiff"),
  # width=8.5,height=6,units = "in",path="figures")
}
bpm_les_plot_list

##### Graphs of DERS vs CBCL by gender ####
outcome_list <- c("yr4_cbcl_int","yr4_cbcl_ext")
cbcl_ders_plot_list <- list()
for (outcome in outcome_list) {
  # Create plot
  cbcl_ders_plot <-
    ggplot(aes(x=yr3_ders_total,y=.data[[outcome]],
               linetype=genderid,
               shape = genderid,
               fill = genderid
    ),data=analysis_data) +
    # geom_point(alpha=.6, size = 2) +
    geom_smooth(method="lm",se=TRUE,color="black") +
    scale_linetype_manual(values = c("cis_boy"="31",
                                     "cis_girl"="11",
                                     "gd"="solid")) +
    scale_shape_manual(values=c(21,22,23)) +
    scale_fill_manual(values=c("grey40","grey85","black")) +
    scale_colour_grey(start=0.9,end=0) +
    scale_x_continuous(expand = c(0,0),
                       breaks=seq(20,130,10),
                       limits = c(25,131)) +
    # scale_y_continuous(expand = c(0,0),
    #                    breaks=seq(0,60,10),
    #                    limits = c(0,52)) +
    theme_classic() +
    guides(shape = guide_legend(override.aes = list(size = 2)))
  cbcl_ders_plot_list[[outcome]] <- cbcl_ders_plot
  # Save plot
  # ggsave(paste0("ders_vs_",outcome,"_bygender.tiff"),
  # width=8.3,height=6,units = "in",path="figures")
}
cbcl_ders_plot_list

##### Graphs of DERS vs CBCL by sex ####
# Note: 22 subjects with sex as NA (but with gender data) will not be plotted
outcome_list <- c("yr4_cbcl_int","yr4_cbcl_ext")
cbcl_ders_plot_list <- list()
for (outcome in outcome_list) {
  # Create plot
  cbcl_ders_plot <-
    ggplot(aes(x=yr3_ders_total,y=.data[[outcome]],
               linetype=sex,
               shape = sex,
               fill = sex
    ),data=analysis_data[-which(is.na(analysis_data$sex)),]) +
    # geom_point(alpha=.6, size = 2) +
    geom_smooth(method="lm",se=TRUE,color="black") +
    scale_linetype_manual(values = c("male"="31",
                                     "female"="11")) +
    scale_shape_manual(values=c(21,22,23)) +
    scale_fill_manual(values=c("grey40","grey85")) +
    scale_colour_grey(start=0.9,end=0) +
    scale_x_continuous(expand = c(0,0),
                       breaks=seq(20,130,10),
                       limits = c(25,131)) +
    # scale_y_continuous(expand = c(0,0),
    #                    breaks=seq(0,60,10),
    #                    limits = c(0,52)) +
    theme_classic() +
    guides(shape = guide_legend(override.aes = list(size = 2)))
  cbcl_ders_plot_list[[outcome]] <- cbcl_ders_plot
  # Save plot
  # ggsave(paste0("ders_vs_",outcome,"_bysex.tiff"),
  # width=8.3,height=6,units = "in",path="figures")
}
cbcl_ders_plot_list

##### Graphs of DERS vs BPM by gender ####
outcome_list <- c("yr4_bpm_int","yr4_bpm_ext")
bpm_ders_plot_list <- list()
for (outcome in outcome_list) {
  # Create plot
  bpm_ders_plot <- ggplot(aes(x=yr3_ders_total,y=.data[[outcome]],
                              linetype=genderid,
                              shape = genderid,
                              fill = genderid
  ),data=analysis_data) +
    # geom_point(alpha=.6, size=2) +
    geom_smooth(method="lm",se=TRUE,color="black") +
    scale_linetype_manual(values = c("cis_boy"="31",
                                     "cis_girl"="11",
                                     "gd"="solid")) +
    scale_shape_manual(values=c(21,22,23)) +
    scale_fill_manual(values=c("grey40","grey85","black")) +
    scale_colour_grey(start=0.9,end=0) +
    scale_x_continuous(expand = c(0,0),
                       breaks=seq(20,130,10),
                       limits = c(25,131)) +
    # scale_y_continuous(expand = c(0,0),
    #                    breaks=seq(0,15,by=5),
    #                    limits=c(0,15)) +
    theme_classic() +
    guides(shape = guide_legend(override.aes = list(size = 2)))
  bpm_ders_plot_list[[outcome]] <- bpm_ders_plot
  # Save plot
  # ggsave(paste0("ders_vs_",outcome,"_bygender.tiff"),
  #        width=8.3,height=6,units = "in",path="figures")
}
bpm_ders_plot_list

##### Graphs of DERS vs BPM by sex ####
# Note: 22 subjects with sex as NA (but with gender data) will not be plotted
outcome_list <- c("yr4_bpm_int","yr4_bpm_ext")
bpm_ders_plot_list <- list()
for (outcome in outcome_list) {
  # Create plot
  bpm_ders_plot <- ggplot(aes(x=yr3_ders_total,y=.data[[outcome]],
                              linetype=sex,
                              shape = sex,
                              fill = sex
  ),data=analysis_data[-which(is.na(analysis_data$sex)),]) +
    # geom_point(alpha=.6, size=2) +
    geom_smooth(method="lm",se=TRUE,color="black") +
    scale_linetype_manual(values = c("male"="31",
                                     "female"="11")) +
    scale_shape_manual(values=c(21,22,23)) +
    scale_fill_manual(values=c("grey40","grey85")) +
    scale_colour_grey(start=0.9,end=0) +
    scale_x_continuous(expand = c(0,0),
                       breaks=seq(20,130,10),
                       limits = c(25,131)) +
    # scale_y_continuous(expand = c(0,0),
    #                    breaks=seq(0,15,by=5),
    #                    limits=c(0,15)) +
    theme_classic() +
    guides(shape = guide_legend(override.aes = list(size = 2)))
  bpm_ders_plot_list[[outcome]] <- bpm_ders_plot
  # Save plot
  # ggsave(paste0("ders_vs_",outcome,"_bysex.tiff"),
  #        width=8.3,height=6,units = "in",path="figures")
}
bpm_ders_plot_list

## STEP ONE: BASIC GROUP DIFFERENCES AND REGRESSION ####

### Kruskal-Wallis (non-parametric version of one-way ANOVA) and Dunn test ie ####
### pairwise Wilcoxon tests to determine whether variables differ based on gender

#### Age (year 4)
##### cis boy, cis girl, or nb: chisq=3.874, p=0.144
kruskal.test(yr4_age ~ genderid, data = analysis_data) 
##### cis boy, cis girl, nb, trans boy, or trans girl: chisq=4.510, p=0.341
kruskal.test(yr4_age ~ gender_details, data = analysis_data)

#### LES (year 3) 
##### cis boy, cis girl, or nb: chisq=23.298, p<0.001
kruskal.test(yr3_total_bad_le ~ genderid, data = analysis_data)
dunnTest(analysis_data$yr3_total_bad_le, analysis_data$genderid, method = "bh")
##### cis boy, cis girl, nb, trans boy, or trans girl: chisq=24.474, p<0.001
kruskal.test(yr3_total_bad_le ~ gender_details, data = analysis_data)
dunnTest(analysis_data$yr3_total_bad_le, analysis_data$gender_details, method = "bh")

#### DERS (year 3) 
##### cis boy, cis girl, or nb: chisq=56.241, p<0.001
kruskal.test(yr3_ders_total ~ genderid, data = analysis_data)
dunnTest(analysis_data$yr3_ders_total, analysis_data$genderid, method = "bh")
##### cis boy, cis girl, nb, trans boy, or trans girl: chisq=57.209, p<0.001
kruskal.test(yr3_ders_total ~ gender_details, data = analysis_data)
dunnTest(analysis_data$yr3_ders_total, analysis_data$gender_details, method = "bh")

#### BPM internalizing (year 4) 
##### cis boy, cis girl, or nb: chisq=459.5, p<0.001
kruskal.test(yr4_bpm_int ~ genderid, data = analysis_data)
dunnTest(analysis_data$yr4_bpm_int, analysis_data$genderid, method = "bh")
##### cis boy, cis girl, nb, trans boy, or trans girl: chisq=465.81, p<0.001
kruskal.test(yr4_bpm_int ~ gender_details, data = analysis_data)
dunnTest(analysis_data$yr4_bpm_int, analysis_data$gender_details, method = "bh")

#### CBCL internalizing (year 4) 
##### cis boy, cis girl, or nb: chisq=127.55, p<0.001
kruskal.test(yr4_cbcl_int ~ genderid, data = analysis_data)
dunnTest(analysis_data$yr4_cbcl_int, analysis_data$genderid, method = "bh")
##### cis boy, cis girl, nb, trans boy, or trans girl: chisq=130.26, p<0.001
kruskal.test(yr4_cbcl_int ~ gender_details, data = analysis_data)
dunnTest(analysis_data$yr4_cbcl_int, analysis_data$gender_details, method = "bh")

#### BPM externalizing (year 4)
##### cis boy, cis girl, or nb: chisq=47.97, p<0.001
kruskal.test(yr4_bpm_ext ~ genderid, data = analysis_data)
dunnTest(analysis_data$yr4_bpm_ext, analysis_data$genderid, method = "bh")
##### cis boy, cis girl, nb, trans boy, or trans girl: chisq=48.716, p<0.001
kruskal.test(yr4_bpm_ext ~ gender_details, data = analysis_data)
dunnTest(analysis_data$yr4_bpm_ext, analysis_data$gender_details, method = "bh")

#### CBCL externalizing (year 4)
##### cis boy, cis girl, or nb: chisq=38.625, p<0.001
kruskal.test(yr4_cbcl_ext ~ genderid, data = analysis_data)
dunnTest(analysis_data$yr4_cbcl_ext, analysis_data$genderid, method = "bh")
##### cis boy, cis girl, nb, trans boy, or trans girl: chisq=42.347, p<0.001
kruskal.test(yr4_cbcl_ext ~ gender_details, data = analysis_data)
dunnTest(analysis_data$yr4_cbcl_ext, analysis_data$gender_details, method = "bh")


### Mann-Whitney U (non-parametric version of two-sample t-test) to determine ####
### whether variables differ based on sex 

#### Age (year 4): W=1791327, p=0.158
wilcox.test(yr4_age ~ sex, data = analysis_data)

#### LES (year 3): W=1655890, p=0.006
wilcox.test(yr3_total_bad_le ~ sex, data = analysis_data)

#### DERS (year 3): W=1937138, p<0.001
wilcox.test(yr3_ders_total ~ sex, data = analysis_data)

#### BPM internalizing (year 4): W=1148045, p<0.001
wilcox.test(yr4_bpm_int ~ sex, data = analysis_data)

#### CBCL internalizing (year 4): W=1465864, p=0.029
wilcox.test(yr4_cbcl_int ~ sex, data = analysis_data)

#### BPM externalizing (year 4): W=1594377, p<0.001
wilcox.test(yr4_bpm_ext ~ sex, data = analysis_data)

#### CBCL externalizing (year 4): W=1859323, p=0.025
wilcox.test(yr4_cbcl_ext ~ sex, data = analysis_data)


### Mixed effect linear regression to determine whether DERS differs based ####
### on LES, using age as fixed effect covariate and site as random intercept
#### DERS ~ LES + age + (1|site) ####
ders_les_age_reg <- lmer(Z_yr4_ders_total ~ Z_yr3_total_bad_le + 
                           Z_yr4_age + 
                           (1|site),
                         data=analysis_data)
summary(ders_les_age_reg) # LES but not age are sig associated with DERS
BIC(ders_les_age_reg)
rsq(ders_les_age_reg,adj=TRUE) # full model = 0.017

### Mixed effect linear regression to determine whether CBCL or BPM differ ####
### based on LES and/or DERS, using age as fixed effect covariate and site as random 
### intercept
#### Checking that all pairs of variables are related (with covariates of
#### age and site) to determine whether mediation analysis is reasonable

##### BPM internalizing ~ LES + age + (1|site) ####
bpm_int_les_age_reg <- lmer(Z_yr4_bpm_int ~ Z_yr3_total_bad_le + 
                              Z_yr4_age + 
                              (1|site),
                            data=analysis_data)
summary(bpm_int_les_age_reg) # LES and age are sig associated with BPM int

##### CBCL internalizing ~ LES + age + (1|site) ####
cbcl_int_les_age_reg <- lmer(Z_yr4_cbcl_int ~ Z_yr3_total_bad_le + 
                               Z_yr4_age + 
                               (1|site),
                             data=analysis_data)
summary(cbcl_int_les_age_reg) # LES but not age are sig associated with CBCL int

##### BPM externalizing ~ LES + age + (1|site) ####
bpm_ext_les_age_reg <- lmer(Z_yr4_bpm_ext ~ Z_yr3_total_bad_le + 
                              Z_yr4_age + 
                              (1|site),
                            data=analysis_data)
summary(bpm_ext_les_age_reg) # LES but not age are sig associated with BPM ext

##### CBCL externalizing ~ LES + age + (1|site) ####
cbcl_ext_les_age_reg <- lmer(Z_yr4_cbcl_ext ~ Z_yr3_total_bad_le + 
                               Z_yr4_age + 
                               (1|site),
                             data=analysis_data)
summary(cbcl_ext_les_age_reg) # LES and age are sig associated with CBCL ext

##### BPM internalizing ~ DERS + age + (1|site) ####
bpm_int_ders_age_reg <- lmer(Z_yr4_bpm_int ~ Z_yr3_ders_total + 
                               Z_yr4_age + 
                               (1|site),
                             data=analysis_data)
summary(bpm_int_ders_age_reg) # DERS and age are sig associated with BPM int

##### CBCL internalizing ~ DERS + age + (1|site) ####
cbcl_int_ders_age_reg <- lmer(Z_yr4_cbcl_int ~ Z_yr3_ders_total + 
                                Z_yr4_age + 
                                (1|site),
                              data=analysis_data)
summary(cbcl_int_ders_age_reg) # DERS but not age are sig associated with CBCL int

##### BPM externalizing ~ DERS + age + (1|site) ####
bpm_ext_ders_age_reg <- lmer(Z_yr4_bpm_ext ~ Z_yr3_ders_total + 
                               Z_yr4_age + 
                               (1|site),
                             data=analysis_data)
summary(bpm_ext_ders_age_reg) # DERS but not age are sig associated with BPM ext

##### CBCL externalizing ~ DERS + age + (1|site) ####
cbcl_ext_ders_age_reg <- lmer(Z_yr4_cbcl_ext ~ Z_yr3_ders_total + 
                                Z_yr4_age + 
                                (1|site),
                              data=analysis_data)
summary(cbcl_ext_ders_age_reg) # DERS and age are sig associated with CBCL ext


#### Second include both LES and DERS in same model of psychopathology #### 
#### symptoms to get variance explained for full linear regression model which
#### does not include sex or gender as a covariate

##### BPM internalizing ~ LES + DERS + age + (1|site) ####
bpm_int_les_ders_age_reg <- lmer(Z_yr4_bpm_int ~ Z_yr3_total_bad_le + Z_yr3_ders_total + 
                                   Z_yr4_age + 
                                   (1|site),
                                 data=analysis_data)
summary(bpm_int_les_ders_age_reg) # LES, DERS, and age all sig associated with BPM int
round(confint(bpm_int_les_ders_age_reg),2)
BIC(bpm_int_les_ders_age_reg)
rsq(bpm_int_les_ders_age_reg, adj=TRUE) # full model = 0.052

##### CBCL internalizing ~ LES + DERS + age + (1|site) ####
cbcl_int_les_ders_age_reg <- lmer(Z_yr4_cbcl_int ~ Z_yr3_total_bad_le + Z_yr3_ders_total + 
                                    Z_yr4_age + 
                                    (1|site),
                                  data=analysis_data)
summary(cbcl_int_les_ders_age_reg) # LES and DERS but not age sig associated with CBCL int
round(confint(cbcl_int_les_ders_age_reg),2)
BIC(cbcl_int_les_ders_age_reg)
rsq(cbcl_int_les_ders_age_reg, adj=TRUE) # full model = 0.168

##### BPM externalizing ~ LES + DERS + age + (1|site) ####
bpm_ext_les_ders_age_reg <- lmer(Z_yr4_bpm_ext ~ Z_yr3_total_bad_le + Z_yr3_ders_total + 
                                   Z_yr4_age + 
                                   (1|site),
                                 data=analysis_data)
summary(bpm_ext_les_ders_age_reg) # LES and DERS but not age sig associated with BPM ext
round(confint(bpm_ext_les_ders_age_reg),2)
BIC(bpm_ext_les_ders_age_reg)
rsq(bpm_ext_les_ders_age_reg, adj=TRUE) # full model = 0.055

##### CBCL externalizing ~ LES + DERS + age + (1|site) ####
cbcl_ext_les_ders_age_reg <- lmer(Z_yr4_cbcl_ext ~ Z_yr3_total_bad_le + Z_yr3_ders_total + 
                                    Z_yr4_age + 
                                    (1|site),
                                  data=analysis_data)
summary(cbcl_ext_les_ders_age_reg) # LES, DERS, and age all sig associated with CBCL ext
round(confint(cbcl_ext_les_ders_age_reg),2)
BIC(cbcl_ext_les_ders_age_reg)
rsq(cbcl_ext_les_ders_age_reg, adj=TRUE) # full model = 0.195


## STEP TWO: MODERATING EFFECTS OF GENDER OR SEX ####
### Mixed effect linear regression to determine whether gender moderates ####
### relationship between DERS and LES, use age as fixed effect covariate and
### site as random intercept

#### DERS ~ LES*gender + age + (1|site) ####
ders_les_gendercisboy_reg <- lmer(Z_yr4_ders_total ~ Z_yr3_total_bad_le*genderid_refcisboy + 
                                    Z_yr4_age + 
                                    (1|site),
                                  data=analysis_data)
summary(ders_les_gendercisboy_reg) 
round(confint(ders_les_gendercisboy_reg),2)
BIC(ders_les_gendercisboy_reg)
anova(ders_les_gendercisboy_reg) # sig main effect LES and gender but not age or interaction
rsq(ders_les_gendercisboy_reg,adj=TRUE) # full model = 0.027

### Mixed effect linear regression to determine whether sex moderates ####
### relationship between DERS and LES, use age as fixed effect covariate and
### site as random intercept
#### DERS ~ LES*sex + age + (1|site) ####
ders_les_sex_reg <- lmer(Z_yr4_ders_total ~ Z_yr3_total_bad_le*sex + 
                           Z_yr4_age + 
                           (1|site),
                         data=analysis_data)
summary(ders_les_sex_reg) # sig main effect LES and sex but not age or interaction
round(confint(ders_les_sex_reg),2)
BIC(ders_les_sex_reg)
rsq(ders_les_sex_reg, adj=TRUE) # full model = 0.017

### Mixed effect linear regression to determine whether gender moderates ####
### relationship between LES, DERS, and CBCL or BPM using age as fixed effect
### covariate and site as random intercept

#### BPM internalizing ~ LES*gender + DERS*gender + age + (1|site) ####
bpm_int_les_gendercisboy_reg <- 
  lmer(Z_yr4_bpm_int ~ Z_yr3_total_bad_le*genderid_refcisboy + Z_yr3_ders_total*genderid_refcisboy +
         Z_yr4_age + 
         (1|site),  
       data=analysis_data, REML=FALSE)
summary(bpm_int_les_gendercisboy_reg)
# sig main effect LES, gender, DERS, age, and interaction between LES and gender but not interaction between DERS and gender
round(confint(bpm_int_les_gendercisboy_reg),2)
anova(bpm_int_les_gendercisboy_reg) # sig omnibus interaction between LES and gender
BIC(bpm_int_les_gendercisboy_reg)
rsq(bpm_int_les_gendercisboy_reg,adj=TRUE) # full model = 0.178

# to get betas for gd vs cis girl
bpm_int_les_gendercisgirl_reg <- 
  lmer(Z_yr4_bpm_int ~ Z_yr3_total_bad_le*genderid_refcisgirl + Z_yr3_ders_total*genderid_refcisgirl +
         Z_yr4_age + 
         (1|site),  
       data=analysis_data, REML=FALSE)
summary(bpm_int_les_gendercisgirl_reg)
round(confint(bpm_int_les_gendercisgirl_reg),2)
anova(bpm_int_les_gendercisgirl_reg)
# sig main effect LES, gender, DERS, age, and interaction between LES and gender but not interaction between DERS and gender
BIC(bpm_int_les_gendercisgirl_reg)
rsq(bpm_int_les_gendercisgirl_reg,adj=TRUE) # full model = 0.178

#### CBCL internalizing ~ LES*gender + DERS*gender + age + (1|site) ####
cbcl_int_les_gendercisboy_reg <- 
  lmer(Z_yr4_cbcl_int ~ Z_yr3_total_bad_le*genderid_refcisboy + Z_yr3_ders_total*genderid_refcisboy +
         Z_yr4_age + 
         (1|site),  
       data=analysis_data, REML=FALSE)
summary(cbcl_int_les_gendercisboy_reg)
round(confint(cbcl_int_les_gendercisboy_reg),2)
anova(cbcl_int_les_gendercisboy_reg) 
# sig main effect LES, gender, DERS, and interaction between DERS and gender but not LES and gender or age [BUT NOT SIG GENDER*DERS IF T SCORE]
BIC(cbcl_int_les_gendercisboy_reg)
rsq(cbcl_int_les_gendercisboy_reg,adj=TRUE) # full model = 0.221

# to get betas for gd vs cis girl
cbcl_int_les_gendercisgirl_reg <- 
  lmer(Z_yr4_cbcl_int ~ Z_yr3_total_bad_le*genderid_refcisgirl + Z_yr3_ders_total*genderid_refcisgirl +
         Z_yr4_age + 
         (1|site),  
       data=analysis_data, REML=FALSE)
summary(cbcl_int_les_gendercisgirl_reg)
round(confint(cbcl_int_les_gendercisgirl_reg),2)
anova(cbcl_int_les_gendercisgirl_reg) 
# sig main effect LES, gender, DERS, and interaction between DERS and gender but not LES and gender or age [BUT NOT SIG GENDER*DERS IF T SCORE]
BIC(cbcl_int_les_gendercisgirl_reg)
rsq(cbcl_int_les_gendercisgirl_reg,adj=TRUE) # full model = 0.221

#### BPM externalizing ~ LES*gender + DERS*gender + age + (1|site) ####
bpm_ext_les_gendercisboy_reg <- 
  lmer(Z_yr4_bpm_ext ~ Z_yr3_total_bad_le*genderid_refcisboy + Z_yr3_ders_total*genderid_refcisboy +
         Z_yr4_age + 
         (1|site),  
       data=analysis_data, REML=FALSE)
summary(bpm_ext_les_gendercisboy_reg)
round(confint(bpm_ext_les_gendercisboy_reg),2)
anova(bpm_ext_les_gendercisboy_reg) # sig main effects LES, gender, and DERS but not age or either interaction
BIC(bpm_ext_les_gendercisboy_reg)
rsq(bpm_ext_les_gendercisboy_reg,adj=TRUE) # full model = 0.066

# to get betas for gd vs cis girl
bpm_ext_les_gendercisgirl_reg <- 
  lmer(Z_yr4_bpm_ext ~ Z_yr3_total_bad_le*genderid_refcisgirl + Z_yr3_ders_total*genderid_refcisgirl +
         Z_yr4_age + 
         (1|site),  
       data=analysis_data, REML=FALSE)
summary(bpm_ext_les_gendercisgirl_reg)
round(confint(bpm_ext_les_gendercisgirl_reg),2)
anova(bpm_ext_les_gendercisgirl_reg) # sig main effects LES, gender, and DERS but not age or either interaction
BIC(bpm_ext_les_gendercisgirl_reg)
rsq(bpm_ext_les_gendercisgirl_reg,adj=TRUE) # full model = 0.066

#### CBCL externalizing ~ LES*gender + DERS*gender + age + (1|site) ####
cbcl_ext_les_gendercisboy_reg <- 
  lmer(Z_yr4_cbcl_ext ~ Z_yr3_total_bad_le*genderid_refcisboy + Z_yr3_ders_total*genderid_refcisboy +
         Z_yr4_age + 
         (1|site),  
       data=analysis_data, REML=FALSE)
summary(cbcl_ext_les_gendercisboy_reg)
round(confint(cbcl_ext_les_gendercisboy_reg),2)
anova(cbcl_ext_les_gendercisboy_reg) # sig main effect LES, gender, DERS, and age but neither interaction
BIC(cbcl_ext_les_gendercisboy_reg)
rsq(cbcl_ext_les_gendercisboy_reg,adj=TRUE) # full model = 0.198

# to get betas for gd vs cis girl
cbcl_ext_les_gendercisgirl_reg <- 
  lmer(Z_yr4_cbcl_ext ~ Z_yr3_total_bad_le*genderid_refcisgirl + Z_yr3_ders_total*genderid_refcisgirl +
         Z_yr4_age + 
         (1|site),  
       data=analysis_data, REML=FALSE)
summary(cbcl_ext_les_gendercisgirl_reg)
round(confint(cbcl_ext_les_gendercisgirl_reg),2)
anova(cbcl_ext_les_gendercisgirl_reg) # sig main effect LES, gender, DERS, and age but neither interaction
BIC(cbcl_ext_les_gendercisgirl_reg)
rsq(cbcl_ext_les_gendercisgirl_reg,adj=TRUE) # full model = 0.198

### Mixed effect linear regression to determine whether sex moderates ####
### relationship between LES, DERS, and CBCL or BPM using age as fixed effect
### covariate and site as random intercept

#### BPM internalizing ~ LES*sex + DERS*sex + age + (1|site) ####
bpm_int_les_sex_reg <- 
  lmer(Z_yr4_bpm_int ~ Z_yr3_total_bad_le*sex + Z_yr3_ders_total*sex + 
         Z_yr4_age + 
         (1|site),
       data=analysis_data, REML=FALSE)
summary(bpm_int_les_sex_reg)
# sig main effect LES, sex, DERS, age, and interaction between LES and sex but not interaction between DERS and sex
round(confint(bpm_int_les_sex_reg),2)
anova(bpm_int_les_sex_reg) #sig omnibus test of interaction between LES and sex
BIC(bpm_int_les_sex_reg)
rsq(bpm_int_les_sex_reg,adj=TRUE) # full model = 0.138

#### CBCL internalizing ~ LES*sex + DERS*sex + age + (1|site) ####
cbcl_int_les_sex_reg <- 
  lmer(Z_yr4_cbcl_int ~ Z_yr3_total_bad_le*sex + Z_yr3_ders_total*sex + 
         Z_yr4_age + 
         (1|site),
       data=analysis_data, REML=FALSE)
summary(cbcl_int_les_sex_reg)
anova(cbcl_int_les_sex_reg)
# sig main effect LES, sex, DERS, interaction between LES and sex, interaction between DERS and sex but not age [BUT NEITHER IXN SIG IF T SCORE]
round(confint(cbcl_int_les_sex_reg),2)
BIC(cbcl_int_les_sex_reg)
rsq(cbcl_int_les_sex_reg,adj=TRUE) # full model = 0.200

#### BPM externalizing ~ LES*sex + DERS*sex + age + (1|site) ####
bpm_ext_les_sex_reg <-
  lmer(Z_yr4_bpm_ext ~ Z_yr3_total_bad_le*sex + Z_yr3_ders_total*sex + 
         Z_yr4_age + 
         (1|site),
       data=analysis_data, REML=FALSE)
summary(bpm_ext_les_sex_reg) # sig main effect LES, sex and DERS but not age or either interaction
round(confint(bpm_ext_les_sex_reg),2)
BIC(bpm_ext_les_sex_reg)
rsq(bpm_ext_les_sex_reg,adj=TRUE) # full model = 0.061

#### CBCL externalizing ~ LES*sex + DERS*sex + age + (1|site) ####
cbcl_ext_les_sex_reg <-
  lmer(Z_yr4_cbcl_ext ~ Z_yr3_total_bad_le*sex + Z_yr3_ders_total*sex + 
         Z_yr4_age + 
         (1|site),
       data=analysis_data, REML=FALSE)
summary(cbcl_ext_les_sex_reg) # sig main effect LES, DERS and age but not sex or either interaction
round(confint(cbcl_ext_les_sex_reg),2)
BIC(cbcl_ext_les_sex_reg)
rsq(cbcl_ext_les_sex_reg,adj=TRUE) # full model = 0.196


## STEP THREE: MEDIATING EFFECT OF ER ON CBCL OR BPM ~ LES #### 

### Simple mediation model for bpm internalizing ####
bpmint_model <-
  ' # direct effect
        Z_yr4_bpm_int~ c*Z_yr3_total_bad_le + Z_yr4_age
      # mediator
        Z_yr3_ders_total ~ a*Z_yr3_total_bad_le  + Z_yr4_age
      # indirect effect
        Z_yr4_bpm_int~ b*Z_yr3_ders_total
      # indirect effect (a*b)
        ab := a*b
      # total effect
        total := c + (a*b)'
bpmint_model <- sem(bpmint_model,
                    data = analysis_data,
                    meanstructure = TRUE,
                    se = "robust.cluster",
                    cluster = "site")
summary(bpmint_model, fit.measures=T,
        standardized=F, ci=TRUE, rsquare=TRUE)
parameterEstimates(bpmint_model, boot.ci.type = "bca.simple")
# Estimates and p-values: a = 0.084, p<0.001; b = 0.115, p<0.001; c = 0.187, p<0.001; ab = 0.010, p<0.001

### Simple mediation model for CBCL internalizing ####
cbclint_model <-
  ' # direct effect
      Z_yr4_cbcl_int~ c*Z_yr3_total_bad_le + Z_yr4_age
    # mediator
      Z_yr3_ders_total ~ a*Z_yr3_total_bad_le  + Z_yr4_age
    # indirect effect
      Z_yr4_cbcl_int~ b*Z_yr3_ders_total
    # indirect effect (a*b)
      ab := a*b
    # total effect
      total := c + (a*b)'
cbclint_model <- sem(cbclint_model, 
                     data = analysis_data, 
                     meanstructure = TRUE,
                     se = "robust.cluster",
                     cluster = "site")
summary(cbclint_model, fit.measures=T, 
        standardized=F, ci=TRUE, rsquare=TRUE)
parameterEstimates(cbclint_model, boot.ci.type = "bca.simple")
# Estimates and p-values: a = 0.084, p<0.001; b = 0.373, p<0.001; c = 0.106, p<0.001; ab = 0.031, p<0.001

### Simple mediation model for bpm externalizing ####
bpmext_model <-
  ' # direct effect
        Z_yr4_bpm_ext~ c*Z_yr3_total_bad_le + Z_yr4_age
      # mediator
        Z_yr3_ders_total ~ a*Z_yr3_total_bad_le  + Z_yr4_age
      # indirect effect
        Z_yr4_bpm_ext~ b*Z_yr3_ders_total
      # indirect effect (a*b)
        ab := a*b
      # total effect
        total := c + (a*b)'
bpmext_model <- sem(bpmext_model, 
                    data = analysis_data, 
                    meanstructure = TRUE,
                    se = "robust.cluster",
                    cluster = "site")
summary(bpmext_model, fit.measures=T, 
        standardized=F, ci=TRUE, rsquare=TRUE)
parameterEstimates(bpmext_model, boot.ci.type = "bca.simple")
# Estimates and p-values: a = 0.084, p<0.001; b = 0.160, p<0.001; c = 0.140, p<0.001; ab = 0.014, p<0.001

### Simple mediation model for CBCL externalizing ####
cbclext_model <-
  ' # direct effect
      Z_yr4_cbcl_ext~ c*Z_yr3_total_bad_le + Z_yr4_age
    # mediator
      Z_yr3_ders_total ~ a*Z_yr3_total_bad_le  + Z_yr4_age
    # indirect effect
      Z_yr4_cbcl_ext~ b*Z_yr3_ders_total
    # indirect effect (a*b)
      ab := a*b
    # total effect
      total := c + (a*b)'
cbclext_model <- sem(cbclext_model, 
                     data = analysis_data, 
                     meanstructure = TRUE,
                     se = "robust.cluster",
                     cluster = "site")
summary(cbclext_model, fit.measures=T, 
        standardized=F, ci=TRUE, rsquare=TRUE)
parameterEstimates(cbclext_model, boot.ci.type = "bca.simple")
# Estimates and p-values: a = 0.084, p<0.001; b = 0.383, p<0.001; c = 0.099, p<0.001; ab = 0.032, p<0.001


## STEP FOUR: MODERATING EFFECT OF GENDER OR SEX ON MEDIATION ####

### Moderated mediation model (Hayes model 15) to test whether gender ####
### moderates mediating effect of DERS on relationship between LES and BPM
### internalizing 
bpm_int_gender_model15 <- PROCESS(
  analysis_data,
  y = "Z_yr4_bpm_int",
  x = "Z_yr3_total_bad_le",
  meds = c("Z_yr3_ders_total"),
  mods = c("genderid_refcisboy"),
  covs = c("Z_yr4_age"),
  clusters = "site",
  # hlm.re.m = "site",
  # hlm.re.y = "site",
  # mod.path = c("x-y"),
  mod.path = c("x-y","m-y"),
  cov.path = c("both"),
  nsim = 1000,
  seed = 1234,
  # center = FALSE,
  # std = FALSE,
  digits = 5)
# all interactions: F=5.73, p<0.001
# LES*gender interaction: F=10.64, p<0.001
# DERS*gender interaction: F=0.36, p = 0.695
# c' (interaction sig): cis boy = 0.080, p<0.001, cis girl = 0.219, p<0.001, gd = 0.252, p<0.001
# b (interaction not sig): cis boy = 0.115, p<0.001, cis girl = 0.132, p<0.001, gd = 0.166, p = 0.007
# ab: cis boy = 0.010, p<0.001, cis girl = 0.011, p=0.001, gd = 0.014, p = 0.088

# are the conditional direct effects [c'] of X on Y significant different for
# different groups? Use Z test to find out. 
# Z = (beta1 - beta2 / (sqrt(SE1^2 + SE2^2)))
# cis boy (beta1) vs cis girl (beta2): 
#     Z = (0.07960 - 0.21889)/(sqrt((0.02275^2)+(0.02319)^2)) = -4.287696
# cis boy (beta1) vs gd (beta2):
#     Z = (0.07960 - 0.25246)/(sqrt((0.02275^2)+(0.05919)^2)) = -2.726004
# cis girl (beta1) vs gd (beta2):
#     Z = (0.21889 - 0.25246)/(sqrt((0.02319^2)+(0.05919)^2)) = -0.5280736
# to go from Z score to p-value, find probability of being outside absolute value
# of Z score (because don't know if beta1 is smaller or larger than beta2) and 
# then multiply that by 2 because two-tailed test. Can use default settings of
# mean = 0 and sd = 1 in pnorm function because that is true of Z scores
# cis boy vs cis girl: Z = -4.287696, so pnorm(-abs(-4.287696))*2 = 0.00001805359
# cis boy vs gd: Z = -2.726004, so pnorm(-abs(-2.726004))*2 = 0.006410622
# cis girl vs gd: Z = -0.5280736, so pnorm(-abs(-0.5280736))*2 = 0.5974483
# Finally, we need to fdr correct for multiple tests:
# p.adjust(c(0.00001805359,0.006410622,0.5974483),method="fdr")
# So final p-values rounded to three places are:
# cis boy vs cis girl: p < .001
# cis boy vs gd: p = 0.010
# cis girl vs gd: 0.597

### Moderated mediation model (Hayes model 15) to test whether gender ####
### moderates mediating effect of DERS on relationship between LES and CBCL
### internalizing  for results and figure 4B
cbcl_int_gender_model15 <- PROCESS(
  analysis_data,
  y = "Z_yr4_cbcl_int",
  x = "Z_yr3_total_bad_le",
  meds = c("Z_yr3_ders_total"),
  mods = c("genderid_refcisboy"),
  covs = c("Z_yr4_age"),
  clusters = "site",
  # hlm.re.m = "site",
  # hlm.re.y = "site",
  # mod.path = c("x-y"),
  mod.path = c("x-y","m-y"),
  cov.path = c("both"),
  nsim = 1000,
  seed = 1234,
  # center = FALSE,
  # std = FALSE,
  digits = 5)
# all interactions: F=7.04, p<0.001
# LES*gender interaction: F=2.09, p = 0.124
# DERS*gender interaction: F=10.51, p<0.001
# c' (interaction not sig): cis boy = 0.062, p = 0.003, cis girl = 0.097, p<0.001, gd = 0.173, p = 0.001
# b (interaction sig): cis boy = 0.323, p<0.001, cis girl = 0.413, p<0.001, gd = 0.572, p<0.001
# ab: cis boy = 0.027, p<0.001, cis girl = 0.035, p<0.001, gd = 0.048, p = 0.001

# are the effects [b] of M on Y significant different for
# different groups? Use Z test to find out. 
# Z = (beta1 - beta2 / (sqrt(SE1^2 + SE2^2)))
# cis boy (beta1) vs cis girl (beta2): 
#     Z = (0.32607 - 0.41279)/(sqrt((0.01944^2)+(0.02260)^2)) = -2.90903
# cis boy (beta1) vs gd (beta2):
#     Z = (0.32607 - 0.57162)/(sqrt((0.01944^2)+(0.05697)^2)) = -4.079211
# cis girl (beta1) vs gd (beta2):
#     Z = (0.41279 - 0.57162)/(sqrt((0.02260^2)+(0.05697)^2)) = -2.591493
# to go from Z score to p-value, find probability of being outside absolute value
# of Z score (because don't know if beta1 is smaller or larger than beta2) and 
# then multiply that by 2 because two-tailed test. Can use default settings of
# mean = 0 and sd = 1 in pnorm function because that is true of Z scores
# cis boy vs cis girl: Z = -2.90903, so pnorm(-abs(-2.90903))*2 = 0.003625521
# cis boy vs gd: Z = -4.079211, so pnorm(-abs(-4.079211))*2 = 0.00004518881
# cis girl vs gd: Z = -2.591493, so pnorm(-abs(-2.591493))*2 = 0.009556049
# Finally, we need to fdr correct for multiple tests:
# p.adjust(c(0.003625521,0.00004518881,0.009556049),method="fdr")
# So final p-values rounded to three places are:
# cis boy vs cis girl: p = 0.005
# cis boy vs gd: p < 0.001
# cis girl vs gd: p = 0.010

# are the indirect effects [ab] of X on Y through M significant different for
# different groups? Use Z test to find out. 
# Z = (beta1 - beta2 / (sqrt(SE1^2 + SE2^2)))
# cis boy (beta1) vs cis girl (beta2): 
#     Z = (0.02725 - 0.03450)/(sqrt((0.00649^2)+(0.00842)^2)) = -0.6819728
# cis boy (beta1) vs gd (beta2):
#     Z = (0.02725 - 0.04777)/(sqrt((0.00649^2)+(0.01491)^2)) = -1.261896
# cis girl (beta1) vs gd (beta2):
#     Z = (0.03450 - 0.04777)/(sqrt((0.00842^2)+(0.01491)^2)) = -0.7749713
# to go from Z score to p-value, find probability of being outside absolute value
# of Z score (because don't know if beta1 is smaller or larger than beta2) and 
# then multiply that by 2 because two-tailed test. Can use default settings of
# mean = 0 and sd = 1 in pnorm function because that is true of Z scores
# cis boy vs cis girl: Z = -0.6819728, so pnorm(-abs(-0.6819728))*2 = 0.4952562
# cis boy vs gd: Z = -1.261896, so pnorm(-abs(-1.261896))*2 = 0.2069862
# cis girl vs gd: Z = -0.7749713, so pnorm(-abs(-0.7749713))*2 = 0.4383566
# Finally, we need to fdr correct for multiple tests:
# p.adjust(c(0.4952562,0.2069862,0.4383566),method="fdr")
# So final p-values rounded to three places are:
# cis boy vs cis girl: p = 0.495 (if a is the same and b is sig diff for all groups, why is this not sig diff?)
# cis boy vs gd: p = 0.495 (if a is the same and b is sig diff for all groups, why is this not sig diff?)
# cis girl vs gd: p = 0.495 (if a is the same and b is sig diff for all groups, why is this not sig diff?)

### Moderated mediation model (Hayes model 15) to test whether gender ####
### moderates mediating effect of DERS on relationship between LES and BPM
### externalizing  
bpm_ext_gender_model15 <- PROCESS(
  analysis_data,
  y = "Z_yr4_bpm_ext",
  x = "Z_yr3_total_bad_le",
  meds = c("Z_yr3_ders_total"),
  mods = c("genderid_refcisboy"),
  covs = c("Z_yr4_age"),
  clusters = "site",
  # hlm.re.m = "site",
  # hlm.re.y = "site",
  # mod.path = c("x-y"),
  mod.path = c("x-y","m-y"),
  cov.path = c("both"),
  nsim = 1000,
  seed = 1234,
  # center = FALSE,
  # std = FALSE,
  digits = 5)
# all interactions: F=0.77, p = 0.544
# LES*gender interaction: F=1.43, p = 0.239
# DERS*gender interaction: F=0.09, p = 0.915
# c' (interaction not sig): cis boy = 0.109, p<0.001, cis girl = 0.161, p<0.001, gd = 0.092, p = 0.126
# b (interaction not sig): cis boy = 0.168, p<0.001, cis girl = 0.159, p<0.001, gd = 0.144, p = 0.023
# ab: cis boy = 0.014, p<0.001, cis girl = 0.013, p<0.001, gd = 0.012, p = 0.058

### Moderated mediation model (Hayes model 15) to test whether gender ####
### moderates mediating effect of DERS on relationship between LES and CBCL
### externalizing  
cbcl_ext_gender_model15 <- PROCESS(
  analysis_data,
  y = "Z_yr4_cbcl_ext",
  x = "Z_yr3_total_bad_le",
  meds = c("Z_yr3_ders_total"),
  mods = c("genderid_refcisboy"),
  covs = c("Z_yr4_age"),
  clusters = "site",
  # hlm.re.m = "site",
  # hlm.re.y = "site",
  # mod.path = c("x-y"),
  mod.path = c("x-y","m-y"),
  cov.path = c("both"),
  nsim = 1000,
  seed = 1234,
  # center = FALSE,
  # std = FALSE,
  digits = 5)
# all interactions: F=1.79, p = 0.129
# LES*gender interaction: F=2.24, p = 0.107
# DERS*gender interaction: F=1.79, p = 0.129
# c' (interaction not sig): cis boy = 0.130, p<0.001, cis girl = 0.072, p<0.001, gd = 0.074, p = 0.155
# b (interaction not sig): cis boy = 0.394, p<0.001, cis girl = 0.353, p<0.001, gd = 0.401, p<0.001
# ab: cis boy = 0.033, p<0.001, cis girl = 0.030, p<0.001, gd = 0.034, p<0.001

### Moderated mediation model (Hayes model 15) to test whether sex ####
### moderates mediating effect of DERS on relationship between LES and BPM
### internalizing 
bpm_int_sex_model15 <- PROCESS(
  analysis_data,
  y = "Z_yr4_bpm_int",
  x = "Z_yr3_total_bad_le",
  meds = c("Z_yr3_ders_total"),
  mods = c("sex"),
  covs = c("Z_yr4_age"),
  clusters = "site",
  # hlm.re.m = "site",
  # hlm.re.y = "site",
  # mod.path = c("x-y"),
  mod.path = c("x-y","m-y"),
  cov.path = c("both"),
  nsim = 1000,
  seed = 1234,
  # center = FALSE,
  # std = FALSE,
  digits = 5)
# all interactions: F=14.86, p<0.001
# LES*sex interaction: F=25.71, p<0.001
# DERS*sex interaction: F=2.45, p = 0.117
# c' (interaction sig): female = 0.244, p<0.001, male = 0.081, p<0.001
# b (interaction not sig): female = 0.167, p<0.001, male = 0.117, p<0.001
# ab: female = 0.014, p<0.001, male = 0.010, p<0.001

### Moderated mediation model (Hayes model 15) to test whether sex ####
### moderates mediating effect of DERS on relationship between LES and CBCL
### internalizing  
cbcl_int_sex_model15 <- PROCESS(
  analysis_data,
  y = "Z_yr4_cbcl_int",
  x = "Z_yr3_total_bad_le",
  meds = c("Z_yr3_ders_total"),
  mods = c("sex"),
  covs = c("Z_yr4_age"),
  clusters = "site",
  # hlm.re.m = "site",
  # hlm.re.y = "site",
  # mod.path = c("x-y"),
  mod.path = c("x-y","m-y"),
  cov.path = c("both"),
  nsim = 1000,
  seed = 1234,
  # center = FALSE,
  # std = FALSE,
  digits = 5)
# all interactions: F=12.88, p<0.001
# LES*sex interaction: F=3.77, p = 0.052
# DERS*sex interaction: F=20.31, p<0.001
# c' (interaction not sig): female = 0.117, p<0.001, male = 0.061, p = 0.004
# b (interaction sig): female = 0.456, p<0.001, male = 0.325, p<0.001
# ab: female = 0.039, p<0.001, male = 0.028, p<0.001

### Moderated mediation model (Hayes model 15) to test whether sex ####
### moderates mediating effect of DERS on relationship between LES and BPM
### externalizing 
bpm_ext_sex_model15 <- PROCESS(
  analysis_data,
  y = "Z_yr4_bpm_ext",
  x = "Z_yr3_total_bad_le",
  meds = c("Z_yr3_ders_total"),
  mods = c("sex"),
  covs = c("Z_yr4_age"),
  clusters = "site",
  # hlm.re.m = "site",
  # hlm.re.y = "site",
  # mod.path = c("x-y"),
  mod.path = c("x-y","m-y"),
  cov.path = c("both"),
  nsim = 1000,
  seed = 1234,
  # center = FALSE,
  # std = FALSE,
  digits = 5)
# all interactions: F=1.10, p = 0.334
# LES*sex interaction: F=2.02, p = 0.156
# DERS*sex interaction: F=0.29, p = 0.588
# c' (interaction not sig): female = 0.155, p<0.001, male = 0.110, p<0.001
# b (interaction not sig): female = 0.157, p<0.001, male = 0.174, p<0.001
# ab: female = 0.013, p<0.001, male = 0.015, p<0.001

### Moderated mediation model (Hayes model 15) to test whether sex ####
### moderates mediating effect of DERS on relationship between LES and CBCL
### externalizing  
cbcl_ext_sex_model15 <- PROCESS(
  analysis_data,
  y = "Z_yr4_cbcl_ext",
  x = "Z_yr3_total_bad_le",
  meds = c("Z_yr3_ders_total"),
  mods = c("sex"),
  covs = c("Z_yr4_age"),
  clusters = "site",
  # hlm.re.m = "site",
  # hlm.re.y = "site",
  # mod.path = c("x-y"),
  mod.path = c("x-y","m-y"),
  cov.path = c("both"),
  nsim = 1000,
  seed = 1234,
  # center = FALSE,
  # std = FALSE,
  digits = 5)
# all interactions: F=2.64, p = 0.071
# LES*sex interaction: F=2.92, p = 0.087
# DERS*sex interaction: F=1.92, p = 0.166
# c' (interaction not sig): female = 0.079, p<0.001, male = 0.126, p<0.001
# b (interaction not sig): female = 0.360, p<0.001, male = 0.398, p<0.001
# ab: female = 0.031, p<0.001, male = 0.034, p<0.001
