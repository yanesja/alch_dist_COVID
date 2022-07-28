---
title: "alch_dist_COVID"
output: html_document
---
## Distress Associated with Alcohol Use/Misuse During the COVID-19 Pandemic

Author(s): Yanes, J. A.

The following code loads, prepares, and analyzes data associated with the NIH/NIMH COVID-19 and Mental Health study.

Hypotheses, statistical models, and data handling procedures were preregistered 06/11/2021. Preregistration information provided below.

Yanes, J. A., Atlas, L. Y., Ramchandani, V., & Chung, J. (2021, June 11). Distress and Alcohol Use/Misuse During the COVID 19 Pandemic. Retrieved from osf.io/4asj9

```{r, include=FALSE}
# ==============================
# initial setup
# 1. clear workspace
# 2. libraries
# 3. load workspace
# 4. load raw reference data (ie sample cohort n=3,655)
# 5. load, combine raw baseline data
# 6. load, combine raw repeated measures (RM) data
# 7. load, combine raw end-of-study (ES) data
# 8. cleanup


# ==============================
# 1. clear workspace

# remove vars
rm(list=ls(all.names=TRUE))


# ---------------
# 2. libraries

# name libraries to load (not working)
libs_to_load = c('arm', 'brms', 'bayestestR','corrplot', 
                 'correlation', 'ggplot2', 'ggpubr', 'lme4', 
                 'lmerTest', 'palettetown', 'sjmisc', 
                 'stringr', 'tidyverse', 'sjPlot')

# install and/or load libraries
for (lib in libs_to_load) {
  
  # check if library already installed; if not, install
  if(lib %in% rownames(installed.packages()) == FALSE) {install.packages(lib)}
  
  # load library
  library(lib, character.only=TRUE, quietly=TRUE, verbose=FALSE)
  
}


# cleanup
rm(lib, libs_to_load)

# ---------------
# # 3. load workspace
# 
# # uncomment, run to load data, models, plots, etc.
# load ('alch_dist_COVID.RData')


# ---------------
# 4. load reference data
# reference data (ie n=3,655)
datum0 = read.csv('../data/BL/final-COVID-19-R-baseline.csv',
                     header=TRUE)


# ---------------
# 5. load, combine baseline (BL) data

# BL measures to load
mesrs_to_load = c('AUDIT',
                  'demographics',
                  'DSMXC',
                  'clinical_history')

# loop thru BL measures; load/combine data
# match to SUBJECT_NUMBER in datum0
datum1 = datum0
for (mesr in mesrs_to_load) {
  
  # load BL measure data
  datum_BL = read.csv(paste0('../data/BL/final-', mesr,
                              '-baseline.csv'),
                       header=TRUE)

  # merge w/ datum
  datum1 = merge(datum1, datum_BL,
                      by.x='SUBJECT_NUMBER',
                      by.y='SUBJECT_NUMBER',
                      all.x=TRUE, all.y=FALSE)
}

# cleanup
rm(datum_BL, mesrs_to_load, mesr)


# ---------------
# 6. load, combine repeated measures (RM) data

# RM measures to load
mesrs_to_load = c('DSMXC',
                 'KESSLER5',
                 'COVID-19-R')

# loop thru RM measures; load/combine data
# match to SUBJECT_NUMBER in datum_ref
datum2 = datum1
for (mesr in mesrs_to_load) {
  # create temp DF
  datum_RM = read.csv(paste0('../data/RM/final-', mesr,
                             '-all_time_points.csv'),
                      stringsAsFactors=FALSE,
                      header=TRUE)

  # remove errant rows
  # match 'Week 24' in INTERVAL_NAME
  datum_RM = datum_RM[!grepl('Week 24',
                             datum_RM$INTERVAL_NAME), ]

  # match 'NA' in INTERVAL_NAME
  datum_RM = datum_RM[!grepl('^$',
                             datum_RM$INTERVAL_NAME), ]

  # match NA in INTERVAL_WEEK
  datum_RM = datum_RM[grepl('',
                             datum_RM$VISIT_WEEK), ]

  # reset temp DF, from long-to-wide
  datum_RM_wide = reshape(datum_RM, idvar='SUBJECT_NUMBER',
                          timevar='VISIT_WEEK',
                          direction='wide')

  # merge w/ datum
  datum2 = merge(datum2, datum_RM_wide,
                    by.x='SUBJECT_NUMBER',
                    by.y='SUBJECT_NUMBER',
                    all.x=TRUE, all.y=FALSE)
}

# cleanup
rm(datum_RM, datum_RM_wide, mesr, mesrs_to_load)
rm(tdatum)


# ---------------
# 7. load, combine raw end-of-study (ES) data

# ES measures to load
mesrs_to_load = c('AUDIT')

# loop thru ES measures; load/combine data
# match to SUBJECT_NUMBER in datum_ref
datum3 = datum2
for (mesr in mesrs_to_load) {
  # load BL measure data
  datum_ES = read.csv(paste0('../data/ES/final-', mesr,
                             '-end_of_study.csv'),
                      header=TRUE)

  # merge w/ datum
  datum3 = merge(datum3, datum_ES,
                 by.x='SUBJECT_NUMBER',
                 by.y='SUBJECT_NUMBER',
                 all.x=TRUE, all.y=FALSE)
}

# cleanup
rm(datum_ES, mesrs_to_load, mesr)


# ---------------
# 8. cleanup

# create final dataset, remove not needed datasets
datum = datum3; rm(datum0, datum1, datum2, datum3)


```

```{r, include=FALSE}
# ==============================
# data preparation
# 1. score COVID RM, remove non-responders
# 2. score KESSLER RM, remove non-responders
# 3. score ALCH RM, remove non-responders
# 4. score AUDIT BL
# 5. score AUDIT ES
# 6. compute AUDIT change
# 7. compute 'substance use treatment' (subx) group
# 8. compute 'mental health treatment' (menx) group
# 9. combine race=asian and race=pacific islander


# ==============================
# 1. score COVID RM

# create temp datum, match subject number and distress items
cols = grep('SUBJECT_NUMBER|DISTRESS_PAST_WEEK.', colnames(datum))
tdatum = datum[cols]
tdatum$COVID.countNA = rowSums(is.na(tdatum)) # count NAs
cols = grep('SUBJECT_NUMBER|countNA', colnames(tdatum))
tdatum = tdatum[cols] # remove unneeded cols

# merge w/ datum
datum = merge(datum, tdatum,
               by.x='SUBJECT_NUMBER',
               by.y='SUBJECT_NUMBER',
               all.x=TRUE, all.y=FALSE)

# drop rows w/ less than 2 obs across weeks
datum = subset(datum, COVID.countNA <= 11)

# cleanup
rm(cols, tdatum)


# ---------------
# 2. score KESSLER RM

# create weeks var
weeks_to_score = seq(0,24,2)

# loop thru weeks to score
for (week in weeks_to_score) {
  # determine start/end cols
  start = grep(paste0('NERVOUS.', week,'$'),
               colnames(datum))
  end = grep(paste0('CANT_CHEER_UP.', week,'$'),
             colnames(datum))
  
  # extract data
  tdatum = datum[, start:end]
  
  # compute sum across cols
  datum[, paste0('KESSLER.', week)] = rowSums(tdatum[1:5], 
                                             na.rm=FALSE)
}

# create temp datum, match subject number and distress items
cols = grep('SUBJECT_NUMBER|KESSLER.', colnames(datum))
tdatum = datum[cols]
tdatum$KESSLER.countNA = rowSums(is.na(tdatum)) # count NAs
cols = grep('SUBJECT_NUMBER|countNA', colnames(tdatum))
tdatum = tdatum[cols] # remove unneeded cols

# merge w/ datum
datum = merge(datum, tdatum,
               by.x='SUBJECT_NUMBER',
               by.y='SUBJECT_NUMBER',
               all.x=TRUE, all.y=FALSE)

# drop rows w/ less than 2 obs across weeks
datum = subset(datum, KESSLER.countNA <= 11) # keep <= 11 NAs

# cleanup
rm(start, end, week, weeks_to_score, cols, tdatum)


# _______________
# 3. score ALCH RM

# create temp datum, match subject number and alcohol items
cols = grep('SUBJECT_NUMBER|MORE_4_DRINK_PER_DAY.', colnames(datum))
tdatum = datum[cols]
tdatum$ALCH.countNA = rowSums(is.na(tdatum)) # count NAs
cols = grep('SUBJECT_NUMBER|countNA', colnames(tdatum))
tdatum = tdatum[cols] # remove unneeded cols

# merge w/ datum
datum = merge(datum, tdatum,
               by.x='SUBJECT_NUMBER',
               by.y='SUBJECT_NUMBER',
               all.x=TRUE, all.y=FALSE)

# drop rows w/ less than 2 obs across weeks
datum = subset(datum, ALCH.countNA <= 11) # keep <= 11 NAs

# cleanup
rm(cols, tdatum)


# _______________
# 4. score AUDIT BL

# determine start/end cols
IDs = grep('SUBJECT_NUMBER',
             colnames(datum))
start = grep('AUDIT_ALCOHOL_FREQ_DRINK.x',
             colnames(datum))
end = grep('AUDIT_PPL_CONCERN_DRINK.x',
           colnames(datum))

# extract data
tdatum = datum[, c(IDs, start:end)]

# compute sum across cols
tdatum$AUDIT_BL = rowSums(tdatum[, 2:11], na.rm=TRUE)

# create col w/ AUDIT BL group label
tdatum = tdatum %>% 
  mutate(AUDIT_BL_GROUP = 
           case_when(AUDIT_BL == 0 ~
                       'Abstainer', 
                     AUDIT_BL > 0 & AUDIT_BL <= 7 ~
                       'Low-Risk', 
                     AUDIT_BL > 7 & AUDIT_BL <= 14 ~
                       'High-Risk',
                     AUDIT_BL > 14 ~
                       'AUD',
                     AUDIT_BL == 'NA' ~
                       'NA'))

# remove unneeded cols
cols = grep('SUBJECT_NUMBER|AUDIT_BL|AUDIT_BL_GROUP', 
            colnames(tdatum))
tdatum = tdatum[cols]

# merge w/ datum
datum = merge(datum, tdatum,
               by.x='SUBJECT_NUMBER',
               by.y='SUBJECT_NUMBER',
               all.x=TRUE, all.y=FALSE)

# sanity check 1: chi-square test for independence across 
# AUDIT items. this tests the question: were there 
# associations between NAs and AUDIT items.

# create table: Non-NAs vs. NAs in AUDIT
tdatum = datum[, c(IDs, start:end)]
tdatum_NAs = as.data.frame(is.na(tdatum))
tdatum_NAs = as.data.frame(lapply(tdatum_NAs, table))
tdatum_NAs = tdatum_NAs[seq(2, 20, 2)]
tdatum_NAs = as.data.frame(t(tdatum_NAs))
colnames(tdatum_NAs) = c('Non-NAs', 'NAs')

# chi-square test for independence across AUDIT items.
# this tests the question: were there associations between
# AUDIT items and NA responses.
model.x2.AUDIT = chisq.test(tdatum_NAs) # p < 0.001

# rerun chi-square test from above w/o items 1, 2.
# becuase of wording, it's possible respondents treated
# these itmes differently (ie responded w/ NA to item 2
# with greater frequency).
tdatum_NAs.3to10 = tdatum_NAs[-c(1:2),]
model.x2.AUDIT.8to10 = chisq.test(tdatum_NAs.3to10) # p = 0.44

# visualize residuals from chi-square test
contrib = 100*model.x2.AUDIT$residuals^2/model.x2.AUDIT$statistic
contrib = round(contrib, 2)
png(file='../output/plot_chsqr_AUDIT.png',
    width=300, height=300, units='mm', res=300) # init plot
corrplot(contrib, is.cor = FALSE)
dev.off() # close plot tool

# cleanup
rm(contrib, tdatum_NAs, tdatum_NAs.3to10)

# sanity check 2: correlation between AUDIT (items 1-3)
# and DSMXC alcohol item. this test the association
# between a well-known use indicator (AUDIT, items 1-3) and
# main outcome variable in subsequent analyses (DSMXC, 
# alcohol item).

# extract data
tdatum_small = tdatum[1:4]

# compute sum across cols
tdatum_small$AUDIT_USE = rowSums(tdatum[,2:4], na.rm=TRUE)

# remove unneeded cols
cols = grep('SUBJECT_NUMBER|AUDIT_USE', 
            colnames(tdatum_small))
tdatum_small = tdatum_small[cols]

# merge w/ datum
datum = merge(datum, tdatum_small,
               by.x='SUBJECT_NUMBER',
               by.y='SUBJECT_NUMBER',
               all.x=TRUE, all.y=FALSE)

# correlation analysis 
# AUDIT time=0 (items 1-3) & DSMXC time=0
cor.test(datum$AUDIT_USE, datum$MORE_4_DRINK_PER_DAY.0, 
    use='complete.obs')

# cleanup
rm(cols, end, IDs, start, tdatum, tdatum_small)


# ---------------
# 5. scoring AUDIT ES 

# determine start/end cols
IDs = grep('SUBJECT_NUMBER',
             colnames(datum))
start = grep('AUDIT_ALCOHOL_FREQ_DRINK.y',
             colnames(datum))
end = grep('AUDIT_PPL_CONCERN_DRINK.y',
           colnames(datum))

# extract data
tdatum = datum[, c(IDs, start:end)]

# compute sum across cols
tdatum$AUDIT_ES = rowSums(tdatum[, 2:11], na.rm=TRUE)

# create col w/ ES AUDIT group labels
tdatum = tdatum %>% 
  mutate(AUDIT_ES_GROUP = 
           case_when(AUDIT_ES == 0 ~
                       'Abstainer', 
                     AUDIT_ES > 0 & AUDIT_ES <= 7 ~
                       'Low-Risk', 
                     AUDIT_ES > 7 & AUDIT_ES <= 14 ~
                       'High-Risk',
                     AUDIT_ES > 14 ~
                       'AUD',
                     AUDIT_ES == 'NA' ~
                       'NA'))

# remove unneeded cols
cols = grep('SUBJECT_NUMBER|AUDIT_ES|AUDIT_ES_GROUP', 
            colnames(tdatum))
tdatum = tdatum[cols]

# merge w/ datum
datum = merge(datum, tdatum,
               by.x='SUBJECT_NUMBER',
               by.y='SUBJECT_NUMBER',
               all.x=TRUE, all.y=FALSE)

# cleanup
rm(cols, end, IDs, start, tdatum)


# ---------------
# 6. compute AUDIT change

# extract data
tdatum = datum[, c('SUBJECT_NUMBER', 'AUDIT_BL', 'AUDIT_ES')]

# compute AUDIT CHANGE (AUDIT ES - AUDIT BL)
tdatum$AUDIT_CHANGE = tdatum$AUDIT_ES - tdatum$AUDIT_BL

# create col w/ AUDIT CHANGE group labels
tdatum = tdatum %>% 
  mutate(AUDIT_CHANGE_GROUP = 
           case_when(AUDIT_CHANGE == 0 ~
                       'No Change', 
                     AUDIT_CHANGE > 0  ~
                       'Increase', 
                     AUDIT_CHANGE < 0  ~
                       'Decrease'))

# remove unneeded cols
cols = grep('SUBJECT_NUMBER|AUDIT_CHANGE|AUDIT_CHANGE_GROUP', 
            colnames(tdatum))
tdatum = tdatum[cols]

# merge w/ datum
datum = merge(datum, tdatum,
               by.x='SUBJECT_NUMBER',
               by.y='SUBJECT_NUMBER',
               all.x=TRUE, all.y=FALSE)

# cleanup
rm(cols, tdatum)


# ---------------
# 7. compute 'substance use treatment' (subx) group

# create temp data; use RM subx treatment items
tdatum = datum[, c(grep('SUBJECT_NUMBER|SUBSTANCE_USE_TX.',
             colnames(datum)))]

# find/replace vals
tdatum[tdatum=='3'] = NA

# compute means across cols
cols = grep('SUBSTANCE_USE_TX', colnames(tdatum))
tdatum$SUBX_RM.mean = rowMeans(tdatum[, cols],
                           na.rm = TRUE)

# determine group labels, dummy code
tdatum = tdatum %>% 
  mutate(SUBX_RM = 
           case_when(SUBX_RM.mean == 0 ~ 0,
                     SUBX_RM.mean > 0 ~ 1))

# remove unneeded cols
cols = grep('SUBJECT_NUMBER|SUBX_RM', colnames(tdatum))
tdatum = tdatum[cols]

# merge w/ datum
datum = merge(datum, tdatum,
               by.x='SUBJECT_NUMBER',
               by.y='SUBJECT_NUMBER',
               all.x=TRUE, all.y=FALSE)

# create col w/ better name
datum$SUBX_BL = datum$TREATMENT_ALCOHOL_DRUG

# cleanup
rm(cols, tdatum)


# ---------------
# 8. compute 'mental health treatment' (menx) group

# create temp data; use RM subx treatment items
tdatum = datum[, c(grep('SUBJECT_NUMBER|MENTAL_HLTH_TX.',
             colnames(datum)))]

# find/replace vals
tdatum[tdatum=='3'] = NA

# compute means across cols
cols = grep('MENTAL_HLTH_TX', colnames(tdatum))
tdatum$MENX_RM.mean = rowMeans(tdatum[, cols],
                           na.rm = TRUE)

# determine group labels, dummy code
tdatum = tdatum %>% 
  mutate(MENX_RM = 
           case_when(MENX_RM.mean == 0 ~ 0,
                     MENX_RM.mean > 0 ~ 1))

# remove unneeded cols
cols = grep('SUBJECT_NUMBER|MENX_RM', colnames(tdatum))
tdatum = tdatum[cols]

# merge w/ datum
datum = merge(datum, tdatum,
               by.x='SUBJECT_NUMBER',
               by.y='SUBJECT_NUMBER',
               all.x=TRUE, all.y=FALSE)

# create col w/ better name
datum$MENX_BL = datum$TREATMENT_MENTAL_HEALTH


# ---------------
# 9. combine race=asian and race=pacific islander

# select cols
cols = grep('SUBJECT_NUMBER|RACE', colnames(datum))
tdatum = datum[cols]

# compute col
tdatum$RACE_API = rowSums(tdatum[c('RACE_ASIAN', 
                                   'RACE_PACIFIC_ISLANDER')], 
                          na.rm=TRUE)

# remove unneeded cols
cols = grep('SUBJECT_NUMBER|API', colnames(tdatum))
tdatum = tdatum[cols]

# merge w/ datum
datum = merge(datum, tdatum,
               by.x='SUBJECT_NUMBER',
               by.y='SUBJECT_NUMBER',
               all.x=TRUE, all.y=FALSE)

# cleanup
rm(cols, tdatum)


```

```{r, include=FALSE}
# ==============================
# assemble analysis dataset (ie distress only)
# 1. select/reshape/merg COVID distress data
# 2. select/reshape/merg KESSLER distress data
# 3. merge distress data
# 4. compute multilevel corr coef
# 
# ==============================
# 1. select/reshape/merg COVID distress data

# create temp data
cols = grep('SUBJECT_NUMBER|DISTRESS_PAST_WEEK.[[:digit:]]',
            colnames(datum))
tdatum_COVID = datum[cols]

# reshape data from wide-to-long
ltdatum_COVID = gather(tdatum_COVID,
                      WEEK, 
                      DISTRESS_COVID, 
                      names(tdatum_COVID)[2:14], 
                      factor_key=TRUE)

# remove text from week col
ltdatum_COVID$WEEK = str_remove_all(ltdatum_COVID$WEEK,
                                    'COVID_19_ADULT_DISTRESS_PAST_WEEK.')

# cleanup
rm(cols, tdatum_COVID)


# ---------------
# 2. select/reshape/merg KESSLER distress data

# create temp data
cols = grep('SUBJECT_NUMBER|KESSLER.[[:digit:]]',
            colnames(datum))
tdatum_KESSLER = datum[cols]

# reshape data from wide-to-long
ltdatum_KESSLER = gather(tdatum_KESSLER,
                         WEEK,
                         DISTRESS_KESSLER,
                         names(tdatum_KESSLER)[2:14],
                         factor_key=TRUE)

# remove text from week col
ltdatum_KESSLER$WEEK = str_remove_all(ltdatum_KESSLER$WEEK,
                                      'KESSLER.')

# cleanup
rm(cols, tdatum_KESSLER)


# ---------------
# 3. merge distress data

# merge COVID distress data w/ KESSLER distress data
ldatum = merge(ltdatum_KESSLER, ltdatum_COVID,
               by=c('SUBJECT_NUMBER', 'WEEK'))

# cleanup
rm(ltdatum_KESSLER, ltdatum_COVID)
```


## multilevel correlation
```{r}
# vars: KESSLER distress, GENERAL distress
correlation(ldatum, partial=FALSE, multilevel=TRUE)

# bivar scatter: KESSLER distress, GENERAL distress
# construct main plot
png(file='../output/plt_scttr_cdist_x_kdist.png',
    width=100, height=100, units='mm', res=300)
ggplot(ldatum, aes(x=DISTRESS_COVID, y=DISTRESS_KESSLER)) +
  geom_jitter(height=.25, width=.25,
              colour='gray', alpha=.5) + 
  geom_smooth(method=lm, size=.66, color='black',
              level=.99, formula=y~x, se=TRUE) +
  labs(x='General Distress',
       y='Psychological Distress') +
    theme_classic()
dev.off() # close plot

```

```{r, include=FALSE}
# ==============================
# assemble analysis dataset (ie remaining measures)
# 1. select/merge BL/ES cols to ldatum
# 2. select/merge RM cols to ldatum

# ==============================
# 1. select/merge BL/ES cols to ldatum

# select BL/ES cols
tdatum = datum[, c('AUDIT_BL',
                   'AUDIT_BL_GROUP',
                   'SEX_MALE',
                   'INCOME',
                   'ETHNICITY_NOTLATINO',
                   'ETHNICITY_LATINO',
                   'RACE_AMERICAN_INDIAN',
                   'RACE_API',
                   'RACE_AFRICAN_BLACK',
                   'RACE_CAUCASIAN',
                   'RACE_MULTIPLE',
                   'SUBX_RM',
                   'SUBX_BL',
                   'MENX_BL',
                   'MENX_RM',
                   'SUBJECT_NUMBER')]

# merge BL/ES cols to ldatum
ldatum = merge(ldatum, tdatum,
               by.x='SUBJECT_NUMBER',
               by.y='SUBJECT_NUMBER',
               all.x = TRUE, all.y = FALSE)

# cleanup
rm(tdatum)


# ---------------
# 2. select/reshape RM vars, merge w/ ldatum

# name vars 
vars = c('MORE_4_DRINK_PER_DAY',
         'TOBACCO', 
         'SELF_MEDICATED',
         'VISIT_ABSOLUTE_DAY')

# loop thru vars
for (var in vars) {
  
  # select cols
  cols = grep(paste0('SUBJECT_NUMBER|', var, '.[[:digit:]]*$'),
            colnames(datum))
  tdatum = datum[cols]
  
  # reshape from wide-to-long
  ltdatum = gather(tdatum,
                   WEEK,
                   !!var,
                   names(tdatum)[2:14],
                   factor_key=TRUE)
  
  # remove text from week col
  ltdatum$WEEK = str_remove_all(ltdatum$WEEK,
                                paste0(var, '.'))
  # merge w/ ldatum
  ldatum = merge(ldatum, ltdatum,
                 by=c('SUBJECT_NUMBER', 'WEEK'))
}

# create cols w/ better names
ldatum$ALCH = ldatum$MORE_4_DRINK_PER_DAY
ldatum$TBCO = ldatum$TOBACCO
ldatum$DRGS_OTHR = ldatum$SELF_MEDICATED
ldatum$DAY = ldatum$VISIT_ABSOLUTE_DAY

# cleanup
rm(cols, var, vars, tdatum, ltdatum)
```

```{r, include=FALSE}
# ==============================
# centering
# 1. grand-mean center BL/ES vars
# 2. grand-mean center & cluster-mean center RM vars


# ==============================
# 1. grand-mean center BL/ES vars
# creates between-subs var

# name IV cols to center
# name vars to center
cols_to_cntr = c('AUDIT_BL', 
                 'AUDIT_CHANGE',
                 'ETHNICITY_NOTLATINO',
                 'ETHNICITY_LATINO',
                 'INCOME',
                 'RACE_AFRICAN_BLACK',
                 'RACE_AMERICAN_INDIAN',
                 'RACE_API',
                 'RACE_CAUCASIAN',
                 'RACE_MULTIPLE',
                 'SUBX_RM',
                 'SUBX_BL',
                 'MENX_BL',
                 'MENX_RM',
                 'SEX_MALE')

# loop thru cols
for (col in cols_to_cntr) {

  # compute grand mean, subtract from obs, create new col
  ldatum[, paste0(col, '.GMC')] = ldatum[[col]] - 
    mean(ldatum[[col]], na.rm = TRUE)
}

# cleanup
rm(col, cols_to_cntr)


# ---------------
# 2. grand-mean center & cluster-mean center RM vars
# creates between-subs var & within-subs var

# 2A. cluster mean centering (CMC) RM KESSLER DISTRESS
ldatum = ldatum %>%
  # by (1) grand mean centering
  mutate(DISTRESS_KESSLER.GMC = DISTRESS_KESSLER - 
           mean(DISTRESS_KESSLER, na.rm = TRUE)) %>% # GMC = grand mean center
  # by (2) cluster mean centering
  group_by(SUBJECT_NUMBER) %>%
  mutate(DISTRESS_KESSLER.CM = mean(DISTRESS_KESSLER,
                                    na.rm = TRUE), # CM = cluster mean
         DISTRESS_KESSLER.CMC = DISTRESS_KESSLER - 
           DISTRESS_KESSLER.CM) %>% # CMC = cluster mean center
  ungroup %>%
  # grand mean center of cluster means
  mutate(DISTRESS_KESSLER.GMCCM = DISTRESS_KESSLER.CM -
           mean(DISTRESS_KESSLER.CM, na.rm = TRUE))

# 2B. cluster mean centering (CMC) RM COVID DISTRESS
ldatum = ldatum %>%
  # by (1) grand mean centering (GMC)
  mutate(DISTRESS_COVID.GMC = DISTRESS_COVID - 
           mean(DISTRESS_COVID, na.rm=TRUE)) %>% 
  # by (2) cluster mean centering (CMC)
  group_by(SUBJECT_NUMBER) %>%
  mutate(DISTRESS_COVID.CM = mean(DISTRESS_COVID, 
                                  na.rm=TRUE),
         DISTRESS_COVID.CMC = DISTRESS_COVID - 
           DISTRESS_COVID.CM) %>%
  ungroup %>%
  # grand mean center of cluster means (GMCCM)
  mutate(DISTRESS_COVID.GMCCM = DISTRESS_COVID.CM - 
           mean(DISTRESS_COVID.CM, na.rm=TRUE))

# 2C. cluster mean centering (CMC) RM VISIT ABSOLUTE STUDY DAY
ldatum = ldatum %>%
  # by (1) grand mean centering
  mutate(DAY.GMC = DAY - 
           mean(DAY, na.rm = TRUE)) %>% # GMC = grand mean center
  # by (2) cluster mean centering
  group_by(SUBJECT_NUMBER) %>%
  mutate(DAY.CM = mean(DAY,
                       na.rm = TRUE), # CM = cluster mean
         DAY.CMC = DAY - 
           DAY.CM) %>% # CMC = cluster mean center
  ungroup %>%
  # grand mean center of cluster means
  mutate(DAY.GMCCM = DAY.CM -
           mean(DAY.CM, na.rm = TRUE)) 

```

## modeling
_1. hypo cdist1: ALCH ~ COVID * BL AUDIT_                  \
_2. hypo kdist1: ALCH ~ KESSLER * BL AUDIT_                \
_3. hypo cdist2: ALCH ~ COVID * BL AUDIT * SEX_            \
_4. hypo kdist2: ALCH ~ KESSLER * BL AUDIT * SEX_          \
_5. hypo cdist3: ALCH ~ COVID * BL AUDIT * INCOME_         \
_6. hypo kdist3: ALCH ~ KESSLER * BL AUDIT * INCOME_       \
_7. hypo cdist4: ALCH ~ COVID * BL AUDIT * ETHNICITY_      \
_8. hypo kdist4: ALCH ~ KESSLER * BL AUDIT * ETHNICITY_    \
_9. hypo cdist5: ALCH ~ COVID * BL AUDIT * RACE_           \
_10. hypo cdist5: ALCH ~ KESSLER * BL AUDIT * RACE_        \
_11. hypo cdist7: ALCH ~ COVID * SUBX_                     \
_12. hypo kdist7: ALCH ~ KESSLER * SUBX_                   \
_13. hypo cdist6: ALCH ~ COVID * MENX_                     \
_14. hypo kdist6: ALCH ~ KESSLER * MENX_                   \
_15. hypo cdist8: TBCO ~ COVID_                            \
_16. hypo kdist8: TBCO ~ KESSLER_                          \
_17. hypo cdist9: DRGS_OTHR ~ COVID_                       \
_18. hypo kdist9: DRGS_OTHR ~ KESSLER_                     \
_19. hypo cdist10: ALCH ~ COVID * BL AUDIT * DAY_          \
_20. hypo kdist10: ALCH ~ KESSLER * BL AUDIT * DAY_        \
\
\
_1. hypo cdist1: ALCH ~ COVID * BL AUDIT_
```{r, message=FALSE}
# run/report freq model
model.cdist1 = lmer(ALCH ~ DISTRESS_COVID.CMC *
                      DISTRESS_COVID.GMCCM *
                      AUDIT_BL.GMC +
                      (DISTRESS_COVID.CMC|SUBJECT_NUMBER),
                    data=ldatum)
sjPlot::tab_model(model.cdist1) # model ok

# # run bayes model
# begin = paste('model.cdist1.bayes begin =', Sys.time())
# model.cdist1.bayes = brm(ALCH ~ DISTRESS_COVID.CMC *
#                            DISTRESS_COVID.GMCCM *
#                            AUDIT_BL.GMC +
#                            (DISTRESS_COVID.CMC|SUBJECT_NUMBER),
#                          data=ldatum,
#                          family=gaussian(),
#                          iter=2000,
#                          chains=4,
#                          cores=4)
# end = paste('model.cdist1.bayes end =', Sys.time())
# print(begin)
# print(end)

# # report bayes model
# describe_posterior(model.cdist1.bayes,
#                    effects='fixed', ci=0.95,
#                    component='all',
#                    test=c('p_direction', 'p_significance', 'rope'),
#                    rope_range='default',
#                    rope_ci=1,
#                    centrality='all')

# run/report sub-model (abstainers only)
# this is for visualizations only
model.cdist1.AB = lmer(ALCH ~ DISTRESS_COVID.CMC *
                         DISTRESS_COVID.GMCCM +
                         (DISTRESS_COVID.CMC|SUBJECT_NUMBER),
                       data=subset(ldatum, AUDIT_BL_GROUP == 'Abstainer'))
sjPlot::tab_model(model.cdist1.AB) # model ok

# run/report submodel (low-risk only)
# this is for visualizations only
model.cdist1.LR = lmer(ALCH ~ DISTRESS_COVID.CMC *
                         DISTRESS_COVID.GMCCM +
                         (DISTRESS_COVID.CMC|SUBJECT_NUMBER),
                       data=subset(ldatum, AUDIT_BL_GROUP == 'Low-Risk'))
sjPlot::tab_model(model.cdist1.LR) # model ok

# run/report submodel (high-risk only)
# this is for visualizations only
model.cdist1.HR = lmer(ALCH ~ DISTRESS_COVID.CMC *
                         DISTRESS_COVID.GMCCM +
                         (DISTRESS_COVID.CMC|SUBJECT_NUMBER),
                       data=subset(ldatum, AUDIT_BL_GROUP == 'High-Risk'))
sjPlot::tab_model(model.cdist1.HR) # model ok

# run/report submodel (AUD only)
# this is for visualizations only
model.cdist1.AUD = lmer(ALCH ~ DISTRESS_COVID.CMC *
                          DISTRESS_COVID.GMCCM +
                          (DISTRESS_COVID.CMC|SUBJECT_NUMBER),
                        data=subset(ldatum, AUDIT_BL_GROUP == 'AUD'))
sjPlot::tab_model(model.cdist1.AUD) # model ok


```
\
\
_2. hypo kdist1: ALCH ~ KESSLER * BL AUDIT_
```{r, message=FALSE}
# run/report model
model.kdist1.A = lmer(ALCH ~ DISTRESS_KESSLER.CMC *
                        DISTRESS_KESSLER.GMCCM *
                        AUDIT_BL.GMC +
                        (DISTRESS_KESSLER.CMC|SUBJECT_NUMBER),
                      data=ldatum) # model conv fail

# rerun/rereport model w/ different optimizer if necessary
# https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html
model.kdist1.A_sc = update(model.kdist1.A, data=ldatum)
ss = getME(model.kdist1.A_sc, c("theta","fixef"))
model.kdist1.B = update(model.kdist1.A,
                        start=ss,
                        control=lmerControl(optimizer="bobyqa",
                                            optCtrl=list(maxfun=2e5)))
sjPlot::tab_model(model.kdist1.B) # model ok

# # run bayes model
# begin = paste('model.kdist1.bayes begin =', Sys.time())
# model.kdist1.bayes = brm(ALCH ~ DISTRESS_KESSLER.CMC *
#                            DISTRESS_KESSLER.GMCCM *
#                            AUDIT_BL.GMC +
#                            (DISTRESS_KESSLER.CMC|SUBJECT_NUMBER),
#                          data=ldatum,
#                          family=gaussian(),
#                          iter=2000,
#                          chains=4,
#                          cores=4)
# end = paste('model.kdist1.bayes end =', Sys.time())
# print(begin)
# print(end)

# # report bayes model
# describe_posterior(model.kdist1.bayes,
#                    effects='fixed', ci=0.95,
#                    component='all',
#                    test=c('p_direction', 'p_significance', 'rope'),
#                    rope_range='default',
#                    rope_ci=1,
#                    centrality='all')

# run/report model (abstainers only)
model.kdist1.AB = lmer(ALCH ~ DISTRESS_KESSLER.CMC *
                         DISTRESS_KESSLER.GMCCM +
                         (DISTRESS_KESSLER.CMC|SUBJECT_NUMBER),
                       data=subset(ldatum, AUDIT_BL_GROUP == 'Abstainer'))
sjPlot::tab_model(model.kdist1.AB) # model ok

# run/report model (low-risk only)
model.kdist1.LR = lmer(ALCH ~ DISTRESS_KESSLER.CMC *
                         DISTRESS_KESSLER.GMCCM +
                         (DISTRESS_KESSLER.CMC|SUBJECT_NUMBER),
                       data=subset(ldatum, AUDIT_BL_GROUP == 'Low-Risk'))
sjPlot::tab_model(model.kdist1.LR) # model ok

# run/report model (high-risk only)
model.kdist1.HR.A = lmer(ALCH ~ DISTRESS_KESSLER.CMC *
                         DISTRESS_KESSLER.GMCCM +
                         (DISTRESS_KESSLER.CMC|SUBJECT_NUMBER),
                       data=subset(ldatum, AUDIT_BL_GROUP == 'High-Risk')) # model conv fail

# rerun/rereport model w/ different optimizer if necessary
# https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html
model.kdist1.HR.A_sc = update(model.kdist1.HR.A, data=ldatum)
ss = getME(model.kdist1.HR.A_sc, c("theta","fixef"))
model.kdist1.HR.B = update(model.kdist1.HR.A,
                        start=ss,
                        control=lmerControl(optimizer="bobyqa",
                                            optCtrl=list(maxfun=2e5)))
sjPlot::tab_model(model.kdist1.HR.B) # model ok

# run/report model (AUD only)
model.kdist1.AUD = lmer(ALCH ~ DISTRESS_KESSLER.CMC *
                      DISTRESS_KESSLER.GMCCM +
                      (DISTRESS_KESSLER.CMC|SUBJECT_NUMBER),
                    data=subset(ldatum, AUDIT_BL_GROUP == 'AUD'))
sjPlot::tab_model(model.kdist1.AUD) # model ok

```
\
\
_3. hypo cdist2: ALCH ~ COVID BL AUDIT * SEX_
```{r, message=FALSE}
model.cdist2 = lmer(ALCH ~ DISTRESS_COVID.CMC *
                      DISTRESS_COVID.GMCCM *
                      AUDIT_BL.GMC *
                      SEX_MALE.GMC +
                      (DISTRESS_COVID.CMC|SUBJECT_NUMBER),
                    data=ldatum)
sjPlot::tab_model(model.cdist2) # model ok


```
\
\
_4. hypo kdist2: ALCH ~ KESSLER * BL AUDIT * SEX_
```{r, message=FALSE}
# run/report model
model.kdist2 = lmer(ALCH ~ DISTRESS_KESSLER.CMC *
                      DISTRESS_KESSLER.GMCCM *
                      AUDIT_BL.GMC *
                      SEX_MALE.GMC +
                      (DISTRESS_KESSLER.CMC|SUBJECT_NUMBER),
                    data=ldatum)
sjPlot::tab_model(model.kdist2) # model ok


```
\
\
_5. hypo cdist3: ALCH ~ COVID * BL AUDIT * INCOME_
```{r, message=FALSE}
model.cdist3.A = lmer(ALCH ~ DISTRESS_COVID.CMC *
                        DISTRESS_COVID.GMCCM *
                        AUDIT_BL.GMC *
                        INCOME.GMC +
                        (DISTRESS_COVID.CMC|SUBJECT_NUMBER),
                      data=ldatum) # model conv fail

# rerun/rereport model w/ different optimizer if necessary
# https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html
model.cdist3.A_sc = update(model.cdist3.A, data=ldatum)
ss = getME(model.cdist3.A_sc, c("theta","fixef"))
model.cdist3.B = update(model.cdist3.A,
                        start=ss,
                        control=lmerControl(optimizer="bobyqa",
                                            optCtrl=list(maxfun=2e5)))
sjPlot::tab_model(model.cdist3.B) # model ok


```
\
\
_6. _hypo kdist3: ALCH ~ KESSLER * BL AUDIT * INCOME_
```{r, message=FALSE}
model.kdist3 = lmer(ALCH ~ DISTRESS_KESSLER.CMC *
                      DISTRESS_KESSLER.GMCCM *
                      AUDIT_BL.GMC *
                      INCOME.GMC +
                      (DISTRESS_KESSLER.CMC|SUBJECT_NUMBER),
                    data=ldatum) # model ok
sjPlot::tab_model(model.kdist3)


```
\
\
_7. hypo cdist4: ALCH ~ COVID * BL AUDIT * ETHNICITY_
```{r, message=FALSE}
model.cdist4 = lmer(ALCH ~ DISTRESS_COVID.CMC *
                      DISTRESS_COVID.GMCCM *
                      AUDIT_BL.GMC *
                      ETHNICITY_LATINO.GMC +
                      (DISTRESS_COVID.CMC|SUBJECT_NUMBER),
                    data=ldatum)
sjPlot::tab_model(model.cdist4) # model ok


```
\
\
_8. hypo kdist4: ALCH ~ KESSLER * BL AUDIT * ETHNICITY_
```{r, message=FALSE}
model.kdist4.A = lmer(ALCH ~ DISTRESS_KESSLER.CMC *
                      DISTRESS_KESSLER.GMCCM *
                      AUDIT_BL.GMC *
                      ETHNICITY_LATINO.GMC +
                      (DISTRESS_KESSLER.CMC|SUBJECT_NUMBER),
                    data=ldatum) # model conv fail

# rerun/rereport model w/ different optimizer if necessary
# https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html
model.kdist4.A_sc = update(model.kdist4.A, data=ldatum)
ss = getME(model.kdist4.A_sc, c("theta","fixef"))
model.kdist4.B = update(model.kdist4.A,
                        start=ss,
                        control=lmerControl(optimizer="bobyqa",
                                            optCtrl=list(maxfun=2e5)))


```
\
\
_9. hypo cdist5: ALCH ~ COVID * BL AUDIT * RACE_
```{r, message=FALSE}
# race: asian/pacific island (A), black (B), american indian (I), white (W)
model.cdist5 = lmer(ALCH ~ DISTRESS_COVID.CMC *
                      DISTRESS_COVID.GMCCM *
                      AUDIT_BL.GMC +
                      RACE_API.GMC +
                      DISTRESS_COVID.CMC * RACE_API.GMC +
                      DISTRESS_COVID.GMCCM * RACE_API.GMC + 
                      AUDIT_BL.GMC * RACE_API.GMC + 
                      (DISTRESS_COVID.CMC * DISTRESS_COVID.GMCCM) * RACE_API.GMC +
                      (DISTRESS_COVID.CMC * AUDIT_BL.GMC) * RACE_API.GMC + 
                      (DISTRESS_COVID.GMCCM * AUDIT_BL.GMC) * RACE_API.GMC +
                      (DISTRESS_COVID.CMC * DISTRESS_COVID.GMCCM * AUDIT_BL.GMC) * RACE_API.GMC +
                      RACE_AFRICAN_BLACK.GMC +
                      DISTRESS_COVID.CMC * RACE_AFRICAN_BLACK.GMC +
                      DISTRESS_COVID.GMCCM * RACE_AFRICAN_BLACK.GMC + 
                      AUDIT_BL.GMC * RACE_AFRICAN_BLACK.GMC + 
                      (DISTRESS_COVID.CMC * DISTRESS_COVID.GMCCM) * RACE_AFRICAN_BLACK.GMC +
                      (DISTRESS_COVID.CMC * AUDIT_BL.GMC) * RACE_AFRICAN_BLACK.GMC + 
                      (DISTRESS_COVID.GMCCM * AUDIT_BL.GMC) * RACE_AFRICAN_BLACK.GMC +
                      (DISTRESS_COVID.CMC * DISTRESS_COVID.GMCCM * AUDIT_BL.GMC) * RACE_AFRICAN_BLACK.GMC +
                      RACE_AMERICAN_INDIAN.GMC +
                      DISTRESS_COVID.CMC * RACE_AMERICAN_INDIAN.GMC +
                      DISTRESS_COVID.GMCCM * RACE_AMERICAN_INDIAN.GMC + 
                      AUDIT_BL.GMC * RACE_AMERICAN_INDIAN.GMC + 
                      (DISTRESS_COVID.CMC * DISTRESS_COVID.GMCCM) * RACE_AMERICAN_INDIAN.GMC +
                      (DISTRESS_COVID.CMC * AUDIT_BL.GMC) * RACE_AMERICAN_INDIAN.GMC + 
                      (DISTRESS_COVID.GMCCM * AUDIT_BL.GMC) * RACE_AMERICAN_INDIAN.GMC +
                      (DISTRESS_COVID.CMC * DISTRESS_COVID.GMCCM * AUDIT_BL.GMC) * RACE_AMERICAN_INDIAN.GMC +
                      RACE_CAUCASIAN.GMC +
                      DISTRESS_COVID.CMC * RACE_CAUCASIAN.GMC +
                      DISTRESS_COVID.GMCCM * RACE_CAUCASIAN.GMC + 
                      AUDIT_BL.GMC * RACE_CAUCASIAN.GMC + 
                      (DISTRESS_COVID.CMC * DISTRESS_COVID.GMCCM) * RACE_CAUCASIAN.GMC +
                      (DISTRESS_COVID.CMC * AUDIT_BL.GMC) * RACE_CAUCASIAN.GMC + 
                      (DISTRESS_COVID.GMCCM * AUDIT_BL.GMC) * RACE_CAUCASIAN.GMC +
                      (DISTRESS_COVID.CMC * DISTRESS_COVID.GMCCM * AUDIT_BL.GMC) * RACE_CAUCASIAN.GMC +
                      (DISTRESS_COVID.CMC|SUBJECT_NUMBER),
                    data=ldatum)
sjPlot::tab_model(model.cdist5) # model ok


```
\
\
_10. hypo kdist5: ALCH ~ KESSLER * BL AUDIT * RACE_
```{r, message=FALSE}
model.kdist5.A = lmer(ALCH ~ DISTRESS_KESSLER.CMC *
                      DISTRESS_KESSLER.GMCCM *
                      AUDIT_BL.GMC +
                      RACE_API.GMC +
                      DISTRESS_KESSLER.CMC * RACE_API.GMC +
                      DISTRESS_KESSLER.GMCCM * RACE_API.GMC + 
                      AUDIT_BL.GMC * RACE_API.GMC + 
                      (DISTRESS_KESSLER.CMC * DISTRESS_KESSLER.GMCCM) * RACE_API.GMC +
                      (DISTRESS_KESSLER.CMC * AUDIT_BL.GMC) * RACE_API.GMC + 
                      (DISTRESS_KESSLER.GMCCM * AUDIT_BL.GMC) * RACE_API.GMC +
                      (DISTRESS_KESSLER.CMC * DISTRESS_KESSLER.GMCCM * AUDIT_BL.GMC) * RACE_API.GMC +
                      RACE_AFRICAN_BLACK.GMC +
                      DISTRESS_KESSLER.CMC * RACE_AFRICAN_BLACK.GMC +
                      DISTRESS_KESSLER.GMCCM * RACE_AFRICAN_BLACK.GMC + 
                      AUDIT_BL.GMC * RACE_AFRICAN_BLACK.GMC + 
                      (DISTRESS_KESSLER.CMC * DISTRESS_KESSLER.GMCCM) * RACE_AFRICAN_BLACK.GMC +
                      (DISTRESS_KESSLER.CMC * AUDIT_BL.GMC) * RACE_AFRICAN_BLACK.GMC + 
                      (DISTRESS_KESSLER.GMCCM * AUDIT_BL.GMC) * RACE_AFRICAN_BLACK.GMC +
                      (DISTRESS_KESSLER.CMC * DISTRESS_KESSLER.GMCCM * AUDIT_BL.GMC) * RACE_AFRICAN_BLACK.GMC +
                      RACE_AMERICAN_INDIAN.GMC +
                      DISTRESS_KESSLER.CMC * RACE_AMERICAN_INDIAN.GMC +
                      DISTRESS_KESSLER.GMCCM * RACE_AMERICAN_INDIAN.GMC + 
                      AUDIT_BL.GMC * RACE_AMERICAN_INDIAN.GMC + 
                      (DISTRESS_KESSLER.CMC * DISTRESS_KESSLER.GMCCM) * RACE_AMERICAN_INDIAN.GMC +
                      (DISTRESS_KESSLER.CMC * AUDIT_BL.GMC) * RACE_AMERICAN_INDIAN.GMC + 
                      (DISTRESS_KESSLER.GMCCM * AUDIT_BL.GMC) * RACE_AMERICAN_INDIAN.GMC +
                      (DISTRESS_KESSLER.CMC * DISTRESS_KESSLER.GMCCM * AUDIT_BL.GMC) * RACE_AMERICAN_INDIAN.GMC +
                      RACE_CAUCASIAN.GMC +
                      DISTRESS_KESSLER.CMC * RACE_CAUCASIAN.GMC +
                      DISTRESS_KESSLER.GMCCM * RACE_CAUCASIAN.GMC + 
                      AUDIT_BL.GMC * RACE_CAUCASIAN.GMC + 
                      (DISTRESS_KESSLER.CMC * DISTRESS_KESSLER.GMCCM) * RACE_CAUCASIAN.GMC +
                      (DISTRESS_KESSLER.CMC * AUDIT_BL.GMC) * RACE_CAUCASIAN.GMC + 
                      (DISTRESS_KESSLER.GMCCM * AUDIT_BL.GMC) * RACE_CAUCASIAN.GMC +
                      (DISTRESS_KESSLER.CMC * DISTRESS_KESSLER.GMCCM * AUDIT_BL.GMC) * RACE_CAUCASIAN.GMC +
                      (DISTRESS_KESSLER.CMC|SUBJECT_NUMBER),
                    data=ldatum)
sjPlot::tab_model(model.kdist5.A) # model conv fail

# rerun/rereport model w/ different optimizer if necessary
# https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html
model.kdist5.A_sc = update(model.kdist5.A, data=ldatum)
ss = getME(model.kdist5.A_sc, c("theta","fixef"))
model.kdist5.B = update(model.kdist5.A,
                        start=ss,
                        control=lmerControl(optimizer="bobyqa",
                                            optCtrl=list(maxfun=2e5))) 
sjPlot::tab_model(model.kdist5.B) # model ok


```
\
\
_11. hypo cdist6: ALCH ~ COVID * BL AUDIT * SUBX_
```{r, message=FALSE}
# run/report model
model.cdist6 = lmer(ALCH ~ DISTRESS_COVID.CMC *
                      DISTRESS_COVID.GMCCM *
                      AUDIT_BL.GMC *
                      SUBX_RM.GMC *
                      SUBX_BL.GMC *
                      (DISTRESS_COVID.CMC|SUBJECT_NUMBER),
                      data=ldatum)
sjPlot::tab_model(model.cdist6) # model ok


```
\
\
_12. hypo kdist6: ALCH ~ KESSLER * BL AUDIT * SUBX_
```{r, message=FALSE}
# run/report model
model.kdist6 = lmer(ALCH ~ DISTRESS_KESSLER.CMC *
                      DISTRESS_KESSLER.GMCCM *
                      AUDIT_BL.GMC *
                      SUBX_RM.GMC *
                      SUBX_BL.GMC +
                      (DISTRESS_KESSLER.CMC|SUBJECT_NUMBER),
                    data=ldatum)
sjPlot::tab_model(model.kdist6) # model ok


```
\
\
_13. hypo cdist7: ALCH ~ COVID * BL AUDIT * MENX_
```{r, message=FALSE}
# run/report model
model.cdist7 = lmer(ALCH ~ DISTRESS_COVID.CMC *
                      DISTRESS_COVID.GMCCM *
                      AUDIT_BL.GMC *
                      MENX_RM.GMC *
                      MENX_BL.GMC +
                      (DISTRESS_COVID.CMC|SUBJECT_NUMBER),
                    data=ldatum)
sjPlot::tab_model(model.cdist7) # model ok


```
\
\
_14. hypo kdist7: ALCH ~ KESSLER * BL AUDIT * MENX_
```{r, message=FALSE}
# run/report model
model.kdist7.A = lmer(ALCH ~ DISTRESS_KESSLER.CMC *
                      DISTRESS_KESSLER.GMCCM *
                      AUDIT_BL.GMC *
                      MENX_RM.GMC *
                      MENX_BL.GMC +
                      (DISTRESS_KESSLER.CMC|SUBJECT_NUMBER),
                    data=ldatum) # model conv fail

# rerun/rereport model w/ different optimizer if necessary
# https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html
model.kdist7.A_sc = update(model.kdist7.A, data=ldatum)
ss = getME(model.kdist7.A_sc, c("theta","fixef"))
model.kdist7.B = update(model.kdist7.A,
                        start=ss,
                        control=lmerControl(optimizer="bobyqa",
                                            optCtrl=list(maxfun=2e5)))
sjPlot::tab_model(model.kdist7.B) # model ok


```
\
\
_15. hypo cdist8: TBCO ~ COVID_
```{r, message=FALSE}
model.cdist8 = lmer(TOBACCO ~ DISTRESS_COVID.CMC *
                      DISTRESS_COVID.GMCCM +
                      (DISTRESS_COVID.CMC|SUBJECT_NUMBER),
                    data=ldatum)
sjPlot::tab_model(model.cdist8) # model ok


```
\
\
_16. hypo kdist8: TBCO ~ KESSLER_
```{r, message=FALSE}
model.kdist8 = lmer(TOBACCO ~ DISTRESS_KESSLER.CMC *
                      DISTRESS_KESSLER.GMCCM +
                      (DISTRESS_KESSLER.CMC|SUBJECT_NUMBER),
                    data=ldatum)
sjPlot::tab_model(model.kdist8) # model ok


```
\
\
_17. hypo cdist9: DRGS_OTHR ~ COVID_
```{r, message=FALSE}
model.cdist9.A = lmer(DRGS_OTHR ~ DISTRESS_COVID.CMC *
                      DISTRESS_COVID.GMCCM +
                      (DISTRESS_COVID.CMC|SUBJECT_NUMBER),
                    data=ldatum) # model conv fail

# rerun/rereport model w/ different optimizer if necessary
# https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html
model.cdist9.A_sc = update(model.cdist9.A, data=ldatum)
ss = getME(model.cdist9.A_sc, c("theta","fixef"))
model.cdist9.B = update(model.cdist9.A,
                        start=ss,
                        control=lmerControl(optimizer="bobyqa",
                                            optCtrl=list(maxfun=2e5)))
sjPlot::tab_model(model.cdist9.B) # model ok


```
\
\
_18. hypo kdist9: DRGS_OTHR ~ KESSLER_
```{r, message=FALSE}
model.kdist9 = lmer(DRGS_OTHR ~ DISTRESS_KESSLER.CMC *
                      DISTRESS_KESSLER.GMCCM +
                      (DISTRESS_KESSLER.CMC|SUBJECT_NUMBER),
                    data=ldatum)
sjPlot::tab_model(model.kdist9) # model ok


```
\
\
_19. hypo cdist10: ALCH ~ COVID * BL AUDIT * DAY_
```{r, message=FALSE}
model.cdist10.A = lmer(ALCH ~ DISTRESS_COVID.CMC *
                      DISTRESS_COVID.GMCCM *
                      AUDIT_BL.GMC *
                      DAY.CMC +
                      (DISTRESS_COVID.CMC|SUBJECT_NUMBER),
                    data=ldatum)
# sjPlot::tab_model(model.cdist10) # model conv fail

# rerun/rereport model w/ different optimizer if necessary
# https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html
model.cdist10.A_sc = update(model.cdist10.A, data=ldatum)
ss = getME(model.cdist10.A_sc, c("theta","fixef"))
model.cdist10.B = update(model.cdist10.A,
                        start=ss,
                        control=lmerControl(optimizer="bobyqa",
                                            optCtrl=list(maxfun=2e5)))
sjPlot::tab_model(model.cdist10.B) # model conv fail

```
\
\
_20. hypo kdist10: ALCH ~ KESSLER * BL AUDIT * DAY_
```{r, message=FALSE}
model.kdist10.A = lmer(ALCH ~ DISTRESS_KESSLER.CMC *
                      DISTRESS_KESSLER.GMCCM *
                      AUDIT_BL.GMC *
                      DAY.CMC +
                      (DISTRESS_KESSLER.CMC|SUBJECT_NUMBER),
                    data=ldatum)
# sjPlot::tab_model(model.kdist10) # model conv fail

# rerun/rereport model w/ different optimizer if necessary
# https://rstudio-pubs-static.s3.amazonaws.com/33653_57fc7b8e5d484c909b615d8633c01d51.html
model.kdist10.A_sc = update(model.kdist10.A, data=ldatum)
ss = getME(model.kdist10.A_sc, c("theta","fixef"))
model.kdist10.B = update(model.kdist10.A,
                        start=ss,
                        control=lmerControl(optimizer="bobyqa",
                                            optCtrl=list(maxfun=2e5)))
sjPlot::tab_model(model.kdist10.B) 
```

## plotting
_1. hist: AUDIT BL_                                        \
_2. bivar scatter: AUDIT BL x ALC_                         \
_3. bivar scatter panel: COVID (W)/COVID (B) X ALCH_       \
_4. bivar scatter panel: KESSLER (W)/KESSLER (B) X ALCH_   \
_5. bivar scatter panel: COVID x ALCH x AUDIT BL_          \
_6. bivar scatter panel: KESSLER x ALCH x AUDIT BL_        \
_7. bivar scatter: COVID x ALCH x AUDIT CHANGE GROUP_      \
_8. bivar scatter: KESSLER x ALCH by AUDIT CHANGE GROUP_   \
\
\
_1. hist: AUDIT BL_
```{r, message=FALSE, warning=FALSE}
# compute relevant vals
mean.AUDIT = round(mean(datum$AUDIT_BL), 2)

# contruct main plot
# png(file='../output/plt_hstgrm_AUDIT_BL.png',
#     width=125, height=100, units='mm', res=300)
ggplot(data=datum, aes(x=AUDIT_BL)) +
  geom_histogram(color='black', fill='gray', position='dodge') +
  labs(x='Baseline Alcohol Problems (AUDIT)', y='Count') +
  geom_text(aes(label=paste0('mean = ', mean.AUDIT), x=35, y=900),
            size=3.5, col='black', hjust=1) +
  theme_classic() +
  theme(legend.position='none')
# dev.off() # close plot


```
\
\
_2. bivar scatter: AUDIT BL x ALCH_
```{r, message=FALSE, warning=FALSE}
# select cols
cols = grep('SUBJECT_NUMBER|MORE_4_DRINK_PER_DAY', 
            colnames(datum))
tdatum = datum[cols]

# compute avg alch consumption across weeks
tdatum$ALCH.AVG = rowMeans(tdatum[, -(1)], na.rm = TRUE)

# remove unneeded cols
cols = grep('SUBJECT_NUMBER|ALCH.AVG', 
            colnames(tdatum))
tdatum = tdatum[cols]

# merge w/ datum (if ALCH.AVG isn't already there)
if ('ALCH.AVG' %in% colnames(datum)) {
cat('ALCH.AVG col exits!')
} else {
  datum = merge(datum, tdatum,
               by.x='SUBJECT_NUMBER',
               by.y='SUBJECT_NUMBER',
               all.x=TRUE, all.y=FALSE)
}

# cleanup
rm(cols, tdatum)

# extract vals from model
coef = round(summary(model.kdist1.B)$coefficients[c('AUDIT_BL.GMC'), 1], 2)
p = round(summary(model.kdist1.B)$coefficients[c('AUDIT_BL.GMC'), 5], 3)
n = sapply(ranef(model.kdist1.B), nrow)
n = as.numeric(str_extract(n, "[[:digit:]]+"))

# replace pval w/ label when applicable
if (p < 0.001) {
p = '<0.001'
}

# construct main plot
# png(file='../output/plt_scttr_AUDIT_x_drinks.png',
#     width=100, height=100, units='mm', res=300)
ggplot(datum, aes(x=AUDIT_BL, y=ALCH.AVG)) +
  geom_jitter(height=.1, width=.1,
              colour='gray', alpha=.5) +
  geom_smooth(method=lm, size=.66, color='black',
              level=.99, formula=y~x, se=FALSE) +
  geom_text(aes(label=paste0('B = ', coef), x=0, y=4.5),
            size=2.5, col='black', hjust=0) +
  geom_text(aes(label=paste0('p = ', p), x=5, y=4.5),
            size=2.5, col='black', hjust=0) +
  geom_text(aes(label=paste0('n = ', n), x=0, y=4.35),
            size=2.5, col='black', hjust=0) +
  coord_cartesian(xlim=c(-0.05, 32.05), ylim=c(-0.05, 4.5)) +
  scale_x_continuous(breaks=c(0, 16, 32)) +
  scale_y_continuous(breaks=c(0, 2, 4)) +
  labs(x='Basline Alcohol Problems (AUDIT)',
       y='Alcohol Use (DSMXC)') +
  theme_classic()
# dev.off() # close plot

# cleanup
rm(coef, n, p)


```
\
\
_3. bivar scatter panel: COVID (W)/COVID (B) X ALCH_
```{r, message=FALSE, warning=FALSE}
# COVID (W)
# extract vals from relevant model
coef.cdist.W = round(summary(model.cdist1)$coefficients[c('DISTRESS_COVID.CMC'), 1], 2)
p.cdist.W = round(summary(model.cdist1)$coefficients[c('DISTRESS_COVID.CMC'), 5], 3)
inter.cdist.W = parse_number(toString(fixef(model.cdist1)['(Intercept)']))
slope.cdist.W = parse_number(toString(fixef(model.cdist1)['DISTRESS_COVID.CMC']))

# create temp dataframe to plot lines
lines.cdist.W = data.frame(X=c(min(ldatum$DISTRESS_COVID, na.rm=T), 
                               max(ldatum$DISTRESS_COVID, na.rm=T)),
                           Y=c(inter.cdist.W + slope.cdist.W * min(ldatum$DISTRESS_COVID, na.rm=T),
                               inter.cdist.W + slope.cdist.W * max(ldatum$DISTRESS_COVID, na.rm=T)))

# replace pval w/ label when applicable
if (p.cdist.W < 0.001) {
p.cdist.W = '<0.001'
}

# construct subplot
cdist.W = ggplot(ldatum, aes(x=DISTRESS_COVID,
                   y=ALCH,
                   group=SUBJECT_NUMBER)) +
  geom_smooth(method='lm', se=FALSE,
              size=.33, color='gray75') +
  geom_line(aes(x=X, y=Y), data=lines.cdist.W, linetype='longdash', 
            colour='black', size=0.75, inherit.aes=F) +
  geom_text(aes(label=paste0('B = ', coef.cdist.W), x=8.5, y=4.25),
            size=3.5, col='black', hjust=1) +
  geom_text(aes(label=paste0('p = ', p.cdist.W), x=10, y=4.25),
            size=3.5, col='black', hjust=1) +
  coord_cartesian(xlim=c(.75, 10.25), ylim=c(0, 4.25)) +
  scale_x_continuous(breaks=c(1, 5.5, 10)) +
  scale_y_continuous(breaks=c(0, 2, 4)) +
  labs(x='Recent Distress (COVID-19 Survey) (W)',
       y='Alcohol Use (DSMXC)') +
  theme_classic()


# ---------------
# COVID (B)
# extract vals from relevant model
coef.cdist.B = round(summary(model.cdist1)$coefficients[c('DISTRESS_COVID.GMCCM'), 1], 2)
p.cdist.B = round(summary(model.cdist1)$coefficients[c('DISTRESS_COVID.GMCCM'), 5], 3)
inter.cdist.B = parse_number(toString(fixef(model.cdist1)['(Intercept)']))
slope.cdist.B = parse_number(toString(fixef(model.cdist1)['DISTRESS_COVID.GMCCM']))

# create temp dataframe to plot lines
lines.cdist.B = data.frame(X=c(min(ldatum$DISTRESS_COVID, na.rm=T), 
                             max(ldatum$DISTRESS_COVID, na.rm=T)),
                         Y=c(inter.cdist.B + slope.cdist.B * min(ldatum$DISTRESS_COVID, na.rm=T),
                             inter.cdist.B + slope.cdist.B * max(ldatum$DISTRESS_COVID, na.rm=T)))

# replace pval w/ label when applicable
if (p.cdist.B < 0.001) {
p.cdist.B = '<0.001'
}

# construct subplot
cdist.B = ggplot(ldatum, aes(x=DISTRESS_COVID,
                   y=ALCH,
                   group=SUBJECT_NUMBER)) +
  geom_smooth(method='lm', se=FALSE,
              size=.33, color='gray75') +
  geom_line(aes(x=X, y=Y), data=lines.cdist.B, linetype='longdash', 
            colour='black', size = 0.75, inherit.aes=F) +
  geom_text(aes(label=paste0('B = ', coef.cdist.B), x=8.5, y=4.25),
            size=3.5, col='black', hjust=1) +
  geom_text(aes(label=paste0('p = ', p.cdist.B), x=10, y=4.25),
            size=3.5, col='black', hjust=1) +
  coord_cartesian(xlim=c(.75, 10.25), ylim=c(0, 4.25)) +
  scale_x_continuous(breaks=c(1, 5.5, 10)) +
  scale_y_continuous(breaks=c(0, 2, 4)) +
  labs(x='Recent Distress (COVID-19 Survey) (B)',
       y='Alcohol Use (DSMXC)') +
  theme_classic()


# ---------------
# panel: COVID (W, B)
# combine/export subplots as panel
png(file='../output/plt_pnnl_cdist_x_drinks.png',
    width=300, height=150, units='mm', res=300) # init plot
ggarrange(cdist.W, cdist.B, ncol=2, nrow=1)
dev.off() # close plotting tool

# ---------------
# cleanup
rm(cdist.W, cdist.B)
rm(list=ls(pattern='p.|coef.|inter.|slope.|lines.'))


```
\
\
_4. bivar scatter panel: KESSLER (W)/KESSLER (B) X ALCH_
```{r, message=FALSE, warning=FALSE}
# KESSLER (W)
# extract vals from relevant model
coef.kdist.W = round(summary(model.kdist1.B)$coefficients[c('DISTRESS_KESSLER.CMC'), 1], 2)
p.kdist.W = round(summary(model.kdist1.B)$coefficients[c('DISTRESS_KESSLER.CMC'), 5], 3)
inter.kdist.W = parse_number(toString(fixef(model.kdist1.B)['(Intercept)']))
slope.kdist.W = parse_number(toString(fixef(model.kdist1.B)['DISTRESS_KESSLER.CMC']))

# create temp dataframe to plot lines
lines.kdist.W = data.frame(X=c(min(ldatum$DISTRESS_KESSLER, na.rm=T), 
                               max(ldatum$DISTRESS_KESSLER, na.rm=T)),
                           Y=c(inter.kdist.W + slope.kdist.W * min(ldatum$DISTRESS_KESSLER, na.rm=T),
                               inter.kdist.W + slope.kdist.W * max(ldatum$DISTRESS_KESSLER, na.rm=T)))

# replace pval w/ label when applicable
if (p.kdist.W < 0.001) {
p.kdist.W = '<0.001'
}

# construct subplot
kdist.W = ggplot(ldatum, aes(x=DISTRESS_KESSLER,
                   y=ALCH,
                   group=SUBJECT_NUMBER)) +
  geom_smooth(method='lm', se=FALSE,
              size=.33, color='gray75') +
  geom_line(aes(x=X, y=Y), data=lines.kdist.W, linetype='longdash', 
            colour='black', size = 0.75, inherit.aes=F) +
  geom_text(aes(label=paste0('B = ', coef.kdist.W), x=16, y=4.25),
            size=3.5, col='black', hjust=1) +
  geom_text(aes(label=paste0('p = ', p.kdist.W), x=20, y=4.25),
            size=3.5, col='black', hjust=1) +
  coord_cartesian(xlim=c(-.05, 20.05), ylim=c(-.05, 4.35)) +
  scale_x_continuous(breaks=c(0, 10, 20)) +
  scale_y_continuous(breaks=c(0, 2, 4)) +
  labs(x='Psychological Distress (KESSLER 5) (W)',
       y='Alcohol Use (DSMXC)') +
  theme_classic()

# ---------------
# KESSLER (B)
# extract vals from relevant model
coef.kdist.B = round(summary(model.kdist1.B)$coefficients[c('DISTRESS_KESSLER.GMCCM'), 1], 2)
p.kdist.B = round(summary(model.kdist1.B)$coefficients[c('DISTRESS_KESSLER.GMCCM'), 5], 3)
inter.kdist.B = parse_number(toString(fixef(model.kdist1.B)['(Intercept)']))
slope.kdist.B = parse_number(toString(fixef(model.kdist1.B)['DISTRESS_KESSLER.GMCCM']))

# create temp dataframe to plot lines
lines.kdist.B = data.frame(X=c(min(ldatum$DISTRESS_KESSLER, na.rm=T), 
                               max(ldatum$DISTRESS_KESSLER, na.rm=T)),
                           Y=c(inter.kdist.B + slope.kdist.B * min(ldatum$DISTRESS_KESSLER, na.rm=T),
                               inter.kdist.B + slope.kdist.B * max(ldatum$DISTRESS_KESSLER, na.rm=T)))

# replace pval w/ label when applicable
if (p.kdist.B < 0.001) {
p.kdist.B = '<0.001'
}

# construct subplot
kdist.B = ggplot(ldatum, aes(x=DISTRESS_KESSLER,
                   y=ALCH,
                   group=SUBJECT_NUMBER)) +
  geom_smooth(method='lm', se=FALSE,
              size=.33, color='gray75') +
  geom_line(aes(x=X, y=Y), data=lines.kdist.B, linetype='longdash', 
            colour='black', size = 0.75, inherit.aes=F) +
  geom_text(aes(label=paste0('B = ', coef.kdist.B), x=16, y=4.25),
            size=3.5, col='black', hjust=1) +
  geom_text(aes(label=paste0('p = ', p.kdist.B), x=20, y=4.25),
            size=3.5, col='black', hjust=1) +
  coord_cartesian(xlim=c(-.05, 20.05), ylim=c(-.05, 4.35)) +
  scale_x_continuous(breaks=c(0, 10, 20)) +
  scale_y_continuous(breaks=c(0, 2, 4)) +
  labs(x='Psychological Distress (KESSLER 5) (B)',
       y='Alcohol Use (DSMXC)') +
  theme_classic()

# ---------------
# panel: KESSLER (W, B)
# combine/export subplots as panel
png(file='../output/plt_pnnl_kdist_x_drinks.png',
    width=300, height=150, units='mm', res=300) # init plot
ggarrange(kdist.W, kdist.B, ncol=2, nrow=1)
dev.off() # close plotting tool

# ---------------
# cleanup
rm(kdist.W, kdist.B)
rm(list=ls(pattern='p.|coef.|inter.|slope.|lines.'))


```
\
\
_5. bivar scatter panel: COVID x ALCH x AUDIT BL_
```{r, message=FALSE, warning=FALSE}
# create col palette
poke1 = ichooseyou(pokemon='squirtle', spread=6)
poke2 = ichooseyou(pokemon='wartortle', spread=6)
poke3 = ichooseyou(pokemon='blastoise', spread=6)

# create new col to enhance gradient vis
ldatum$AUDIT_BL.log.1 = log(ldatum$AUDIT_BL+1)

# ---------------
# COVID x ALCH x AUDIT
# extract vals from relevant model
coef.cdist = round(summary(model.cdist1)$
                     coefficients[c('DISTRESS_COVID.CMC:AUDIT_BL.GMC'), 1], 2)
p.cdist = round(summary(model.cdist1)$
                  coefficients[c('DISTRESS_COVID.CMC:AUDIT_BL.GMC'), 5], 3)
n.cdist = sapply(ranef(model.cdist1), nrow)
n.cdist = as.numeric(str_extract(n.cdist, "[[:digit:]]+"))
inter.cdist = parse_number(toString(fixef(model.cdist1)['(Intercept)']))
slope.cdist = parse_number(toString(fixef(model.cdist1)['DISTRESS_COVID.CMC']))

# create temp dataframe to plot lines
lines.cdist = data.frame(X=c(min(ldatum$DISTRESS_COVID, na.rm=T), 
                             max(ldatum$DISTRESS_COVID, na.rm=T)),
                         Y=c(inter.cdist + slope.cdist * min(ldatum$DISTRESS_COVID, na.rm=T),
                             inter.cdist + slope.cdist * max(ldatum$DISTRESS_COVID, na.rm=T)))

# replace pval w/ label when applicable
if (p.cdist < 0.001) {
  p.cdist = '<0.001'
}

# construct main plot
F1 = ggplot(ldatum, aes(x=DISTRESS_COVID,
                        y=ALCH,
                        group=SUBJECT_NUMBER,
                        col=AUDIT_BL.log.1)) +
  geom_smooth(method='lm', se=FALSE, size=.33) +
  scale_colour_gradient2(low=poke3[3],
                         mid=poke1[1],
                         high=poke2[5],
                         midpoint=1,
                         guide='colourbar',
                         na.value=NA,
                         limits = c(-0.05, 3.05)) +
  geom_line(aes(x=X, y=Y), data=lines.cdist, linetype='longdash', 
            colour='black', size = 0.75, inherit.aes=F) +
  geom_text(aes(label=paste0('B = ', coef.cdist), x=8.5, y=4.5),
            size=3.5, col='black', hjust=1) +
  geom_text(aes(label=paste0('p = ', p.cdist), x=10, y=4.5),
            size=3.5, col='black', hjust=1) +
  geom_text(aes(label=paste0('n = ', n.cdist), x=10, y=4.35),
            size=3.5, col='black', hjust=1) +
  coord_cartesian(xlim=c(.75, 10.25), ylim=c(0, 4.50)) +
  scale_x_continuous(breaks=c(1, 5.5, 10)) +
  scale_y_continuous(breaks=c(0, 2, 4)) +
  labs(x='General Distress (W)',
       y='Alcohol Misuse',
       col='Alcohol Problems') +
  theme_classic() +
  theme(legend.position='top')


# ---------------
# subplot 1: AUDIT abstainers
# extract vals from relevant model
coef.cdist.AB = round(summary(model.cdist1.AB)$
                        coefficients[c('DISTRESS_COVID.CMC'), 1], 2)
p.cdist.AB = round(summary(model.cdist1.AB)$
                     coefficients[c('DISTRESS_COVID.CMC'), 5], 3)
n.cdist.AB = sapply(ranef(model.cdist1.AB), nrow)
n.cdist.AB = as.numeric(str_extract(n.cdist.AB, "[[:digit:]]+"))
inter.cdist.AB = parse_number(toString(fixef(model.cdist1.AB)['(Intercept)']))
slope.cdist.AB = parse_number(toString(fixef(model.cdist1.AB)['DISTRESS_COVID.CMC']))

# create temp dataframe to plot lines
lines.cdist.AB = data.frame(X=c(min(ldatum$DISTRESS_COVID, na.rm=T), 
                                max(ldatum$DISTRESS_COVID, na.rm=T)),
                            Y=c(inter.cdist.AB + slope.cdist.AB * min(ldatum$DISTRESS_COVID, na.rm=T),
                                inter.cdist.AB + slope.cdist.AB * max(ldatum$DISTRESS_COVID, na.rm=T)))

# replace pval w/ label when applicable
if (p.cdist.AB < 0.001) {
  p.cdist.AB = '<0.001'
}

# construct subplot
AB = ggplot(subset(ldatum, AUDIT_BL_GROUP == 'Abstainer'),
            aes(x=DISTRESS_COVID,
                y=ALCH,
                group=SUBJECT_NUMBER)) +
  geom_smooth(method='lm', se=FALSE,
              size=.33, color='#55788D') +
  geom_line(aes(x=X, y=Y), data=lines.cdist.AB, linetype='longdash', 
            colour='black', size = 0.75, inherit.aes = FALSE ) +
  geom_text(aes(label=paste0('B = ', coef.cdist.AB), x=7, y=4.5),
            size=3.5, col='black', hjust=1) +
  geom_text(aes(label=paste0('p = ', p.cdist.AB), x=10, y=4.5),
            size=3.5, col='black', hjust=1) +
  geom_text(aes(label=paste0('n = ', n.cdist.AB), x=10, y=4.25),
            size=3.5, col='black', hjust=1) +
  coord_cartesian(xlim=c(.75, 10.25), ylim=c(0, 4.50)) +
  scale_x_continuous(breaks=c(1, 5.5, 10)) +
  scale_y_continuous(breaks=c(0, 2, 4)) +
  labs(title='Abstainers',
       x='General Distress (W)',
       y='Alcohol Misuse') +
  theme_classic() +
  theme(legend.position='none',
        plot.title = element_text(hjust = 0.5))


# ---------------
# subset 2: AUDIT low-risk
# extract vals from relevant model
coef.cdist.LR = round(summary(model.cdist1.LR)$
                        coefficients[c('DISTRESS_COVID.CMC'), 1], 2)
p.cdist.LR = round(summary(model.cdist1.LR)$
                     coefficients[c('DISTRESS_COVID.CMC'), 5], 3)
n.cdist.LR = sapply(ranef(model.cdist1.LR), nrow)
n.cdist.LR = as.numeric(str_extract(n.cdist.LR, "[[:digit:]]+"))
inter.cdist.LR = parse_number(toString(fixef(model.cdist1.LR)['(Intercept)']))
slope.cdist.LR = parse_number(toString(fixef(model.cdist1.LR)['DISTRESS_COVID.CMC']))

# create temp dataframe to plot lines
lines.cdist.LR = data.frame(X=c(min(ldatum$DISTRESS_COVID, na.rm=T), 
                                max(ldatum$DISTRESS_COVID, na.rm=T)),
                            Y=c(inter.cdist.LR + slope.cdist.LR * min(ldatum$DISTRESS_COVID, na.rm=T),
                                inter.cdist.LR + slope.cdist.LR * max(ldatum$DISTRESS_COVID, na.rm=T)))

# replace pval w/ label when applicable
if (p.cdist.LR < 0.001) {
  p.cdist.LR = '<0.001'
}

# construct subplot
LR = ggplot(subset(ldatum, AUDIT_BL_GROUP == 'Low-Risk'),
            aes(x=DISTRESS_COVID,
                y=ALCH,
                group=SUBJECT_NUMBER)) +
  geom_smooth(method='lm', se=FALSE,
              size=.33, color='#62AB97') +
  geom_line(aes(x=X, y=Y), data=lines.cdist.LR, linetype='longdash', 
            colour='black', size = 0.75, inherit.aes = FALSE ) +
  geom_text(aes(label=paste0('B = ', coef.cdist.LR), x=7, y=4.5),
            size=3.5, col='black', hjust=1) +
  geom_text(aes(label=paste0('p = ', p.cdist.LR), x=10, y=4.5),
            size=3.5, col='black', hjust=1) +
  geom_text(aes(label=paste0('n = ', n.cdist.LR), x=10, y=4.25),
            size=3.5, col='black', hjust=1) +
  coord_cartesian(xlim=c(.75, 10.25), ylim=c(0, 4.50)) +
  scale_x_continuous(breaks=c(1, 5.5, 10)) +
  scale_y_continuous(breaks=c(0, 2, 4)) +
  labs(title='Low Risk',
       x='General Distress (W)',
       y='Alcohol Misuse') +
  theme_classic() +
  theme(legend.position='none',
        plot.title = element_text(hjust = 0.5))


# ---------------
# subset 3: AUDIT high-risk
# extract vals from relevant model
coef.cdist.HR = round(summary(model.cdist1.HR)$
                        coefficients[c('DISTRESS_COVID.CMC'), 1], 2)
p.cdist.HR = round(summary(model.cdist1.HR)$
                     coefficients[c('DISTRESS_COVID.CMC'), 5], 3)
n.cdist.HR = sapply(ranef(model.cdist1.HR), nrow)
n.cdist.HR = as.numeric(str_extract(n.cdist.HR, "[[:digit:]]+"))
inter.cdist.HR = parse_number(toString(fixef(model.cdist1.HR)['(Intercept)']))
slope.cdist.HR = parse_number(toString(fixef(model.cdist1.HR)['DISTRESS_COVID.CMC']))

# create temp dataframe to plot lines
lines.cdist.HR = data.frame(X=c(min(ldatum$DISTRESS_COVID, na.rm=T), 
                                max(ldatum$DISTRESS_COVID, na.rm=T)),
                            Y=c(inter.cdist.HR + slope.cdist.HR * min(ldatum$DISTRESS_COVID, na.rm=T),
                                inter.cdist.HR + slope.cdist.HR * max(ldatum$DISTRESS_COVID, na.rm=T)))

# replace pval w/ label when applicable
if (p.cdist.HR < 0.001) {
  p.cdist.HR = '<0.001'
}

# construct subplot
HR = ggplot(subset(ldatum, AUDIT_BL_GROUP == 'High-Risk'),
            aes(x=DISTRESS_COVID,
                y=ALCH,
                group=SUBJECT_NUMBER)) +
  geom_smooth(method='lm', se=FALSE,
              size=.33, color='#B7CD90') +
  geom_line(aes(x=X, y=Y), data=lines.cdist.HR, linetype='longdash', 
            colour='black', size = 0.75, inherit.aes = FALSE ) +
  geom_text(aes(label=paste0('B = ', coef.cdist.HR), x=7, y=4.5),
            size=3.5, col='black', hjust=1) +
  geom_text(aes(label=paste0('p = ', p.cdist.HR), x=10, y=4.5),
            size=3.5, col='black', hjust=1) +
  geom_text(aes(label=paste0('n = ', n.cdist.HR), x=10, y=4.25),
            size=3.5, col='black', hjust=1) +
  coord_cartesian(xlim=c(.75, 10.25), ylim=c(0, 4.50)) +
  scale_x_continuous(breaks=c(1, 5.5, 10)) +
  scale_y_continuous(breaks=c(0, 2, 4)) +
  labs(title='High Risk',
       x='General Distress (W)',
       y='Alcohol Misuse') +
  theme_classic() +
  theme(legend.position='none',
        plot.title = element_text(hjust = 0.5))


# ---------------
# subset 4: AUDIT AUD
# extract vals from relevant model
coef.cdist.AUD = round(summary(model.cdist1.AUD)$
                         coefficients[c('DISTRESS_COVID.CMC'), 1], 2)
p.cdist.AUD = round(summary(model.cdist1.AUD)$
                      coefficients[c('DISTRESS_COVID.CMC'), 5], 3)
n.cdist.AUD = sapply(ranef(model.cdist1.AUD), nrow)
n.cdist.AUD = as.numeric(str_extract(n.cdist.AUD, "[[:digit:]]+"))
inter.cdist.AUD = parse_number(toString(fixef(model.cdist1.AUD)['(Intercept)']))
slope.cdist.AUD = parse_number(toString(fixef(model.cdist1.AUD)['DISTRESS_COVID.CMC']))

# create temp dataframe to plot lines
lines.cdist.AUD = data.frame(X=c(min(ldatum$DISTRESS_COVID, na.rm=T), 
                                max(ldatum$DISTRESS_COVID, na.rm=T)),
                            Y=c(inter.cdist.AUD + slope.cdist.AUD * min(ldatum$DISTRESS_COVID, na.rm=T),
                                inter.cdist.AUD + slope.cdist.AUD * max(ldatum$DISTRESS_COVID, na.rm=T)))

# replace pval w/ label when applicable
if (p.cdist.AUD < 0.001) {
  p.cdist.AUD = '<0.001'
}

# construct subplot
AUD = ggplot(subset(ldatum, AUDIT_BL_GROUP == 'AUD'),
             aes(x=DISTRESS_COVID,
                 y=ALCH,
                 group=SUBJECT_NUMBER)) +
  geom_smooth(method='lm', se=FALSE,
              size=.33, color='#DEDF8B') +
  geom_line(aes(x=X, y=Y), data=lines.cdist.AUD, linetype='longdash', 
            colour='black', size = 0.75, inherit.aes = FALSE ) +
  geom_text(aes(label=paste0('B = ', coef.cdist.AUD), x=7, y=4.5),
            size=3.5, col='black', hjust=1) +
  geom_text(aes(label=paste0('p = ', p.cdist.AUD), x=10, y=4.5),
            size=3.5, col='black', hjust=1) +
  geom_text(aes(label=paste0('n = ', n.cdist.AUD), x=10, y=4.25),
            size=3.5, col='black', hjust=1) +
  coord_cartesian(xlim=c(.75, 10.25), ylim=c(0, 4.50)) +
  scale_x_continuous(breaks=c(1, 5.5, 10)) +
  scale_y_continuous(breaks=c(0, 2, 4)) +
  labs(title='AUD',
       x='General Distress (W)',
       y='Alcohol Misuse') +
  theme_classic() +
  theme(legend.position='none',
        plot.title = element_text(hjust = 0.5))


# ---------------
# combine subplots as panel
F2 = ggarrange(HR, AUD, AB, LR,
               ncol=2, nrow=2)

# combine/export main plot + subplots as panel
png(file='../output/plt_pnnl_cdist_x_drinks_x_AUDIT_gendist.png',
    width=300, height=150, units='mm', res=300) # init plot
ggarrange(F1, F2,
          ncol=2, nrow=1,
          labels=c('A', 'B'))
dev.off() # close plot tool


# ---------------
# cleanup
rm(HR, AUD, AB, LR, F1, F2)
rm(list=ls(pattern='p.|coef.|n.|inter.|slope.|lines.'))


```
\
\
_6. bivar scatter panel: KESSLER x ALCH x AUDIT BL_        \
```{r, message=FALSE, warning=FALSE}
# create col palette
poke1 = ichooseyou(pokemon='squirtle', spread=6)
poke2 = ichooseyou(pokemon='wartortle', spread=6)
poke3 = ichooseyou(pokemon='blastoise', spread=6)

# create new col to enhance gradient vis
ldatum$AUDIT_BL.log.1 = log(ldatum$AUDIT_BL+1)


# ---------------
# KESSLER x DRINKS x AUDIT
# extract vals from relevant model
coef.kdist = round(summary(model.kdist1.B)$
                       coefficients[c('DISTRESS_KESSLER.CMC:AUDIT_BL.GMC'), 1], 2)
p.kdist = round(summary(model.kdist1.B)$
                    coefficients[c('DISTRESS_KESSLER.CMC:AUDIT_BL.GMC'), 5], 3)
n.kdist = sapply(ranef(model.kdist1.B), nrow)
n.kdist = as.numeric(str_extract(n.kdist, "[[:digit:]]+"))
inter.kdist = parse_number(toString(fixef(model.kdist1.B)['(Intercept)']))
slope.kdist = parse_number(toString(fixef(model.kdist1.B)['DISTRESS_KESSLER.CMC']))

# create temp dataframe to plot lines
lines.kdist = data.frame(X=c(min(ldatum$DISTRESS_KESSLER, na.rm=T), 
                             max(ldatum$DISTRESS_KESSLER, na.rm=T)),
                         Y=c(inter.kdist + slope.kdist * min(ldatum$DISTRESS_KESSLER, na.rm=T),
                             inter.kdist + slope.kdist * max(ldatum$DISTRESS_KESSLER, na.rm=T)))

# replace pval w/ label when applicable
if (p.kdist < 0.001) {
p.kdist = '<0.001'
}

# construct main plot
F1 = ggplot(ldatum, aes(x=DISTRESS_KESSLER,
                        y=ALCH,
                        group=SUBJECT_NUMBER,
                        col=AUDIT_BL.log.1)) +
  geom_smooth(method='lm', se=FALSE, size=.33) +
  scale_colour_gradient2(low=poke3[3],
                         mid=poke1[1],
                         high=poke2[5],
                         midpoint=1,
                         guide='colourbar',
                         na.value=NA,
                         limits = c(-0.05, 3.05)) +
  geom_line(aes(x=X, y=Y), data=lines.kdist, linetype='longdash', 
            colour='black', size = 0.75, inherit.aes = FALSE ) +
  geom_text(aes(label=paste0('B = ', coef.kdist), x=16, y=4.5),
            size=3.5, col='black', hjust=1) +
  geom_text(aes(label=paste0('p = ', p.kdist), x=20, y=4.5),
            size=3.5, col='black', hjust=1) +
  geom_text(aes(label=paste0('n = ', n.kdist), x=20, y=4.35),
            size=3.5, col='black', hjust=1) +
  coord_cartesian(xlim=c(-0.05, 20.05), ylim=c(-0.05, 4.50)) +
  scale_x_continuous(breaks=c(0, 10, 20)) +
  scale_y_continuous(breaks=c(0, 2, 4)) +
  labs(x='Psychological Distress (W)',
       y='Alcohol Misuse',
       col='Alcohol Problems') +
  theme_classic() +
  theme(legend.position='top')


# ---------------
# subplot 1: AUDIT abstainers
# extract vals from relevant model
coef.kdist.AB = round(summary(model.kdist1.AB)$
                       coefficients[c('DISTRESS_KESSLER.CMC'), 1], 2)
p.kdist.AB = round(summary(model.kdist1.AB)$
                    coefficients[c('DISTRESS_KESSLER.CMC'), 5], 3)
n.kdist.AB = sapply(ranef(model.kdist1.AB), nrow)
n.kdist.AB = as.numeric(str_extract(n.kdist.AB, "[[:digit:]]+"))
inter.kdist.AB = parse_number(toString(fixef(model.kdist1.AB)['(Intercept)']))
slope.kdist.AB = parse_number(toString(fixef(model.kdist1.AB)['DISTRESS_KESSLER.CMC']))

# create temp dataframe to plot lines
lines.kdist.AB = data.frame(X=c(min(ldatum$DISTRESS_KESSLER, na.rm=T), 
                                max(ldatum$DISTRESS_KESSLER, na.rm=T)),
                            Y=c(inter.kdist.AB + slope.kdist.AB * min(ldatum$DISTRESS_KESSLER, na.rm=T),
                                inter.kdist.AB + slope.kdist.AB * max(ldatum$DISTRESS_KESSLER, na.rm=T)))

# replace pval w/ label when applicable
if (p.kdist.AB < 0.001) {
p.kdist.AB = '<0.001'
}

# construct subplot
AB = ggplot(subset(ldatum, AUDIT_BL_GROUP == 'Abstainer'),
            aes(x=DISTRESS_KESSLER,
                y=ALCH,
                group=SUBJECT_NUMBER)) +
  geom_smooth(method='lm', se=FALSE,
              size=.33, color='#55788D') +
  geom_line(aes(x=X, y=Y), data=lines.kdist.AB, linetype='longdash', 
            colour='black', size = 0.75, inherit.aes = FALSE ) +
  geom_text(aes(label=paste0('B = ', coef.kdist.AB), x=13.5, y=4.5),
            size=3.5, col='black', hjust=1) +
  geom_text(aes(label=paste0('p = ', p.kdist.AB), x=20, y=4.5),
            size=3.5, col='black', hjust=1) +
  geom_text(aes(label=paste0('n = ', n.kdist.AB), x=20, y=4.25),
            size=3.5, col='black', hjust=1) +
  coord_cartesian(xlim=c(-0.05, 20.05), ylim=c(-0.05, 4.50)) +
  scale_x_continuous(breaks=c(0, 10, 20)) +
  scale_y_continuous(breaks=c(0, 2, 4)) +
  labs(title='Abstainers',
       x='Psychological Distress (W)',
       y='Alcohol Misuse') +
  theme_classic() +
  theme(legend.position='none',
        plot.title = element_text(hjust = 0.5))


# ---------------
# subset 2: AUDIT low-risk
# extract vals from relevant model
coef.kdist.LR = round(summary(model.kdist1.LR)$
                       coefficients[c('DISTRESS_KESSLER.CMC'), 1], 2)
p.kdist.LR = round(summary(model.kdist1.LR)$
                    coefficients[c('DISTRESS_KESSLER.CMC'), 5], 3)
n.kdist.LR = sapply(ranef(model.kdist1.LR), nrow)
n.kdist.LR = as.numeric(str_extract(n.kdist.LR, "[[:digit:]]+"))
inter.kdist.LR = parse_number(toString(fixef(model.kdist1.LR)['(Intercept)']))
slope.kdist.LR = parse_number(toString(fixef(model.kdist1.LR)['DISTRESS_KESSLER.CMC']))

# create temp dataframe to plot lines
lines.kdist.LR = data.frame(X=c(min(ldatum$DISTRESS_KESSLER, na.rm=T), 
                             max(ldatum$DISTRESS_KESSLER, na.rm=T)),
                         Y=c(inter.kdist.LR + slope.kdist.LR * min(ldatum$DISTRESS_KESSLER, na.rm=T),
                             inter.kdist.LR + slope.kdist.LR * max(ldatum$DISTRESS_KESSLER, na.rm=T)))

# replace pval w/ label when applicable
if (p.kdist.LR < 0.001) {
p.kdist.LR = '<0.001'
}

# construct subplot
LR = ggplot(subset(ldatum, AUDIT_BL_GROUP == 'Low-Risk'),
            aes(x=DISTRESS_KESSLER,
                y=ALCH,
                group=SUBJECT_NUMBER)) +
  geom_smooth(method='lm', se=FALSE,
              size=.33, color='#62AB97') +
  geom_line(aes(x=X, y=Y), data=lines.kdist.AB, linetype='longdash', 
            colour='black', size = 0.75, inherit.aes = FALSE ) +
  geom_text(aes(label=paste0('B = ', coef.kdist.LR), x=13.5, y=4.5),
            size=3.5, col='black', hjust=1) +
  geom_text(aes(label=paste0('p = ', p.kdist.LR), x=20, y=4.5),
            size=3.5, col='black', hjust=1) +
  geom_text(aes(label=paste0('n = ', n.kdist.LR), x=20, y=4.25),
            size=3.5, col='black', hjust=1) +
  coord_cartesian(xlim=c(-0.05, 20.05), ylim=c(-0.05, 4.50)) +
  scale_x_continuous(breaks=c(0, 10, 20)) +
  scale_y_continuous(breaks=c(0, 2, 4)) +
  labs(title='Low Risk',
       x='Psychological Distress (W)',
       y='Alcohol Misuse') +
  theme_classic() +
  theme(legend.position='none',
        plot.title = element_text(hjust = 0.5))


# ---------------
# subset 3: AUDIT high-risk
# extract vals from relevant model
coef.kdist.HR = round(summary(model.kdist1.HR.B)$
                       coefficients[c('DISTRESS_KESSLER.CMC'), 1], 2)
p.kdist.HR = round(summary(model.kdist1.HR.B)$
                    coefficients[c('DISTRESS_KESSLER.CMC'), 5], 3)
n.kdist.HR = sapply(ranef(model.kdist1.HR.B), nrow)
n.kdist.HR = as.numeric(str_extract(n.kdist.HR, "[[:digit:]]+"))
inter.kdist.HR.B = parse_number(toString(fixef(model.kdist1.HR.B)['(Intercept)']))
slope.kdist.HR.B = parse_number(toString(fixef(model.kdist1.HR.B)['DISTRESS_KESSLER.CMC']))

# create temp dataframe to plot lines
lines.kdist.HR.B = data.frame(X=c(min(ldatum$DISTRESS_KESSLER, na.rm=T), 
                             max(ldatum$DISTRESS_KESSLER, na.rm=T)),
                         Y=c(inter.kdist.HR.B + slope.kdist.HR.B * min(ldatum$DISTRESS_KESSLER, na.rm=T),
                             inter.kdist.HR.B + slope.kdist.HR.B * max(ldatum$DISTRESS_KESSLER, na.rm=T)))

# replace pval w/ label when applicable
if (p.kdist.HR < 0.001) {
p.kdist.HR = '<0.001'
}

# construct subplot
HR = ggplot(subset(ldatum, AUDIT_BL_GROUP == 'High-Risk'),
            aes(x=DISTRESS_KESSLER,
                y=ALCH,
                group=SUBJECT_NUMBER)) +
  geom_smooth(method='lm', se=FALSE,
              size=.33, color='#B7CD90') +
  geom_line(aes(x=X, y=Y), data=lines.kdist.HR.B, linetype='longdash', 
            colour='black', size = 0.75, inherit.aes = FALSE ) +
  geom_text(aes(label=paste0('B = ', coef.kdist.HR), x=13.5, y=4.5),
            size=3.5, col='black', hjust=1) +
  geom_text(aes(label=paste0('p = ', p.kdist.HR), x=20, y=4.5),
            size=3.5, col='black', hjust=1) +
  geom_text(aes(label=paste0('n = ', n.kdist.HR), x=20, y=4.25),
            size=3.5, col='black', hjust=1) +
  coord_cartesian(xlim=c(-0.05, 20.05), ylim=c(-0.05, 4.50)) +
  scale_x_continuous(breaks=c(0, 10, 20)) +
  scale_y_continuous(breaks=c(0, 2, 4)) +
  labs(title='High Risk',
       x='Psychological Distress (W)',
       y='Alcohol Miuse') +
  theme_classic() +
  theme(legend.position='none',
        plot.title = element_text(hjust = 0.5))


# ---------------
# subset 4: AUDIT AUD
# extract vals from relevant model
coef.kdist.AUD = round(summary(model.kdist1.AUD)$
                       coefficients[c('DISTRESS_KESSLER.CMC'), 1], 2)
p.kdist.AUD = round(summary(model.kdist1.AUD)$
                    coefficients[c('DISTRESS_KESSLER.CMC'), 5], 3)
n.kdist.AUD = sapply(ranef(model.kdist1.AUD), nrow)
n.kdist.AUD = as.numeric(str_extract(n.kdist.AUD, "[[:digit:]]+"))
inter.kdist.AUD = parse_number(toString(fixef(model.kdist1.AUD)['(Intercept)']))
slope.kdist.AUD = parse_number(toString(fixef(model.kdist1.AUD)['DISTRESS_KESSLER.CMC']))

# create temp dataframe to plot lines
lines.kdist.AUD = data.frame(X=c(min(ldatum$DISTRESS_KESSLER, na.rm=T), 
                             max(ldatum$DISTRESS_KESSLER, na.rm=T)),
                         Y=c(inter.kdist.AUD + slope.kdist.AUD * min(ldatum$DISTRESS_KESSLER, na.rm=T),
                             inter.kdist.AUD + slope.kdist.AUD * max(ldatum$DISTRESS_KESSLER, na.rm=T)))

# replace pval w/ label when applicable
if (p.kdist.AUD < 0.001) {
p.kdist.AUD = '<0.001'
}

# construct subplot
AUD = ggplot(subset(ldatum, AUDIT_BL_GROUP == 'AUD'),
            aes(x=DISTRESS_KESSLER,
                y=ALCH,
                group=SUBJECT_NUMBER)) +
  geom_smooth(method='lm', se=FALSE,
              size=.33, color='#DEDF8B') +
  geom_line(aes(x=X, y=Y), data=lines.kdist.AUD, linetype='longdash', 
            colour='black', size = 0.75, inherit.aes = FALSE ) +
  geom_text(aes(label=paste0('B = ', coef.kdist.AUD), x=13.5, y=4.5),
            size=3.5, col='black', hjust=1) +
  geom_text(aes(label=paste0('p = ', p.kdist.AUD), x=20, y=4.5),
            size=3.5, col='black', hjust=1) +
  geom_text(aes(label=paste0('n = ', n.kdist.AUD), x=20, y=4.25),
            size=3.5, col='black', hjust=1) +
  coord_cartesian(xlim=c(-0.05, 20.05), ylim=c(-0.05, 4.50)) +
  scale_x_continuous(breaks=c(0, 10, 20)) +
  scale_y_continuous(breaks=c(0, 2, 4)) +
  labs(title='AUD',
       x='Psychological Distress (W)',
       y='Alcohol Misuse') +
  theme_classic() +
  theme(legend.position='none',
        plot.title = element_text(hjust = 0.5))


# ---------------
# combine subplots as panel
F2 = ggarrange(HR, AUD, AB, LR,
               ncol=2, nrow=2)

# combine/export main plot + subplots as panel
png(file='../output/plt_pnnl_kdist_x_drinks_x_AUDIT_psychdist.png',
    width=300, height=150, units='mm', res=300) # init plot
ggarrange(F1, F2, ncol=2, nrow=1)
dev.off() # close plot tool

# ---------------
# cleanup
rm(HR, AUD, AB, LR, F1, F2)
rm(list=ls(pattern='p.|coef.|n.|inter.|slope.|lines.'))


```
\
\
_7. bivar scatter: COVID x ALCH x AUDIT CHANGE GROUP_      \
```{r, message=FALSE, warning=FALSE}
# AUDIT CHANGE group: no change
CN = ggplot(subset(ldatum, AUDIT_CHANGE_GROUP == 'No Change'),
            aes(x=DISTRESS_COVID,
                y=ALCH,
                group=SUBJECT_NUMBER)) +
  geom_smooth(method='lm', se=FALSE,
              size=.33, color='gray75') +
  geom_abline(aes(intercept=parse_number(toString(fixef(model.cdist6.CN)['(Intercept)'])),
                  slope=parse_number(toString(fixef(model.cdist6.CN)['DISTRESS_COVID.CMC']))),
              size=1, linetype='solid') +
  coord_cartesian(xlim=c(.75, 10.25), ylim=c(0, 4.50)) +
  scale_x_continuous(breaks=c(1, 5.5, 10)) +
  scale_y_continuous(breaks=c(0, 2, 4)) +
  labs(title='No Change',
       x='Recent Distress (COVID-19) (W)',
       y='Alcohol Use (DSMXC)') +
  theme_classic() +
  theme(legend.position='none',
        plot.title = element_text(hjust = 0.5))

# AUDIT CHANGE group: decrease
CD = ggplot(subset(ldatum, AUDIT_CHANGE_GROUP == 'Decrease'),
            aes(x=DISTRESS_COVID,
                y=ALCH,
                group=SUBJECT_NUMBER)) +
  geom_smooth(method='lm', se=FALSE,
              size=.33, color='gray75') +
  geom_abline(aes(intercept=parse_number(toString(fixef(model.cdist6.CD.B)['(Intercept)'])),
                  slope=parse_number(toString(fixef(model.cdist6.CD.B)['DISTRESS_COVID.CMC']))),
              size=1, linetype='solid') +
  coord_cartesian(xlim=c(.75, 10.25), ylim=c(0, 4.50)) +
  scale_x_continuous(breaks=c(1, 5.5, 10)) +
  scale_y_continuous(breaks=c(0, 2, 4)) +
  labs(title='Decrease',
       x='Recent Distress (COVID-19) (W)',
       y='Alcohol Use (DSMXC)') +
  theme_classic() +
  theme(legend.position='none',
        plot.title = element_text(hjust = 0.5))

# AUDIT CHANGE group: increase
CI = ggplot(subset(ldatum, AUDIT_CHANGE_GROUP == 'Increase'),
            aes(x=DISTRESS_COVID,
                y=ALCH,
                group=SUBJECT_NUMBER)) +
  geom_smooth(method='lm', se=FALSE,
              size=.33, color='gray75') +
  geom_abline(aes(intercept=parse_number(toString(fixef(model.cdist6.CI)['(Intercept)'])),
                  slope=parse_number(toString(fixef(model.cdist6.CI)['DISTRESS_COVID.CMC']))),
              size=1, linetype='solid') +
  coord_cartesian(xlim=c(.75, 10.25), ylim=c(0, 4.50)) +
  scale_x_continuous(breaks=c(1, 5.5, 10)) +
  scale_y_continuous(breaks=c(0, 2, 4)) +
  labs(title='Increase',
       x='Recent Distress (COVID-19) (W)',
       y='Alcohol Use (DSMXC)') +
  theme_classic() +
  theme(legend.position='none',
        plot.title = element_text(hjust = 0.5))

# combine AUDIT CHANGE groups, export as panel
# png(file='../output/plt_pnnl_alch_cdist_x_alch_x_AUDITchange.png',
#     width=225, height=75,
#     units='mm', res=300) # init plot
ggarrange(CN, CD, CI, ncol=3, nrow=1,
          labels=c('A', 'B', 'C'))
# dev.off() # close plot tool
```
\
\
_8. bivar scatter: KESSLER x ALCH by AUDIT CHANGE GROUP_   \
```{r, message=FALSE, warning=FALSE}
# AUDIT CHANGE group: no change
CN = ggplot(subset(ldatum, AUDIT_CHANGE_GROUP == 'No Change'),
            aes(x=DISTRESS_KESSLER,
                y=ALCH,
                group=SUBJECT_NUMBER)) +
  geom_smooth(method='lm', se=FALSE,
              size=.33, color='gray75') +
  geom_abline(aes(intercept=parse_number(toString(fixef(model.kdist6.CN)['(Intercept)'])),
                  slope=parse_number(toString(fixef(model.kdist6.CN)['DISTRESS_KESSLER.CMC']))),
              size=1, linetype='solid') +
  coord_cartesian(xlim=c(.75, 20.25), ylim=c(0, 4.50)) +
  scale_x_continuous(breaks=c(0, 10, 20)) +
  scale_y_continuous(breaks=c(0, 2, 4)) +
  labs(title='No Change',
       x='Psychological Distress (KESSLER 5) (W)',
       y='Alcohol Use (DSMXC)') +
  theme_classic() +
  theme(legend.position='none',
        plot.title = element_text(hjust = 0.5))

# AUDIT CHANGE group: decrease
CD = ggplot(subset(ldatum, AUDIT_CHANGE_GROUP == 'Decrease'),
            aes(x=DISTRESS_KESSLER,
                y=ALCH,
                group=SUBJECT_NUMBER)) +
  geom_smooth(method='lm', se=FALSE,
              size=.33, color='gray75') +
  geom_abline(aes(intercept=parse_number(toString(fixef(model.kdist6.CD.B)['(Intercept)'])),
                  slope=parse_number(toString(fixef(model.kdist6.CD.B)['DISTRESS_KESSLER.CMC']))),
              size=1, linetype='solid') +
  coord_cartesian(xlim=c(.75, 20.25), ylim=c(0, 4.50)) +
  scale_x_continuous(breaks=c(0, 10, 20)) +
  scale_y_continuous(breaks=c(0, 2, 4)) +
  labs(title='Decrease',
       x='Psychological Distress (KESSLER 5) (W)',
       y='Alcohol Use (DSMXC)') +
  theme_classic() +
  theme(legend.position='none',
        plot.title = element_text(hjust = 0.5))

# AUDIT CHANGE group: increase
CI = ggplot(subset(ldatum, AUDIT_CHANGE_GROUP == 'Increase'),
            aes(x=DISTRESS_KESSLER,
                y=ALCH,
                group=SUBJECT_NUMBER)) +
  geom_smooth(method='lm', se=FALSE,
              size=.33, color='gray75') +
  geom_abline(aes(intercept=parse_number(toString(fixef(model.kdist6.CI)['(Intercept)'])),
                  slope=parse_number(toString(fixef(model.kdist6.CI)['DISTRESS_KESSLER.CMC']))),
              size=1, linetype='solid') +
  coord_cartesian(xlim=c(.75, 20.25), ylim=c(0, 4.50)) +
  scale_x_continuous(breaks=c(0, 10, 20)) +
  scale_y_continuous(breaks=c(0, 2, 4)) +
  labs(title='Increase',
       x='Psychological Distress (KESSLER 5) (W)',
       y='Alcohol Use (DSMXC)') +
  theme_classic() +
  theme(legend.position='none',
        plot.title = element_text(hjust = 0.5))

# combine AUDIT CHANGE groups, export as panel
# png(file='../output/plt_pnnl_alch_kdist_x_alch_x_AUDITchange.png',
#     width=225, height=75,
#     units='mm', res=300) # init plot
ggarrange(CN, CD, CI, ncol=3, nrow=1,
          labels=c('A', 'B', 'C'))
# dev.off() # close plot tool
```
