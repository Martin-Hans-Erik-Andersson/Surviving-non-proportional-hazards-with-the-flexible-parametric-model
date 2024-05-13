library(tidyverse)
library(knitr)
library(haven)
library(dplyr)
library(survival)
library(flexsurv)
library(rstpm2)
library(ggfortify)
library(broom)
library(ggsurvfit)
library(lubridate)
library(gtsummary)




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# Cleaning dermatology data.


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -






# Real data

rawData <- read_sas("data/anapop.sas7bdat")  %>%
  select(-BIOLOGIC, - Methotrexate, -death_date, -emig_date, -birth_date, -outVTE, -eos_date) %>%
  rename(BIOLOGIC = BIOLOGIC1, 
         Methotrexate = Methotrexate1, 
         birth_date = birth_date1, 
         death_date = death_date1, 
         emig_date = emig_date1, 
         outVTE = outVTE1,
         eos_date = eos_date1)


# Scrambled Data

rawData2 <- read_sas("data/final2.sas7bdat") %>%
 rename(BIOLOGIC = BIOLOGIC1, 
         Methotrexate = Methotrexate1, 
         birth_date = birth_date1, 
         death_date = death_date1, 
         emig_date = emig_date1, 
         outVTE = outVTE1,
         eos_date = eos_date1)




# Duplicating any person who received both methods of treatment.
dermatologyData <- rawData %>%
  mutate(Test = if_else(Methotrexate %>% is.na() | BIOLOGIC %>% is.na(), 0, 1 ))


dermatologyData <- rbind(dermatologyData,
                         dermatologyData %>% 
                      filter(Test == 1) %>%
                      mutate(Test = 2))






# Add which treatment each observation received.

dermatologyData <- dermatologyData %>%
  mutate(Treatment = case_when(Test == 2 ~ "Biological",
                               
                               
                               Test == 1 ~ "Methotrexate",
                               
                               
                               Test == 0 ~ if_else(Methotrexate %>% is.na(), 
                                                   "Biological", 
                                                   "Methotrexate"))
         )





# Add what type of time each observation has (Death, emigration, end of study or VTE).


# Format the data for simplicity


temp1 <- dermatologyData %>%
  filter(Test == 1) %>%
  mutate(death_date = ifelse(death_date - Methotrexate >= 0, 
                               death_date - Methotrexate, 
                               NA),
         emig_date = ifelse(emig_date - Methotrexate >= 0, 
                              emig_date - Methotrexate,
                              NA),
         eos_date = ifelse(eos_date - Methotrexate >= 0,
                             eos_date - Methotrexate,
                             NA),
         outVTE = ifelse(outVTE - Methotrexate >= 0,
                          outVTE - Methotrexate,
                          NA)) %>%
  mutate(Methotrexate = 0,
         BIOLOGIC = NA)


temp2 <- dermatologyData %>%
  filter(Test == 2) %>%
  mutate(death_date = ifelse(death_date - BIOLOGIC >= 0, 
                             death_date - BIOLOGIC, 
                             NA),
         emig_date = ifelse(emig_date - BIOLOGIC >= 0, 
                            emig_date - BIOLOGIC,
                            NA),
         eos_date = ifelse(eos_date - BIOLOGIC >= 0,
                           eos_date - BIOLOGIC,
                           NA),
         outVTE = ifelse(outVTE - BIOLOGIC >= 0,
                         outVTE - BIOLOGIC,
                         NA)) %>%
  mutate(BIOLOGIC = 0,
         Methotrexate = NA)


dermatologyData <- dermatologyData %>%
  filter(Test == 0) %>%
  full_join(temp1) %>%
  full_join(temp2)



remove(temp1, temp2)



dermatologyData <- dermatologyData %>%
  mutate(timeType = case_when(death_date == pmin(death_date, emig_date, eos_date, outVTE, na.rm = TRUE) ~ "Death",
                              emig_date == pmin(death_date, emig_date, eos_date, outVTE, na.rm = TRUE) ~ "Emigration",
                              eos_date == pmin(death_date, emig_date, eos_date, outVTE, na.rm = TRUE) ~ "End of study",
                              outVTE == pmin(death_date, emig_date, eos_date, outVTE, na.rm = TRUE) ~ "VTE"))




# Add the row showing if an observation was censored or not

dermatologyData <- dermatologyData %>%
  mutate(status = case_when(timeType != "VTE" ~ 0,
                              timeType == "VTE" ~ 1))





# Reformat the age column to be in years.

dermatologyData <- dermatologyData %>% 
  mutate(age =(-birth_date) %/% 365)


# Add the time column

dermatologyData <- dermatologyData %>%
  mutate(time = case_when(timeType == "Death" ~ death_date,
                           timeType == "Emigration" ~ emig_date,
                           timeType == "End of study" ~ eos_date,
                           timeType == "VTE" ~ outVTE ))




# Change time from days to years.

dermatologyData <- dermatologyData %>%
  mutate(time = time / 365.25)



# Add age groups

dermatologyData <- dermatologyData %>%
  mutate(ageGroup = case_when(age >= 75 ~ '75+',
                              age >= 60 ~ '60-74',
                              age >= 45 ~ '45-59',
                              age >= 0 ~ '0-44'))





# Format the tibble to be easier to read
# Add the duplicate column for any patient that was doubled to make the original data frame easier to use.
# Need to know this so they can be removed for visualizing age for example.

dermatologyData <- dermatologyData %>%
  select(Treatment, time, timeType, status, eos_date, age, ageGroup, AD, sex:History_BIOLOGIC, Test, -Smoking) %>%
  rename(Duplicate = Test) %>%
  mutate(Duplicate = if_else(Duplicate == 2, 1, 0))




# Remove 0 times

dermatologyData <- dermatologyData %>%
  filter(time != 0)



# Write the data to csv file.

write_csv(dermatologyData, file = "data/cleaned_data.csv")





# - - - - - - - - - - - 


# Base hazard dataframe


# - - - - - - - - - - -


# Define a data frame to be used for the base hazard rate.
# To be used later to make graphs.


baseNewdata <- data.frame(VTE = 0,
                          ageGroup = '0-44',
                          Treatment = 'Methotrexate',
                          AD = 0, 
                          AD_Drugs = 0, 
                          Anticoagulant = 0, 
                          Asthma = 0, 
                          Atrial_fibrillitation = 0, 
                          CD = 0, 
                          COPD = 0, 
                          Cancer = 0, 
                          Diabetes = 0, 
                          Fractures = 0, 
                          Heart_failure = 0, 
                          History_BIOLOGIC = 0, 
                          History_Methotrexate = 0, 
                          Hyperlipidemia = 0, 
                          Hypertension = 0, 
                          MACE = 0, 
                          Melanoma_skin_cancer = 0, 
                          Obesity = 0,
                          PSO = 0, 
                          RA = 0, 
                          Rhinitis = 0, 
                          Skin_cancer = 0, 
                          Stroke = 0, 
                          UC = 0, 
                          Uveit = 0, 
                          sex = 1 ) %>%
  as_tibble()






# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -




# Formula




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -








# Define the regression formula to be used.

Formula = Surv(time, status) ~ 
  ageGroup +
  Treatment +
  AD + 
  AD_Drugs + 
  Anticoagulant + 
  Asthma + 
  Atrial_fibrillitation + 
  CD + 
  COPD + 
  Cancer + 
  Diabetes + 
  Fractures + 
  Heart_failure + 
  History_BIOLOGIC + 
  History_Methotrexate + 
  Hyperlipidemia + 
  Hypertension + 
  MACE + 
  Melanoma_skin_cancer + 
  Obesity + 
  PSO+ 
  RA+ 
  Rhinitis + 
  Skin_cancer + 
  Stroke + 
  UC + 
  Uveit +
  VTE +
  sex







# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -




# Kaplan-Meier Plots




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



# Make a new data set with some covariates removed to make some things easier to program.



dermatologyData2 <- dermatologyData %>% 
  select(-age, 
         -timeType, 
         -Duplicate, 
         -eos_date,
         -index_date)








# Kaplan-Meier estimate without covariates

dermatologyKM <- survfit2(Surv(time, status) ~ 1, data = dermatologyData) %>% 
  ggsurvfit() +
  labs(
    x = "Years",
    y = "Overall survival probability"
  ) + 
  ylim(0, 1) +
  add_confidence_interval() +
  add_risktable()


# Kaplan-Meier for treatments

dermatologyKMtreatment <- survfit2(Surv(time, status) ~ Treatment, data = dermatologyData) %>% 
  ggsurvfit() +
  labs(
    x = "Years",
    y = "Overall survival probability"
  ) + 
  ylim(0, 1) +
  add_confidence_interval() +
  add_risktable()


# Kaplan-Meier for age group

dermatologyKMageGroup <- survfit2(Surv(time, status) ~ ageGroup, data = dermatologyData) %>% 
  ggsurvfit() +
  labs(
    x = "Years",
    y = "Overall survival probability"
  ) + 
  ylim(0, 1) +
  add_confidence_interval() +
  add_risktable()



# Define a function that takes the string of a covariate (for example "AD") as argument and returns a Kaplan-Meier graph for that variable.


KaplanMeierPlotFunction <- function(Covariate){
  
  
  form <- reformulate(termlabels = Covariate, response = "Surv(time, status)")
  
  
  survfit2(formula = form, data = dermatologyData) %>% 
    ggsurvfit() +
    labs(
      x = "Years",
      y = "Overall survival probability"
    ) + 
    ylim(0, 1) +
    add_confidence_interval() +
    add_risktable()
  
}






# - - - - - - - - - - - - - - 


# Kaplan-Meier graphs Without the risktable (Ended up using these in the thesis since they take up less space).


# Kaplan-Meier without covariates

dermatologyKM2 <- survfit2(Surv(time, status) ~ 1, data = dermatologyData) %>% 
  ggsurvfit() +
  labs(
    x = "Years",
    y = "S(t)"
  ) + 
  ylim(0, 1) +
  add_confidence_interval()


# Kaplan-Meier for treatments

dermatologyKMtreatment2 <- survfit2(Surv(time, status) ~ Treatment, data = dermatologyData) %>% 
  ggsurvfit() +
  labs(
    x = "Years",
    y = "S(t)"
  ) + 
  ylim(0, 1) +
  add_confidence_interval() +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00"))



# Kaplan-Meier for age group

dermatologyKMageGroup2 <- survfit2(Surv(time, status) ~ ageGroup, data = dermatologyData) %>% 
  ggsurvfit() +
  labs(
    x = "Years",
    y = "S(t)"
  ) + 
  ylim(0, 1) +
  add_confidence_interval() +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  scale_fill_manual(values = c("#0072B2", "#D55E00"))




# Kaplan-Meier plotting function.


KaplanMeierPlotFunction2 <- function(Covariate){
  
  
  form <- reformulate(termlabels = Covariate, response = "Surv(time, status)")
  
  
  survfit2(formula = form, data = dermatologyData) %>% 
    ggsurvfit() +
    labs(
      x = "Years",
      y = "S(t)"
    ) + 
    ylim(0, 1) +
    add_confidence_interval() +
    add_legend_title() +
    scale_color_manual(values = c("#0072B2", "#D55E00")) +
    scale_fill_manual(values = c("#0072B2", "#D55E00"))
  
}





#- - - - - - - - - - - - - - - - -



# Kaplan-Meier plots using the base plot function instead of ggplot.





KaplanMeierPlotFunction3 <- function(Covariate){
  
  
  form <- reformulate(termlabels = Covariate, response = "Surv(time, status)")
  
  plot(survfit2(formula = form, data = dermatologyData),
       col = c("#0072B2", "#D55E00"),
       xlab = "Years",
       ylab = "Overall survival")
  legend("topleft", legend = c("Base Hazard", Covariate), lty = 1, col = c("#0072B2", "#D55E00"), bty = "n")
}



AgeGroupKM <- function(){
  plot(survfit2(Surv(time, status) ~ ageGroup, data = dermatologyData),
       col = c("#0072B2", "#D55E00", "#F0E442", "#009E73"),
       xlab = "Years",
       ylab = "Overall survival")
  legend("topleft", legend = c("0-44", "45-59", "60-75", "75+" ), lty = 1, col = c("#0072B2", "#D55E00", "#F0E442", "#009E73"), bty = "n")
}







# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -











# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -











# Cox model










# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# Define a new data set with some covariates removed to make programming easier.

dermatologyData2 <- dermatologyData %>% 
  select(-age, 
         -timeType, 
         -Duplicate, 
         -eos_date,
         -index_date)




# Create a Cox model for the data with all covariates included.

coxDermatology <- coxph(Surv(time, status) ~ ., data = dermatologyData2)




# Schoenfeld test

# Perform a Schoenfeld test for our Cox regression model.

SchoenfeldDermatology <- cox.zph(coxDermatology)









# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -








# Parametric flexible model







# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -




# Define the formula to be used

Formula = Surv(time, status) ~ 
  ageGroup +
  Treatment +
  AD + 
  AD_Drugs + 
  Anticoagulant + 
  Asthma + 
  Atrial_fibrillitation + 
  CD + 
  COPD + 
  Cancer + 
  Diabetes + 
  Fractures + 
  Heart_failure + 
  History_BIOLOGIC + 
  History_Methotrexate + 
  Hyperlipidemia + 
  Hypertension + 
  MACE + 
  Melanoma_skin_cancer + 
  Obesity + 
  PSO+ 
  RA+ 
  Rhinitis + 
  Skin_cancer + 
  Stroke + 
  UC + 
  Uveit +
  VTE +
  sex




# Create the model

# The model takes a relatively long time to run so it has been written as a comment to make sure it doesnt run
# every time the pdf has to be knitted. The model is saved to a file in a later step and then loaded to save time.


#flexibleDermatology <- stpm2(formula = Formula, 
#                             # Data to be used
#                             data = dermatologyData,
#                             # List the Covariates that should be time dependent.
#                             tvc = list(VTE = 5),
#                             smooth.formula = ~nsx(log(time), df = 5), 
#                             # Degrees of freedom translates to number of knots.
#                             df = 5)





# Save the model to a file.


#saveRDS(flexibleDermatology, "data/flex_dermatology.rds")


# Then load it.


flexibleDermatology <- readRDS("data/flex_dermatology.rds")



# A function that takes a covariate in string form to create a graph for the estimated hazard curve.


HazardPlotFunction <- function(Covariate, limits = c(0,0.04)){
  
  # The "0-vector" to be used for the base hazard rate.
  baseNewdata <- data.frame(VTE = 0,
                            ageGroup = '0-44',
                            Treatment = 'Methotrexate',
                            AD = 0, 
                            AD_Drugs = 0, 
                            Anticoagulant = 0, 
                            Asthma = 0, 
                            Atrial_fibrillitation = 0, 
                            CD = 0, 
                            COPD = 0, 
                            Cancer = 0, 
                            Diabetes = 0, 
                            Fractures = 0, 
                            Heart_failure = 0, 
                            History_BIOLOGIC = 0, 
                            History_Methotrexate = 0, 
                            Hyperlipidemia = 0, 
                            Hypertension = 0, 
                            MACE = 0, 
                            Melanoma_skin_cancer = 0, 
                            Obesity = 0,
                            PSO = 0, 
                            RA = 0, 
                            Rhinitis = 0, 
                            Skin_cancer = 0, 
                            Stroke = 0, 
                            UC = 0, 
                            Uveit = 0, 
                            sex = 1 ) %>%
    as_tibble()
  
  

  
  newdata2 <- baseNewdata %>% 
    mutate({{Covariate}} := 1)
  
  
  plot(flexibleDermatology,
       newdata = baseNewdata,
       xlab = "Years",
       ylab = "Hazard",
       type = "hazard",
       ci = FALSE,
       ylim = limits,
       rug = FALSE)
  lines(flexibleDermatology, 
        newdata = newdata2, 
        col=2, 
        lty=2, 
        type="hazard")
  legend("topright", legend = c("Base Hazard", str_c(Covariate)), lty = 1:2, bty = "n")
}




# Some covariates need to have different arguments.

# Treatment Hazard


TreatmentHazard <-  function(){
  plot(flexibleDermatology,
       newdata = baseNewdata,
       xlab = "Years",
       ylab = "Hazard",
       type = "hazard",
       ci = FALSE,
       ylim = c(0, 0.04),
       rug = FALSE)
  lines(flexibleDermatology, newdata = baseNewdata %>% mutate(Treatment = "Biological"), col=2, lty=2, type="hazard")
  legend("topright", legend = c("Base Hazard", "Methotrexate"), lty = 1:2, bty = "n")
}



# Age group hazard


AgegroupHazard <- function(){
plot(flexibleDermatology,
     newdata = baseNewdata,
     xlab = "Years",
     ylab = "Hazard",
     type = "hazard",
     ci = FALSE,
     ylim = c(0, 0.04),
     rug = FALSE)
lines(flexibleDermatology, newdata = baseNewdata %>% mutate(ageGroup = "45-59"), col=2, lty=2, type="hazard")
lines(flexibleDermatology, newdata = baseNewdata %>% mutate(ageGroup = "60-74"), col=2, lty=2, type="hazard")
lines(flexibleDermatology, newdata = baseNewdata %>% mutate(ageGroup = "75+"), col=2, lty=2, type="hazard")
legend("topright", legend = c("Base Hazard", "45-59", "60-74", "75+"), lty = 1:4, bty = "n")
}



# VTE Hazard

HistoryOfVTEHazard <-function(){
plot(flexibleDermatology,
     newdata = baseNewdata,
     xlab = "Years",
     ylab = "Hazard",
     type = "haz",
     ci = FALSE,
     ylim = c(0, 0.04),
     rug = FALSE)
lines(flexibleDermatology, newdata = baseNewdata %>% mutate(VTE  = 1), col=2, lty=2, type="hazard")
legend("topright", legend = c("Base Hazard", "History of VTE"), lty = 1:2, bty = "n")
}


# Sex Hazard

SexHazard <-function(){
  plot(flexibleDermatology,
       newdata = baseNewdata,
       xlab = "Years",
       ylab = "Hazard",
       type = "haz",
       ci = FALSE,
       ylim = c(0, 0.04),
       rug = FALSE)
  lines(flexibleDermatology, newdata = baseNewdata %>% mutate(sex  = 2), col=2, lty=2, type="hazard")
  legend("topright", legend = c("Base Hazard", "Female hazard rate"), lty = 1:2, bty = "n")
}




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



# Color-blind friendly colors.


cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")




# A function that takes the name of a covariate in string form to create a log-minus-log survival curve for that covariate.


loglogPlotFunction2 <- function(Covariate){

  
  form <- reformulate(termlabels = Covariate, response = "Surv(time, status)")
  
  plot(survfit2(formula = form, data = dermatologyData),
       col = c("#0072B2", "#D55E00"),
       fun = "cloglog",
       xlab = "log(years)",
       ylab = "log-minus-log survival")
  legend("topleft", legend = c("Base Hazard", Covariate), lty = 1, col = c("#0072B2", "#D55E00"), bty = "n")
}


# Age Group has four levels and need to be defined differently from the others.


loglogAgeGroup2 <- function(){
  plot(survfit2(Surv(time, status) ~ ageGroup, data = dermatologyData),
       col = c("#0072B2", "#D55E00", "#F0E442", "#009E73"),
       fun = "cloglog",
       xlab = "log(years)",
       ylab = "log-minus-log survival")
  legend("topleft", legend = c("0-44", "45-59", "60-75", "75+" ), lty = 1, col = c("#0072B2", "#D55E00", "#F0E442", "#009E73"), bty = "n")
}



