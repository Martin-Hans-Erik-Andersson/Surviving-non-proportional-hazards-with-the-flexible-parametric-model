library(tidyverse)
library(knitr)
library(haven)
library(dplyr)
library(survival)
library(flexsurv)
library(rstpm2)
library(ggfortify)
library(broom)







# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# Data cleaning


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



# Loading and cleaning the data set.


strippedData <- read_sas("data/strip.sas7bdat")




strippedData <- read_sas("data/strip.sas7bdat") %>%
  rename(BIOLOGIC = BIOLOGIC1, 
         Methotrexate = Methotrexate1, 
         death_date = death_date1, 
         emig_date = emig_date1, 
         outVTE = outVTE1,
         eos_date = eos_date1)




# Duplicating any person who received both methods of treatment.
BaseSimulationData <- strippedData %>%
  mutate(Test = if_else(Methotrexate %>% is.na() | BIOLOGIC %>% is.na(), 0, 1 ))


BaseSimulationData <- rbind(BaseSimulationData,
                         BaseSimulationData %>% 
                           filter(Test == 1) %>%
                           mutate(Test = 2))






# Add which treatment each observation received.

BaseSimulationData <- BaseSimulationData %>%
  mutate(Treatment = case_when(Test == 2 ~ "Biological",
                               
                               
                               Test == 1 ~ "Methotrexate",
                               
                               
                               Test == 0 ~ if_else(Methotrexate %>% is.na(), 
                                                   "Biological", 
                                                   "Methotrexate"))
  )





# Add what type of time each observation has (Death, emigration, end of study or VTE).


# Format the data for simplicity


temp1 <- BaseSimulationData %>%
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


temp2 <- BaseSimulationData %>%
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


BaseSimulationData <- BaseSimulationData %>%
  filter(Test == 0) %>%
  full_join(temp1) %>%
  full_join(temp2)



remove(temp1, temp2)



BaseSimulationData <- BaseSimulationData %>%
  mutate(timeType = case_when(death_date == pmin(death_date, emig_date, eos_date, outVTE, na.rm = TRUE) ~ "Death",
                              emig_date == pmin(death_date, emig_date, eos_date, outVTE, na.rm = TRUE) ~ "Emigration",
                              eos_date == pmin(death_date, emig_date, eos_date, outVTE, na.rm = TRUE) ~ "End of study",
                              outVTE == pmin(death_date, emig_date, eos_date, outVTE, na.rm = TRUE) ~ "VTE"))




# Add the row showing if an observation was censored or not

BaseSimulationData <- BaseSimulationData %>%
  mutate(status = case_when(timeType != "VTE" ~ 0,
                            timeType == "VTE" ~ 1))





# Reformat the age column to be in years.

#BaseSimulationData <- BaseSimulationData %>% 
#  mutate(age =(-birth_date) %/% 365)


# Add the time column

BaseSimulationData <- BaseSimulationData %>%
  mutate(time = case_when(timeType == "Death" ~ death_date,
                          timeType == "Emigration" ~ emig_date,
                          timeType == "End of study" ~ eos_date,
                          timeType == "VTE" ~ outVTE ))




# Change time from days to years.

BaseSimulationData <- BaseSimulationData %>%
  mutate(time = time / 365.25,
         eos_date = eos_date / 365.25)



# Add age groups

BaseSimulationData <- BaseSimulationData %>%
  mutate(ageGroup = case_when(age >= 75 ~ '75+',
                              age >= 60 ~ '60-74',
                              age >= 45 ~ '45-59',
                              age >= 0 ~ '0-44'))





# Format the tibble to be easier to read
# Add the duplicate column for any patient that was doubled to make the original data frame easier to use.
# Need to know this so they can be removed for visualizing age for example.

BaseSimulationData <- BaseSimulationData %>%
  select(Treatment, time, timeType, status, eos_date, age, ageGroup, Test) %>%
  rename(Duplicate = Test) %>%
  mutate(Duplicate = if_else(Duplicate == 2, 1, 0))




# Remove 0 times

BaseSimulationData <- BaseSimulationData %>%
  filter(time != 0)




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# Proportional hazards simulation


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -




SimulationFunction1 <- function(DATA = BaseSimulationData, N = 500, m = 45000, Rate1, Rate2, AgeFactor){
  
  
  
  
  # Create the dataframe that will contain the final results.
  
  Results <- c(Index = list(1:N), 
               CoxAIC = list(1:N), 
               FlexAIC = list(1:N), 
               CoxBeta = list(1:N), 
               FlexBeta = list(1:N), 
               SchoenfeldP = list(1:N),
               CoxConcordance = list(1:N),
               UncensoredPercent = list(1:N))
  
  # Results$Index[[y]] gives the value the "Index" column and row y.  
  
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  

  
  
  # THis is the function for simulating times.
  
  SimulationFunction <- function(TIBBLE, 
                                 meanM = Rate1, 
                                 meanB = Rate2, 
                                 ageGroupFactor = AgeFactor){
    
    
    TIBBLE %>%
      mutate(ageFactor = case_when(ageGroup == '0-44' ~ 1,
                                   ageGroup == '45-59' ~ ageGroupFactor,
                                   ageGroup == '60-74' ~ ageGroupFactor * 2,
                                   ageGroup == '75+' ~ ageGroupFactor * 3)) %>%
      mutate(meanLambda = case_when(Treatment == "Methotrexate" ~ meanM,
                                    Treatment == "Biological" ~ meanB)) %>%
      mutate(meanLambda = meanLambda *  1 / ageFactor) %>%
      mutate(simTime = rexp(n = n(), rate = 1 / meanLambda)) %>%
      mutate(outputTime = case_when(status == 0 ~ if_else(simTime > time, time, simTime),
                                    status == 1 ~ if_else(simTime > eos_date, eos_date, simTime))) %>%
      mutate(status = case_when(
        outputTime == simTime ~ 1,
        outputTime != simTime ~ 0)) %>%
      filter(outputTime != 0) %>%
      return()
    
    
  }
  

  
  
  
  # Create a for loop that repeats n number of times.
  
  
  for (i in 1:N){
    
    # Create a plasmode dataset by sampling from the original data with 45000 out of 48591 bootstrap with replacement.
    PlasmodeData <- DATA %>%
      slice_sample(n = 45000 , replace = TRUE)
    
    
    NewData <- SimulationFunction(PlasmodeData) %>%
      filter(time != 0)
    
    
    # Fit cox model.
    
    cox <- coxph(Surv(outputTime, status) ~ Treatment + ageGroup, data = NewData)
    
    
    
    
    coxAIC <- AIC(cox)
    
    
    
    # Fit flexible parametric model.
    
    
    fpma <- stpm2(Surv(outputTime, status) ~ Treatment + ageGroup, 
                  data=NewData %>% filter(time != 0) %>% mutate(Treatment = as.factor(Treatment)), df=5)
    
    
    flexAIC <- AIC(fpma)
    
    
    
    
    
    
    
    # Proportionality test
    
    SchoenfeldCox<- cox.zph(cox)
    

    # Insert results into the results tibble.
    
    
    Results$CoxAIC[i] <- coxAIC
    Results$FlexAIC[i] <- flexAIC
    Results$CoxBeta[i] <- cox$coefficients[1] %>% exp()
    Results$FlexBeta[i] <- coefficients(fpma)[2] %>% exp()
    Results$SchoenfeldP[i] <- SchoenfeldCox$table[1,3]
    Results$CoxConcordance[i] <- cox$concordance[6] 
    Results$UncensoredPercent[i] <- NewData %>%
      summarise(mean = mean(status)) %>%
      pull()
  }
  
  
  
  
  
  # Return results tibble
  
  
  Results %>%
    as_tibble()%>%
    return()

}




# Run the simulation 500 times.


set.seed(2024)



ProportionalSimulationResults <- SimulationFunction1(N = 500, Rate1 = 220 * 1.8, Rate2 = 191 * 1.8, AgeFactor =  2)



# Save results


write_csv(ProportionalSimulationResults, "data/SimulationResults1.csv")






# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# Non Proportional hazards simulation


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -








SimulationFunction2 <- function(DATA = BaseSimulationData, N = 500, m = 45000, Rate11, Rate12, Rate21, Rate22, AgeFactor, Mid, DF = 5){
  
  
  
  
  # Create the dataframe that will contain the final results.
  
  Results <- c(Index = list(1:N), 
               CoxAIC = list(1:N), 
               FlexAIC = list(1:N), 
               CoxBeta = list(1:N), 
               FlexBeta = list(1:N), 
               SchoenfeldP = list(1:N),
               CoxConcordance = list(1:N),
               UncensoredPercent = list(1:N))
  
  # Results$Index[[y]] gives the value the "Index" column and row y.  
  
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  
  
  
  # THis is the function for simulating times.
  
  
  
  SimulationFunction <- function(TIBBLE, 
                                  meanM1 = Rate11, 
                                  meanM2 = Rate12, 
                                  meanB1 = Rate21, 
                                  meanB2 = Rate22,
                                  ageGroupFactor = AgeFactor,
                                  Midpoint = Mid){
    
    
    
    TIBBLE %>%
      
  
      
      mutate(ageFactor = case_when(ageGroup == '0-44' ~ 1,
                                   ageGroup == '45-59' ~ ageGroupFactor,
                                   ageGroup == '60-74' ~ ageGroupFactor * 2,
                                   ageGroup == '75+' ~ ageGroupFactor * 3)) %>%
      
      
      mutate(meanLambda1 = case_when(Treatment == "Methotrexate" ~ meanM1,
                                     Treatment == "Biological" ~ meanB1),
             meanLambda2 = if_else(meanLambda1 == meanM1, meanM2, meanB2)) %>%
      
      
      mutate(meanLambda1 = meanLambda1 * 1 / ageFactor,
             meanLambda2 = meanLambda2 * 1 / ageFactor) %>%
      
      
      mutate(simTime1 = rexp(n = n(), rate = 1 / meanLambda1),
             simTime2 = rexp(n = n(), rate = 1 / meanLambda2),
             simTime = if_else(simTime1 < Midpoint, simTime1, simTime2)) %>%
      
      
      mutate(outputTime = case_when(status == 0 ~ if_else(simTime > time, time, simTime),
                                    status == 1 ~ if_else(simTime > eos_date, eos_date, simTime))) %>%
      
      
      mutate(status = case_when(
        outputTime == simTime ~ 1,
        outputTime != simTime ~ 0)) %>%
      filter(outputTime != 0) %>%
      return()
  }
  
  
  
  
  
  # Create a for loop that repeats n number of times.
  
  
  for (i in 1:N){
    
    # Create a plasmode dataset by sampling from the original data with 45000 out of 48591 bootstrap with replacement.
    PlasmodeData <- DATA %>%
      slice_sample(n = 45000 , replace = TRUE)
    
    
    NewData <- SimulationFunction(PlasmodeData) %>%
      filter(time != 0) %>%
      mutate(Treatment = as.factor(Treatment),
             ageGroup = as.factor(ageGroup))
    
    
    # Fit cox model.
    
    cox <- coxph(Surv(outputTime, status) ~ Treatment + ageGroup, data = NewData)
    
    
    
    
    coxAIC <- AIC(cox)
    
    
    
    # Fit flexible parametric model.
    
    
    fpma <- stpm2(Surv(outputTime, status) ~ Treatment + ageGroup, 
                  data=NewData, 
                  df=DF)
    
    fpma2 <- stpm2(Surv(outputTime, status) ~ Treatment + ageGroup, 
                  data=NewData, 
                  tvc = list(Treatment = 5),
                  smooth.formula = ~ nsx(log(time), df = DF),  
                  df=DF)
    
    
    flexAIC <- AIC(fpma2)
    
    
    
    
    
    
    
    # Proportionality test
    
    SchoenfeldCox <- cox.zph(cox)
    
    
    # Insert results into the results tibble.
    
    
    Results$CoxAIC[i] <- coxAIC
    Results$FlexAIC[i] <- flexAIC
    Results$CoxBeta[i] <- cox$coefficients[1] %>% exp()
    Results$FlexBeta[i] <- coefficients(fpma)[2] %>% exp()
    Results$SchoenfeldP[i] <- SchoenfeldCox$table[1,3]
    Results$CoxConcordance[i] <- cox$concordance[6] 
    Results$UncensoredPercent[i] <- NewData %>%
      summarise(mean = mean(status)) %>%
      pull()
  }
  
  
  
  
  
  # Return results tibble
  
  
  Results %>%
    as_tibble()%>%
    return()
  
}



# Run the simulation

set.seed(2024)


constant = 430


NonProportionalSimulationResults <- SimulationFunction2(N = 500, 
                                                        Rate11 = 2 * constant, 
                                                        Rate12 = 1 * constant, 
                                                        Rate21 = 1 * constant, 
                                                        Rate22 = 2 * constant, 
                                                        AgeFactor =  2, 
                                                        Mid = 2,
                                                        DF = 5)



# Save results

write_csv(NonProportionalSimulationResults, "data/SimulationResults2.csv")











