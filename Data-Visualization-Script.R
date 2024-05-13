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

# Ideas

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -




# Colorblind friendly colors to add to graphs.


cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# Colorblind friendly package

# library(viridis)





covariates <- c("ageGroup",
                "Treatment",
                "AD", 
                "AD_Drugs", 
                "Anticoagulant", 
                "Asthma", 
                "Atrial fibrillation", 
                "CD",
                "COPD",
                "Cancer",
                "Diabetes",
                "fractures",
                "heart failure",
                "History BIOLOGIC",
                "History Methotrexate",
                "Hyperlipidemia",
                "Hypertension",
                "MACE",
                "Melanoma skin cancer",
                "Obesity",
                "PSO",
                "RA",
                "Rhinitis",
                "Skin cancer",
                "Stroke",
                "UC",
                "Uveitis",
                "VTE",
                "sex")


description <- c("The persons age group. The original data did not define any age groups. We have chosen to use the age groups 0-44, 45-59, 60-74 and 75+.",
                 'The type of medicine that the patient filled a prescription for. Either "Methotrexate" or "Biological".',
                 "History of atopical dermatitis, also known as eczema, a chronic disease that causes inflammation/redness/irritation of the skin",
                 "History of drugs for the aforementioned atopical dermatitis",
                 "History of anticoagulant drugs which makes your blood less likely to clot.",
                 "History of asthma which causes inflammation in the throat which in turn can cause coughing, wheezing and shortness of breath.",
                 "History of atrial fibrillation which is irregular heart beats.",
                 "History of Crohns disease, a disease that causes inflammation of the bowel.",
                 "History of chronical obstructive pulmonary disease which is a chronic inflammatory lung disease that obstructs airflow from the lungs.",
                 "History of cancer.",
                 "History of diabetes.",
                 "History of bone fractures.",
                 "History of heart failure.",
                 "History of biological medicine before enrolling in the trial.",
                 "History of being treated with methotrexate before enrolling in the trial.",
                 "History of hyperlipidemia (High blood fat).",
                 "History of hypertension (High blood pressure).",
                 "History of MACE (major adverse cardiovascular events) in which medicine causes major problems with the cardiovascular system (Heart and blood)",
                 "History of skin cancer with melanoma, that being cancer that starts in the melanocytes that create the pigment that gives the skin its color.",
                 "History of obesity.",
                 "History of psoriasis which is signified by an overactive immune system that causes skin cells to multiply at an overly fast rate.",
                 "History of reumatoid arthritis which is an autoimmune disease which causes your immune system to attack itself and causing inflammation.",
                 "History of rhinitis which is a reaction in your nose which leads to congestion, runny nose, sneezing and itching.",
                 "History of skin cancer without the aforementioned melanoma.",
                 "History of stroke which is caused by cell death in the brain due to poor blood flow.",
                 "History of ulcerous cholitis which is a chronic condition where the large intestine and rectum become inflamed.",
                 "History of uveitis which is an inflammation in the uvea of the eye.",
                 "History of venous tromboembolism. Venous thromboembolism is either a blood clot that is not close to the surface of the skin or in the arteries connected to the lungs.",
                 "Sex of patient, 1 for male, 2 for female.")


CovariateTable <- list(Covariate = covariates, Description = description) %>%
  as_tibble()



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Illustrating concepts

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# How to illustrate arbitrary function.

base <-
  ggplot() +
  xlim(0, 5)

base + geom_function(fun = ~ 2*.x)

base + 
  geom_function(aes(colour = "Hello"),fun = dexp, args = list(rate = 1)) + 
  geom_function(aes(colour = "Goodbye"), fun = dexp, args = list(rate = 3))

base + geom_function(fun = dexp, args = list(rate = 2))


# Proportional Hazards

prophazardIllustration <- base +
  geom_function(color = "#56B4E9",aes(colour = "Treatment 1"), fun = ~ 1000*exp(-(.x+6))) +
  geom_function(color = "#D55E00",aes(colour = "Treatment 2"), fun = ~ 1000*2*exp(-(.x+6))) +
  geom_hline(yintercept =  0) +
  geom_vline(xintercept =  0) +
  labs(y = "Hazard rate", x = "Years")



# Non proportional hazards

nonprophazardIllustration <- base +
  geom_function(color = "#56B4E9", aes(colour = "Conservative"), fun = ~ 1000*exp(-(.x+5.8)) + 0.2) +
  geom_function(color = "#D55E00", aes(colour = "Aggressive"), fun = ~ 1000*2*exp(-(.x+6))) +
  geom_hline(yintercept =  0) +
  geom_vline(xintercept =  0) +
  labs(y = "Hazard rate", x = "Years")






# Crossing Piecewise Hazards

CrossingPiecewiseHazards <- base +
  geom_function(color = "#56B4E9", aes(colour = "Treatment 1"), fun = ~ ifelse(.x < 2, 1, 2)) +
  geom_function(color = "#D55E00", aes(colour = "Treatment 2"), fun = ~ ifelse(.x < 2, 2, 1)) +
  geom_hline(yintercept =  0) +
  geom_vline(xintercept =  0) +
  labs(y = "Hazard rate", x = "Years")




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Loading data

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -






# Load data that was cleaned in Hazard-script

dermatologyData <- read_csv("data/cleaned_data.csv")









# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# Visualization


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -








# Average number of people with each covariate

MeanTable <- dermatologyData %>% 
  filter(Duplicate == 0) %>%
  select(age:History_BIOLOGIC, -Duplicate) %>%
  mutate(sex = sex - 1) %>%
  rename(femaleRatio = sex) %>%
  summarise(across(where(is.numeric), mean)) %>%
  mutate(across(where(is.numeric), round, 3)) %>%
  pivot_longer(cols = age:History_BIOLOGIC, names_to = "Covariate", values_to = "mean_value") %>%
  mutate("Mean value" =ifelse(Covariate != "age", str_c(mean_value * 100, " %"), mean_value)) %>%
  select("Covariate", "Mean value")



# Ratio of people with each covariate by age group.

AgeGroupRelativeFrequencyTable <- dermatologyData %>%
  filter(Duplicate == 0) %>%
  select(ageGroup:History_BIOLOGIC) %>%
  mutate(sex = sex - 1) %>%
  group_by(ageGroup) %>%
  summarise(across(where(is.numeric), sum),
            n = n()) %>%
  mutate_at(vars(-ageGroup), ~ . / n) %>%
  select(-n) %>%
  ungroup() %>%
  pivot_longer(names_to = "covariate", cols = AD:History_BIOLOGIC, values_to = "relative") %>%
  mutate(relative = relative %>% round(4)) %>%
  mutate(relative = str_c(relative * 100, " %")) %>%
  rename("Relative frequency" = relative) %>%
  pivot_wider(names_from = ageGroup, values_from = 'Relative frequency')
  
  
# Ratio of people on each treatment by age group.

TreatmentByAgegroupTable <- dermatologyData %>%
  group_by(ageGroup, Treatment) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = Treatment, values_from = n) %>%
  mutate(total = Biological + Methotrexate) %>%
  mutate(Biological = Biological / total,
         Methotrexate = Methotrexate / total) %>%
  select(-total) %>%
  mutate(Biological = Biological %>% round(2),
         Methotrexate = Methotrexate %>% round(2)) %>%
  mutate(Biological = str_c(Biological * 100, " %"),
         Methotrexate = str_c(Methotrexate * 100, " %")) %>%
  ungroup() %>%
  pivot_longer(names_to = "Treatment", cols = Biological:Methotrexate) %>%
  pivot_wider(names_from = ageGroup, values_from = value)
  


# Visualizing number of people on each treatment

TreatmentVisualization <- dermatologyData %>%
  select(Treatment) %>%
  group_by(Treatment) %>%
  summarise(n())


# Visualizing ages

AgeVisualization <- dermatologyData %>% 
  filter(Duplicate == 0) %>%
  ggplot() +
  geom_histogram(fill = "#00AFBB", color = "#00AFBB", aes(x = age), binwidth = 2)


# Visualizing survival times

# Time of deaths

DeathTimes <- dermatologyData %>%
  filter(timeType == "Death") %>%
  ggplot() +
  geom_histogram(fill = "#00AFBB", color = "#00AFBB", aes(x = time), binwidth = 0.1) + 
  ylab("Count") +
  xlab("Years")


# End of study times

eosTimes <- dermatologyData %>%
  filter(timeType == "End of study") %>%
  ggplot() +
  geom_histogram(fill = "#00AFBB", color = "#00AFBB", aes(x = time), binwidth = 0.1) + 
  ylab("Count") +
  xlab("Years")


# Emigration times

EmigrationTimes <- dermatologyData %>%
  filter(timeType == "Emigration") %>%
  ggplot() +
  geom_histogram(fill = "#00AFBB", color = "#00AFBB", aes(x = time), binwidth = 0.1) + 
  ylab("Count") +
  xlab("Years")


# Time of outcomes of venous thromoboembolism


VTETimes <- dermatologyData %>%
  filter(timeType == "VTE") %>%
  ggplot() +
  geom_histogram(fill = "#00AFBB",color = "#00AFBB", aes(x = time), binwidth = 0.1) + 
  ylab("Count") +
  xlab("Years")







#Number of people in each age group


agegrouptibble <- dermatologyData %>%
  filter(Duplicate == 0) %>%
  group_by(ageGroup) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(percentage = 100 * n / 48591)



AgeGroupVisualization <- dermatologyData %>%
  ggplot() +
  geom_bar(fill = "#00AFBB", color = "#00AFBB", aes(x = ageGroup))



# Box plot of age groups

AgeBoxPlot <- dermatologyData %>%
  filter(Duplicate == 0) %>%
  ggplot() +
  geom_boxplot(fill = "#00AFBB", aes(y = age, x = ageGroup))




# Number of outcomes of VTE by age group.

a1 <-dermatologyData %>%
  filter(timeType == "VTE",
         Duplicate == 0) %>%
  group_by(ageGroup) %>%
  summarise(n_VTE = n()) %>%
  full_join(agegrouptibble) %>%
  mutate(normalized_n_VTE = 1 / (percentage / 100),
         proportion_VTE = 100 * n_VTE / n)

# Number of outcomes of emigration by age group


a2 <- dermatologyData %>%
  filter(timeType == "Emigration",
         Duplicate == 0) %>%
  group_by(ageGroup) %>%
  summarise(n_emigration = n()) %>%
  full_join(agegrouptibble) %>%
  mutate(normalized_n_emigration = 1 / (percentage / 100),
         proportion_emigration = 100 * n_emigration / n)


# Number of outcomes of death by age group

a3 <- dermatologyData %>%
  filter(timeType == "Death",
         Duplicate == 0) %>%
  group_by(ageGroup) %>%
  summarise(n_death = n()) %>%
  full_join(agegrouptibble) %>%
  mutate(normalized_n_death = 1 / (percentage / 100),
         proportion_death = 100 * n_death / n)

# Number of outcomes of end of study by age group


a4 <- dermatologyData %>%
  filter(timeType == "End of study",
         Duplicate == 0) %>%
  group_by(ageGroup) %>%
  summarise(n_eos = n()) %>%
  full_join(agegrouptibble) %>%
  mutate(normalized_n_eos = 1 / (percentage / 100),
         proportion_eos = 100 * n_eos / n)


# Merge the previous tibbles and remove temp variables.

AgeGroupOutcomes <- a1 %>%
  full_join(a2) %>%
  full_join(a3) %>%
  full_join(a4)



remove(a1, a2, a3, a4)




# Clean the merged tibble

AgeGroupOutcomes <- AgeGroupOutcomes %>%
  select(ageGroup, n, proportion_VTE, proportion_emigration, proportion_death, proportion_eos) %>%
  mutate(across(3:6, round, 2)) %>%
  mutate(proportion_VTE = str_c(proportion_VTE, " %"),
         proportion_emigration = str_c(proportion_emigration, " %"),
         proportion_death = str_c(proportion_death, " %"),
         proportion_eos = str_c(proportion_eos, " %"))






# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# Simulation


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



SimulationResults1 <- read_csv("data/SimulationResults1.csv")

SimulationResults2 <- read_csv("data/SimulationResults2.csv")






# Create table with mean and 95% confidence interval for the proportional simulation results.


a1 <- SimulationResults1 %>%
  summarise(mean = mean(CoxConcordance),
            lowQuantile = quantile(CoxConcordance, 0.025),
            highQuantile = quantile(CoxConcordance, 0.975)) %>%
  mutate(Type = "Cox Concordance")

a2 <- SimulationResults1 %>%
  summarise(mean = mean(SchoenfeldP),
            lowQuantile = quantile(SchoenfeldP, 0.025),
            highQuantile = quantile(SchoenfeldP, 0.975)) %>%
  mutate(Type = "Schoenfeld p value")

a3 <- SimulationResults1 %>%
  summarise(mean = mean(CoxAIC),
            lowQuantile = quantile(CoxAIC, 0.025),
            highQuantile = quantile(CoxAIC, 0.975)) %>%
  mutate(Type = "Cox AIC")

a4 <- SimulationResults1 %>%
  summarise(mean = mean(FlexAIC),
            lowQuantile = quantile(FlexAIC, 0.025),
            highQuantile = quantile(FlexAIC, 0.975)) %>%
  mutate(Type = "FPM AIC")


a5 <- SimulationResults1 %>%
  mutate(difference = CoxBeta - FlexBeta) %>%
  summarise(mean = mean(difference),
            lowQuantile = quantile(difference, 0.025),
            highQuantile = quantile(difference, 0.975)) %>%
  mutate(Type = "Difference in beta estimate")

SimulationMeanProportional <- full_join(a1, a2) %>%
  full_join(a3) %>%
  full_join(a4) %>%
  full_join(a5) %>%
  relocate(Type, lowQuantile, mean, highQuantile)

remove(a1, a2, a3, a4, a5)





# Create the same table but for the non proportional simulation



a1 <- SimulationResults2 %>%
  summarise(mean = mean(CoxConcordance),
            lowQuantile = quantile(CoxConcordance, 0.025),
            highQuantile = quantile(CoxConcordance, 0.975)) %>%
  mutate(Type = "Cox Concordance")

a2 <- SimulationResults2 %>%
  summarise(mean = mean(SchoenfeldP),
            lowQuantile = quantile(SchoenfeldP, 0.025),
            highQuantile = quantile(SchoenfeldP, 0.975)) %>%
  mutate(Type = "Schoenfeld p value")

a3 <- SimulationResults2 %>%
  summarise(mean = mean(CoxAIC),
            lowQuantile = quantile(CoxAIC, 0.025),
            highQuantile = quantile(CoxAIC, 0.975)) %>%
  mutate(Type = "Cox AIC")

a4 <- SimulationResults2 %>%
  summarise(mean = mean(FlexAIC),
            lowQuantile = quantile(FlexAIC, 0.025),
            highQuantile = quantile(FlexAIC, 0.975)) %>%
  mutate(Type = "FPM AIC")


a5 <- SimulationResults2 %>%
  mutate(difference = CoxBeta - FlexBeta) %>%
  summarise(mean = mean(difference),
            lowQuantile = quantile(difference, 0.025),
            highQuantile = quantile(difference, 0.975)) %>%
  mutate(Type = "Difference in beta estimate")

SimulationMeanNonProportional <- full_join(a1, a2) %>%
  full_join(a3) %>%
  full_join(a4) %>%
  full_join(a5) %>%
  relocate(Type, lowQuantile, mean, highQuantile)




# Make tables with probability of type-I and type-II errors for the simulation.


a <- SimulationResults1 %>%
  select(SchoenfeldP) %>%
  mutate(FivePercentSignificance = if_else(SchoenfeldP <= 0.05, 1, 0),
         OnePercentSignificance = if_else(SchoenfeldP <= 0.01, 1, 0)) %>%
  summarise(FivePercentSignificance = mean(FivePercentSignificance),
            OnePercentSignificance = mean(OnePercentSignificance)) %>%
  mutate('Probability of type I error (5% significance)' = str_c(FivePercentSignificance * 100, " %"),
         'Probability of type I error (1% significance)' = str_c(OnePercentSignificance* 100, " %")) %>%
  mutate(Hazards = "Proportional")


SimulationSchoenfeldTable <-  SimulationResults2 %>%
  select(SchoenfeldP) %>%
  mutate(FivePercentSignificance = if_else(SchoenfeldP <= 0.05, 1, 0),
         OnePercentSignificance = if_else(SchoenfeldP <= 0.01, 1, 0)) %>%
  summarise(FivePercentSignificance = mean(FivePercentSignificance),
            OnePercentSignificance = mean(OnePercentSignificance)) %>%
  mutate(Hazards = "Non-Proportional") %>%  
  mutate('Probability of type II error (5% significance)' = str_c((1 - FivePercentSignificance) * 100, " %"),
         'Probability of type II error (1% significance)' = str_c((1 - OnePercentSignificance) * 100, " %")) %>%
  full_join(a) %>%
  select(-FivePercentSignificance, -OnePercentSignificance) %>%
  relocate(Hazards, 'Probability of type I error (5% significance)', 'Probability of type I error (1% significance)')

remove(a)
