
## PROSTATE CANCER DATA ANALYSIS

 ## SETUP & REPRODUCIBILITY


# set.seed() ensures that any random processes (if used later)
# produce identical results each time the script is run.
# This is essential for reproducibility.

set.seed(123)

##  LOAD LIBRARIES

# Load required packages at the beginning so that any missing
# package errors appear immediately rather than mid-analysis.

library(tidyverse)        # data manipulation + visualisation
library(lubridate)        # date handling
library(broom)            # tidy model outputs
library(pROC)             # ROC curves and AUC
library(gtsummary)        # publication-ready tables
library(pscl)             # pseudo R-squared
library(car)              # VIF (multicollinearity check)
library(ResourceSelection) # Hosmer-Lemeshow test
library(lmtest)           # likelihood ratio test
library(kableExtra)       # formatted tables
library(randomForest)     # random forest 

## DATA IMPORT

# The dataset is imported as a data frame.
# Missing values are explicitly defined to ensure correct handling.

data_raw <- read.csv(
  "prostate_cancer_data_2025.csv",
  stringsAsFactors = FALSE,
  na.strings = c("NA", "", "N/A")
)

# Confirm the file loaded correctly
cat("Rows loaded:", nrow(data_raw), "| Columns:", ncol(data_raw), "\n")


## DATA INSPECTION AND UNDERSTANDING


# Inspect dataset to understand:
# - variable types
# - presence of missing data
# - general distributions

str(data_raw)        # structure and variable types
head(data_raw)       # first few rows
summary(data_raw)    # summary statistics
colSums(is.na(data_raw))  # missing values per variable


# Understanding the dataset prevents errors such as:
# - treating categorical variables as numeric
# - ignoring missing data
# - misinterpreting variable meaning


## DATA CLEANING & PREPARATION

# Prepare the dataset for modelling by:
# - converting variables into correct formats
# - applying inclusion criteria
# - defining reference categories

data_clean <- data_raw %>%
  mutate(
    
    # Convert Sex to factor
    Sex = factor(Sex),
    
    # Ethnicity with White as reference group
    # Justification: largest and most commonly used baseline
    Ethnicity = factor(Ethnicity,
                       levels = c("White","Black","South Asian","East Asian","Other")),
    
    # Smoking status (Never = baseline, no exposure)
    Smoking_Status = factor(Smoking_Status,
                            levels = c("Never","Former","Current")),
    
    # Family history (No = baseline risk)
    Family_History = factor(Family_History,
                            levels = c("No","Yes")),
    
    # Convert binary variables into labelled factors
    Diabetes = factor(Diabetes, levels=c(0,1), labels=c("No","Yes")),
    CVD = factor(CVD, levels=c(0,1), labels=c("No","Yes")),
    
    # Convert date variables into proper date format
    Date_of_Birth = dmy(Date_of_Birth),
    PSA_Test_Date = dmy(PSA_Test_Date),
    Diagnosis_Date = dmy(Diagnosis_Date),
    
    # Convert outcome variable to numeric (0/1)
    Diagnosed = as.integer(as.logical(Diagnosed))
  )

# Apply inclusion criteria
analysis_data <- data_clean %>%
  filter(Sex == "Male") %>%        # PSA only relevant for males
  filter(Age >= 44 & Age <= 94) %>% # screening age group
  filter(!is.na(PSA))              # PSA required for model


# These criteria ensure clinical relevance and avoid bias
# from inappropriate or missing data.


nrow(analysis_data)
mean(analysis_data$Diagnosed)

# Check Data size is suitable for Model

# Outcome Balance (Events vs Non-events)
table(analysis_data$Diagnosed)
prop.table(table(analysis_data$Diagnosed))

# Number of predictors in your model
num_predictors <- length(coef(model3)) - 1  # remove intercept

# Number of events (cancer cases)
events <- sum(analysis_data$Diagnosed == 1)

events

# EPV calculation
EPV <- events / num_predictors
EPV
## DESCRIPTIVE STATISTICS

# Summarise population characteristics and compare
# diagnosed vs non-diagnosed groups.

summary_table <- analysis_data %>%
  select(Age, PSA, BMI, Ethnicity,
         Smoking_Status, Family_History, Diagnosed) %>%
  mutate(Diagnosed = factor(Diagnosed,
                            levels=c(0,1), labels=c("No Cancer","Cancer"))) %>%
  tbl_summary(by = Diagnosed) %>%
  add_p() %>%        # statistical tests between groups
  bold_labels()

summary_table

# --- PSA and diagnosis rates by ethnicity --------------------


# To quantify how PSA levels and diagnosis rates vary by ethnicity before modelling. 

part1_table <- analysis_data %>%
  group_by(Ethnicity) %>%
  summarise(
    N              = n(),
    Median_PSA     = round(median(PSA, na.rm = TRUE), 2),
    IQR_PSA        = paste0(round(quantile(PSA, 0.25), 2),
                            " - ",
                            round(quantile(PSA, 0.75), 2)),
    Mean_PSA       = round(mean(PSA, na.rm = TRUE), 2),
    SD_PSA         = round(sd(PSA, na.rm = TRUE), 2),
    n_Diagnosed    = sum(Diagnosed),
    Diagnosis_Rate = paste0(round(100 * mean(Diagnosed), 1), "%")
  )

kable(part1_table,
      caption = "Table 2. PSA Levels and Diagnosis Rates by Ethnicity",
      col.names = c("Ethnicity", "N", "Median PSA",
                    "IQR PSA", "Mean PSA", "SD PSA",
                    "n Diagnosed", "Diagnosis Rate (%)")) %>%
  kable_styling(bootstrap_options = c("striped", "hover"),
                full_width = FALSE)


# To Identify potential confounders and important predictors.

## EXPLORATORY DATA ANALYSIS (EDA)

# Visualisation helps identify patterns and relationships.
# EDA reveals trends, outliers, and supports model selection.

# PSA vs Age
ggplot(analysis_data, aes(Age, PSA, color = factor(Diagnosed))) +
  geom_point(alpha = 0.3) +
  scale_color_manual(values = c("0" = "green", "1" = "red"),
                     labels = c("0" = "No Cancer", "1" = "Cancer")) +
  labs(title = "PSA vs Age by Diagnosis",
       color = "Diagnosis") +
  theme_minimal()


# Boxplot: PSA by Ethnicity and Diagnosis
ggplot(analysis_data, aes(x = Ethnicity, y = PSA, fill = factor(Diagnosed))) +
  geom_boxplot() +
  scale_fill_manual(values = c("0"="green", "1"="red"), labels = c("No", "Yes")) +
  labs(title="PSA Distribution by Ethnicity and Diagnosis", fill="Diagnosed") +
  theme_minimal()



# PSA distribution by diagnosis 
ggplot(analysis_data, aes(x = factor(Diagnosed), 
                          y = PSA, fill = factor(Diagnosed))) + 
  geom_boxplot() + scale_fill_manual(values=c("0"="green","1"="red"), 
                                     labels=c("No Cancer","Cancer")) + 
  labs(title="PSA Levels by Diagnosis", x="Diagnosis", y="PSA") + 
  theme_minimal()


## LOGISTIC REGRESSION

#	-Appropriate for binary outcomes 
#	-Produces interpretable results (odds ratios) 


# Model 1: PSA Only

# Baseline model using PSA only
# To establishe baseline predictive performance of PSA alone.

model1 <- glm(Diagnosed ~ PSA,
              data=analysis_data, family=binomial)

prob1 <- predict(model1, type="response")

roc1 <- roc(analysis_data$Diagnosed, prob1)
auc(roc1)


# Model 2: Adjusted Model

# Adjusted for key confounders

model2 <- glm(Diagnosed ~ PSA + Age + Family_History + Ethnicity,
              data=analysis_data, family=binomial)

prob2 <- predict(model2, type="response")

roc2 <- roc(analysis_data$Diagnosed, prob2)
auc(roc2)

# Compare models
lrtest(model1, model2)



## MODEL DIAGNOSTICS


# Multicollinearity check
vif(model2)

# Model fit comparison
AIC(model1, model2)

# Calibration check
hoslem.test(analysis_data$Diagnosed, prob2, g=10)

## SENSITIVITY, SPECIFICITY & OPTIMAL PSA THRESHOLD

# These metrics quantify the model's practical diagnostic value.

psa_threshold <- 4.  #Performance at the standard 4.0 ng/mL threshold

pred <- analysis_data$PSA >= psa_threshold
actual <- analysis_data$Diagnosed == 1

with(list(
  TP = sum(pred & actual),
  FP = sum(pred & !actual),
  TN = sum(!pred & !actual),
  FN = sum(!pred & actual)
), round(100 * c(
  Sensitivity = TP/(TP+FN),
  Specificity = TN/(TN+FP),
  PPV = TP/(TP+FP),
  NPV = TN/(TN+FN)
), 1))

# Optimal threshold via Youden Index

# The Youden Index (J = Sensitivity + Specificity - 1)
# identifies the predicted probability threshold that
# maximises the balance between sensitivity and specificity.

coords(roc2, "best",
       ret = c("threshold", "sensitivity", "specificity"),
       quiet = TRUE) |> round(3)


## INTERACTION ANALYSIS

# To Test whether PSA effect differs by ethnicity

model_interaction <- glm(
  Diagnosed ~ PSA * Ethnicity + Age + Family_History,
  data=analysis_data, family=binomial)

# To identify effect modification.

lrtest(model2, model_interaction)

## EXTENDED MODEL (MODEL 3)

# Add additional predictors

model3 <- glm(
  Diagnosed ~ PSA + Age + Family_History + Ethnicity +
    BMI + Smoking_Status + Diabetes + CVD + Townsend_Index,
  data=analysis_data, family=binomial)


prob3 <- predict(model3, type="response")

roc3 <- roc(analysis_data$Diagnosed, prob3)
auc(roc3)



# Compare Model 2 vs Model 3

# To know  best performed model based on statistical and practical criteria.

lrtest(model2, model3)
AIC(model1, model2, model3)


pR2(model2)
pR2(model3)

auc(roc3) - auc(roc2)

# regression table

Reg_table <- tbl_regression(model3, exponentiate=TRUE) %>%
  bold_labels()

Reg_table

## Random Forest Model (Advanced Model)

# Why Random Forest?
# Captures non-linear relationships 
# Handles interactions automatically 
# Often improves prediction accuracy

# Prepare data (VERY IMPORTANT)

rf_data <- analysis_data

rf_data$Diagnosed <- factor(rf_data$Diagnosed, levels = c(0,1),
                            labels = c("No","Yes"))

# Build Random Forest Model
set.seed(123)

rf <- randomForest(
  Diagnosed ~ PSA + Age + Family_History + Ethnicity,
  data = rf_data,
  ntree = 500,
  importance = TRUE
)

# Predict probabilities

rf_prob <- predict(rf, type = "prob")[,2]  # probability of "Yes"

# ROC + AUC (compare with model3)


roc_rf <- roc(rf_data$Diagnosed, rf_prob)
auc(roc_rf)


## MODEL COMPARISON

# Predictions

# Logistic models
prob1 <- predict(model1, type = "response")
prob2 <- predict(model2, type = "response")
prob3 <- predict(model3, type = "response")


roc1 <- roc(analysis_data$Diagnosed, prob1)
roc2 <- roc(analysis_data$Diagnosed, prob2)
roc3 <- roc(analysis_data$Diagnosed, prob3)


# Plot comparison
plot(roc1, col="blue", lwd=2, main="Model Comparison")
lines(roc2, col="red", lwd=2)
lines(roc3, col="darkgreen", lwd=2)
lines(roc_rf, col="yellow", lwd=2)

legend("bottomright",
       legend=c("Model 1 (PSA)",
                "Model 2 (Model2)",
                "Model 3 (Model3)",
                "Model 4 (Random Forest)"),
       col=c("blue","red","darkgreen","yellow"),
       lwd=2)

 
auc(roc1)
auc(roc2)
auc(roc3)
auc(roc_rf)



