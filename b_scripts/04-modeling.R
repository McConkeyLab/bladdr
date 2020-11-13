
# Load Packages -----------------------------------------------------------
library(tidyverse)
library(readxl)
library(faraway)
library(broom)
library(car)

# Load Data ---------------------------------------------------------------

train <- read_excel("./b_datasets/training_data.xlsx")

# How's My Driving? -------------------------------------------------------
## bunch of functions to make comparisons between models easier and faster
r_values <- function(model, n) {
    chi_model <- model$null.deviance - model$deviance
    chi_df <- model$df.null - model$df.residual
    chisq_prob <- 1 - pchisq(chi_model, chi_df)

    R2_hl <- chi_model/model$null.deviance
    R2_cs <- 1 - exp((model$deviance - model$null.deviance)/n)
    R2_n <- R2_cs/(1 - (exp(-(model$null.deviance/n))))

    df <- rbind(R2_hl, R2_cs, R2_n, chisq_prob)

    rownames(df) <- c("Hosmer and Lemeshow’s R^2",
                      "Cox and Snell's R^2",
                      "Nagelkerke’s R^2",
                      "Chi-Squared Statistic")
    colnames(df) <- "Value"
    df
}

OR <- function(model) {
    exp(model$coefficients)
}

conf_ints <- function(model) {
    df <- data.frame(exp(confint(model)))
    names <- rownames(df)
    colnames(df) <- c("a", "b")
    find_t <- apply(df, 1, function(df) between(1, df["a"], df["b"]))
    df <- df %>% `colnames<-`(., c("2.5%", "97.5%")) %>%
        mutate(., Usable = case_when(as.vector(find_t) == F ~ "Y",
                                     TRUE ~ "N")) %>%
        `rownames<-`(., names)
    df
}

quick_plot <- function(model, train, title) {
    predicted_data <- data.frame(
        prob_success = model$fitted.values,
        success = train$success)

    predicted_data <- predicted_data[
        order(predicted_data$prob_success, decreasing = FALSE),]
    predicted_data$rank <- 1:nrow(predicted_data)

    ggplot(data = predicted_data, aes(x = rank, y = prob_success)) +
        geom_point(aes(color = success), alpha=1) +
        ylab("Probability of Success") + xlab(element_blank()) +
        labs(title = title)
}

## this function helps you find problem samples within the training set
## n = number of samples, k = number of predictors
case_diagnositcs <- function(model, train, k) {
   train %>% mutate(., predicted_prob = fitted(model),
                                    standard_residuals = rstandard(model),
                                    student_residuals = rstudent(model),
                                    df_beta = dfbeta(model),
                                    df_fits = dffits(model),
                                    leverage = hatvalues(model),
                                    expected_leverage = (k + 1)/nrow(train),
                                    cooks_distance = cooks.distance(model))
}

## dont wan't standardized residuals to be abs(>3), 1% of samples to have
## stand. resid > 2.5, or more than 5% to be > 2


## does NOT find cases based off leverage
problem_samples <- function(train_plus) {
    std_residuals <- abs(train_plus$standard_residuals) > 2.0
    stud_residuals <- abs(train_plus$student_residuals) > 2.0
    df_beta <- abs(train_plus$df_beta) > 1.0
        ## difference in model when each case is either included or not
    df_fit <- abs(train_plus$df_fits) > 1.0
        ## difference in predicted value for case when that case is included
        ## vs not
    cooks <- train_plus$cooks_distance > 1.0
        ## measure of overall influence of sample on model
    find_t <- std_residuals | stud_residuals | df_beta | df_fit |
        cooks
    find_t <- find_t[,1]

    train_plus[find_t,]

}

## to note with this function - model1 needs to be the model with more
## variables. P < 0.05 means second model is significant improvement
## over first
compare_models <- function(model1, model2) {
    model_chi <- model1$deviance - model2$deviance
    chi_df <- model1$df.residual - model2$df.residual
    chisq_prob <- 1 - pchisq(model_chi, chi_df)
    df <- rbind(model_chi, chi_df, chisq_prob)
    rownames(df) <- c("Model Chi", "Chi DF Residuals", "Chi P Value")
    df

}

log_data <- function(model, train, prob) {
    data <- select_if(train, is.numeric)
    variables <- colnames(data)
    data <- data %>%
        mutate(.,logit = log(prob/(1 - prob))) %>%
        pivot_longer(!logit, names_to = "variables")
}

log_view <- function(log_data) {
    ggplot(log_data, aes(logit, value)) +
               geom_point(size = 0.5, alpha = 0.5) +
               geom_smooth(method = "loess") +
               facet_wrap(~variables, scales = "free_y")
}

pred_values <- function(model, train, k, p) {
    case_diagnositcs(model, train, k = 9) %>%
        mutate(., probability = predict(model, type = "response")) %>%
        select(., colnames(train), probability) %>%
        mutate(., call = ifelse(probability > p, 1, 0))
}

pred_value_test <- function(model, testdata, p) {
    testdata %>%
        mutate(., probability = predict(model, newdata = testdata, type = "response")) %>%
        select(., colnames(testdata), probability) %>%
        mutate(., call = ifelse(probability > p, 1, 0))
}
## most accurate when prob set to 0.15???


# Model 1 -----------------------------------------------------------------

log_model <- glm(success_binary ~ nano_conc + protein + salt + tape_conc +
                 conc_200 + dv200 + input_rna + tech + kit,
                 data = train, family = "binomial")
summary(log_model)

## from summary: intercept, protein, tape_conc, conc_200, dv200
    ## matches what conf_ints says
        ## not kit though?
## AIC == 396.14


# Model 1 - QC ------------------------------------------------------------

r_values(log_model, nrow(train))
## p < 0.05

OR(log_model)
conf_ints(log_model)
        ## if conf int. crosses 1 then its useless --> cant use that variable?

quick_plot(log_model, train, "first model")
## its okay

## checking for multicollinearity(idk how to spell it)
vif(log_model)
## tape_conc and conc_200 (to be expected, theyre pretty much the same)
        ## deciding to keep tape_conc because thats the value we use to
        ## determine rna dilutions


problem <- case_diagnositcs(log_model, train, k = 9) %>% problem_samples(.)

train_add <- pred_values(log_model, train, 9, 0.5)

sum_log <- train_add %>%
    group_by(., success, call) %>%
    summarise(., n())

find_wrong <- train_add$success_binary == train_add$call

bad_calls <- train_add[!find_wrong,]

## 24% of calls are wrong (using 0.5 cut off)

probs <- predict(log_model, type = "response")

log_view(log_data(log_model, train, probs))



# Model 2 -----------------------------------------------------------------
## basically model 1, but with bad variables removed
log_model2 <- glm(success_binary ~ nano_conc + protein + salt + tape_conc +
                                 dv200 + kit + input_rna,
                         data = train, family = "binomial")
summary(log_model2)
## AIC == 401.42; went up, interesting

r_values(log_model2, nrow(train))
## p < 0.05

OR(log_model2)
conf_ints(log_model2)
quick_plot(log_model2, train, "second model - no conc_200")
vif(log_model2)

## testing linearity of the logit
probs2 <- predict(log_model2, type = "response")

log_view(log_data(log_model2, train, probs2))

problem2 <- case_diagnositcs(log_model2, train, k = 7) %>% problem_samples(.)

train_add2 <- pred_values(log_model2, train, 7, 0.5)

sum_log2 <- train_add2 %>%
    group_by(., success, call) %>%
    summarise(., n())

find_wrong2 <- train_add2$success_binary == train_add2$call

bad_calls2 <- train_add2[!find_wrong2,]
##27% miscalled



# Testing -----------------------------------------------------------------
binary_complete <- read_excel("./b_datasets/binary_complete.xlsx")

test_add <- pred_value_test(log_model2, binary_complete, 0.15)

test_sum_log <- test_add %>%
    group_by(., success, call) %>%
    summarise(., n())

test_find_wrong <- test_add$success_binary == test_add$call
table(test_find_wrong) %>% .[2]/(nrow(binary_complete)) * 100
    ## gives % of correct calls

test_bad_calls <- test_add[!test_find_wrong,]

