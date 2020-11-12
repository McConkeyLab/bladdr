

# Load Packages -----------------------------------------------------------
library(tidyverse)
library(readxl)
library(faraway)
library(broom)
library(car)

# Load Data ---------------------------------------------------------------

train <- read_excel("./b_datasets/training_data.xlsx")

# How's My Driving? -------------------------------------------------------
r_values <- function(model, n) {
    chi_model <- model$null.deviance - model$deviance
    chi_df <- model$df.null - model$df.residual
    chisq_prob <- 1 - pchisq(chi_model, chi_df)

    R2_hl <- chi_model/model$null.deviance
    R2_cs <- 1 - exp((model$deviance - model$null.deviance)/n)
    R2_n <- R2_cs/(1 - (exp(-(log_model$null.deviance/n))))

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
case_diagnositcs <- function(model, train, n, k) {
   train %>% mutate(., predicted_prob = fitted(model),
                                    standard_residuals = rstandard(model),
                                    student_residuals = rstudent(model),
                                    df_beta = dfbeta(model),
                                    df_fits = dffits(model),
                                    leverage = hatvalues(model),
                                    expected_leverage = (k + 1)/n,
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


# Model 1 -----------------------------------------------------------------

log_model <- glm(success_binary ~ nano_conc + protein + salt + tape_conc +
                 conc_200 + dv200 + input_rna + tech + kit,
                 data = train, family = "binomial")
summary(log_model1)

## from summary: intercept, protein, tape_conc, conc_200, dv200
    ## matches what conf_ints says
        ## not kit though?
## AIC == 352.88

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



# Model 2 -----------------------------------------------------------------
log_model2 <- glm(success_binary ~ nano_conc + protein + salt + tape_conc +
                                 dv200 + tech + input_rna + kit,
                         data = train, family = "binomial")
summary(log_model2)
## AIC == 364.27; went up, interesting

r_values(log_model2, nrow(train))
## p < 0.05

OR(log_model2)
conf_ints(log_model2)
quick_plot(log_model2, train, "second model - no conc_200")
vif(log_model2)

## testing linearity of the logit
probs <- predict(log_model2)


## when you run this as a model, you want the new variables to NOT be
## significant -- that would violate the assumption of linearity

test_model <- glm(success_binary ~ nano_conc + protein + salt + tape_conc +
                          dv200 + tech + input_rna + kit + nano_conc_log +
                          protein_log + salt_log + tape_conc_log + dv200_log,
                  data = logit_frame, family = "binomial")
summary(test_model)


## ROC cut off of like 75%






