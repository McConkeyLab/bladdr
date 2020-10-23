
# Load Data ---------------------------------------------------------------
train <- read.table("training_data")


# Load Packages -----------------------------------------------------------
library(tidyverse)
library(faraway)
library(broom)

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

## from summary: intercept, protein, salt, tape_conc, conc_200, dv200, kit
    ## matches what conf_ints says
## AIC == 446.73

r_values(log_model, 386)
## p < 0.05

OR(log_model)

conf_ints(log_model)

quick_plot(log_model, train, "first model")
## its okay


# Model 2 -----------------------------------------------------------------

log_model2 <- glm(success_binary ~ protein + salt + tape_conc + conc_200
                  + dv200 + kit,
                 data = train, family = "binomial")

## AIC == 448.77

conf_ints(log_model2)

quick_plot(log_model2, train, "second model")

train2 <- case_diagnositcs(log_model2, train)





removedata <- train_model2.1$cooks_distance > 1.0
train2.2 <- train[!removedata,]

log_model2.2 <- glm(success_binary ~ kit +  protein + salt,
                    data = train2.2, family = "binomial")

train_model2.2 <- train2.2 %>% mutate(., predicted_prob = fitted(log_model2.2),
                                   standard_residuals = rstandard(log_model2.2),
                                   student_residuals = rstudent(log_model2.2),
                                   df_beta = dfbeta(log_model2.2),
                                   df_fits = dffits(log_model2.2),
                                   leverage = hatvalues(log_model2.2),
                                   cooks_distance = cooks.distance(log_model2.2))

OR2.2 <- exp(log_model2.2$coefficients)
CI2.2 <- exp(confint(log_model2.2))






## lines 5-24 stolen from
## https://github.com/StatQuest/logistic_regression_demo/
##    blob/master/logistic_regression_demo.R

ll.null <- log_model$null.deviance/-2
ll.proposed <- log_model$deviance/-2
## McFadden's Pseudo R^2 = [ LL(Null) - LL(Proposed) ] / LL(Null)
(ll.null - ll.proposed) / ll.null

## The p-value for the R^2
1 - pchisq(2*(ll.proposed - ll.null), df=(length(log_model$coefficients)-1))

predicted.data <- data.frame(
    probability.of.success = log_model2.2$fitted.values,
    success = train2.2$success)

predicted.data <- predicted.data[
    order(predicted.data$probability.of.success, decreasing = FALSE),]
predicted.data$rank <- 1:nrow(predicted.data)

## Lastly, we can plot the predicted probabilities for each sample having
## heart disease and color by whether or not they actually had heart disease
ggplot(data = predicted.data, aes(x = rank, y = probability.of.success)) +
    geom_point(aes(color = success), alpha=1) +
    ylab("Probability of Success") + xlab(element_blank()) +
    labs(title = "Probability of Library Success using (Incomplete) Training Data")



ggsave("log_modeling.pdf")


