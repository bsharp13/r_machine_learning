library(tidyverse)
library(caret)
library(ranger)

#-------------------------------------------------------------------------------

# Introcution
#-------------------------------------------------------------------------------
# When it comes to machine learning in R, no tool is more useful than Caret's 
# train() function. In my time using it, I've encountered very few problems, but 
# I've recently run in to some difficulties creating probability trees using the 
# ranger method. The simplest viable solution to this issue was to continue using 
# pure classification trees and aggregate the predictions of the individual 
# decision trees. The mean response value would serve as a probability estimate 
# (e.g., if 450 of the 500 predictions for a particular observation were "1", the
# probability for that observation would be 90%).

# Such predictions would not be pure probability trees, but I had had enough of 
# trying to force probability trees. Rather than blindly pursuing my hunch, 
# however, I decided to run this simulation in an attempt to measure the 
# difference in predictions using these two methods. Each simulation fits two 
# models (one probability forest and one classification forest) using the ranger 
# method on a set of newly generated data. The function is designed to output 
# the mean, median, or max difference in these predictive methods.


#-------------------------------------------------------------------------------


# Simulation
#-------------------------------------------------------------------------------

run <- function(seed, output = 'mean') {
  
  # Generate data
  set.seed(seed)
  df <- data_frame(
    x1 = rnorm(1000),
    x2 = rnorm(1000),
    x3 = rnorm(1000),
    x4 = sample(0:1, 1000, prob = c(0.3, 0.7), replace = TRUE)
  ) %>% mutate(
    y = as.numeric(
      round(0.5 + (x1 / 3) + (x2 / 3) - (0.2 * x4) + (x1 * x4)) <= 0
    )
  )
  
  # Partition Data
  set.seed(seed)
  train_row <- createDataPartition(df$y, p = 0.8, list = FALSE)
  training <- df[train_row, ]
  testing <- df[-train_row, ]
  
  # Model fitting
  set.seed(seed)
  prob_trees <- ranger(
    y ~ .,
    data = training,
    num.trees = 200,
    mtry = 2,
    probability = TRUE
  )
  
  binary_trees <- ranger(
    y ~ .,
    data = training,
    num.trees = 200,
    mtry = 2
  )
  
  # Model assessment
  prob_tree_pred <- predict(prob_trees, testing)$predictions[,1]
  binary_tree_pred <- 
    rowMeans(predict(binary_trees, testing, predict.all = TRUE)$predictions)
  
  # Calculate output metric difference in predictions
  if (output == 'mean') {
    result <- mean(prob_tree_pred - binary_tree_pred)
  } else if (output == 'median') {
    result <- median(prob_tree_pred - binary_tree_pred)
  } else if (output == 'max') {
    result <- max(prob_tree_pred - binary_tree_pred)
  }
  
  return(result)
  
}

# Run simulations
n_reps <- 1000
mean_diff <- sapply(1:n_reps, function(i) run(i, output = 'mean'))
mean_diff <- sapply(1:n_reps, function(i) run(i, output = 'median'))
max_diff <- sapply(1:n_reps, function(i) run(i, output = 'max'))






