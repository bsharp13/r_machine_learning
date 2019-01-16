
library(magrittr)
library(tidyverse)
library(gridExtra)
library(caret)

#-------------------------------------------------------------------------------

# Preparation
#-------------------------------------------------------------------------------

# Read in data
mouse <- read_csv('~/Downloads/Data_Cortex_Nuclear.csv')

X <- mouse %>% select(DYRK1A_N:CaNA_N)
y <- as.factor(as.numeric(mouse$Genotype != 'Control'))

#-------------------------------------------------------------------------------

# EDA
#-------------------------------------------------------------------------------
# Final class distribution
table(y)

# Create custom histogram function to plot X variables
my_histogram <- function(data, var, bins = 25) {
  
  width <- (max(data[, var], na.rm = T) - min(data[, var], na.rm = T)) / bins
  ggplot(data, aes_string(x = var)) +
    geom_histogram(binwidth = width, fill = '#099E6A', color = '#FFFFFF') +
    ggtitle(var) +
    theme_bw()
  
}

# Plot X values to identify strange distributions; most are approximately
# normally or lognormally distributed. GluR4_N may have outliers. Note warnings
# about missing values
for (i in seq_along(X)) print(my_histogram(X, colnames(X)[i]))
rm(i)
#-------------------------------------------------------------------------------

# Preprocessing & Test/Train Split
#-------------------------------------------------------------------------------
# Check for missing data in X
colSums(is.na(X))

# Many variables have 3 missing values; check row sums to see if any rows are 
# missing observations for more than half of the variables in the data
remove <- which(rowSums(is.na(X)) > (0.5 * ncol(X)))

# Discard those rows
X <- X[-remove, ]
y <- y[-remove]


# Observe missing-rates in other rows; no other rows are missing values in more
# than five variables, so imputation for the remainder of missing data should
# be fine
table(rowSums(is.na(X)))

# Identify variables with more than 3% of values missing. Perform chi-square
# tests to determine whether those missing values are significantly related
# to the outcome variable
missing <- colMeans(is.na(X))[colMeans(is.na(X)) > 0.03]
remove <- rep(FALSE, length(missing))

for (i in 1:length(missing)) {
  test <- chisq.test(is.na(X[, names(missing)[i]]), y)
  if (test$p.value < 0.05) remove[i] <- TRUE
}

# Remove variables with missingness significantly related to the outcome
# variable. We could keep these as potential predictors, but we will be running
# PCA on these data, so categorical predictors cannot be included
X %<>% select(-c(names(missing)[remove]))

# Check for outliers 
check <- X %>% select(GluR4_N)
no_outlier <- X %>% select(GluR4_N) %>% filter(GluR4_N < 0.2)

plot1 <- check %>%
  ggplot(aes(y = GluR4_N)) +
  geom_boxplot()

plot2 <- no_outlier %>% 
  ggplot(aes(y = GluR4_N)) +
  geom_boxplot()

plot3 <- check %>% 
  my_histogram('GluR4_N')

plot4 <- no_outlier %>% 
  my_histogram('GluR4_N')

grid.arrange(
  arrangeGrob(plot1, plot3, top = 'With Outliers'), 
  arrangeGrob(plot2, plot4, top = 'Without Outliers'),
  ncol = 2, 
  nrow = 1
)

# Outliers confirmed; check for relationship with outlier status and genotype
chisq.test((X$GluR4_N > 0.2), y)$p.value

# Relationship isn't significant; remove outliers
remove <- which(X$GluR4_N > 0.2)
X <- X[-remove, ]
y <- y[-remove]

# Test/Train split
set.seed(314)
training <- createDataPartition(y, p = 0.7, list = FALSE)

X_train <- X[training, ]
X_test <- X[-training, ]

y_train <- y[training]
y_test <- y[-training]

#-------------------------------------------------------------------------------

# Dimensionality Reduction
#-------------------------------------------------------------------------------
# Center, scale, and impute
prep <- preProcess(X_train, method = c('center', 'scale', 'medianImpute'))

X_train_pre <- predict(prep, X_train)
X_test_pre  <- predict(prep, X_test)

# Perform PCA; 28 principal components required to capture 95% of variance
set.seed(314)
pca <- prcomp(X_train_pre)
summary(pca)

# Final PCA
pca_final <- prcomp(X_train_pre, rank = 28)
summary(pca_final)

# Fit PCA to data
X_train_pca <- as.data.frame(predict(pca_final, X_train_pre))
X_test_pca  <- as.data.frame(predict(pca_final, X_test_pre))

#-------------------------------------------------------------------------------

# Model Type Selection
#-------------------------------------------------------------------------------
# Create custom model training function
my_train <- function(x, y, method) {
  
  if (method == 'logreg') {
    mod <- train(
      x = x,
      y = y,
      method = 'glm',
      family = 'binomial',
      tuneLength = 5,
      trControl = trainControl(method = 'cv', number = 10, verboseIter = TRUE)
    )
  } else {
    mod <- train(
      x = x,
      y = y,
      method = method,
      tuneLength = 5,
      trControl = trainControl(method = 'cv', number = 10, verboseIter = TRUE)
    )
  }
  
  return(mod)
  
}

# Train models
set.seed(314)
fit_logreg <- my_train(X_train_pca, y_train, 'logreg')

set.seed(314)
fit_bayes <- my_train(X_train_pca, y_train, 'nb')

set.seed(314)
fit_rf <- my_train(X_train_pca, y_train, 'ranger')

set.seed(314)
fit_gbt <- my_train(X_train_pca, y_train, 'xgbTree')

set.seed(314)
fit_svm <- my_train(X_train_pca, y_train, 'svmLinear')

# Assess Performance; random forest performs best, with gradient boosted trees
# at a close second
confusionMatrix(predict(fit_logreg, X_test_pca), y_test)
confusionMatrix(predict(fit_bayes, X_test_pca), y_test)
confusionMatrix(predict(fit_rf, X_test_pca), y_test)
confusionMatrix(predict(fit_gbt, X_test_pca), y_test)
confusionMatrix(predict(fit_svm, X_test_pca), y_test)

#-------------------------------------------------------------------------------

# Parameter Tuning
#-------------------------------------------------------------------------------
# Create training grid; original train model used mtry of 2
my_grid <- expand.grid(
  mtry = seq(2, 8, by = 2),
  splitrule = c('gini', 'extratrees'),
  min.node.size = c(3, 5, 7, 9)
)

# Train model using custom grid
set.seed(314)
final_rf <- train(
  x = X_train_pca, 
  y = y_train,
  method = 'ranger',
  tuneGrid = my_grid,
  trControl = trainControl(method = 'oob', verboseIter = TRUE)
)

#-------------------------------------------------------------------------------

# Results
#-------------------------------------------------------------------------------
confusionMatrix(predict(final_rf, X_test_pca), y_test)


