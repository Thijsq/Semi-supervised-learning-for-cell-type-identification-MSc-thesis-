####################
# Author: Thijs Quast
# Contact: thijsquast@gmail.com
####################

# Loading training and test data with common genes
# Both datasets overlap in 12333 genes, the last column is the cell type
# Cells are in the rows, genes are in the columns

library(dplyr)
shekhar <- read.csv("shekhar_normal_downsampled.csv")

start_time <- Sys.time()
macosko <- read.csv("macosko_normal.csv")
end_time <- Sys.time()
running_time <- end_time - start_time
running_time

library(glmnet)

# Cross validation grouped LASSO regression to find optimal lambda paramter, also measure running time #

start_time <- Sys.time()
model_cv <- cv.glmnet(as.matrix(shekhar[,-ncol(shekhar)]), as.matrix(shekhar[, ncol(shekhar)]), 
                     alpha = 1, family = "multinomial", lambda = seq(0, 1, 0.01), type.multinomial = "grouped")
end_time <- Sys.time()
running_time <- end_time - start_time
running_time

# Results: CV for grouped lasso, optimal lambda = 0.01


# Preprocess data for grouped LASSO regression
covariates <- shekhar[,-ncol(shekhar)]
response <- shekhar[,ncol(shekhar)]

# Train grouped LASSO regression on original data representation
start_time <- Sys.time()
model_lasso <- glmnet(as.matrix(covariates), as.matrix(response), 
                      alpha = 1, family = "multinomial", lambda = 0.01, type.multinomial = "grouped")

# Predict using trained model
# type = "response" in order to obtain probabilities
predictions <- predict(model_lasso, newx = as.matrix(macosko[,-ncol(macosko)]), type = "response")
predictions <- as.matrix(predictions[,,1])
predictions


end_time <- Sys.time()
running_time <- end_time - start_time
running_time


# Obtain maximum probability per prediction
# This is used for the self learning algorithm
predictions_max <- c()
for (i in 1:nrow(macosko)){
  predictions_max[i] <- max(predictions[i,])
}

# Plot histogram of maximum probability per predictions
hist(predictions_max, breaks = 300)



# Self labelling implementation 

train <- shekhar
test <- macosko

# Just make a copy of test data, need lateron
test2 <- test


# Supervised_model
# This is the grouped LASSO regression in supervised learning format, so no self-learning
supervised_model <- glmnet(as.matrix(train[,-ncol(train)]), as.matrix(train[,ncol(train)]), alpha = 1, 
                           family = "multinomial", lambda = 0.01, type.multinomial = "grouped")

# Predict classes, obtain confusion matrix and accuracy score
supervised_predictions <- predict(supervised_model, newx = as.matrix(test2[,-ncol(test2)]), type = "class")
supervised_confusion <- table(test2$cell_type, supervised_predictions)
accuracy_supervised <- sum(diag(supervised_confusion))/sum(supervised_confusion)
supervised_confusion

# Store all predictions in a dataframe, also store all class predictions as a number, this is lateron used
# When computing classification reports (In Python)
supervised_predictions <- as.data.frame(supervised_predictions)
colnames(supervised_predictions) <- c("classes")
supervised_predictions$nr <- NA
supervised_predictions$nr[supervised_predictions$classes == "amacrine"] = 1
supervised_predictions$nr[supervised_predictions$classes == "bipolar"] = 2
supervised_predictions$nr[supervised_predictions$classes == "cones"] = 3
supervised_predictions$nr[supervised_predictions$classes == "muller"] = 4
supervised_predictions$nr[supervised_predictions$classes == "rods"] = 5

# Store predictions file
#write.csv(supervised_predictions, "lasso_original_mouse.csv", row.names = FALSE)



# Here is where the self-learning algorithm starts
# 1:10 is iterations 1 to 10, value of 10 can be modified
for (iteration in 1:10){
  
  print(paste("Model is currently at iteration:", iteration))
  
  rownames(test) <- c(1:nrow(test))
  rownames(train) <- c(1:nrow(train))
  
  # Train model
  model <- glmnet(x = as.matrix(train[,-ncol(train)]), y = as.matrix(train[,ncol(train)]),
                  family = "multinomial", alpha = 1, lambda = 0.01, type.multinomial = "grouped")
  
  # Predict probabilities
  predictions <- predict(model, as.matrix(test[,-ncol(test)]), type="response")
  predictions <- predictions[,,1]
  #predictions[is.nan(predictions)] <- 1
  
  # Store predictions in dataframe
  predictions <- as.data.frame(predictions)
  
  # Obtain maximum probability per prediction, needed for self-learning
  predictions_max <- c()
  
  for (i in 1:nrow(test)){
    predictions_max[i] <- max(predictions[i,])
  }
  hist(predictions_max, breaks = sqrt(length(predictions_max)))
  
  # Obtain predicted class per prediction
  predicted_classes <- apply(predictions, 1, which.max)
  
  # Store predicted classes as class names
  predicted_classes[predicted_classes == 1] = "amacrine"
  predicted_classes[predicted_classes == 2] = "bipolar"
  predicted_classes[predicted_classes == 3] = "cones"
  predicted_classes[predicted_classes == 4] = "muller"
  predicted_classes[predicted_classes == 5] = "rods"
  
  # Obtain indices of confident predictions,
  # those which exceed the confidence threshold of 95%
  confident_indices <- which(predictions_max > 0.95)
  
  # And if statement, if there are no more confident predictions, stop the algorithm
  if (length(confident_indices) == 0){
    print(paste("Model has converged after number of iterations:", iteration-1))
    break
  } else {
    
    # Store predicted classes of confident predictions
    confident_classes <- predicted_classes[confident_indices]
    
    # Store confident predictions
    predictions_to_train <- test[confident_indices,]
    predictions_to_train$cell_type <- predicted_classes[confident_indices]
    
    # Remove confident predictions from the 'test' data
    test <- test[-confident_indices,]
    # Append confident predictions to the training data. This is self-learning
    train <- as.data.frame(rbind(train, predictions_to_train))
  }
}

# Semi supervised model
# Store model that comes out of self-learning algorithm in final_model
final_model <- model

# Predict on test2, which is a copy of the original test data
final_predictions <- predict(final_model, as.matrix(test2[,-ncol(test2)]), type="class")

# Obtain confusion matrix
semi_supervised_confusion <- table(test2$cell_type, final_predictions)

# Store predictions made by self-learning algorithm in a dataframe
final_predictions <- as.data.frame(final_predictions)
colnames(final_predictions) <- c("classes")
final_predictions$nr <- NA
final_predictions$nr[final_predictions$classes == "amacrine"] = 1
final_predictions$nr[final_predictions$classes == "bipolar"] = 2
final_predictions$nr[final_predictions$classes == "cones"] = 3
final_predictions$nr[final_predictions$classes == "muller"] = 4
final_predictions$nr[final_predictions$classes == "rods"] = 5

# Save as csv file
#write.csv(final_predictions, "ssl_lasso_original_mouse.csv", row.names = FALSE)



# Accuracy scores
accuracy_semi_supervised <- sum(diag(semi_supervised_confusion))/sum(semi_supervised_confusion)


# Results
supervised_confusion
accuracy_supervised

semi_supervised_confusion
accuracy_semi_supervised


### descriptive statistics ###
# genes
count_zeros <- function(x){sum(x==0)}
zeros <- apply(shekhar[, -ncol(shekhar)], 2, count_zeros)
non_zeros <- nrow(shekhar) - zeros

# cells
zeros <- apply(shekhar[, -ncol(shekhar)], 1, count_zeros)
non_zeros <- (ncol(shekhar)-1) - zeros
