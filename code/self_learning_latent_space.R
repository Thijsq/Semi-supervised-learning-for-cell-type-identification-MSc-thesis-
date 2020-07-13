####################
# Author: Thijs Quast
# Contact: thijsquast@gmail.com
####################

# Import latent dimensionality created by scVI autoencoder
latent_mouse <- read.csv("latent_integrated_mouse_downsampled.csv")

# Drop first column
latent_mouse <- latent_mouse[,-1]

# Import cell type labels for both datasets
shekhar_labels <- read.csv("shekhar_labels_downsampled.csv")
shekhar_labels <- shekhar_labels[,2]
macosko_labels <- read.csv("macosko_labels.csv")
macosko_labels <- macosko_labels[,2]



# Self learning implementation 

# Now transform the integretad latent dimensionality back to training and test data
train <- latent_mouse[1:length(shekhar_labels),]
train$labels <- shekhar_labels



test <- latent_mouse[(length(shekhar_labels) + 1):nrow(latent_mouse),]

# check if dimensions match
nrow(train) + nrow(test) == nrow(latent_mouse)

# Create copy of test data
test2 <- test
test2$labels <- macosko_labels



##### Grouped Lasso regression ########
library(glmnet)
# Run cross validation to find optimal lambda value, uncomment model_cv below
#model_cv <- cv.glmnet(as.matrix(train[,-31]), as.matrix(train[,31]), alpha = 1, family = "multinomial", 
#                     lambda = seq(0, 1, 0.01), type.multinomial = "grouped")

# Run supervised model
supervised_model <- glmnet(as.matrix(train[,-31]), as.matrix(train[,31]), alpha = 1, family = "multinomial", lambda = 0, type.multinomial = "grouped")
predictions <- predict(supervised_model, as.matrix(test2[,-31]), type = "response")
predictions <- predictions[,,1]

# Obtain maximum probability per prediction
predictions_max <- c()

for (i in 1:nrow(test)){
  predictions_max[i] <- max(predictions[i,])
}
hist(predictions_max, breaks = sqrt(nrow(test)))



# Predict classes for supervised model
supervised_model <- glmnet(as.matrix(train[,-31]), as.matrix(train[,31]), alpha = 1, family = "multinomial", lambda = 0, type.multinomial = "grouped")
supervised_predictions <- predict(supervised_model, as.matrix(test2[,-31]), type = "class")

# Generate confusion matrix for grouped LASSO regression in latent dimensionality
supervised_confusion <- table(test2$labels, supervised_predictions)
supervised_confusion

# Store predictions, and numbered labels
supervised_predictions <- as.data.frame(supervised_predictions)
colnames(supervised_predictions) <- c("classes")
supervised_predictions$nr <- NA
supervised_predictions$nr[supervised_predictions$classes == "amacrine"] = 1
supervised_predictions$nr[supervised_predictions$classes == "bipolar"] = 2
supervised_predictions$nr[supervised_predictions$classes == "cones"] = 3
supervised_predictions$nr[supervised_predictions$classes == "muller"] = 4
supervised_predictions$nr[supervised_predictions$classes == "rods"] = 5

# Save as csv file
#write.csv(supervised_predictions, "lasso_latent_mouse.csv", row.names = FALSE)


for (iteration in 1:10){
  
  rownames(test) <- c(1:nrow(test))
  rownames(train) <- c(1:nrow(train))
  
  # train model
  model <- glmnet(as.matrix(train[,-31]), as.matrix(train[,31]), alpha = 1, family = "multinomial", 
                  lambda = 0, type.multinomial = "grouped")
  
  predictions <- predict(model, as.matrix(test), type="response")
  predictions <- predictions[,,1]
  predictions <- as.data.frame(predictions)
  predictions
  
  # Obtain maximum predictions
  predictions_max <- c()
  
  for (i in 1:nrow(test)){
    predictions_max[i] <- max(predictions[i,])
  }
  
  # Plot histogram
  hist(predictions_max, breaks = sqrt(nrow(test)))
  
  predicted_classes <- apply(predictions, 1, which.max)
  predicted_classes[predicted_classes == 1] = "amacrine"
  predicted_classes[predicted_classes == 2] = "bipolar"
  predicted_classes[predicted_classes == 3] = "cones"
  predicted_classes[predicted_classes == 4] = "muller"
  predicted_classes[predicted_classes == 5] = "rods"
  
  # Confident indices, above 99% threshold
  confident_indices <- which(predictions_max > 0.99)
  
  print(paste("Model is currently at iteration:", iteration))
  
  if (length(confident_indices) == 0){
    print(paste("Model has converged after number of iterations:", iteration-1))
    break
  } else {
    
    # Confident predictions
    confident_classes <- predicted_classes[confident_indices]
    
    
    predictions_to_train <- test[confident_indices,]
    predictions_to_train$labels <- predicted_classes[confident_indices]
    
    # Remove confident predictions from test data 
    # Append training data with confident predictions
    test <- test[-confident_indices,]
    train <- as.data.frame(rbind(train, predictions_to_train))
  }
}

# final model
final_model <- model
final_predictions <- predict(final_model, as.matrix(test2[,-31]), type="class")
semi_supervised_confusion <- table(test2$labels, final_predictions)


# Store final predictions
final_predictions <- as.data.frame(final_predictions)
colnames(final_predictions) <- c("classes")
final_predictions$nr <- NA
final_predictions$nr[final_predictions$classes == "amacrine"] = 1
final_predictions$nr[final_predictions$classes == "bipolar"] = 2
final_predictions$nr[final_predictions$classes == "cones"] = 3
final_predictions$nr[final_predictions$classes == "muller"] = 4
final_predictions$nr[final_predictions$classes == "rods"] = 5


#write.csv(final_predictions, "ssl_lasso_latent_mouse.csv", row.names = FALSE)



accuracy_supervised <- sum(diag(supervised_confusion))/sum(supervised_confusion)
accuracy_semi_supervised <- sum(diag(semi_supervised_confusion))/sum(semi_supervised_confusion)

accuracy_supervised
accuracy_semi_supervised

