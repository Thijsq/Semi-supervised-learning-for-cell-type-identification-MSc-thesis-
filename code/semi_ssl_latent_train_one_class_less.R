####################
# Author: Thijs Quast
# Contact: thijsquast@gmail.com
####################

#
# This is the scenario where one class is missing in the training data
#


# For this file, manually remove one class from the test data, 
# merge both datasets and transition them into latent space.
# This is then the input file for latent_mouse

# Code comments for the self-learning implementation can be found in other self-learning 
# code files shared in the folder

latent_mouse <- read.csv("integrated_latent_shekhar4celltypes.csv")
latent_mouse <- latent_mouse[,-1]

shekhar_labels <- read.csv("shekhar_4labels.csv")
shekhar_labels <- shekhar_labels[,2]
macosko_labels <- read.csv("macosko_labels.csv")
macosko_labels <- macosko_labels[,2]



# Self labelling implementation 2

train <- latent_mouse[1:length(shekhar_labels),]
train$labels <- shekhar_labels



test <- latent_mouse[(length(shekhar_labels) + 1):nrow(latent_mouse),]

# check if dimensions match
nrow(train) + nrow(test) == nrow(latent_mouse)

test2 <- test
test2$labels <- macosko_labels



##### Lasso regression ########
library(glmnet)
#model_cv <- cv.glmnet(as.matrix(train[,-31]), as.matrix(train[,31]), alpha = 1, family = "multinomial", 
#                      lambda = seq(0, 1, 0.01), type.multinomial = "grouped")

# optimal lambda = 0

supervised_model <- glmnet(as.matrix(train[,-31]), as.matrix(train[,31]), alpha = 1, 
                           family = "multinomial", lambda = 0, type.multinomial = "grouped")
predictions <- predict(supervised_model, as.matrix(test2[,-31]), type = "response")
predictions <- predictions[,,1]

predicted_classes <- rep(NA, nrow(predictions))

predictions_max <- c()

for (i in 1:nrow(test)){
  predictions_max[i] <- max(predictions[i,])
}
hist(predictions_max, breaks = sqrt(nrow(test)))


predicted_classes <- apply(predictions, 1, which.max)
predicted_classes[predicted_classes == 1] = "bipolar"
predicted_classes[predicted_classes == 2] = "cones"
predicted_classes[predicted_classes == 3] = "muller"
predicted_classes[predicted_classes == 4] = "rods"

# When the predictions are less than 70% sure, I consider it an unrecognized cell type
unsure_predictions <- which(predictions_max < 0.70)
predicted_classes[unsure_predictions] = "unrecognized celltype"


confusion1 <- table(macosko_labels, factor(predicted_classes, 
                                           levels = c("unrecognized celltype",
                                                      "bipolar", "cones",
                                                      "muller", "rods")))
accuracy_supervised <- sum(diag(confusion1))/sum(confusion1)

supervised_predictions <- as.data.frame(predicted_classes)
colnames(supervised_predictions) <- c("classes")
supervised_predictions$nr <- NA
supervised_predictions$nr[supervised_predictions$classes == "unrecognized celltype"] = 1
supervised_predictions$nr[supervised_predictions$classes == "bipolar"] = 2
supervised_predictions$nr[supervised_predictions$classes == "cones"] = 3
supervised_predictions$nr[supervised_predictions$classes == "muller"] = 4
supervised_predictions$nr[supervised_predictions$classes == "rods"] = 5

write.csv(supervised_predictions, "lasso_latent_shekhar4.csv", row.names = FALSE)


for (iteration in 1:10){
  
  rownames(test) <- c(1:nrow(test))
  rownames(train) <- c(1:nrow(train))
  
  model <- glmnet(as.matrix(train[,-31]), as.matrix(train[,31]), alpha = 1, family = "multinomial", 
                  lambda = 0, type.multinomial = "grouped")
  
  predictions <- predict(model, as.matrix(test), type="response")
  predictions <- predictions[,,1]
  predictions <- as.data.frame(predictions)
  predictions
  
  predictions_max <- c()
  
  for (i in 1:nrow(test)){
    predictions_max[i] <- max(predictions[i,])
  }
  hist(predictions_max, breaks = sqrt(nrow(test)))
  
  predicted_classes <- apply(predictions, 1, which.max)
  predicted_classes[predicted_classes == 1] = "bipolar"
  predicted_classes[predicted_classes == 2] = "cones"
  predicted_classes[predicted_classes == 3] = "muller"
  predicted_classes[predicted_classes == 4] = "rods"
  
  confident_indices <- which(predictions_max > 0.99)
  
  print(paste("Model is currently at iteration:", iteration))
  
  if (length(confident_indices) == 0){
    print(paste("Model has converged after number of iterations:", iteration-1))
    break
  } else {
    
    
    confident_classes <- predicted_classes[confident_indices]
    
    predictions_to_train <- test[confident_indices,]
    predictions_to_train$labels <- predicted_classes[confident_indices]
    
    test <- test[-confident_indices,]
    train <- as.data.frame(rbind(train, predictions_to_train))
  }
}

final_model <- model
predictions <- predict(final_model, as.matrix(test2[,-31]), type="response")
predictions <- predictions[,,1]

predicted_classes <- rep(NA, nrow(predictions))

predictions_max <- c()
for (i in 1:nrow(test)){
  predictions_max[i] <- max(predictions[i,])
}
hist(predictions_max)

predicted_classes <- apply(predictions, 1, which.max)
predicted_classes[predicted_classes == 1] = "bipolar"
predicted_classes[predicted_classes == 2] = "cones"
predicted_classes[predicted_classes == 3] = "muller"
predicted_classes[predicted_classes == 4] = "rods"

unsure_predictions <- which(predictions_max < 0.70)
predicted_classes[unsure_predictions] = "unrecognized celltype"

confusion2 <- table(macosko_labels, factor(predicted_classes, 
                                           levels = c("unrecognized celltype",
                                                      "bipolar", "cones",
                                                      "muller", "rods")))


final_predictions <- as.data.frame(predicted_classes)
colnames(final_predictions) <- c("classes")
final_predictions$nr <- NA
final_predictions$nr[final_predictions$classes == "unrecognized celltype"] = 1
final_predictions$nr[final_predictions$classes == "bipolar"] = 2
final_predictions$nr[final_predictions$classes == "cones"] = 3
final_predictions$nr[final_predictions$classes == "muller"] = 4
final_predictions$nr[final_predictions$classes == "rods"] = 5

write.csv(final_predictions, "ssl_lasso_latent_shekhar4.csv", row.names = FALSE)
confusion1
confusion2


