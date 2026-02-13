# Project 1 template

# important -- this makes sure our runs are consistent and reproducible
set.seed(0)

file <- "TCGA_breast_cancer_ERpositive_vs_ERnegative_PAM50.tsv"
first10 <- c('NAT1','BIRC5','BAG1','BCL2','BLVRA','CCNB1','CCNE1','CDC6','CDC20','CDH3')
nfold=5
  
header <- scan(file, nlines = 1, sep="\t", what = character())
data <- read.table(file, skip = 2, header = FALSE, sep = "\t", quote = "", check.names=FALSE)

header[1] <- "gene_id"
names(data) <- header

header2 <- scan(file, skip = 1, nlines = 1, sep="\t", what = character())

Positive <- data[data$gene_id %in% first10,header2=='Positive']
Negative <- data[data$gene_id %in% first10,header2=='Negative']


  # split each cancer type samples into nfold groups
  Positive_groups <- split(sample(colnames(Positive)), 1+(seq_along(colnames(Positive)) %% nfold))
  Negative_groups <- split(sample(colnames(Negative)), 1+(seq_along(colnames(Negative)) %% nfold))
  
  centroid_result <- array()
  glm_result <- array()
 
  # iterate from 1 to nfold groups -- to choose test group
  for (test_group in 1:nfold) {
   
    # return all samples in the chosen test group
    testPositive <- Positive[,colnames(Positive) %in% unlist(Positive_groups[test_group])]
    testNegative <- Negative[,colnames(Negative) %in% unlist(Negative_groups[test_group])]
    
    # return all samples *not* in the chosen test group 
    trainingPositive <- Positive[,!(colnames(Positive) %in% unlist(Positive_groups[test_group]))]
    trainingNegative <- Negative[,!(colnames(Negative) %in% unlist(Negative_groups[test_group]))]

    # compute centroid for each cancer type -- mean for each gene based on all samples
    # note -- rows are gene
    centroidPositive <- rowMeans(trainingPositive)
    centroidNegative <- rowMeans(trainingNegative)

    #kNNC
    # For each sample in the test set decide whether it will be classified
    # distance from centroid Positive: sum(abs(x-centroidPositive))
    # distance from centroid Negative: sum(abs(x-centroidNegative))
    # distance is a sum of distances over all genes
    # misclassification if when the distance is greater from centroid associated with known result
    centroid_misclassifiedPositive <- sum(sapply(testPositive, function(x) { sum(abs(x-centroidPositive))>sum(abs(x-centroidNegative)) }))
    centroid_misclassifiedNegative <- sum(sapply(testNegative, function(x) { sum(abs(x-centroidPositive))<sum(abs(x-centroidNegative)) }))

    # GLM
    # prepare/transform training data
    glm_training_data <- rbind(cbind(data.frame(t(trainingPositive)), cancer=0), cbind(data.frame(t(trainingNegative)), cancer=1))
    # transform test data
    glm_test_data_positive <- data.frame(t(testPositive))
    glm_test_data_negative <- data.frame(t(testNegative))
    # fit model
    glm_model <- glm(cancer~., data=glm_training_data, family=binomial)
    # make predictions
    glm_p_positive <- predict(glm_model, newdata=glm_test_data_positive, type="response")
    glm_p_negative <- predict(glm_model, newdata=glm_test_data_negative, type="response")
    # calculate misclassified
    glm_misclassifiedPositive <- sum(ifelse(glm_p_positive<0.5,0,1))
    glm_misclassifiedNegative <- sum(ifelse(glm_p_negative>=0.5,0,1))
    
    centroid_result[test_group] <- (centroid_misclassifiedPositive+centroid_misclassifiedNegative)/(ncol(testPositive)+ncol(testNegative))
    glm_result[test_group] <- (glm_misclassifiedPositive+glm_misclassifiedNegative)/(ncol(testPositive)+ncol(testNegative))
 }
 
 c(mean(centroid_result), sd(centroid_result))
 c(mean(glm_result), sd(glm_result))


#x<-data.frame(three=cross_validation(3), five=cross_validation(5), ten=cross_validation(10))
#rownames(x) <- c('mean','sd')
#x

