library(glmnet)

# important -- this makes sure our runs are consistent and reproducible
set.seed(0)

file <- "TCGA_breast_cancer_ERstatus_allGenes.txt"
nfold <- 5
sd_threashold <- 1
  
header <- scan(file, nlines = 1, sep="\t", what = character())
data <- read.table(file, skip = 2, header = FALSE, sep = "\t", quote = "", check.names=FALSE)

header[1] <- "gene_id"
names(data) <- header

header2 <- scan(file, skip = 1, nlines = 1, sep="\t", what = character())

# cleanup - remove genes with sd < 1
# compute sd for each gene
data_sd<-sapply(seq(nrow(data)), function(x) { as.numeric(sd(data[x,-1])) })

# add gene names to the sd list
data_sd_names<-cbind(data.frame(data_sd),data[,1])

# create an "include" list of all those genes where sd > threshold
include_list <- data_sd_names[data_sd_names[,1]>sd_threashold,2]

Positive <- data[data$gene_id %in% include_list,header2=='Positive']
Negative <- data[data$gene_id %in% include_list,header2=='Negative']

# define function cross_valid so we can rerun the cross validation with various parameters
cross_validation <- function (nfold, alg="centroid") {

  Positive_groups <- split(sample(colnames(Positive)), 1+(seq_along(colnames(Positive)) %% nfold))
  Negative_groups <- split(sample(colnames(Negative)), 1+(seq_along(colnames(Negative)) %% nfold))
  
  result <- array()
  
  for (test_group in 1:nfold) {
    
    testPos <- Positive[,colnames(Positive) %in% unlist(Positive_groups[test_group])]
    testNeg <- Negative[,colnames(Negative) %in% unlist(Negative_groups[test_group])]
    
    trainingPos <- Positive[,!(colnames(Positive) %in% unlist(Positive_groups[test_group]))]
    trainingNeg <- Negative[,!(colnames(Negative) %in% unlist(Negative_groups[test_group]))]
    
    # Feature selection -- 
    
    # compute t-statistic for each row
    training_t_stat<-data.frame(sapply(seq(nrow(trainingPos)), function(x) { abs(as.numeric(t.test(trainingPos[x,], trainingNeg[x,])$statistic)) }))
    
    # add gene id column
    training_t_stat_geneid<-cbind(training_t_stat,rownames(trainingPos))
    colnames(training_t_stat_geneid) <- c('t','id')
    
    # pick top 50 based on t-statistic
    selected_genes <- head(training_t_stat_geneid[order(-training_t_stat_geneid$t),],n=top_num)[,2]
    
    # narrow down the list of genes based on t-statistic
    testPos <- testPos[rownames(testPos) %in% selected_genes,]
    testNeg <- testNeg[rownames(testNeg) %in% selected_genes,]
    trainingPos <- trainingPos[rownames(trainingPos) %in% selected_genes,]
    trainingNeg <- trainingNeg[rownames(trainingNeg) %in% selected_genes,]
    
    if (alg == "centroid") {
      centroidPos <- rowMeans(trainingPos)
      centroidNeg <- rowMeans(trainingNeg)
      
      misclassifiedPos <- sum(sapply(testPos, function(x) { sqrt(sum((x-centroidPos)^2))-sqrt(sum((x-centroidNeg)^2))>0 }))
      misclassifiedNeg <- sum(sapply(testNeg, function(x) { sqrt(sum((x-centroidPos)^2))-sqrt(sum((x-centroidNeg)^2))<0 }))
    }
    
    if (alg == "GLM") {
      # format training data
      training_data <- rbind(cbind(data.frame(t(trainingPos)), cancer=0), cbind(data.frame(t(trainingNeg)), cancer=1))
      
      # format test data
      test_data_pos <- data.frame(t(testPos))
      test_data_neg <- data.frame(t(testNeg))
      
      # fit model
      model <- glm(cancer~., data=training_data, family=binomial)
      
      # make predictions
      p_positive <- predict(model, newdata=test_data_pos, type="response")
      p_negative <- predict(model, newdata=test_data_neg, type="response")
      
      # calculate misclassified
      misclassifiedPos <- sum(ifelse(p_positive < 0.5, 0, 1))
      misclassifiedNeg <- sum(ifelse(p_negative >= 0.5, 0, 1))
    }
    
    if (alg == "Lasso") {
      # format training data 
      # glmnet requires a matrix 'x' (genes) and a vector 'y' (cancer labels)
      x_train <- rbind(t(trainingPos), t(trainingNeg))
      y_train <- c(rep(0, ncol(trainingPos)), rep(1, ncol(trainingNeg)))
      
      # format test data (must be converted to matrix for glmnet)
      test_data_pos <- as.matrix(t(testPos))
      test_data_neg <- as.matrix(t(testNeg))
      
      # fit model
      # cv.glmnet automatically solves the convergence/error issues
      # alpha=1 means Lasso
      model <- cv.glmnet(as.matrix(x_train), y_train, family="binomial", alpha=1)
      
      # make predictions
      # use s="lambda.min" to use the best penalty found
      p_positive <- predict(model, newx=test_data_pos, s="lambda.min", type="response")
      p_negative <- predict(model, newx=test_data_neg, s="lambda.min", type="response")
      
      # calculate misclassified
      misclassifiedPos <- sum(ifelse(p_positive < 0.5, 0, 1))
      misclassifiedNeg <- sum(ifelse(p_negative >= 0.5, 0, 1))
    }
    
    result[test_group] <- (misclassifiedPos+misclassifiedNeg)/(ncol(testPos)+ncol(testNeg))
  }
 
  paste0(round(mean(result),4)," sd=(",round(sd(result),4),")")
}

# Update top_num to 50
top_num <- 50
results_50_centroid <- cross_validation(5, alg="centroid")
results_50_glm <- cross_validation(5, alg="GLM")
results_50_lasso <- cross_validation(5, alg="Lasso")

# Update top_num to 100
top_num <- 100
results_100_centroid <- cross_validation(5, alg="centroid")
results_100_glm <- cross_validation(5, alg="GLM")
results_100_lasso <- cross_validation(5, alg="Lasso")

# Create a final comparison table
fs_summary <- data.frame(
  "Genes" = c("Top 50", "Top 50", "Top 50", "Top 100", "Top 100", "Top 100"),
  "Alg"   = c("Centroid", "GLM", "Lasso", "Centroid", "GLM", "Lasso"),
  "Result" = c(results_50_centroid, results_50_glm, results_50_lasso, results_100_centroid, results_100_glm, results_100_lasso)
)

kable(fs_summary)

