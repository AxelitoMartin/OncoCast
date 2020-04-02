####################################################################
##################### ONCOCAST ANALYSIS  ###########################
####################################################################

#' OncoCast
#'
#' This functions let's the user select one or multiple statistical learning algorithms (penalized regression and generalized boosted regression).
#' This is intended for survival data, and the mehods can handle left truncated survival data with a high number of various of predictors.
#' The inputed data should be a data frame with columns representing the variables of interest while the rows should correspond to patients.
#' All methods selected will all return a list of length equal to the number of crossvalidation performed,
#' including the predicted risk score at each cross validation for all the patients falling in the test set.
#' @param data Data frame with variables as columns and patients as rows. Must have no missing data and should contain only the outcome and the predictors to be used.
#' We recommend the time variables to use the month unit. Moreover note that these methods cannot handle categorical variables. Dummy binary variables must be created
#' prior to input.
#' @param family Name of the distribution of the outcome. Options are "cox", "binomial". Default is "cox".
#' @param formula A formula with the names of the variables to be used in the data frame provided in the first argument.
#'  eg : Surv(time,status)~. for cox; y ~ . for binomial (Note all the variable available will be used regardless of the right
#'  side of the formula).
#' @param method Character vector of the names of the method(s) to be used, options are : 1) LASSO ("LASSO") 2) Ridge ("RIDGE")
#'  3) Elastic Net ("ENET"), 4) Random Forest ("RF"), 5) Support vector machine ("SVM") 6) Boosted Forest ("GBM") 7) Neural network ("NN"). Note that GBM requires training in order to perform optimally. Arguments in that objectives are listed below.
#'  Default is ENET.
#' @param runs Number of cross validation iterations to be performed. Default is 100. We highly recommend performing at least 50.
#'
# #' @param sampling The method use for sampling, options are bootstrapping ("boot") and cross-validation ("cv").
# #' Default is cross-validation.
#'
#' @param cores If you wish to run this function in parallel set the number of cores to be used to be greater than 1. Default is 1.
#' CAUTION : Overloading your computer can lead to serious issues, please check how many cores are available on your machine
#' before selecting an option!
#' @param pathResults String of where the users wishes to output the results. Note that paths are relative in this context.
#' Default is current directory.
#' @param studyType String that will be the prefix to the name of the outputed results. Default is empty.
#' @param save Boolean value : Default is TRUE, the results will be saved with the specified name in the specified path. If FALSE the results
#' will be returned directly from the function and won't be saved. Be sure to save them in an object in your environment.
#' @param nonPenCol Name of variables you do not with to penalize (available only for LASSO, RIDGE and ENET methods). Default is NULL.
#' @param nTree Argument required for RF and GBM. Designates the number of trees to be built in each forest.
#' Default is 500.
#' @param interactions For GBM only. Integer specifying the maximum depth of each tree (i.e., the highest level of
#' variable interactions allowed). A value of 1 implies an additive model, a value
#' of 2 implies a model with up to 2-way interactions, etc. Default is 1 and 2 (as a vector c(1,2)).
#' @param shrinkage For GBM only. a shrinkage parameter applied to each tree in the expansion. Also known as
#' the learning rate or step-size reduction; 0.001 to 0.1 usually work, but a smaller
#' learning rate typically requires more trees. Default are c(0.001,0.01).
#' @param min.node For GBM only. Integer specifying the minimum number of observations in the terminal nodes
#' of the trees. Note that this is the actual number of observations, not the total
#' weight. Beware that nodes must be small if n is small. Default are c(10,20).
#' @param rf_gbm.save For RF, GBM, SVM and NN only. In order to perform proper validation we must save the entire fit. This will require more memory space to save.
#' If you plan to perform validation set to TRUE. Default is FALSE.
#' @param cv.folds Number of internal cross-validations to be performed in GBM.
#' @param mtry Number of features to include in each tree (for random forest). Default is 1/3 of all features.
#' @param replace Boolean to sample with or without replacement for the samples (for random forest). Default is T.
#' @param sample.fraction Fraction of samples to be used in random forest (default is 1).
#' @param max.depth Maximum of depth of trees to be grown. Default trees are grown as much as possible.
#' @param rf.node The minimal size of terminal nodes for random forest. Default is 5 (recommended for regression trees).
#' @param out.ties phcpe argument to calculate the concordance index. If out.ties is set to FALSE,
#' pairs of observations tied on covariates will be used to calculate the CPE. Otherwise, they will not be used.
#' @param epsilon.svm epsilon sequence for SVM. Default is seq(0,0.2,0.02)
#' @param cost.svm cost sequence for SVM. Default is 2^(2:9)
#' @param layers a vector input for neural network method, each entry representing the number of nodes to be used at the given hidden layer.
#' The number of hidden layers is determined by the length of the vector. Default is NULL in which case a single hidden layer will be used with diverse neuron numbers (recommended).
#' @param norm.nn Boolean value to decide if the data for the nerual network should be normalized. Default is true (recommended).
#' @return CI : For each iteration the concordance index of the generated model will be calculated on the testing set
#' @return fit : For LASSO, RIDGE and ENET methods return the refitted cox proportional hazard model with the beta coefficients found
#' in the penalized regression model.
#' @return The formula used to perform the analysis.
#' @return predicted : The predicted risk for all the patients in the testing set.
#' @return means : The mean value of the predictors that were not shrunken to zero in the penalized regression method.
#' @return method : The name of the method that was used to generate the output.
#' @return data : The data used to fit the model (available only in the first element of the list).
#' @keywords Selection, penalized regression
#' @export
#' @examples library(OncoCast)
#' test <- OncoCast(data=survData,formula = Surv(time,status)~.,
#' method = "LASSO",runs = 50,
#' save = F,nonPenCol = NULL,cores =2)

OncoCast <- function(data,family = "cox",formula, method = c("ENET"),
                     runs = 100,cores = 1,
                     pathResults = ".",studyType = "",save = T,
                     nonPenCol = NULL,
                     nTree=500,interactions=c(1,2),
                     shrinkage=c(0.001,0.01),min.node=c(10,20),rf_gbm.save = F,
                     out.ties=F,cv.folds=5,rf.node=5,mtry = floor(ncol(data)/3),replace = T,
                     sample.fraction=1,max.depth = NULL,
                     epsilon.svm = seq(0,0.2,0.02), cost.svm = 2^(2:9),
                     layers = NULL,norm.nn = T){

  # Missingness
  if(anyNA(data)){
    stop("ERROR : Missing data is not allowed at this time, please remove or impute missing data.")
  }

  # check arguments
  if(!(family %in% c("cox","binomial","gaussian"))) stop("ERROR: Family must be one of: 'cox', 'binomial', 'gaussian'.")
  if( anyNA(match(method,c("LASSO","ENET","RIDGE","RF","GBM","SVM","NN"))) ){stop("ERROR : The method you have selected is not available.")}


  if(runs < 0){stop("The number of cross-validation/bootstraps MUST be positive.")}
  else if(runs < 50){warning("We do not recommend using a number of cross-validation/bootstraps lower than 50.")}

  # perform data cleaning #
  dums <- apply(data,2,function(x){anyNA(as.numeric(as.character(x)))})
  data.save = F
  if(sum(dums) > 0){
    tmp <- data %>%
      select(which(dums)) %>%
      fastDummies::dummy_cols(remove_first_dummy = T) %>%
      select(-one_of(names(which(dums))))
    data <- as.data.frame(cbind(
      data %>% select(-one_of(names(which(dums)))),
      tmp
    ) %>% mutate_all(as.character) %>%
      mutate_all(as.numeric)
    )
    data.save = T
    warning("Character variables were transformed to dummy numeric variables. If you didn't have any character variables make sure all columns in your input data are numeric. The transformed data will be saved as part of the output.")
  }

  #
  print("Data check performed, ready for analysis.")

  ### Prepare cores
  cl <- makeCluster(cores,outfile = paste0(studyType,"_OncoCast_log.txt"))
  registerDoParallel(cl)

  ##### generate empty output objects #####
  LASSO <- NULL
  ENET <- NULL
  RIDGE <- NULL
  RF <- NULL
  GBM <- NULL
  SVM <- NULL
  NN <- NULL

  if(family == "cox"){
    # appropriate formula
    survFormula <- as.formula(formula)
    survResponse <- survFormula[[2]]

    if(!length(as.list(survResponse)) %in% c(3,4)){
      stop("ERROR : Response must be a 'survival' object with 'Surv(time, status)' or 'Surv(time1, time2, status)'.")
    }
    ### reprocess data
    if(length(as.list(survResponse)) == 3){
      colnames(data)[match(as.list(survResponse)[2:3],colnames(data))] <- c("time","status")
      LT = FALSE
      timevars <- match(c("time","status"),colnames(data))
      data <- data[,c(timevars,
                      c(1:ncol(data))[-timevars])]
    }
    if(length(as.list(survResponse)) == 4){
      colnames(data)[match(as.list(survResponse)[2:4],colnames(data))] <- c("time1","time2","status")
      LT = TRUE
      timevars <- match(c("time1","time2","status"),colnames(data))
      data <- data[,c(timevars,
                      c(1:ncol(data))[-timevars])]
    }


    ### max number of iterations in penalized ###
    iter.max = 10^4


    ########## LASSO #############
    final.lasso <- list()
    if("LASSO" %in% method) {
      print("LASSO SELECTED")
      LASSO <- foreach(run=1:runs) %dopar% {

        ### BUILD TRAINING AND TESTING SET ###
        set.seed(run)
        cat(paste("Run : ", run,"\n",sep=""), file=paste0(studyType,"LASSO_log.txt"), append=TRUE)
        # split data
        rm.samples <- sample(1:nrow(data), ceiling(nrow(data)*1/3),replace = FALSE)
        train <- data[-rm.samples,]
        test <- data[rm.samples,]

        # make new survival objects
        if(LT) {
          trainSurv <- with(train,Surv(time1,time2,status))
          testSurv <- with(test,Surv(time1,time2,status))
          if(is.null(nonPenCol)){
            opt <- try(optL1(trainSurv,data = train,penalized = train[,4:ncol(train)],fold = 5,trace=FALSE))}
          else{
            noPen.index <- match(nonPenCol,colnames(train))
            noPen.index.pen <-c(1:3,noPen.index)
            opt <- try(optL1(trainSurv,data = train,penalized = train[,-noPen.index.pen],
                             unpenalized = train[,nonPenCol],fold = 5,trace=FALSE))
          }
        }

        else {
          trainSurv <- with(train,Surv(time,status))
          testSurv <- with(test,Surv(time,status))
          if(is.null(nonPenCol)){
            opt <- try(optL1(trainSurv,data = train,penalized = as.matrix(train[,3:ncol(train)]),fold = 5,trace=FALSE))}
          else{
            noPen.index <- match(nonPenCol,colnames(train))
            noPen.index.pen <-c(1,2,noPen.index)
            opt <- try(optL1(trainSurv,data = train,penalized = train[,-noPen.index.pen],
                             unpenalized = train[,nonPenCol],fold = 5,trace=FALSE))
          }
        }


        if(typeof(opt) == "list" && length(coefficients(opt$fullfit)) != 0){
          optimal.coefs <- coefficients(opt$fullfit)
          coefs.left <- names(coefficients(opt$fullfit))
          if(LT) lasso.formula <- as.formula(paste("Surv(time1,time2,status) ~ ",paste(coefs.left, collapse= "+")))
          if(!LT) lasso.formula <- as.formula(paste("Surv(time,status) ~ ",paste(coefs.left, collapse= "+")))
          # Refit model with the corresponding coefficients
          lasso.fit <- coxph(lasso.formula, data=train,init=optimal.coefs,iter=0)
          # save what you need to save
          CI <- as.numeric(concordance(coxph(testSurv ~predict(lasso.fit, newdata=test),
                                             test))$concordance)
          CPE <- as.numeric(phcpe(coxph(testSurv ~predict(lasso.fit, newdata=test)),out.ties=out.ties ))
          predicted <- rep(NA,nrow(data))
          names(predicted) <- rownames(data)
          predicted[match(names(predict(lasso.fit, newdata=test)),names(predicted))] <- as.numeric(predict(lasso.fit, newdata=test))

          if(LT) {coefs <- rep(NA,ncol(data)-3); names(coefs) <- colnames(data)[4:ncol(data)]}
          if(!LT) {coefs <- rep(NA,ncol(data)-2); names(coefs) <- colnames(data)[3:ncol(data)]}
          coefs[match(coefs.left,names(coefs))] <- summary(lasso.fit)$coefficients[,1]

          final.lasso$formula <- formula
          final.lasso$CPE <- CPE
          final.lasso$CI <- CI
          final.lasso$fit <- coefs
          final.lasso$predicted <- predicted
          final.lasso$means <- lasso.fit$means
          final.lasso$method <- "LASSO"
          final.lasso$family <- family
        }
        else{
          final.lasso <- NULL
        }
        return(final.lasso)
      }

      if(save){
        save(LASSO,file = paste0(pathResults,"/",studyType,"_LASSO.Rdata"))}
    }

    ##### RIDGE #####

    final.ridge <- list()
    if("RIDGE" %in% method) {
      print("RIDGE SELECTED")
      RIDGE <- foreach(run=1:runs) %dopar% {

        ### BUILD TRAINING AND TESTING SET ###
        set.seed(run)
        cat(paste("Run : ", run,"\n",sep=""), file=paste0(studyType,"RIDGE_log.txt"), append=TRUE)
        # split data
        rm.samples <- sample(1:nrow(data), ceiling(nrow(data)*1/3),replace = FALSE)
        train <- data[-rm.samples,]
        test <- data[rm.samples,]


        # make new survival objects
        if(LT) {
          trainSurv <- with(train,Surv(time1,time2,status))
          testSurv <- with(test,Surv(time1,time2,status))
          if(is.null(nonPenCol)){
            opt <- try(optL1(trainSurv,data = train,penalized = train[,4:ncol(train)],fold = 5,trace=FALSE))}
          else{
            noPen.index <- match(nonPenCol,colnames(train))
            noPen.index.pen <-c(1:3,noPen.index)
            opt <- try(optL1(trainSurv,data = train,penalized = train[,-noPen.index.pen],
                             unpenalized = train[,nonPenCol],fold = 5,trace=FALSE))
          }
        }

        else {
          trainSurv <- with(train,Surv(time,status))
          testSurv <- with(test,Surv(time,status))
          if(is.null(nonPenCol)){
            opt <- try(optL1(trainSurv,data = train,penalized = train[,3:ncol(train)],fold = 5,trace=FALSE))}
          else{
            noPen.index <- match(nonPenCol,colnames(train))
            noPen.index.pen <-c(1,2,noPen.index)
            opt <- try(optL1(trainSurv,data = train,penalized = train[,-noPen.index.pen],
                             unpenalized = train[,nonPenCol],fold = 5,trace=FALSE))
          }
        }

        if(typeof(opt) == "list" && length(coefficients(opt$fullfit)) != 0){
          # get optimal coefficients
          optimal.coefs <- coefficients(opt$fullfit)
          coefs.left <- names(coefficients(opt$fullfit))
          if(LT) {ridge.formula <- as.formula(paste("Surv(time1,time2,status) ~ ",paste(coefs.left, collapse= "+")))}
          else {ridge.formula <- as.formula(paste("Surv(time,status) ~ ",paste(coefs.left, collapse= "+")))}

          # Refit model with the corresponding coefficients
          ridge.fit <- coxph(ridge.formula, data=train,init=optimal.coefs,iter=0)

          # save what you need to save
          CPE <- as.numeric(phcpe(coxph(testSurv ~predict(ridge.fit, newdata=test)),out.ties=out.ties))
          CI <- as.numeric(concordance(coxph(testSurv ~predict(ridge.fit, newdata=test),
                                             test))$concordance)
          predicted <- rep(NA,nrow(data))
          names(predicted) <- rownames(data)
          predicted[match(names(predict(ridge.fit, newdata=test)),names(predicted))] <- as.numeric(predict(ridge.fit, newdata=test))

          if(LT) {coefs <- rep(NA,ncol(data)-3); names(coefs) <- colnames(data)[4:ncol(data)]}
          if(!LT) {coefs <- rep(NA,ncol(data)-2); names(coefs) <- colnames(data)[3:ncol(data)]}
          coefs[match(coefs.left,names(coefs))] <- summary(ridge.fit)$coefficients[,1]

          final.ridge$CPE <- CPE
          final.ridge$CI <- CI
          final.ridge$fit <- coefs
          final.ridge$predicted <- predicted
          final.ridge$means <- ridge.fit$means
          final.ridge$method <- "RIDGE"
          final.ridge$formula <- formula
          final.ridge$family <- family
        }
        else{
          final.ridge <- NULL
        }

        return(final.ridge)
      }
      if(save){
        save(RIDGE,file = paste0(pathResults,"/",studyType,"_RIDGE.Rdata"))}
    }

    #### LASSO ENET #####
    final.enet <- list()
    if("ENET" %in% method) {
      print("ENET SELECTED")
      ENET <- foreach(run=1:runs) %dopar% {

        ### BUILD TRAINING AND TESTING SET ###
        set.seed(run)
        cat(paste("Run : ", run,"\n",sep=""), file=paste0(studyType,"ENET_log.txt"), append=TRUE)

        rm.samples <- sample(1:nrow(data), ceiling(nrow(data)*1/3),replace = FALSE)
        train <- data[-rm.samples,]
        test <- data[rm.samples,]

        opt <- NULL

        if(LT) {
          trainSurv <- with(train,Surv(time1,time2,status))
          testSurv <- with(test,Surv(time1,time2,status))
          noPen.index <- match(nonPenCol,colnames(train))
          noPen.index.pen <-c(1:3,noPen.index)}
        else{
          trainSurv <- with(train,Surv(time,status))
          testSurv <- with(test,Surv(time,status))
          noPen.index <- match(nonPenCol,colnames(train))
          noPen.index.pen <-c(1:2,noPen.index)}

        lambda1s <- c(0,0.01,0.05,0.1,0.2,0.5,1,2,5,10)
        lambda2s <- c(0,0.01,0.05,0.1,0.2,0.5,1,2,5,10)

        # lambda1s <- c(0.05,0.1)
        # lambda2s <- c(0.05,0.1)

        lambdaGrid <- expand.grid(.lambda1s = lambda1s,
                                  .lambda2s = lambda2s)

        allfits <- apply(lambdaGrid,1,function(x){

          if(is.null(nonPenCol)){
            opt <- try(cvl(trainSurv,data = train,
                           penalized = train[,-na.omit(noPen.index.pen)],
                           unpenalized = ~0,fold=5,lambda1 = as.numeric(x[1]),lambda2 = as.numeric(x[2]),
                           model = "cox",maxiter = iter.max),
                       silent=T)
          }
          else{
            opt <- try(cvl(trainSurv,data = train,
                           penalized = train[,-noPen.index.pen],
                           unpenalized = train[,nonPenCol],fold=5,lambda1 = as.numeric(x[1]),lambda2 = as.numeric(x[2]),
                           model = "cox",maxiter = iter.max),
                       silent=T)
          }


          if(typeof(opt) == "list" && length(coefficients(opt$fullfit)) != 0 && opt$fullfit@converged){
            optimal.coefs <- coefficients(opt$fullfit)
            coefs.left <- names(coefficients(opt$fullfit))
            if(LT) lasso.formula <- as.formula(paste("Surv(time1,time2,status) ~ ",paste(coefs.left, collapse= "+")))
            if(!LT) lasso.formula <- as.formula(paste("Surv(time,status) ~ ",paste(coefs.left, collapse= "+")))
            # Refit model with the corresponding coefficients
            lasso.fit <- coxph(lasso.formula, data=train,init=optimal.coefs,iter=0)
            # save what you need to save
            CPE <- as.numeric(phcpe(coxph(testSurv ~predict(lasso.fit, newdata=test)),out.ties=out.ties))
            CI <- as.numeric(concordance(coxph(testSurv ~predict(lasso.fit, newdata=test),
                                               test))$concordance)
            predicted <- rep(NA,nrow(data))
            names(predicted) <- rownames(data)
            predicted[match(names(predict(lasso.fit, newdata=test)),names(predicted))] <- as.numeric(predict(lasso.fit, newdata=test))

            if(LT) {coefs <- rep(NA,ncol(data)-3); names(coefs) <- colnames(data)[4:ncol(data)]}
            if(!LT) {coefs <- rep(NA,ncol(data)-2); names(coefs) <- colnames(data)[3:ncol(data)]}
            coefs[match(coefs.left,names(coefs))] <- summary(lasso.fit)$coefficients[,1]

            return(list("CPE"=CPE,"CI"=CI,"fit"=coefs,"predicted"=predicted,"means"=lasso.fit$means))
          }

          else{return(list("CPE"=NA,"CI"=NA,"fit"=NULL,"predicted"=NULL,"means"=NULL))}

        })

        if(sum(is.na(vapply(allfits,"[[","CPE",FUN.VALUE = numeric(1)))) == nrow(lambdaGrid)){
          final.enet <- NULL
          return(final.enet)
        }

        else{
          bestRun <- which.max(vapply(allfits,"[[","CPE",FUN.VALUE = numeric(1)))

          final.enet$CPE <- allfits[[bestRun]]$CPE
          final.enet$CI <- allfits[[bestRun]]$CI
          final.enet$fit <- allfits[[bestRun]]$fit
          final.enet$predicted <- allfits[[bestRun]]$predicted
          final.enet$means <- allfits[[bestRun]]$means
          final.enet$method <- "ENET"
          final.enet$formula <- formula
          final.enet$family <- family

          return(final.enet)
        }
      }
      if(save){
        save(ENET,file = paste0(pathResults,"/",studyType,"_ENET.Rdata"))
        ENET<- NULL}
    }



    ##########
    ### RF ###
    ##########
    if("RF" %in% method){
      final.rf <- list()
      print("RF SELECTED")
      RF <- foreach(run=1:runs) %dopar% {

        set.seed(run)
        cat(paste("Run : ", run,"\n",sep=""), file=paste0(studyType,"RF_log.txt"), append=TRUE)

        rm.samples <- sample(1:nrow(data), ceiling(nrow(data)*1/3),replace = FALSE)
        train <- data[-rm.samples,]
        test <- data[rm.samples,]

        if(LT){
          fit <- coxph(Surv(time1,time2,status)~1,data=train)
          train <- train[,-match(c("time1","time2","status"),colnames(train))]
          testSurv <- with(test,Surv(time1,time2,status))
        }
        else{
          fit <- coxph(Surv(time,status)~1,data=train)
          train <- train[,-match(c("time","status"),colnames(train))]
          testSurv <- with(test,Surv(time,status))
        }
        train$y <- residuals(fit,type="martingale")


        if(!is.null(max.depth)) rfGrid <- expand.grid(.mtry = mtry,
                                                      .n.trees = nTree,
                                                      .n.minobsinnode = rf.node,
                                                      .sample.fraction = sample.fraction,
                                                      .max.depth = max.depth)
        else rfGrid <- expand.grid(.mtry = mtry,
                                   .n.trees = nTree,
                                   .n.minobsinnode = rf.node,
                                   .sample.fraction = sample.fraction)

        BestPerf <- apply(rfGrid,1,function(x){
          set.seed(21071993)

          if(!is.null(max.depth)) rf <- ranger(formula = y~., data = train, num.trees = x[2],replace = replace,
                                               importance = "impurity",mtry = x[1],
                                               min.node.size = x[3],
                                               sample.fraction = x[4],
                                               max.depth = x[5])
          else rf <- ranger(formula = y~., data = train, num.trees = x[2],replace = replace,
                            importance = "impurity",mtry = x[1],
                            min.node.size = x[3],
                            sample.fraction = x[4])

          if(is.character(rf)){
            return(list("CPE"=0,"CI"=0,"infl"=NA,"predicted"=NA,
                        "bestTreeForPrediction"=NA))
          }

          predicted<- predict(rf,test)$predictions

          CPE <- as.numeric(phcpe(coxph(testSurv ~predicted,data=test),out.ties=out.ties))
          CI <- as.numeric(concordance(coxph(testSurv ~ predicted,
                                             test))$concordance)
          return(list("CPE"=CPE,"CI"=CI,"infl"=importance(rf),"predicted"=predicted,"rf"=rf))

        })


        if(sum(vapply(BestPerf,"[[","CPE",FUN.VALUE = numeric(1)),na.rm = T)==0){return(NULL)}
        index <- which.max(vapply(BestPerf,"[[","CPE",FUN.VALUE = numeric(1)))

        rf <- BestPerf[[index]]$rf
        final.rf$method <- "RF"
        final.rf$Vars <- importance(rf)

        predicted <- rep(NA,nrow(data))
        names(predicted) <- rownames(data)
        predicted[match(rownames(test),names(predicted))] <- predict(rf,test)$predictions
        final.rf$CPE <-  BestPerf[[index]]$CPE
        final.rf$CI <- BestPerf[[index]]$CI
        final.rf$predicted <- predicted
        final.rf$formula <- formula
        final.rf$family <- family

        if(rf_gbm.save){
          final.rf$RF <- rf
        }

        # if(tune.rf == T){
        #   train.task = makeRegrTask(data = train, target = "y")
        #   # Tuning
        #   res = tuneRanger(train.task, num.trees = nTree,  #measure = list(mse),
        #                    parameters = list(replace = FALSE, respect.unordered.factors = "order"),
        #                    num.threads = 1, iters = 50, save.file.path = NULL,show.info = F,
        #                    build.final.model = F)
        #
        #   if(!is.null(res$recommended.pars[1])) mtry <- as.numeric(res$recommended.pars[1])
        #   if(!is.null(res$recommended.pars[2])) rf.node <- as.numeric(res$recommended.pars[2])
        #   if(!is.null(res$recommended.pars[3])) sample.fraction <- as.numeric(res$recommended.pars[3])
        #   else sample.fraction <- 2/3
        #
        #   rf <- ranger(formula = y~., data = train, num.trees = nTree,replace = F,
        #                importance = "impurity",mtry = mtry,
        #                min.node.size = rf.node,
        #                sample.fraction = sample.fraction) #floor(ncol(train)/3)
        # }
        # else{
        #   rf <- ranger(formula = y~., data = train, num.trees = nTree,replace = F,
        #                importance = "impurity",mtry = mtry,
        #                min.node.size = rf.node)
        # }

        # if(typeof(rf) != "character"){
        #   final.rf$method <- "RF"
        #   final.rf$Vars <- importance(rf)
        #
        #   predicted <- rep(NA,nrow(data))
        #   names(predicted) <- rownames(data)
        #   predicted[match(rownames(test),names(predicted))] <- predict(rf,test)$predictions
        #   final.rf$CPE <-  as.numeric(phcpe(coxph(testSurv ~predict(rf,test)$predictions),out.ties=out.ties))
        #   final.rf$CI <- as.numeric(concordance(coxph(testSurv ~predict(rf,test)$predictions,
        #                                               test))$concordance)
        #   #as.numeric(mean((predict(rf,test)-test$y)^2))
        #   final.rf$predicted <- predicted
        #   final.rf$formula <- formula
        #
        #   if(rf_gbm.save){
        #     final.rf$RF <- rf
        #   }
        #
        # }
        # else{
        #   final.rf <- NULL
        # }


        return(final.rf)
      }
      if(save){save(RF,file = paste0(pathResults,"/",studyType,"_RF.Rdata"))
        RF <- NULL}
    }


    ###########
    ### SVM ###
    ###########

    if("SVM" %in% method) {
      print("SVM SELECTED")
      final.svm <- list()
      SVM <- foreach(run=1:runs) %dopar% {

        set.seed(run)
        cat(paste("Run : ", run,"\n",sep=""), file=paste0(studyType,"SVM_log.txt"), append=TRUE)

        rm.samples <- sample(1:nrow(data), ceiling(nrow(data)*1/3),replace = FALSE)
        train <- data[-rm.samples,]
        test <- data[rm.samples,]

        if(LT){
          fit <- coxph(Surv(time1,time2,status)~1,data=train)
          train <- train[,-match(c("time1","time2","status"),colnames(train))]
          testSurv <- with(test,Surv(time1,time2,status))
        }
        else{
          fit <- coxph(Surv(time,status)~1,data=train)
          train <- train[,-match(c("time","status"),colnames(train))]
          testSurv <- with(test,Surv(time,status))
        }
        train$y <- residuals(fit,type="martingale")


        SVM.tune <- tune(svm, y~.,  data = train,
                         ranges = list(epsilon = epsilon.svm, cost = cost.svm))

        SVM <- svm(y~.,  data = train,
                   epsilon = as.numeric(SVM.tune$best.parameters[1]),
                   cost = as.numeric(SVM.tune$best.parameters[2]))
        cat(typeof(SVM)!= "character", file=paste0(studyType,"SVM_log.txt"), append=TRUE)

        # apply(SVM$SV,2,mean)

        if(typeof(SVM)!= "character"){
          predicted <- rep(NA,nrow(data))
          names(predicted) <- rownames(data)
          predicted[match(rownames(test),names(predicted))] <- predict(SVM,newdata = test)

          CPE <- as.numeric(phcpe(coxph(testSurv ~predict(SVM,newdata = test)),out.ties=out.ties))
          final.svm$CPE <- CPE
          CI <- as.numeric(concordance(coxph(testSurv ~predict(SVM,newdata = test),
                                             test))$concordance)
          final.svm$CI <- CI

          # coefs <- as.data.frame(apply(SVM$SV,2,mean))
          final.svm$Vars <- apply(SVM$SV,2,function(x){median(abs(x)) }) #as.data.frame(
          final.svm$predicted <- predicted
          final.svm$formula <- formula
          final.svm$method <- "SVM"
          final.svm$family <- family
          if(rf_gbm.save) final.svm$SVM <- SVM
        }

        else{
          final.svm <- NULL
        }

        return(final.svm)
      }
      if(save){
        save(SVM,file = paste0(pathResults,studyType,"_SVM.Rdata"))
        SVM <- NULL}
    }

    ###########
    ### GBM ###
    ###########
    if("GBM" %in% method){
      final.gbm <- list()
      print("GBM SELECTED")
      GBM <- foreach(run=1:runs) %dopar% {

        gbmGrid <- expand.grid(.interaction.depth = interactions,
                               .n.trees = nTree, .shrinkage = shrinkage,
                               .n.minobsinnode = min.node)

        set.seed(run)
        cat(paste("Run : ", run,"\n",sep=""), file=paste0(studyType,"GBM_log.txt"), append=TRUE)

        rm.samples <- sample(1:nrow(data), ceiling(nrow(data)*1/3),replace = FALSE)
        train <- data[-rm.samples,]
        test <- data[rm.samples,]

        BestPerf <- apply(gbmGrid,1,function(x){
          set.seed(21071993)
          if(!LT){
            GBM <- try(withTimeout(gbm(Surv(time,status)~.,data = train,distribution = "coxph",
                                       interaction.depth = x[1],n.trees=as.numeric(x[2]),shrinkage = x[3],n.minobsinnode = x[4],
                                       bag.fraction = 0.5,train.fraction = 1,cv.folds = cv.folds,n.cores = 1),timeout = 500 ),silent=T)

            if(is.character(GBM)){
              return(list("CPE"=0,"CI"=0,"infl"=NA,"predicted"=NA,
                          "bestTreeForPrediction"=NA))
            }

            bestTreeForPrediction <- gbm.perf.noprint(GBM)
            predicted<- predict(GBM,newdata=test,n.trees = bestTreeForPrediction,type="response")
            CPE <- as.numeric(phcpe(coxph(Surv(time,status) ~predicted,data=test),out.ties=out.ties))
            CI <- as.numeric(concordance(coxph(Surv(time,status) ~ predicted,
                                               test))$concordance)
            return(list("CPE"=CPE,"CI"=CI,"infl"=relative.influence.noprint(GBM),"predicted"=predicted,
                        "bestTreeForPrediction"=bestTreeForPrediction,"GBM"=GBM))
          }


          if(LT){

            fit <- coxph(Surv(time1,time2,status)~1,data=train)
            temp <- train
            temp$y <- residuals(fit,type="martingale")
            temp <- temp[,-match(c("time1","time2","status"),colnames(temp))]
            testSurv <- with(test,Surv(time1,time2,status))
            # test$y <- residuals(coxph(Surv(time1,time2,status)~1,data=test),type="martingale")

            GBM <- try(gbm(y~.,data = temp,distribution = "gaussian",
                           interaction.depth = x[1],n.trees=as.numeric(x[2]),shrinkage = x[3],n.minobsinnode = x[4],
                           bag.fraction = 0.5,train.fraction = 1,cv.folds = 5,n.cores = 1),silent=T)

            if(is.character(GBM)){
              return(list("CPE"=0,"CI"=0,"infl"=NA,"predicted"=NA,
                          "bestTreeForPrediction"=NA))
            }

            bestTreeForPrediction <- gbm.perf.noprint(GBM)
            predicted<- predict(GBM,newdata=test,n.trees = bestTreeForPrediction,type="response")
            CPE <- as.numeric(phcpe(coxph(Surv(time1,time2,status) ~predicted,data=test),out.ties=out.ties))
            CI <- as.numeric(concordance(coxph(testSurv ~ predicted,
                                               test))$concordance)

            return(list("CPE"=CPE,"CI"=CI,"infl"=relative.influence.noprint(GBM),"predicted"=predicted,
                        "bestTreeForPrediction"=bestTreeForPrediction,"GBM"=GBM))

          }
          gc()
        })

        # if(!LT){
        if(sum(vapply(BestPerf,"[[","CPE",FUN.VALUE = numeric(1)),na.rm = T)==0){return(NULL)}
        index <- which.max(vapply(BestPerf,"[[","CPE",FUN.VALUE = numeric(1)))

        CI <- BestPerf[[index]]$CI
        final.gbm$CI <- CI
        final.gbm$CPE <- BestPerf[[index]]$CPE


        final.gbm$method <- "GBM"
        final.gbm$Vars <- BestPerf[[index]]$infl

        predicted <- rep(NA,nrow(data))
        names(predicted) <- rownames(data)
        predicted[match(rownames(test),names(predicted))] <- BestPerf[[index]]$predicted

        final.gbm$predicted <- predicted
        final.gbm$bestTreeForPrediction <- BestPerf[[index]]$bestTreeForPrediction
        final.gbm$params <- gbmGrid[index,]
        final.gbm$formula <- formula
        final.gbm$family <- family

        if(rf_gbm.save) {
          BestPerf[[index]]$GBM$data <- NULL
          BestPerf[[index]]$GBM$valid.error <- NULL
          BestPerf[[index]]$GBM$oobag.improve <- NULL
          BestPerf[[index]]$GBM$bag.fraction <- NULL
          BestPerf[[index]]$GBM$interaction.depth <- NULL
          BestPerf[[index]]$GBM$nTrain <- NULL
          BestPerf[[index]]$GBM$response.name <- NULL
          BestPerf[[index]]$GBM$shrinkage <- NULL
          BestPerf[[index]]$GBM$var.monotone <- NULL
          BestPerf[[index]]$GBM$verbose <- NULL
          BestPerf[[index]]$GBM$Terms <- NULL
          BestPerf[[index]]$GBM$cv.error <- NULL
          BestPerf[[index]]$GBM$cv.folds <- NULL
          BestPerf[[index]]$GBM$call <- NULL
          BestPerf[[index]]$GBM$cv.fitted <- NULL
          BestPerf[[index]]$GBM$trees <- BestPerf[[index]]$GBM$trees[1:BestPerf[[index]]$bestTreeForPrediction]
          final.gbm$GBM <- BestPerf[[index]]$GBM
        }
        return(final.gbm)
      }

      if(save){save(GBM,file = paste0(pathResults,"/",studyType,"_GBM.Rdata"))
        GBM <- NULL}
    }


    ##########
    ### NN ###
    ##########
    if("NN" %in% method){
      final.nn <- list()
      print("NN SELECTED")

      if(norm.nn){
        if(LT) rm <- 1:3
        if(!LT) rm <- 1:2
        maxs <- apply(data[,-rm], 2, max)
        mins <- apply(data[,-rm], 2, min)

        data[,-rm] <- as.data.frame(scale(data[-rm], center = mins, scale = maxs - mins))
        data.save = T
      }

      NN <- foreach(run=1:runs) %dopar% {

        set.seed(run)
        cat(paste("Run : ", run,"\n",sep=""), file=paste0(studyType,"NN_log.txt"), append=TRUE)

        rm.samples <- sample(1:nrow(data), ceiling(nrow(data)*1/3),replace = FALSE)
        train <- data[-rm.samples,]
        test <- data[rm.samples,]

        if(LT){
          fit <- coxph(Surv(time1,time2,status)~1,data=train)
          train <- train[,-match(c("time1","time2","status"),colnames(train))]
          testSurv <- with(test,Surv(time1,time2,status))
        }
        else{
          fit <- coxph(Surv(time,status)~1,data=train)
          train <- train[,-match(c("time","status"),colnames(train))]
          testSurv <- with(test,Surv(time,status))
        }
        train$y <- residuals(fit,type="martingale")
        covs <- colnames(train)
        f <- as.formula(paste("y ~", paste(covs[!covs %in% "y"], collapse = " + ")))

        if(is.null(layers)){
          n <- ncol(train)
          sub.layers <- round(c(n*4/5,n*3/4,n*2/3,n/2,n/3,n/4,n/5))}
        else sub.layers <- layers
        out <- lapply(sub.layers,function(x){
          nn <- neuralnet(f,data=train,hidden=x,linear.output=T)
          pred <- predict(nn,test)

          predicted <- rep(NA,nrow(data))
          names(predicted) <- rownames(data)
          predicted[match(rownames(test),names(predicted))] <- pred

          CPE <- as.numeric(phcpe(coxph(testSurv ~pred),out.ties=out.ties))
          CI <- as.numeric(concordance(coxph(testSurv ~pred,
                                             test))$concordance)
          return(list(nnet=nn,CPE=CPE,CI=CI))
        })
        index <- which.max(vapply(out, "[[", "CPE", FUN.VALUE = numeric(1)))
        nn <- out[[index]]$nnet
        pred <- predict(nn,test)

        predicted <- rep(NA,nrow(data))
        names(predicted) <- rownames(data)
        predicted[match(rownames(test),names(predicted))] <- pred

        CPE <- as.numeric(phcpe(coxph(testSurv ~pred),out.ties=out.ties))
        final.nn$CPE <- CPE
        CI <- as.numeric(concordance(coxph(testSurv ~pred,
                                           test))$concordance)
        final.nn$CI <- CI

        Vars <- garson(nn)$data
        Vars <- Vars[match(colnames(train),Vars$x_names),]
        Vars <- Vars[complete.cases(Vars),1]
        names(Vars) <- colnames(train)[-match("y",colnames(train))]
        final.nn$Vars <- Vars
        final.nn$predicted <- predicted
        final.nn$formula <- formula
        final.nn$method <- "NN"
        final.nn$layer <- sub.layers[index]
        final.nn$family <- family

        if(rf_gbm.save) final.nn$NN <- nn

        return(final.nn)
      }

      if(save){save(NN,file = paste0(pathResults,"/",studyType,"_NN.Rdata"))
        NN <- NULL}
    }
  }


  #################
  # binomial onco #
  #################

  if(family %in% c("binomial","gaussian")){

    # appropriate formula
    temp.formula <- as.formula(formula)
    if(!length(as.list(temp.formula)) == 3){
      stop("ERROR : Response must have three components only. eg: 'y ~.'.")
    }

    colnames(data)[match(as.list(temp.formula)[2],colnames(data))] <- "y"
    data <- data[,c(match("y",colnames(data)),
                    c(1:ncol(data))[-match("y",colnames(data))])]


    ########## LASSO #############
    final.lasso <- list()
    if("LASSO" %in% method) {
      print("LASSO SELECTED")
      LASSO <- foreach(run=1:runs) %dopar% {

        ### BUILD TRAINING AND TESTING SET ###
        set.seed(run)
        cat(paste("Run : ", run,"\n",sep=""), file=paste0(studyType,"LASSO_log.txt"), append=TRUE)
        # split data
        rm.samples <- sample(1:nrow(data), ceiling(nrow(data)*1/3),replace = FALSE)
        train <- data[-rm.samples,]
        test <- data[rm.samples,]

        cv.fit <- cv.glmnet(x = as.matrix(train[,-match("y",colnames(train))]), y = train[,"y"],family = family,
                            nfolds = 5,alpha=1)
        fit <- glmnet(x = as.matrix(train[,-match("y",colnames(train))]), y = train[,"y"],family = family,alpha=1)

        if(family == "gaussian"){
          pred <- predict(cv.fit,newx = as.matrix(test[,-match("y",colnames(test))]), s="lambda.min")
          predicted <- rep(NA,nrow(data))
          names(predicted) <- rownames(data)
          predicted[match(rownames(pred),names(predicted))] <- as.numeric(pred)
          mse <- mean((test$y - pred)^2)
        }

        if(family == "binomial"){
          pred <- predict(cv.fit,newx = as.matrix(test[,-match("y",colnames(test))]), s="lambda.min",type="class")
          predicted <- rep(NA,nrow(data))
          names(predicted) <- rownames(data)
          predicted[match(rownames(pred),names(predicted))] <- as.numeric(pred)
          CI <- sum(as.numeric(pred) == test$y) / nrow(test)
        }
        Coefficients <- coef(fit, s = cv.fit$lambda.min)
        coefs <- Coefficients[-1,1]

        final.lasso$formula <- formula
        final.lasso$fit <- coefs
        final.lasso$predicted <- predicted
        # final.lasso$means <- lasso.fit$means --> need to think of how to do validation ...
        final.lasso$method <- "LASSO"
        if(family == "binomial") final.lasso$CI <- CI
        if(family == "gaussian") final.lasso$MSE <- MSE
        final.lasso$family <- family
        return(final.lasso)
      }
      if(save){
        save(LASSO,file = paste0(pathResults,"/",studyType,"_LASSO.Rdata"))}

    }


    ########## RIDGE #############
    final.ridge <- list()
    if("RIDGE" %in% method) {
      print("RIDGE SELECTED")
      RIDGE <- foreach(run=1:runs) %dopar% {

        ### BUILD TRAINING AND TESTING SET ###
        set.seed(run)
        cat(paste("Run : ", run,"\n",sep=""), file=paste0(studyType,"RIDGE_log.txt"), append=TRUE)
        # split data
        rm.samples <- sample(1:nrow(data), ceiling(nrow(data)*1/3),replace = FALSE)
        train <- data[-rm.samples,]
        test <- data[rm.samples,]

        cv.fit <- cv.glmnet(x = as.matrix(train[,-match("y",colnames(train))]), y = train[,"y"],family = family,
                            nfolds = 5,alpha=0)
        fit <- glmnet(x = as.matrix(train[,-match("y",colnames(train))]), y = train[,"y"],family = family,alpha=0)

        if(family == "gaussian"){
          pred <- predict(cv.fit,newx = as.matrix(test[,-match("y",colnames(test))]), s="lambda.min")
          predicted <- rep(NA,nrow(data))
          names(predicted) <- rownames(data)
          predicted[match(rownames(pred),names(predicted))] <- as.numeric(pred)
          mse <- mean((test$y - pred)^2)
        }

        if(family == "binomial"){
          pred <- predict(cv.fit,newx = as.matrix(test[,-match("y",colnames(test))]), s="lambda.min",type="class")
          predicted <- rep(NA,nrow(data))
          names(predicted) <- rownames(data)
          predicted[match(rownames(pred),names(predicted))] <- as.numeric(pred)
          CI <- sum(as.numeric(pred) == test$y) / nrow(test)
        }
        Coefficients <- coef(fit, s = cv.fit$lambda.min)
        coefs <- Coefficients[-1,1]

        final.ridge$formula <- formula
        final.ridge$fit <- coefs
        final.ridge$predicted <- predicted
        # final.ridge$means <- ridge.fit$means --> need to think of how to do validation ...
        final.ridge$method <- "RIDGE"
        if(family == "binomial") final.ridge$CI <- CI
        if(family == "gaussian") final.ridge$MSE <- MSE
        final.ridge$family <- family
        return(final.ridge)
      }
      if(save){
        save(RIDGE,file = paste0(pathResults,"/",studyType,"_RIDGE.Rdata"))}

    }


    ########## ENET #############
    final.enet <- list()
    if("ENET" %in% method) {
      print("ENET SELECTED")
      ENET <- foreach(run=1:runs) %dopar% {

        ### BUILD TRAINING AND TESTING SET ###
        set.seed(run)
        cat(paste("Run : ", run,"\n",sep=""), file=paste0(studyType,"ENET_log.txt"), append=TRUE)
        # split data
        rm.samples <- sample(1:nrow(data), ceiling(nrow(data)*1/3),replace = FALSE)
        train <- data[-rm.samples,]
        test <- data[rm.samples,]

        alphas <- seq(0,1,0.05)

        allCV <- lapply(1:length(alphas),function(x){
          cv.fit <- cv.glmnet(x = as.matrix(train[,-match("y",colnames(train))]), y = train[,"y"],family = family,
                              nfolds = 5,alpha=alphas[x])
          fit <- glmnet(x = as.matrix(train[,-match("y",colnames(train))]), y = train[,"y"],family = family,alpha=alphas[x])

          if(family == "gaussian"){
            pred <- predict(cv.fit,newx = as.matrix(test[,-match("y",colnames(test))]), s="lambda.min")
            mse <- mean((test$y - pred)^2)
            return(list("cv.fit"=cv.fit,"fit"=fit,"mse"=mse,"predicted"=pred))
          }

          if(family == "binomial"){
            pred <- predict(cv.fit,newx = as.matrix(test[,-match("y",colnames(test))]), s="lambda.min",type="class")
            CI <- sum(as.numeric(pred) == test$y) / nrow(test)
            return(list("cv.fit"=cv.fit,"fit"=fit,"CI"=CI,"predicted"=pred))
          }
        })

        if(family == "gaussian"){
          best.run <- which.min(sapply(allCV,"[[","mse"))
          cv.fit <- allCV[[best.run]]$cv.fit
          fit <- allCV[[best.run]]$fit
          mse <- allCV[[best.run]]$mse
          # predicted <- allCV[[best.run]]$predicted
        }

        if(family == "binomial"){
          best.run <- which.max(sapply(allCV,"[[","CI"))
          cv.fit <- allCV[[best.run]]$cv.fit
          fit <- allCV[[best.run]]$fit
          CI <- allCV[[best.run]]$CI
          # predicted <- allCV[[best.run]]$predicted
        }

        predicted <- rep(NA,nrow(data))
        names(predicted) <- rownames(data)
        predicted[match(rownames(allCV[[best.run]]$predicted),names(predicted))] <- as.numeric(allCV[[best.run]]$predicted)

        Coefficients <- coef(fit, s = cv.fit$lambda.min)
        coefs <- Coefficients[-1,1]

        final.enet$formula <- formula
        final.enet$fit <- coefs
        final.enet$predicted <- predicted
        # final.lasso$means <- lasso.fit$means --> need to think of how to do validation ...
        final.enet$method <- "LASSO"
        if(family == "binomial") final.enet$CI <- CI
        if(family == "gaussian") final.enet$MSE <- MSE
        final.enet$family <- family
        return(final.enet)
      }
      if(save){
        save(ENET,file = paste0(pathResults,"/",studyType,"_ENET.Rdata"))}

    }


    ########## GBM #############
    final.gbm <- list()
    if("GBM" %in% method) {
      print("GBM SELECTED")
      GBM <- foreach(run=1:runs) %dopar% {

        ### BUILD TRAINING AND TESTING SET ###
        set.seed(run)
        cat(paste("Run : ", run,"\n",sep=""), file=paste0(studyType,"GBM_log.txt"), append=TRUE)
        # split data
        rm.samples <- sample(1:nrow(data), ceiling(nrow(data)*1/3),replace = FALSE)
        train <- data[-rm.samples,]
        test <- data[rm.samples,]

        gbmGrid <- expand.grid(.interaction.depth = interactions,
                               .n.trees = nTree, .shrinkage = shrinkage,
                               .n.minobsinnode = min.node)

        BestPerf <- apply(gbmGrid,1,function(x){
          set.seed(21071993)

          if(family == "binomial"){
            GBM <- try(withTimeout(gbm(y~.,data = train,distribution = "bernoulli",
                                       interaction.depth = x[1],n.trees=as.numeric(x[2]),shrinkage = x[3],n.minobsinnode = x[4],
                                       bag.fraction = 0.5,train.fraction = 1,cv.folds = cv.folds,n.cores = 1),timeout = 500 ),silent=T)

            if(is.character(GBM)){
              return(list("CI"=0,"infl"=NA,"predicted"=NA,
                          "bestTreeForPrediction"=NA))
            }

            bestTreeForPrediction <- gbm.perf.noprint(GBM)
            predicted<- predict(GBM,newdata=test,n.trees = bestTreeForPrediction,type="response")
            CI <- sum(round(predicted) == test$y) / nrow(test)
            return(list("CI"=CI,"infl"=relative.influence.noprint(GBM),"predicted"=predicted,
                        "bestTreeForPrediction"=bestTreeForPrediction,"GBM"=GBM))
          }
          gc()
        })

        if(sum(vapply(BestPerf,"[[","CI",FUN.VALUE = numeric(1)),na.rm = T)==0){return(NULL)}
        index <- which.max(vapply(BestPerf,"[[","CI",FUN.VALUE = numeric(1)))

        CI <- BestPerf[[index]]$CI
        final.gbm$CI <- CI

        final.gbm$method <- "GBM"
        final.gbm$Vars <- BestPerf[[index]]$infl

        predicted <- rep(NA,nrow(data))
        names(predicted) <- rownames(data)
        predicted[match(rownames(test),names(predicted))] <- BestPerf[[index]]$predicted

        final.gbm$predicted <- predicted
        final.gbm$bestTreeForPrediction <- BestPerf[[index]]$bestTreeForPrediction
        final.gbm$params <- gbmGrid[index,]
        final.gbm$formula <- formula
        final.gbm$family <- family

        if(rf_gbm.save) {
          BestPerf[[index]]$GBM$data <- NULL
          BestPerf[[index]]$GBM$valid.error <- NULL
          BestPerf[[index]]$GBM$oobag.improve <- NULL
          BestPerf[[index]]$GBM$bag.fraction <- NULL
          BestPerf[[index]]$GBM$interaction.depth <- NULL
          BestPerf[[index]]$GBM$nTrain <- NULL
          BestPerf[[index]]$GBM$response.name <- NULL
          BestPerf[[index]]$GBM$shrinkage <- NULL
          BestPerf[[index]]$GBM$var.monotone <- NULL
          BestPerf[[index]]$GBM$verbose <- NULL
          BestPerf[[index]]$GBM$Terms <- NULL
          BestPerf[[index]]$GBM$cv.error <- NULL
          BestPerf[[index]]$GBM$cv.folds <- NULL
          BestPerf[[index]]$GBM$call <- NULL
          BestPerf[[index]]$GBM$cv.fitted <- NULL
          BestPerf[[index]]$GBM$trees <- BestPerf[[index]]$GBM$trees[1:BestPerf[[index]]$bestTreeForPrediction]
          final.gbm$GBM <- BestPerf[[index]]$GBM
        }
        return(final.gbm)
      }

      if(save){save(GBM,file = paste0(pathResults,"/",studyType,"_GBM.Rdata"))
        GBM <- NULL}
    }


  }


  ################
  stopCluster(cl)
  ################

  ##################################

  OUTPUT <- list()
  if(!is.null(LASSO)){OUTPUT$LASSO <- LASSO}
  if(!is.null(RIDGE)){OUTPUT$RIDGE <- RIDGE}
  if(!is.null(ENET)){OUTPUT$ENET <- ENET}
  if(!is.null(RF)){OUTPUT$RF <- RF}
  if(!is.null(GBM)){OUTPUT$GBM <- GBM}
  if(!is.null(SVM)){OUTPUT$SVM <- SVM}
  if(!is.null(NN)){OUTPUT$NN <- NN}

  if(data.save == T) {OUTPUT$data <-data}
  if(!is.null(OUTPUT)){
    return(OUTPUT)}
  else{return(0)}

}
