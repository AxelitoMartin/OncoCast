#' predIncomingSurv
#'
#' This function let's user predict the genetic risk score of new/incoming patients. This function
#' makes use of the output of the OncoCast function and an inputed data set similar to the one use to generate
#' those results. The user can retrieve the genetic risk score of those new patients and their predicted survival curves.
#'
#' @param OC_object Output of the OncoCast function.
#' @param new.data New data set containing the information of the incoming patients. Should be a dataframe
#' with patients as rows and features as columns.
#' @param surv.print A numeric vector indicating the patients for which the user wishes to print the predicted
#' survival curves.
#' @param riskRefit The refitted cox proportional hazard model with the risk score as predictor.
#'
#' @return data.out : The data frame inputted in the function with an additional column giving the predicted risk
#' score of the incoming patients.
#' @return RiskHist : A histogram of the distribution of the risk scores of patients in the given dataset.
#' @return IncKM : An interactive Kaplan-Meier plot of the selected patients in the surv.print argument.
#' @export
#' @examples library(OncoCast)
#' test <- OncoCast(data=survData,formula = Surv(time,status)~.,
#' method = "LASSO",runs = 50,
#' save = F,nonPenCol = NULL,cores =2)
#' results <- getResults_OC(OC_object=test$LASSO,data=survData,
#' numGroups=5,cuts=c(0.2,0.4,0.6,0.8),
#' geneList="NULL",mut.data=T)
#' new.data <- as.data.frame(matrix(rbinom(5*20,1,0.5),nrow=20,ncol = 5))
#' colnames(new.data) <- c("ImpCov1","ImpCov2","ImpCov3","ImpCov4","Cov7")
#' rownames(new.data) <- paste0("Incoming",1:20)
#' Incoming <- predIncomingSurv(test$LASSO,new.data=new.data,
#' surv.print = c(5,10,15),riskRefit = results$RiskRefit)


predIncomingSurv <- function(OC_object,new.data,surv.print= NULL,riskRefit){

  OC_object <- Filter(Negate(is.null), OC_object)
  # get all information needed from the oncocast object
  # 1. risk
  final.pred <- sapply(OC_object,"[[","predicted")
  ori.risk <- apply(final.pred,1,function(x){
    mean(as.numeric(x),na.rm = TRUE)
  })


  ### FOR PENALIZED REG ###

  if(OC_object[[1]]$method %in% c("LASSO","RIDGE","ENET")){
    # 2. Fits
    LassoFits <- t(sapply(OC_object,"[[","fit"))
    LassoFits[is.na(LassoFits)] <- 0

    ################################
    features <- colnames(LassoFits)

    dums <- apply(new.data,2,function(x){anyNA(as.numeric(as.character(x)))})
    if(sum(dums)){
      tmp <- new.data %>%
        select(which(dums)) %>%
        fastDummies::dummy_cols(remove_first_dummy = T) %>%
        select(-one_of(names(which(dums))))
      new.data <- as.data.frame(cbind(
        new.data %>% select(-one_of(names(which(dums)))),
        tmp
      ) %>% mutate_all(as.character) %>%
        mutate_all(as.numeric)
      )
      warning("Character variables were transformed to dummy numeric variables. If you didn't have any character variables make sure all columns in your input data are numeric. The transformed data will be saved as part of the output.")
    }

    if(!all(is.na(match(colnames(new.data),features)))){
      matched.genes <- c(na.omit(match(colnames(new.data),features)))
      new.dat <- new.data[,which(!is.na(match(colnames(new.data),features)))]

      ## ADD ALL MISSING GENES TO BE ALL zero ##
      missing <- features[which(is.na(match(features,colnames(new.dat))))]
      to.add <- as.data.frame(matrix(0L,nrow=nrow(new.dat),ncol=length(missing)))
      colnames(to.add) <- missing
      rownames(to.add) <- rownames(new.dat)
      new.dat <- as.data.frame(cbind(new.dat,to.add))

      new.dat <- new.dat[,match(features,colnames(new.dat))]

      #############################################

      all.pred <- lapply(1:nrow(LassoFits),function(x){
        ### Subset to the coefs of that cv ###
        coefs <- LassoFits[x,LassoFits[x,] != 0]
        new.temp <- select(new.dat,names(coefs))

        ## substract mean mutation rate of TRAINING SET !!!###
        new.x <- new.temp - rep(OC_object[[x]]$means[match(names(coefs),names(OC_object[[x]]$means))], each = nrow(new.temp))
        cal.risk.test <- drop(as.matrix(new.x) %*% coefs)
        return(cal.risk.test)
      })
    }
    else{
      stop("No gene overlapped be sure they are correctly matched.")
    }
  }

  ### For GBM ###
  if(OC_object[[1]]$method  %in% c("GBM","RF","SVM")){

    if(OC_object[[1]]$method == "GBM") features <- OC_object[[1]]$GBM$var.names
    if(OC_object[[1]]$method == "RF") features <- OC_object[[1]]$RF$forest$independent.variable.names
    if(OC_object[[1]]$method == "SVM") features <- names(OC_object[[1]]$Vars)

    dums <- apply(new.data,2,function(x){anyNA(as.numeric(as.character(x)))})
    if(sum(dums)){
      tmp <- new.data %>%
        select(which(dums)) %>%
        fastDummies::dummy_cols(remove_first_dummy = T) %>%
        select(-one_of(names(which(dums))))
      new.data <- as.data.frame(cbind(
        new.data %>% select(-one_of(names(which(dums)))),
        tmp
      ) %>% mutate_all(as.character) %>%
        mutate_all(as.numeric)
      )
      warning("Character variables were transformed to dummy numeric variables. If you didn't have any character variables make sure all columns in your input data are numeric. The transformed data will be saved as part of the output.")
    }

    if(!all(is.na(match(colnames(new.data),features)))){
      matched.genes <- c(na.omit(match(colnames(new.data),features)))
      new.dat <- new.data[,which(!is.na(match(colnames(new.data),features)))]

      ## ADD ALL MISSING GENES TO BE ALL zero ##
      missing <- features[which(is.na(match(features,colnames(new.dat))))]
      to.add <- as.data.frame(matrix(0L,nrow=nrow(new.dat),ncol=length(missing)))
      colnames(to.add) <- missing
      rownames(to.add) <- rownames(new.dat)
      new.dat <- as.data.frame(cbind(new.dat,to.add))

      new.dat <- new.dat[,match(features,colnames(new.dat))]

      if(OC_object[[1]]$method == "GBM") {
        all.pred <- lapply(OC_object,function(x){
          predict(x$GBM,newdata=new.dat,
                  n.trees = x$bestTreeForPrediction,
                  type="response")
        })}

      if(OC_object[[1]]$method == "RF") {
        all.pred <- lapply(OC_object,function(x){
          predict(x$RF,new.dat)$predictions
        })}

      if(OC_object[[1]]$method == "SVM") {
        all.pred <- lapply(OC_object,function(x){
          predict(x$SVM,new.dat)
        })}

    }
    else{
      stop("No gene overlapped be sure they are correctly matched.")
    }
  }

  all.pred <- do.call("cbind",all.pred)
  Risk <- apply(all.pred,1,mean)
  names(Risk) <- rownames(new.dat)
  # Risk.all <- as.matrix(coefs) %*% as.matrix(t(new.dat))
  # Risk <- apply(Risk.all,2,mean)
  #new.data$Risk <- Risk
  ##########################################
  ori.risk.range <- range(ori.risk)
  new.data$OncoCastRiskScore <- rescale(Risk, to = c(0, 10), from = ori.risk.range) #WithOriginal
  #new.data$rescaledRisk <- rescale(new.data$Risk, to = c(0, 10), from = range(new.data$Risk, na.rm = TRUE, finite = TRUE))
  RiskHistogram.new <- ggplot(new.data, aes(x = OncoCastRiskScore, y = ..density..)) +
    geom_histogram(show.legend = FALSE, aes(fill=..x..),
                   breaks=seq(min(new.data$OncoCastRiskScore,na.rm = T), max(new.data$OncoCastRiskScore,na.rm = T))) +#, by=20/nrow(new.data))) +
    geom_density(show.legend = FALSE) +
    theme_minimal() +
    labs(x = "Average risk score", y = "Density") +
    scale_fill_gradient(high = "red", low = "green")

  #return(list("RiskHistogram.new"=RiskHistogram.new,"out.data"=new.data))


  ####################################################
  ## Creat survival curves for patients of interest ##
  ####################################################

  if(!is.null(surv.print)){
    mut <- new.data[surv.print,]
    colnames(mut)[ncol(mut)] <- "RiskScore"

    allSurvs <- data.frame(nrow= 5)

    for(j in 1:nrow(mut)){
      survival.probs <- as.data.frame(matrix(nrow=6,ncol=15))
      rownames(survival.probs) <- c("Patient","Surv","Lower","Upper","Time","OncoRiskScore")
      surv.temp <- survfit(riskRefit, newdata = mut[j,])
      for(i in 1:ncol(survival.probs)){
        survival.probs[,i] <- try(c(rownames(mut)[j],as.numeric(summary(surv.temp, times = (i*3-3))$surv),
                                    round(summary(surv.temp, times = (i*3-3))$lower,digits=2),
                                    round(summary(surv.temp, times = (i*3-3))$upper,digits=2),
                                    i*3-3,as.numeric(mut$RiskScore[j])),silent=T)
      }
      allSurvs <- cbind(allSurvs,survival.probs)
    }
    allSurvs <- allSurvs[,-1]

    a <- list(
      autotick = FALSE,
      dtick = 6,
      tickcolor = toRGB("black")
    )

    t.survival.probs <- as.data.frame(t(allSurvs))
    for(k in 2:ncol(t.survival.probs)){
      t.survival.probs[,k] <- as.numeric(as.character(t.survival.probs[,k]))
    }
    y <- list(
      title = "Survival Probability"
    )
    IndSurvKM <- plot_ly(t.survival.probs, x = ~Time, y = ~Surv, name = ~Patient, type = 'scatter',
                         mode = 'lines+markers',hoverinfo="hovertext",color = ~Patient,
                         hovertext = ~paste("Genetic Risk Score :",round(OncoRiskScore,digits=3))
    ) %>% layout(yaxis = y,xaxis = ~a) %>%
      layout(xaxis = list(title = paste0("Time (Months)"), showgrid = TRUE),showlegend = FALSE) %>%
      add_ribbons(data = t.survival.probs,
                  ymin = ~Lower,
                  ymax = ~Upper,
                  line = list(color = 'rgba(7, 164, 181, 0.05)'),
                  fillcolor = 'rgba(7, 164, 181, 0.2)',
                  name = "Confidence Interval") #%>% layout(showlegend = FALSE)

  }
  else{IndSurvKM = NULL}
  return(list("data.out" = new.data,"RiskHist"=RiskHistogram.new,"IncKM" = IndSurvKM,"survivalEst"=t.survival.probs))
}
