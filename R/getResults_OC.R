##################################################################
##################### ONCOCAST RESULTS ###########################
##################################################################

#' getResults_OC
#'
#' This functions let's the user study the output of the OncoCast function. This function takes as input
#' one of the (or the one) objects returned from the different machine learning algorithms chosen previously.
#' Only one such object can be inputted at a time in the getResults_OC function.
#' @param OC_object A list object outputed by the OncoCast function.
#' @param data A dataframe that corresponds to the data used to generate the OncoCast output.
#' @param numGroups The number of groups to be made when stratifying by risk groups. Options are 2,3 and 4 (for now implementing
#' broader version). Default is 2.
#' @param cuts Numeric vector of the points in the distribution of risk scores where groups will be splitted. Needs to be of length
#' numGroups - 1. eg : c(0.25,0.5,0.75) when numgroups is 4. Default is 0.5.
#' @param geneList Optional character vector of features of potential higher interest. Default is NULL, which
#' leads to using the 5 most frequently selected features.
#' @param mut.data Boolean argument indicating if the user is using mutation predictors (binary data). Default is FALSE.
#' @param plotQuant A numeric entry between 0-1 that defines what proportion of patients will be represented on the Kaplan-Meier plot.
#' Particularly useful when a lot of patients with long survival are censored. Default is 1 (all patients are plotted).
#' @return CPE Summary of the distribution of the concordance probability estimate (see phcpe function) accross all runs. (Recommended)
#' @return CI Summary of the distribution of the concordance index accross all runs. (Depreciated)
#' @return inflPlot Bar plot of frequency of the 20 most selected features.
#' @return topHits Character vector of the top 10 most selected features.
#' @return average.risk Average predicted risk score for each patient in the data set.
#' @return RiskScore The average.risk output rescaled between 0-10.
#' @return data.out The data that was used for the analysis.
#' @return selectInflPlot Volcano plot of the selection frequency against the average mean coefficient of each feature accross all runs.
#' Note that this plot is interactive.
#' @return RiskRefit Refitted cox proportional hazard model with the predicted average risk score as the continuous covariate.
#' @return RiskHistogram Histogram of the density distribution of the average predicted risk score. Note it has been rescaled from 0 to 10
#' for simplicity.
#' @return Fits Data frame reporting the coefficients found for each feature at each run.
#' @return time.type Time unit used. Options are Days or Months.
#' @return RiskScoreSummary Distribution summary of the average predicted risk score.
#' @return KM_Plot Kaplan-Meier plot stratified by risk group.
#' @return SurvSum Summary of survival metrics per risk group.
#' @return RiskGroup The risk group assigned to each patients (ordered in the order of the input dataset).
#' @return mut_Plot Mutation distribution by features bar plot.
#' @return PieChart Interactive pie chart of the mutation dsitribution using either the most frequently selected features or
#' the manually inputted gene list. Each of the pies represent one of the risk groups.
#' @return GenesUsed Character vector of the features used to make the pie charts.
#' @keywords Results
#' @export
#' @examples library(OncoCast)
#' test <- OncoCast(data=survData,formula = Surv(time,status)~.,
#' method = "LASSO",runs = 50,
#' save = F,nonPenCol = NULL,cores =2)
#' results <- getResults_OC(OC_object=test$LASSO,data=survData,
#' numGroups=5,cuts=c(0.2,0.4,0.6,0.8),
#' geneList=NULL,mut.data=T)



getResults_OC <- function(OC_object,data,cuts=NULL,geneList=NULL,mut.data=F,plotQuant=1,...){

  OC_object <- Filter(Negate(is.null), OC_object)
  method <- OC_object[[1]]$method
  formula <- OC_object[[1]]$formula

  dums <- apply(data,2,function(x){anyNA(as.numeric(as.character(x)))})
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
    warning("Character variables were transformed to dummy numeric variables. If you didn't have any character variables make sure all columns in your input data are numeric. The transformed data will be saved as part of the output.")
  }

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

  return(outputSurv(OC_object,data,method,geneList,numGroups,cuts,plotQuant,mut.data,LT,...))

}


#############################################################



#############################################################

#' outputSurv
#'
#' This functions let's the user study the basic output of the OncoCast function. This function takes as input
#' one of the (or the one) objects returned from the different machine learning algorithms chosen previously.
#' Only one such object can be inputted everytime in the outputSummary function.
#' @param OC_object A list object outputed by the VariableSelection function.
#' @param data A dataframe that corresponds to the data used to generate the OncoCast output.
#' @param method Method used to generate the OC_object (e.g.: "LASSO").
#' @param geneList Optional argument enabling the user to use a particular list of features of interest.
#' @param numGroups Integer argument for the number of groups to be made in stratified Kaplan-Meier
#' plot. Recommend to note make groups smaller than 10 patients.
#' @param cuts Vector argument of decimal numbers between 0 and 1 representing the quantiles where the groups will be formed.
#' @param plotQuant Decimal number between 0 and 1, for the proportion of patients to be shown on the Kaplan-Meier plot.
#' @param mut.data Boolean indicating if the set of predictors are binary variables.
#' @return ciSummary Summary of the distribution of the concordance index accross all runs.
#' @return inflPlot Bar plot of frequency of the 20 most selected features.
#' @return topHits Character vector of the top 10 most selected features.
#' @return average.risk Average predicted risk score for each patient in the data set.
#' @return scaled.risk The average.risk output rescaled between 0-10.
#' @return data.out The data that was used for the analysis.
#' @return selectInflPlot Volcano plot of the selection frequency against the average mean coefficient of each feature accross all runs.
#' Note that this plot is interactive.
#' @return RiskRefit Refitted cox proportional hazard model with the predicted average risk score as the continuous covariate.
#' @return RiskHistogram Histogram of the density distribution of the average predicted risk score. Note it has been rescaled from 0 to 10
#' for simplicity.
#' @return Fits Data frame reporting the coefficients found for each feature at each run.
#' @return time.type Time unit used. Options are Days or Months.
#' @return RiskScoreSummary Distribution summary of the average predicted risk score.
#' @keywords Results
#' @export
#' @examples library(OncoCast2)
#' test <- OncoCast(data=survData,formula = Surv(time,status)~.,
#' family = "cox",method = c("LASSO"),runs = 50,
#' sampling = "cv",save = F,nonPenCol = NULL,cores =2)
#' OC_object <- test$LASSO
#' data <- survData
#' numGroups <- 5
#' cuts <- c(0.2,0.4,0.6,0.8)
#' out.test <- outputSurv(OC_object,data,method="GBM",
#'                        geneList="NULL",numGroups,cuts,mut.data=T)


outputSurv <- function(OC_object,data,method,geneList=NULL,numGroups=2,cuts=0.5,plotQuant=1,mut.data=F,LT,...){

  args <- list()
  OC_object <- Filter(Negate(is.null), OC_object)
  MD <- 12

  ConcordanceIndex <- as.data.frame(as.vector(unlist(vapply(OC_object, "[[", "CI",FUN.VALUE = numeric(1)))))
  summary.CI <- round(as.data.frame(c(quantile(ConcordanceIndex[,1],c(0.1,0.25,0.5,0.75,0.9),na.rm = T))),digits = 2)
  colnames(summary.CI) <- "Concordance Index"
  rownames(summary.CI) <- c("Lower 10%","1st Quarter","Median","3rd Quarter","Upper 10%")
  CI.BP <- as.data.frame(t(summary.CI))

  CPEIndex <- as.data.frame(as.vector(unlist(vapply(OC_object, "[[", "CPE",FUN.VALUE = numeric(1)))))
  summary.CPE <- round(as.data.frame(c(quantile(CPEIndex[,1],c(0.1,0.25,0.5,0.75,0.9),na.rm = T))),digits = 2)
  colnames(summary.CPE) <- "Concordance probability estimate"
  rownames(summary.CPE) <- c("Lower 10%","1st Quarter","Median","3rd Quarter","Upper 10%")
  CPE <- as.data.frame(t(summary.CPE))

  if(method %in% c("LASSO","RIDGE","ENET")){
    allCoefs <- t(sapply(OC_object,"[[","fit"))
    allCoefs[is.na(allCoefs)] <- 0

    meanCoefs <- apply(allCoefs,2,function(x){mean(x,na.rm = TRUE)})
    selectFreq <- apply(allCoefs,2,function(x){
      length(which(x!=0))/length(x)
    })

    if(length(selectFreq[selectFreq > 0.5]) > 2) {
      topHits <- names(selectFreq[order(selectFreq,decreasing = T)])[1:sum(selectFreq > 0.5)] }

    else topHits <- names(selectFreq[order(selectFreq,decreasing = T)])[1:10]
    ## get mu freq
    if(LT) data.temp <- data[,-c(1:3)]
    if(!LT) data.temp <- data[,-c(1:2)]
    MutationFrequency <- apply(data.temp,2,function(x){
      sum(x)/length(x)
    })

    resultsAll <- as.data.frame(cbind(meanCoefs,selectFreq,MutationFrequency))
    colnames(resultsAll) <- c("MeanCoefficient","SelectionFrequency","MutationFrequency")
    rownames(resultsAll) <- names(meanCoefs)
    resultsAll <- resultsAll[complete.cases(resultsAll),]
    resultsAll$GeneName <- rownames(resultsAll)
    resultsAll$MutationFrequency2 <- cut(resultsAll$MutationFrequency, c(0,0.10,0.20,0.40))

    if(length(geneList)!=0){
      m <- resultsAll[match(geneList,rownames(resultsAll)), ]

      a <- list(
        x = m$MeanCoefficient,
        y = m$SelectionFrequency,
        text = rownames(m),
        xref = "x",
        yref = "y",
        showarrow = TRUE,
        arrowhead = 7,
        ax = 20,
        ay = -40
      )

      selectInflPlot <- plot_ly(data = resultsAll, x = ~MeanCoefficient, y = ~SelectionFrequency,
                                text = ~paste('Gene :',GeneName,
                                              '</br> Hazard Ratio :',round(exp(MeanCoefficient),digits=2)),
                                mode = "markers",size = ~MutationFrequency,color = ~MeanCoefficient) %>%
        layout(title ="Volcano Plot",annotations = a)
    }

    else{
      selectInflPlot <- plot_ly(data = resultsAll, x = ~MeanCoefficient, y = ~SelectionFrequency,
                                text = ~paste('Gene :',GeneName,
                                              '</br> Hazard Ratio :',round(exp(MeanCoefficient),digits=2)),
                                mode = "markers",size = ~MutationFrequency,color = ~MeanCoefficient) %>%
        layout(title ="Volcano Plot")
    }

  }

  if(method %in% c("GBM","RF","SVM","NN")) {
    #### Variables ####
    if(LT) Variables <- colnames(data)[-c(1:3)]
    if(!LT) Variables <- colnames(data)[-c(1:2)]

    imp <- sapply(OC_object,"[[","Vars")
    # if(method == "NN") rownames(imp) <- Variables
    mean.imp <- apply(scale(imp),1,mean)

    # for 10 top medians make boxplot
    topHits <- sort(mean.imp,decreasing = T)[1:pmin(length(Variables),20)]
    topHitsMat <- t(imp[match(names(topHits),rownames(imp)),])
    topHitsMat.melt <- melt(topHitsMat)[,-1]
    colnames(topHitsMat.melt) <- c("Gene","Rank")

    selectInflPlot <- ggplot(topHitsMat.melt, aes(x=Gene,y = Rank)) +
      geom_boxplot(outlier.colour = "red", outlier.shape = 1,colour = "#3366FF") +
      #geom_jitter(width = 0.05) +
      theme_grey() +
      labs(title = "Variable use in the forest", subtitle = "Using genetic data only",y = "Importance") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))

    allCoefs <- mean.imp
    resultsAll <- NULL
  }

  final.pred <- sapply(OC_object,"[[","predicted")
  average.risk <- apply(final.pred,1,function(x){
    mean(as.numeric(x),na.rm = TRUE)
  })
  average.risk[which(is.na(average.risk))] <- NA
  to <- c(0,10)
  from <- range(average.risk, na.rm = TRUE, finite = TRUE)
  RiskScore <- (as.numeric(average.risk)-from[1])/diff(from)*diff(to)+to[1]

  summary.RiskScore <- round(as.data.frame(c(quantile(RiskScore,c(0.1,0.25,0.33,0.5,0.66,0.75,0.9),na.rm = TRUE))),digits = 2)
  colnames(summary.RiskScore) <- "Risk Score"
  rownames(summary.RiskScore) <- c("Lower 10%","1st Quarter","1st Tertile","Median","2nd Tertile","3rd Quarter","Upper 10%")

  # refit coxph model with average risk as covariate
  meanRS <- mean(RiskScore)
  if(LT) refit.risk <- coxph(Surv(data$time1,data$time2,data$status)~RiskScore)
  if(!LT) refit.risk <- coxph(Surv(data$time,data$status)~RiskScore)
  Risk <- as.data.frame(RiskScore)

  RiskHistogram <- ggplot(Risk, aes(x = RiskScore, y = ..density..)) +
    geom_histogram(show.legend = FALSE, aes(fill=..x..),
                   breaks=seq(min(Risk$RiskScore,na.rm = T), max(Risk$RiskScore,na.rm = T), by=0.25)) +
    geom_density(show.legend = FALSE) +
    theme_minimal() +
    labs(x = "Average risk score", y = "Density") +
    scale_fill_gradient(high = "red", low = "green")


  ##### univariate volcano plot #####
  if(LT == F){
    uni <- as.data.frame(t(apply(data[,-c(1:2)],2,function(x){
      fit <- coxph(Surv(data$time,data$status)~x)
      out <- as.numeric(summary(fit)$coefficients[c(1,5)])
      return(out)
    })))
    colnames(uni) <- c("Coefficient","Pvalue")
    uni$Feature <- colnames(data)[-c(1:2)]
  }

  if(LT == T){
    uni <- as.data.frame(t(apply(data[,-c(1:3)],2,function(x){
      fit <- coxph(Surv(data$time1,data$time2,data$status)~x)
      out <- as.numeric(summary(fit)$coefficients[c(1,5)])
      return(out)
    })))
    colnames(uni) <- c("Coefficient","Pvalue")
    uni$Feature <- colnames(data)[-c(1:3)]
  }


  uniVolcano <- plot_ly(data = uni, x = ~Coefficient, y = ~-log10(Pvalue),
                        text = ~paste('Feature :',Feature,
                                      '<br> Coefficient :',round(exp(Coefficient),digits=3)),
                        mode = "markers") %>%
    layout(title ="Volcano Plot")


  ##### heatmap #####
  average.risk.order <- order(average.risk,decreasing = F)
  risk.ordered <- average.risk[average.risk.order]

  if(LT) uniques <- apply(data[,-c(1:3)],2,function(x){length(unique(x))})
  else uniques <- apply(data[,-c(1:2)],2,function(x){length(unique(x))})

  ## for binary ##
  heatmap.sorted.bin <- NULL

  if(sum(uniques == 2) > 2){
    genes <- names(uniques[which(uniques == 2)])
    if(method %in% c("LASSO","RIDGE","ENET") ) keep <- names(sort(selectFreq[match(genes,names(selectFreq))],decreasing = T)[1:pmin(15,length(genes))])
    if(method %in% c("GBM","RF","SVM","NN")) keep <- names(sort(mean.imp[match(genes,names(mean.imp))],decreasing = T)[1:pmin(15,length(genes))])

    if(!is.null(geneList)){
      if(length(geneList) < 15){
        keep <- c(geneList,keep[-c((15-length(geneList)):15)])
      }
    }

    data.sub <- data[,match(keep,colnames(data))]
    data.ordered <- t(as.matrix(data.sub[average.risk.order,]))
    heatmap.sorted.bin <- pheatmap(data.ordered,cluster_rows = F,cluster_cols = F,show_rownames=T,show_colnames = F,
                                   color = c("cyan","black"),border_color = "black",silent =T,legend=F,
                                   fontsize = 15)
  }

  ## for continuous ##
  heatmap.sorted.cont <- NULL

  if(sum(uniques > 2) > 2){
    genes <- names(uniques[which(uniques > 2)])
    if(method %in% c("LASSO","RIDGE","ENET") ) keep <- names(sort(selectFreq[match(genes,names(selectFreq))],decreasing = T)[1:pmin(15,length(genes))])
    if(method %in% c("GBM","RF","SVM","NN")) keep <- names(sort(mean.imp[match(genes,names(mean.imp))],decreasing = T)[1:pmin(15,length(genes))])

    if(!is.null(geneList)){
      if(length(geneList) < 15){
        keep <- c(geneList,keep[-c((15-length(geneList)):15)])
      }
    }

    data.sub <- data[,match(keep,colnames(data))]
    data.ordered <- t(as.matrix(data.sub[average.risk.order,]))
    data.ordered.new <- t(apply(data.ordered,1,function(x){
      to <- c(0,1)
      from <- range(x, na.rm = TRUE, finite = TRUE)
      new <- (as.numeric(x)-from[1])/diff(from)*diff(to)+to[1]
      return(new)
    }))
    heatmap.sorted.cont <- pheatmap(data.ordered.new,cluster_rows = F,cluster_cols = F,show_rownames=T,show_colnames = F,silent =T)
  }



  ######### RISK STRATIFICATION ############
  if(is.null(cuts)){
    # apply kmeans and take smallest #
    dists <- c()
    set.seed(21071993)
    for(i in 2:5){
      temp <- kmeans(average.risk,centers = i)
      dists[i] <- temp$tot.withinss + 2*i*nrow(temp$centers)
    }
    numGroups <- which.min(dists)
    temp <- kmeans(average.risk,centers = numGroups)
    riskGroupTemp <- temp$cluster
    qts <- c()
    count <- 1
    for(i in order(temp$centers)) {
      qts[count] <- as.numeric(average.risk[names(which.max(average.risk[riskGroupTemp == i]))])
      count = count + 1
    }
    riskGroup <- c()
    map <- order(temp$centers)
    names(map) <- as.character(1:numGroups)
    riskGroup <- names(map[match(riskGroupTemp,map)])
  }
  else {
    numGroups <- length(cuts)+1
    qts <- quantile(average.risk,cuts)
    riskGroup <- c()
    for(i in 1:length(average.risk)){
      temp <- average.risk[i] - qts
      if(length(match(names(which.max(temp[temp<0])),names(qts))) == 1) riskGroup[i] <- match(names(which.max(temp[temp<0])),names(qts))
      else riskGroup[i] <- numGroups
    }
    if(numGroups == 2){riskGroup[is.na(riskGroup)] <- 1}
  }

  data$RiskGroup <- riskGroup
  data$RiskGroup <- factor(data$RiskGroup, levels = c(1:numGroups) )

  if(LT == TRUE) {
    fit0 <- coxph(Surv(time1,time2,status) ~ RiskGroup,data=data,
                  na.action=na.exclude)
    if(max(data$time2) > 1000){
      timeType = "Days"
      intercept = 5*365}
    else{
      timeType = "Months"
      intercept = 5*12}
  }
  if(LT == FALSE) {
    fit0 <- coxph(Surv(time,status) ~ RiskGroup,data=data,
                  na.action=na.exclude)
    if(max(data$time) > 1000){
      timeType = "Days"
      intercept = 5*365}
    else{
      timeType = "Months"
      intercept = 5*12}
  }

  log.test.pval <- as.vector(summary(fit0)[10][[1]])[3]
  CI <- as.numeric(as.vector(summary(fit0)[14])[[1]][1])
  if(LT) limit <- as.numeric(quantile(data$time2,plotQuant))
  if(!LT) limit <- as.numeric(quantile(data$time,plotQuant))

  # if(!exists("surv.median.line")) surv.median.line <- "hv"
  if(is.null(args[['surv.median.line']])) surv.median.line <- "hv"
  if(LT) {KM <- ggsurvplot(survfit(Surv(time1,time2,status) ~ RiskGroup,data=data, conf.type = "log-log"),conf.int  = TRUE,surv.median.line = surv.median.line,
                           data = data,xlim=c(0,limit),break.time.by = 6,risk.table = T) + xlab("Time (Months)") +
    labs(title = paste("Kaplan Meier Plot (p-value : " ,round(log.test.pval,digits =5),")",sep=""))}
  if(!LT){KM <- ggsurvplot(survfit(Surv(time,status) ~ RiskGroup,data=data, conf.type = "log-log"),conf.int  = TRUE,surv.median.line = surv.median.line,
                           data = data,xlim=c(0,limit),break.time.by = 6,risk.table = T) + xlab("Time (Months)") +
    labs(title = paste("Kaplan Meier Plot (p-value : " ,round(log.test.pval,digits =5),")",sep=""))}



  survivalGroup <- as.data.frame(matrix(nrow=numGroups,ncol=4))
  rownames(survivalGroup) <- 1:numGroups
  colnames(survivalGroup) <- c("MedianOS","95%CI","1Ysurvival","3Ysurvival")
  # for each group find closest value to median
  if(timeType == "Months"){YR1 <- 1*12;YR3 <- 3*12}
  if(timeType == "Days"){YR1 <- 1*365;YR3 <- 3*365}
  for(i in 1:numGroups){
    if(LT == TRUE){NewObject <- with(data[data$RiskGroup == i,],Surv(time1,time2,status))}
    if(LT == FALSE){NewObject <- with(data[data$RiskGroup == i,],Surv(time,status))}
    Fit <- survfit(NewObject ~ 1,data=data[data$RiskGroup == i,], conf.type = "log-log")
    # med.index <- which.min(abs(Fit$surv-0.5))
    YR1.index <- which.min(abs(Fit$time-YR1))
    YR3.index <- which.min(abs(Fit$time-YR3))
    survivalGroup[i,] <- c(as.numeric(round(summary(Fit)$table[7],digits=2)),
                           paste0("(",as.numeric(round(summary(Fit)$table[8],digits=2)),",",
                                  as.numeric(round(summary(Fit)$table[9],digits=2)),")"),
                           paste0(round(Fit$surv[YR1.index],digits=2)," (",
                                  round(Fit$lower[YR1.index],digits=2),",",
                                  round(Fit$upper[YR1.index],digits=2),")"),
                           paste0(round(Fit$surv[YR3.index],digits=2)," (",
                                  round(Fit$lower[YR3.index],digits=2),",",
                                  round(Fit$upper[YR3.index],digits=2),")"))
  }


  ################## Mutational profiles #################
  if(method %in% c("GBM","RF","SVM","NN")){topHits <- names(topHits)}
  if(mut.data){

    sub.dat <- data[,match(c(topHits),colnames(data))]
    #mut.table <- group_by(sub.dat,RiskGroup) %>% summarise(Mutations = n())
    mutDistrib <- do.call("cbind",apply(sub.dat,2,function(x){
      temp <- as.data.frame(cbind(x,riskGroup))
      temp2 <- group_by(temp,riskGroup) %>% summarise(Mutations = sum(as.numeric(as.character(x))))
      return(temp2[,2]/nrow(data))
    }))
    colnames(mutDistrib) <- colnames(sub.dat)
    mutDistrib$Risk <- as.character(1:numGroups)
    rownames(mutDistrib) <- as.character(1:numGroups)
    melted.prop <- melt(mutDistrib)
    colnames(melted.prop) <- c("Risk","Gene","Proportion")

    mutplot <- ggplot(melted.prop, aes(x = Gene, y = Proportion, fill = Risk)) +
      geom_bar(stat = "identity", position=position_dodge()) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      labs(title = "Proportion of mutations for each genes per risk and method", subtitle = "Using genetic data only")

  }

  else {
    mutplot <- NULL
  }

  return(list("CPE"=CPE,"CI" = CI.BP,"risk.raw"=average.risk,"scaled.risk"=RiskScore,
              "RiskHistogram"=RiskHistogram,"RiskScoreSummary"=as.data.frame(t(summary.RiskScore)),
              "RiskRefit"=refit.risk,"rawCuts"= as.numeric(qts),
              "uniVolcano"=uniVolcano,"topHits" = topHits,
              "selectInflPlot" = selectInflPlot,"Fits"=allCoefs,
              "heatmap.sorted.bin"=heatmap.sorted.bin,"heatmap.sorted.cont"=heatmap.sorted.cont,
              "RiskGroup"=riskGroup,"KM"=KM,"Strat_pval"=log.test.pval,"survivalTable"=survivalGroup,
              "mut_Plot" = mutplot,"resultsAll" = resultsAll))
}
