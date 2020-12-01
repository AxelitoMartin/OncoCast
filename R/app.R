#' app UI functions
#'
#' This file contains the ui functions for running the internal application of
#' OncoCast.
#' @return Web-based application
#' @keywords application
#' @export
#' @examples library(OncoCast)
#' @import
#' shiny
#' shinydashboard
#' survival
#' gbm
#' ranger
#' foreach
#' doParallel
#' ggplot2
#' reshape2
#' scales
#' pheatmap
#' survminer
#' rmarkdown
#' knitr
#' CPE
#' mlr
#' dtplyr
#' fastDummies
#' NeuralNetTools
#' @importFrom dplyr select filter mutate group_by rename summarise arrange
#' @importFrom neuralnet neuralnet
#' @importFrom e1071 tune svm
#' @importFrom plotly plot_ly layout toRGB add_ribbons
#' @importFrom penalized optL1 optL2 predict

ui <- shinydashboard::dashboardPage(
  shinydashboard::dashboardHeader(title = "Prognostic machine learning for survival analysis",titleWidth = 500),

  shinydashboard::dashboardSidebar(width = 200,
                   selectInput("choice", label = "App type:",
                               choices = c("Create OncoCast run","Load OncoCast run")),
                   shinydashboard::sidebarMenuOutput("menu")
  ),

  shinydashboard::dashboardBody(
    shinydashboard::tabItems(


      ### intro ###
      shinydashboard::tabItem(tabName = "intro",
              h1("Introduction to OncoCast"),
              mainPanel(width = 14,
                        tags$b(tags$u(h2("Development context"))),

                        h4("OncoCast is an analytical tool for prediction and risk stratification of survival
                   statistical problems, includind delayed entry problems. This framework relies both on
                   classical and more advanced machine learning methods and cross-validation to generate an
                   ensemble learning framework. This enables an unbiased and optimized assesment of the predicted risk.
                   More information about OncoCast can be found ",a("here",href="https://axelitomartin.github.io/software.html"),
                           " and be found as an R package with additional functionalities on ",a("github",href="https://github.com/AxelitoMartin/OncoCast"),".",
                           "We orginally developed this framework as a tool to analyse genomic data for cancer patients,
                   our findings were published in this ",
                           a('manuscript', href = 'https://ascopubs.org/doi/10.1200/PO.18.00307'),"."),

                        tags$b(tags$u(h2("Usage"))),

                        h4("This online application enables users to either perform a simple OncoCast run by
                   uploading their own data, or load a pre-existing run and exploring its properties.
                   Given the latter implies an a-priori knowledge of the method (most likely R users)
                   we will not cover it in this abstract but if you wish to read more about it I will
                   refer you to the ",
                           a("vignette",href = "https://axelitomartin.github.io/docs/OncoCast_vignette.html"),".",
                           "These options can be toggled in the scroll down at the top of the tab list (top left).",
                           "To perform an OncoCast run select the following tab (OncoCast run). There you will be
                   asked to choose the model of your choice, we recommend penalized regression methods (LASSO,
                   RIDGE or ENET) for uncorrelated and sparse data, and tree based methods (RF, GBM) otherwise.
                   Next the user will be asked to select the number of 'runs', this is the number
                   of cross-validations that will be performed. A greater number of runs will lead to a more
                   robust run at the cost of computing time.
                   In order for the application to remain relatively simple we ask that users upload the data
                   in a comma separated value ('.csv'). Moreover it is necessary to clearly state the survival variables as ",
                           tags$b("time"), " and ", tags$b("status")," (or ",tags$b("time1,time2,status"),"for
                           left-truncated data).",
                           "Also not that all categorical variables will be recoded as dummy binary variables due to
                           the complexitity of the method. Example datasets are provided below.
                           The following tabs will enable the user to explore the results and perform
                            risk group stratification, validation with external datasets and individual survival predictions
                           in the following tabs."),

                        h3("Download example data to run OncoCast:"),
                        downloadLink('downloadAppData', 'Download'),
                        h3('Download example validation data:'),
                        downloadLink('downloadValAppData', 'Download'),

                        h4("Contact axel.steph.martin@gmail.com for any questions or comments.")
              )

      ),

      ### Run oncocast ###
      shinydashboard::tabItem(tabName = "run",
              h1("Choose OncoCast parameters"),
              mainPanel(
                radioButtons("method", "Select a method:",
                             c("LASSO" = "LASSO",
                               "RIDGE" = "RIDGE",
                               "ENET" = "ENET",
                               "RF" = "RF",
                               "GBM" = "GBM",
                               "SVM"="SVM",
                               "NN"="NN")),

                sliderInput("runs", label = "Select number of runs", min = 30,
                            max = 200, value = 50,step = 10),
                conditionalPanel(
                  condition = 'input.method == "RF" || input.method == "GBM"', #c("RF","GBM")
                  sliderInput("nTree", label = "Select number of trees per run", min = 100,
                              max = 1000, value = 500,step = 100)
                ),
                fileInput("data.run", "Choose a time data CSV file (make sure you set all the other parameters first):",
                          accept = c(
                            "text/csv",
                            "text/comma-separated-values,text/plain",
                            ".csv")
                ),
                radioButtons("LT", "Is the data left-truncated:",
                             c("No" = "0",
                               "Yes" = "1")),
                htmlOutput("launchRun")
              )),

      ### Load the file up your choosing ###
      shinydashboard::tabItem(tabName = "load",
              h1("Choose OncoCast file"),
              mainPanel(
                fileInput("file", label = ""),
                htmlOutput("Load"),
                fileInput("data", "Choose the corresponding CSV file",
                          accept = c(
                            "text/csv",
                            "text/comma-separated-values,text/plain",
                            ".csv")
                ),
                htmlOutput("dataLoad")
              )),


      # First tab content
      shinydashboard::tabItem(tabName = "fixed",
              sidebarLayout(
                sidebarPanel(
                  width=12,
                  h2("Prognostic Performance"),
                  textInput("cuts", label = h3("Enter quantiles for stratification (e.g.: 25,50,75).
                                               By default kmeans clustering will be used."), value = "")#,
                  # submitButton("Submit")
                ),
                mainPanel(width=12,
                          htmlOutput("RiskHeader"),
                          plotOutput("RiskHistogram"),
                          tableOutput("RiskSummary"),
                          htmlOutput("CIHeader"),
                          tableOutput("CPE"),
                          htmlOutput("RefitHeader"),
                          tableOutput("RefitRisk")
                          #tableOutput("ClinRefit")
                )
              )
      ),


      shinydashboard::tabItem(tabName = "strat",
              sidebarLayout(
                sidebarPanel(width = 12,
                             h2("Risk Group Stratification")#,
                             # downloadButton("report", "Generate report")
                ),
                mainPanel(width = 12,
                          htmlOutput("predRiskText"),
                          htmlOutput("KMText"),
                          plotOutput("KM",height="500px"),
                          tableOutput("SurvSum")
                )
              )
      ),

      shinydashboard::tabItem(tabName = "gene",
              h2("Exploratory interactive gene plots"),
              sidebarLayout(
                sidebarPanel(width = 12,
                             textInput("GeneListRisk",
                                       "Find gene(s) : ",
                                       value = "")#,
                             # submitButton("Submit")
                ),
                mainPanel(
                  htmlOutput("VolcanoHeader"),width = 12,
                  plotly::plotlyOutput("effectPlot"),
                  htmlOutput("heatmaphead"),
                  plotOutput("binMap"),
                  plotOutput("contMap")
                )
              )
      ),

      shinydashboard::tabItem(tabName = "valid",
              h2("Validation of the model on an exterior dataset"),
              sidebarLayout(
                sidebarPanel(width = 12,
                             fileInput("newdata", "Upload a new file, including the outcome (must be named the same as in the oncocast run)",
                                       accept = c(
                                         "text/csv",
                                         "text/comma-separated-values,text/plain",
                                         ".csv")
                             ),
                             textInput("custFormula",
                                       "If your validation file has a different survival formula then the original data (e.g: time1,time2,status): ",
                                       value = "")#,
                             # submitButton("Submit")

                ),
                mainPanel(width = 12,
                          plotOutput("newRisk"),
                          plotOutput("KM.val")
                )
              )
      ),



      shinydashboard::tabItem(tabName = "patient",
              h2("Predict the survival of individual patients"),
              sidebarLayout(
                sidebarPanel(width = 12,
                             fileInput("newdataind", "Upload a new file (Note that only the 10 first patients will be printed)",
                                       accept = c(
                                         "text/csv",
                                         "text/comma-separated-values,text/plain",
                                         ".csv")
                             )
                ),
                mainPanel(width = 12,
                          plotOutput("newRisk.ind"),
                          plotly::plotlyOutput("KM.ind")#,
                          # tableOutput("survprob.ind")
                )
              )
      ),


      shinydashboard::tabItem(tabName = "DL",
              h2("Download summary report"),
              sidebarLayout(
                sidebarPanel(width = 12,
                             radioButtons('typeDL',"Choose format (note some features will be lost when converted to PDF):",
                                          c("HTML"="HTML",
                                            "PDF (not currently available)"="PDF"))
                ),
                mainPanel(width = 12,
                          downloadLink('report', 'Download'),

                          h3("You can also download the OncoCast run raw output here:"),
                          downloadLink('downloadOncoCastrun', 'Download Run')
                )
              )
      )


    )
  )
)

##################################


#' app server function
#'
#' This file contains the  server functions for running the internal application of
#' OncoCast.
#' @param input input for server
#' @param output output for server
#' @return Web-based application
#' @keywords application
#' @export
#' @examples library(OncoCast)
#' @import
#' shiny
#' shinydashboard
#' survival
#' gbm
#' ranger
#' foreach
#' doParallel
#' ggplot2
#' reshape2
#' scales
#' pheatmap
#' survminer
#' rmarkdown
#' knitr
#' CPE
#' mlr
#' fastDummies
#' NeuralNetTools
#' @importFrom dplyr select filter mutate group_by rename summarise arrange
#' @importFrom neuralnet neuralnet
#' @importFrom e1071 tune svm
#' @importFrom plotly plot_ly layout toRGB add_ribbons
#' @importFrom penalized optL1 optL2 predict

server <- function(input, output) {

  options(shiny.maxRequestSize=500*1024^2)


  output$menu <- shinydashboard::renderMenu({

    if(input$choice=="Load OncoCast run"){
      my_list = list(
        shinydashboard::menuItem("Abstract", tabName = "intro", icon = icon("atlas")),
        shinydashboard::menuItem("Loading Preview", tabName = "load", icon = icon("arrow-alt-circle-down")),
        shinydashboard::menuItem("Dashboard", tabName = "fixed", icon = icon("bar-chart")),
        shinydashboard::menuItem("Risk Group Stratification", tabName = "strat", icon = icon("scissors")),
        shinydashboard::menuItem("Gene View", tabName = "gene", icon = icon("chart-bar")),
        shinydashboard::menuItem("Validation", tabName = "valid", icon = icon("check-circle")),
        shinydashboard::menuItem("Patient View", tabName = "patient", icon = icon("user-circle-o")),
        shinydashboard::menuItem("Download Summary", tabName = "DL",icon = icon("download")))
    }

    if(input$choice=="Create OncoCast run"){
      my_list = list(
        shinydashboard::menuItem("Abstract", tabName = "intro", icon = icon("atlas")),
        shinydashboard::menuItem("OncoCast Run", tabName = "run", icon = icon("gears")),
        shinydashboard::menuItem("Dashboard", tabName = "fixed", icon = icon("bar-chart")),
        shinydashboard::menuItem("Risk Group Stratification", tabName = "strat", icon = icon("scissors")),
        shinydashboard::menuItem("Gene View", tabName = "gene", icon = icon("chart-bar")),
        shinydashboard::menuItem("Validation", tabName = "valid", icon = icon("check-circle")),
        shinydashboard::menuItem("Patient View", tabName = "patient", icon = icon("user-circle-o")),
        shinydashboard::menuItem("Download Summary", tabName = "DL",icon = icon("download"))
      )
    }
    shinydashboard::sidebarMenu(my_list)

  })

  ##### Loading corresponding data as csv ... #####
  dat <- reactive({
    if(input$choice=="Load OncoCast run") temp <- read.csv(input$data$datapath,row.names = 1)
    else temp <- read.csv(input$data.run$datapath,row.names = 1)

    dums <- apply(temp,2,function(x){anyNA(as.numeric(as.character(x)))})
    if(sum(dums) > 0){
      tmp <- temp %>%
        select(which(dums)) %>%
        fastDummies::dummy_cols(remove_first_dummy = T) %>%
        select(-one_of(names(which(dums))))
      data <- as.data.frame(cbind(
        temp %>% select(-one_of(names(which(dums)))),
        tmp
      ) %>% mutate_all(as.character) %>%
        mutate_all(as.numeric)
      )
      return(data)
    }
    return(temp)
  })

  ##### Make sure file was loaded correctly and that all#####

  OC_object <- reactive({
    if(input$choice=="Load OncoCast run"){
      OC_object <-load_object(input$file$datapath)
      OC_object <- Filter(Negate(is.null), OC_object)
    }

    else{
      if(input$LT == "0") out <- OncoCast(dat(),formula= Surv(time,status)~., method = input$method,
                                          runs = as.numeric(input$runs),cores = 1,save = F,
                                          nonPenCol = NULL,nTree=as.numeric(input$nTree),rf_gbm.save = T,
                                          epsilon.svm = seq(0,0.2,0.1), cost.svm = 2^(2:4))
      else out <- OncoCast(dat(),formula= Surv(time1,time2,status)~., method = input$method,
                           runs = as.numeric(input$runs),cores = 1,save = F,
                           nonPenCol = NULL,nTree=as.numeric(input$nTree),rf_gbm.save = T,
                           epsilon.svm = seq(0,0.2,0.1), cost.svm = 2^(2:4))
      OC_object <- out[[1]]
      # if(!is.null(out[[2]])) dat() <- out[[2]] #reactive({out[[2]]})
      out <- NULL
    }

    return(OC_object)})

  output$launchRun <- renderText({ paste("<h4> <u> <font color=\"red\"><b>",
                                         "Once you have loaded the data and selected all the parameters, please select the following tab (Dashboard).","</b></font> </u> </h4>",
                                         "Note that depending on your sample size and the method selected the run
                                        time may take some time, please be patient."
  )})

  #####
  observeEvent(input$file,if(!is.null(input$file)){
    method <- OC_object()[[1]]$method

    if(method %in% c("ENET","LASSO","RIDGE","GBM","RF","SVM")){
      output$Load <- renderText({ paste("<h2> <u> <font color=\"black\"><b>",
                                        "LOAD SUCCESFUL. You selected a",method,"OncoCast run",
                                        "</b></font> </u> </h2>")
      })
    }
  }
  else{output$Load <- renderText({ paste("<h2> <u> <font color=\"red\"><b>",
                                         "LOAD FAILURE! Please choose another file.",
                                         "</b></font> </u> </h2>") })}
  )


  ##### Run get results and store them under appropriate name outputs #####

  results <- reactive({
    return(Progno.tab(OC_object(),data=dat()))
  })


  observeEvent(input$data,if(!is.null(input$data)){
    # output$dataLoad <- renderText({ paste(colnames(dat()))})
    # output$dataLoad <- renderTable({ dat()})
    if(is.list(results())){
      output$dataLoad <- renderText({ paste("<h2> <u> <font color=\"black\"><b>",
                                            "Corresponding data read in.",
                                            "</b></font> </u> </h2>")})}

    else{output$dataLoad <- renderText({ paste("<h2> <u> <font color=\"red\"><b>",
                                               "The data loaded does not correspond the oncocast run loaded.",
                                               "</b></font> </u> </h2>")})}
  })

  ##### Now that everything is looking good proceed #####

  ##### Tab2 #####
  # headers
  output$VariableHeader <- renderText({ paste("<h3> <u> <font color=\"black\"><b>","Prognostic performance", "</b></font> </u> </h3>") })
  output$RiskHeader <- renderText({paste("<h4> <u> <font color=\"black\"><b>","Histogram of predicted risk score", "</b></font> </u> </h4>")})
  output$CIHeader <- renderText({paste("<h4> <u> <font color=\"black\"><b>","Concordance probability estimate in predicting overall survival", "</b></font> </u> </h4>")})
  output$RefitHeader <- renderText({paste("<h4> <u> <font color=\"black\"><b>","Cox regression estimates and significance P-values", "</b></font> </u> </h4>")})

  # results

  output$RiskHistogram <- renderPlot({results()$RiskHistogram})
  output$RiskSummary <- renderTable({results()$RiskScoreSummary},rownames = TRUE)
  output$CPE <- renderTable({results()$CPE},rownames = TRUE)
  output$RefitRisk <- renderTable({
    temp <- summary(results()$RiskRefit)$coefficients
    colnames(temp) <- c("Coefficient","Hazard Ratio","S.E.","Z-value","P-value")
    temp[1,] <- round(temp[1,],digits=4)
    if(temp[,5] < 0.00001) temp[,5] <- "<0.00001"
    return(temp)
  },rownames = TRUE)


  #### TAB3 ####
  # headers
  output$predRiskText <- renderText({ paste("<h3> <u> <font color=\"black\"><b>","Risk group stratification", "</b></font> </u> </h3>") })
  output$KMText <- renderText({ paste("<h4> <u> <font color=\"black\"><b>","Kaplan-Meier plot of overall survival",
                                      "</b></font> </u> </h4>") })


  # results
  strat <- reactive({riskStrat(OC_object(),dat(),
                               cuts = input$cuts,
                               results()$risk.raw,results()$LT,plotQuant=1)})
  output$KM <- renderPlot({strat()$KM})
  output$SurvSum <- renderTable({strat()$survivalTable},rownames = T)


  #### TAB4 ####
  output$VolcanoHeader <- renderText({ paste("<h3> <u> <font color=\"black\"><b>","Gene effect & selection",
                                             "</b></font> </u> </h3>") })

  gene.plot <- reactive({gene.view(OC_object(),dat(),
                                   geneList = gsub(" ","",unlist(strsplit(input$GeneListRisk, split ="," ))),
                                   results()$LT,results()$risk.raw)})
  output$effectPlot <- plotly::renderPlotly({gene.plot()$selectInflPlot})

  output$heatmaphead <- renderText({ paste("<h3> <u> <font color=\"black\"><b>","Distribution of features by risk",
                                           "</b></font> </u> </h3>") })

  output$binMap <- renderPlot({gene.plot()$heatmap.sorted.bin})

  output$contMap <- renderPlot({gene.plot()$heatmap.sorted.cont})


  ##### VALIDATION ####
  in.dat <- reactive({read.csv(input$newdata$datapath,row.names = 1)})

  v.results <- reactive({

    validate_app(OC_object = OC_object(),
                 LassoFits=gene.plot()$Fits,
                 ori.risk=results()$risk.raw,
                 qts=strat()$rawCuts,
                 in.data= in.dat(),
                 formula=input$custFormula)#OC_object()[[1]]$formula)

  })

  output$newRisk <- renderPlot({v.results()$RiskHistogram.new})
  output$KM.val <- renderPlot({v.results()$KM})



  ##### Individual prediction ####
  in.dat.pred <- reactive({read.csv(input$newdataind$datapath,row.names = 1)})

  ind.results <- reactive({

    predIncomingSurv_app(OC_object = OC_object(),
                         new.data=in.dat.pred(),
                         surv.print= 1:min(10,nrow(in.dat.pred)),
                         riskRefit=results()$RiskRefit)

  })

  output$newRisk.ind <- renderPlot({ind.results()$RiskHist})
  output$KM.ind <- plotly::renderPlotly({ind.results()$IncKM})
  # output$survprob.ind <- renderTable(({ind.results()$survivalEst}))



  ##### DOWNLOAD REPORT ######
  output$report <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = "report.html",
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path("../inst/", "HTML_report.Rmd") #tempdir()
      file.copy("HTML_report.Rmd", tempReport, overwrite = TRUE)

      # possible validation parameters #
      if(is.character(try(v.results()$RiskHistogram.new,silent = T))){
        ValHist <- NA
        ValKM <- NA}
      else{
        ValHist <- v.results()$RiskHistogram.new
        ValKM <- v.results()$KM
      }
      if(is.character(try(ind.results()$IncKM))) IndKM <- NA
      else IndKM=ind.results()$IncKM

      # Set up parameters to pass to Rmd document
      params <- list(
        # data stuff #
        "n"=nrow(dat()),"n.events"=sum(dat()$status),"n.features"=ncol(dat())-2,
        "method"=OC_object()[[1]]$method, "n.CV"=input$runs,
        # summary page #
        "riskHist"=results()$RiskHistogram,"riskSum"=results()$RiskScoreSummary,
        "CPE"=results()$CPE,"riskTable"=summary(results()$RiskRefit)$coefficients,
        "cuts"=unlist(strsplit(input$cuts, split ="," )),"KM"=strat()$KM,
        "survTable"=strat()$survivalTable,"effectPlot"=gene.plot()$selectInflPlot,
        "binMap"=gene.plot()$heatmap.sorted.bin,"contMap"=gene.plot()$heatmap.sorted.cont,
        # optional validation stuff #
        "ValHist"=ValHist,"ValKM"=ValKM,
        "IndKM"=IndKM
        # "ValHist"=v.results()$RiskHistogram.new,"ValKM"=v.results()$KM,
        # "IndKM"=ind.results()$IncKM
      )
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  )

  # oncocast run #
  OncoRun <- reactiveValues()
  observe({
    if(!is.character(try(OC_object(),silent = T)))
      isolate(
        OncoRun <<- OC_object()
      )
  })

  output$downloadOncoCastrun <- downloadHandler(
    filename = function() {
      "OncoCastRun.Rdata"
    },
    content = function(file) {
      save(OncoRun,file = file)
    }
  )

  # data to download #
  output$downloadAppData <- downloadHandler(
    filename = function() {
      "appTest_data.csv"
    },
    content = function(file) {
      data <- read.csv("../inst/extdata/appTest_data.csv",row.names = 1)
      write.csv(data, file)
    }
  )
  output$downloadValAppData <- downloadHandler(
    filename = function() {
      "appTest_valdata.csv"
    },
    content = function(file) {
      data <- read.csv("../inst/extdata/appTest_valdata.csv",row.names = 1)
      write.csv(data, file)
    }
  )
}

shinyApp(ui, server)
