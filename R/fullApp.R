library(shiny)
library(shinydashboard)

# source("./R/appFunctions.R")

ui <- dashboardPage(
  dashboardHeader(title = "Prognostic models in metastatic lung adenocarcinoma (BETA)",titleWidth = 400),

  dashboardSidebar(width = 200,
                   selectInput("choice", label = "App type:",
                               choices = c("Create OncoCast run","Load OncoCast run")),
                   # sidebarMenu(
                   sidebarMenuOutput("menu")
                   # )

  ),

  dashboardBody(
    tabItems(

      ### Load the file up your choosing ###
      tabItem(tabName = "load",
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

      ### Run oncocast ###
      tabItem(tabName = "run",
              h1("Choose OncoCast parameters"),
              mainPanel(
                radioButtons("dist", "Distribution type:",
                             c("Normal" = "norm",
                               "Uniform" = "unif",
                               "Log-normal" = "lnorm",
                               "Exponential" = "exp"))
              )),

      # First tab content
      tabItem(tabName = "fixed",
              sidebarLayout(
                sidebarPanel(
                  width=12,
                  h2("Prognostic Performance"),
                  textInput("cuts", label = h3("Enter quantiles for stratification (e.g.: 25,50,75)"), value = "50")#,
                  # submitButton("Submit")
                ),
                mainPanel(width=12,
                          htmlOutput("RiskHeader"),
                          plotOutput("RiskHistogram"),
                          tableOutput("RiskSummary"),
                          htmlOutput("CIHeader"),
                          tableOutput("CI"),
                          htmlOutput("RefitHeader"),
                          tableOutput("RefitRisk")
                          #tableOutput("ClinRefit")
                )
              )
      ),


      tabItem(tabName = "strat",
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

      tabItem(tabName = "gene",
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
                  plotlyOutput("effectPlot"),
                  htmlOutput("heatmaphead"),
                  plotOutput("binMap"),
                  plotOutput("contMap")
                )
              )
      ),

      tabItem(tabName = "valid",
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



      tabItem(tabName = "patient",
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
                          plotlyOutput("KM.ind")#,
                          # tableOutput("survprob.ind")
                )
              )
      )


    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

  options(shiny.maxRequestSize=500*1024^2)


  output$menu <- renderMenu({

    if(input$choice=="Load OncoCast run"){
      my_list = list(menuItem("Loading Preview", tabName = "load", icon = icon("arrow-alt-circle-down")),
                     menuItem("Dashboard", tabName = "fixed", icon = icon("bar-chart")),
                     menuItem("Risk Group Stratification", tabName = "strat", icon = icon("scissors")),
                     menuItem("Gene View", tabName = "gene", icon = icon("chart-bar")),
                     menuItem("Validation", tabName = "valid", icon = icon("check-circle")),
                     menuItem("Patient View", tabName = "patient", icon = icon("user-circle-o")))
    }

    if(input$choice=="Create OncoCast run"){
      my_list = list(
        menuItem("OncoCast Run", tabName = "run", icon = icon("gears")),
        menuItem("Dashboard", tabName = "fixed", icon = icon("bar-chart")),
        menuItem("Risk Group Stratification", tabName = "strat", icon = icon("scissors")),
        menuItem("Gene View", tabName = "gene", icon = icon("chart-bar")),
        menuItem("Validation", tabName = "valid", icon = icon("check-circle")),
        menuItem("Patient View", tabName = "patient", icon = icon("user-circle-o"))
      )
    }
    sidebarMenu(my_list)

  })

  ##### Make sure file was loaded correctly and that all#####

  OC_object <- reactive({
    OC_object <-load_object(input$file$datapath)
    OC_object <- Filter(Negate(is.null), OC_object)
    return(OC_object)})


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


  ##### Loading corresponding data as csv ... #####
  ##### Run get results and store them under appropriate name outputs #####

  dat <- reactive({read.csv(input$data$datapath,row.names = 1)})
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
  output$CIHeader <- renderText({paste("<h4> <u> <font color=\"black\"><b>","Concordance index in predicting overall survival", "</b></font> </u> </h4>")})
  output$RefitHeader <- renderText({paste("<h4> <u> <font color=\"black\"><b>","Cox regression estimates and significance P-values", "</b></font> </u> </h4>")})

  # results

  output$RiskHistogram <- renderPlot({results()$RiskHistogram})
  output$RiskSummary <- renderTable({results()$RiskScoreSummary})
  output$CI <- renderTable({results()$CI})
  output$RefitRisk <- renderTable({
    temp <- summary(results()$RiskRefit)$coefficients
    colnames(temp) <- c("Coefficient","Hazard Ratio","S.E.","Z-value","P-value")
    temp[1,] <- round(temp[1,],digits=4)
    if(temp[,5] < 0.00001) temp[,5] <- "<0.00001"
    return(temp)
  })


  #### TAB3 ####
  # headers
  output$predRiskText <- renderText({ paste("<h3> <u> <font color=\"black\"><b>","Risk group stratification", "</b></font> </u> </h3>") })
  output$KMText <- renderText({ paste("<h4> <u> <font color=\"black\"><b>","Kaplan-Meier plot of overall survival",
                                      "</b></font> </u> </h4>") })


  # results
  strat <- reactive({riskStrat(OC_object(),dat(),
                               cuts = unlist(strsplit(input$cuts, split ="," )),
                               results()$risk.raw,results()$LT,plotQuant=1)})
  output$KM <- renderPlot({strat()$KM})
  output$SurvSum <- renderTable({strat()$survivalTable})


  #### TAB4 ####
  output$VolcanoHeader <- renderText({ paste("<h3> <u> <font color=\"black\"><b>","Gene effect & selection",
                                             "</b></font> </u> </h3>") })

  gene.plot <- reactive({gene.view(OC_object(),dat(),
                                   geneList = unlist(strsplit(input$GeneListRisk, split ="," )),
                                   results()$LT,results()$risk.raw)})
  output$effectPlot <- renderPlotly({gene.plot()$selectInflPlot})

  output$heatmaphead <- renderText({ paste("<h3> <u> <font color=\"black\"><b>","Distribution of features by risk",
                                           "</b></font> </u> </h3>") })

  output$binMap <- renderPlot({gene.plot()$heatmap.sorted.bin})

  output$contMap <- renderPlot({gene.plot()$heatmap.sorted.cont})


  ##### VALIDATION ####
  in.dat <- reactive({read.csv(input$newdata$datapath,row.names = 1)})

  v.results <- reactive({

    validate.app(OC_object = OC_object(),
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

    predIncomingSurv.app(OC_object = OC_object(),
                         new.data=in.dat.pred(),
                         surv.print= 1:min(10,nrow(in.dat.pred)),
                         riskRefit=results()$RiskRefit)

  })

  output$newRisk.ind <- renderPlot({ind.results()$RiskHist})
  output$KM.ind <- renderPlotly({ind.results()$IncKM})
  # output$survprob.ind <- renderTable(({ind.results()$survivalEst}))



  ##### DOWNLOAD REPORT ######
  output$report <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = "report.html",
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "report.Rmd")
      file.copy("report.Rmd", tempReport, overwrite = TRUE)

      # Set up parameters to pass to Rmd document
      params <- list("riskHist"=results()$RiskHistogram,"riskSum"=results()$RiskScoreSummary,
                     "CI"=results()$CI)
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  )


}

runApp(shinyApp(ui, server))
