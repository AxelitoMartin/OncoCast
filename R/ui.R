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
                radioButtons("method", "Select a method:",
                             c("LASSO" = "LASSO",
                               "ENET" = "ENET",
                               "RF" = "RF",
                               "GBM" = "GBM")),

                sliderInput("runs", label = "Select number of runs", min = 25,
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
                )
                # ,
                # htmlOutput("dataLoad")
                #,
                # submitButton("Run OncoCast")
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
