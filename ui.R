library(shiny)
library(shinyBS)
library(shinydashboard)
library(shinyjs)
library(shinycssloaders)
library(plotly)

source("addDeps.R")
source("ui-about.R")
source("modals.R")

commonStyles <- list(
  includeCSS('www/ecsu.css'),
  HTML("<link href='https://fonts.googleapis.com/css?family=Courgette' rel='stylesheet' type='text/css'>"),
  HTML("<link rel='stylesheet' href='https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css' integrity='sha384-BVYiiSIFeK1dGmJRAkycuHAHRg32OmUcww7on3RYdg4Va+PmSTsz/K68vbdEjh4u' crossorigin='anonymous'>"),
  tags$style(HTML("
                               
                  /* color table selections */
                  table.dataTable tr.selected td, table.dataTable td.selected {
                  background-color: maroon !important;
                  color: white
                  }
                  
              .blue-button {
              #  position:relative;
              #  bottom: 12px;
              #  height: 35px;
                #width: 66px;
                #color: white;
                background-image: linear-gradient(#04519b,#044687 60%,#033769);
                color:#FFF;
              }

                  /* 'hide' tab */
                  #sidebar, .navbar {
                  background-image: linear-gradient(#04519b,#044687 60%,#033769);
                  color:#FFF;
                  }
                  
                  td, td {
                    font-size:80%;
                  }
 
                  .navbar-header > .navbar-brand {
                  font-family:Courgette;
                  color:#FFF;
                  }
                  
                  .navbar-nav > li > a[data-value='hide'] {color: maroon !important;}
                  
                  /* Top level tabs */
                  .navbar-default .navbar-nav > .active > a, .navbar-default .navbar-nav > .active > a:focus, .navbar-default .navbar-nav > .active > a:hover {
                  color: white;
                  background-color: maroon
                  }
                  
                  
                  
                  .navbar-default .navbar-nav > li > a:hover {
                  color: darkred;
                  }
                  
                  
                  .navbar-default .navbar-nav > li > a {
                  color: white;
                  }
                  
                  
                  /* 2nd level (Disease menu) tabs */
                  .navbar-default .navbar-nav > li > ul > li[class='active'] > a {
                  color: white;
                  background-color: maroon
                  }
                  
                  .navbar-default .navbar-nav > li > ul > li > a {
                  color: darkblue;
                  }
                  
                  
                  .navbar-default .navbar-nav > .active > li > ul > li > a, .navbar-default .navbar-nav > .active > li > ul > li > a:focus, .navbar-default .navbar-nav > .active > li > ul > li > a:hover {
                  color: white;
                  background-color: maroon
                  }
                  
                  .navbar-default .navbar-nav > li >ul > li > a:hover {
                  color: darkred;
                  }
                  
                #summaryRow {
                  position:sticky;
                  position:-webkit-sticky;
                  top:0;
                  background:white;
                  z-index:1000;
                  }                  

                  "))  
  
)

# common header for all pages
commonHeader <- list(
  
  welcomeModal, 
  filterModal,
  
  # summary row
  fluidRow(id = "summaryRow",
    #shiny::column(width = 2,
    #              radioButtons("rbDiseaseLimits", "Article Limits:",
    #                           c("Cancer-related only" = "cancer", "All articles" = "none"), selected = "cancer")),

    shiny::column(width=2,
                  actionButton("btnNewSearch", "New Gene Search", class = "blue-button")
    ),
    shiny::column(width =10 ,
        htmlOutput("summaryHeader")              
    ),
    br()

  ),
  
  
  hr(style = "padding: 0px; margin: 0px"),
 
  fluidRow(shiny::column(width = 12,
                         conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                          HTML("<div class=\"progress\" style=\"position: fixed;  width: 100%; height:25px; !important\">
                                               <div class=\"progress-bar progress-bar-striped active\" role=\"progressbar\" aria-valuenow=\"100\" aria-valuemin=\"0\" aria-valuemax=\"100\" style=\"opacity: 1; width:100%\" !important>
                                               <span id=\"bar-text\"><b>Loading, please wait...</b></span>
                                               </div></div>")
                                          ))
           
                         ),
  br(), commonStyles
)


# add a tab panel with given title, table, and graph,
# may want to consider setting click = "CancerGraph_click", etc on graph
addTabPanel <- function(title, tableId, graphId = NULL) {
  
  tabPanel(title, 
           fluidRow(div(style = "height: 300px",
             shiny::column(width = 5,
                           withSpinner(DT::dataTableOutput(tableId), type = 3,
                                       color.background = "white")
             ),
             if (!is.null(graphId)) {
               shiny::column(width = 7,
                             if (graphId == "cancerGraph") {
                               withSpinner(plotOutput(graphId), type = 3, 
                                           color.background = "white")  
                             } else {
                               withSpinner(plotlyOutput(graphId), type = 3, 
                                           color.background = "white")  
                             }
                             
                             
               )
             }
           ))
  )
}

articlesPanel <- function() {
  tabPanel('Articles',         
           useShinyjs(),  
           br(),
           fluidRow(
             shiny::column(width = 3,
                           bsButton(inputId = "btnPubTator", label = "Load/Refresh PubTator Results", style = "info")
             ),
             shiny::column(width = 3,
                           div(align = "right",
                               bsButton("btnPubTatorGo", label = "View Results in PubTator", style = "danger")
                           )
             )
           ),
           
           fluidRow(
             shiny::column(id = "colPubs", width = 2, 
                            DT::dataTableOutput("articleTable")),
             shiny::column(id = "colPubs", width = 10,
                                uiOutput("articles")
             ),
             div(id = "pageBottom", style ="visibility:hidden",
                 a(id = "pageDownLink", href = "#pageBottom", "click")
             )
           )
  )
}


addDownloadsTabPanel <- function(title) {
  tabPanel(title,
           #    fluidRow(
           #      downloadButton("downloadCancerTypesData", "Download Cancer Types Table")
           #    ),
           
           
           sidebarLayout(
             sidebarPanel(
               
               fluidRow(
                 downloadButton("downloadCancerTypesData", "Download Cancer Types Table")
               ),
               
               fluidRow(
                 downloadButton("downloadDrugTreatmentsData", "Download Drugs Table")
               ),
               
               fluidRow(
                 downloadButton("downloadMutationsData", "Download Mutations Table")
               ),
               
               fluidRow(
                 downloadButton("downloadGenesData", "Download Genes Table")
               )
             ),
             mainPanel(
               h3("To download one of the tables created from your inquiry, click on the corresponding button."),
               h3("Ex: Click on the \"Download Mutations Table\" button to acquire a .csv file of your gene mutations search.")
             )
           )
  )
}



logPanel <- function() {
  tabPanel("Log", verbatimTextOutput("log"))
}




shinyUI(


  navbarPage(title = 'Cancer Publication Portal',
             id = "headerNavBarPage", 
    tabPanel("Home",
          commonHeader, 
          fluidRow(column(style='border-right: 1px solid',width = 12,
          tabsetPanel(id = "MainPage",
            #addTabPanel('Cancer Types', "cancerSummaryTable"),
            #addTabPanel('Treatments', "paResults"),
            addTabPanel("Cancer Types", "diseaseResults", "cancerGraph"),
            addTabPanel('Drugs', 'chemResults', "chemGraph"),
            addTabPanel("Mutations", "mutationResults", "mutGraph"),
            addTabPanel("Genes", "geneResults"),
            addDownloadsTabPanel("Download"),
            tabPanel("Articles", articlesPanel())

          )))#, # end tabsetPanel and 1st row
          #fluidRow(column(width = 12, articlesPanel() ) )# end second row (articles panel)
          ), # end Portal Panel
      tabAbout,
      logPanel()
  ) # end navbarPage
) # end shinyUI
               

