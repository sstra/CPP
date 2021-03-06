welcomeModal <-  bsModal("welcomeModal",HTML("<i>Cancer Publication Portal</i>"), trigger = "btnNewSearch", size = "large",
      p(strong("Instructions:"), "This is a beta version of a Cancer Publication Portal for summarizing and searching cancer-related literature.",
        "To start, select a gene and click Search to summarize cancer publications for that gene. When the results are displayed, you can click on a row in any 
        table to further filter the results. Filters can be removed by removing the term from the 
        drop down box."),
      p("Feedback is welcome, and can be sent to dancikg@easternct.edu"),
      
      hr(class = "blue-button", style="height: 2px"),
      
      fluidRow(
        column(2, style="padding-right:0px",
          selectizeInput("geneInput", label = "Select a gene", choices = NULL)
        ),
        column(2,style = "vertical-align:middle; padding-left:0px",
           HTML("<label class = 'control-label' style='visibility:hidden'>Hello</label>"),
           div(
                actionButton("btnGeneSearch", "Summarize Cancer Articles", class = "blue-button")
           )
        )
      )
        
  )


filterModal <- bsModal("filterModal", "Remove filters", "btnRemoveFilters",
    fluidRow(column(12,
     HTML("<p>To remove a filter simply delete the term from the dropdown menus below. Changes take effect immediately.</p><br>")
    )),
    fluidRow(
      shiny::column(width=3,
                  selectInput("filterDisease", "Disease Filters", choices = NULL, multiple = TRUE, selectize = TRUE)
      ),
      shiny::column(width=3,
                  selectInput("filterChem", "Chem Filters", choices = NULL, multiple = TRUE, selectize = TRUE)
      ),
      shiny::column(width=3,
                  selectInput("filterMutations", "Mutation Filters", choices = NULL, multiple = TRUE, selectize = TRUE)
      ),    
      shiny::column(width=3,
                  selectInput("filterGenes", "Additional Gene Filters", choices = NULL, multiple = TRUE, selectize = TRUE)
      )
    ), size = "large"
)

