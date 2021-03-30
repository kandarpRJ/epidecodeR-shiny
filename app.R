options(shiny.maxRequestSize = 100*1024^2)
library(shiny)
library(shinybusy)
library(DT)
library(shinyBS)

options(repos = BiocManager::repositories())

source("./epidecodeR.R")

ui <- fluidPage(
  add_busy_bar(color = "red", height = "8px"),
  titlePanel("epidecodeR | an online functional exploration tool for epigenetic and epitranscriptomic regulation"),
  sidebarLayout(
    sidebarPanel(
      fluidRow(
        h3("Download demo files here: ", downloadButton("demo1files", "demo1 files"), downloadButton("demo2files", "demo2 files"))
      ),
      fluidRow(
        fileInput("events", h3("Events input*")),
      ),
      # conditionalPanel(
      #   condition = "output.dataformat1 == 'bed'",
      #   
      #   fluidRow(
      #     selectInput("mol", h3("Molecule type"),
      #                 choices = list("DNA" = "dna", "RNA" = "rna"), selected = "rna")
      #   ),
      #   fluidRow(
      #     fileInput("gtf_file", h3("GTF input"))
      #   ),
      #   fluidRow(
      #     textInput("id_type", h3("ID type to select from GTF"), 
      #               value = "")
      #   ),
      #   conditionalPanel(
      #     condition = "input.mol == 'dna' && output.dataformat2 == 'gtf'",
      #     
      #     fluidRow(
      #       sliderInput("boundaries", h3("Genomic boundaries of genes"),
      #                   min = 0, max = 5000, value = 0)
      #     )
      #   ),
      # ),
      
      fluidRow(
        fileInput("deg", h3("DEG list file*"))
      ),
      fluidRow(
        textInput("pval", h3("P value"), 
                  value = "0.05")
      ),
      fluidRow(
        selectInput("param", h3("Select group"), 
                    choices = list("Two groups [0:1+]" = 1, "Three groups (0:[1-N]:N+)" = 2,
                                   "Four groups (0:1:[2-N]:N+)" = 3), selected = 1)
      ),
      conditionalPanel(
        condition = "input.param == 2",
        
        fluidRow(
          sliderInput("ints2", h3("Interval range for variable group [1-N]"),
                      min = 1, max = 8, value = 4)
        )
      ),
      conditionalPanel(
        condition = "input.param == 3",
        
        fluidRow(
          sliderInput("ints3", h3("Interval range for variable group [2-N]"),
                      dragRange = FALSE, min = 2, max = 8, value = 5)
        )
      ),
      fluidRow(
        selectInput("type", h3("Plot type"),
                    choices = list("Theoretical" = "t", "Emperical" = "e", "Both" = "both"), selected = "both")
      ),
      fluidRow(
        sliderInput("range", h3("X-axis range"),
                    min = -20, max = 20, value = c(-2,2))
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Renders", fluidPage(textOutput("message"),
                                      tags$head(tags$style("#message{color: red;
                                                                   font-size: 20px;
                                                                   font-style: italic;
                                                                  }"
                                                                        )
                                      ),
                                      tableOutput("paramvals"),
                                      fluidRow(
                                        column(6, uiOutput("dbutton1"), 
                                               plotOutput("outgraph")),
                                        # conditionalPanel(
                                        #   condition = "input.type == 'e' || input.type=='both'",
                                        #   
                                        #   column(6, uiOutput("dbutton2"), 
                                        #          plotOutput("plottestgraph"))),
                                        column(6, uiOutput("dbutton2"), 
                                               plotOutput("plottestgraph")),
                                        uiOutput("tablelist"),
                                        dataTableOutput("outputtable")
                                      )
                                    )
                 ),
        tabPanel("Reference manual", fluidPage(mainPanel(
          htmlOutput("inc")
        )))
      )
    )
  )
)

server <- function(input, output, session) {
  getPage<-function() {
    return(includeHTML("shinyepidecodeR.html"))
  }
  output$inc<-renderUI({getPage()})
  #efiledf<-read.table("eventcounts.txt", header = TRUE, row.names = NULL, stringsAsFactors = F, sep = "\t", fill = TRUE)
  #degfiledf<-read.table("deg.txt", header = TRUE, row.names = NULL, stringsAsFactors = F, sep = "\t", fill = TRUE)
  
  output$demo1files<-downloadHandler(
    filename = 'demo1.zip',
    content <- function(file) {
      file.copy("demo1.zip", file)
    },
    contentType = "application/zip")
  
  output$demo2files<-downloadHandler(
    filename = 'demo2.zip',
    content = function(file) {
      file.copy("demo2.zip", file)
    },
    contentType = "application/zip")
  
  obj<-reactive({
    eventsfile<-input$events
    if (! isEmpty(eventsfile) && tolower(tools::file_ext(eventsfile$datapath)) == "bed") { 
      
      output$dataformat1 <- renderText("bed")
      outputOptions(output, "dataformat1", suspendWhenHidden = FALSE)
      
    } else {output$dataformat1<-renderText(tolower(tools::file_ext(eventsfile$datapath)))}
    
    if (!is.null(input$gtf_file)) {
      output$dataformat2 <- renderText("gtf")
      outputOptions(output, "dataformat2", suspendWhenHidden = FALSE)
      file.rename(input$gtf_file$datapath, gsub(".gz", ".gtf.gz", input$gtf_file$datapath))
      gtffile<-gsub(".gz", ".gtf.gz", input$gtf_file$datapath)
    } else {
      gtffile=input$gtf_file
    }
    output$message<-NULL
    degfile<-input$deg
    idtype<-input$id_type
    req(input$events)
    req(input$deg)
    output$interval<-renderText(input$ints)
    inputints=NULL
    if (input$param==2) {
      inputints<-c(1, input$ints2)
    }
    if (input$param==3) {
      inputints<-c(2, input$ints3)
    }
    tc<-tryCatch({anobj<-epidecodeR(events = eventsfile$datapath, gtf_file = gtffile, boundaries = input$boundaries,
                  id_type = input$id_type, deg = degfile$datapath, pval = as.numeric(input$pval), param = input$param, 
                  ints = inputints)},
             error=function(cond) {
               output$message<-renderText("Plot rendering did not complete as expected!\nPlease check inputs! [Hint: Are the id types matching between events and DEG list?]")
             })
  })
  output$class<-renderText(class(obj())[1])
  output$tablelist<-renderUI({
    req(obj())
    if (class(obj())[1]=="epidecodeR") {
      l<-names(obj()@grptables)
      selectInput("ntab", h5("Gene list for selected group"), choices = l, selected = as.character(l[1]))
    }
  })
  
  output$outputtable<-renderDataTable(server = FALSE,{
    if (class(obj())[1] == "epidecodeR" & !is.null(input$ntab)) {
      tablist<-(obj()@grptables[[input$ntab]])
      rownames(tablist)<-NULL
      tablist
    }},
  extensions = c('Buttons'), 
  options = list(
    dom = 'Bfrtip',
    buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
  ))
  
  output$dbutton1<-renderUI({
    req(obj())
    if (class(obj())[1]=="epidecodeR") {
      downloadButton('downb1')
    }
  })
  
  output$dbutton2<-renderUI({
    req(obj())
    if (class(obj())[1]=="epidecodeR") {
      downloadButton('downb2')
    }
  })
  
  output$outgraph<-renderPlot({
    req(obj())
    if (class(obj())[1]=="epidecodeR") {
      tryCatch({makeplot(anobj = obj(), type = input$type, lim = input$range)})      
    }
  },
  res = 100
  #height = function() {session$clientData$output_outgraph_width}
  )
  
  output$plottestgraph<-renderPlot({
    req(obj())
    if (class(obj())[1]=="epidecodeR") {
      plot.test(obj())      
    }

  },
  res = 100
  #height = function() {session$clientData$output_outgraph_width}
  )
  
  output$downb1 = downloadHandler(
    req(obj()),
    filename = 'cdfplot.png',
    content = function(file) {
      ggsave(file, plot = makeplot(anobj = obj(), type = input$type, lim = input$range), device = "png", dpi = "print")
    })
  
  output$downb2 = downloadHandler(
    req(obj()),
    filename = 'plottest.png',
    content = function(file) {
      ggsave(file, plot = plot.test(obj()), device = "png", dpi = "print")
    })
  
  session$onSessionEnded(function() {
    cat("Session Ended\n")
    unlink(eventsfile)
    unlink(gtffile)
    unlink(degfile)
  })
}



shinyApp(ui = ui, server = server)
