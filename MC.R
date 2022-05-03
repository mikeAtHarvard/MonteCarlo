library(readr)
library(tidyverse)
library(ggplot2)
library(mc2d)
library(boot)
library(scales)
library(magrittr)
library(shiny)

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#Shiny app specific code below

library(shiny)

# Define UI for dataset viewer app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Shiny App for Sensitivity Analysis"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      # Input: Select a file ----
      width = 3,
      fileInput("inputfile", "Choose CSV File",
                multiple = TRUE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      actionButton("startButton", "Start")

    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      width = 9,
      tabsetPanel(type = "tabs",
                  tabPanel("Introduction", ""),
                  tabPanel("Table", tableOutput(outputId = "TableView")),
                  tabPanel("Summary", plotOutput(outputId = "SummaryPlot")),
                  tabPanel("Overview", plotOutput(outputId = "MCPlot")),
                  tabPanel("Tornado", plotOutput(outputId = "TornadoPlot")),
                  tabPanel("Inputs", plotOutput(outputId = "InputsPlot")),
                  tabPanel("Outputs", plotOutput(outputId = "OutputsPlot")),
                  tabPanel("Sensitivity", plotOutput(outputId = "SensitivityPlot"))
      )
    )
  )
)

# Define server logic to summarize and view selected dataset ----
server <- function(input, output) {

    observeEvent(input$startButton, {
      #Using <<- assigns the variable to be global rather than <- which makes a variable local
      dm <<- as.matrix(read.csv(input$inputfile$datapath, sep=","), rownames.force = TRUE)
      mclbl <<- list()
      for (r in 1:nrow(dm)) {
        mclbl[r] <- dm[r,1]
        for (c in 1:ncol(dm)) {
        }
      }
    
      #Create the Monte Carlo model
      
      mcstring <<- list()
      mctext <<- "mc("
      outputstr <<- "output <- "
      
      modMC <- mcmodel({
        
        set.seed(NULL)
        
        for (r in 1:nrow(dm)) {
          if (dm[r,8] ==  "norm") { mcstring[r] = paste(rawToChar(as.raw(97+r-1)), " <- mcstoc(r", dm[r,8], ", mean=as.numeric(dm[", r, ",6]), sd=as.numeric(dm[", r, ",7])", sep="") }
          if (dm[r,8] ==  "pert") { mcstring[r] = paste(rawToChar(as.raw(97+r-1)), " <- mcstoc(r", dm[r,8], ", min=as.numeric(dm[", r, ",2]), mode=as.numeric(dm[", r, ",4]), max=as.numeric(dm[", r, ",5]), shape=4", sep="") }
          if (dm[r,9] ==  "yes") { mcstring[r] = paste(mcstring[r], ", rtrunc=TRUE, linf=as.numeric(dm[", r, ",10]), lsup=as.numeric(dm[", r, ",11])", sep="") }
          mcstring[r] = paste(mcstring[r], ")", sep="")
          eval(parse(text=toString(mcstring[r])))
          if (r!= nrow(dm)) {
            mctext = paste(mctext,  rawToChar(as.raw(97+r-1)), ", ", sep="")
          } else {
            mctext = paste(mctext, rawToChar(as.raw(97+r-1)), ", output)", sep="")
            outputstr = paste(outputstr, dm[r,12], sep="")
          }
        }
        
        eval(parse(text=toString(outputstr)))
        eval(parse(text=toString(mctext)))
        
      })
      
      evalmcmod(modMC, nsv=1000, nsu=1000)
      m <<- evalmcmod(modMC, nsv=1000, nsu=1000)
      
      #These lines would be used in a standalone R file
      #plot(m)
      #plot(tornado(m, (nrow(dm) + 1), 'complete.obs', 'spearman', c(0.025, 0.975)))
      #plot(tornado(m, 7, "complete.obs", "pearson", c(0.025, 0.975)))
      #plot(tornado(m, 7, "complete.obs", "kendall", c(0.025, 0.975)))
      
      #Un-Monte Carlo the model and convert to a dataframe
      unmcModMC <- unmc(m)
      dframe <- data.frame(unmcModMC)
      
      #To Output the results into an external file
      #write.table(dframe, "Output.csv", sep = ",", row.names = FALSE)
      
      #Plot the output
      r1 <<- ggplot(dframe, aes(x=output)) + stat_ecdf(pad=FALSE) +
        labs(c="Output", y="Cumulative Frequency")
      r2 <<- ggplot(dframe, aes(x=output)) + geom_histogram(aes(y=..density..), alpha=0.5, color="black", fill="white") +
        geom_vline(aes(xintercept=mean(output))) +
        geom_text(aes(x=5, y=0.4, label=round(mean(output), digits=2))) +
        labs(x="Output", y="Density")
      #This line would be used in a standalone R file
      #multiplot(r1, r2, cols=2)
      
      #Histogram plots for the inputs
      p <- list()
      pstr <<- "multiplot("
      for (r in 1:nrow(dm)) {
        p[r] = paste("p", r, " <<- ", "ggplot(dframe, aes(x=", rawToChar(as.raw(97+r-1)), ")) + geom_histogram(aes(y=..density..), alpha=0.5, color='black', fill='white') + geom_density(alpha=0.2) +labs(x='", mclbl[r], "', y='Density') + scale_y_continuous(labels=NULL)", sep="")
        if (r!= nrow(dm)) {
          pstr = paste(pstr, "p", r, ", ", sep="")
        } else {
          pstr = paste(pstr,"p", r, ", cols=", nrow(dm) %/% 3, ")", sep="")
          glblpstr <<- pstr
        }
        eval(parse(text=toString(p[r])))
      }
      for (r in 1:nrow(dm)) {
        eval(parse(text=toString(p[r])))
      }
      #This line would be used in a standalone R file
      #eval(parse(text=toString(pstr)))
      
      #Scatter plots for the outputs
      q <- list()
      qstr <<- "multiplot("
      for (r in 1:nrow(dm)) {
        q[r] = paste("q", r, " <<- ", "ggplot(dframe, aes(x=", rawToChar(as.raw(97+r-1)), ", y=output)) + geom_point() + geom_density_2d() +labs(x='", mclbl[r], "', y='Output')", sep="")
        ifelse ((r != nrow(dm)), (qstr = paste(qstr, "q", r, ", ", sep="")), (qstr = paste(qstr,"q", r, ", cols=", nrow(dm) %/% 3, ")", sep="")))
        glblqstr <<- qstr
      }
      for (r in 1:nrow(dm)) {
        eval(parse(text=toString(q[r])))
      }
      #This line would be used in a stand alone R file
      #eval(parse(text=toString(qstr)))
      
    }, once = TRUE)
  
    tmpdm <- reactive ({
    req(input$inputfile, file.exists(input$inputfile$datapath))
    read.csv(input$inputfile$datapath, sep=",")
    })
      
    output$TableView <- renderTable({
      req(input$inputfile, file.exists(input$inputfile$datapath))
      dm <<- read.csv(input$inputfile$datapath, sep=",")
      return(dm)
    })
    
    output$SummaryPlot <- renderPlot({
      multiplot(r1, r2, cols=2)
    })
    
    output$MCPlot <- renderPlot({
      plot(m)
    })
    
    output$TornadoPlot <- renderPlot({
      plot(tornado(m, (nrow(dm) + 1), 'complete.obs', 'spearman', c(0.025, 0.975)))
    })
    
    output$InputsPlot <- renderPlot({
      eval(parse(text=toString(glblpstr)))
    })
    
    output$OutputsPlot <- renderPlot({
      eval(parse(text=toString(glblqstr)))
    })
    
    output$SensitivityPlot <- renderPlot({
      multiplot(r1, r2, cols=2)
    })
    
}
    
# Create Shiny app ----
shinyApp(ui, server)
