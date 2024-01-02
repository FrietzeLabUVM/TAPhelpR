# library(shiny)
# library(DT)
# library(RColorBrewer)
# library(colourpicker)
# library(magrittr)
# source("setup.R")


#' @import DT RColorBrewer colourpicker
#' @rawNamespace import(shiny, except = c(runExample, renderDataTable, dataTableOutput))
#'
ucsc_track_ui = function(){
  fluidPage(
    tags$head(
      tags$link(
        rel = "stylesheet",
        type = "text/css",
        href = "my.css")
    ),
    conditionalPanel(
      condition = "false == true",
      colourInput("col", "Select colour", "purple")
    ),
    tabsetPanel(id = "tabsetPanel",
                tabPanel("File Selection", tags$p(),
                         sidebarLayout(
                           sidebarPanel(width = 3,
                                        actionButton("refreshFiles", "Start Over"),
                                        tags$hr(),
                                        actionButton("removeSelected", "Remove Selected"),
                                        actionButton("limitToSelected", "Limit To Selection"),
                                        actionButton("limitToVisible", "Limit To Filtered"),
                                        tags$hr(),
                                        radioButtons("fileType", label = "File Type", choices = c("bigWig", "bigBed"), selected = "bigWig"),
                                        tags$hr(),
                                        uiOutput("colorBy"),
                                        selectInput("colorStyle", label = "Color Style", selected = "Set1", choices = rownames(RColorBrewer::brewer.pal.info)),
                                        checkboxInput("freeColor", "Free Color"),
                                        uiOutput("colorAssign")
                           ),
                           mainPanel(width = 9,

                                     DT::dataTableOutput(outputId = "bwTable")
                           )
                         )
                ),
                tabPanel("Track Configuration", br(),
                         sidebarLayout(
                           sidebarPanel(width = 3,
                                        sliderInput("trackSize", label = "Track Size", min = 16, max = 128, value = 32),
                                        numericInput("viewLimitsMin", label = "View Limits Min", value = 0),
                                        numericInput("viewLimitsMax", label = "View Limits Max", value = 50),
                                        checkboxInput("autoScale", "Auto Scale"),
                                        checkboxInput("alwaysZero", "Auto Include Zero"),
                                        selectInput("visibility", label = "Visibility",
                                                    choices = names(UCSC_VIS), selected = "full"),
                                        selectInput("windowFun", "Windowing Function", choices = UCSC_WIN, selected = "mean"),
                                        sliderInput("smoothWin", "Smoothing Window", value = 0, min = 0, max = 16, step = 1),
                                        checkboxInput("showZero", "Line At Zero", value = TRUE),
                                        checkboxInput("showYref", "Include Custom Y-line", value = FALSE),
                                        numericInput("numYref", "Custom Y-line", value = 0, step = .5)

                           ),
                           mainPanel(width = 9,
                                     tags$h3("Use configuration by link, download file, or copy-paste."),
                                     tags$hr(),
                                     tags$h4("Links: (temporary)"),
                                     uiOutput("ucscLink"),
                                     uiOutput("cfgLink"),
                                     tags$hr(),
                                     tags$h4("Download:"),
                                     downloadButton(outputId = "dlTracks", label = "Download Tracks"),
                                     tags$hr(),
                                     tags$h4("Configuration:"),
                                     verbatimTextOutput("cfgText")
                                     # DT::dataTableOutput(outputId = "cfgTable")
                           )
                         )
                )
    )
  )
}
