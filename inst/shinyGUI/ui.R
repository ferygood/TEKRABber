fluidPage(
    titlePanel(
        "TEKRABber",
    ),
    tabsetPanel(
        tabPanel(
            "Expression",
            sidebarLayout(
                sidebarPanel(
                    selectInput("gene", "Gene name", choices = geneChoices),
                    selectInput("te", "Transposable element", choices = teChoices)
                ),
                mainPanel(
                    h4("Log2 expression of selected gene and transposable element (TE):"),
                    fluidRow(
                        splitLayout(cellWidths = c("48%", "48%"), plotOutput("boxplotGene"), plotOutput("boxplotTE"))
                    ),
                    HTML("<br/>"),
                    h4("Differentially expressed selected gene and TE table:"),
                    tableOutput("DETable"),
                    HTML("<br/>"),
                    h4("Correlation between seleted gene and TE:"),
                    tableOutput("corrTable")
                )
            )
        ),
        tabPanel(
            "Correlation",
            sidebarLayout(
                sidebarPanel(
                    selectInput("geneCorr", "Gene name", choices = appRef$geneName),
                    selectInput("teCorr", "Transposable element", choices = appRef$teName)
                ),
                mainPanel(
                    h4("Distribution of Coefficients with adjusted p-value: "),
                    fluidRow(
                        splitLayout(
                            cellWidths = c("48%", "48%"),
                            plotOutput("scatterRef", click = "scatterRefClick", brush = brushOpts(id = "scatterRefBrush")),
                            plotOutput("scatterCompare", click = "scatterCompareClick", brush = brushOpts(id = "scatterCompareBrush"))
                        )
                    ),
                    HTML("<br/>"),
                    h4("Data points near to your selection:"),
                    fluidRow(
                        splitLayout(cellWidths = c("48%", "48%"), verbatimTextOutput("pointNearRef"), verbatimTextOutput("pointNearCompare"))
                    ),
                    HTML("<br/>"),
                    h4("Data points in your selected area:"),
                    fluidRow(
                        splitLayout(cellWidths = c("48%", "48%"), verbatimTextOutput("pointBrushRef"), verbatimTextOutput("pointBrushCompare"))
                    )
                )
            )
        )
    )
)