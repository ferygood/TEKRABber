#' Visualize TEKRABber results with shiny app
#' @description To help user explore their results using TEKRABber, it visualize the results
#' using a self-written shiny app with two tabs, including the expression and correlation of
#' genes and TEs.
#' @usage appTEKRABber(DEresult, corrRef, corrCompare, metadata)
#' @param DEresult the output variable from using DEgeneTE()
#' @param corrRef the correlation result of your reference species using corrOthologTE()
#' @param corrCompare the correlation result of your compare species using corrOrthologTE()
#' @param metadata the same metadata you use for DEgeneTE()
#' @return an app can display differentially expressed genes/TE and the correlation results
#' @export
#' @examples 
#' data(ctInputDE)
#' geneInputDE <- ctInputDE$gene
#' teInputDE <- ctInputDE$te
#' 
#' metaExp <- data.frame(experiment = c(rep("control", 5), rep("treatment", 5)))
#' rownames(metaExp) <- colnames(geneInputDE)
#' metaExp$experiment <- factor(metaExp$experiment, levels = c("control", "treatment"))
#' 
#' resultDE <- DEgeneTE(
#'   geneTable = geneInputDE,
#'   teTable = teInputDE,
#'   metadata = metaExp,
#'   contrastVector = c("experiment", "control", "treatment"),
#'   expDesign = FALSE
#' )
#' 
#' #library(SummarizedExperiment)
#' #data(ctCorr)
#' #geneConCorrInput <- assay(ctCorr$geneCorr[, ctCorr$geneCorr$experiment=="control"])
#' #teConCorrInput <- assay(ctCorr$teCorr[, ctCorr$teCorr$experiment=="control"])
#' #geneTreatCorrInput <- assay(ctCorr$geneCorr[, ctCorr$geneCorr$experiment=="treatment"])
#' #teTreatCorrInput <- assay(ctCorr$teCorr[, ctCorr$teCorr$experiment=="treatment"])
#' 
#' #controlCorr <- corrOrthologTE(
#' #  geneInput = geneConCorrInput,
#' #  teInput = teConCorrInput,
#' #  corrMethod = "pearson",
#' #  padjMethod = "fdr",
#' #  filename = "controlCorrResult.csv"
#' #)
#' 
#' #treatmentCorr <- corrOrthologTE(
#' #  geneInput = geneTreatCorrInput,
#' #  teInput = teTreatCorrInput,
#' #  corrMethod = "pearson",
#' #  padjMethod = "fdr",
#' #  filename = "treatmentCorrResult.csv"
#' #)
#' 
#' #appTEKRABber(
#' #  DEresult = resultDE,
#' #  corrRef = controlCorr,
#' #  corrCompare = treatmentCorr,
#' #  metadata = metaExp
#' #)
#' 
appTEKRABber <- function(DEresult, corrRef, corrCompare, metadata) {
    requireNamespace("shiny")
    requireNamespace("ggplot2")
    requireNamespace("ggpubr")
    geneInput <- DEresult$normalized_gene_counts
    geneRes <- DEresult$gene_res
    teInput <- DEresult$normalized_te_counts
    teRes <- DEresult$te_res
    
    rownames(geneInput) <- gsub("[()]", "", rownames(geneInput))
    rownames(teInput) <- gsub("[()]", "", rownames(teInput))
    rownames(geneRes) <- gsub("[()]", "", rownames(geneRes))
    rownames(teRes) <- gsub("[()]", "", rownames(teRes))
    corrRef$geneName <- gsub("[()]", "", corrRef$geneName)
    corrRef$teName <- gsub("[()]", "", corrRef$teName)
    corrCompare$geneName <- gsub("[()]", "", corrCompare$geneName)
    corrCompare$teName <- gsub("[()]", "", corrCompare$teName)
    
    # for expression tab: avoid subscript out of bounds with empty expression
    geneChoices <- intersect(corrRef$geneName, rownames(geneInput))
    geneChoices <- intersect(corrCompare$geneName, geneChoices)
    teChoices <- intersect(corrRef$teName, rownames(teInput))
    teChoices <- intersect(corrCompare$teName, teChoices)
    
    ui <- fluidPage(
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
                        selectInput("geneCorr", "Gene name", choices = corrRef$geneName),
                        selectInput("teCorr", "Transposable element", choices = corrRef$teName)
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
    
    server <- function(input, output, session) {
        
        ## Expression tab
        expData <- reactive({
            req(input$gene)
            req(input$te)
            list("gene" = input$gene, "te" = input$te)
        })
        
        output$boxplotGene <- renderPlot({
            exp <- expData()
            df.gene <- data.frame(geneInput[exp$gene, ], metadata[, 1])
            colnames(df.gene) <- c("expression", "group")
            df.gene <- df.gene %>%
                mutate(expression = log2(expression))
            
            ggplot(df.gene, aes(x = group, y = expression, fill = group)) +
                geom_boxplot(outlier.shape = NA) +
                geom_jitter(width = 0.2) +
                theme(text = element_text(size = 15), legend.position = "bottom") +
                xlab("") +
                ylab("expression") +
                labs(subtitle = exp$gene)
        })
        
        output$boxplotTE <- renderPlot({
            exp <- expData()
            df.te <- data.frame(teInput[exp$te, ], metadata[, 1])
            
            colnames(df.te) <- c("expression", "group")
            
            df.te <- df.te %>%
                mutate(expression = log2(expression))
            
            ggplot(df.te, aes(x = group, y = expression, fill = group)) +
                geom_boxplot(outlier.shape = NA) +
                geom_jitter(width = 0.2) +
                theme(text = element_text(size = 15), legend.position = "bottom") +
                xlab("") +
                ylab("expression") +
                labs(subtitle = exp$te)
        })
        
        output$DETable <- renderTable(
            {
                exp <- expData()
                df.ref <- geneRes[exp$gene, ]
                df.compare <- teRes[exp$te, ]
                df <- rbind(df.ref, df.compare)
                df
            },
            rownames = TRUE
        )
        
        output$corrTable <- renderTable({
            exp <- expData()
            df.ref <- corrRef %>%
                filter(geneName == exp$gene) %>%
                filter(teName == exp$te)
            
            df.compare <- corrCompare %>%
                filter(geneName == exp$gene) %>%
                filter(teName == exp$te)
            df <- rbind(df.ref, df.compare)
            df$design <- c("reference", "comparison")
            df <- df[, c(6, 1, 2, 3, 4, 5)]
            df
        })
        
        ## Correlation
        scatterRefPlot <- ggplot(corrRef, mapping = aes(x = padj, y = coef)) +
            geom_point(data = corrRef, color = "black", size = 2) +
            theme(text = element_text(size = 15)) +
            xlab("adjusted p-value") +
            ylab("correlation coefficient") +
            labs(title = "Reference") +
            geom_vline(xintercept = 0.05, linetype="dashed", color = "blue") +
            geom_vline(xintercept = 0.01, linetype="dashed", color = "orange")
        
        scatterComparePlot <- ggplot(corrCompare, mapping = aes(x = padj, y = coef)) +
            geom_point(data = corrCompare, color = "black", size = 2) +
            theme(text = element_text(size = 15)) +
            xlab("adjusted p-value") +
            ylab("correlation coefficient") +
            labs(title = "Comparison") +
            geom_vline(xintercept = 0.05, linetype="dashed", color = "blue") +
            geom_vline(xintercept = 0.01, linetype="dashed", color = "orange")
        
        corrGene <- reactive({
            input$geneCorr
        })
        corrTE <- reactive({
            input$teCorr
        })
        
        output$scatterRef <- renderPlot({
            scatterRefPlot +
                geom_point(
                    data = corrRef %>%
                        filter(geneName == input$geneCorr & teName == input$teCorr),
                    color = "red", size = 2.5
                )
        })
        
        output$scatterCompare <- renderPlot({
            scatterComparePlot +
                geom_point(
                    data = corrCompare %>%
                        filter(geneName == input$geneCorr & teName == input$teCorr),
                    color = "red", size = 2.5
                )
        })
        
        output$pointNearRef <- renderPrint({
            nearPoints(corrRef[, c(1, 2, 4, 5)], input$scatterRefClick)
        })
        
        output$pointBrushRef <- renderPrint({
            brushedPoints(corrRef[, c(1, 2, 4, 5)], input$scatterRefBrush)
        })
        
        output$pointNearCompare <- renderPrint({
            nearPoints(corrCompare[, c(1, 2, 4, 5)], input$scatterCompareClick)
        })
        
        output$pointBrushCompare <- renderPrint({
            brushedPoints(corrCompare[, c(1, 2, 4, 5)], input$scatterCompareBrush)
        })
    }
    
    shinyApp(ui, server)
}