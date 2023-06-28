library(shiny)
library(plotly)
library(gridlayout)
library(bslib)
library(ggplot2)
library(ggpubr)

runMyApp <- function(corrRef, corrCompare=NULL, DEobject) {
    
    # metadata subset reference and comparison
    metadata <- data.frame(DEobject$gene_dds[[1]])
    unique_value <- unique(metadata[,1])
    ref_indices <- which(metadata[,1] == unique_value[1])
    compare_indices <- which(metadata[,1] == unique_value[2])
    
    # normalized expression
    geneRef <- data.frame(DEobject$normalized_gene_counts)[, ref_indices]
    geneCompare <- data.frame(DEobject$normalized_gene_counts)[, compare_indices]
    teRef <- data.frame(DEobject$normalized_te_counts)[, ref_indices]
    teCompare <- data.frame(DEobject$normalized_te_counts)[, compare_indices]
    
    # log2FC information
    geneLFC <- data.frame(DEobject$gene_res)
    teLFC <- data.frame(DEobject$te_res)
    
    ui <- grid_page(
        layout = c(
            "header header header",
            "area1  area2  area3",
            "area1  area4  area5",
            ".      area6  area7"
        ),
        row_sizes = c(
            "100px",
            "0.79fr",
            "0.79fr",
            "0.79fr"
        ),
        col_sizes = c(
            "250px",
            "0.79fr",
            "0.79fr"
        ),
        gap_size = "1rem",
        grid_card_text(
            area = "header",
            content = "TEKRABber",
            alignment = "start",
            is_title = TRUE
        ),
        grid_card(
            area = "area1",
            card_header("parameters"),
            card_body_fill(
                textInput(
                    inputId = "geneIdInput",
                    label = "Ensembl Gene ID",
                    value = "ENSG00000162571"
                ),
                textInput(
                    inputId = "teInput",
                    label = "Transposable Elements",
                    value = "AluJb"
                ),
                selectInput(
                    inputId = "mySelectInput",
                    label = "Log2FC",
                    choices = list("1.2" = 1.2, "1.5" = 1.5, "2.0" = 2.0)
                ),
                selectInput(
                    inputId = "mySelectInput",
                    label = "adjusted p-value",
                    choices = list("0.05" = 0.05, "0.01" = 0.01)
                ),
                actionButton(inputId = "myButton", label = "Submit ")
            )
        ),
        grid_card(
            area = "area2",
            card_body_fill(plotlyOutput(outputId = "scatterPlot"))
        ),
        grid_card(
            area = "area4",
            card_body_fill(plotOutput(outputId = "coefPlotRef"))
        ),
        grid_card(
            area = "area5",
            card_body_fill(plotOutput(outputId = "coefPlotCompare"))
        ),
        
        grid_card(
            area = "area6",
            card_body_fill(plotOutput(outputId = "deGenePlot"))
        ),
        
        grid_card(
            area = "area7",
            card_body_fill(plotOutput(outputId = "deTEPlot"))
        )
        
        
    )
    
    server <- function(input, output) {
        
        output$scatterPlot <- renderPlotly({
            scatter_data <- plot_ly(corrRef, x = ~pvalue, y = ~coef, text = ~geneName, type = "scatter", mode = "markers")
            scatter_data <- scatter_data %>% 
                layout(shapes=list(
                    list(type="line", y0 = -1, y1 = 1, x0 = 0.05, x1 = 0.05, line=list(color="red"))))
        })
        
        observeEvent(input$myButton, {
            ensemblID <- input$geneIdInput
            teName <- input$teInput
            
            geneRef_select <- data.frame(t(geneRef[ensemblID, ]))
            geneCompare_select <- data.frame(t(geneCompare[ensemblID, ]))
            
            teRef_select <- data.frame(t(teRef[teName, ]))
            teCompare_select <- data.frame(t(teCompare[teName, ]))
            
            df_ref <- cbind(geneRef_select, teRef_select)
            colnames(df_ref) <- c("gene", "TE")
            
            df_compare <- cbind(geneCompare_select, teCompare_select)
            colnames(df_compare) <- c("gene", "TE")
            
            df_ref$group <- "reference"
            df_compare$group <- "comparison"
            df_all <- rbind(df_ref, df_compare)
            
            # render ref coef plot
            output$coefPlotRef <- renderPlot({
                
                coef <- corrRef[corrRef$geneName==ensemblID & corrRef$teName==teName, "coef"]
                coef <- round(coef, 4)
                pvalue <- corrRef[corrRef$geneName==ensemblID & corrRef$teName==teName, "pvalue"]
                pvalue <- round(pvalue, 4)
                
                ggplot(df_ref, aes(x=gene, y=TE)) +
                    geom_point(colour="black", shape=21, size=3, fill="#45A9EC") +
                    labs(x=ensemblID, y=teName) +
                    theme_bw() +
                    ggtitle("Reference Correlation") +
                    annotate(
                        "text", x=max(df_ref$gene)-3, y=max(df_ref$TE)-2, 
                        label = paste0("Coefficient: ", coef, "\npvalue: ", pvalue)
                    )
            })
            
            # render compare coef plot
            output$coefPlotCompare <- renderPlot({
                
                pcor <- cor(df_compare$gene, df_compare$TE, method="pearson")
                
                ggplot(df_ref, aes(x=gene, y=TE)) +
                    geom_point(colour="black", shape=21, size=3, fill="#8BE748") +
                    labs(x=ensemblID, y=teName) +
                    theme_bw() +
                    ggtitle("Comparison Correlation") +
                    annotate(
                        "text",
                        x = max(df_compare$gene), y = max(df_compare$TE),
                        label = paste("Coefficient", round(pcor,2)),
                        hjust=1, vjust=1
                    )
            })
            
            # render expression plot
            output$deGenePlot <- renderPlot({
                
                lfc2 <- round(geneLFC[ensemblID, "log2FoldChange"], 4)
                pvalue <- round(geneLFC[ensemblID, "pvalue"], 4)
                
                plot <- ggviolin(df_all, x="group", y="gene", fill="group",
                         palette=c("#45A9EC", "#8BE748"), add="boxplot", 
                         add.params=list(fill="white")) +
                    ylab(ensemblID)
                
                plot <- plot +
                    geom_text(
                        aes(x = Inf, y = Inf, label = paste0("lfc2: ", lfc2, "\npvalue: ", pvalue)),
                        hjust = 1, vjust = 1
                    )
                
                plot
            })
            
            # render TE plot
            output$deTEPlot <- renderPlot({
                
                lfc2 <- round(teLFC[teName, "log2FoldChange"], 4)
                pvalue <- round(teLFC[teName, "pvalue"], 4)
                
                plot <- ggviolin(df_all, x="group", y="TE", fill="group",
                         palette=c("#45A9EC", "#8BE748"), add="boxplot", 
                         add.params=list(fill="white")) +
                    ylab(teName)
                
                # Add annotation
                plot <- plot +
                    geom_text(
                        aes(x = Inf, y = Inf, label = paste0("lfc2: ", lfc2, "\npvalue: ", pvalue)),
                        hjust = 1, vjust = 1
                    )
                
                plot
                
            })
            
        })
        
    }
    
    shinyApp(ui, server)
}

# Example usage:
testDEobject <- testDE



runMyApp(hmCorrResult, DEobject=testDEobject)
