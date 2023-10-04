#' appTEKRABber
#'
#' @param corrRef correlation results for reference using corrOrtholgScale()
#' @param corrCompare correlation results for comparison using 
#' corrOrthologScale()
#' @param DEobject DE object using DEgeneTE() 
#'
#' @return provide an interactive shinyapp
#' @importFrom plotly renderPlotly plot_ly
#' @importFrom ggpubr ggviolin
#' @importFrom gridlayout grid_page grid_card grid_card_text
#' @import bslib
#' @import shiny
#' @import ggplot2
#' @export
#'
#' @examples
#' 
#' 
#' 
appTEKRABber <- function(corrRef, corrCompare, DEobject) {
    
    # metadata subset reference and comparison
    metadata <- data.frame(DEobject$gene_dds[[1]])
    unique_value <- unique(metadata[,1])
    ref_indices <- which(metadata[,1] == unique_value[1])
    compare_indices <- which(metadata[,1] == unique_value[2])
    
    # normalized expression
    geneRef <- log(data.frame(DEobject$normalized_gene_counts)[, ref_indices] + 1)
    geneCompare <- log(data.frame(DEobject$normalized_gene_counts)[, compare_indices] + 1)
    teRef <- log(data.frame(DEobject$normalized_te_counts)[, ref_indices] + 1)
    teCompare <- log(data.frame(DEobject$normalized_te_counts)[, compare_indices] + 1)
    
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
            card_body(
                selectizeInput(
                    inputId = "geneIdInput",
                    label = "Gene Name",
                    choices = unique(corrRef$geneName)
                ),
                selectizeInput(
                    inputId = "teInput",
                    label = "Transposable Elements",
                    choices = unique(corrRef$teName)
                ),
                actionButton(inputId = "myButton", label = "Submit ")
            )
        ),
        grid_card(
            card_header("Distribution of Gene:TE in reference"),
            area = "area2",
            card_body(plotlyOutput(outputId = "scatterPlotRef"))
        ),
        grid_card(
            card_header("Distribution of Gene:TE in comparison"),
            area = "area3",
            card_body(plotlyOutput(outputId = "scatterPlotCompare"))
        ),
        grid_card(
            card_header("Reference Correlation"),
            area = "area4",
            card_body(plotOutput(outputId = "coefPlotRef"))
        ),
        grid_card(
            card_header("Comparison Correlation"),
            area = "area5",
            card_body(plotOutput(outputId = "coefPlotCompare"))
        ),
        grid_card(
            card_header("Gene expression (Log normalized)"),
            area = "area6",
            card_body(plotOutput(outputId = "deGenePlot"))
        ),
        grid_card(
            card_header("TE expression (Log normalized)"),
            area = "area7",
            card_body(plotOutput(outputId = "deTEPlot"))
        )
        
        
    )
    
    server <- function(input, output) {
        
        
        output$scatterPlotRef <- renderPlotly({
            corrRef_sig <- corrRef[corrRef$pvalue < 0.05, ]
            corrRef_sig$pair <- paste0(corrRef_sig$geneName, " : ", corrRef_sig$teName)
            
            plot_ly(data = corrRef_sig, x = ~pvalue, y = ~coef, 
                    text = ~pair, 
                    marker = list(size=10, color="#45A9EC", line=list(color="black", width=1)),
                    type = "scatter",
                    mode = "markers")
            
        })
        
        
        output$scatterPlotCompare <- renderPlotly({
            corrCompare_sig <- corrCompare[corrCompare$pvalue < 0.05, ]
            corrCompare_sig$pair <- paste0(corrCompare_sig$geneName, " : ", corrCompare_sig$teName)
            
            plot_ly(data = corrCompare_sig, x = ~pvalue, y = ~coef, 
                    text = ~pair, 
                    marker = list(size=10, color="#8BE748", line=list(color="black", width=1)),
                    type = "scatter", 
                    mode = "markers")
            
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
                pvalue <- sprintf("%.2e", pvalue)
                
                ggplot(df_ref, aes(x=gene, y=TE)) +
                    geom_point(colour="black", shape=21, size=3, fill="#45A9EC") +
                    labs(x=ensemblID, y=teName) +
                    geom_smooth(method = "lm") +
                    theme_bw() +
                    ggtitle(paste0("Coefficient: ", coef, "\npvalue: ", pvalue))
            })
            
            # render compare coef plot
            output$coefPlotCompare <- renderPlot({
                
                coef <- corrCompare[corrCompare$geneName==ensemblID & corrCompare$teName==teName, "coef"]
                coef <- round(coef, 4)
                pvalue <- corrCompare[corrCompare$geneName==ensemblID & corrCompare$teName==teName, "pvalue"]
                pvalue <- sprintf("%.2e", pvalue)
                
                ggplot(df_compare, aes(x=gene, y=TE)) +
                    geom_point(colour="black", shape=21, size=3, fill="#8BE748") +
                    labs(x=ensemblID, y=teName) +
                    geom_smooth(method = "lm") +
                    theme_bw() +
                    ggtitle(paste0("Coefficient: ", coef, "\npvalue: ", pvalue))
            })
            
            # render expression plot
            output$deGenePlot <- renderPlot({
                
                lfc2 <- round(geneLFC[ensemblID, "log2FoldChange"], 4)
                pvalue <- sprintf("%.2e", geneLFC[ensemblID, "pvalue"])
                
                ggviolin(df_all, x="group", y="gene", fill="group",
                         palette=c("#45A9EC", "#8BE748"), add="boxplot", 
                         add.params=list(fill="white")) +
                    ylab(ensemblID) +
                    xlab("") +
                    ggtitle(paste0("LFC2: ", lfc2, "\npvalue: ", pvalue)) +
                    theme(legend.position = "none")
                
            })
            
            # render TE plot
            output$deTEPlot <- renderPlot({
                
                lfc2 <- round(teLFC[teName, "log2FoldChange"], 4)
                pvalue <- sprintf("%.2e", teLFC[teName, "pvalue"])
                
                ggviolin(df_all, x="group", y="TE", fill="group",
                         palette=c("#45A9EC", "#8BE748"), add="boxplot", 
                         add.params=list(fill="white")) +
                    ylab(teName) +
                    xlab("") +
                    ggtitle(paste0("LFC2: ", lfc2, "\npvalue: ", pvalue)) +
                    theme(legend.position = "none")
                    
            })
            
        })
        
    }
    
    shinyApp(ui, server)
}
