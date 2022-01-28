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