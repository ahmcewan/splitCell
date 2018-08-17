library(shiny)
library(reshape2)
library(xlsx)
library(DT)
library(dplyr)
library(readr)
options(shiny.sanitize.errors = FALSE)
#Table used to select control wells
toggleTable <- matrix(" ", nrow = 8, ncol = 12, dimnames = list(c(LETTERS[1:8]), seq.int(1, 12, 1)))


localMaxima <- function(x) {
  b <- diff(c(-Inf, x)) > 0L
  rle(b)$lengths
  b <- cumsum(rle(b)$lengths)
  b <- b[seq.int(1L, length(b), 2L)]
  if (x[[1]] == x[[2]]) {
    b <- b[-1]
  }
  return(b)
}
#Table used to select control wells
toggleTable <- matrix(" ", nrow = 8, ncol = 12, dimnames = list(c(LETTERS[1:8]), seq.int(1, 12, 1)))


localMaxima <- function(x) {
  b <- diff(c(-Inf, x)) > 0L
  rle(b)$lengths
  b <- cumsum(rle(b)$lengths)
  b <- b[seq.int(1L, length(b), 2L)]
  if (x[[1]] == x[[2]]) {
    b <- b[-1]
  }
  return(b)
}

localMinima <- function(x) {
  a <- diff(c(Inf, x)) > 0L
  rle(a)$lengths
  a <- cumsum(rle(a)$lengths)
  a <- a[seq.int(1L, length(a), 2L)]
  if (x[[1]] == x[[2]]) {
    a <- a[-1]
  }
  return(a)
}

#ui


ui <- shinyUI(fluidPage(
  mainPanel(
    fileInput("files",label = h4("Upload Data"), multiple = TRUE, accept = c(".csv")),
    fluidRow(column(3, uiOutput("edu")),
             column(3, uiOutput("dapi")),
             column(3, uiOutput("foci")),
             column(3, uiOutput("dapi_int"))),
    DT::dataTableOutput("userChoiceTbl", width = "50%"),
    tags$b("Cells Selected:"),
    verbatimTextOutput("selectedInfo"), 
    fluidRow(column(6,h4("Density Plot of Mean Edu Intensity"),plotOutput("edu_plot")),
             column(6,h4("Density Plot of Sum Dapi Intensity"),plotOutput("dapi_plot"))),
    fluidRow(column(6,uiOutput("slider_edu")),
             column(2,uiOutput("slider_sub")),
             column(2,uiOutput("slider_dapi")),
             column(2,uiOutput("slider_cut"))),
    fluidRow(column(12,h4("Proportion of Cells",tableOutput("population_stats")))),
    downloadButton("downloadData", "Download")
    
  )
))


#server
#load data 
server <- function(input, output, session){
  dataset <- reactive({
    validate(
      need(input$files != "", "")
    )
    path_list <- as.list(input$files$datapath)
    tbl_list <- lapply(input$files$datapath, read.table, header=TRUE, sep=",")
    df <- do.call(rbind, tbl_list)
    df <- df[order(df$WellName),]
    col <- colnames(df)
    updateSelectInput(session, inputId ="edu", choices =c("select" , col ))
    updateSelectInput(session, inputId ="dapi", choices =c("select" , col ))
    updateSelectInput(session, inputId ="foci", choices =c("select" , col ))
    updateSelectInput(session, inputId ="dapi_int", choices =c("select" , col ))
    list(df=data.frame(df),
         path_list=path_list)
    
  })
  
  #update selectInput wells
  output$edu <- renderUI({
    selectInput("edu","Edu Mean",choices = as.vector(dataset()$col))
  })
  output$dapi <- renderUI({
    selectInput("dapi","Dapi sum",choices = as.vector(dataset()$col))
  })
  output$foci <- renderUI({
    selectInput("foci","Num spots",choices = as.vector(dataset()$col))
  })
  output$dapi_int <- renderUI({
    selectInput("dapi_int","H2AX",choices = as.vector(dataset()$col))
  })
  output$userChoiceTbl <- DT::renderDataTable({
    datatable(toggleTable,options = list(dom = 't',ordering = F),selection = list(target = 'cell'),class = 'cell-border compact') %>% formatStyle(1:12, cursor = 'pointer')
  })
  
  
  cont_wells <- reactive({
    validate(
      need(input$userChoiceTbl_cells_selected != "", "Please select control cells if applicable")
    )
    wells <- input$userChoiceTbl_cells_selected
    let <- LETTERS[wells[,1]]
    num <- sort(wells[,2])
    coordinates <- paste0(let,num)
    list(coordinates = coordinates)
  })
  
  #Outputs cells selected in datatable used to choose control wells  
  output$selectedInfo <- renderPrint({
    (cont_wells()$coordinates)
  })
  
  reshape_cols <- function(x, cols=cols){
    x$id <- with(x, ave(rep(1, nrow(x)), WellName, FUN = seq_along))
    x <- dcast(data = x,formula = id~WellName,fun.aggregate = sum,value.var =cols,fill=-1)
    x[x == -1] <-  NA
    x <- x[-1]
    dat <- x[,order(as.character(sub("[0-9]","",names(x))),as.integer(sub("[A-Z]", "", names(x))))]
    return(dat)
  }
  
  subset_file <- reactive({
    validate(
      need(input$edu != "", ""),
      need(input$dapi != "", ""),
      need(input$foci != "", ""),
      need(input$dapi_int != "", "")
    )
    if(input$edu %in% c("","select") | input$dapi %in% c("","select") | input$foci %in% c("","select") | input$dapi_int %in% c("","select")){
      return(NULL)
    }else{
      rawData = dataset()$df
      dat.cols <- c("WellName", input$edu, input$dapi, input$foci, input$dapi_int)
      dat.all <- dataset()$df[,colnames(dataset()$df) %in% dat.cols]
      dat.cont <- dat.all
      dat.all.count <- reshape_cols(dat.cont, cols=input$dapi_int)
      if (length(input$userChoiceTbl_cells_selected)!=0){
        dat.cont <- dat.cont[rawData$WellName %in% cont_wells()$coordinates,]
      }
    }
    list(dat.cont=dat.cont, 
         dat.all=dat.all,
         dat.all.count=dat.all.count)
  })
  
  dat_g1 <- reactive({
    validate(
      need(input$edu != "", ""),
      need(input$dapi != "", ""),
      need(input$foci != "", ""),
      need(input$dapi_int != "", "")
    )
    if (input$edu %in% c("","select") | input$dapi %in% c("","select") | input$foci %in% c("","select")| input$dapi_int %in% c("","select")){
      return(NULL)
    }else{
    if(nrow(subset_file()$dat.cont)==0){
    edu_values <- subset_file()$dat.all[[input$edu]]
    }else{
    edu_values <- subset_file()$dat.cont[[input$edu]]
    }
    edu_dens <- density(edu_values, na.rm=TRUE)
    dens <- edu_dens$y
    #add 'max' localMaxima?
    edu.max.x <-localMaxima(dens)[1]
    edu.min.x <-localMinima(dens)
    edu.min.thres <-edu.min.x[min(which(edu.min.x > edu.max.x))]
    edu.min.x <- edu_dens$x[edu.min.thres]
    

    list(edu_values=edu_values, 
         dens=dens,
         edu_dens=edu_dens,
         edu.min.x=edu.min.x,
         edu.max.x=edu.max.x,
         edu.min.thres=edu.min.thres)
}
  })
  
  output$slider_edu <- renderUI({
    validate(
      need(input$edu != "", ""),
      need(input$dapi != "", ""),
      need(input$foci != "", ""),
      need(input$dapi_int != "", "")
    )
    if (input$edu %in% c("","select") | input$dapi %in% c("","select") | input$foci %in% c("","select")| input$dapi_int %in% c("","select")){
      return(NULL)
    }else{
    sliderInput("eduSlider", "Edu Threshold", min=0, max=150, value=dat_g1()$edu.min.thres)
    }
  })
  
  dat_edu <- reactive({
    if (input$edu %in% c("","select") | input$dapi %in% c("","select") | input$foci %in% c("","select")| input$dapi_int %in% c("","select")){
      return(NULL)
    }else{
    validate(
      need(input$edu != "", ""),
      need(input$dapi != "", ""),
      need(input$foci != "", ""),
      need(input$dapi_int != "", "")
    )
      if(nrow(subset_file()$dat.cont)==0){
        dat.cont <- subset_file()$dat.all
      }else{
        dat.cont <- subset_file()$dat.cont
      }  
    #create edu neg and edu positive populations
    edu_threshold <- dat_g1()$edu_dens$x[input$eduSlider]
    edu_neg_pooled_dapi <- dat.cont[dat_g1()$edu_values < edu_threshold,]
    edu_pos_pooled_dapi <- dat.cont[dat_g1()$edu_values > edu_threshold,]
    
    #create edu neg and edu positive populations of all data based on control thresholds
    edu_neg_pooled_dapi_all <- subset_file()$dat.all[dat_g1()$edu_values < edu_threshold,]
    edu_pos_pooled_dapi_all <- subset_file()$dat.all[dat_g1()$edu_values >edu_threshold,]
    
    list(edu_neg_pooled_dapi=edu_neg_pooled_dapi, 
         edu_pos_pooled_dapi=edu_pos_pooled_dapi,
         edu_neg_pooled_dapi_all=edu_neg_pooled_dapi_all,
         edu_pos_pooled_dapi_all=edu_pos_pooled_dapi_all)
    }
  })
  
  dapi_threshold <- reactive({
    if (input$edu %in% c("","select") | input$dapi %in% c("","select") | input$foci %in% c("","select")| input$dapi_int %in% c("","select")){
      return(NULL)
    }else{
    validate(
      need(input$edu != "", ""),
      need(input$dapi != "", ""),
      need(input$foci != "", ""),
      need(input$dapi_int != "", "")
    )
    #get G1 and S/G2 threshold using dapi measurments
    edu_neg <- dat_edu()$edu_neg_pooled_dapi
    dapi_dens <- density(edu_neg[,input$dapi],na.rm=TRUE)
    maxPeaks <- localMaxima(dapi_dens$y)
    minValley <- localMinima(dapi_dens$y)
    peak1 <- maxPeaks[match(max(dapi_dens$y[maxPeaks]),dapi_dens$y[maxPeaks])]
    peak2 <- maxPeaks[match(max(dapi_dens$y[maxPeaks]),dapi_dens$y[maxPeaks]) + 1]
    valley <- minValley[minValley < peak2 & minValley > peak1]
    valley0 <- peak1 - (valley-peak1)
    cutoff <- peak2 + (peak2-valley)
    
    list(maxPeaks=maxPeaks,
         dapi_dens=dapi_dens,
         peak1=peak1,
         peak2=peak2, 
         minValley=minValley,
         valley=valley,
         valley0=valley0,
         cutoff=cutoff)
    }
  })
  
  output$slider_sub <- renderUI({
    validate(
      need(input$edu != "", ""),
      need(input$dapi != "", ""),
      need(input$foci != "", ""),
      need(input$dapi_int != "", "")
    )
    if (input$edu %in% c("","select") | input$dapi %in% c("","select") | input$foci %in% c("","select")| input$dapi_int %in% c("","select")){
      return(NULL)
    }else{
      validate(
        need(input$edu != "", ""),
        need(input$dapi != "", ""),
        need(input$foci != "", ""),
        need(input$dapi_int != "", "")
      )
      sliderInput("dapiSub", "SubG1 Threshold", min=dapi_threshold()$valley0-50, max=dapi_threshold()$peak1, value=dapi_threshold()$valley0)
    }
  })
  
  output$slider_dapi <- renderUI({
    validate(
      need(input$edu != "", ""),
      need(input$dapi != "", ""),
      need(input$foci != "", ""),
      need(input$dapi_int != "", "")
    )
    if (input$edu %in% c("","select") | input$dapi %in% c("","select") | input$foci %in% c("","select")| input$dapi_int %in% c("","select")){
      return(NULL)
    }else{
    validate(
      need(input$edu != "", ""),
      need(input$dapi != "", ""),
      need(input$foci != "", ""),
      need(input$dapi_int != "", "")
    )
    sliderInput("dapiSlider", "Dapi Threshold", min=dapi_threshold()$peak1, max=dapi_threshold()$peak2, value=dapi_threshold()$valley)
    }
  })
  
  
  output$slider_cut <- renderUI({
    validate(
      need(input$edu != "", ""),
      need(input$dapi != "", ""),
      need(input$foci != "", ""),
      need(input$dapi_int != "", "")
    )
    if (input$edu %in% c("","select") | input$dapi %in% c("","select") | input$foci %in% c("","select")| input$dapi_int %in% c("","select")){
      return(NULL)
    }else{
    sliderInput("cutSlider", "Dapi cutoff", min=dapi_threshold()$peak1+50, max=dapi_threshold()$peak2+100, value=dapi_threshold()$cutoff)
    }
  })
  
  
  downloads <- reactive({
    validate(
      need(input$edu != "", ""),
      need(input$dapi != "", ""),
      need(input$foci != "", ""),
      need(input$dapi_int != "", "")
    )
    if (input$edu %in% c("","select") | input$dapi %in% c("","select") | input$foci %in% c("","select")| input$dapi_int %in% c("","select")){
      return(NULL)
    }else{
    #get foci of edu + cells--edu pos
    edu_pos_dapi <- dat_edu()$edu_pos_pooled_dapi_all[colnames(dat_edu()$edu_pos_pooled_dapi_all) %in% c("WellName", input$dapi_int)]
    edu_pos_int <- reshape_cols(edu_pos_dapi,cols=input$dapi_int)
    #edu neg
    edu_neg_dapi <- dat_edu()$edu_neg_pooled_dapi_all[colnames(dat_edu()$edu_neg_pooled_dapi_all) %in% c("WellName", input$dapi_int)]
    edu_neg_int <- reshape_cols(edu_neg_dapi,cols=input$dapi_int)
    
    #split data based on G1 and S/G2 threshold
    cell_cycle_threshold <- dapi_threshold()$dapi_dens$x[input$dapiSlider]
    sub_g1_threshold <- dapi_threshold()$dapi_dens$x[input$dapiSub]
    dapi_cutoff <- dapi_threshold()$dapi_dens$x[input$cutSlider]
    dapi_values <- dat_edu()$edu_neg_pooled_dapi[[input$dapi]]
    dapi_values <- dapi_values[dapi_values < dapi_cutoff]
    dat.sub <- dat_edu()$edu_neg_pooled_dapi_all[dapi_values < sub_g1_threshold,]
    dat.g1 <- dat_edu()$edu_neg_pooled_dapi_all[dapi_values < cell_cycle_threshold & dapi_values > sub_g1_threshold,]
    dat.S <- dat_edu()$edu_neg_pooled_dapi_all[dapi_values > cell_cycle_threshold,]
    edu_neg_subG1_int <- reshape_cols(dat.sub, cols=input$dapi_int)
    edu_neg_g1_int <- reshape_cols(dat.g1, cols=input$dapi_int)
    edu_neg_S_int <- reshape_cols(dat.S, cols=input$dapi_int)
    dat.g1 <- reshape_cols(dat.g1,cols=input$foci)
    dat.S <- reshape_cols(dat.S,cols=input$foci)
    
  
    #edu positive foci
    edu_pos_foci <- dat_edu()$edu_pos_pooled_dapi_all
    edu_pos_foci <- reshape_cols(edu_pos_foci,cols=input$foci)
    
    
    list(edu_neg_subG1_int=edu_neg_subG1_int,
         dat.g1=dat.g1, 
         dat.S=dat.S,
         edu_pos_int=edu_pos_int,
         edu_neg_int=edu_neg_int,
         edu_pos_foci=edu_pos_foci,
         edu_neg_g1_int=edu_neg_g1_int,
         edu_neg_S_int=edu_neg_S_int)
    }
  })
  
  populations <- reactive({
    validate(
      need(input$edu != "", ""),
      need(input$dapi != "", ""),
      need(input$foci != "", ""),
      need(input$dapi_int != "", "")
    )
    if (input$edu %in% c("","select") | input$dapi %in% c("","select") | input$foci %in% c("","select")| input$dapi_int %in% c("","select")){
      return(NULL)
    }else{
      if (length(dataset()$path_list) > 1){
        total_cells <- apply(subset_file()$dat.all.count, 2, function(x) length(na.omit(x)))
        count_subG1 <- apply(downloads()$edu_neg_subG1_int, 2, function(x) length(na.omit(x)))
        count_G1 <- apply(downloads()$dat.g1, 2, function(x) length(na.omit(x)))
        count_G2 <- apply(downloads()$dat.S, 2, function(x) length(na.omit(x)))
        count_edu_neg <- apply(downloads()$edu_neg_int, 2, function(x) length(na.omit(x)))
        count_edu_pos <- apply(downloads()$edu_pos_int, 2, function(x) length(na.omit(x)))
      }else{
        total_cells <- length(subset_file()$dat.all.count)
        count_subG1 <- length(downloads()$edu_neg_subG1_int)
        count_G1 <- length(downloads()$dat.g1)
        count_G2 <- length(downloads()$dat.S)
        count_edu_neg <- length(downloads()$edu_neg_subG1_int)
        count_edu_pos <- length(downloads()$count_edu_pos)
        
      }
      
      pop_table <- data.frame(rbind(count_subG1,count_G1, count_G2, count_edu_neg,count_edu_pos,total_cells))
      rownames(pop_table) <- c("subG1","G1", "G2","total edu-","total edu+","total")
      list(pop_table=pop_table)
    }
  })
  
  
  output$population_stats <- renderTable({
    rownames=TRUE
    validate(
      need(input$edu != "", ""),
      need(input$dapi != "", ""),
      need(input$foci != "", ""),
      need(input$dapi_int != "", "")
    )
    if (input$edu %in% c("","select") | input$dapi %in% c("","select") | input$foci %in% c("","select")| input$dapi_int %in% c("","select")){
      return(NULL)
    }else{
      populations()$pop_table
    }
  },
  rownames = TRUE)
  
  output$rawData <- renderDataTable({
    dat_edu()$edu_pos_pooled_dapi_all
  },options = list(pageLength = 10))
  
  output$edu_plot <- renderPlot({
    validate(
      need(input$edu != "", ""),
      need(input$dapi != "", ""),
      need(input$foci != "", ""),
      need(input$dapi_int != "", "")
    )
    if (input$edu %in% c("","select") | input$dapi %in% c("","select") | input$foci %in% c("","select")| input$dapi_int %in% c("","select")){
      return(NULL)
    }else{
      plot(dat_g1()$dens,type="l")
      abline(v=input$eduSlider, col="red")
    }
  })
  
  output$dapi_plot <- renderPlot({
    validate(
      need(input$edu != "", ""),
      need(input$dapi != "", ""),
      need(input$foci != "", ""),
      need(input$dapi_int != "", "")
    )
    if (input$edu %in% c("","select") | input$dapi %in% c("","select") | input$foci %in% c("","select")| input$dapi_int %in% c("","select")){
      return(NULL)
    }else{
      plot(dapi_threshold()$dapi_dens$y, type="l")
      abline(v=input$dapiSub, col="blue")
      abline(v=input$dapiSlider, col="red")
      abline(v=input$cutSlider, col="green")
    }
  })
  
  
  ###download data###
  ### remove #id column
  output$downloadData <- downloadHandler(
    filename = function(file) {
      paste("data-", Sys.Date(), ".xlsx", sep="")
    },
    content = function(con) {
      write.xlsx2(populations()$pop_table, con, sheetName="population_props", row.names=TRUE)
      write.xlsx2(downloads()$dat.g1, con, sheetName="edu_neg_G1_foci", row.names=FALSE, append = TRUE)
      write.xlsx2(downloads()$dat.S, con, sheetName="edu_neg_G2_foci", row.names=FALSE,append = TRUE)
      write.xlsx2(downloads()$edu_pos_foci, con, sheetName="edu_pos_foci", row.names=FALSE,append = TRUE)
      write.xlsx2(downloads()$edu_pos_int, con, sheetName="edu_pos_nuclear_int", row.names=FALSE, append = TRUE)
      write.xlsx2(downloads()$edu_neg_g1_int, con, sheetName="edu_neg_G1_int", row.names=FALSE,append = TRUE)
      write.xlsx2(downloads()$edu_neg_S_int, con, sheetName="edu_neg_G2_int", row.names=FALSE,append = TRUE)
    }
  )
  
  ###/download data###
  
  
}  


# Run the application 
shinyApp(ui = ui, server = server)