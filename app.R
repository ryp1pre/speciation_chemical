#
# Speciation Of Polyprotic Acids And Bases.
# With ggplot plot wrapped in Plotly.
# Calculates solution pH and concentration of species.
#

library(shiny)
library(shinydashboard)
library(tidyr)
library(ggplot2)
library(plotly)
library(DT)


###########################################################################
#                                                                         #
#                       FUNCTIONS DEFINITIONS                             #
#                                                                         #
###########################################################################
# Speciation Function -----------------------------------------------------
speciate <- function(pH = seq(0, 14, 0.2), pKa1=NA, pKa2=NA, pKa3=NA, pKa4=NA) {
  # The function requres the 'tidyr' library.
  # The function calculates the speciation of 
  # polyprotic acids/bases with up to 4 pKa's.
  # The output is a list of two dataframes
  # (1st-wide; 2nd-tall) containing the pH
  # and the mol fraction of the ion species.
  
  pKa <- c(pKa1,pKa2,pKa3,pKa4)
  number_na <- sum(is.na(pKa))
  
  if (number_na == 3) {
    # pKa1 provided
    Ka1 <- 10**(-pKa1)
    D <- (10**(-pH)) + Ka1
    
    HA <-  (10**(-pH))/D
    A <-   Ka1/D
    
    dfr <- data.frame(pH, HA, A)
    dfr <- round(dfr,3)
    dfr_tall <- gather(dfr,"ion", "fraq", HA:A)
    dfr_list <- list(dfr, dfr_tall)
    return(dfr_list)
    
  } else if (number_na == 2) {
    # pKa's 1-2 provided
    Ka1 <- 10**(-pKa1)
    Ka2 <- 10**(-pKa2)
    D <- ((10**(-pH))**2 + ((10**(-pH))*Ka1) + (Ka1*Ka2))
    
    H2A <- ((10**(-pH))**2)/D
    HA <-  (10**(-pH)*Ka1)/D
    A <-   (Ka1*Ka2)/D
    
    dfr <- data.frame(pH, H2A, HA, A)
    dfr <- round(dfr,3)
    dfr_tall <- gather(dfr,"ion", "fraq", H2A:A)
    dfr_list <- list(dfr, dfr_tall)
    return(dfr_list)
    
  } else if (number_na == 1) {
    # pKa's 1-3 provided
    Ka1 <- 10**(-pKa1)
    Ka2 <- 10**(-pKa2)
    Ka3 <- 10**(-pKa3)
    D <- ((10**(-pH))**3 + ((10**(-pH))**2)*Ka1) + (10**(-pH))*Ka1*Ka2 + (Ka1*Ka2*Ka3)
    
    H3A <- ((10**(-pH))**3)/D
    H2A <- (((10**(-pH))**2)*Ka1)/D
    HA <-  (10**(-pH)*Ka1*Ka2)/D
    A <-   (Ka1*Ka2*Ka3)/D
    
    dfr <- data.frame(pH, H3A, H2A, HA, A)
    dfr <- round(dfr,3)
    dfr_tall <- gather(dfr,"ion", "fraq", H3A:A)
    dfr_list <- list(dfr, dfr_tall)
    return(dfr_list)
    
  } else {
    # all pKa's (1-4) provided
    Ka1 <- 10**(-pKa1)
    Ka2 <- 10**(-pKa2)
    Ka3 <- 10**(-pKa3)
    Ka4 <- 10**(-pKa4)
    D <- ((10**(-pH))**4 + ((10**(-pH))**3)*Ka1) + ((10**(-pH))**2)*Ka1*Ka2 + (10**(-pH))*Ka1*Ka2*Ka3 + (Ka1*Ka2*Ka3*Ka4)
    
    H4A <- ((10**(-pH))**4)/D
    H3A <- ((10**(-pH))**3)*Ka1/D
    H2A <- (((10**(-pH))**2)*Ka1*Ka2)/D
    HA <-  (10**(-pH)*Ka1*Ka2*Ka3)/D
    A <-   (Ka1*Ka2*Ka3*Ka4)/D
    
    dfr <- data.frame(pH, H4A, H3A, H2A, HA, A)
    dfr <- round(dfr,3)
    dfr_tall <- gather(dfr,"ion", "fraq", H4A:A)
    dfr_list <- list(dfr, dfr_tall)
    return(dfr_list)
  }
} # end of Speciation Function ----------------------------------------------


# pH Function -------------------------------------------------------------
ph_solution <- function(C, pKa1=NA, pKa2=NA, pKa3=NA, pKa4=NA) {
  # The function calculates the solution pH of 
  # polyprotic acids/bases with up to 4 pKa's.
  # The 'C' argument is the molar concentration.
  
  pKa <- c(pKa1,pKa2,pKa3,pKa4)
  number_na <- sum(is.na(pKa))
  
  
  if (number_na == 3) {
    # Mono-protic acid 
    Ka1 <- 10**(-pKa1)
    
    H <- (-Ka1 + sqrt(Ka1^2 + 4*Ka1*C))/2
    pH <- -log10(H)
    
  } else if (number_na == 2) {
    # Di-protic acid 
    Ka1 <- 10**(-pKa1)
    Ka2 <- 10**(-pKa2)
    
    H.step1 <- (-Ka1 + (Ka1^2 + 4*Ka1*C)^(1/2))/2
    H.step2 <- 0.5*(sqrt(H.step1^2 + 6*H.step1*Ka2 + Ka2^2) - H.step1 - Ka2)
    pH <- -log10(H.step1 + H.step2)
    
  } else if (number_na == 1) {
    # Tri-protic acid 
    Ka1 <- 10**(-pKa1)
    Ka2 <- 10**(-pKa2)
    Ka3 <- 10**(-pKa3)
    
    H.step1 <- (-Ka1 + (Ka1^2 + 4*Ka1*C)^(1/2))/2
    H.step2 <- 0.5*(sqrt(H.step1^2 + 6*H.step1*Ka2 + Ka2^2) - H.step1 - Ka2)
    H.step3 <- 0.5*(sqrt((- H.step1 - H.step2 - Ka3)^2 + 4*H.step2*Ka3) - H.step1 - H.step2 - Ka3)
    pH <- -log10(H.step1 + H.step2 + H.step3)
    
  } else {
    # Tetra-protic acid 
    Ka1 <- 10**(-pKa1)
    Ka2 <- 10**(-pKa2)
    Ka3 <- 10**(-pKa3)
    Ka4 <- 10**(-pKa4)
    
    H.step1 <- (-Ka1 + (Ka1^2 + 4*Ka1*C)^(1/2))/2
    H.step2 <- 0.5*(sqrt(H.step1^2 + 6*H.step1*Ka2 + Ka2^2) - H.step1 - Ka2)
    H.step3 <- 0.5*(sqrt(4*H.step2*Ka3 + (H.step1 + H.step2)^2 + 2*(H.step1 + H.step2)*Ka3 + Ka3^2) - H.step1 + H.step2 - Ka3)
    H.step4 <- 0.5*(sqrt((- H.step1 - H.step2 - H.step3 - Ka4)^2 + 4*H.step3*Ka4) - H.step1 - H.step2 - H.step3 - Ka4)
    pH <- -log10(H.step1 + H.step2 + H.step3 + H.step4)
  }
} # end of pH Function
# End of Function Definitions ///////////////////////////////////////////////

###########################################################################
#                                                                         #
#                               MAIN                                      #
#                                                                         #
###########################################################################

# HEADER
header <- dashboardHeader(title = "Speciation Of Polyprotic Acids And Bases", titleWidth = 400, disable = FALSE)

# SIDEBAR
sidebar <- dashboardSidebar(disable = FALSE,
                            h4("Enter up to four pKa values:"),
                            numericInput(inputId = "pka1", label = "pKa1:", value = ""),
                            numericInput(inputId = "pka2", label = "pKa2:", value = ""),
                            numericInput(inputId = "pka3", label = "pKa3:", value = ""),
                            numericInput(inputId = "pka4", label = "pKa4:", value = ""),
                            numericInput(inputId = "conc", label = "Total concentration [mg/mL]:*", value = ""),
                            helpText("* Optional: Only required for species concentration!"),
                            numericInput(inputId = "Mw", label = "Molecular weight:*", value = ""),
                            helpText("* Optional: Only required for solution pH calculation! Total concentration must also be provided!"),
                            actionButton(inputId = "calc", label = "Calculate")
)

# BODY
body <- dashboardBody(
  # add custom CSS to make the title background area the same
  # color as the rest of the header.
  tags$head(tags$style(HTML('
        .skin-blue .main-header .logo {
        background-color: #3c8dbc;
        }
        .skin-blue .main-header .logo:hover {
        background-color: #3c8dbc;
        }
        '))),
  
  
  fluidRow(
    box(width = 10, title = "Graph", status = "primary",solidHeader = TRUE, collapsible = TRUE,
        plotlyOutput("distPlot"),
        textOutput("ph"),
        tags$head(tags$style("#ph{color: red;
                            font-size: 20px;
                            font-style: italic;
                             "
                            )
                 )
    ),
    
    box(width = 10, title = "Mol Fraction Table", status = "primary",solidHeader = TRUE, collapsible = TRUE,
        DT::dataTableOutput("table_fraq")
    ),
    
    box(width = 10, title = "Concentration Table", status = "primary",solidHeader = TRUE, collapsible = TRUE,
      DT::dataTableOutput("table_conc")
    )
    
  )
)
  

ui <- dashboardPage(header, sidebar, body, 
      skin = "blue", tags$head(tags$style(HTML("
     .skin-blue .main-sidebar {background-color:  black;}")))) # CHANGES COLOR OF SIDE PANEL

server <- function(input, output) {
  observeEvent(input$calc, {
    
    pKa1 <- as.numeric(input$pka1)
    pKa2 <- as.numeric(input$pka2)
    pKa3 <- as.numeric(input$pka3)
    pKa4 <- as.numeric(input$pka4)
    

# Speciation function call ------------------------------------------------
  dfr_list <- speciate(pKa1=pKa1, pKa2=pKa2, pKa3=pKa3, pKa4=pKa4)

# Fraction table ----------------------------------------------------------
    output$table_fraq <- DT::renderDataTable(
      datatable(dfr_list[[1]], extensions = 'Buttons', 
                class="cell-border stripe", 
                rownames = FALSE, 
                options = list(dom = "Blfrtip",
                               buttond = list("copy", list(extend = "collection",
                                                           buttons = c("csv", "excel", "pdf"), 
                                                           text = "Download")), pageLength=10, autoWidth = TRUE, searchHighlight = TRUE, filter = "top")))

# Plot data ---------------------------------------------------------------
    output$distPlot <- renderPlotly({ggplotly
      (ggplot(data = dfr_list[[2]], aes(x = pH, y = fraq, color = ion)) +
        geom_line(size=1.5) +
          labs(title = "Speciation Plot", x="pH", y="Mol Fraction") +
          scale_x_continuous(breaks = round(seq(min(dfr_list[[2]]$pH), max(dfr_list[[2]]$pH), by = 1),1))
      )
    })

# Convert fraction to concentration ---------------------------------------
      total_conc <- as.numeric(input$conc)
      pKa <- c(pKa1,pKa2,pKa3,pKa4)
      number_na <- sum(is.na(pKa))

      if (number_na == 3) {
        Conc_HA <- dfr_list[[1]]$HA*total_conc
        Conc_A <- dfr_list[[1]]$A*total_conc
        pH <- dfr_list[[1]]$pH
        dfr_conc <- data.frame(pH,Conc_HA,Conc_A)
      } else if (number_na == 2) {
        Conc_H2A <- dfr_list[[1]]$H2A*total_conc
        Conc_HA <- dfr_list[[1]]$HA*total_conc
        Conc_A <- dfr_list[[1]]$A*total_conc
        pH <- dfr_list[[1]]$pH
        dfr_conc <- data.frame(pH,Conc_H2A,Conc_HA,Conc_A)
      } else if (number_na == 1) {
        Conc_H3A <- dfr_list[[1]]$H3A*total_conc
        Conc_H2A <- dfr_list[[1]]$H2A*total_conc
        Conc_HA <- dfr_list[[1]]$HA*total_conc
        Conc_A <- dfr_list[[1]]$A*total_conc
        pH <- dfr_list[[1]]$pH
        dfr_conc <- data.frame(pH,Conc_H3A,Conc_H2A,Conc_HA,Conc_A)
      } else {
        Conc_H4A <- dfr_list[[1]]$H4A*total_conc
        Conc_H3A <- dfr_list[[1]]$H3A*total_conc
        Conc_H2A <- dfr_list[[1]]$H2A*total_conc
        Conc_HA <- dfr_list[[1]]$HA*total_conc
        Conc_A <- dfr_list[[1]]$A*total_conc
        pH <- dfr_list[[1]]$pH
        dfr_conc <- data.frame(pH,Conc_H4A, Conc_H3A,Conc_H2A,Conc_HA,Conc_A)
      } # end of "if"

# Concentration table -----------------------------------------------------
      output$table_conc <- DT::renderDataTable(
        datatable(dfr_conc, extensions = 'Buttons',
                  class="cell-border stripe",
                  rownames = FALSE,
                  options = list(dom = "Blfrtip",
                  buttond = list("copy", list(extend = "collection",
                  buttons = c("csv", "excel", "pdf"),
                  text = "Download")), pageLength=10, autoWidth = TRUE,
                  searchHighlight = TRUE, filter = "top")))
      
# Solution pH -------------------------------------------------------------      
      output$ph <- renderText({
        C <- isolate(input$conc/input$Mw)
        # pH function call ------------------------------------------------
        pH <- ph_solution(C=C, pKa1=pKa1, pKa2=pKa2, pKa3=pKa3, pKa4=pKa4)
        paste("Solution pH = ", round(pH, 3))
      })
    
  }) # end of observeEvent for input$calc 
}
shinyApp(ui, server)

